"""基序富集分析：用 FIMO 扫描高/低活性组，统计 TF 基序富集/耗尽。

方法（对应步骤 2）:
    1. 从 JASPAR2024 vertebrates MEME 文件提取人类关键 TF 白名单子库
    2. 高/低活性组序列转 fasta，用 fimo (p-value < 1e-4) 扫描
    3. 对每个基序构造 2x2 列联表，Fisher 精确检验 + BH-FDR 校正
    4. 筛选 enriched_in_high (FDR<0.01, OR>1.5) / depleted_in_high (FDR<0.01, OR<1/1.5)
"""

from __future__ import annotations

import logging
import math
import re
import subprocess
from pathlib import Path
from typing import Iterable, Optional

import numpy as np

logger = logging.getLogger("motif_enrichment")

# 精选人类关键 TF 白名单（按细胞系分组；ID 为 JASPAR2024 中的 ID）
# K562 (红系/髓系) / HepG2 (肝) / SK-N-SH (神经) + 广泛表达的管家因子
MOTIF_WHITELIST: dict[str, list[str]] = {
    "K562": [
        "MA0035.5",   # GATA1
        "MA0036.4",   # GATA2
        "MA0140.3",   # GATA1::TAL1
        "MA0493.3",   # KLF1
        "MA0841.2",   # NFE2
        "MA0501.2",   # MAF::NFE2
        "MA0091.2",   # TAL1::TCF3
        "MA0002.3",   # Runx1
        "MA0080.7",   # Spi1 (PU.1)
        "MA0475.3",   # FLI1
        "MA0100.4",   # MYB
    ],
    "HepG2": [
        "MA0114.5",   # HNF4A
        "MA0046.3",   # HNF1A
        "MA0102.5",   # CEBPA
        "MA0466.4",   # CEBPB
        "MA0148.5",   # FOXA1
        "MA0047.4",   # FOXA2
        "MA0679.3",   # ONECUT1
    ],
    "SK-N-SH": [
        "MA0669.1",   # NEUROG2
        "MA1100.3",   # ASCL1
        "MA1109.2",   # NEUROD1
        "MA0138.3",   # REST
        "MA0442.3",   # SOX10
        "MA0069.1",   # PAX6
        "MA1648.2",   # TCF12
    ],
    "common": [
        "MA0079.5",   # SP1
        "MA0139.2",   # CTCF
        "MA0062.4",   # GABPA
        "MA0098.4",   # ETS1
        "MA0473.4",   # ELF1
        "MA0506.3",   # Nrf1
        "MA0150.3",   # Nfe2l2
        "MA0095.4",   # Yy1
    ],
}

FIMO_PVALUE = 1e-4


def get_whitelist_for(cell_line: str) -> list[str]:
    """返回某细胞系 + 公共因子的白名单 ID 列表。"""
    ids = list(MOTIF_WHITELIST.get(cell_line, []))
    ids += list(MOTIF_WHITELIST.get("common", []))
    # 去重保序
    seen = set()
    out = []
    for i in ids:
        if i not in seen:
            seen.add(i)
            out.append(i)
    return out


def parse_jaspar_meme(meme_path: str | Path) -> dict[str, dict]:
    """解析 JASPAR MEME v4 文件 -> {motif_id: {name, block}}。"""
    meme_path = Path(meme_path)
    motifs: dict[str, dict] = {}
    current_id: Optional[str] = None
    current_name: Optional[str] = None
    lines = meme_path.read_text(encoding="utf-8").splitlines()

    i = 0
    n = len(lines)
    while i < n:
        line = lines[i]
        if line.startswith("MOTIF "):
            parts = line.split()
            current_id = parts[1]
            current_name = parts[2] if len(parts) > 2 else current_id
            # 收集 block 直到下一个 MOTIF 行（含 URL 行）
            j = i + 1
            block = []
            while j < n and not lines[j].startswith("MOTIF "):
                block.append(lines[j])
                j += 1
            motifs[current_id] = {"name": current_name, "block": block}
            i = j
        else:
            i += 1
    return motifs


def extract_whitelist_meme(
    meme_path: str | Path,
    motif_ids: Iterable[str],
    out_path: str | Path,
) -> int:
    """从 JASPAR 文件提取白名单基序，写出 MEME v4 子文件。返回写出数量。"""
    motifs = parse_jaspar_meme(meme_path)
    id_set = set(motif_ids)
    missing = id_set - set(motifs.keys())
    if missing:
        logger.warning("白名单中不存在于 JASPAR 文件的 ID: %s", sorted(missing))

    out_path = Path(out_path)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with open(out_path, "w", encoding="utf-8") as f:
        f.write("MEME version 4\n\n")
        for mid in motif_ids:
            if mid in motifs:
                f.write(f"MOTIF {mid} {motifs[mid]['name']}\n")
                f.write("\n".join(motifs[mid]["block"]))
                f.write("\n\n")
    written = sum(1 for m in motif_ids if m in motifs)
    logger.info("已写出 %d 个白名单基序 -> %s", written, out_path)
    return written


def sequences_to_fasta(
    seqs: Iterable[str],
    out_path: str | Path,
    prefix: str = "seq",
) -> int:
    """序列列表 -> fasta 文件。返回条数。"""
    out_path = Path(out_path)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    count = 0
    with open(out_path, "w", encoding="utf-8") as f:
        for i, s in enumerate(seqs):
            f.write(f">{prefix}_{i}\n{s}\n")
            count += 1
    return count


def run_fimo(
    fimo_bin: str,
    motif_meme: str | Path,
    seq_fasta: str | Path,
    out_dir: str | Path,
    pvalue: float = FIMO_PVALUE,
    threads: int = 4,
) -> list[dict]:
    """运行 fimo 扫描，返回 [{motif_id, sequence_name}] 命中列表（p-value 阈值内）。"""
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    cmd = [
        fimo_bin,
        "--thresh", str(pvalue),
        "--text",
        "--verbosity", "1",
        # MPRA 片段方向随机，扫双链以保证基序可落在任一条链上
        str(motif_meme),
        str(seq_fasta),
    ]
    logger.info("运行 fimo: %s", " ".join(cmd))
    proc = subprocess.run(cmd, capture_output=True, text=True, timeout=3600)
    if proc.returncode != 0:
        raise RuntimeError(f"fimo 失败 rc={proc.returncode}\n{proc.stderr[-2000:]}")
    hits: list[dict] = []
    for line in proc.stdout.splitlines():
        if line.startswith("#") or line.startswith("motif_id"):
            continue
        parts = line.split("\t")
        # fimo --text 列: motif_id, motif_alt_id, sequence_name, start, stop, strand, score, p-value, q-value, matched_sequence
        if len(parts) < 10:
            continue
        hits.append({"motif_id": parts[0], "sequence_name": parts[2]})
    logger.info("fimo 完成: %d 条命中记录", len(hits))
    return hits


def bh_fdr(p_values: np.ndarray) -> np.ndarray:
    """Benjamini-Hochberg FDR 校正（返回与输入等长的 q-values）。"""
    p = np.asarray(p_values, dtype=np.float64)
    n = len(p)
    if n == 0:
        return np.array([], dtype=np.float64)
    order = np.argsort(p)
    ranked = p[order]
    q = ranked * n / np.arange(1, n + 1)
    # 单调性保证
    q = np.minimum.accumulate(q[::-1])[::-1]
    q = np.clip(q, 0.0, 1.0)
    out = np.empty_like(q)
    out[order] = q
    return out


def enrich_analysis(
    high_hits: set[str],
    low_hits: set[str],
    n_high: int,
    n_low: int,
) -> dict:
    """对单个基序做富集分析。

    2x2 列联表:
              有匹配  无匹配
        高组   a       b
        低组   c       d

    Returns:
        {odds_ratio, p_value, n_high_hit, n_low_hit, ...}
    """
    a = len(high_hits)  # 高组命中数（去重后的序列数）
    c = len(low_hits)
    b = n_high - a
    d = n_low - c
    table = [[a, b], [c, d]]
    # Fisher 精确检验（双侧）
    from scipy.stats import fisher_exact

    odds_ratio, p_value = fisher_exact(table, alternative="two-sided")
    return {
        "odds_ratio": float(odds_ratio),
        "p_value": float(p_value),
        "n_high_hit": int(a),
        "n_low_hit": int(c),
        "frac_high": float(a / n_high) if n_high else 0.0,
        "frac_low": float(c / n_low) if n_low else 0.0,
    }


def analyze_motifs(
    motif_ids: list[str],
    motif_names: dict[str, str],
    motif_meme: str | Path,
    high_fasta: str | Path,
    low_fasta: str | Path,
    fimo_bin: str,
    work_dir: str | Path,
    n_high: int,
    n_low: int,
    pvalue: float = FIMO_PVALUE,
) -> list[dict]:
    """对高/低组 fasta 各跑一次 fimo，统计每个基序的富集。

    Returns:
        rows: [{motif_id, motif_name, odds_ratio, p_value, fdr, direction, ...}]
    """
    work_dir = Path(work_dir)
    work_dir.mkdir(parents=True, exist_ok=True)

    high_hits = run_fimo(fimo_bin, motif_meme, high_fasta, work_dir / "fimo_high", pvalue)
    low_hits = run_fimo(fimo_bin, motif_meme, low_fasta, work_dir / "fimo_low", pvalue)

    # 按 motif 分组（sequence_name 无重复需求；只需集合判断）
    from collections import defaultdict

    high_by_motif: dict[str, set[str]] = defaultdict(set)
    low_by_motif: dict[str, set[str]] = defaultdict(set)
    for h in high_hits:
        high_by_motif[h["motif_id"]].add(h["sequence_name"])
    for h in low_hits:
        low_by_motif[h["motif_id"]].add(h["sequence_name"])

    rows: list[dict] = []
    for mid in motif_ids:
        if mid not in high_by_motif and mid not in low_by_motif:
            # 无任何命中 -> OR 无法计算，跳过
            continue
        res = enrich_analysis(high_by_motif.get(mid, set()), low_by_motif.get(mid, set()), n_high, n_low)
        rows.append({"motif_id": mid, "motif_name": motif_names.get(mid, mid), **res})

    if not rows:
        return rows

    # BH-FDR
    pvals = np.array([r["p_value"] for r in rows])
    fdrs = bh_fdr(pvals)
    for r, q in zip(rows, fdrs):
        r["fdr"] = float(q)

    # 方向
    for r in rows:
        if r["fdr"] < 0.01 and r["odds_ratio"] > 1.5:
            r["direction"] = "enriched_in_high"
        elif r["fdr"] < 0.01 and r["odds_ratio"] < 1 / 1.5:
            r["direction"] = "depleted_in_high"
        else:
            r["direction"] = "not_significant"

    rows.sort(key=lambda r: (-r["odds_ratio"], r["p_value"]))
    return rows
