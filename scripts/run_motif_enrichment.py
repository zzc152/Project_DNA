"""MPRA 知识库 · 步骤 2：基序富集分析（FIMO + Fisher + BH-FDR）。

对每个细胞系:
    - 从 JASPAR2024 提取人类关键 TF 白名单子库
    - 高/低活性组序列转 fasta，用 fimo (p<1e-4) 扫描
    - Fisher 精确检验 + BH-FDR，筛选 enriched/depleted
    - 保存 data/processed/motif_enrichment_<CELL>.csv

用法（项目根目录）:
    python scripts/run_motif_enrichment.py
"""

from __future__ import annotations

import argparse
import csv
import logging
import sys
from pathlib import Path

import h5py
import numpy as np
import yaml

ROOT = Path(__file__).resolve().parent.parent
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from src.knowledge.motif_enrichment import (  # noqa: E402
    analyze_motifs,
    extract_whitelist_meme,
    get_whitelist_for,
    parse_jaspar_meme,
    sequences_to_fasta,
)
from src.knowledge.mpra_loader import onehot_to_seq  # noqa: E402

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s | %(levelname)s | %(message)s",
    datefmt="%H:%M:%S",
)
logger = logging.getLogger("run_motif_enrichment")

CSV_COLUMNS = [
    "motif_id", "motif_name", "odds_ratio", "p_value", "fdr",
    "n_high_hit", "n_low_hit", "frac_high", "frac_low", "direction",
]


def load_split_indices(outdir: Path, cell_line: str) -> tuple[np.ndarray, np.ndarray]:
    high_path = outdir / f"high_{cell_line}_indices.npy"
    low_path = outdir / f"low_{cell_line}_indices.npy"
    if not high_path.exists() or not low_path.exists():
        raise FileNotFoundError(f"缺少步骤1输出: {high_path} / {low_path}")
    return np.load(high_path), np.load(low_path)


def main() -> None:
    parser = argparse.ArgumentParser(description="MPRA 基序富集分析")
    parser.add_argument("--config", default="configs/default_config.yaml")
    parser.add_argument("--outdir", default="data/processed")
    parser.add_argument("--workdir", default="data/processed/fimo_work")
    parser.add_argument("--cells", default=None, help="逗号分隔的细胞系子集，如 K562")
    args = parser.parse_args()

    cfg_path = ROOT / args.config
    with open(cfg_path, "r", encoding="utf-8") as f:
        cfg = yaml.safe_load(f)
    mpra = cfg.get("mpra", {})
    if not mpra:
        logger.error("配置缺少 mpra 段")
        sys.exit(1)

    h5_path = mpra["h5_path"]
    x_key = mpra.get("x_key", "x_train")
    y_key = mpra.get("y_key", "y_train")
    cell_lines = mpra.get("cell_lines", ["K562", "HepG2", "SK-N-SH"])
    if args.cells:
        cell_lines = [c.strip() for c in args.cells.split(",") if c.strip()]
    fimo_bin = mpra["fimo"]
    jaspar_meme = mpra["jaspar_meme"]
    pvalue = mpra.get("fimo_pvalue", 1e-4)

    outdir = ROOT / args.outdir
    workdir = ROOT / args.workdir
    workdir.mkdir(parents=True, exist_ok=True)

    # 解析 JASPAR 名称映射
    motifs = parse_jaspar_meme(jaspar_meme)
    motif_names = {mid: m["name"] for mid, m in motifs.items()}

    # 加载 one-hot（内存约 0.5GB，安全）
    logger.info("加载 MPRA 数据: %s", h5_path)
    with h5py.File(h5_path, "r") as f:
        x = f[x_key][...]

    for cell in cell_lines:
        logger.info("############ %s ############", cell)
        high_idx, low_idx = load_split_indices(outdir, cell)

        # 1. 白名单子库
        wl_ids = get_whitelist_for(cell)
        motif_meme = workdir / f"motifs_{cell}.meme"
        written = extract_whitelist_meme(jaspar_meme, wl_ids, motif_meme)
        if written == 0:
            logger.error("%s 无白名单基序可写", cell)
            sys.exit(1)
        wl_ids = [m for m in wl_ids if m in motifs]

        # 2. 序列 -> fasta
        high_seqs = onehot_to_seq(x[high_idx])
        low_seqs = onehot_to_seq(x[low_idx])
        high_fasta = workdir / f"high_{cell}.fasta"
        low_fasta = workdir / f"low_{cell}.fasta"
        sequences_to_fasta(high_seqs, high_fasta, prefix="H")
        sequences_to_fasta(low_seqs, low_fasta, prefix="L")
        logger.info("%s: 高 %d 条 / 低 %d 条 fasta 已生成", cell, len(high_seqs), len(low_seqs))

        # 3. 富集分析
        rows = analyze_motifs(
            wl_ids, motif_names, motif_meme,
            high_fasta, low_fasta,
            fimo_bin, workdir / f"fimo_{cell}",
            n_high=len(high_seqs), n_low=len(low_seqs),
            pvalue=pvalue,
        )

        # 4. 保存 CSV
        csv_path = outdir / f"motif_enrichment_{cell}.csv"
        with open(csv_path, "w", newline="", encoding="utf-8") as f:
            writer = csv.DictWriter(f, fieldnames=CSV_COLUMNS)
            writer.writeheader()
            for r in rows:
                writer.writerow({k: r.get(k, "") for k in CSV_COLUMNS})
        logger.info("已保存: %s (%d 行)", csv_path, len(rows))

        # 5. 摘要
        enriched = [r for r in rows if r["direction"] == "enriched_in_high"]
        depleted = [r for r in rows if r["direction"] == "depleted_in_high"]
        logger.info("--- %s Top5 富集 (enriched_in_high) ---", cell)
        for r in sorted(enriched, key=lambda r: -r["odds_ratio"])[:5]:
            logger.info("  %-12s OR=%.2f p=%.2e fdr=%.2e  (高%.1f%% 低%.1f%%)",
                        r["motif_name"], r["odds_ratio"], r["p_value"], r["fdr"],
                        r["frac_high"] * 100, r["frac_low"] * 100)
        logger.info("--- %s Top5 耗尽 (depleted_in_high) ---", cell)
        for r in sorted(depleted, key=lambda r: r["odds_ratio"])[:5]:
            logger.info("  %-12s OR=%.2f p=%.2e fdr=%.2e  (高%.1f%% 低%.1f%%)",
                        r["motif_name"], r["odds_ratio"], r["p_value"], r["fdr"],
                        r["frac_high"] * 100, r["frac_low"] * 100)
        if not enriched and not depleted:
            logger.info("  (无显著富集/耗尽的基序)")

    logger.info("=== 全部完成 ===")


if __name__ == "__main__":
    main()
