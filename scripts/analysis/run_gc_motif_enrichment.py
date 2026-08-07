#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""run_gc_motif_enrichment.py — GC 含量与 motif 富集关系分析（新增 claim 层）。

核心思想（用户 2026-08-07 提出）:
    motif 有 GC 倾向（GC-rich motif / AT-rich motif），
    高 GC 序列更可能富集 GC-rich motif —— 这是可统计验证的 claim。

分析流程（每个细胞系独立）:
    1. 加载 MPRA one-hot 序列，计算每条序列 GC 含量
    2. 按 GC 分位数划分 GC-rich 组（GC >= p75）与 GC-poor 组（GC <= p25）
    3. 两组转 fasta，fimo 扫描白名单 motif
    4. 每个 motif 做 2x2 Fisher（GC-rich 组命中 vs GC-poor 组命中）+ BH-FDR
    5. 关联 motif 自身 PWM GC 特征，输出关系 claim

输出:
    data/processed/gc_motif_enrichment_<CELL>.csv   — 富集统计
    data/processed/gc_motif_claims_<CELL>.jsonl     — 可入库 claim

用法（项目根目录）:
    python scripts/analysis/run_gc_motif_enrichment.py
"""

from __future__ import annotations

import argparse
import csv
import json
import logging
import sys
from pathlib import Path

import h5py
import numpy as np
import yaml

ROOT = Path(__file__).resolve().parents[2]  # 项目根（scripts/analysis/ -> 上两级）
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from src.knowledge.gc_analysis import gc_content_from_onehot  # noqa: E402
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
logger = logging.getLogger("gc_motif_enrichment")

CSV_COLUMNS = [
    "motif_id", "motif_name", "motif_gc", "gc_class",
    "odds_ratio", "p_value", "fdr",
    "n_high_hit", "n_low_hit", "frac_high", "frac_low", "direction",
]

# GC 组分位数
GC_RICH_Q = 0.75
GC_POOR_Q = 0.25
MAX_PER_GROUP = 20000
SEED = 42

CLAIM_CN = {
    "gc_rich": "GC-rich motif（PWM中GC比例{gc:.2f}）",
    "at_rich": "AT-rich motif（PWM中GC比例{gc:.2f}）",
    "neutral": "中性GC偏好motif（PWM中GC比例{gc:.2f}）",
}


def build_claim(row: dict, cell: str, n_total: int) -> dict:
    """由富集统计行生成可入库 claim。"""
    motif_gc = row["motif_gc"]
    cls = row["gc_class"]
    desc = CLAIM_CN.get(cls, "motif").format(gc=motif_gc)

    or_v = row["odds_ratio"]
    if row["direction"] == "enriched_in_high":  # 此处 high = GC-rich 组
        rel = "显著富集"
    elif row["direction"] == "depleted_in_high":
        rel = "显著减少（耗尽）"
    else:
        rel = "无显著差异"

    frac_rich = row.get("frac_high", 0.0)   # analyze_motifs 用 frac_high/frac_low 命名
    frac_poor = row.get("frac_low", 0.0)
    claim = (
        f"在{cell}的MPRA序列中，{row['motif_name']}（{row['motif_id']}，{desc}）"
        f"在GC含量高（≥{GC_RICH_Q:.0%}分位）的序列中{rel}"
        f"（OR={or_v:.2f}, FDR={row['fdr']:.2e}），"
        f"GC-rich组命中率{frac_rich*100:.1f}% vs GC-poor组{frac_poor*100:.1f}%，"
        f"提示{row['motif_name']}结合位点与序列GC含量相关。"
    )
    return {
        "claim": claim,
        "claim_type": "association",
        "entities": {
            "tf": row["motif_name"],
            "motif": row["motif_id"],
            "cell_line": cell,
            "gc_class": cls,
            "motif_gc": round(motif_gc, 4),
        },
        "evidence": {
            "source": "Gosai et al. 2022 MPRA + FIMO scan of JASPAR2024 motifs (p<1e-4), "
                      "GC quantile split, Fisher exact + BH-FDR",
            "direction": "enriched_in_gc_rich" if row["direction"] == "enriched_in_high"
                         else "depleted_in_gc_rich" if row["direction"] == "depleted_in_high"
                         else "not_significant",
            "p_value": row["p_value"],
            "effect_size": row["odds_ratio"],
        },
        "metadata": {
            "n_total": n_total,
            "gc_rich_q": GC_RICH_Q,
            "gc_poor_q": GC_POOR_Q,
        },
    }


def main() -> None:
    parser = argparse.ArgumentParser(description="GC 含量与 motif 富集关系分析")
    parser.add_argument("--config", default="configs/default_config.yaml")
    parser.add_argument("--outdir", default="data/processed")
    parser.add_argument("--workdir", default="data/processed/fimo_gc_work")
    parser.add_argument("--cells", default=None, help="逗号分隔细胞系，如 K562")
    args = parser.parse_args()

    with open(ROOT / args.config, "r", encoding="utf-8") as f:
        cfg = yaml.safe_load(f)
    mpra = cfg.get("mpra", {})
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

    motifs = parse_jaspar_meme(jaspar_meme)
    motif_names = {mid: m["name"] for mid, m in motifs.items()}

    # 关联 motif GC 特征（来自步骤2.5 输出）
    gc_feat = {}
    gc_path = outdir / "motif_gc_features.json"
    if gc_path.exists():
        for r in json.loads(gc_path.read_text(encoding="utf-8")):
            gc_feat[r["motif_id"]] = (r["motif_gc"], r["gc_class"])
        logger.info("加载 motif GC 特征 %d 条", len(gc_feat))

    logger.info("加载 MPRA one-hot: %s", h5_path)
    with h5py.File(h5_path, "r") as f:
        x = f[x_key][...]
    logger.info("x shape=%s", x.shape)

    for cell in cell_lines:
        logger.info("############ %s ############", cell)

        # 1. GC 含量 + 分组
        gc_frac = gc_content_from_onehot(x)
        gc_rich_q = np.quantile(gc_frac, GC_RICH_Q)
        gc_poor_q = np.quantile(gc_frac, GC_POOR_Q)
        gc_rich_idx = np.where(gc_frac >= gc_rich_q)[0]
        gc_poor_idx = np.where(gc_frac <= gc_poor_q)[0]
        logger.info(
            "GC 分组: rich>=%.3f (%d 条) / poor<=%.3f (%d 条)",
            gc_rich_q, len(gc_rich_idx), gc_poor_q, len(gc_poor_idx),
        )
        # 上限采样
        rng = np.random.default_rng(SEED)
        if len(gc_rich_idx) > MAX_PER_GROUP:
            gc_rich_idx = rng.choice(gc_rich_idx, MAX_PER_GROUP, replace=False)
        if len(gc_poor_idx) > MAX_PER_GROUP:
            gc_poor_idx = rng.choice(gc_poor_idx, MAX_PER_GROUP, replace=False)

        # 2. 白名单 + fasta
        wl_ids = get_whitelist_for(cell)
        motif_meme = workdir / f"motifs_{cell}.meme"
        extract_whitelist_meme(jaspar_meme, wl_ids, motif_meme)
        wl_ids = [m for m in wl_ids if m in motifs]

        rich_seqs = onehot_to_seq(x[gc_rich_idx])
        poor_seqs = onehot_to_seq(x[gc_poor_idx])
        rich_fasta = workdir / f"gcrich_{cell}.fasta"
        poor_fasta = workdir / f"gcpoor_{cell}.fasta"
        sequences_to_fasta(rich_seqs, rich_fasta, prefix="GR")
        sequences_to_fasta(poor_seqs, poor_fasta, prefix="GP")
        logger.info("%s: GC-rich %d 条 / GC-poor %d 条", cell, len(rich_seqs), len(poor_seqs))

        # 3. 富集分析（复用 analyze_motifs，high=GC-rich）
        rows = analyze_motifs(
            wl_ids, motif_names, motif_meme,
            rich_fasta, poor_fasta,
            fimo_bin, workdir / f"fimo_gc_{cell}",
            n_high=len(rich_seqs), n_low=len(poor_seqs),
            pvalue=pvalue,
        )

        # 4. 合并 motif GC 特征
        for r in rows:
            mgc, mcls = gc_feat.get(r["motif_id"], (float("nan"), "neutral"))
            r["motif_gc"] = mgc
            r["gc_class"] = mcls

        # 5. 保存 CSV
        csv_path = outdir / f"gc_motif_enrichment_{cell}.csv"
        with open(csv_path, "w", newline="", encoding="utf-8") as f:
            writer = csv.DictWriter(f, fieldnames=CSV_COLUMNS)
            writer.writeheader()
            for r in rows:
                writer.writerow({k: r.get(k, "") for k in CSV_COLUMNS})
        logger.info("已保存: %s (%d 行)", csv_path.name, len(rows))

        # 6. 生成 claim（仅显著富集/耗尽 + 有 GC 分类的）
        n_total = len(gc_frac)
        claims = []
        for r in rows:
            if r["direction"] == "not_significant":
                continue
            claims.append(build_claim(r, cell, n_total))
        claim_path = outdir / f"gc_motif_claims_{cell}.jsonl"
        with open(claim_path, "w", encoding="utf-8") as f:
            for c in claims:
                f.write(json.dumps(c, ensure_ascii=False) + "\n")
        logger.info("已生成 claim: %s (%d 条)", claim_path.name, len(claims))

        # 7. 摘要：GC 分类 × 富集方向交叉表
        logger.info("--- %s GC class × enrichment 交叉表 ---", cell)
        cross = {}
        for r in rows:
            key = (r["gc_class"], r["direction"])
            cross[key] = cross.get(key, 0) + 1
        for (cls, direc), n in sorted(cross.items()):
            logger.info("  %-10s %-18s %d", cls, direc, n)

    logger.info("=== 全部完成 ===")


if __name__ == "__main__":
    main()
