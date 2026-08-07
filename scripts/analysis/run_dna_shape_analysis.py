"""Module 1 · DNA 形状分析（步骤 3b）：MGW/ProT/Roll/HelT 与增强子活性。

对每个细胞系:
    - 复用 split_mpra_groups 的高/低活性索引（或按 p90/p10 重算，seed=42）
    - 对高/低组各 20000 条 200nt 序列预测 4 种 DNAshape 形状
    - 折叠为序列级特征（均值/波动/极值）
    - 逐特征 Welch t 检验 + Cohen's d，判定方向（high=富集于高活性）
    - 保存 shape_analysis_<CELL>.csv + shape_summary.json

用法（项目根目录）:
    python scripts/analysis/run_dna_shape_analysis.py
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

ROOT = Path(__file__).resolve().parent.parent
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from src.knowledge.dna_shape import (  # noqa: E402
    build_lookup_arrays,
    compare_high_low,
    load_shape_table,
    predict_shapes,
    sequence_shape_features,
)
from src.knowledge.mpra_loader import split_high_low  # noqa: E402

logging.basicConfig(level=logging.INFO, format="%(asctime)s | %(levelname)s | %(message)s", datefmt="%H:%M:%S")
logger = logging.getLogger("run_dna_shape_analysis")

CSV_COLUMNS = [
    "shape", "feature", "direction", "cohens_d", "p_value",
    "mean_high", "mean_low", "sd_high", "sd_low", "n_high", "n_low",
]

SHAPE_LABELS = {
    "MGW": "小沟宽度",
    "ProT": "螺旋桨扭转",
    "Roll": "Roll 弯曲",
    "HelT": "螺旋扭转",
}


def load_split_indices(outdir: Path, cell: str, y_col: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """优先复用 split_mpra_groups 的索引文件，否则重算。"""
    hi_path = outdir / f"high_{cell}_indices.npy"
    lo_path = outdir / f"low_{cell}_indices.npy"
    if hi_path.exists() and lo_path.exists():
        hi = np.load(hi_path)
        lo = np.load(lo_path)
        logger.info("%s: 复用索引 high=%d low=%d", cell, len(hi), len(lo))
        return hi, lo
    hi, lo = split_high_low(y_col, high_quantile=0.90, low_quantile=0.10, max_per_group=20000, seed=42)
    logger.info("%s: 重算索引 high=%d low=%d", cell, len(hi), len(lo))
    return hi, lo


def main() -> None:
    parser = argparse.ArgumentParser(description="MPRA DNA 形状分析")
    parser.add_argument("--config", default="configs/default_config.yaml")
    parser.add_argument("--outdir", default="data/processed")
    parser.add_argument("--cells", default=None, help="逗号分隔的细胞系子集，如 K562")
    parser.add_argument("--shape-table", default="data/processed/dna_shape_querytable.tsv")
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

    outdir = ROOT / args.outdir
    outdir.mkdir(parents=True, exist_ok=True)

    # 查询表
    table_path = ROOT / args.shape_table
    if not table_path.exists():
        logger.error("查询表不存在: %s\n请先运行 scripts/extract/extract_dna_shape_table.py", table_path)
        sys.exit(1)
    table = load_shape_table(table_path)
    lookup, col_idx = build_lookup_arrays(table)
    logger.info("查询数组: %s", lookup.shape)

    logger.info("加载 MPRA 数据: %s", h5_path)
    with h5py.File(h5_path, "r") as f:
        x = f[x_key][...]
        y = f[y_key][...]
    x = np.asarray(x)
    y = np.asarray(y)
    logger.info("x=%s y=%s", x.shape, y.shape)

    summary_all = {"cell_lines": cell_lines, "results": {}}
    for ci, cell in enumerate(cell_lines):
        logger.info("############ %s ############", cell)
        y_col = y[:, ci]
        hi_idx, lo_idx = load_split_indices(outdir, cell, y_col)

        x_high = x[hi_idx]
        x_low = x[lo_idx]
        logger.info("预测形状: high=%d low=%d", x_high.shape[0], x_low.shape[0])
        mgw_h, prot_h, roll_h, helt_h = predict_shapes(x_high, lookup, col_idx)
        mgw_l, prot_l, roll_l, helt_l = predict_shapes(x_low, lookup, col_idx)
        logger.info("形状矩阵: MGW/ProT %s  Roll/HelT %s", mgw_h.shape, roll_h.shape)

        feats_high = sequence_shape_features(mgw_h, prot_h, roll_h, helt_h)
        feats_low = sequence_shape_features(mgw_l, prot_l, roll_l, helt_l)
        rows = compare_high_low(feats_high, feats_low, cell)

        # 保存 CSV
        csv_path = outdir / f"shape_analysis_{cell}.csv"
        with open(csv_path, "w", newline="", encoding="utf-8") as f:
            writer = csv.DictWriter(f, fieldnames=CSV_COLUMNS)
            writer.writeheader()
            for r in rows:
                writer.writerow({k: r.get(k, "") for k in CSV_COLUMNS})
        logger.info("已保存: %s (%d 特征)", csv_path, len(rows))

        sig = [r for r in rows if r["p_value"] < 0.01 and abs(r["cohens_d"]) >= 0.2]
        summary_all["results"][cell] = {
            "n_high": rows[0]["n_high"] if rows else 0,
            "n_low": rows[0]["n_low"] if rows else 0,
            "significant": [
                {"feature": r["feature"], "direction": r["direction"],
                 "cohens_d": r["cohens_d"], "p_value": r["p_value"]}
                for r in sig
            ],
        }
        for r in sig[:6]:
            logger.info("  %-10s %s  d=%+.3f  p=%.2e  (high=%.3f low=%.3f)",
                        r["feature"], SHAPE_LABELS.get(r["shape"], r["shape"]),
                        r["cohens_d"], r["p_value"], r["mean_high"], r["mean_low"])

    summary_path = outdir / "shape_summary.json"
    with open(summary_path, "w", encoding="utf-8") as f:
        json.dump(summary_all, f, ensure_ascii=False, indent=2)
    logger.info("汇总已保存: %s", summary_path)


if __name__ == "__main__":
    main()
