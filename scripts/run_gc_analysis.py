"""MPRA 知识库 · 步骤 3：GC 含量与活性关系分析。

对每个细胞系:
    - 计算全部序列的 GC 比例，按 5% GC 分桶
    - 统计每桶平均活性（mean/std/median/n）
    - 判定最优 GC 区间（平均活性最高）与抑制 GC 区间（最低）
    - 保存 gc_analysis_<CELL>.csv + gc_summary.json

用法（项目根目录）:
    python scripts/run_gc_analysis.py
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

from src.knowledge.gc_analysis import bin_gc_activity, find_optimal_suppressive, gc_content_from_onehot  # noqa: E402

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s | %(levelname)s | %(message)s",
    datefmt="%H:%M:%S",
)
logger = logging.getLogger("run_gc_analysis")

CSV_COLUMNS = ["gc_bin", "gc_lo", "gc_hi", "n", "mean", "std", "median", "reliable"]


def main() -> None:
    parser = argparse.ArgumentParser(description="MPRA GC 含量分析")
    parser.add_argument("--config", default="configs/default_config.yaml")
    parser.add_argument("--outdir", default="data/processed")
    parser.add_argument("--bin-width", type=float, default=0.05)
    parser.add_argument("--min-count", type=int, default=100)
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

    outdir = ROOT / args.outdir
    outdir.mkdir(parents=True, exist_ok=True)

    logger.info("加载 MPRA 数据: %s", h5_path)
    with h5py.File(h5_path, "r") as f:
        x = f[x_key][...]
        y = f[y_key][...]
    x = np.asarray(x)
    y = np.asarray(y)
    logger.info("x=%s y=%s", x.shape, y.shape)

    # GC 比例只算一次（与细胞系无关）
    gc_frac = gc_content_from_onehot(x)
    logger.info("GC 比例: mean=%.4f std=%.4f", gc_frac.mean(), gc_frac.std())

    summary_all = {"cell_lines": cell_lines, "bin_width": args.bin_width, "min_count": args.min_count, "results": {}}
    for ci, cell in enumerate(cell_lines):
        logger.info("############ %s ############", cell)
        y_col = y[:, ci]
        rows = bin_gc_activity(gc_frac, y_col, bin_width=args.bin_width, min_count=args.min_count)
        verdict = find_optimal_suppressive(rows, min_count=args.min_count)

        # 保存 CSV
        csv_path = outdir / f"gc_analysis_{cell}.csv"
        with open(csv_path, "w", newline="", encoding="utf-8") as f:
            writer = csv.DictWriter(f, fieldnames=CSV_COLUMNS)
            writer.writeheader()
            for r in rows:
                writer.writerow({k: r.get(k, "") for k in CSV_COLUMNS})
        logger.info("已保存: %s (%d 桶)", csv_path, len(rows))

        summary_all["results"][cell] = {
            "gc_mean": round(float(gc_frac.mean()), 4),
            "gc_std": round(float(gc_frac.std()), 4),
            "optimal": verdict["optimal"],
            "suppressive": verdict["suppressive"],
            "diff_ratio": verdict["diff_ratio"],
        }

        opt = verdict["optimal"]
        sup = verdict["suppressive"]
        if opt:
            logger.info("最优 GC 区间: %s (mean=%.3f, n=%d)", opt["gc_bin"], opt["mean"], opt["n"])
        else:
            logger.info("最优 GC 区间: (无显著差异)")
        if sup:
            logger.info("抑制 GC 区间: %s (mean=%.3f, n=%d)", sup["gc_bin"], sup["mean"], sup["n"])
        else:
            logger.info("抑制 GC 区间: (无显著差异)")

        # 打印 top/bottom 桶明细
        ranked = sorted([r for r in rows if r["reliable"]], key=lambda r: -r["mean"])
        logger.info("Top3 高活性 GC 桶:")
        for r in ranked[:3]:
            logger.info("  %-10s mean=%+.3f n=%d", r["gc_bin"], r["mean"], r["n"])
        logger.info("Bottom3 低活性 GC 桶:")
        for r in ranked[-3:]:
            logger.info("  %-10s mean=%+.3f n=%d", r["gc_bin"], r["mean"], r["n"])

    summary_path = outdir / "gc_summary.json"
    with open(summary_path, "w", encoding="utf-8") as f:
        json.dump(summary_all, f, ensure_ascii=False, indent=2)
    logger.info("已保存汇总: %s", summary_path)
    logger.info("=== 全部完成 ===")


if __name__ == "__main__":
    main()
