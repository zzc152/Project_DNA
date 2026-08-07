"""MPRA 知识库 · 步骤 1：加载数据并划分高/低活性组。

对每个细胞系 (K562, HepG2, SK-N-SH):
    - 高活性组: 活性 >= p90（至多 20000 条，随机采样）
    - 低活性组: 活性 <= p10（至多 20000 条，随机采样）
保存索引到 data/processed/high_<CELL>_indices.npy / low_<CELL>_indices.npy，
并验证高低组活性差异显著。

用法（项目根目录）:
    python scripts/analysis/split_mpra_groups.py
    python scripts/analysis/split_mpra_groups.py --config configs/default_config.yaml
"""

from __future__ import annotations

import argparse
import json
import logging
import sys
from pathlib import Path

import numpy as np
import yaml

ROOT = Path(__file__).resolve().parent.parent
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from src.knowledge.mpra_loader import load_mpra_data, split_high_low, validate_split  # noqa: E402

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s | %(levelname)s | %(message)s",
    datefmt="%H:%M:%S",
)
logger = logging.getLogger("split_mpra_groups")


def main() -> None:
    parser = argparse.ArgumentParser(description="划分 MPRA 高/低活性组")
    parser.add_argument("--config", default="configs/default_config.yaml", help="配置文件")
    parser.add_argument("--outdir", default="data/processed", help="输出目录（相对项目根目录）")
    args = parser.parse_args()

    cfg_path = ROOT / args.config
    if not cfg_path.exists():
        logger.error("配置文件不存在: %s", cfg_path)
        sys.exit(1)
    with open(cfg_path, "r", encoding="utf-8") as f:
        cfg = yaml.safe_load(f)
    mpra_cfg = cfg.get("mpra", {})
    if not mpra_cfg:
        logger.error("配置缺少 mpra 段")
        sys.exit(1)

    h5_path = mpra_cfg.get("h5_path")
    cell_lines = mpra_cfg.get("cell_lines", ["K562", "HepG2", "SK-N-SH"])
    x_key = mpra_cfg.get("x_key", "x_train")
    y_key = mpra_cfg.get("y_key", "y_train")
    high_q = mpra_cfg.get("high_quantile", 0.90)
    low_q = mpra_cfg.get("low_quantile", 0.10)
    max_per = mpra_cfg.get("max_per_group", 20000)
    seed = mpra_cfg.get("split_seed", 42)
    outdir = ROOT / args.outdir
    outdir.mkdir(parents=True, exist_ok=True)

    # 1. 加载数据
    logger.info("加载 MPRA 数据: %s", h5_path)
    x, y = load_mpra_data(h5_path, x_key, y_key)
    if y.shape[1] != len(cell_lines):
        logger.warning("y 列数 %d 与 cell_lines 数 %d 不一致，按列数截断", y.shape[1], len(cell_lines))
        cell_lines = cell_lines[: y.shape[1]]

    # 2. 逐细胞系划分
    summary = {}
    for i, cell in enumerate(cell_lines):
        logger.info("===== %s (列 %d) =====", cell, i)
        y_col = y[:, i]
        high_idx, low_idx = split_high_low(
            y_col,
            high_quantile=high_q,
            low_quantile=low_q,
            max_per_group=max_per,
            seed=seed,
        )
        stat = validate_split(y_col, high_idx, low_idx)

        high_path = outdir / f"high_{cell}_indices.npy"
        low_path = outdir / f"low_{cell}_indices.npy"
        np.save(high_path, high_idx)
        np.save(low_path, low_idx)
        logger.info("已保存: %s (%d), %s (%d)", high_path.name, len(high_idx), low_path.name, len(low_idx))

        summary[cell] = {
            "n_high": stat["n_high"],
            "n_low": stat["n_low"],
            "high_mean": stat["high_mean"],
            "low_mean": stat["low_mean"],
            "high_median": stat["high_median"],
            "low_median": stat["low_median"],
            "p_value": stat["p_value"],
            "high_file": str(high_path),
            "low_file": str(low_path),
        }

    # 3. 汇总
    summary_path = outdir / "split_summary.json"
    with open(summary_path, "w", encoding="utf-8") as f:
        json.dump(summary, f, ensure_ascii=False, indent=2)
    logger.info("汇总已写入: %s", summary_path)
    logger.info("=== 结果摘要 ===")
    for cell, s in summary.items():
        logger.info(
            "%s: 高 %d (mean=%.3f) | 低 %d (mean=%.3f) | 差异 p=%.2e",
            cell, s["n_high"], s["high_mean"], s["n_low"], s["low_mean"], s["p_value"],
        )


if __name__ == "__main__":
    main()
