"""MPRA 数据加载与高/低活性组划分工具。

数据来源: Gosai et al. MPRA（640,029 条 200nt 序列，3 细胞系活性）。
y_train 3 列顺序: K562 (erythroid precursors), HepG2 (hepatocytes), SK-N-SH (neuroblastoma)。

模块职责:
    - load_mpra_data     : 从 h5 加载 one-hot 序列与三细胞系活性
    - onehot_to_seq      : one-hot -> ACGT 字符串（供 FIMO fasta 使用）
    - split_high_low     : 按分位数划分高/低活性组索引
    - validate_split     : 验证高低组活性差异显著
"""

from __future__ import annotations

import logging
from pathlib import Path
from typing import Optional, Tuple

import h5py
import numpy as np

logger = logging.getLogger("mpra_loader")

# one-hot 编码的碱基顺序（与参考代码 x.permute(0,2,1).argmax(dim=1) 一致）
BASES = "ACGT"
BASE_TO_IDX = {b: i for i, b in enumerate(BASES)}


def load_mpra_data(
    h5_path: str | Path,
    x_key: str = "x_train",
    y_key: str = "y_train",
) -> Tuple[np.ndarray, np.ndarray]:
    """从 h5 文件加载 one-hot 序列与活性矩阵。

    Returns:
        x_onehot: (N, 200, 4) int8，A/C/G/T one-hot
        y: (N, 3) float32，三细胞系活性 [K562, HepG2, SK-N-SH]
    """
    h5_path = Path(h5_path)
    if not h5_path.exists():
        raise FileNotFoundError(f"h5 文件不存在: {h5_path}")

    with h5py.File(h5_path, "r") as f:
        if x_key not in f or y_key not in f:
            raise KeyError(f"h5 缺少数据集 {x_key}/{y_key}，现有: {list(f.keys())}")
        x = f[x_key][...]
        y = f[y_key][...]

    x = np.asarray(x)
    y = np.asarray(y)
    logger.info("加载完成: x=%s (%s), y=%s", x.shape, x.dtype, y.shape)
    if x.ndim != 3 or x.shape[2] != 4:
        raise ValueError(f"x 应为 (N,200,4)，实际 {x.shape}")
    if y.ndim != 2 or y.shape[1] != 3:
        raise ValueError(f"y 应为 (N,3)，实际 {y.shape}")
    if x.shape[0] != y.shape[0]:
        raise ValueError(f"x/y 样本数不一致: {x.shape[0]} vs {y.shape[0]}")
    return x, y


def onehot_to_seq(x: np.ndarray) -> list[str]:
    """将 (N, 200, 4) one-hot 转为 ACGT 字符串列表。"""
    if x.ndim != 3 or x.shape[2] != 4:
        raise ValueError(f"x 应为 (N,200,4)，实际 {x.shape}")
    idx = np.argmax(x, axis=2).astype(np.int8)  # (N, 200)
    # 防御: one-hot 全零列（argmax 会取 0=A），罕见但无害
    seqs = ["".join(BASES[i] for i in row) for row in idx]
    return seqs


def split_high_low(
    y_col: np.ndarray,
    high_quantile: float = 0.90,
    low_quantile: float = 0.10,
    max_per_group: int = 20000,
    seed: int = 42,
) -> Tuple[np.ndarray, np.ndarray]:
    """按活性分位数划分高/低活性组索引。

    - 高活性组: y >= p_high 的样本，随机采样至多 max_per_group 条
    - 低活性组: y <= p_low 的样本，随机采样至多 max_per_group 条

    Returns:
        high_idx, low_idx: int 数组（原始样本索引）
    """
    y_col = np.asarray(y_col, dtype=np.float64).ravel()
    p_high = np.quantile(y_col, high_quantile)
    p_low = np.quantile(y_col, low_quantile)

    rng = np.random.default_rng(seed)
    high_all = np.where(y_col >= p_high)[0]
    low_all = np.where(y_col <= p_low)[0]

    if len(high_all) > max_per_group:
        high_idx = rng.choice(high_all, size=max_per_group, replace=False)
    else:
        high_idx = high_all
    if len(low_all) > max_per_group:
        low_idx = rng.choice(low_all, size=max_per_group, replace=False)
    else:
        low_idx = low_all

    high_idx = np.sort(np.asarray(high_idx, dtype=np.int64))
    low_idx = np.sort(np.asarray(low_idx, dtype=np.int64))

    logger.info(
        "划分: p90=%.3f p10=%.3f | 高活性候选 %d -> 取 %d | 低活性候选 %d -> 取 %d",
        p_high, p_low, len(high_all), len(high_idx), len(low_all), len(low_idx),
    )
    return high_idx, low_idx


def validate_split(
    y_col: np.ndarray,
    high_idx: np.ndarray,
    low_idx: np.ndarray,
) -> dict:
    """验证高低活性组活性差异（均值/中位数 + 独立样本 t 检验）。"""
    from scipy import stats

    high_vals = y_col[high_idx]
    low_vals = y_col[low_idx]
    t_stat, p_val = stats.ttest_ind(high_vals, low_vals, equal_var=False)

    result = {
        "n_high": int(len(high_vals)),
        "n_low": int(len(low_vals)),
        "high_mean": float(np.mean(high_vals)),
        "low_mean": float(np.mean(low_vals)),
        "high_median": float(np.median(high_vals)),
        "low_median": float(np.median(low_vals)),
        "t_stat": float(t_stat),
        "p_value": float(p_val),
    }
    logger.info(
        "验证: 高组 mean=%.3f median=%.3f | 低组 mean=%.3f median=%.3f | t=%.1f p=%.2e",
        result["high_mean"], result["high_median"],
        result["low_mean"], result["low_median"],
        result["t_stat"], result["p_value"],
    )
    return result
