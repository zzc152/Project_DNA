"""GC 含量与活性关系分析（步骤 3）。

对每个细胞系:
    - 计算每条 200nt 序列的 GC 比例（one-hot 中 C/G 列占比）
    - 按 5% GC 分桶（0-5%, 5-10%, ..., 95-100%）
    - 统计每桶平均活性/标准差/样本数
    - 判定最优 GC 区间（平均活性最高）与抑制 GC 区间（平均活性最低）

输出: gc_analysis_<CELL>.csv + gc_summary.json
"""

from __future__ import annotations

import logging
import math
from pathlib import Path
from typing import Iterable, Optional

import numpy as np

logger = logging.getLogger("gc_analysis")

# BASES = "ACGT": A=0, C=1, G=2, T=3。GC = 列 1(C) + 列 2(G)
GC_COLS = (1, 2)


def gc_content_from_onehot(x: np.ndarray) -> np.ndarray:
    """(N, 200, 4) one-hot -> (N,) GC 比例数组（0~1）。"""
    x = np.asarray(x)
    if x.ndim != 3 or x.shape[2] != 4:
        raise ValueError(f"x 应为 (N,200,4)，实际 {x.shape}")
    gc_count = x[:, :, GC_COLS[0]].astype(np.float32) + x[:, :, GC_COLS[1]].astype(np.float32)
    gc_frac = gc_count.sum(axis=1) / x.shape[1]  # (N,)
    return gc_frac


def gc_content_from_seq(seqs: Iterable[str]) -> np.ndarray:
    """ACGT 字符串列表 -> (N,) GC 比例数组（备用，供纯字符串输入）。"""
    fracs = []
    for s in seqs:
        s = s.upper()
        if len(s) == 0:
            fracs.append(0.0)
            continue
        fracs.append((s.count("G") + s.count("C")) / len(s))
    return np.asarray(fracs, dtype=np.float64)


def bin_gc_activity(
    gc_frac: np.ndarray,
    y_col: np.ndarray,
    bin_width: float = 0.05,
    min_count: int = 100,
) -> list[dict]:
    """按 GC 分桶统计平均活性。

    Args:
        gc_frac: (N,) GC 比例
        y_col: (N,) 活性值
        bin_width: 桶宽，默认 0.05（5%）
        min_count: 桶内样本数低于该值的桶标记为不可靠

    Returns:
        按桶排序的 [{gc_bin, gc_lo, gc_hi, n, mean, std, median}] 列表
    """
    gc_frac = np.asarray(gc_frac, dtype=np.float64)
    y_col = np.asarray(y_col, dtype=np.float64).ravel()
    if gc_frac.shape != y_col.shape:
        raise ValueError(f"gc_frac 与 y_col 长度不一致: {gc_frac.shape} vs {y_col.shape}")

    nbins = int(math.ceil(1.0 / bin_width))
    rows = []
    for b in range(nbins):
        lo = b * bin_width
        hi = (b + 1) * bin_width
        # 最后一桶含 1.0（如 [0.95, 1.0] 闭区间）
        if b == nbins - 1:
            mask = (gc_frac >= lo) & (gc_frac <= 1.0 + 1e-9)
        else:
            mask = (gc_frac >= lo) & (gc_frac < hi)
        vals = y_col[mask]
        n = int(mask.sum())
        if n == 0:
            rows.append({
                "gc_bin": f"{lo:.0%}-{hi:.0%}",
                "gc_lo": round(lo, 4), "gc_hi": round(hi, 4),
                "n": 0, "mean": None, "std": None, "median": None,
                "reliable": False,
            })
            continue
        rows.append({
            "gc_bin": f"{lo:.0%}-{hi:.0%}",
            "gc_lo": round(lo, 4), "gc_hi": round(hi, 4),
            "n": n,
            "mean": float(np.mean(vals)),
            "std": float(np.std(vals)),
            "median": float(np.median(vals)),
            "reliable": n >= min_count,
        })
    return rows


def find_optimal_suppressive(
    rows: list[dict],
    min_count: int = 100,
    min_diff_frac: float = 0.2,
) -> dict:
    """在可靠桶中判定最优/抑制 GC 区间。

    Args:
        rows: bin_gc_activity 输出
        min_count: 桶最小样本数（可靠性阈值）
        min_diff_frac: 最优与抑制桶均值差 / 全数据活性标准差，低于该值则视为无显著差异

    Returns:
        {"optimal": {...} | None, "suppressive": {...} | None, "diff_ratio": float}
    """
    reliable = [r for r in rows if r["reliable"] and r["n"] >= min_count]
    if len(reliable) < 2:
        return {"optimal": None, "suppressive": None, "diff_ratio": 0.0}

    # 全局活性标准差（用所有桶的合并方差近似）
    all_means = np.asarray([r["mean"] for r in reliable], dtype=np.float64)
    all_stds = np.asarray([r["std"] for r in reliable], dtype=np.float64)
    all_n = np.asarray([r["n"] for r in reliable], dtype=np.float64)
    total_n = all_n.sum()
    global_var = (all_stds ** 2 * all_n).sum() / total_n
    global_std = float(np.sqrt(global_var)) if global_var > 0 else 1.0

    best = max(reliable, key=lambda r: r["mean"])
    worst = min(reliable, key=lambda r: r["mean"])
    diff_ratio = (best["mean"] - worst["mean"]) / global_std if global_std > 0 else 0.0

    optimal = None
    if best["mean"] > worst["mean"] and diff_ratio >= min_diff_frac:
        optimal = {k: best[k] for k in ("gc_bin", "gc_lo", "gc_hi", "n", "mean", "median")}
        optimal["mean_minus_baseline"] = round(best["mean"] - worst["mean"], 4)
    suppressive = None
    if worst["mean"] < best["mean"] and diff_ratio >= min_diff_frac:
        suppressive = {k: worst[k] for k in ("gc_bin", "gc_lo", "gc_hi", "n", "mean", "median")}
        suppressive["mean_minus_baseline"] = round(worst["mean"] - best["mean"], 4)

    return {
        "optimal": optimal,
        "suppressive": suppressive,
        "diff_ratio": round(diff_ratio, 4),
        "global_std": round(global_std, 4),
    }


def summarize_gc(
    gc_frac: np.ndarray,
    y_col: np.ndarray,
    bin_width: float = 0.05,
    min_count: int = 100,
) -> dict:
    """一站式：分桶 + 最优/抑制区间判定，返回汇总 dict。"""
    rows = bin_gc_activity(gc_frac, y_col, bin_width=bin_width, min_count=min_count)
    result = find_optimal_suppressive(rows, min_count=min_count)
    result["bins"] = rows
    return result
