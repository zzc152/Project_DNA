"""DNA 形状预测与高/低活性组比较（Module 1 · 统计事实）。

基于标准 DNAshape 五核苷酸查询表（Zhou et al. 2013, NAR）预测 4 种核心形状:
    - MGW  : minor groove width（小沟宽度，每个碱基对一个值）
    - ProT : propeller twist（螺旋桨扭转，每个碱基对一个值）
    - Roll : roll（每个 step 一个值，5-mer 覆盖取平均）
    - HelT : helix twist（每个 step 一个值，5-mer 覆盖取平均）

查询语义（与 DNAshape/DNAshapeR 一致）:
    - 位置 j (2<=j<=N-3) 的 MGW/ProT: 5-mer seq[j-2:j+3] 的中央值
    - step j (1<=j<=N-3) 的 Roll/HelT: (pentamer[j-1].{roll1,twist1}
      + pentamer[j-2].{roll2,twist2}) / 2
    - 查询表中缺失的 5-mer 用其反向互补条目（值直接采用）

模块职责:
    - load_shape_table     : 加载查询表 tsv -> {5mer: {shape: value}}
    - build_lookup_arrays  : 构建 1024 全覆盖向量化查询数组
    - predict_shapes       : 一批 one-hot -> MGW/ProT (N,200) + Roll/HelT (N,199)
    - sequence_shape_features : 序列级形状统计特征
    - compare_high_low     : 高/低活性组逐特征 Welch t 检验 + Cohen's d
"""

from __future__ import annotations

import logging
from pathlib import Path

import numpy as np

logger = logging.getLogger("dna_shape")

BASES = "ACGT"
BASE_TO_IDX = {b: i for i, b in enumerate(BASES)}
_COMP = str.maketrans("ACGT", "TGCA")

# 查询表每行 91 字段中关键形状值的位置（见 extract_dna_shape_table.py）
SHAPE_COLS = {
    "MGW": 0,
    "ProT": 6,
    "roll1": 15,
    "roll2": 18,
    "twist1": 21,
    "twist2": 24,
}

# 输出用形状名 -> (组特征名前缀)
PER_BP_SHAPES = ("MGW", "ProT")
PER_STEP_SHAPES = ("Roll", "HelT")


def load_shape_table(tsv_path: str | Path) -> dict[str, dict[str, float]]:
    """加载查询表 tsv -> {5mer: {"MGW":.., "ProT":.., "roll1":.., "roll2":.., "twist1":.., "twist2":..}}。

    只保留标准 ACGT 条目与缺失 5-mer 的反向互补解析所需信息（简并 M/Q 条目忽略）。
    """
    tsv_path = Path(tsv_path)
    if not tsv_path.exists():
        raise FileNotFoundError(f"查询表不存在: {tsv_path}")
    table: dict[str, dict[str, float]] = {}
    n_std = 0
    with open(tsv_path, encoding="utf-8") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            parts = line.split()
            if len(parts) != 91:
                continue
            mer = parts[0]
            if not set(mer) <= set("ACGT"):
                continue  # 忽略 M/Q 简并条目
            vals = np.asarray(parts[1:], dtype=np.float64)
            table[mer] = {name: float(vals[col]) for name, col in SHAPE_COLS.items()}
            n_std += 1
    # 用反向互补补齐缺失的标准 5-mer（DNAshapeR 语义：反向互补条目 roll1<->roll2、twist1<->twist2 互换）
    n_rc = 0
    for p in _all_pentamers():
        if p not in table:
            rcp = p.translate(_COMP)[::-1]
            if rcp in table:
                table[p] = dict(table[rcp])
                # 链方向翻转：per-step 参数的 1/2 互换（MGW/ProT 为每 bp 值，不交换）
                table[p]["roll1"], table[p]["roll2"] = table[rcp]["roll2"], table[rcp]["roll1"]
                table[p]["twist1"], table[p]["twist2"] = table[rcp]["twist2"], table[rcp]["twist1"]
                n_rc += 1
    logger.info("查询表加载: %d 标准条目 + %d 反向互补补齐 = %d/1024",
                n_std, n_rc, len(table))
    if len(table) != 1024:
        raise ValueError(f"查询表覆盖不完整: {len(table)}/1024")
    return table


def _all_pentamers() -> list[str]:
    from itertools import product

    return ["".join(p) for p in product("ACGT", repeat=5)]


def build_lookup_arrays(table: dict[str, dict[str, float]]) -> tuple[np.ndarray, dict[str, int]]:
    """构建 (1024, n_shapes) 查询数组 + 形状列索引。

    5-mer 编码: base-4 (A=0,C=1,G=2,T=3)，p0*256+p1*64+p2*16+p3*4+p4。
    """
    shape_names = ["MGW", "ProT", "roll1", "roll2", "twist1", "twist2"]
    lookup = np.zeros((1024, len(shape_names)), dtype=np.float64)
    col_idx = {n: i for i, n in enumerate(shape_names)}
    for mer, vals in table.items():
        e = encode_pentamer(mer)
        for n in shape_names:
            lookup[e, col_idx[n]] = vals[n]
    return lookup, col_idx


def encode_pentamer(mer: str) -> int:
    e = 0
    for ch in mer:
        e = e * 4 + BASE_TO_IDX[ch]
    return e


def onehot_to_idx(x: np.ndarray) -> np.ndarray:
    """(N,200,4) one-hot -> (N,200) 碱基索引 (A=0..T=3)。"""
    if x.ndim != 3 or x.shape[2] != 4:
        raise ValueError(f"x 应为 (N,200,4)，实际 {x.shape}")
    return np.argmax(x, axis=2).astype(np.int16)


def _enc_window(idx: np.ndarray, start: int) -> np.ndarray:
    """(N,200) idx -> (N,) 窗口 [start, start+5) 的 5-mer 编码。"""
    w = idx[:, start : start + 5].astype(np.int32)
    return w[:, 0] * 256 + w[:, 1] * 64 + w[:, 2] * 16 + w[:, 3] * 4 + w[:, 4]


def predict_shapes(
    x: np.ndarray, lookup: np.ndarray, col_idx: dict[str, int]
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """预测一批序列的 4 种形状。

    Args:
        x: (N,200,4) one-hot
        lookup: (1024, n_shapes)
        col_idx: 形状列索引

    Returns:
        (mgw, prot, roll, helt):
            mgw/prot: (N, 200) 每碱基对值（边界位置为 NaN）
            roll/helt: (N, 199) 每 step 值（step j 位于碱基 j 与 j+1 之间，边界为 NaN）
    """
    idx = onehot_to_idx(x)
    n = x.shape[0]
    L = x.shape[1]
    mgw = np.full((n, L), np.nan, dtype=np.float32)
    prot = np.full((n, L), np.nan, dtype=np.float32)
    for j in range(2, L - 2):
        e = _enc_window(idx, j - 2)
        mgw[:, j] = lookup[e, col_idx["MGW"]]
        prot[:, j] = lookup[e, col_idx["ProT"]]

    # step 数 = L-1（step j 位于 bp j 与 j+1 之间）；DNAshapeR 语义：
    #   - j=1       : 只用 pentamer seq[0:5] 的 roll1/twist1
    #   - j=2..L-4  : (pentamer[j-1].roll1 + pentamer[j-2].roll2) / 2
    #   - j=L-3     : 只用 pentamer seq[L-5:L] 的 roll2/twist2（f1 越界为 null）
    roll = np.full((n, L - 1), np.nan, dtype=np.float32)
    helt = np.full((n, L - 1), np.nan, dtype=np.float32)
    for j in range(1, L - 2):  # step j=1..L-3 可预测
        if j == 1:
            e1 = _enc_window(idx, 0)  # pentamer seq[0:5]
            roll[:, 1] = lookup[e1, col_idx["roll1"]]
            helt[:, 1] = lookup[e1, col_idx["twist1"]]
        elif j == L - 3:
            e2 = _enc_window(idx, L - 5)  # pentamer seq[L-5:L]
            roll[:, L - 3] = lookup[e2, col_idx["roll2"]]
            helt[:, L - 3] = lookup[e2, col_idx["twist2"]]
        else:
            e1 = _enc_window(idx, j - 1)  # pentamer seq[j-1:j+4]
            e2 = _enc_window(idx, j - 2)  # pentamer seq[j-2:j+3]
            roll[:, j] = 0.5 * (lookup[e1, col_idx["roll1"]] + lookup[e2, col_idx["roll2"]])
            helt[:, j] = 0.5 * (lookup[e1, col_idx["twist1"]] + lookup[e2, col_idx["twist2"]])
    return mgw, prot, roll, helt


# ---------------------------------------------------------------------------
# 序列级特征
# ---------------------------------------------------------------------------

def sequence_shape_features(mgw, prot, roll, helt) -> dict[str, np.ndarray]:
    """将 (N,200)/(N,199) 形状矩阵折叠为序列级特征。

    Returns: {feature_name: (N,) 数组}
        MGW_mean/ProT_mean/Roll_mean/HelT_mean   : 全序列均值
        MGW_std/ProT_std/Roll_std/HelT_std       : 全序列波动
        MGW_min/MGW_max/ProT_min/ProT_max        : 极值（MGW/ProT 每 bp）
    """
    feats: dict[str, np.ndarray] = {}
    for name, arr in (("MGW", mgw), ("ProT", prot)):
        valid = ~np.isnan(arr)
        with np.errstate(invalid="ignore"):
            feats[f"{name}_mean"] = np.where(valid.any(axis=1), np.nanmean(arr, axis=1), np.nan)
            feats[f"{name}_std"] = np.where(valid.any(axis=1), np.nanstd(arr, axis=1), np.nan)
            feats[f"{name}_min"] = np.where(valid.any(axis=1), np.nanmin(arr, axis=1), np.nan)
            feats[f"{name}_max"] = np.where(valid.any(axis=1), np.nanmax(arr, axis=1), np.nan)
    for name, arr in (("Roll", roll), ("HelT", helt)):
        valid = ~np.isnan(arr)
        with np.errstate(invalid="ignore"):
            feats[f"{name}_mean"] = np.where(valid.any(axis=1), np.nanmean(arr, axis=1), np.nan)
            feats[f"{name}_std"] = np.where(valid.any(axis=1), np.nanstd(arr, axis=1), np.nan)
    return feats


# ---------------------------------------------------------------------------
# 组间比较
# ---------------------------------------------------------------------------

def compare_high_low(
    high_feats: dict[str, np.ndarray],
    low_feats: dict[str, np.ndarray],
    cell_line: str,
) -> list[dict]:
    """逐特征比较高/低活性组。

    Returns: [{shape, feature, direction, cohens_d, p_value, mean_high,
               mean_low, sd_high, sd_low, n_high, n_low}]，按 |cohens_d| 降序。
    """
    from scipy import stats

    rows: list[dict] = []
    for feat, h in high_feats.items():
        l = low_feats[feat]
        mh, ml = float(np.nanmean(h)), float(np.nanmean(l))
        sh, sl = float(np.nanstd(h, ddof=1)), float(np.nanstd(l, ddof=1))
        nh, nl = int((~np.isnan(h)).sum()), int((~np.isnan(l)).sum())
        if nh < 10 or nl < 10 or sh == 0 and sl == 0:
            continue
        # Welch t 检验
        t_stat, p_val = stats.ttest_ind(h[~np.isnan(h)], l[~np.isnan(l)], equal_var=False)
        # Cohen's d（池化标准差）
        pooled_sd = np.sqrt(((nh - 1) * sh ** 2 + (nl - 1) * sl ** 2) / (nh + nl - 2))
        d = (mh - ml) / pooled_sd if pooled_sd > 0 else 0.0
        direction = "high" if d > 0 else "low"
        rows.append({
            "shape": feat.split("_")[0],
            "feature": feat,
            "direction": direction,
            "cohens_d": round(float(d), 4),
            "p_value": float(p_val),
            "mean_high": round(mh, 4),
            "mean_low": round(ml, 4),
            "sd_high": round(sh, 4),
            "sd_low": round(sl, 4),
            "n_high": nh,
            "n_low": nl,
        })
    rows.sort(key=lambda r: -abs(r["cohens_d"]))
    logger.info("%s: %d 个形状特征完成比较", cell_line, len(rows))
    return rows
