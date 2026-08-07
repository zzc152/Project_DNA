#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""motif_gc_features.py — 计算白名单 motif 自身的 GC 含量特征（从 PWM）。

原理：
    JASPAR PWM 每列是 4 碱基概率（A,C,G,T），motif 序列水平的 GC 倾向 =
    每列 C+G 概率，按该列信息量（保守度）加权平均 → motif_gc ∈ [0,1]。

    GC 倾向性判定（通用阈值）:
        motif_gc >= 0.60  → GC-rich motif（如 SP1 GC-box、KLF 家族）
        motif_gc <= 0.40  → AT-rich motif（如 GATA1、REST）
        之间             → 中性 motif

输出:
    data/processed/motif_gc_features.csv
    data/processed/motif_gc_features.json  （含每 motif 判定结果）

用法（项目根目录）:
    python scripts/analysis/run_motif_gc_features.py
"""

from __future__ import annotations

import argparse
import csv
import json
import logging
import math
import sys
from pathlib import Path

import numpy as np
import yaml

ROOT = Path(__file__).resolve().parents[2]  # 项目根（scripts/analysis/ -> 上两级）
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from src.knowledge.motif_enrichment import parse_jaspar_meme  # noqa: E402

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s | %(levelname)s | %(message)s",
    datefmt="%H:%M:%S",
)
logger = logging.getLogger("motif_gc_features")

# GC-rich / AT-rich 判定阈值
GC_RICH_THRESHOLD = 0.60
AT_RICH_THRESHOLD = 0.40

CSV_COLUMNS = [
    "motif_id", "motif_name", "w", "motif_gc",
    "gc_class", "gc_class_detail",
]


def column_info(pwm_col: np.ndarray) -> float:
    """单列 PWM 信息量（bits）：2 + Σ p·log2(p)。"""
    p = np.asarray(pwm_col, dtype=np.float64)
    p = p[p > 0]
    if len(p) < 4:
        return 0.0
    return 2.0 + float(np.sum(p * np.log2(p)))


def pwm_to_gc(block: list[str]) -> tuple[float, int]:
    """从 MEME letter-probability matrix 块计算加权 GC 含量。

    Returns:
        (motif_gc, w): 加权 GC 比例 + motif 宽度
    """
    rows = []
    for line in block:
        s = line.strip()
        if not s or s.startswith("letter-probability") or s.startswith("URL") or s.startswith("URL"):
            continue
        try:
            vals = [float(v) for v in s.split()]
        except ValueError:
            continue
        if len(vals) == 4:
            rows.append(vals)
    if not rows:
        return float("nan"), 0

    pwm = np.asarray(rows)  # (w, 4) 行=位置，列序 A,C,G,T
    gc_col = pwm[:, 1] + pwm[:, 2]  # 每个位置的 C+G 概率，形状 (w,)
    # 每个位置的保守度（信息量）：对每行（4 碱基概率）计算
    infos = np.asarray([column_info(row) for row in pwm])  # 形状 (w,)
    if infos.sum() == 0:
        return float(np.mean(gc_col)), pwm.shape[0]
    # 信息量加权（保守位置权重大）
    motif_gc = float(np.sum(gc_col * infos) / np.sum(infos))
    return motif_gc, pwm.shape[0]


def classify_gc(motif_gc: float) -> str:
    if motif_gc >= GC_RICH_THRESHOLD:
        return "gc_rich"
    if motif_gc <= AT_RICH_THRESHOLD:
        return "at_rich"
    return "neutral"


def main() -> None:
    parser = argparse.ArgumentParser(description="计算白名单 motif 的 GC 含量特征")
    parser.add_argument("--config", default="configs/default_config.yaml")
    parser.add_argument("--outdir", default="data/processed")
    args = parser.parse_args()

    with open(ROOT / args.config, "r", encoding="utf-8") as f:
        cfg = yaml.safe_load(f)
    mpra = cfg.get("mpra", {})
    jaspar_meme = mpra["jaspar_meme"]

    outdir = ROOT / args.outdir
    outdir.mkdir(parents=True, exist_ok=True)

    motifs = parse_jaspar_meme(jaspar_meme)
    logger.info("JASPAR 全库 %d 个 motif", len(motifs))

    # 白名单交集：三个细胞系 + common
    from src.knowledge.motif_enrichment import MOTIF_WHITELIST

    wl_ids = []
    for group in MOTIF_WHITELIST.values():
        for mid in group:
            if mid not in wl_ids:
                wl_ids.append(mid)
    wl_ids = [m for m in wl_ids if m in motifs]
    logger.info("白名单可解析 %d 个 motif", len(wl_ids))

    rows = []
    for mid in wl_ids:
        info = motifs[mid]
        motif_gc, w = pwm_to_gc(info["block"])
        cls = classify_gc(motif_gc)
        rows.append({
            "motif_id": mid,
            "motif_name": info["name"],
            "w": w,
            "motif_gc": round(motif_gc, 4),
            "gc_class": cls,
            "gc_class_detail": (
                f"GC-rich motif (PWM GC={motif_gc:.2f})" if cls == "gc_rich"
                else f"AT-rich motif (PWM GC={motif_gc:.2f})" if cls == "at_rich"
                else f"neutral motif (PWM GC={motif_gc:.2f})"
            ),
        })

    rows.sort(key=lambda r: -r["motif_gc"])
    csv_path = outdir / "motif_gc_features.csv"
    with open(csv_path, "w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=CSV_COLUMNS)
        writer.writeheader()
        for r in rows:
            writer.writerow(r)

    json_path = outdir / "motif_gc_features.json"
    with open(json_path, "w", encoding="utf-8") as f:
        json.dump(rows, f, ensure_ascii=False, indent=2)

    logger.info("已保存: %s / %s", csv_path.name, json_path.name)
    logger.info("=== GC 特征排序（从高到低）===")
    for r in rows:
        logger.info("  %-10s %-14s GC=%.2f  %s", r["motif_id"], r["motif_name"], r["motif_gc"], r["gc_class"])


if __name__ == "__main__":
    main()
