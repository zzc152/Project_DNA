"""MPRA 知识三元组构建（步骤 4）。

从步骤 2（基序富集）与步骤 3（GC 分析）结果构建知识三元组，
schema: {
    head, relation, tail,
    source, statistic, cell_line, description
}

关系类型:
    - enriched_in_high   : TF 基序在高活性组显著富集 (head=TF, tail=cell_line)
    - depleted_in_high   : TF 基序在高活性组显著耗尽 (head=TF, tail=cell_line)
    - optimal_range_for  : GC 区间为活性最优区间 (head=GC区间, tail=cell_line)
    - suppressive_range_for : GC 区间为活性抑制区间 (head=GC区间, tail=cell_line)
"""

from __future__ import annotations

import csv
import json
import logging
from pathlib import Path
from typing import Iterable, Optional

logger = logging.getLogger("triplet_builder")

# 步骤2 阈值（与 motif_enrichment.py 保持一致）
ENRICH_FDR = 0.01
ENRICH_OR_MIN = 1.5
DEPLETE_OR_MAX = 1 / 1.5

SOURCE_MPRA = "Gosai et al. 2022 MPRA (640,029 x 200bp, K562/HepG2/SK-N-SH)"
SOURCE_FIMO = "FIMO scan of JASPAR2024 vertebrate motifs, p<1e-4"
SOURCE_GC = "GC content analysis of MPRA 200bp sequences"


def motif_triplets_from_csv(
    csv_path: str | Path,
    cell_line: str,
    enrich_fdr: float = ENRICH_FDR,
    enrich_or_min: float = ENRICH_OR_MIN,
    deplete_or_max: float = DEPLETE_OR_MAX,
) -> list[dict]:
    """从步骤2的 motif_enrichment_<CELL>.csv 生成三元组。

    Returns:
        [{head, relation, tail, source, statistic, cell_line, description}]
    """
    csv_path = Path(csv_path)
    if not csv_path.exists():
        raise FileNotFoundError(f"缺少步骤2输出: {csv_path}")

    triplets: list[dict] = []
    with open(csv_path, "r", encoding="utf-8") as f:
        reader = csv.DictReader(f)
        for row in reader:
            motif_id = row.get("motif_id", "")
            motif_name = row.get("motif_name", "")
            try:
                or_val = float(row.get("odds_ratio", "nan"))
                fdr = float(row.get("fdr", "1"))
                frac_high = float(row.get("frac_high", "0"))
                frac_low = float(row.get("frac_low", "0"))
            except ValueError:
                continue
            direction = row.get("direction", "not_significant")

            if direction == "enriched_in_high":
                rel = "enriched_in_high"
            elif direction == "depleted_in_high":
                rel = "depleted_in_high"
            else:
                continue

            head = motif_name or motif_id
            statistic = {
                "motif_id": motif_id,
                "odds_ratio": round(or_val, 4),
                "fdr": fdr,
                "frac_high": round(frac_high, 5),
                "frac_low": round(frac_low, 5),
            }
            if rel == "enriched_in_high":
                desc = (f"{motif_name} motif ({motif_id}) 在高活性组显著富集: "
                        f"OR={or_val:.2f}, FDR={fdr:.2e} "
                        f"(高活性组{frac_high*100:.1f}% vs 低活性组{frac_low*100:.1f}%)")
            else:
                desc = (f"{motif_name} motif ({motif_id}) 在高活性组显著耗尽: "
                        f"OR={or_val:.2f}, FDR={fdr:.2e} "
                        f"(高活性组{frac_high*100:.1f}% vs 低活性组{frac_low*100:.1f}%)")

            triplets.append({
                "head": head,
                "relation": rel,
                "tail": cell_line,
                "source": SOURCE_FIMO,
                "statistic": statistic,
                "cell_line": cell_line,
                "description": desc,
            })
    logger.info("%s: 生成 %d 条基序三元组", cell_line, len(triplets))
    return triplets


def gc_triplets_from_summary(
    summary_path: str | Path,
    cell_line: str,
) -> list[dict]:
    """从步骤3的 gc_summary.json 生成 GC 区间三元组。

    Returns:
        [{head, relation, tail, source, statistic, cell_line, description}]
    """
    summary_path = Path(summary_path)
    if not summary_path.exists():
        raise FileNotFoundError(f"缺少步骤3输出: {summary_path}")

    with open(summary_path, "r", encoding="utf-8") as f:
        summary = json.load(f)

    triplets: list[dict] = []
    results = summary.get("results", {})
    cell_result = results.get(cell_line, {})
    if not cell_result:
        logger.warning("%s 在 gc_summary.json 中无结果", cell_line)
        return triplets

    optimal = cell_result.get("optimal")
    suppressive = cell_result.get("suppressive")
    diff_ratio = cell_result.get("diff_ratio", 0.0)

    if optimal:
        bin_label = optimal.get("gc_bin", "")
        statistic = {
            "gc_bin": bin_label,
            "gc_lo": optimal.get("gc_lo"),
            "gc_hi": optimal.get("gc_hi"),
            "mean_activity": optimal.get("mean"),
            "n_sequences": optimal.get("n"),
            "diff_ratio_std": diff_ratio,
        }
        triplets.append({
            "head": f"GC {bin_label}",
            "relation": "optimal_range_for",
            "tail": cell_line,
            "source": SOURCE_GC,
            "statistic": statistic,
            "cell_line": cell_line,
            "description": (f"GC 含量 {bin_label} 区间在 {cell_line} 中平均活性最高 "
                            f"({optimal.get('mean'):.3f}, n={optimal.get('n')})，"
                            f"相对最低桶差 {diff_ratio:.2f} 倍活性标准差"),
        })

    if suppressive:
        bin_label = suppressive.get("gc_bin", "")
        statistic = {
            "gc_bin": bin_label,
            "gc_lo": suppressive.get("gc_lo"),
            "gc_hi": suppressive.get("gc_hi"),
            "mean_activity": suppressive.get("mean"),
            "n_sequences": suppressive.get("n"),
            "diff_ratio_std": diff_ratio,
        }
        triplets.append({
            "head": f"GC {bin_label}",
            "relation": "suppressive_range_for",
            "tail": cell_line,
            "source": SOURCE_GC,
            "statistic": statistic,
            "cell_line": cell_line,
            "description": (f"GC 含量 {bin_label} 区间在 {cell_line} 中平均活性最低 "
                            f"({suppressive.get('mean'):.3f}, n={suppressive.get('n')})，"
                            f"相对最高桶差 {diff_ratio:.2f} 倍活性标准差"),
        })

    logger.info("%s: 生成 %d 条 GC 三元组", cell_line, len(triplets))
    return triplets
