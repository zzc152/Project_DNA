"""[DEPRECATED] 升级知识库为"行动型"三元组（可指导模型生成）。

⚠️ 已废弃：本脚本产出旧版三元组格式。
项目现统一使用 claim 版本知识库（knowledge_base_stat.jsonl，8 键 schema），
由 build_stat_claims.py + build_literature_claims.py 构建。
旧文件已移至 data/backup/。请勿再使用本脚本。

输入:
    data/processed/knowledge_base_mpra.jsonl  (40 条 MPRA 统计三元组)
    data/processed/knowledge_base_full.jsonl  (3946 条文献三元组)
    data/processed/motif_enrichment_<CELL>.csv
    data/processed/gc_analysis_<CELL>.csv
    data/processed/split_summary.json

输出:
    data/processed/knowledge_base_actionable.jsonl

行动型三元组 schema:
    {
        "head": "motif:FLI1" | "gc_range:0.80-0.85" | "tf:TP53" | ...,
        "relation": "increases_activity" | "decreases_activity" |
                    "optimal_for" | "suppressive_for" |
                    "regulates" | "binds_motif" | "associated_with",
        "tail": "cell_line:K562" | "gene:TP53" | "disease:ALZHEIMER'S DISEASE" | ...,
        "effect_size": {...},          # 数值证据（MPRA 统计有，文献为空）
        "confidence": {...},           # p_value/fdr/sample_size 或 evidence_type
        "recommendation": "...",       # 可直接指导生成的行动建议
        "source": "...",               # 数据来源
        "knowledge_type": "motif_activity" | "gc_activity" | "literature"
    }

用法（项目根目录）:
    python scripts/upgrade_knowledge_base.py
"""

from __future__ import annotations

import argparse
import csv
import json
import logging
import math
import sys
from collections import Counter
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s | %(levelname)s | %(message)s",
    datefmt="%H:%M:%S",
)
logger = logging.getLogger("upgrade_knowledge_base")

DEFAULT_CELLS = ["K562", "HepG2", "SK-N-SH"]
SOURCE_MPRA = "Gosai et al. 2022 MPRA (640,029 x 200bp, K562/HepG2/SK-N-SH)"
SOURCE_FIMO = "FIMO scan of JASPAR2024 vertebrate motifs, p<1e-4"
SOURCE_GC = "GC content analysis of MPRA 200bp sequences"
SOURCE_LIT = "PubMed literature extraction (Qwen2.5-7B-Instruct)"


# ---------------------------------------------------------------- 类型 1: 基序-活性效应
def _motif_actionable_triplets(
    csv_path: Path, cell: str, split: dict
) -> list[dict]:
    """从 motif_enrichment CSV 生成行动型基序-活性三元组。"""
    out: list[dict] = []
    sample_size = (split.get("n_high", 0) or 0) + (split.get("n_low", 0) or 0)
    with open(csv_path, "r", encoding="utf-8") as f:
        reader = csv.DictReader(f)
        for row in reader:
            direction = row.get("direction", "")
            if direction not in ("enriched_in_high", "depleted_in_high"):
                continue
            name = row.get("motif_name", "").strip() or row.get("motif_id", "").strip()
            try:
                or_val = float(row.get("odds_ratio", "nan"))
                p_val = float(row.get("p_value", "1"))
                fdr = float(row.get("fdr", "1"))
                frac_high = float(row.get("frac_high", "0"))
                frac_low = float(row.get("frac_low", "0"))
            except ValueError:
                continue

            relation = ("increases_activity" if direction == "enriched_in_high"
                        else "decreases_activity")
            effect_size = {
                "type": "odds_ratio",
                "log2_odds_ratio": round(math.log2(or_val), 3) if or_val > 0 else None,
                "value": round(or_val, 3),
                "frac_high": round(frac_high, 5),
                "frac_low": round(frac_low, 5),
                "fold_high_vs_low": round(frac_high / frac_low, 2) if frac_low > 0 else None,
            }
            if relation == "increases_activity":
                recommendation = (
                    f"To maximize {cell} activity, include at least one {name} motif "
                    f"(JASPAR {row.get('motif_id')}). Sequences containing it are "
                    f"{or_val:.1f}x more likely to be in the top-activity group "
                    f"({frac_high*100:.1f}% of high vs {frac_low*100:.1f}% of low)."
                )
            else:
                recommendation = (
                    f"Avoid {name} motif in {cell}: sequences containing it are "
                    f"{1/or_val:.1f}x less likely to be in the top-activity group "
                    f"({frac_high*100:.1f}% of high vs {frac_low*100:.1f}% of low)."
                )
            out.append({
                "head": f"motif:{name}",
                "relation": relation,
                "tail": f"cell_line:{cell}",
                "effect_size": effect_size,
                "confidence": {
                    "p_value": p_val,
                    "fdr": fdr,
                    "sample_size": sample_size,
                },
                "recommendation": recommendation,
                "source": SOURCE_FIMO,
                "knowledge_type": "motif_activity",
            })
    return out


# ---------------------------------------------------------------- 类型 2: GC 含量-活性偏好
def _gc_actionable_triplets(
    gc_csv: Path, gc_summary: dict, cell: str
) -> list[dict]:
    """从 gc_analysis CSV + gc_summary 生成行动型 GC-活性三元组。"""
    # 读取所有 bin 用于 percentile 计算
    bins: list[dict] = []
    with open(gc_csv, "r", encoding="utf-8") as f:
        for row in csv.DictReader(f):
            try:
                mean = float(row.get("mean", "nan"))
                n = int(float(row.get("n", "0")))
            except ValueError:
                continue
            if not math.isnan(mean) and n > 0:
                bins.append({"bin": row.get("gc_bin", ""), "mean": mean, "n": n})

    if not bins:
        return []
    # 全局均值 = n 加权平均
    total_n = sum(b["n"] for b in bins)
    global_mean = sum(b["mean"] * b["n"] for b in bins) / total_n if total_n else 0.0
    # percentile: mean 在所有 bin 中的百分位（升序）
    sorted_means = sorted(b["mean"] for b in bins)
    n_bins = len(sorted_means)

    cell_result = gc_summary.get("results", {}).get(cell, {})
    out: list[dict] = []

    def _percentile(mean: float) -> int:
        idx = sum(1 for m in sorted_means if m <= mean) - 1
        return int(round(idx / max(n_bins - 1, 1) * 100))

    for key, relation, template in (
        ("optimal", "optimal_for",
         "Design sequences with GC content {lo_pct:.0f}%-{hi_pct:.0f}% for highest average "
         "activity in {cell} (mean {mean:.2f}, {pct}th percentile vs global {gmean:.2f})."),
        ("suppressive", "suppressive_for",
         "Avoid GC content {lo_pct:.0f}%-{hi_pct:.0f}% in {cell}: mean activity {mean:.2f} "
         "is {pct}th percentile, far below the global mean {gmean:.2f}."),
    ):
        seg = cell_result.get(key)
        if not seg:
            continue
        lo, hi = seg.get("gc_lo", 0), seg.get("gc_hi", 0)
        mean = seg.get("mean", 0.0)
        n_seq = seg.get("n", 0)
        pct = _percentile(mean)
        bin_label = seg.get("gc_bin", f"{lo*100:.0f}%-{hi*100:.0f}%")
        out.append({
            "head": f"gc_range:{lo:.2f}-{hi:.2f}",
            "relation": relation,
            "tail": f"cell_line:{cell}",
            "effect_size": {
                "type": "mean_activity",
                "mean_activity": round(mean, 3),
                "global_mean": round(global_mean, 3),
                "percentile": pct,
                "median": seg.get("median"),
                "diff_ratio_std": cell_result.get("diff_ratio"),
            },
            "confidence": {
                "sample_size": n_seq,
                "n_bins": n_bins,
            },
            "recommendation": template.format(
                lo_pct=lo * 100, hi_pct=hi * 100, cell=cell, mean=mean, pct=pct,
                gmean=global_mean,
            ),
            "source": SOURCE_GC,
            "knowledge_type": "gc_activity",
        })
    return out


# ---------------------------------------------------------------- 文献库: 定性行动型
_LIT_REL_RECOMMEND = {
    "regulates": (
        "{head} regulates {tail} (literature evidence). Include {head} binding "
        "sites when designing sequences to control {tail} expression."
    ),
    "binds_motif": (
        "{head} binds motif {tail} (literature evidence). Including this motif in "
        "a candidate sequence creates a binding opportunity for {head}."
    ),
    "associated_with": (
        "{head} is associated with {tail} (literature evidence). Relevant if the "
        "target phenotype involves {tail}."
    ),
}


# 这些关系属于 MPRA 统计库（已由类型1/类型2 生成行动型版本），文献库升级时跳过
_MPRA_STAT_RELS = {"enriched_in_high", "depleted_in_high",
                   "optimal_range_for", "suppressive_range_for"}


def _literature_actionable_triplets(path: Path) -> list[dict]:
    """文献三元组 -> 定性行动型（无数值效应量，标注 evidence_type）。

    跳过 MPRA 统计关系（enriched_in_high 等），避免与类型1/2 重复。
    """
    out: list[dict] = []
    with open(path, "r", encoding="utf-8") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            try:
                t = json.loads(line)
            except json.JSONDecodeError:
                continue
            head, rel, tail = t.get("head", ""), t.get("relation", ""), t.get("tail", "")
            if not head or not tail:
                continue
            if rel in _MPRA_STAT_RELS:
                continue  # 已由类型1/类型2 生成行动型版本
            template = _LIT_REL_RECOMMEND.get(rel)
            recommendation = (
                template.format(head=head, tail=tail) if template
                else f"{head} {rel} {tail} (literature evidence)."
            )
            out.append({
                "head": head,
                "relation": rel,
                "tail": tail,
                "effect_size": None,
                "confidence": {"evidence_type": "literature", "source": "PubMed"},
                "recommendation": recommendation,
                "source": SOURCE_LIT,
                "knowledge_type": "literature",
            })
    return out


def main() -> None:
    parser = argparse.ArgumentParser(description="升级知识库为行动型三元组")
    parser.add_argument("--outdir", default="data/processed")
    parser.add_argument("--output", default="data/processed/knowledge_base_actionable.jsonl")
    parser.add_argument("--literature", default="data/processed/knowledge_base_full.jsonl")
    args = parser.parse_args()

    outdir = ROOT / args.outdir
    output_path = ROOT / args.output

    split_path = outdir / "split_summary.json"
    split = {}
    if split_path.exists():
        split = json.load(open(split_path, encoding="utf-8"))

    gc_summary_path = outdir / "gc_summary.json"
    gc_summary = json.load(open(gc_summary_path, encoding="utf-8")) if gc_summary_path.exists() else {}

    triplets: list[dict] = []
    n_motif = n_gc = 0

    # 类型 1 + 类型 2（MPRA 统计）
    for cell in DEFAULT_CELLS:
        motif_csv = outdir / f"motif_enrichment_{cell}.csv"
        gc_csv = outdir / f"gc_analysis_{cell}.csv"
        if motif_csv.exists():
            ts = _motif_actionable_triplets(motif_csv, cell, split.get(cell, {}))
            n_motif += len(ts)
            triplets.extend(ts)
        if gc_csv.exists():
            ts = _gc_actionable_triplets(gc_csv, gc_summary, cell)
            n_gc += len(ts)
            triplets.extend(ts)

    # 文献库（定性）
    lit_path = ROOT / args.literature
    lit_ts = _literature_actionable_triplets(lit_path) if lit_path.exists() else []
    triplets.extend(lit_ts)

    output_path.parent.mkdir(parents=True, exist_ok=True)
    with open(output_path, "w", encoding="utf-8") as f:
        for t in triplets:
            f.write(json.dumps(t, ensure_ascii=False) + "\n")

    logger.info("=" * 60)
    logger.info("行动型知识库构建统计:")
    logger.info("  类型1 基序-活性效应  : %d", n_motif)
    logger.info("  类型2 GC-活性偏好    : %d", n_gc)
    logger.info("  文献定性三元组       : %d", len(lit_ts))
    logger.info("  三元组总数           : %d", len(triplets))
    rel_counter: Counter = Counter(t["relation"] for t in triplets)
    logger.info("  关系类型分布         : %s", dict(rel_counter))
    logger.info("  输出文件             : %s", output_path)
    logger.info("=" * 60)

    # 展示样例
    for t in triplets[:5]:
        logger.info("样例: %s", json.dumps(t, ensure_ascii=False)[:400])


if __name__ == "__main__":
    main()
