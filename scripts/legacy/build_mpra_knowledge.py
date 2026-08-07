"""[DEPRECATED] MPRA 知识库 · 步骤 4：构建 knowledge_base.jsonl（MPRA 部分）。

⚠️ 已废弃：本脚本产出旧版三元组格式。
项目现统一使用 claim 版本知识库（knowledge_base_stat.jsonl，8 键 schema），
由 scripts/build/build_stat_claims.py + scripts/build/build_literature_claims.py 构建。
旧文件已移至 data/backup/。请勿再使用本脚本。

从步骤 2（motif_enrichment_<CELL>.csv）与步骤 3（gc_summary.json）结果
生成知识三元组，输出 data/processed/knowledge_base_mpra.jsonl。

三元组 schema: {head, relation, tail, source, statistic, cell_line, description}

用法（项目根目录）:
    python scripts/legacy/build_mpra_knowledge.py
    python scripts/legacy/build_mpra_knowledge.py --merge   # 与 PubMed 抽取知识库合并
"""

from __future__ import annotations

import argparse
import json
import logging
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from src.knowledge.triplet_builder import gc_triplets_from_summary, motif_triplets_from_csv  # noqa: E402

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s | %(levelname)s | %(message)s",
    datefmt="%H:%M:%S",
)
logger = logging.getLogger("build_mpra_knowledge")

DEFAULT_CELLS = ["K562", "HepG2", "SK-N-SH"]


def main() -> None:
    parser = argparse.ArgumentParser(description="构建 MPRA 知识三元组")
    parser.add_argument("--outdir", default="data/processed")
    parser.add_argument("--output", default="data/processed/knowledge_base_mpra.jsonl")
    parser.add_argument("--merge", action="store_true",
                        help="与 PubMed 抽取知识库 knowledge_base.jsonl 合并")
    parser.add_argument("--cells", default=None, help="逗号分隔的细胞系子集")
    args = parser.parse_args()

    outdir = ROOT / args.outdir
    cells = DEFAULT_CELLS if not args.cells else [c.strip() for c in args.cells.split(",") if c.strip()]

    triplets: list[dict] = []
    for cell in cells:
        motif_csv = outdir / f"motif_enrichment_{cell}.csv"
        gc_summary = outdir / "gc_summary.json"
        if not motif_csv.exists():
            logger.error("缺少 %s（请先运行步骤2）", motif_csv)
            sys.exit(1)
        if not gc_summary.exists():
            logger.error("缺少 %s（请先运行步骤3）", gc_summary)
            sys.exit(1)
        triplets.extend(motif_triplets_from_csv(motif_csv, cell))
        triplets.extend(gc_triplets_from_summary(gc_summary, cell))

    output_path = ROOT / args.output
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with open(output_path, "w", encoding="utf-8") as f:
        for t in triplets:
            f.write(json.dumps(t, ensure_ascii=False) + "\n")

    logger.info("=" * 60)
    logger.info("MPRA 知识库构建统计:")
    logger.info("  三元组总数        : %d", len(triplets))
    rel_counter: dict[str, int] = {}
    for t in triplets:
        rel_counter[t["relation"]] = rel_counter.get(t["relation"], 0) + 1
    logger.info("  关系类型分布      : %s", rel_counter)
    logger.info("  输出文件          : %s", output_path)
    logger.info("=" * 60)
    for t in triplets[:8]:
        logger.info("  (%s, %s, %s)", t["head"], t["relation"], t["tail"])

    if args.merge:
        pubmed_kb = ROOT / "data/processed/knowledge_base.jsonl"
        if not pubmed_kb.exists():
            logger.warning("PubMed 知识库 %s 不存在，跳过合并", pubmed_kb)
            return
        merged = ROOT / "data/processed/knowledge_base_full.jsonl"
        n_pubmed = 0
        with open(merged, "w", encoding="utf-8") as fout:
            for src in (pubmed_kb, output_path):
                with open(src, "r", encoding="utf-8") as fin:
                    for line in fin:
                        line = line.strip()
                        if not line:
                            continue
                        fout.write(line + "\n")
                        if src == pubmed_kb:
                            n_pubmed += 1
        logger.info("已合并: %s (PubMed %d 条 + MPRA %d 条 = %d 条)",
                    merged, n_pubmed, len(triplets), n_pubmed + len(triplets))


if __name__ == "__main__":
    main()
