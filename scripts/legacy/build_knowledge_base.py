"""[DEPRECATED] 从原始抽取结果构建知识三元组 knowledge_base.jsonl。

⚠️ 已废弃：本脚本产出旧版 head/relation/tail 三元组格式。
项目现统一使用 claim 版本知识库（knowledge_base_stat.jsonl，8 键 schema），
由 scripts/build/build_stat_claims.py + scripts/build/build_literature_claims.py 构建。
旧文件已移至 data/backup/。请勿再使用本脚本产出新知识。

用法（在项目根目录运行）:
    python scripts/legacy/build_knowledge_base.py                          # 默认输入
    python scripts/legacy/build_knowledge_base.py --input data/processed/raw_extractions.jsonl
"""

import argparse
import json
import logging
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
logger = logging.getLogger("build_knowledge_base")


def normalize(ent) -> str:
    """实体标准化：去除空白、统一大小写。"""
    s = str(ent).strip()
    return s.upper() if s else s


def build_triples(record: dict) -> set[tuple[str, str, str]]:
    """根据抽取结果构造 (head, relation, tail) 三元组。

    规则:
        tf -regulates-> gene
        tf -binds_motif-> motif
        tf -associated_with-> disease
        gene -associated_with-> disease
    """
    triples: set[tuple[str, str, str]] = set()
    tfs = {normalize(e) for e in record.get("tf", []) if e}
    genes = {normalize(e) for e in record.get("gene", []) if e}
    motifs = {normalize(e) for e in record.get("motif", []) if e}
    diseases = {normalize(e) for e in record.get("disease", []) if e}

    for tf in tfs:
        for g in genes:
            triples.add((tf, "regulates", g))
        for m in motifs:
            triples.add((tf, "binds_motif", m))
        for d in diseases:
            triples.add((tf, "associated_with", d))
    for g in genes:
        for d in diseases:
            triples.add((g, "associated_with", d))

    return triples


def load_flagged_entities(flagged_path: Path) -> dict[str, set[str]]:
    """加载 flagged 清单（review_extractions.py 产物），返回 pmid -> 需排除的实体集合。

    只采用 reason="not_in_abstract"（实体不在原文 = 疑似幻觉）的实体；
    黑名单实体本就不该进知识库，由本模块内置过滤处理。
    """
    result: dict[str, set[str]] = {}
    if not flagged_path.exists():
        logger.warning("flagged 清单 %s 不存在，跳过幻觉过滤", flagged_path)
        return result
    with open(flagged_path, "r", encoding="utf-8") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            try:
                item = json.loads(line)
            except json.JSONDecodeError:
                continue
            if item.get("reason") != "not_in_abstract":
                continue
            pmid = str(item.get("pmid"))
            ent = str(item.get("entity", "")).strip()
            if pmid and ent:
                result.setdefault(pmid, set()).add(normalize(ent))
    logger.info("  幻觉过滤: %d 个 pmid 的 %d 个实体将被排除",
                len(result), sum(len(v) for v in result.values()))
    return result


# 泛化词黑名单: 命中的实体不进知识库（与 review_extractions.BLACKLIST_WORDS 一致）
BLACKLIST_WORDS = [
    "MASTER REGULATORS", "TRANSCRIPTION FACTORS", "TFS", "REGULATORS",
    "ACTIVATORS", "SILENCERS", "BHLH", "HOMEODOMAIN", "HD",
    "TRANSFER RNAS", "DNA SEQUENCES", "BINDING SITES", "COORDINATOR",
    "COMPOSITE MOTIFS", "MOTIFS", "DISEASE", "DISEASES",
    "AUTOIMMUNE DISEASES", "DISEASE RISK ALLELES", "TN5",
    "BASIC HELIX-LOOP-HELIX",
]
_SPECIAL_EXACT = {"DISEASE", "DISEASES"}


def is_blacklisted(entity_norm: str) -> bool:
    """实体是否命中泛化词黑名单（"disease"/"diseases" 仅精确匹配）。"""
    for w in BLACKLIST_WORDS:
        if w in _SPECIAL_EXACT:
            if entity_norm == w:
                return True
        elif w in entity_norm:
            return True
    return False


def main() -> None:
    parser = argparse.ArgumentParser(description="构建知识三元组")
    parser.add_argument("--input", default="data/processed/raw_extractions.jsonl", help="原始抽取结果")
    parser.add_argument("--output", default="data/processed/knowledge_base.jsonl", help="输出三元组文件")
    parser.add_argument("--flagged", default=None,
                        help="review_extractions.py 的 flagged 清单（含幻觉实体，可选）")
    args = parser.parse_args()

    input_path = Path(args.input)
    if not input_path.exists():
        logger.error("输入文件不存在: %s（请先运行 extract_knowledge.py）", input_path)
        sys.exit(1)
    output_path = Path(args.output)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    records: list[dict] = []
    with open(input_path, "r", encoding="utf-8") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            try:
                records.append(json.loads(line))
            except json.JSONDecodeError:
                logger.warning("跳过无法解析的行: %s", line[:100])

    flagged_entities = load_flagged_entities(Path(args.flagged)) if args.flagged else {}

    all_triples: set[tuple[str, str, str]] = set()
    relation_counter: Counter = Counter()
    pmid_with_triples = 0
    skipped_blacklist = 0
    for rec in records:
        if not rec.get("parsed"):
            continue
        pmid = str(rec.get("pmid", ""))
        # 幻觉过滤: 该 pmid 下被 flagged 的实体直接剔除
        exclude = flagged_entities.get(pmid, set())
        for field in ("tf", "gene", "motif", "disease"):
            rec[field] = [e for e in rec.get(field, [])
                          if normalize(e) not in exclude
                          and not is_blacklisted(normalize(e))]
        triples = build_triples(rec)
        if triples:
            pmid_with_triples += 1
        for t in triples:
            all_triples.add(t)
            relation_counter[t[1]] += 1

    # 按 (head, relation, tail) 排序后写入
    sorted_triples = sorted(all_triples)
    with open(output_path, "w", encoding="utf-8") as f:
        for head, rel, tail in sorted_triples:
            obj = {"head": head, "relation": rel, "tail": tail}
            f.write(json.dumps(obj, ensure_ascii=False) + "\n")

    logger.info("=" * 60)
    logger.info("知识库构建统计:")
    logger.info("  抽取记录总数      : %d", len(records))
    logger.info("  产生三元组的记录  : %d", pmid_with_triples)
    logger.info("  去重后三元组总数  : %d", len(all_triples))
    logger.info("  关系类型分布      : %s", dict(relation_counter))
    logger.info("  输出文件          : %s", output_path)
    logger.info("=" * 60)

    # 展示少量样例
    logger.info("三元组样例（前 10 条）:")
    for head, rel, tail in sorted_triples[:10]:
        logger.info("  (%s, %s, %s)", head, rel, tail)


if __name__ == "__main__":
    main()
