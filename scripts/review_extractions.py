"""复查全量抽取结果，定位幻觉与错误模式（复查迭代 Step A）。

用法（在项目根目录 /workspace/zzc/BioDesign-Agent 下运行）:
    python scripts/review_extractions.py
    python scripts/review_extractions.py --input data/processed/raw_extractions.jsonl

检查项:
    1) 一致性检查（幻觉检测）: 每个实体（tf/gene/motif/disease）归一化后必须是
       对应摘要原文的子串，否则标记为 "not_in_abstract"（疑似幻觉）。
    2) 黑名单检查: 实体命中黑名单泛化词（master regulators/bHLH/Coordinator 等）
       → 标记为 "blacklist"（即使原文出现也不该作为命名实体输出）。
    3) 错误模式统计: 解析失败条目、全空抽取条目、各字段实体频率 TOP 榜。

输出:
    - 控制台复查报告
    - data/processed/flagged_extractions.jsonl —— 含疑似问题的条目清单，
      供 Step B（模型二次校验）与 Step C（提示词定点修改后重抽）使用。
"""

import argparse
import json
import logging
import re
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
logger = logging.getLogger("review_extractions")

# 与 prompts.py HARD BLACKLIST 保持一致 + 历次测试暴露的泛化词。
# 命中即认为"不是命名实体"，不管原文是否出现。
BLACKLIST_WORDS = [
    "master regulators",
    "transcription factors",
    "tfs",
    "regulators",
    "activators",
    "silencers",
    "bhlh",
    "homeodomain",
    "hd",
    "transfer rnas",
    "dna sequences",
    "binding sites",
    "coordinator",
    "composite motifs",
    "motifs",
    "disease",  # 仅精确匹配（见 hit_blacklist），避免误伤具体疾病名
    "autoimmune diseases",
    "disease risk alleles",
    "tn5",
    "basic helix-loop-helix",
]

# "disease"/"diseases" 只在实体本身恰为该泛化词时才算黑名单，
# 不拦截 "Alzheimer's disease" 这类带具体病名的实体。
_SPECIAL_EXACT = {"DISEASE", "DISEASES"}

ENTITY_FIELDS = ("tf", "gene", "motif", "disease")

_TAG_RE = re.compile(r"<[^>]+>")  # HTML 标签（如 <sup>1-3</sup>）

# 实体别名表: 归一化后的实体名 -> 别名集合。
# 用于消解书写差异: TP53/P53、AP1/AP-1/activator protein 1、NFY/NF-Y、E-box/E box 等。
ALIAS_MAP: dict[str, set[str]] = {
    "TP53": {"P53"},
    "P53": {"TP53"},
    "AP1": {"AP-1", "ACTIVATOR PROTEIN 1"},
    "AP-1": {"AP1", "ACTIVATOR PROTEIN 1"},
    "NFY": {"NF-Y", "NF-YA", "NF-YB"},
    "NF-Y": {"NFY"},
    "HNF1B": {"HNF1A/B", "HNF1"},
    "HNF1": {"HNF1A/B"},
    "SP1": set(),  # 大小写已由 upper() 统一，无需别名
}


def normalize(s: str) -> str:
    """归一化: 去 HTML 标签（删除标签本身，保留内部文字）、大写、压缩空白、去首尾标点。"""
    s = _TAG_RE.sub("", str(s))  # <sub>n</sub> -> n（不保留空格，避免 (AU)n 被拆散）
    s = s.upper()
    s = re.sub(r"\s+", " ", s).strip()
    return s.strip(" .,;:!?\"'()[]{}<>-—–")


def compact(s: str) -> str:
    """紧凑形式: 去掉所有连字符/空格/点号，消除书写差异（E-BOX vs E BOX）。"""
    return re.sub(r"[\s.\-–—]", "", s)


def variants_of(ent_norm: str) -> set[str]:
    """生成实体全部匹配变体: 原形 + 紧凑形 + 别名 + 别名紧凑形。"""
    vs = {ent_norm, compact(ent_norm)}
    for alt in ALIAS_MAP.get(ent_norm, set()):
        vs.add(alt)
        vs.add(compact(alt))
    return vs


def in_abstract(ent_norm: str, abs_norm: str, abs_compact: str) -> bool:
    """实体是否出现在摘要中（支持别名与书写差异）。"""
    for v in variants_of(ent_norm):
        if v and v in abs_norm:
            return True
        cv = compact(v)
        if cv and cv in abs_compact:
            return True
    # 宽松回退: 实体核心词（去括号/去通用修饰后缀）全部出现在摘要中即算匹配。
    # 解决 "NFY binding site" vs 原文 "NFY binding"、"cleft lip with or cleft palate (CL/P)"
    # vs "cleft lip with or without cleft palate (CL/P)" 这类补充词/丢词误报。
    return core_tokens_in_abstract(ent_norm, abs_norm)


# 通用修饰后缀: 作为实体尾部时允许松弛匹配（不算核心语义词）
_GENERIC_SUFFIX = {
    "BINDING SITE", "BINDING SITES", "SITE", "SITES",
    "MOTIF", "MOTIFS", "REGION", "REGIONS",
    "ELEMENT", "ELEMENTS", "ENHANCER", "ENHANCERS",
    "PROMOTER", "PROMOTERS", "BINDING", "FACTOR", "FACTORS",
    "PROTEIN", "PROTEINS", "TFS", "TF", "REPEATS", "REPEAT",
}


def core_tokens_in_abstract(ent_norm: str, abs_norm: str) -> bool:
    """宽松匹配: 实体的核心词（去括号内容/通用后缀）是否全部出现在摘要中。"""
    # 去括号内容（如 (CL/P)）
    s = re.sub(r"\([^)]*\)", "", ent_norm)
    tokens = [t for t in re.split(r"[\s\-–—/:]+", s.strip()) if t]
    if not tokens:
        return False
    # 从尾部去掉通用修饰词（FACTOR/BINDING/SITE/MOTIF 等）
    while tokens and tokens[-1] in _GENERIC_SUFFIX:
        tokens.pop()
    if not tokens:
        return False
    # 所有核心词必须都在摘要中出现
    for t in tokens:
        if t not in abs_norm:
            return False
    return True


def hit_blacklist(entity_norm: str) -> bool:
    """实体（已归一化）是否命中黑名单泛化词。

    "disease"/"diseases" 做精确匹配（避免误伤 "Alzheimer's disease"
    这类具体疾病名），其余黑名单词按子串匹配。
    """
    for w in BLACKLIST_WORDS:
        wu = w.upper()
        if wu in _SPECIAL_EXACT:
            if entity_norm == wu:
                return True
        elif wu in entity_norm:
            return True
    return False


def load_jsonl(path: Path) -> list[dict]:
    items: list[dict] = []
    with open(path, "r", encoding="utf-8") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            try:
                items.append(json.loads(line))
            except json.JSONDecodeError:
                logger.warning("跳过无法解析的行: %s", line[:100])
    return items


def main() -> None:
    parser = argparse.ArgumentParser(description="复查抽取结果，定位幻觉与错误模式")
    parser.add_argument("--input", default="data/processed/raw_extractions.jsonl")
    parser.add_argument("--abstracts", default="data/raw/abstracts.jsonl", help="原始摘要（一致性检查用）")
    parser.add_argument("--flagged", default="data/processed/flagged_extractions.jsonl", help="问题清单输出")
    args = parser.parse_args()

    input_path = Path(args.input)
    abs_path = Path(args.abstracts)
    if not input_path.exists() or not abs_path.exists():
        logger.error("输入文件不存在（请先运行 extract_knowledge.py 并在项目根目录执行）")
        sys.exit(1)

    records = load_jsonl(input_path)
    abstracts = load_jsonl(abs_path)
    abs_by_pmid = {str(a.get("pmid")): a.get("abstract", "") for a in abstracts}
    logger.info("共加载 %d 条抽取记录、%d 条摘要", len(records), len(abs_by_pmid))

    n_parsed = sum(1 for r in records if r.get("parsed"))
    n_failed = len(records) - n_parsed

    # ---------- 检查 ----------
    not_in_abs: Counter = Counter()          # (entity, field) -> count
    blacklist_hits: Counter = Counter()      # (entity, field) -> count
    field_freq: dict[str, Counter] = {f: Counter() for f in ENTITY_FIELDS}
    flagged: list[dict] = []
    empty_records: list[str] = []
    empty_rel_text: list[str] = []
    total_entities = 0

    for rec in records:
        pmid = str(rec.get("pmid"))
        abstract = abs_by_pmid.get(pmid, "")
        abs_norm = normalize(abstract)
        abs_compact = compact(abs_norm)

        if not rec.get("parsed"):
            flagged.append({
                "pmid": pmid, "reason": "parse_failed",
                "raw_output": str(rec.get("raw_output", ""))[:500],
            })
            continue

        all_empty = True
        for field in ENTITY_FIELDS:
            for ent in rec.get(field, []) or []:
                ent = str(ent)
                if not ent:
                    continue
                all_empty = False
                total_entities += 1
                ent_norm = normalize(ent)
                field_freq[field][ent_norm] += 1

                # 1) 黑名单检查
                if hit_blacklist(ent_norm):
                    blacklist_hits[(ent, field)] += 1
                    flagged.append({
                        "pmid": pmid, "field": field, "entity": ent,
                        "reason": "blacklist",
                        "abstract": abstract[:800],
                    })
                    continue

                # 2) 一致性检查（幻觉检测）——支持别名/紧凑形匹配
                if not ent_norm or not in_abstract(ent_norm, abs_norm, abs_compact):
                    not_in_abs[(ent, field)] += 1
                    flagged.append({
                        "pmid": pmid, "field": field, "entity": ent,
                        "reason": "not_in_abstract",
                        "abstract": abstract[:800],
                    })

        if all_empty:
            empty_records.append(pmid)
            empty_rel_text.append(str(rec.get("relation", "")).strip())

    # ---------- 报告 ----------
    print("=" * 64)
    print("复查报告: %s" % input_path)
    print("=" * 64)
    print("总记录数      : %d" % len(records))
    print("解析成功      : %d (%.1f%%)" % (n_parsed, n_parsed / len(records) * 100 if records else 0))
    print("解析失败      : %d" % n_failed)
    print("实体总数      : %d" % total_entities)
    print("全空抽取条目  : %d" % len(empty_records))
    n_with_reason = sum(1 for t in empty_rel_text if t)
    if empty_records:
        print("  全空但有说明  : %d (%.0f%%)" % (n_with_reason, n_with_reason / len(empty_records) * 100))
        print("  空条目 pmid   : %s" % ", ".join(empty_records[:30]))
        print("  说明抽样      :")
        shown = 0
        for t in empty_rel_text:
            if t:
                print("      %s" % t[:150])
                shown += 1
                if shown >= 5:
                    break

    print("-" * 64)
    print("疑似幻觉（实体不在原文）: %d 个" % sum(not_in_abs.values()))
    for (ent, field), c in not_in_abs.most_common(30):
        print("  [%s] %s  ×%d" % (field, ent, c))

    print("-" * 64)
    print("黑名单命中（泛化词当实体）: %d 个" % sum(blacklist_hits.values()))
    for (ent, field), c in blacklist_hits.most_common(30):
        print("  [%s] %s  ×%d" % (field, ent, c))

    print("-" * 64)
    for field in ENTITY_FIELDS:
        print("字段 %-7s 实体频率 TOP15:" % field)
        for ent, c in field_freq[field].most_common(15):
            print("    %-40s ×%d" % (ent, c))

    # ---------- 写 flagged 清单 ----------
    flagged_path = Path(args.flagged)
    flagged_path.parent.mkdir(parents=True, exist_ok=True)
    with open(flagged_path, "w", encoding="utf-8") as f:
        for item in flagged:
            f.write(json.dumps(item, ensure_ascii=False) + "\n")
    print("=" * 64)
    print("问题清单已写入: %s（共 %d 条，供模型二次校验/重抽）" % (flagged_path, len(flagged)))


if __name__ == "__main__":
    main()
