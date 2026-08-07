"""文献知识源（Source 2）· 科学结论 claims 构建（写入 knowledge_base_stat.jsonl）。

从 Qwen 抽取结果 activity_extractions.jsonl（findings 列表）生成统一 schema
的 claim 版本知识，与统计知识（scripts/build/build_stat_claims.py 产出）写入同一文件
knowledge_base_stat.jsonl（追加，保留原有统计 claims）。

每条 finding 生成 1 条 association/mechanistic claim；对方向明确的
（increases/required_for/decreases）额外生成 1 条 design_rule claim。
同一 (cell_line, regulatory_element, factor, effect) 跨多篇文献时合并
证据（pmids + evidence_sentences）。

Claim schema（与统计部分一致）:
{
    "claim": str,                    # 完整、可独立理解的科学结论句子
    "claim_type": str,               # association | mechanistic | design_rule
    "entities": {...},               # factor / cell_line / regulatory_element / effect / pmids
    "evidence": {
        "source": str,               # 文献来源（PubMed PMIDs）
        "direction": str             # increases / required_for / decreases ...
    },
    "reasoning_chain": [str, ...],   # 文献证据链
    "experimental_context": {
        "assay": "literature",
        "cell_type": str,
        "source_type": "PubMed"
    },
    "confidence": float              # 0-1
}

用法（项目根目录）:
    python scripts/build/build_literature_claims.py
    python scripts/build/build_literature_claims.py --input data/processed/activity_extractions.jsonl
    python scripts/build/build_literature_claims.py --output data/processed/knowledge_base_stat.jsonl
"""

from __future__ import annotations

import argparse
import json
import logging
import re
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s | %(levelname)s | %(message)s",
    datefmt="%H:%M:%S",
)
logger = logging.getLogger("build_literature_claims")

# 与统计部分保持一致的中文细胞系描述
CELL_CN = {
    "K562": "K562（红系前体细胞）",
    "HepG2": "HepG2（肝癌细胞）",
    "SK-N-SH": "SK-N-SH（神经母细胞瘤细胞）",
    "not_specified": "未指定细胞系",
}

# Qwen 抽取中常见的 cell_line 别名/变体 → 规范名
CELL_ALIASES = {
    "SH-SY5Y": "SK-N-SH",          # SH-SY5Y 是 SK-N-SH 的亚系
    "neuroblastoma": "SK-N-SH",
    "HepG2-NTCP": "HepG2",         # HepG2 的稳定转染亚系
    "k562": "K562",
    "hepg2": "HepG2",
    "sk-n-sh": "SK-N-SH",
}

# Qwen 抽取中 effect 的同义词/变体 → 受控词表
EFFECT_SYNONYMS = {
    "increases": "increases",
    "increased": "increases",
    "enhances": "increases",
    "promotes": "increases",
    "upregulates": "increases",
    "upregulated": "increases",
    "activates": "increases",
    "transcriptionally activates": "increases",
    "assists": "increases",
    "required for": "required_for",
    "required_for": "required_for",
    "is required for": "required_for",
    "needed for": "required_for",
    "essential for": "required_for",
    "decreases": "decreases",
    "decreased": "decreases",
    "reduces": "decreases",
    "reduced": "decreases",
    "represses": "decreases",
    "suppresses": "decreases",
    "inhibits": "decreases",
    "modulates": "modulates",
    "regulates": "modulates",
    "correlates_with": "modulates",
    "correlates with": "modulates",
    "binds_to": "modulates",
    "no_effect": "no_effect",
    "no effect": "no_effect",
}

# effect 受控词表 → 中文描述 / 方向标签
EFFECT_CN = {
    "increases": ("增强（increases）", "increases_activity"),
    "decreases": ("降低（decreases）", "decreases_activity"),
    "required_for": ("所必需的（required_for）", "required_for_activity"),
    "modulates": ("调节（modulates，方向依赖上下文）", "modulates_activity"),
    "no_effect": ("无显著影响（no_effect）", "no_effect"),
}

# effect → 自然句式（嵌入 claim，{element} 为调控元件名）
EFFECT_PHRASE = {
    "increases": "增强（increases）{element}活性",
    "decreases": "降低（decreases）{element}活性",
    "required_for": "对{element}活性是必需的（required_for）",
    "modulates": "调节（modulates）{element}活性（方向依赖上下文）",
    "no_effect": "对{element}活性无显著影响（no_effect）",
}

EFFECT_CONF = {
    "increases": 0.90,
    "required_for": 0.92,
    "decreases": 0.90,
    "modulates": 0.85,
    "no_effect": 0.80,
}

# 因子类型启发式推断（用于 entities.factor_type 与 design_rule 措辞）
_TF_RE = re.compile(r"^[A-Z][A-Z0-9]{1,7}(?:::[A-Z0-9]+)?$")
_SEQ_RE = re.compile(r"motif|GC|shape|spacing|position|sequence|nucleotid", re.I)
_EPI_RE = re.compile(r"H3K|methylation|chromatin|nucleosome|accessib|acetylation", re.I)


def infer_factor_type(factor: str) -> str:
    """粗略推断因子类型：TF / motif / sequence_feature / epigenetic / other。"""
    if not factor:
        return "other"
    if _EPI_RE.search(factor):
        return "epigenetic"
    if _TF_RE.match(factor.strip()) and len(factor.strip()) <= 12:
        return "TF"
    if _SEQ_RE.search(factor):
        return "sequence_feature"
    return "other"


def normalize_finding(f: dict) -> dict | None:
    """清洗单条 finding：字段规范化、丢弃缺失关键字段的条目。"""
    factor = str(f.get("factor") or "").strip()
    effect_raw = str(f.get("effect") or "").strip().lower()
    # effect 同义词归一化到受控词表；无法归一化则丢弃该条
    effect = EFFECT_SYNONYMS.get(effect_raw, "").strip()
    if not factor or effect not in EFFECT_CN:
        return None
    cell_line = str(f.get("cell_line") or "not_specified").strip()
    # cell_line 别名归一化；无法识别则回退 not_specified
    if cell_line in CELL_ALIASES:
        cell_line = CELL_ALIASES[cell_line]
    elif cell_line not in CELL_CN and cell_line.lower() not in {k.lower() for k in CELL_CN}:
        cell_line = "not_specified"
    element = str(f.get("regulatory_element") or "regulatory element").strip().lower()
    if not element:
        element = "regulatory element"
    return {
        "cell_line": cell_line,
        "regulatory_element": element,
        "factor": factor,
        "effect": effect,
        "mechanism": str(f.get("mechanism") or "").strip(),
        "evidence_sentence": str(f.get("evidence_sentence") or "").strip(),
    }


def load_extractions(path: Path) -> list[dict]:
    """读取抽取结果，展平 findings 并按 (cell, element, factor, effect) 合并。"""
    merged: dict[tuple, dict] = {}
    with open(path, "r", encoding="utf-8") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            try:
                rec = json.loads(line)
            except json.JSONDecodeError:
                continue
            pmid = str(rec.get("pmid") or "")
            for finding in rec.get("findings") or []:
                nf = normalize_finding(finding)
                if nf is None:
                    continue
                key = (nf["cell_line"], nf["regulatory_element"], nf["factor"], nf["effect"])
                if key not in merged:
                    nf["pmids"] = [pmid] if pmid else []
                    nf["evidence_sentences"] = [nf.pop("evidence_sentence")] if nf["evidence_sentence"] else []
                    nf["mechanisms"] = [nf.pop("mechanism")] if nf["mechanism"] else []
                    merged[key] = nf
                else:
                    e = merged[key]
                    if pmid and pmid not in e["pmids"]:
                        e["pmids"].append(pmid)
                    if nf["evidence_sentence"] and nf["evidence_sentence"] not in e["evidence_sentences"]:
                        e["evidence_sentences"].append(nf["evidence_sentence"])
                    if nf["mechanism"] and nf["mechanism"] not in e["mechanisms"]:
                        e["mechanisms"].append(nf["mechanism"])
    return list(merged.values())


def fmt_pmids(pmids: list[str]) -> str:
    return "、".join(f"PMID {p}" for p in pmids) if pmids else "文献未注明PMID"


def build_association_claim(nf: dict) -> dict:
    """association/mechanistic claim：factor 对 cell 中 element 活性的影响。"""
    cell_cn = CELL_CN[nf["cell_line"]]
    _, direction = EFFECT_CN[nf["effect"]]
    factor_type = infer_factor_type(nf["factor"])
    mechanisms = nf.get("mechanisms") or []
    sentences = nf.get("evidence_sentences") or []

    if mechanisms:
        claim_type = "mechanistic"
        mech_part = "；机制上，" + mechanisms[0] if mechanisms[0] else ""
    else:
        claim_type = "association"
        mech_part = ""

    effect_phrase = EFFECT_PHRASE[nf["effect"]].format(element=nf["regulatory_element"])
    claim = (
        f"文献证据表明，在{cell_cn}的实验中，因子{nf['factor']}{effect_phrase}"
        f"{mech_part}（来源：{fmt_pmids(nf['pmids'])}）。"
    )

    chain = [
        f"来源: {fmt_pmids(nf['pmids'])}",
    ]
    if sentences:
        chain.append("摘要原文证据: " + sentences[0])
    if mechanisms:
        chain.append("文献机制描述: " + mechanisms[0])
    chain.append(f"因子类型判定: {factor_type}；细胞系: {nf['cell_line']}；元件: {nf['regulatory_element']}")

    evidence = {
        "source": f"PubMed {'/'.join(nf['pmids'])}" if nf["pmids"] else "PubMed",
        "direction": direction,
        "p_value": None,
        "effect_size": None,
    }

    return {
        "claim": claim,
        "claim_type": claim_type,
        "entities": {
            "factor": nf["factor"],
            "factor_type": factor_type,
            "cell_line": nf["cell_line"],
            "regulatory_element": nf["regulatory_element"],
            "effect": nf["effect"],
            "pmids": nf["pmids"],
        },
        "evidence": evidence,
        "reasoning_chain": chain,
        "experimental_context": {
            "assay": "literature",
            "cell_type": nf["cell_line"],
            "source_type": "PubMed",
        },
        "confidence": EFFECT_CONF.get(nf["effect"], 0.85),
    }


def build_design_rule_claim(nf: dict) -> dict | None:
    """design_rule claim：对方向明确的因子给出设计建议。"""
    if nf["effect"] not in ("increases", "required_for", "decreases"):
        return None
    cell_cn = CELL_CN[nf["cell_line"]]
    if nf["effect"] in ("increases", "required_for"):
        action = "应包含"
        reason = "因其与高活性正相关"
    else:
        action = "应避免引入"
        reason = "因其与高活性负相关"
    factor_type = infer_factor_type(nf["factor"])

    claim = (
        f"设计建议（文献增强）：在设计{cell_cn}的{nf['regulatory_element']}序列时，{action}"
        f"{nf['factor']}相关特征，{reason}（来源：{fmt_pmids(nf['pmids'])}）。"
    )

    chain = [
        f"由文献关联推导: {fmt_pmids(nf['pmids'])}",
        f"原始发现: {EFFECT_PHRASE[nf['effect']].format(element=nf['regulatory_element'])}",
        "需功能实验验证（文献结论的工程化转写）",
    ]

    return {
        "claim": claim,
        "claim_type": "design_rule",
        "entities": {
            "factor": nf["factor"],
            "factor_type": factor_type,
            "cell_line": nf["cell_line"],
            "regulatory_element": nf["regulatory_element"],
            "effect": nf["effect"],
            "pmids": nf["pmids"],
        },
        "evidence": {
            "source": f"PubMed {'/'.join(nf['pmids'])}" if nf["pmids"] else "PubMed",
            "direction": EFFECT_CN[nf["effect"]][1],
            "p_value": None,
            "effect_size": None,
        },
        "reasoning_chain": chain,
        "experimental_context": {
            "assay": "literature",
            "cell_type": nf["cell_line"],
            "source_type": "PubMed",
        },
        "confidence": EFFECT_CONF.get(nf["effect"], 0.85) - 0.05,
    }


def main() -> None:
    parser = argparse.ArgumentParser(description="构建文献知识源 claims 并追加到 knowledge_base_stat.jsonl")
    parser.add_argument("--input", default="data/processed/activity_extractions.jsonl")
    parser.add_argument("--output", default="data/processed/knowledge_base_stat.jsonl")
    parser.add_argument("--no-design-rules", action="store_true", help="不生成 design_rule claims")
    args = parser.parse_args()

    input_path = ROOT / args.input
    if not input_path.exists():
        logger.error("输入文件不存在: %s（请先运行 extract_activity_knowledge.py）", input_path)
        sys.exit(1)

    findings = load_extractions(input_path)
    logger.info("合并后独立 finding: %d 条", len(findings))
    if not findings:
        logger.error("没有可用的 findings，退出")
        sys.exit(1)

    claims: list[dict] = []
    for nf in findings:
        claims.append(build_association_claim(nf))
        if not args.no_design_rules:
            dr = build_design_rule_claim(nf)
            if dr:
                claims.append(dr)

    # 写入：追加到既有 knowledge_base_stat.jsonl（保留统计 claims）
    output_path = ROOT / args.output
    existed = output_path.exists() and output_path.stat().st_size > 0
    n_before = 0
    if existed:
        with open(output_path, "r", encoding="utf-8") as f:
            n_before = sum(1 for line in f if line.strip())

    with open(output_path, "a", encoding="utf-8") as f:
        for c in claims:
            f.write(json.dumps(c, ensure_ascii=False) + "\n")

    # 统计
    from collections import Counter
    type_counter = Counter(c["claim_type"] for c in claims)
    cell_counter = Counter(c["entities"]["cell_line"] for c in claims)
    effect_counter = Counter(c["entities"]["effect"] for c in claims)
    logger.info("=" * 70)
    logger.info("文献知识源 claims 构建统计:")
    logger.info("  新增 claims      : %d", len(claims))
    logger.info("  类型分布         : %s", dict(type_counter))
    logger.info("  细胞系分布       : %s", dict(cell_counter))
    logger.info("  效应方向分布     : %s", dict(effect_counter))
    logger.info("  输出文件         : %s（原 %d 条 → 现 %d 条）", output_path, n_before, n_before + len(claims))
    logger.info("=" * 70)
    for c in claims[:6]:
        logger.info("  · %s", c["claim"])


if __name__ == "__main__":
    main()
