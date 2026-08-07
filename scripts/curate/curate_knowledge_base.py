# -*- coding: utf-8 -*-
"""
知识库清洗与去重脚本（curate step）

输入: data/processed/knowledge_base_stat.jsonl (Module 1 原始输出, 633 条)
输出:
  - data/processed/knowledge_base_clean.jsonl  清洗后的知识库（新文件，不覆盖原文件）
  - data/processed/curation_report.jsonl        逐条清洗动作日志（复核清单）
  - data/processed/curation_summary.json        汇总统计

清洗内容:
  [1] 去重
      - 精确字符串去重（claim 完全一致）→ 保留 confidence 最高
      - 实体指纹去重（按记录类型分派: 文献/TF基序/DNA shape/GC/其他，
        指纹 + 证据来源(PMID) + 方向一致才去重；不同 PMID 是独立证据不删除）
      - 命名变体去重（仅文献记录: 同 PMID + 同方向 + factor 互为大小写/连字符变体 → 合并，
        如 PGC1A-En1 vs Pgc1a-En1）
  [2] 清洗
      - 必填字段检查: claim / claim_type / entities / evidence.source / confidence
      - 实体规范化: 基因名大小写（仅限确认的基因列表，避免误改非基因实体）
      - claim_type 一致性校验（association/design_rule/mechanistic）
      - 证据合理性: p_value ∈ [0,1]
  [3] 方向语义修复（人工复核确认的抽取错误）
      - idx=206/207: DNA methylation/promoter 机制为抑制(represses LEF1) 却标 increases → 改 decreases
      - idx=474: Osa "required to constrain"(负调控) 标 required_for → 改 decreases + 上下文标注
      - idx=476: Osa no_effect 保留 + 上下文标注（区分调控层级, 不合并）
  [4] 上下文标注增强
      - 对"方向看似矛盾实为不同调控层级"的记录补充 context 字段
      - 应用: Osa(474/476), MYCN 分靶基因组(361/390/394, 371/393/410)
  [4.5] 互补关系识别（context-dependent 原则）
      - 同 PMID 内同一 factor/元件出现不同效应 → 先挖底层逻辑（关系对象/测量层级不同）
      - 能解释 → 标注 context 区分 biological axis（互补非矛盾, 如 enhancer/41617689, rs2294510/38365720）
      - 白名单未覆盖 → 输出 complementary_unclassified 供人工复核
  [5] 一致性验证（输出到 report）
      - 同指纹方向冲突检测（排除已标注 context 的互补组）
      - 统计部分 MPRA 记录: p_value/effect_size 有效性

用法:
  python scripts/curate/curate_knowledge_base.py
  python scripts/curate/curate_knowledge_base.py --input data/processed/knowledge_base_stat.jsonl --output data/processed/knowledge_base_clean.jsonl
"""
from __future__ import annotations

import argparse
import json
import logging
import re
import sys
from collections import defaultdict
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent.parent
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

logging.basicConfig(level=logging.INFO, format="%(levelname)s %(message)s")
logger = logging.getLogger("curate")

# ---------------- 配置 ----------------
VALID_CLAIM_TYPES = {"association", "design_rule", "mechanistic"}
REQUIRED_FIELDS = ["claim", "claim_type", "entities", "evidence", "confidence"]
REQUIRED_EVIDENCE = ["source"]

# 真基因名规范化映射（人工确认的基因，勿加非基因实体）
GENE_NORMALIZATION = {
    "Bclaf1": "BCLAF1",
    "FoxO1": "FOXO1",
    "Sp3": "SP3",
    "FoxA1": "FOXA1",
    "Spi1": "SPI1",
    "Nrf1": "NRF1",
    "Nfe2l2": "NFE2L2",
    "Ets-1": "ETS1",
    "Smyd2": "SMYD2",
    "hemogen": "HEMOGEN",
}

# 方向语义修复（factor, regulatory_element, source 关键词 → 修正)
DIRECTION_FIXES = [
    {
        "match": lambda e, ev: e.get("factor") == "DNA methylation"
        and e.get("regulatory_element") == "promoter"
        and "40918098" in ev.get("source", ""),
        "effect_new": "decreases",
        "direction_new": "decreases_activity",
        "reason": "PMID 40918098: 机制为'represses LEF1 expression'(甲基化↑→LEF1↓), 原标 increases 为抽取错误",
    },
    {
        "match": lambda e, ev: e.get("factor") == "Osa"
        and e.get("regulatory_element") == "enhancer"
        and "37949841" in ev.get("source", ""),
        "effect_new": "decreases",
        "direction_new": "decreases_activity",
        "reason": "PMID 37949841: 'required to constrain the levels of activity' = 负调控(buffer factor), 原标 required_for 语义不准",
    },
]

# 上下文标注增强（direction 看似矛盾实为不同调控层级/靶基因）
# (factor, regulatory_element, source关键词, effect原值) → context 描述
CONTEXT_ENRICHMENT = [
    {
        "match": lambda e, ev: e.get("factor") == "Osa"
        and e.get("regulatory_element") == "enhancer"
        and "37949841" in ev.get("source", ""),
        "context_by_effect": {
            "decreases": "regulatory_axis=enhancer_transcriptional_output; context=brdisc active enhancer (Osa 约束活性水平, 负调控)",
            "no_effect": "regulatory_axis=chromatin_accessibility_transition; context=wing development early vs late (Osa 不影响染色质开放状态变化)",
        },
        "reason": "同 PMID 两条 claim 描述不同调控层级, 非矛盾, 不合并, 仅补充 axis 上下文",
    },
    {
        "match": lambda e, ev: e.get("factor") == "MYCN"
        and e.get("regulatory_element") == "promoter",
        "context_by_effect": {
            "required_for": "target_gene=DLX5; MYCN 直接结合 DLX5 启动子共调控表达 (PMID 41135749)",
            "increases": "target_gene=MYC targets; MYCN 结合增加 H3K4me3/H3K27ac (PMID 37781575/36539767)",
            "decreases": "target_genes=neuronal differentiation genes; MYCN 招募 G9a 抑制 (PMID 37781575)",
        },
        "reason": "MYCN 对不同靶基因方向不同, 补充靶基因上下文避免误判矛盾",
    },
    {
        "match": lambda e, ev: e.get("factor") == "MYCN"
        and e.get("regulatory_element") == "enhancer",
        "context_by_effect": {
            "modulates": "mechanism=MYCN binds enhancers when aberrantly upregulated in neuroblastoma (PMID 36598365)",
            "decreases": "target_genes=neuronal differentiation; MYCN 抑制型增强子 (PMID 37781575)",
            "increases": "target_gene=SMAD9; MYCN 在 enhancer 区域转录调控 SMAD9 (PMID 36539767)",
        },
        "reason": "MYCN/enhancer 不同靶基因/上下文, 补充说明避免误判矛盾",
    },
]

# 互补关系白名单: 同 PMID 内同一 factor/元件出现不同效应（方向看似矛盾），
# 实为 context-dependent（关系对象/测量层级不同）→ 补充 context 区分 biological axis，绝不合并/删除。
# 原则: 生物调控关系通常是 context-dependent 的，先挖掘底层逻辑，再决定保留与标注。
COMPLEMENTARY_PAIRS = [
    {
        "match": lambda e, ev: e.get("factor") == "enhancer"
        and e.get("regulatory_element") == "promoter"
        and "41617689" in ev.get("source", ""),
        "contexts": {
            "required_for": "context=element_combination; promoter 活性对 enhancer 功能是必需但非充分条件 (necessary but not sufficient)",
            "modulates": "context=reciprocal_modulation; 不同 minimal promoter 组合下 enhancer 对 promoter 活性呈双向(reciprocal)调节",
        },
        "reason": "同 PMID 两条 claim: '必需但非充分' vs '组合依赖的双向调节', 属不同测量维度, 互补非矛盾",
    },
    {
        "match": lambda e, ev: e.get("factor") == "rs2294510"
        and e.get("regulatory_element") == "promoter"
        and "38365720" in ev.get("source", ""),
        "contexts": {
            "modulates": "context=allele_specific_effect; rs2294510 改变 TF 结合位点, 效应方向依赖等位基因/细胞上下文",
            "increases": "context=NAFLD_model; 在 HepG2 NAFLD 模型中突变体增强脂质积累表型, 该表型下表现为增强",
        },
        "reason": "同 PMID 两条 claim: '方向依赖上下文' vs '特定表型下增强', 测量层级不同, 互补非矛盾",
    },
]

# 跨 PMID 互补关系: 同 factor/元件在不同文献中方向不同，但底层逻辑可解释为
# context-dependent（细胞类型/测量层级/机制视角不同）→ 补充 context 区分，避免误判为矛盾。
CROSS_PMID_COMPLEMENTARY = [
    {
        "match": lambda e, ev: e.get("factor") == "transcription factor binding motifs"
        and e.get("regulatory_element") == "promoter",
        "contexts": {
            "required_for": "context=motif_functional_necessity; 扰动 TF 结合基序致活性下降, 表明 motif 对 promoter 活性功能性必需 (PMID 41617689)",
            "modulates": "context=variant_driven_modulation; caQTL 变异破坏 TF 结合基序, 通过遗传变异实现调节效应 (PMID 34038741)",
        },
        "reason": "两条 claim 测量层级不同: motif 功能必要性(扰动实验) vs 遗传变异驱动的调节(关联分析), 互补非矛盾",
    },
    {
        "match": lambda e, ev: e.get("factor") == "BRG1"
        and e.get("regulatory_element") == "enhancer",
        "contexts": {
            "required_for": "context=erythroid_enhancers; BRG1 依赖 hemogen 调节红系基因增强子染色质可及性, 该语境下必需 (PMID 35297980)",
            "modulates": "context=basal_activity_inducibility; BRG1 通过影响 cohesin 结合与 enhancer-promoter 互作调节增强子基础活性/可诱导性 (PMID 34260936)",
        },
        "reason": "细胞类型/增强子子集不同: 红系特异(依赖 hemogen) vs 一般增强子基础活性调节, 互补非矛盾",
    },
    {
        "match": lambda e, ev: e.get("factor") == "H3K27ac"
        and e.get("regulatory_element") == "enhancer",
        "contexts": {
            "increases": "context=epigenetic_mark_accumulation; H3K27ac 作为活性增强子表观标志, 其积累与增强活性正相关 (PMID 36698204/37979036)",
            "required_for": "context=functional_necessity; H3K27ac 标记 BRD4 增强转录的区域, 对 BRD4 介导的增强子转录功能必需 (PMID 41984404)",
        },
        "reason": "测量层级不同: 表观标志积累的关联效应 vs 功能必要性(BRD4 介导), 互补非矛盾",
    },
]


# ---------------- 核心逻辑 ----------------


def fingerprint(r: dict) -> tuple:
    """
    实体指纹（按记录类型分派，避免统计记录被空指纹误并）:
      - 文献记录(factor):  LIT / factor / regulatory_element / cell_line / claim_type
      - 统计 TF 基序(tf):  STAT_TF / tf / motif / cell_line / claim_type
      - 统计 DNA shape:    STAT_SHAPE / shape / feature / cell_line / claim_type
      - 统计 GC:           STAT_GC / gc_range 或 (gc_lo,gc_hi) / cell_line / claim_type
      - 其他(家族/画像):   STAT_OTHER / 全实体键值签名 / cell_line / claim_type
    """
    e = r.get("entities", {})
    cell = str(e.get("cell_line") or "")
    ct = r.get("claim_type")
    factor = e.get("factor")
    if factor:  # 文献记录
        return ("LIT", str(factor), str(e.get("regulatory_element") or ""), cell, ct)
    tf = e.get("tf")
    if tf:  # 统计 TF/基序
        return ("STAT_TF", str(tf), str(e.get("motif") or ""), cell, ct)
    shape = e.get("shape")
    if shape:  # 统计 DNA shape
        return ("STAT_SHAPE", str(shape), str(e.get("feature") or ""), cell, ct)
    gc = e.get("gc_range")
    if gc is None and ("gc_lo" in e or "gc_hi" in e):
        gc = (e.get("gc_lo"), e.get("gc_hi"))
    if gc is not None:  # 统计 GC
        return ("STAT_GC", str(gc), cell, ct)
    # 其他（家族/画像/GC分布等）: 全实体键值做签名，避免同 cell 同 claim_type 误并
    sig = tuple(sorted(f"{k}={v}" for k, v in e.items()))
    return ("STAT_OTHER", sig, cell, ct)


def norm_factor(s: str) -> str:
    """实体规范化（仅真基因名，非基因实体原样保留）"""
    if s in GENE_NORMALIZATION:
        return GENE_NORMALIZATION[s]
    return s


def direction_of(r: dict) -> str:
    """归一化方向: 优先 entities.effect, 其次 evidence.direction"""
    e = r.get("entities", {})
    eff = e.get("effect")
    if eff:
        return str(eff).lower()
    d = r.get("evidence", {}).get("direction")
    return str(d).lower() if d else "NA"


def dedupe_exact(rows: list[dict]) -> tuple[list[dict], list[dict]]:
    """精确字符串去重: claim 完全相同 → 保留 confidence 最高"""
    seen: dict[str, dict] = {}
    removed = []
    for r in rows:
        key = r["claim"]
        if key in seen:
            if r.get("confidence", 0) > seen[key].get("confidence", 0):
                removed.append(seen[key])
                seen[key] = r
            else:
                removed.append(r)
        else:
            seen[key] = r
    return list(seen.values()), removed


def dedupe_fingerprint(rows: list[dict]) -> tuple[list[dict], list[dict]]:
    """
    实体指纹去重: 指纹相同 + 证据来源(PMID/MPRA批次)相同 + 方向一致 → 保留 confidence 最高。
    不同 PMID 的记录是独立证据，绝不删除；方向不同（矛盾）的不去重，
    交给互补关系识别/人工。
    """
    groups: dict[tuple, list[dict]] = defaultdict(list)
    for r in rows:
        key = (fingerprint(r), r.get("evidence", {}).get("source", ""))
        groups[key].append(r)
    kept, removed = [], []
    for (fp, src), grp in groups.items():
        if len(grp) < 2:
            kept.extend(grp)
            continue
        by_dir: dict[str, list[dict]] = defaultdict(list)
        for r in grp:
            by_dir[direction_of(r)].append(r)
        for d, sub in by_dir.items():
            if len(sub) > 1:
                best = max(sub, key=lambda x: x.get("confidence", 0))
                kept.append(best)
                removed.extend(x for x in sub if x is not best)
            else:
                kept.extend(sub)
    return kept, removed


def is_naming_variant(a: str, b: str) -> bool:
    """a/b 互为命名变体: 仅大小写/连字符/空格差异（如 PGC1A-En1 vs Pgc1a-En1）"""
    def sig(s: str) -> str:
        return s.lower().replace("-", "").replace(" ", "").replace("_", "")
    return a != b and sig(a) == sig(b)


def dedupe_semantic(rows: list[dict]) -> tuple[list[dict], list[dict]]:
    """
    命名变体去重（仅文献记录）: 同 PMID + 同指纹 + 同方向(effect) + factor 互为命名变体 → 合并。
    已确认案例: PGC1A-En1 vs Pgc1a-En1 (同 PMID 41869877)。
    注意: 同 PMID 内 factor 相同但 effect 不同（如 Osa required_for vs no_effect）是
    context-dependent 互补关系，绝不合并。
    """
    lit = [r for r in rows if r.get("entities", {}).get("factor")]
    stat = [r for r in rows if not r.get("entities", {}).get("factor")]
    groups: dict[tuple, list[dict]] = defaultdict(list)
    for r in lit:
        e = r.get("entities", {})
        f_sig = str(e.get("factor")).lower().replace("-", "").replace(" ", "").replace("_", "")
        key = (
            f_sig,
            str(e.get("regulatory_element") or ""),
            str(e.get("cell_line") or ""),
            r.get("claim_type"),
            r.get("evidence", {}).get("source", ""),
        )
        groups[key].append(r)
    kept, removed = [], []
    for grp in groups.values():
        if len(grp) < 2:
            kept.extend(grp)
            continue
        # 组内再按效应细分: 效应不同绝不合并
        by_eff: dict[str, list[dict]] = defaultdict(list)
        for r in grp:
            by_eff[direction_of(r)].append(r)
        for eff, sub in by_eff.items():
            if len(sub) > 1:
                factors = [str(r.get("entities", {}).get("factor")) for r in sub]
                if all(is_naming_variant(factors[0], f) for f in factors[1:]):
                    best = max(sub, key=lambda x: x.get("confidence", 0))
                    kept.append(best)
                    removed.extend(x for x in sub if x is not best)
                else:
                    kept.extend(sub)
            else:
                kept.extend(sub)
    return kept + stat, removed


def clean_entity_names(r: dict) -> list[str]:
    """实体名称规范化，返回变更说明"""
    actions = []
    e = r.get("entities", {})
    for field in ("factor", "tf"):
        if field in e and isinstance(e[field], str):
            new = norm_factor(e[field])
            if new != e[field]:
                actions.append(f"实体规范化 {field}: {e[field]} → {new}")
                e[field] = new
    return actions


def fix_direction(r: dict) -> list[str]:
    """方向语义修复（人工复核确认的抽取错误）"""
    actions = []
    e = r.get("entities", {})
    ev = r.get("evidence", {})
    for fix in DIRECTION_FIXES:
        if fix["match"](e, ev):
            old_eff = e.get("effect")
            old_dir = ev.get("direction")
            if old_eff != fix["effect_new"]:
                actions.append(f"方向修复 effect: {old_eff} → {fix['effect_new']} ({fix['reason']})")
                e["effect"] = fix["effect_new"]
            if old_dir != fix["direction_new"]:
                actions.append(f"方向修复 direction: {old_dir} → {fix['direction_new']}")
                ev["direction"] = fix["direction_new"]
            # 修复后若方向变 decreases, 同步修正派生 design_rule 的 claim 文本
            if fix["effect_new"] == "decreases" and r.get("claim_type") == "design_rule":
                if "应包含" in r["claim"]:
                    r["claim"] = r["claim"].replace("应包含", "应避免引入").replace("正相关", "负相关")
                    actions.append("design_rule claim 同步: '应包含/正相关' → '应避免引入/负相关'")
    return actions


def enrich_context(r: dict) -> list[str]:
    """上下文标注增强: 对调控层级/靶基因敏感的记录补充 context 字段"""
    actions = []
    e = r.get("entities", {})
    ev = r.get("evidence", {})
    for c in CONTEXT_ENRICHMENT:
        if c["match"](e, ev):
            eff = direction_of(r)
            ctx = c["context_by_effect"].get(eff)
            if ctx and not e.get("context"):
                e["context"] = ctx
                actions.append(f"上下文标注: {ctx[:80]}...")
    return actions


def mark_complementary(rows: list[dict]) -> list[dict]:
    """
    互补关系识别（context-dependent 原则落地）:
    对同 PMID 内同一 factor/元件出现不同效应的记录对，若能在底层逻辑上解释
    （关系对象/测量层级不同），则补充 context 区分 biological axis，而非当作矛盾删除。
    白名单未覆盖的多效应组 → 输出 unclassified 供人工复核。
    返回 report 动作日志。
    """
    actions = []
    groups: dict[tuple, list[dict]] = defaultdict(list)
    for r in rows:
        e = r.get("entities", {})
        if not e.get("factor"):
            continue
        groups[(r.get("evidence", {}).get("source", ""), str(e.get("factor")), str(e.get("regulatory_element") or ""))].append(r)

    for (src, factor, reg), grp in groups.items():
        effs = {direction_of(r) for r in grp}
        if len(effs) < 2:
            continue
        matched = False
        for pair in COMPLEMENTARY_PAIRS:
            sample = {"entities": {"factor": factor, "regulatory_element": reg}, "evidence": {"source": src}}
            if pair["match"](sample["entities"], sample["evidence"]):
                matched = True
                pair_idx = [r.get("_src_idx") for r in grp]
                for r in grp:
                    eff = direction_of(r)
                    ctx = pair["contexts"].get(eff)
                    if ctx and not r.get("entities", {}).get("context"):
                        r["entities"]["context"] = ctx
                        actions.append({
                            "action": "complementary_pair",
                            "src_idx": r["_src_idx"],
                            "factor": factor, "reg": reg, "source": src,
                            "effect": eff,
                            "context": ctx,
                            "reason": pair["reason"],
                            "pair_idx": pair_idx,
                        })
                break
        if not matched:
            actions.append({
                "action": "complementary_unclassified",
                "factor": factor, "reg": reg, "source": src,
                "effects": sorted(effs),
                "src_idx": [r.get("_src_idx") for r in grp],
                "note": "同 PMID 内同 factor 出现不同效应, 白名单未覆盖, 需人工判断是否 context-dependent 互补",
            })
    return actions


def mark_cross_pmid_complementary(rows: list[dict]) -> list[dict]:
    """
    跨 PMID 互补关系标注（context-dependent 原则的扩展）:
    同 factor/元件（同指纹）在不同 PMID 中方向不同，但底层逻辑可解释为
    context-dependent（细胞类型/测量层级/机制视角不同）→ 按 effect 补充 context，
    从而在冲突检测中被豁免。
    返回 report 动作日志。
    """
    actions = []
    groups: dict[tuple, list[dict]] = defaultdict(list)
    for r in rows:
        groups[fingerprint(r)].append(r)
    for fp, grp in groups.items():
        if fp[0] != "LIT" or len(grp) < 2:
            continue
        dirs = {direction_of(r) for r in grp}
        if len(dirs) < 2:
            continue
        e = grp[0].get("entities", {})
        for rule in CROSS_PMID_COMPLEMENTARY:
            if not rule["match"](e, grp[0].get("evidence", {})):
                continue
            for r in grp:
                eff = direction_of(r)
                ctx = rule["contexts"].get(eff)
                if ctx and not r.get("entities", {}).get("context"):
                    r["entities"]["context"] = ctx
                    actions.append({
                        "action": "complementary_cross_pmid",
                        "src_idx": r["_src_idx"],
                        "factor": e.get("factor"), "reg": e.get("regulatory_element"),
                        "effect": eff,
                        "context": ctx,
                        "reason": rule["reason"],
                        "pair_idx": [x["_src_idx"] for x in grp],
                    })
            break
    return actions


def validate(r: dict) -> list[str]:
    """必填字段/一致性/证据合理性校验，返回问题列表"""
    problems = []
    for f in REQUIRED_FIELDS:
        if f not in r or r[f] in (None, ""):
            problems.append(f"缺少必填字段: {f}")
    if r.get("claim_type") not in VALID_CLAIM_TYPES:
        problems.append(f"claim_type 非法: {r.get('claim_type')}")
    ev = r.get("evidence", {})
    for f in REQUIRED_EVIDENCE:
        if f not in ev or not ev[f]:
            problems.append(f"evidence 缺少: {f}")
    p = ev.get("p_value")
    if p is not None:
        try:
            p = float(p)
            if not (0 <= p <= 1):
                problems.append(f"p_value 越界: {p}")
        except (TypeError, ValueError):
            problems.append(f"p_value 非数值: {p}")
    conf = r.get("confidence")
    if conf is not None and not (0 <= conf <= 1):
        problems.append(f"confidence 越界: {conf}")
    return problems


def main():
    ap = argparse.ArgumentParser(description="知识库清洗与去重")
    ap.add_argument("--input", default=str(ROOT / "data/processed/knowledge_base_stat.jsonl"))
    ap.add_argument("--output", default=str(ROOT / "data/processed/knowledge_base_clean.jsonl"))
    ap.add_argument("--report", default=str(ROOT / "data/processed/curation_report.jsonl"))
    ap.add_argument("--summary", default=str(ROOT / "data/processed/curation_summary.json"))
    args = ap.parse_args()

    inp = Path(args.input)
    rows = [json.loads(l) for l in open(inp, encoding="utf-8")]
    logger.info(f"输入: {inp} ({len(rows)} 条)")

    stats = {"input": len(rows)}
    actions_log = []

    # ---- 1. 去重 ----
    rows, removed_exact = dedupe_exact(rows)
    for r in removed_exact:
        actions_log.append({"action": "dedupe_exact", "removed_idx": r.get("_src_idx"), "claim": r["claim"][:80]})
    stats["dedupe_exact_removed"] = len(removed_exact)

    rows, removed_fp = dedupe_fingerprint(rows)
    for r in removed_fp:
        actions_log.append({"action": "dedupe_fingerprint", "claim": r["claim"][:80]})
    stats["dedupe_fingerprint_removed"] = len(removed_fp)

    rows, removed_sem = dedupe_semantic(rows)
    for r in removed_sem:
        actions_log.append({"action": "dedupe_semantic", "claim": r["claim"][:80]})
    stats["dedupe_semantic_removed"] = len(removed_sem)

    # ---- 2. 清洗 + 修复 + 标注 ----
    n_fix, n_ctx = 0, 0
    for i, r in enumerate(rows):
        r["_src_idx"] = i
        acts = []
        acts += clean_entity_names(r)
        acts += fix_direction(r)
        acts += enrich_context(r)
        problems = validate(r)
        if acts:
            n_fix += 1
            actions_log.append({"action": "clean", "src_idx": r["_src_idx"], "changes": acts})
        if problems:
            actions_log.append({"action": "problem", "src_idx": r["_src_idx"], "problems": problems})
        if r.get("entities", {}).get("context"):
            n_ctx += 1
    stats["records_modified"] = n_fix
    stats["records_with_context"] = n_ctx

    # ---- 2.5 互补关系识别（context-dependent 原则）----
    comp_actions = mark_complementary(rows)
    actions_log.extend(comp_actions)
    stats["complementary_pairs_marked"] = len([a for a in comp_actions if a["action"] == "complementary_pair"])
    stats["complementary_unclassified"] = len([a for a in comp_actions if a["action"] == "complementary_unclassified"])
    # 跨 PMID 互补标注
    cross_actions = mark_cross_pmid_complementary(rows)
    actions_log.extend(cross_actions)
    stats["complementary_cross_pmid"] = len(cross_actions)
    # 重新统计带 context 的记录数（互补标注也会加 context）
    stats["records_with_context"] = sum(1 for r in rows if r.get("entities", {}).get("context"))

    # ---- 3. 一致性验证 ----
    # 同指纹方向冲突检测（排除已标注 context 的互补组）
    groups: dict[tuple, list[dict]] = defaultdict(list)
    for r in rows:
        groups[fingerprint(r)].append(r)
    conflicts = []
    for grp in groups.values():
        if len(grp) < 2:
            continue
        dirs = {direction_of(r) for r in grp}
        if len(dirs) > 1:
            # 只有全部带 context 的组豁免（context 已说明为何方向不同）
            if all(r.get("entities", {}).get("context") for r in grp):
                continue
            conflicts.append({
                "fingerprint": list(grp[0].get("entities", {}).keys()),
                "directions": sorted(dirs),
                "src_idx": [r["_src_idx"] for r in grp],
            })
    stats["direction_conflicts_after"] = len(conflicts)
    for c in conflicts:
        actions_log.append({"action": "conflict_remaining", **c})

    # ---- 4. 输出 ----
    out_rows = []
    for r in rows:
        out = {k: v for k, v in r.items() if not k.startswith("_")}
        out_rows.append(out)
    Path(args.output).parent.mkdir(parents=True, exist_ok=True)
    with open(args.output, "w", encoding="utf-8") as f:
        for r in out_rows:
            f.write(json.dumps(r, ensure_ascii=False) + "\n")
    with open(args.report, "w", encoding="utf-8") as f:
        for a in actions_log:
            f.write(json.dumps(a, ensure_ascii=False) + "\n")

    stats["output"] = len(out_rows)
    stats["confidence_mean"] = round(sum(r.get("confidence", 0) for r in out_rows) / len(out_rows), 4) if out_rows else 0
    stats["principles"] = [
        "生物调控关系通常是 context-dependent 的：关系对象（靶基因/调控元件类型）与测量层级"
        "（活性水平/染色质状态/组合上下文）不同，会导致看似矛盾但实际互补的 claims。"
        "处理原则：不简单合并或删除方向相反的记录，先挖掘底层逻辑（调控层级/靶基因/等位基因"
        "或细胞上下文），能解释的保留并标注 context 字段区分 biological axis；"
        "无法自动判断的列入复核清单（complementary_unclassified / conflict_remaining）供人工裁决。"
    ]
    with open(args.summary, "w", encoding="utf-8") as f:
        json.dump(stats, f, ensure_ascii=False, indent=2)

    logger.info(f"输出: {args.output} ({len(out_rows)} 条)")
    logger.info(f"汇总: {json.dumps(stats, ensure_ascii=False)}")
    if conflicts:
        logger.warning(f"仍有 {len(conflicts)} 组方向冲突（见 report）")


if __name__ == "__main__":
    main()
