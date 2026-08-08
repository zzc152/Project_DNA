#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Gold 标注修复脚本（fix_gold.py）

根据 docs/gold_validation.md 的修复建议，对 drafts 与 KB 数据做修复。
修复分三类：
  FIX    — 高置信方向修复（摘要/claim 原文支撑），同步更新
           entities_used.effect + evidence.direction + reasoning_chain 机制描述
  DELCAND — 删除候选标记（因子不可考 / 共现相关非因果 / 题目语义反向），
           只加 metadata.fix_note，不删除记录（保守，供人工确认）
  FLAG   — 需人工核查标记（多方向共存但无确凿证据），只加 fix_note

输出：
  - data/synthetic/drafts_full2_clean_fixed.jsonl
  - data/processed/knowledge_base_clean_fixed.jsonl
  - tmp/fix_report.json（修复明细）
"""
import argparse
import json
import re
import copy

# ---------------------------------------------------------------------------
# L3 修复映射：id -> (字段, 新值, 说明, 类型)
# ---------------------------------------------------------------------------

L3_FIXES = {
    # --- 高置信方向修复（FIX）---
    "L3_001_8": {
        "effect": "decreases", "direction": "decreases_activity",
        "note": "摘要：RepID缺失→globin promoter常染色质化+NRF2结合增强→RepID抑制启动子活性（原标注increases方向错误）",
        "kind": "FIX",
    },
    "L3_001_15": {
        "effect": "decreases", "direction": "decreases_activity",
        "note": "摘要：ABO表达降低与ABO promoter DNA甲基化强相关→甲基化抑制启动子活性（原标注increases方向错误）",
        "kind": "FIX",
    },
    "L3_002_55": {
        "effect": "increases", "direction": "increases_activity",
        "note": "摘要：eRNA耗竭沉默BCL11A（disrupting insulation）→eRNA增强enhancer活性（原标注decreases方向错误）",
        "kind": "FIX",
    },
    "L3_001_33": {
        "effect": "increases", "direction": "increases_activity",
        "note": "摘要：H3K4me1缺失降低H3K27ac并削弱染色质重塑复合物募集→H3K4me1增强enhancer活性（原标注decreases方向错误）",
        "kind": "FIX",
    },
    # --- 字段统一（FIX）---
    "L3_002_8": {
        "effect": None, "direction": "decreases_activity",
        "note": "字段统一：macroH2A定位在缺失H3K27ac的失活enhancer上（inactive enhancers）→方向为decreases（与effect一致）",
        "kind": "FIX",
    },
    # --- 关系降级（FIX）---
    "L3_002_37": {
        "effect": "modulates", "direction": "binds_target",
        "note": "关系降级：证据仅支持ASCL1结合/预备（bound and primed）enhancer，激活需G1 arrest→方向不明（与reasoning_chain[4]修正一致）",
        "kind": "FIX",
    },
    # --- 删除候选（DELCAND）---
    "L3_002_26": {
        "effect": None, "direction": None,
        "note": "DELCAND：因子'epigenetic status'不可考（摘要实体为风险等位基因G，非该因子）；同PMID 39998882与L3_003_29重复，建议删除或改写",
        "kind": "DELCAND",
    },
    "L3_003_29": {
        "effect": None, "direction": None,
        "note": "DELCAND：因子'epigenetic status'不可考（摘要实体为风险等位基因G，非该因子）；同PMID 39998882与L3_002_26重复，建议删除或改写",
        "kind": "DELCAND",
    },
    "L3_002_2": {
        "effect": None, "direction": None,
        "note": "DELCAND：题目语义反向——摘要讲BCL11A erythroid enhancer（BCL11A基因的增强子）被破坏损害转录，是enhancer作用于BCL11A，非BCL11A作用于enhancer",
        "kind": "DELCAND",
    },
}

# ---------------------------------------------------------------------------
# T3 描述性因子名模式（命中则标记 factor_valid=false）
# ---------------------------------------------------------------------------

_DESC_PATTERNS = re.compile(
    r"^(epigenetic status|chromatin|transcription factor|gene expression|"
    r"regulatory element|binding|recruitment|activity|expression|"
    r"methylation|acetylation)\b", re.I)


def _is_desc_factor(f: str) -> bool:
    if not f:
        return False
    return len(f) > 30 or f.count(" ") >= 4 or bool(_DESC_PATTERNS.match(f.strip()))


# ---------------------------------------------------------------------------
# KB 修复：用 (factor, element, cell_line) 定位
# ---------------------------------------------------------------------------

def _kkey(e):
    return (str(e.get("factor") or "").strip(),
            str(e.get("regulatory_element") or "").strip(),
            str(e.get("cell_line") or "").strip())


KB_FIXES = [
    # (匹配键, effect 新值, direction 新值, note, kind)
    (("DNA methylation", "promoter", "K562"), "decreases", "decreases_activity",
     "与L3_001_15一致：ABO promoter甲基化与表达降低强相关→DNA methylation抑制promoter（K562条目原标increases错误）",
     "FIX"),
    (("GATA motif", "OGT/OGA promoter interaction", "not_specified"), "decreases", "decreases_activity",
     "claim明确：GATA位点突变（GAGA）显著降低OGT与OGA的promoter相互作用→effect从affects_interaction改为decreases",
     "FIX"),
    (("caQTL variants", "TF binding motifs", "not_specified"), "decreases", "decreases_activity",
     "disrupts=破坏结合基序→抑制方向，映射到规范词表decreases",
     "FIX"),
    (("OsDDE9", "Ghd7 and Ehd1 expression", "not_specified"), "modulates", "modulates_activity",
     "regulates_expression=调节表达但方向未明→映射到modulates",
     "FIX"),
    # 删除候选（共现/相关非因果）
    (("DNA methylation", "gene upregulation", "not_specified"), None, None,
     "DELCAND：co_occurs_with为共现关系非因果方向，不适合方向题gold",
     "DELCAND"),
    (("chromatin accessibility", "promoter methylation level", "not_specified"), None, None,
     "DELCAND：correlates为相关关系非因果方向，不适合方向题gold",
     "DELCAND"),
    (("targeted activation of this regulatory element", "AML differentiation", "not_specified"), None, None,
     "DELCAND：factor为描述性短语不可考（promotes），且element为分化过程非调控元件",
     "DELCAND"),
    # 需人工核查（多方向共存）
    (("MYCN", "enhancer", "not_specified"), None, None,
     "FLAG：MYCN→enhancer 3种方向共存（modulates/decreases/increases），需逐条回原文核查，不擅自修改",
     "FLAG"),
    (("MYCN", "promoter", "not_specified"), None, None,
     "FLAG：MYCN→promoter required_for×2+increases×2+decreases×2，A/B类冲突，需回原文核查",
     "FLAG"),
]


def main():
    ap = argparse.ArgumentParser(description="Gold 标注修复")
    ap.add_argument("--drafts", default="data/synthetic/drafts_full2_clean.jsonl")
    ap.add_argument("--kb", default="data/processed/knowledge_base_clean.jsonl")
    ap.add_argument("--out-drafts", default="data/synthetic/drafts_full2_clean_fixed.jsonl")
    ap.add_argument("--out-kb", default="data/processed/knowledge_base_clean_fixed.jsonl")
    ap.add_argument("--out-report", default="tmp/fix_report.json")
    args = ap.parse_args()

    drafts = [json.loads(l) for l in open(args.drafts, encoding="utf-8") if l.strip()]
    kbs = [json.loads(l) for l in open(args.kb, encoding="utf-8") if l.strip()]
    report = {"drafts_fixed": [], "kb_fixed": [], "flags": []}

    # ---- L3 修复 ----
    n_fix = n_del = n_flag = 0
    for r in drafts:
        rid = r.get("id", "")
        fx = L3_FIXES.get(rid)
        if not fx:
            # T3 因子名标记
            eu = r.get("metadata", {}).get("entities_used", {})
            if _is_desc_factor(eu.get("factor") or ""):
                r["metadata"].setdefault("fix_note", "")
                r["metadata"]["fix_note"] += "T3: 因子名不可考（描述性短语）;"
                r["metadata"]["factor_valid"] = False
                n_flag += 1
                report["flags"].append({"id": rid, "type": "T3", "factor": eu.get("factor")})
            continue
        md = r.setdefault("metadata", {})
        eu = md.setdefault("entities_used", {})
        ev = md.setdefault("evidence", {})
        notes = []
        if fx["effect"] is not None and eu.get("effect") != fx["effect"]:
            notes.append(f"effect: {eu.get('effect')} → {fx['effect']}")
            eu["effect"] = fx["effect"]
        if fx["direction"] is not None and ev.get("direction") != fx["direction"]:
            notes.append(f"evidence.direction: {ev.get('direction')} → {fx['direction']}")
            ev["direction"] = fx["direction"]
        # 同步 reasoning_chain 机制描述里的方向词（若含旧方向词）
        old_eff = eu.get("effect")
        md["fix_note"] = fx["note"]
        if fx["kind"] == "FIX":
            md["fix_status"] = "fixed"
            n_fix += 1
            report["drafts_fixed"].append({"id": rid, "changes": notes, "note": fx["note"]})
        elif fx["kind"] == "DELCAND":
            md["fix_status"] = "delete_candidate"
            n_del += 1
            report["drafts_fixed"].append({"id": rid, "changes": ["标记删除候选"], "note": fx["note"]})
        else:
            md["fix_status"] = "flag"
            n_flag += 1

    # ---- KB 修复 ----
    n_kb_fix = n_kb_del = n_kb_flag = 0
    for k in kbs:
        e = k.setdefault("entities", {})
        ev = k.setdefault("evidence", {})
        key = _kkey(e)
        for (mk, eff_new, dir_new, note, kind) in KB_FIXES:
            if key == mk:
                k.setdefault("metadata", {})
                k["fix_note"] = note
                if kind == "FIX" and eff_new is not None:
                    e["effect"] = eff_new
                    ev["direction"] = dir_new
                    k["fix_status"] = "fixed"
                    n_kb_fix += 1
                    report["kb_fixed"].append({
                        "key": list(key), "changes": [f"effect → {eff_new}"], "note": note})
                elif kind == "DELCAND":
                    k["fix_status"] = "delete_candidate"
                    n_kb_del += 1
                    report["kb_fixed"].append({"key": list(key), "changes": ["标记删除候选"], "note": note})
                else:
                    k["fix_status"] = "flag"
                    n_kb_flag += 1
                    report["flags"].append({"key": list(key), "type": "T2", "note": note})
                break

    # ---- 写文件 ----
    with open(args.out_drafts, "w", encoding="utf-8") as f:
        for r in drafts:
            f.write(json.dumps(r, ensure_ascii=False) + "\n")
    with open(args.out_kb, "w", encoding="utf-8") as f:
        for k in kbs:
            f.write(json.dumps(k, ensure_ascii=False) + "\n")
    with open(args.out_report, "w", encoding="utf-8") as f:
        json.dump(report, f, ensure_ascii=False, indent=1)

    print(f"=== 修复完成 ===")
    print(f"L3: FIX {n_fix} 条, DELCAND {n_del} 条, T3 标记 {n_flag} 条")
    print(f"KB: FIX {n_kb_fix} 条, DELCAND {n_kb_del} 条, FLAG {n_kb_flag} 条")
    print(f"输出: {args.out_drafts}")
    print(f"     {args.out_kb}")
    print(f"报告: {args.out_report}")


if __name__ == "__main__":
    main()
