#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Gold 人工核查落地脚本（review_gold_manual.py）

在 fix_gold.py 的基础上，落地人工核查裁决（docs/gold_manual_review.md）：

  1. MYCN FLAG 11 条 — 逐条回原文核查：
     - 9 条方向正确，保留（fix_status=resolved，补 target_gene 上下文）
     - 2 条（DLX5 required_for×2，PMID 41135749）证据不足降级为 modulates
       —— 摘要原文只说 "MYCN coregulates DLX5 ... by directly binding to the DLX5
          promoter"，"coregulates" 未指明激活/抑制方向，required_for 是过度推断。
  2. DELCAND 6 条 — 人工确认删除：
     - L3 3 条：L3_002_26 / L3_003_29（因子"epigenetic status"不可考+同 PMID 重复）
               L3_002_2（BCL11A 题目语义反向）
     - KB 3 条：DNA methylation→gene upregulation（co_occurs_with 共现非因果）
               chromatin accessibility→promoter methylation level（correlates 相关非因果）
               targeted activation...→AML differentiation（promotes 描述性因子不可考）

输出：
  - data/synthetic/drafts_full2_clean_final.jsonl（删除 DELCAND 后）
  - data/processed/knowledge_base_clean_final.jsonl（删除 DELCAND + MYCN 裁决后）
  - tmp/manual_review.json（逐条裁决明细）
"""
import argparse
import json

# ---------------------------------------------------------------------------
# MYCN 逐条裁决表：按 (factor, element, effect, pmids) 定位
# 裁决说明见 docs/gold_manual_review.md
# ---------------------------------------------------------------------------

MYCN_REVIEW = [
    # --- MYCN→promoter，DLX5，required_for → modulates（原文"coregulates"方向未明）---
    {"match": {"effect": "required_for", "pmids": ["41135749"], "target": "DLX5"},
     "action": "down_grade", "new_effect": "modulates", "new_direction": "modulates_activity",
     "target_gene": "DLX5",
     "note": "人工核查：摘要原文仅证明MYCN直接结合DLX5启动子并coregulates表达，coregulates未指明激活/抑制方向，required_for属过度推断→降级为modulates（与L3_002_37 ASCL1处理一致）",
     "verdict": "证据不足，降级"},
    # --- MYCN→enhancer，modulates（PMID 36598365，仅"binds"）---
    {"match": {"effect": "modulates", "pmids": ["36598365"]},
     "action": "keep", "target_gene": "neuroblastoma targets",
     "note": "人工核查：摘要原文'MYCN binds to enhancers when aberrantly upregulated'，仅证明结合未指明方向，modulates保守正确（保留）",
     "verdict": "方向正确，保留"},
    # --- MYCN→promoter，increases（PMID 37781575/36539767，MYC targets）---
    {"match": {"effect": "increases", "regulatory_element": "promoter", "pmids": ["37781575", "36539767"]},
     "action": "keep", "target_gene": "MYC targets",
     "note": "人工核查：摘要原文'MYCN-binding increases H3K4me3 and H3K27ac on canonical MYC target promoters'，激活标记增加→increases正确（保留）",
     "verdict": "方向正确，保留"},
    # --- MYCN→enhancer，decreases（PMID 37781575，neuronal differentiation）---
    {"match": {"effect": "decreases", "regulatory_element": "enhancer", "pmids": ["37781575"]},
     "action": "keep", "target_gene": "neuronal differentiation genes",
     "note": "人工核查：摘要原文'MYCN recruits G9a to enhancers to repress neuronal differentiation genes'，G9a为抑制性修饰→decreases正确（保留）",
     "verdict": "方向正确，保留"},
    # --- MYCN→promoter，decreases（PMID 37781575，neuronal differentiation）---
    {"match": {"effect": "decreases", "regulatory_element": "promoter", "pmids": ["37781575"]},
     "action": "keep", "target_gene": "neuronal differentiation genes",
     "note": "人工核查：同摘要后半句'decreases H3K27ac on neuronal differentiation gene enhancers and promoters'，明确覆盖promoter→decreases正确（保留）",
     "verdict": "方向正确，保留"},
    # --- MYCN→enhancer，increases（PMID 36539767，SMAD9）---
    {"match": {"effect": "increases", "regulatory_element": "enhancer", "pmids": ["36539767"]},
     "action": "keep", "target_gene": "SMAD9",
     "note": "人工核查：摘要原文'SMAD9 specifically induced by TFs including MYCN at the enhancer region'，诱导高表达→increases正确（保留）",
     "verdict": "方向正确，保留"},
]

# 待删除的 DELCAND（人工确认）
DELCAND_DRAFTS_IDS = {"L3_002_26", "L3_003_29", "L3_002_2"}
DELCAND_KB_KEYS = {
    ("DNA methylation", "gene upregulation", "not_specified"),
    ("chromatin accessibility", "promoter methylation level", "not_specified"),
    ("targeted activation of this regulatory element", "AML differentiation", "not_specified"),
}


def _kkey(e):
    return (str(e.get("factor") or "").strip(),
            str(e.get("regulatory_element") or "").strip(),
            str(e.get("cell_line") or "").strip())


def _match(rec, spec):
    """按 spec 中的字段匹配 KB 记录（未指定的字段不参与匹配）。"""
    e = rec.get("entities", {})
    if "effect" in spec and e.get("effect") != spec["effect"]:
        return False
    if "regulatory_element" in spec and e.get("regulatory_element") != spec["regulatory_element"]:
        return False
    if "pmids" in spec and e.get("pmids") != spec["pmids"]:
        return False
    if "target" in spec and spec["target"] not in (e.get("context") or ""):
        return False
    return True


def main():
    ap = argparse.ArgumentParser(description="Gold 人工核查落地")
    ap.add_argument("--drafts", default="data/synthetic/drafts_full2_clean_fixed.jsonl")
    ap.add_argument("--kb", default="data/processed/knowledge_base_clean_fixed.jsonl")
    ap.add_argument("--out-drafts", default="data/synthetic/drafts_full2_clean_final.jsonl")
    ap.add_argument("--out-kb", default="data/processed/knowledge_base_clean_final.jsonl")
    ap.add_argument("--out-report", default="tmp/manual_review.json")
    args = ap.parse_args()

    report = {"mycn_review": [], "delcand_deleted": []}

    # ---- 1. drafts：删除 DELCAND ----
    drafts = [json.loads(l) for l in open(args.drafts, encoding="utf-8") if l.strip()]
    n_before = len(drafts)
    kept, deleted = [], []
    for r in drafts:
        rid = r.get("id", "")
        if rid in DELCAND_DRAFTS_IDS:
            deleted.append({"id": rid, "factor": r.get("metadata", {}).get("entities_used", {}).get("factor"),
                            "note": r.get("metadata", {}).get("fix_note")})
        else:
            kept.append(r)
    report["delcand_deleted"].extend(deleted)

    # ---- 2. KB：删除 DELCAND + MYCN 裁决 ----
    kbs = [json.loads(l) for l in open(args.kb, encoding="utf-8") if l.strip()]
    n_kb = len(kbs)
    kb_kept, kb_del, n_mycn = [], [], 0
    for k in kbs:
        e = k.get("entities", {})
        if k.get("fix_status") == "delete_candidate" and _kkey(e) in DELCAND_KB_KEYS:
            kb_del.append({"key": list(_kkey(e)), "note": k.get("fix_note")})
            continue
        # MYCN 逐条裁决
        if e.get("factor") == "MYCN":
            for spec in MYCN_REVIEW:
                if _match(k, spec["match"]):
                    n_mycn += 1
                    k["fix_status"] = "resolved"
                    k.setdefault("metadata", {})
                    k["metadata"]["target_gene"] = spec.get("target_gene", "")
                    if spec["action"] == "down_grade":
                        e["effect"] = spec["new_effect"]
                        k["evidence"]["direction"] = spec["new_direction"]
                    k["fix_note"] = spec["note"]
                    report["mycn_review"].append({
                        "key": list(_kkey(e)), "action": spec["action"],
                        "target_gene": spec.get("target_gene", ""),
                        "verdict": spec["verdict"], "note": spec["note"]})
                    break
        kb_kept.append(k)

    # ---- 3. 写文件 ----
    with open(args.out_drafts, "w", encoding="utf-8") as f:
        for r in kept:
            f.write(json.dumps(r, ensure_ascii=False) + "\n")
    with open(args.out_kb, "w", encoding="utf-8") as f:
        for k in kb_kept:
            f.write(json.dumps(k, ensure_ascii=False) + "\n")
    with open(args.out_report, "w", encoding="utf-8") as f:
        json.dump(report, f, ensure_ascii=False, indent=2)

    # ---- 4. 摘要 ----
    print(f"drafts: {n_before} → {len(kept)}（删除 {len(deleted)} 条 DELCAND）")
    print(f"  L3 DELCAND 删除: {[d['id'] for d in deleted]}")
    print(f"kb: {n_kb} → {len(kb_kept)}（删除 {len(kb_del)} 条 DELCAND）")
    print(f"  KB DELCAND 删除: {[tuple(d['key']) for d in kb_del]}")
    print(f"MYCN 裁决: {n_mycn} 条（9 保留 + 2 降级）")
    for r in report["mycn_review"]:
        print(f"  {tuple(r['key'])} → {r['action']}（{r['verdict']}）")
    print(f"报告: {args.out_report}")


if __name__ == "__main__":
    main()
