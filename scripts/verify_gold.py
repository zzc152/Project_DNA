#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Gold 标注质量验证脚本（verify_gold.py）

对数据集做系统性的 gold 验证，回答"标注的答案是否正确"：

  T1  L3 证据字段内部一致性
      metadata.entities_used.effect  vs  metadata.evidence.direction
      归一化后必须一致；不一致 = 标注矛盾（同一记录两个字段说不同方向）。
  T2  KB（L1/L2 数据源）实体对方向冲突
      同一 (factor, element) 对出现 increases 与 decreases 共存 = 数据冲突。
  T3  L3 题面因子名合法性
      因子名是否是可考实体（过滤描述性短语，如 "epigenetic status"）。
  T4  L1/L2 题面 vs KB 依据一致性（gold 可复现）
      从 KB effect 重新推导的 gold 必须 == 题目实际 gold。
  T5  摘要↔gold 语义一致性（--model 提供时启用，rc 型）
      仅给摘要原文句，让模型独立判断方向，与数据标注 gold 对比。
      一致率高 → gold 可从证据推出（数据质量好）；
      不一致 → 候选问题条目（gold 错 或 摘要信息不足）。

输出：
  - 控制台逐项 PASS/FAIL + 汇总
  - --out 指定的 JSON 报告
  - --audit 指定的人工抽检清单（markdown）
"""
import argparse
import json
import re
import random
import sys
import os
from collections import defaultdict

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(ROOT, "src"))
try:
    from data_synthesis.l3_benchmark_v2 import (
        _gold_direction, _norm_ent, _abs_sentence, _factor_name, _element,
        _cell_line, load_model, generate,
    )
except Exception:
    from l3_benchmark_v2 import (
        _gold_direction, _norm_ent, _abs_sentence, _factor_name, _element,
        _cell_line, load_model, generate,
    )


# ---------------------------------------------------------------------------
# 归一化
# ---------------------------------------------------------------------------

def norm_direction(d: str):
    """把方向字段归一化到 {increases, decreases, modulates, required_for, no_effect, None}。"""
    d = (d or "").lower()
    if "increas" in d:
        return "increases"
    if "decreas" in d:
        return "decreases"
    if "modulat" in d or "bind" in d:
        return "modulates"
    if "required" in d:
        return "required_for"
    if "no_effect" in d or "noeffect" in d:
        return "no_effect"
    return d or None


def norm_ent(s: str) -> str:
    return _norm_ent(s) or ""


# ---------------------------------------------------------------------------
# T1: L3 证据字段内部一致性
# ---------------------------------------------------------------------------

def t1_l3_evidence_consistency(l3_recs) -> dict:
    issues = []
    no_effect = []
    for r in l3_recs:
        md = r.get("metadata", {})
        eu = md.get("entities_used", {})
        ev = md.get("evidence", {})
        eff = eu.get("effect")
        dirv = ev.get("direction")
        if not dirv:
            continue  # 无 evidence 的跳过（L3_003 数据型）
        if not eff:
            no_effect.append(r.get("id"))
            continue
        n_eff, n_dir = norm_direction(eff), norm_direction(dirv)
        if n_eff != n_dir:
            issues.append({
                "id": r.get("id"),
                "factor": eu.get("factor"),
                "element": eu.get("regulatory_element"),
                "effect": eff,
                "evidence_direction": dirv,
                "gold_by_effect": _gold_direction(eff),
                "gold_by_evidence": _gold_direction(dirv),
            })
    return {"issues": issues, "no_effect": no_effect}


# ---------------------------------------------------------------------------
# T2: KB 实体对方向冲突
# ---------------------------------------------------------------------------

def t2_kb_conflict(kb_recs) -> dict:
    """同一 (factor, element) 对方向冲突。

    升级：加入 target_gene 维度（若记录有 metadata.target_gene），
    因为同因子对不同靶基因可以有真实相反的效应（如 MYCN 激活 MYC
    靶基因、抑制神经元分化基因），此时不算标注冲突。
    无 target_gene 的记录退回原 (factor, element) 分组。
    """
    pairs = defaultdict(list)
    for k in kb_recs:
        e = k.get("entities", {})
        key0 = (norm_ent(e.get("factor")), norm_ent(e.get("regulatory_element")))
        if not key0[0] or not key0[1]:
            continue
        tg = (k.get("metadata", {}) or {}).get("target_gene") or ""
        key = (key0[0], key0[1], tg)
        pairs[key].append((e.get("effect"), e.get("cell_line")))
    conflicts = []
    for key, vals in pairs.items():
        effs = [v for v, _ in vals if v]
        has_inc = any("increas" in str(v) or "required" in str(v) for v in effs)
        has_dec = any("decreas" in str(v) for v in effs)
        if has_inc and has_dec:
            conflicts.append({"pair": list(key), "records": vals[:6]})
    return {"conflicts": conflicts, "n_pairs": len(pairs)}


# ---------------------------------------------------------------------------
# T3: L3 因子名合法性
# ---------------------------------------------------------------------------

_DESC_PATTERNS = re.compile(
    r"^(epigenetic status|chromatin|transcription factor|gene expression|"
    r"regulatory element|binding|recruitment|activity|expression|"
    r"methylation|acetylation)\b", re.I)


def t3_factor_name(l3_recs) -> dict:
    issues = []
    for r in l3_recs:
        md = r.get("metadata", {})
        eu = md.get("entities_used", {})
        f = eu.get("factor") or ""
        # 描述性短语或超长
        if (len(f) > 30 or f.count(" ") >= 4 or _DESC_PATTERNS.match(f.strip())) and f:
            issues.append({"id": r.get("id"), "factor": f})
    return {"issues": issues}


# ---------------------------------------------------------------------------
# T4: L1/L2 题面 gold 可复现（依赖 l12_benchmark 的构造）
# ---------------------------------------------------------------------------

def t4_l12_gold_reproducible(kb_recs, l12_scores_fn=None) -> dict:
    """从 KB 重新推导 L1/L2 的 gold 并与题目实际 gold 对比。

    L1/L2 题面是运行时从 KB 构造的，这里无法拿到构造后的逐题 gold，
    只能验证 KB 中每条有效记录能推导出非 E 的 gold，且 effect 合法。
    """
    issues = []
    for k in kb_recs:
        e = k.get("entities", {})
        eff = e.get("effect")
        if not eff or str(eff).lower() in ("none", "null"):
            continue
        g = _gold_direction(eff)
        if g == "E":
            issues.append({
                "factor": e.get("factor"),
                "element": e.get("regulatory_element"),
                "effect": eff,
                "gold": g,
                "msg": "effect 无法映射到 A-D（异常值）",
            })
    return {"issues": issues}


# ---------------------------------------------------------------------------
# T5: 摘要↔gold 语义一致性（模型）
# ---------------------------------------------------------------------------

_T5_PROMPT = (
    "以下是某篇生物文献摘要中的一句话：\n“{abs}”\n\n"
    "问题：根据这句话，因子“{factor}”对“{element}”活性的作用是？\n"
    "只回答一个字母：\n"
    "A. 增强/促进（含必需）\n"
    "B. 减弱/抑制\n"
    "C. 无影响\n"
    "D. 有调节作用但方向不明\n"
    "E. 无法判断（这句话信息不足）\n"
    "请只输出选项字母，不要输出其他内容。"
)

_OPT_RE = re.compile(r"(?<![A-Za-z0-9])([A-E])(?![A-Za-z0-9])")


def t5_semantic_consistency(l3_recs, model, tok, device, limit=63) -> dict:
    rc_recs = [r for r in l3_recs
               if r.get("template_id") in ("L3_001", "L3_002")]
    rng = random.Random(42)
    rng.shuffle(rc_recs)
    rc_recs = rc_recs[:limit]
    results = []
    for r in rc_recs:
        md = r.get("metadata", {})
        eu = md.get("entities_used", {})
        factor = _factor_name(r) or eu.get("factor") or ""
        elem = _element(r) or eu.get("regulatory_element") or ""
        abs_s = _abs_sentence(r)
        if not abs_s or not factor:
            continue
        prompt = _T5_PROMPT.format(abs=abs_s, factor=factor, element=elem)
        out = generate(model, tok, prompt, 16)
        m = _OPT_RE.search(out)
        pred = m.group(1) if m else "?"
        gold = _gold_direction(eu.get("effect"))
        results.append({
            "id": r.get("id"),
            "factor": factor,
            "element": elem,
            "abs": abs_s[:200],
            "gold": gold,
            "pred": pred,
            "match": pred == gold,
        })
    return {"results": results}


# ---------------------------------------------------------------------------
# 人工抽检清单
# ---------------------------------------------------------------------------

def build_audit(l3_recs, kb_recs, n=30, seed=42) -> list:
    rng = random.Random(seed)
    # 混合：L3 方向题（rc/know）+ 异常条目（矛盾/缺 effect）
    rc_recs = [r for r in l3_recs if r.get("template_id") in ("L3_001", "L3_002")]
    rng.shuffle(rc_recs)
    audit = []
    for r in rc_recs[:n]:
        md = r.get("metadata", {})
        eu = md.get("entities_used", {})
        ev = md.get("evidence", {})
        abs_s = _abs_sentence(r)
        audit.append({
            "id": r.get("id"),
            "factor": eu.get("factor"),
            "element": eu.get("regulatory_element"),
            "cell": eu.get("cell_line"),
            "abstract": abs_s,
            "effect": eu.get("effect"),
            "evidence_direction": ev.get("direction"),
            "gold": _gold_direction(eu.get("effect")),
            "source": "L3",
        })
    return audit


def audit_markdown(audit: list) -> str:
    lines = ["| id | 因子 | 元件 | 细胞 | 摘要原文 | 标注 effect | evidence.direction | gold |",
             "|----|------|------|------|----------|-------------|--------------------|------|"]
    for a in audit:
        lines.append(
            f"| {a['id']} | {a['factor'] or ''} | {a['element'] or ''} | {a['cell'] or ''} | "
            f"{str(a['abstract'])[:60]} | {a['effect'] or ''} | {a['evidence_direction'] or ''} | {a['gold']} |")
    return "\n".join(lines)


# ---------------------------------------------------------------------------
# main
# ---------------------------------------------------------------------------

def main():
    ap = argparse.ArgumentParser(description="Gold 标注质量验证")
    ap.add_argument("--drafts", default="data/synthetic/drafts_full2_clean.jsonl")
    ap.add_argument("--kb", default="data/processed/knowledge_base_clean.jsonl")
    ap.add_argument("--model", default=None, help="模型路径（提供则启用 T5 语义一致性）")
    ap.add_argument("--device", default="cuda:0")
    ap.add_argument("--t5-limit", type=int, default=63)
    ap.add_argument("--out", default=None, help="JSON 报告输出路径")
    ap.add_argument("--audit", default=None, help="人工抽检清单输出路径")
    ap.add_argument("--seed", type=int, default=42)
    args = ap.parse_args()

    drafts = [json.loads(l) for l in open(args.drafts, encoding="utf-8") if l.strip()]
    kbs = [json.loads(l) for l in open(args.kb, encoding="utf-8") if l.strip()]
    l3 = [d for d in drafts if d.get("level") == "L3"]
    print(f"数据规模: drafts={len(drafts)} L3={len(l3)} KB={len(kbs)}")

    report = {}

    # T1
    print("\n=== T1 L3 证据字段内部一致性 ===")
    t1 = t1_l3_evidence_consistency(l3)
    print(f"  矛盾条数: {len(t1['issues'])} / {len(l3)}")
    for i in t1["issues"]:
        print(f"    [{i['id']}] effect={i['effect']}(gold {i['gold_by_effect']}) vs "
              f"evidence.direction={i['evidence_direction']}(gold {i['gold_by_evidence']})")
    if t1["no_effect"]:
        print(f"  缺 effect 但带 evidence.direction: {len(t1['no_effect'])} 条 {t1['no_effect'][:8]}")
    report["T1"] = t1

    # T2
    print("\n=== T2 KB 实体对方向冲突 ===")
    t2 = t2_kb_conflict(kbs)
    print(f"  方向冲突对: {len(t2['conflicts'])} / {t2['n_pairs']} 有 effect 的对")
    for c in t2["conflicts"]:
        print(f"    {c['pair']}: {c['records'][:3]}")
    report["T2"] = t2

    # T3
    print("\n=== T3 L3 因子名合法性 ===")
    t3 = t3_factor_name(l3)
    print(f"  可疑因子名: {len(t3['issues'])}")
    for i in t3["issues"]:
        print(f"    [{i['id']}] {i['factor']!r}")
    report["T3"] = t3

    # T4
    print("\n=== T4 KB effect→gold 可映射 ===")
    t4 = t4_l12_gold_reproducible(kbs)
    print(f"  无法映射的 effect: {len(t4['issues'])}")
    for i in t4["issues"][:10]:
        print(f"    {i['factor']}/{i['element']}: {i['effect']!r}")
    report["T4"] = t4

    # T5
    if args.model:
        print(f"\n=== T5 摘要↔gold 语义一致性 (model={args.model}) ===")
        print("  加载模型...")
        model, tok = load_model(args.model, args.device)
        t5 = t5_semantic_consistency(l3, model, tok, args.device, limit=args.t5_limit)
        res = t5["results"]
        n = len(res)
        match = sum(1 for x in res if x["match"])
        print(f"  样本数: {n}, 一致: {match} ({match / n:.1%} 若 n>0)")
        mism = [x for x in res if not x["match"]]
        for x in mism[:15]:
            print(f"    [{x['id']}] gold={x['gold']} pred={x['pred']} factor={x['factor']} "
                  f"elem={x['element']}\n      摘要: {x['abs'][:100]}")
        print(f"  (不一致 {len(mism)} 条)")
        report["T5"] = t5
    else:
        print("\n=== T5 摘要↔gold 语义一致性（跳过：未提供 --model）===")

    # 汇总
    print("\n================ 汇总 ================")
    print(f"T1 证据矛盾: {len(t1['issues'])} 条")
    print(f"T2 KB 方向冲突: {len(t2['conflicts'])} 对")
    print(f"T3 可疑因子名: {len(t3['issues'])} 条")
    print(f"T4 不可映射 effect: {len(t4['issues'])} 条")
    if args.model and "T5" in report:
        res = report["T5"]["results"]
        n = len(res)
        match = sum(1 for x in res if x["match"])
        print(f"T5 摘要↔gold 一致率: {match}/{n} ({match / n:.1%})")

    if args.out:
        with open(args.out, "w", encoding="utf-8") as f:
            json.dump(report, f, ensure_ascii=False, indent=1)
        print(f"\n报告已写入: {args.out}")

    if args.audit:
        audit = build_audit(l3, kbs, n=30, seed=args.seed)
        with open(args.audit, "w", encoding="utf-8") as f:
            f.write(audit_markdown(audit))
        print(f"抽检清单已写入: {args.audit} ({len(audit)} 条)")


if __name__ == "__main__":
    main()
