#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""修复前后三层级评测对比（基线 scores vs final scores）"""
import json, sys, os
from collections import Counter, defaultdict

def load_scores(path):
    recs = []
    for l in open(path, encoding="utf-8"):
        l = l.strip()
        if not l:
            continue
        r = json.loads(l)
        recs.append(r)
    return recs

def report(recs, name=""):
    n = len(recs)
    if n == 0:
        return {"n": 0}
    avg = sum(r.get("total", 0) for r in recs) / n
    acc = sum(1 for r in recs if r.get("acc") == 60) / n
    pass80 = sum(1 for r in recs if r.get("total", 0) >= 80)
    claimed = sum(1 for r in recs if r.get("claimed"))
    confident_acc = sum(1 for r in recs if r.get("claimed") and r.get("acc") == 60) / max(1, claimed)
    overclaim = sum(1 for r in recs if r.get("overclaim"))
    ansE = sum(1 for r in recs if r.get("answer") == "E")
    return {
        "n": n, "avg": round(avg, 1), "acc": round(acc, 3), "pass80": pass80,
        "claimed": claimed, "confident_acc": round(confident_acc, 3),
        "overclaim": overclaim, "ansE_rate": round(ansE / n, 3),
    }

def by_group(recs, key):
    g = defaultdict(list)
    for r in recs:
        g[r.get(key, "?")].append(r)
    return g

def show(title, rows):
    print(f"\n### {title}")
    hdr = ["指标"] + [r[0] for r in rows]
    print("| " + " | ".join(hdr) + " |")
    print("|" + "---|" * len(hdr))
    keys = ["n", "avg", "acc", "pass80", "claimed", "confident_acc", "overclaim", "ansE_rate"]
    for k in keys:
        print("| " + k + " | " + " | ".join(str(r[1].get(k, "")) for r in rows) + " |")

def main():
    base_dir = "tmp"
    # 参数：--l3-32b-base ... 等，或默认路径
    pairs = {
        "L3-32B": (sys.argv[1] if len(sys.argv) > 1 else f"{base_dir}/l3_v21_clean_32b_scores.jsonl",
                   f"{base_dir}/l3_final_32b_scores.jsonl"),
        "L3-7B":  (sys.argv[2] if len(sys.argv) > 2 else f"{base_dir}/l3_v21_clean_7b_scores.jsonl",
                   f"{base_dir}/l3_final_7b_scores.jsonl"),
        "L12-32B":(sys.argv[3] if len(sys.argv) > 3 else f"{base_dir}/l12_full_32b_scores.jsonl",
                   f"{base_dir}/l12_final_32b_scores.jsonl"),
        "L12-7B": (sys.argv[4] if len(sys.argv) > 4 else f"{base_dir}/l12_full_7b_scores.jsonl",
                   f"{base_dir}/l12_final_7b_scores.jsonl"),
    }
    for tag, (base_p, final_p) in pairs.items():
        if not os.path.exists(base_p) or not os.path.exists(final_p):
            print(f"\n## {tag}: 缺文件 base={os.path.exists(base_p)} final={os.path.exists(final_p)}")
            continue
        base = load_scores(base_p)
        final = load_scores(final_p)
        print(f"\n## {tag}  基线: {os.path.basename(base_p)} ({len(base)} 条)  vs  final: {os.path.basename(final_p)} ({len(final)} 条)")
        show("总体", [("基线", report(base)), ("final", report(final))])
        # 分模板
        gb = by_group(base, "template_id")
        gf = by_group(final, "template_id")
        tpls = sorted(set(gb) | set(gf))
        rows = [(t, report(gf.get(t, []))) for t in tpls]
        show("分模板 (final)", rows)
        # 分 qtype
        gb2 = by_group(base, "qtype")
        gf2 = by_group(final, "qtype")
        qts = sorted(set(gb2) | set(gf2))
        rows = [(q, report(gf2.get(q, []))) for q in qts]
        show("分 qtype (final)", rows)

if __name__ == "__main__":
    main()
