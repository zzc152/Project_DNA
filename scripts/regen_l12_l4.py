#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""regen_l12_l4.py — 局部重生成：L1/L2 用新直出逻辑（补核心实体），
L3/L5 保留旧 enhanced 结果，L4 交给 llm_enhancer 断点续跑（模板已更新）。
用法: python scripts/regen_l12_l4.py \
        --drafts data/synthetic/drafts_full2.jsonl \
        --old data/synthetic/enhanced_full.jsonl \
        --out data/synthetic/enhanced_full2.jsonl
"""

import argparse
import json
import sys

sys.path.insert(0, "src")
from data_synthesis.llm_enhancer import DIRECT_ANSWER_RULES  # noqa: E402


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--drafts", required=True)
    ap.add_argument("--old", required=True)
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    drafts = [json.loads(l) for l in open(args.drafts, encoding="utf-8")]
    old = {}
    for l in open(args.old, encoding="utf-8"):
        r = json.loads(l)
        old[r["id"]] = r

    n_new, n_keep, n_l4 = 0, 0, 0
    with open(args.out, "w", encoding="utf-8") as f:
        for d in drafts:
            tid = d["template_id"]
            if d["level"] == "L4":
                n_l4 += 1  # L4 模板已更新，一律留给 llm_enhancer 续跑
                continue
            if d["level"] in ("L1", "L2"):
                fn = DIRECT_ANSWER_RULES.get(tid)
                if fn is None:
                    raise ValueError(f"无直出规则: {tid}")
                rec = dict(d)
                rec["output"] = fn(d)
                rec["direct"] = True
                rec["metadata"] = dict(d["metadata"])
                rec["metadata"]["llm_model"] = "rule"
                n_new += 1
            elif d["id"] in old:
                rec = dict(old[d["id"]])  # L3/L5 保留旧结果（模板未变）
                rec.setdefault("direct", False)
                n_keep += 1
            else:
                raise ValueError(f"缺少旧结果且非直出/L4: {d['id']}")
            f.write(json.dumps(rec, ensure_ascii=False) + "\n")

    print(f"[L1/L2 新直出] {n_new} | [L3/L5 保留] {n_keep} | "
          f"[L4 待续跑] {n_l4} | 共 {len(drafts)}")


if __name__ == "__main__":
    main()
