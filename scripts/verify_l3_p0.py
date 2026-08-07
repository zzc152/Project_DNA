# -*- coding: utf-8 -*-
"""验证 l3_benchmark_p0 新题型（矛盾检测 + evidence-level）的题面构造与评分逻辑。

用法：
  python scripts/verify_l3_p0.py
"""
import json
import random
import sys
from collections import Counter

sys.path.insert(0, "src")
from data_synthesis.l3_benchmark_p0 import (
    build_all, score_p0, _extract_evidence_sent, _evid_type,
)

recs = [json.loads(l) for l in open("data/processed/knowledge_base_clean.jsonl", encoding="utf-8")]

print("=" * 70)
print("1. 题面构造（seed=42 确定性）")
qs = build_all(recs, random.Random(42), 60, 60)
c = Counter((q["qtype"], q["gold"]) for q in qs.values())
print(f"  总计 {len(qs)} 条")
for k, v in sorted(c.items()):
    print(f"    {k[0]} gold={k[1]}: {v}")

# 确定性检查：同 seed 两次构造一致
qs2 = build_all(recs, random.Random(42), 60, 60)
assert [q["id"] for q in qs.values()] == [q["id"] for q in qs2.values()], "seed=42 构造不稳定!"
print("  ✓ seed=42 确定性通过")

# 结构检查
for qid, q in qs.items():
    assert "question" in q and "gold" in q and "opts" in q, f"{qid} 缺字段"
    assert q["gold"] in ("A", "B", "C"), f"{qid} gold 非法: {q['gold']}"
    if q["qtype"] == "contradiction":
        assert "结论 1" in q["question"] and "结论 2" in q["question"]
        assert q["claim_a"] and q["claim_b"]
    if q["qtype"] == "evidence":
        assert "以下是一条文献证据" in q["question"]
        assert q["evidence"]
print("  ✓ 题面结构完整")

# evidence 分类健全性：能推出(A) 必须来自活性证据，不能推出(B) 必须来自表达证据
for qid, q in qs.items():
    if q["qtype"] != "evidence":
        continue
    if q["gold"] == "A":
        assert _evid_type(q["evidence"]) == "act", f"{qid} gold=A 但证据非活性类"
    elif q["gold"] == "B":
        assert _evid_type(q["evidence"]) == "expr", f"{qid} gold=B 但证据非表达类"
print("  ✓ evidence 证据-结论层级匹配（表达≠活性）")

print("=" * 70)
print("2. 评分逻辑")
# 答对 / 答错但声称 / 诚实答C / 空答 四类
cases = [
    # (q, output, 期望)
    ({"id": "T1", "qtype": "contradiction", "question": "结论 1: MYCN increases promoter\n结论 2: MYCN decreases promoter\n选项 A 相互矛盾 B 支持 C 无关", "gold": "A", "opts": "A. 相互矛盾\nB. 相互支持\nC. 无关", "claim_a": "x", "claim_b": "y"},
     "A。MYCN对promoter活性方向相反", 60 + 40),
    ({"id": "T2", "qtype": "contradiction", "question": "结论 1: X increases\n结论 2: X decreases\n选项 A 相互矛盾 B 支持 C 无关", "gold": "A", "opts": "A. 相互矛盾\nB. 相互支持\nC. 无关", "claim_a": "x", "claim_b": "y"},
     "B。两条结论一致", 0 + 20),  # 答错(0) + 无编造(20)，理由未引用实体
    ({"id": "T3", "qtype": "evidence", "question": "证据: RNA-seq 显示表达上调\n结论: X 调节 promoter 活性\n选项 A 能推出 B 不能推出 C 无法判断", "gold": "B", "opts": "A. 能推出\nB. 不能推出\nC. 无法判断", "evidence": "RNA-seq 显示表达上调"},
     "C。信息不足无法判断", 30 + 20),  # 诚实分(30) + 无编造(20)，理由未引用实体
    ({"id": "T4", "qtype": "evidence", "question": "证据: reporter 活性\n结论: X 调节 enhancer 活性\n选项 A 能推出 B 不能推出 C 无法判断", "gold": "A", "opts": "A. 能推出\nB. 不能推出\nC. 无法判断", "evidence": "reporter 活性"},
     "", 0 + 20),  # 空答(0) + 无编造(20)
]
for q, out, expect in cases:
    s = score_p0({"id": q["id"]}, q, out)
    status = "✓" if s["total"] == expect else "✗"
    print(f"  {status} {q['id']} [{q['qtype']}] gold={q['gold']} ans={s['answer']} "
          f"total={s['total']} (期望 {expect}) claimed={s['claimed']}")
    assert s["total"] == expect, f"{q['id']} 评分不符: {s['total']} != {expect}"

print("=" * 70)
print("3. 诚实性字段")
s_over = score_p0({"id": "T5"}, cases[1][0], "B。两条结论一致")
assert s_over["claimed"] and s_over["overclaim"], "T5 应标记 overclaim"
print(f"  ✓ 声称知道但答错 → claimed={s_over['claimed']} overclaim={s_over['overclaim']}")

print("=" * 70)
print("全部通过 ✓")
