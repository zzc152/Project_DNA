# -*- coding: utf-8 -*-
"""验证 l12_benchmark（L1 知识回忆 + L2 知识关联）的题面构造与评分逻辑。

用法：
  python scripts/verify_l12.py
"""
import json
import random
import sys
from collections import Counter

sys.path.insert(0, "src")
from data_synthesis.l12_benchmark import (
    build_all, score_l12, _clean_factor, _kb_effect,
)
from data_synthesis.l3_benchmark_v2 import _balance_directions

recs = [json.loads(l) for l in open("data/processed/knowledge_base_clean.jsonl", encoding="utf-8")]

print("=" * 70)
print("1. 题面构造（seed=42 确定性）")
qs = build_all(recs, random.Random(42), 60, 60)
c = Counter((q["qtype"], q["gold"]) for q in qs.values())
print(f"  总计 {len(qs)} 条")
for k, v in sorted(c.items()):
    print(f"    {k[0]} gold={k[1]}: {v}")

# 确定性检查
qs2 = build_all(recs, random.Random(42), 60, 60)
assert list(qs.keys()) == list(qs2.keys()), "seed=42 构造不稳定!"
print("  ✓ seed=42 确定性通过")

# 结构检查
for qid, q in qs.items():
    assert "question" in q and "gold" in q and "opts" in q, f"{qid} 缺字段"
    assert q["gold"] in ("A", "B", "C", "D", "E"), f"{qid} gold 非法: {q['gold']}"
    assert "结论" not in q["question"] and "摘要" not in q["question"], \
        f"{qid} 泄露了 claim/摘要信息!"
    if q["qtype"] == "L2":
        assert q["cell"], f"{qid} L2 缺细胞系约束"
    else:
        assert "在" not in q["question"].split("问题：")[1].split("细胞中")[0], \
            f"{qid} L1 不应带细胞系"
print("  ✓ 题面结构完整 + 无 claim/摘要泄露 + 信息量递减(L1<L2)")

print("=" * 70)
print("2. 平衡性（A/B 严格各半，防蒙选项）")
a = sum(1 for q in qs.values() if q["gold"] == "A")
b = sum(1 for q in qs.values() if q["gold"] == "B")
print(f"  A={a} B={b} 差={abs(a-b)}")
assert a == b, "A/B 未严格平衡!"
print("  ✓ A/B 平衡通过")

print("=" * 70)
print("3. 因子名过滤")
noise = [r for r in recs if _kb_effect(r) and _clean_factor(r) is None]
print(f"  被过滤的噪声记录数: {len(noise)}")
for r in noise[:3]:
    print(f"    - {_kb_effect(r)} | 因子: {r.get('entities', {}).get('factor')}")
assert all(_clean_factor(r) for r in recs if _kb_effect(r) and r["entities"].get("factor") and len(r["entities"]["factor"]) <= 30 and r["entities"]["factor"].count(" ") < 4)
print("  ✓ 过滤规则正确")

print("=" * 70)
print("4. 评分逻辑")
cases = [
    # (q, output, 期望 total)
    ({"id": "T1", "qtype": "L1", "question": "问题：因子X对enhancer活性的影响方向？\n选项 A增强 B减弱 C无影响 D方向不明 E无法判断", "gold": "A", "opts": "A增强 B减弱 C无影响 D方向不明 E无法判断", "factor": "MYCN", "element": "enhancer", "cell": None, "effect": "increases"},
     "A。MYCN增强enhancer活性", 60 + 40),
    ({"id": "T2", "qtype": "L1", "question": "问题：因子X对enhancer活性的影响方向？", "gold": "A", "opts": "A增强 B减弱 C无影响 D方向不明 E无法判断", "factor": "MYCN", "element": "enhancer", "cell": None, "effect": "increases"},
     "B。MYCN减弱", 0 + 40),
    ({"id": "T3", "qtype": "L2", "question": "问题：在K562细胞中，因子X对enhancer活性的影响方向？", "gold": "B", "opts": "A增强 B减弱 C无影响 D方向不明 E无法判断", "factor": "BCL11A", "element": "enhancer", "cell": "K562", "effect": "decreases"},
     "E。信息不足无法判断", 30 + 20),
    ({"id": "T4", "qtype": "L1", "question": "问题：因子X对enhancer活性的影响方向？", "gold": "A", "opts": "A增强 B减弱 C无影响 D方向不明 E无法判断", "factor": "GATA1", "element": "enhancer", "cell": None, "effect": "increases"},
     "", 0 + 20),
]
for q, out, expect in cases:
    s = score_l12(q["id"], q, out)
    status = "✓" if s["total"] == expect else "✗"
    print(f"  {status} {q['id']} [{q['qtype']}] gold={q['gold']} ans={s['answer']} "
          f"total={s['total']} (期望 {expect}) claimed={s['claimed']}")
    assert s["total"] == expect, f"{q['id']} 评分不符: {s['total']} != {expect}"

print("=" * 70)
print("5. 诚实性字段")
s_over = score_l12("T5", cases[1][0], "B。MYCN减弱")
assert s_over["claimed"] and s_over["overclaim"], "T5 应标记 overclaim"
print(f"  ✓ 声称知道但答错 → claimed={s_over['claimed']} overclaim={s_over['overclaim']}")

print("=" * 70)
print("全部通过 ✓")
