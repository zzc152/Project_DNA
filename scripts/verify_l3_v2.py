# -*- coding: utf-8 -*-
"""验证 l3_benchmark_v2 题面构造与解析。"""
import json
import random
import sys
sys.path.insert(0, "src")
from data_synthesis.l3_benchmark_v2 import (
    build_questions, parse_answer, _clean_data_line, score_v2, _balance_directions,
)

recs = [json.loads(l) for l in open("data/synthetic/drafts_full2.jsonl", encoding="utf-8")
        if json.loads(l).get("level") == "L3"]
rng = random.Random(42)
qs = {d["id"]: build_questions(d, i, rng) for i, d in enumerate(recs)}

# 1. 清洗效果
print("=== 清洗效果 ===")
print(repr(_clean_data_line("Spi1仅在K562富集（OR=1.93），其它细胞系未达阈值")))
print(repr(_clean_data_line("富集基序6个，最强OR=10.41")))

# 2. 类型/金标分布
from collections import Counter
print("\n=== qtype 分布 ===", dict(Counter(q["qtype"] for q in qs.values())))
print("=== gold 分布 ===", dict(Counter(q["gold"] for q in qs.values())))
by_tpl = {}
for d in recs:
    by_tpl.setdefault(d["template_id"], []).append(qs[d["id"]])
for tid, qlist in by_tpl.items():
    print(f"  {tid}: n={len(qlist)} gold={dict(Counter(q['gold'] for q in qlist))} "
          f"qtype={dict(Counter(q['qtype'] for q in qlist))}")

# 3. 解析器测试
print("\n=== parse_answer 测试 ===")
tests = [
    ("A。摘要提到该因子是正向调节因子", "A"),
    ("答案是 B", "B"),
    ("C、无影响", "C"),
    ("我选择D，因为调节方向未明确", "D"),
    ("E 无法判断", "E"),
    ("该因子增强启动子活性", "A"),
    ("结论不成立，摘要未支持", "B"),
    ("无法从题面判断", "E"),
    ("成立", "A"),
]
for out, exp in tests:
    got = parse_answer(out)
    flag = "OK" if got == exp else f"FAIL(exp={exp})"
    print(f"  [{flag}] {out!r} -> {got!r}")

# 4. 评分逻辑测试
print("\n=== score_v2 测试 ===")
d0 = recs[0]  # L3_001_1, rc, gold=A
q0 = qs[d0["id"]]
print("题面类型:", q0["qtype"], "| gold:", q0["gold"])
cases = [
    ("A。SP1是promoter的正向调节因子", 60, 40),
    ("B。SP1抑制promoter", 0, 40),
    ("E。题面信息不足", 0, 20),  # 未引用实体 → 仅"无编造"20 分
]
for out, exp_acc, exp_fid in cases:
    s = score_v2(d0, q0, out)
    ok = "OK" if (s["acc"], s["fidelity"]) == (exp_acc, exp_fid) else f"FAIL exp={(exp_acc,exp_fid)}"
    print(f"  [{ok}] {out!r} -> acc={s['acc']} fid={s['fidelity']} total={s['total']}")

# 5. know 型诚实分测试（找一个 know 型）
dk = next(d for d in recs if qs[d["id"]]["qtype"] == "know")
qk = qs[dk["id"]]
s = score_v2(dk, qk, "E。题面没有给文献信息，无法判断")
print(f"  know型诚实分: gold={qk['gold']} ans=E -> acc={s['acc']} (期望30)")

# 6. 正负方向平衡测试（v2.1）
print("\n=== _balance_directions 平衡测试 ===")
qb = _balance_directions(qs, rng)
gold_b = Counter(q["gold"] for q in qb.values())
print(f"  平衡后 gold 分布: {dict(gold_b)} | 条数 {len(qb)}")
ok = gold_b["A"] == gold_b["B"] and len(qb) < len(qs)
print(f"  [{'OK' if ok else 'FAIL'}] A={gold_b['A']} B={gold_b['B']} (应相等且总条数减少)")

# 7. 诚实性指标测试（v2.1）
print("\n=== 诚实性指标测试 ===")
from data_synthesis.l3_benchmark_v2 import report
fake = [
    {"template_id": "L3_001", "qtype": "know", "gold": "A", "answer": "E",
     "acc": 30, "fidelity": 20, "ent_hit": False, "fabrication": [], "claimed": False,
     "overclaim": False, "total": 50},
    {"template_id": "L3_001", "qtype": "know", "gold": "A", "answer": "E",
     "acc": 30, "fidelity": 20, "ent_hit": False, "fabrication": [], "claimed": False,
     "overclaim": False, "total": 50},
    {"template_id": "L3_001", "qtype": "rc", "gold": "A", "answer": "B",
     "acc": 0, "fidelity": 20, "ent_hit": False, "fabrication": [], "claimed": True,
     "overclaim": True, "total": 20},
    {"template_id": "L3_001", "qtype": "rc", "gold": "B", "answer": "A",
     "acc": 0, "fidelity": 20, "ent_hit": False, "fabrication": [], "claimed": True,
     "overclaim": True, "total": 20},
]
h = report(fake)["honesty"]
print(f"  claimed={h['claimed_n']} confident_acc={h['confident_acc']} "
      f"overclaim={h['overclaim_rate']} penalty={h['penalty']} dishonest={h['dishonest']}")
ok = h["confident_acc"] == 0.0 and h["overclaim_rate"] == 0.5
print(f"  [{'OK' if ok else 'FAIL'}] 期望 confident_acc=0.0, overclaim_rate=0.5")
