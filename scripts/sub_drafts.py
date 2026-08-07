# -*- coding: utf-8 -*-
"""用更新后的模板 + 原 KB 记录重新填充问题记录 instruction，写临时 drafts 供小批量验证。"""
import json
import sys

sys.path.insert(0, "src/data_synthesis")
from template_engine import load_templates, fill_template

tpls = {t["id"]: t for t in load_templates("templates")}
print("模板数:", len(tpls))

kb = [json.loads(l) for l in open("data/processed/knowledge_base_clean.jsonl", encoding="utf-8")]
drafts = [json.loads(l) for l in open("data/synthetic/drafts.jsonl", encoding="utf-8")]
by_id = {d["id"]: d for d in drafts}

ids = [
    "L4_001_1", "L4_001_2", "L4_001_3",
    "L4_003_1", "L4_003_3", "L4_003_4",
    "L5_001_1", "L5_001_2",
    "L3_003_1",
]

out = []
for rid in ids:
    d = by_id[rid]
    t = tpls[d["template_id"]]
    # 原 knowledge 记录：source_claim_ids[0] 是 KB 列表索引
    src = d["metadata"]["source_claim_ids"][0]
    record = kb[src]
    new_inst = fill_template(record, t["instruction_template"])
    d2 = dict(d)
    d2["instruction"] = new_inst
    out.append(d2)

with open("tmp/drafts_sub.jsonl", "w", encoding="utf-8") as f:
    for d in out:
        f.write(json.dumps(d, ensure_ascii=False) + "\n")
print("写出", len(out), "条 -> tmp/drafts_sub.jsonl")

# 对照新旧 instruction（自补越界模板 L4_001）
d = out[0]
print("=== 示例", d["id"], "===")
print("[新]", d["instruction"][:260])
print("[旧]", by_id[d["id"]]["instruction"][:120])
