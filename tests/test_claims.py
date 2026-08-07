# -*- coding: utf-8 -*-
"""临时验证脚本：本地合成数据验证 build_literature_claims 逻辑"""
import json
import sys
from pathlib import Path

sys.path.insert(0, r"d:\code_home\project")

test = [
    {"pmid": "111", "parsed": True, "findings": [
        {"cell_line": "K562", "regulatory_element": "enhancer", "factor": "GATA1",
         "effect": "required_for",
         "mechanism": "GATA1 occupancy is necessary for enhancer function",
         "evidence_sentence": "Deletion of GATA1 sites abolished expression"},
        {"cell_line": "not_specified", "regulatory_element": "promoter",
         "factor": "nucleosome occupancy", "effect": "decreases",
         "mechanism": "High occupancy represses promoter activity",
         "evidence_sentence": "high nucleosome occupancy represses"},
    ]},
    {"pmid": "222", "parsed": True, "findings": [
        {"cell_line": "K562", "regulatory_element": "enhancer", "factor": "GATA1",
         "effect": "required_for", "mechanism": "",
         "evidence_sentence": "GATA1 needed again"},
    ]},
]

tmp = Path(r"d:\code_home\project\_test_extractions.jsonl")
tmp.write_text("\n".join(json.dumps(x, ensure_ascii=False) for x in test), encoding="utf-8")

from scripts.build_literature_claims import (
    load_extractions,
    build_association_claim,
    build_design_rule_claim,
)

fs = load_extractions(tmp)
print("合并 findings 数:", len(fs))
for nf in fs:
    c = build_association_claim(nf)
    print("---", c["claim_type"], "conf=", c["confidence"])
    print(c["claim"])
    dr = build_design_rule_claim(nf)
    if dr:
        print("   DR:", dr["claim"])
tmp.unlink()
print("OK")
