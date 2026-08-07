#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
fix_v3.py — 二轮收尾修正（基于 fix_v2 输出，625 行，新 idx 行号）

处置原则（用户决策：只处理可优化的，其余删除）：
- 修正 3 条（可优化）:
    #43  QTLs modulates → affects（association 降级，claim 与证据文本对齐）
    #212 ASCL1 required_for → binds（证据为 G1 arrest 对 enhancer 必需，ASCL1 仅结合）
    #225 ATF4 chain 因果强度修正句改写为正向自洽表述
- 删除 5 条（Qwen 反推局限、人工已确认成立但验证无法通过）:
    #86  H3K27me3 decreases（loss-of-function 反推）
    #109 BCL11A no_effect（Qwen 语义局限）
    #137 ETO2 decreases（loss-of-function 反推）
    #271 Sp3 required_for（binding≠required_for 反推 + chain）
    #285 ERRα decreases（loss-of-function 反推）

注意：本脚本输入为 fix_v2 输出（625 行），idx 为 fix_v2 后的新行号。
"""
import json
from pathlib import Path

KB = Path("data/processed/knowledge_base_clean.jsonl")

# ---- 二轮删除（新 idx）----
DELETES = {86, 109, 137, 271, 285}

# ---- 二轮修正（新 idx）----
FIXES: dict[int, dict] = {
    43: {
        "claim": "文献证据表明，在未指定细胞系的实验中，遗传变异（QTLs）影响（affects）SPI1、GABPB和STAT3的预测结合位点；机制上，Many detected QTLs influence the predicted binding of myeloid transcription factors, including SPI1, GABPB and STAT3.（来源：PMID 39870618）。",
        "entities": {
            "factor": "QTLs",
            "factor_type": "variant",
            "regulatory_element": "predicted binding sites of SPI1/GABPB/STAT3",
            "effect": "affects",
        },
        "evidence": {"direction": "affects_binding"},
        "reasoning_chain": [
            "来源: PMID 39870618",
            "摘要原文证据: many detected QTLs overlap disease susceptibility loci and influence the predicted binding of myeloid transcription factors, including SPI1, GABPB and STAT3.",
            "文献机制描述: Many detected QTLs influence the predicted binding of myeloid transcription factors, including SPI1, GABPB and STAT3.",
            "因子类型判定: variant；细胞系: not_specified；元件: SPI1/GABPB/STAT3 预测结合位点",
            "association≠mechanism: variant→predicted binding change 不能推出 TF→enhancer activity；结论限定为「影响预测结合位点」而非「调节活性」（error pattern: unsupported specificity / overclaim）",
        ],
    },
    212: {
        "claim": "文献证据表明，在未指定细胞系的实验中，因子ASCL1结合（binds）神经元增强子；机制上，Prolonged G1 arrest is required to activate ASCL1-bound and primed neuronal enhancers to drive neuronal differentiation.（来源：PMID 40452575）。",
        "entities": {
            "factor": "ASCL1",
            "factor_type": "TF",
            "regulatory_element": "neuronal enhancers",
            "effect": "binds",
        },
        "evidence": {"direction": "binds_target"},
        "reasoning_chain": [
            "来源: PMID 40452575",
            "摘要原文证据: Prolonged G1 arrest is further required to activate these ASCL1-bound and primed neuronal enhancers to drive neuronal differentiation.",
            "文献机制描述: ASCL1 binds and primes neuronal enhancers; prolonged G1 arrest is required to activate them.",
            "因子类型判定: TF；细胞系: not_specified；元件: neuronal enhancers",
            "关系粒度修正: 证据仅支持 ASCL1 结合/预备（bound and primed）enhancer；「必需」的是 G1 arrest 而非 ASCL1，故降级为 binds（error pattern: evidence gap / relation overreach）",
        ],
    },
    225: {
        "reasoning_chain": [
            "来源: PMID 38192280",
            "摘要原文证据: ATF4 occupancy increases at the ATF5 promoter site in brains isolated from TMG treated mice suggesting that O-GlcNAc is regulating ATF4 targeted gene expression.",
            "文献机制描述: ATF4 occupancy increases at the ATF5 promoter site in brains isolated from TMG treated mice suggesting that O-GlcNAc is regulating ATF4 targeted gene expression.",
            "因子类型判定: TF；体系: in vivo tissue（mouse brain）；元件: ATF5 promoter",
            "实验体系修正: 证据为小鼠脑组织体内实验，原抽取写 SK-N-SH 细胞系（error pattern: incorrect extraction / context grounding）",
            "结论限定: occupancy 增加支持 ATF4 与 ATF5 promoter 结合（binds），不强推 required_for promoter activity",
        ],
    },
}


def deep_merge(base: dict, patch: dict) -> dict:
    """dict 递归合并（patch 覆盖 base 的同名键；嵌套 dict 递归）。"""
    out = dict(base)
    for k, v in patch.items():
        if isinstance(v, dict) and isinstance(out.get(k), dict):
            out[k] = deep_merge(out[k], v)
        else:
            out[k] = v
    return out


def main() -> None:
    rows = [json.loads(l) for l in KB.open(encoding="utf-8")]
    report = []
    print(f"[fix_v3] 输入 {len(rows)} 条")

    for i, r in enumerate(rows):
        if i in DELETES:
            continue
        fix = FIXES.get(i)
        if fix is None:
            continue
        merged = deep_merge(r, fix)
        rows[i] = merged
        report.append({"idx": i, "new_claim": merged.get("claim"), "changed_fields": list(fix.keys())})

    rows = [r for i, r in enumerate(rows) if i not in DELETES]
    out = "\n".join(json.dumps(r, ensure_ascii=False) for r in rows) + "\n"
    KB.write_text(out, encoding="utf-8")
    print(f"[fix_v3] 输出 {len(rows)} 条（删除 {len(DELETES)}，修正 {len(report)}）")
    for d in report:
        print(f"  - 修正 #{d['idx']}: {d['changed_fields']}")
    print(f"[fix_v3] 删除 idx: {sorted(DELETES)}")


if __name__ == "__main__":
    main()
