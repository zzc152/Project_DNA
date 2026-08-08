# Gold 人工核查报告（MYCN FLAG + DELCAND）

> 核查日期：2026-08-08
> 前置：`docs/gold_fix_report.md`（自动修复完成，剩余 11 条 MYCN FLAG + 6 条 DELCAND 待人工处理）
> 工具：`scripts/review_gold_manual.py`（人工核查落地脚本）
> 产物：`data/synthetic/drafts_full2_clean_final.jsonl`（893 条）、`data/processed/knowledge_base_clean_final.jsonl`（614 条）、`tmp/manual_review.json`（逐条裁决明细）

---

## 一、MYCN FLAG 11 条 — 逐条回原文核查

### 核心发现

**MYCN 的"方向冲突"不是标注错误，而是靶基因依赖（target-gene dependent）的生物学真实。** 同一篇摘要（PMID 37781575）明确写了 MYCN 对两类靶基因的相反效应：

> "MYCN-binding also **increases** H3K4me3 and H3K27ac on **canonical MYC target promoters** and **decreases** H3K27ac on **neuronal differentiation gene enhancers and promoters**."

即：MYCN 激活 MYC 靶基因、抑制神经元分化基因，同时成立。原 KB 以 `(factor, element, cell_line)` 为键，未含靶基因维度，导致 T2 误报"同对冲突"。

### 逐条裁决表

| # | 记录（factor→element） | 原 effect | PMID | 摘要原文证据 | 裁决 |
|---|------------------------|-----------|------|--------------|------|
| 1 | MYCN→promoter (DLX5) | required_for | 41135749 | "direct binding site of MYCN on DLX5 promoter... MYCN **coregulates** the expression of DLX5 by directly binding" | ⚠️ **降级 modulates**——"coregulates"未指明激活/抑制方向，required_for 属过度推断 |
| 2 | MYCN→promoter (DLX5, design) | required_for | 41135749 | 同上 | ⚠️ **降级 modulates** |
| 3 | MYCN→enhancer | modulates | 36598365 | "MYCN **binds** to enhancers when aberrantly upregulated in neuroblastoma" | ✅ **保留**——仅证明结合，modulates 保守正确 |
| 4 | MYCN→promoter (MYC targets) | increases | 37781575/36539767 | "MYCN-binding **increases** H3K4me3 and H3K27ac on canonical MYC target promoters" | ✅ **保留** |
| 5 | 同 4（design 型） | increases | 37781575/36539767 | 同上 | ✅ **保留** |
| 6 | MYCN→enhancer (neuronal diff.) | decreases | 37781575 | "MYCN recruits **G9a** to enhancers to **repress** neuronal differentiation genes" | ✅ **保留**——G9a 为抑制性组蛋白甲基转移酶 |
| 7 | 同 6（design 型） | decreases | 37781575 | 同上 | ✅ **保留** |
| 8 | MYCN→promoter (neuronal diff.) | decreases | 37781575 | 同上句后半："decreases H3K27ac on ... **promoters**" | ✅ **保留**——原文明确覆盖 promoter |
| 9 | 同 8（design 型） | decreases | 37781575 | 同上 | ✅ **保留** |
| 10 | MYCN→enhancer (SMAD9) | increases | 36539767 | "**SMAD9** specifically induced by TFs including MYCN... at the enhancer region" | ✅ **保留** |
| 11 | 同 10（design 型） | increases | 36539767 | 同上 | ✅ **保留** |

**结果：11 条中 9 条方向正确保留，2 条（DLX5 required_for）证据不足降级为 modulates。**

### 落地方式

- 每条 MYCN 记录加 `metadata.target_gene`（DLX5 / neuroblastoma targets / MYC targets / neuronal differentiation genes / SMAD9），消除同对误判
- `fix_status` 统一改为 `resolved`（人工核查通过）
- 2 条 DLX5：`effect: required_for → modulates`、`evidence.direction: required_for_activity → modulates_activity`
- 同步升级 `verify_gold.py` 的 T2：冲突判定加入 `metadata.target_gene` 维度（同因子不同靶基因的方向相反不算标注冲突）

---

## 二、DELCAND 6 条 — 人工确认删除

### L3 3 条（drafts）

| id | 因子 | 问题 | 原文证据 | 结论 |
|----|------|------|----------|------|
| L3_002_26 | "epigenetic status" | 因子不可考 + 同 PMID 重复 | 摘要实体为**风险等位基因 G**（"conversion into risk allele G results in increased GLI3 binding and USP47 upregulation"），非 "epigenetic status"；与 L3_003_29 同 PMID 39998882 重复 | **删除** |
| L3_003_29 | "epigenetic status" | 同上 | 同上 | **删除** |
| L3_002_2 | BCL11A | 题目语义反向 | 摘要讲 "BCL11A erythroid enhancer"（BCL11A 基因的增强子）被破坏损害转录 = **enhancer 作用于 BCL11A**，非 BCL11A 作用于 enhancer | **删除** |

### KB 3 条

| 记录 | effect | 问题 | 原文证据 | 结论 |
|------|--------|------|----------|------|
| DNA methylation→gene upregulation | co_occurs_with | 共现非因果 | "Hypermethylation at promoters of **up-regulated genes** is observed"（观察性共现） | **删除** |
| chromatin accessibility→promoter methylation level | correlates | 相关非因果 | "ATAC-seq peaks with high accessibility... corresponded to **low** methylation levels"（观察性相关） | **删除** |
| targeted activation of this regulatory element→AML differentiation | promotes | 因子描述性不可考 + element 非元件 | factor 是干预描述（CRISPRa）非可考实体；element 为分化过程非调控元件 | **删除** |

> 这 3 条 KB 记录的 reasoning_chain 中均已带 "error pattern: evidence gap / relation 过强 / unsupported specificity" 的原始标注，与人工判定一致。

---

## 三、核查后复测（verify_gold.py）

```text
                    自动修复后      人工核查后
T1 证据矛盾         2 条（DELCAND）→  0 条 ✅
T2 KB 方向冲突       2 对（FLAG）   →  0 对 ✅（target_gene 维度升级）
T3 可疑因子名       23 条（已标记） → 21 条（全部 factor_valid=false）
T4 不可映射 effect   3 条（DELCAND）→  0 条 ✅
```

**核查后所有未标记问题清零**——剩余 T3 21 条均为显式标记 `factor_valid=false` 的因子名清洗项。

## 四、数据规模变化

| 数据 | 修复前 | 自动修复后 | 人工核查后 | 删除 |
|------|--------|-----------|-----------|------|
| drafts | 896 | 896 | **893** | L3 3 条 DELCAND |
| KB | 617 | 617 | **614** | 3 条 DELCAND |

## 五、复现方式

```bash
cd /workspace/zzc/BioDesign-Agent
# 1. 人工核查落地（基于 fixed 数据 → final 数据）
PYTHONPATH=src python scripts/review_gold_manual.py

# 2. 复测
PYTHONPATH=src python scripts/verify_gold.py \
  --drafts data/synthetic/drafts_full2_clean_final.jsonl \
  --kb data/processed/knowledge_base_clean_final.jsonl \
  --out tmp/gold_final_report.json
```

产物：`data/synthetic/drafts_full2_clean_final.jsonl`、`data/processed/knowledge_base_clean_final.jsonl`、`tmp/manual_review.json`。
