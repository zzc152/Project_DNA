# Gold 标注修复报告（fix_gold.py）

> 修复日期：2026-08-08
> 前置：`docs/gold_validation.md`（gold 验证报告，T1-T5 发现 37 处问题）
> 工具：`scripts/fix_gold.py`（可复现修复脚本）
> 修复产物：`data/synthetic/drafts_full2_clean_fixed.jsonl`、`data/processed/knowledge_base_clean_fixed.jsonl`

## 修复摘要

| 项目 | 修复前 | 修复后 | 剩余（均为已标记待处理项） |
|------|--------|--------|--------------------------|
| T1 证据字段矛盾 | 5 条 | **2 条** | L3_002_26 / L3_003_29（DELCAND，因子不可考） |
| T2 KB 方向冲突 | 3 对 | **2 对** | MYCN→promoter / MYCN→enhancer（FLAG，待人工核查） |
| T4 非法 effect | 6 条 | **3 条** | 3 条 DELCAND（共现/相关非因果，待删除） |
| T3 因子名 | 23 条 | 23 条（已标记） | `metadata.factor_valid=false`，待清洗 |

- **FIX（高置信修复）**：L3 6 条 + KB 5 条
- **DELCAND（删除候选标记）**：L3 3 条 + KB 3 条（保留记录、加 `fix_status=delete_candidate`，供人工确认后删除）
- **FLAG（人工核查标记）**：KB 11 条（MYCN 两对全部记录，加 `fix_status=flag`）
- **T3 标记**：全部 drafts 66 条（含 L3 23 条 + L1/L2/L4/L5 原始记录），加 `factor_valid=false`

---

## 一、L3 高置信方向修复（6 条，gold 直接改变）

| id | 修复前 effect (gold) | 修复后 effect (gold) | 摘要证据（修复依据） |
|----|---------------------|---------------------|---------------------|
| L3_001_8 | increases (**A**) | **decreases (B)** | "euchromatinization of the globin promoter and enhanced NRF2 binding **in RepID-depleted cells**" → RepID 缺失→活性↑ → RepID 本身**抑制**启动子 |
| L3_001_15 | increases (**A**) | **decreases (B)** | "reduction of ABO expression is strongly associated with **DNA methylation** of the ABO promoter" → 甲基化与表达降低相关 → **抑制** |
| L3_002_55 | decreases (**B**) | **increases (A)** | "Targeted **depletion** of enhancer RNAs ... **silences** BCL11A" → eRNA 缺失→沉默 → eRNA 本身**增强** enhancer 活性 |
| L3_001_33 | decreases (**B**) | **increases (A)** | "**loss of H3K4me1 reduced** H3K27ac at the β-globin enhancer" → H3K4me1 缺失→活性↓ → H3K4me1 本身**增强** |
| L3_002_8 | —（字段矛盾） | direction: modulates_activity → **decreases_activity** | "macroH2A variants localized at enhancer elements ... **devoid of H3K27ac** ... a role for macroH2A at **inactive enhancers**" → 方向=decreases，与 effect 对齐 |
| L3_002_37 | increases (**A**) | **modulates (D)** | 关系降级："ASCL1-**bound** and **primed** neuronal enhancers"，激活需 G1 arrest → 仅结合/预备，方向不明（与 reasoning_chain[4] 修正逻辑一致） |

> 说明：前 4 条同时被 T5（32B 独立判断）证实——模型从同一摘要独立推出方向与修复后 gold 一致。

## 二、L3 删除候选标记（3 条，DELCAND）

| id | 问题 | 依据 |
|----|------|------|
| L3_002_26 | 因子 "epigenetic status" 不可考 | 摘要实体为风险等位基因 G（SNP），非 "epigenetic status"；与 L3_003_29 同 PMID 39998882 重复 |
| L3_003_29 | 同上 | 同上 |
| L3_002_2 | 题目语义反向 | 摘要讲 "BCL11A erythroid enhancer"（BCL11A 基因的增强子）被破坏损害转录 = enhancer 作用于 BCL11A，非 BCL11A 作用于 enhancer |

> 建议：3 条在下一轮数据清洗时删除（当前保留并标记，避免评测集数量骤变）。

## 三、KB 修复（5 条 FIX + 3 条 DELCAND + 11 条 FLAG）

### FIX（5 条，含 DNA methylation K562 2 条）

| 记录 | 修复前 effect | 修复后 effect | 依据 |
|------|--------------|--------------|------|
| DNA methylation→promoter (K562) ×2 | increases | **decreases** | 与 L3_001_15 同源错误：ABO 启动子甲基化→表达降低；与 not_specified 的 decreases×2 冲突消除 |
| GATA motif→OGT/OGA promoter interaction | affects_interaction | **decreases** | claim 明确："GATA 位点突变（GAGA）显著降低 OGT 和 OGA 的 promoter 相互作用" |
| caQTL variants→TF binding motifs | disrupts | **decreases** | disrupts=破坏结合基序→抑制方向，映射到规范词表 |
| OsDDE9→Ghd7 and Ehd1 expression | regulates_expression | **modulates** | regulates=调节但方向未明→modulates |

### DELCAND（3 条，共现/相关非因果，不适合方向题 gold）

| 记录 | effect | 问题 |
|------|--------|------|
| DNA methylation→gene upregulation | co_occurs_with | 共现非因果 |
| chromatin accessibility→promoter methylation level | correlates | 相关非因果 |
| targeted activation of this regulatory element→AML differentiation | promotes | factor 描述性不可考 + element 为分化过程非元件 |

### FLAG（11 条，MYCN 两对全部记录，待人工核查）

- **MYCN→promoter**（6 条）：required_for×2 + increases×2 + decreases×2 → A/B 类冲突
- **MYCN→enhancer**（5 条）：modulates + decreases×2 + increases×2 → 3 种方向

> 说明：MYCN 作为 TF 对靶位点可能同时存在激活与抑制（招募 G9a 抑制神经元分化基因，如 L3_002_60），需逐条回原文裁决，故仅标记不擅自修改。

## 四、T3 因子名标记（23 条 L3）

已全部加 `metadata.factor_valid=false` + fix_note，典型不可考因子名：

- `acetylating and de-methylating agents`（×4）
- `H3 acetylation at lysine 4 and 9 (K4 and K9)`（×2）
- `epigenetic status`（×2，已 DELCAND）
- `TLR4 downstream NF-κB p65 and c-Jun`（×2）
- `GATA1, TAL1, and CCAAT/enhancer-binding protein (C/EBPβ)`（×2）
- `RET-dependent transcription factors`（×2）
- `adipogenic lineage-specific transcription factors (TFs)`（×2）

> 建议：下一轮清洗时替换为可考实体名或删除对应题目。

## 五、修复前后对比（verify_gold.py 复测）

```text
                   修复前      修复后
T1 证据矛盾         5 条   →   2 条（剩余=DELCAND）
T2 KB 方向冲突       3 对   →   2 对（剩余=FLAG）
T3 可疑因子名       23 条   →  23 条（已标记 factor_valid=false）
T4 不可映射 effect   6 条   →   3 条（剩余=DELCAND）
```

**修复后无未标记的遗留问题**——剩余全部为已显式标记待处理项。

## 六、复现方式

```bash
# 1. 修复（生成 _fixed.jsonl）
cd /workspace/zzc/BioDesign-Agent
PYTHONPATH=src python scripts/fix_gold.py

# 2. 复测
PYTHONPATH=src python scripts/verify_gold.py \
  --drafts data/synthetic/drafts_full2_clean_fixed.jsonl \
  --kb data/processed/knowledge_base_clean_fixed.jsonl \
  --out tmp/gold_fixed_report.json

# 3. T5 语义复测（可选，需 GPU）
PYTHONPATH=src python scripts/verify_gold.py \
  --drafts data/synthetic/drafts_full2_clean_fixed.jsonl \
  --model models/Qwen2.5-32B-AWQ --device cuda:0 --t5-limit 63 \
  --out tmp/gold_fixed_t5.json
```

产物：`data/synthetic/drafts_full2_clean_fixed.jsonl`、`data/processed/knowledge_base_clean_fixed.jsonl`、`tmp/fix_report.json`（逐条修复明细）。
