# Gold 标注质量验证报告（verify_gold.py）

> 验证日期：2025-08-08（session 8）
> 数据源：`data/synthetic/drafts_full2_clean.jsonl`（896 条，L3 192 条）、`data/processed/knowledge_base_clean.jsonl`（617 条）
> 工具：`scripts/verify_gold.py`（T1-T5）
> 模型：Qwen2.5-32B-AWQ（仅 T5 用，seed=42，63 条 rc 型抽样）

## 结论摘要

| 测试 | 内容 | 结果 | 判定 |
|------|------|------|------|
| T1 | L3 证据字段内部一致性（effect vs evidence.direction） | **5/192 条矛盾** | ⚠️ FAIL |
| T2 | KB 实体对方向冲突 | **3/253 对冲突** | ⚠️ FAIL |
| T3 | L3 题面因子名合法性 | **23/192 条可疑** | ⚠️ FAIL |
| T4 | KB effect → gold 可映射性 | **6 条非法 effect** | ⚠️ FAIL |
| T5 | 摘要↔gold 语义一致性（32B 独立判断） | **36/63 (57.1%)** | ❌ FAIL |

**总体判定：gold 标注存在系统性质量问题，用户质疑成立。** 数据不适合直接用于可靠评测，需修复后重新标注。

---

## T1：L3 证据字段内部一致性（5/192 条矛盾）

`entities_used.effect` 与 `evidence.direction` 归一化后不一致：

| id | 因子 | 元件 | effect | evidence.direction | 判定 |
|----|------|------|--------|-------------------|------|
| L3_001_8 | RepID | globin promoter | increases (A) | modulates_activity (D) | **gold 错**：摘要"RepID-depleted cells 中启动子常染色质化+NRF2 结合增强"→ RepID 缺失增强活性 → 应为 B（decreases） |
| L3_002_8 | macroH2A histone variants | enhancer | decreases (B) | modulates_activity (D) | 字段不统一 |
| L3_002_26 | epigenetic status | promoter | increases (A) | modulates_activity (D) | 字段不统一（方向字段保守化） |
| L3_002_37 | ASCL1 | neuronal enhancers | increases (A) | binds_target (D) | gold=A 合理（摘要支持"drive neuronal differentiation"），direction 字段应改为 increases_activity |
| L3_003_29 | epigenetic status | promoter | increases (A) | modulates_activity (D) | 同 L3_002_26 |

**关键发现：L3_001_8 被 T5 独立证伪**——32B 从同一摘要独立判断为 B，与摘要原文语义一致，数据标 A 错误。

---

## T2：KB 实体对方向冲突（3/253 对）

| factor | element | 冲突内容 | 备注 |
|--------|---------|---------|------|
| DNA methylation | promoter | decreases×2 (not_specified) + increases×2 (K562) | **K562 的 increases 被 L3_001_15 证伪**（见 T5），KB 与 L3 数据同源错误 |
| MYCN | promoter | required_for×2 + increases×2 (均 not_specified) | required_for 与 increases 语义接近，可归并为 A 类，不算硬冲突 |
| MYCN | enhancer | modulates + decreases×2 + increases | **最严重**：3 种方向共存 |

---

## T3：L3 因子名合法性（23/192 条可疑）

典型不可考因子名（描述性短语混入实体位）：

- `acetylating and de-methylating agents`（×4 条）
- `H3 acetylation at lysine 4 and 9 (K4 and K9)`（×2）
- `epigenetic status`（×2）
- `TLR4 downstream NF-κB p65 and c-Jun`（×2）
- `GATA1, TAL1, and CCAAT/enhancer-binding protein (C/EBPβ)`（×2）
- `RET-dependent transcription factors`（×2）
- `adipogenic lineage-specific transcription factors (TFs)`（×2）
- 等 23 条

这些"因子"不是单一可考实体，题目可答性存疑（模型无法从摘要中锁定对应实体）。

---

## T4：KB effect → gold 可映射性（6 条非法 effect）

`_gold_direction()` 无法映射到 A-D 的 effect 值：

| factor | element | effect |
|--------|---------|--------|
| targeted activation of this regulatory element | AML differentiation | promotes |
| GATA motif | OGT/OGA promoter interaction | affects_interaction |
| DNA methylation | gene upregulation | co_occurs_with |
| chromatin accessibility | promoter methylation level | correlates |
| OsDDE9 | Ghd7 and Ehd1 expression | regulates_expression |
| caQTL variants | TF binding motifs | disrupts |

这些值不在规范方向词表（increases/decreases/modulates/required_for/no_effect）内，gold 无法推导。

---

## T5：摘要↔gold 语义一致性（32B，36/63 = 57.1%）

**设计**：只给摘要原文句（`摘要原文证据` 字段），32B 独立判断方向（A-E 五选），与数据 gold 对比。seed=42 抽 63 条 rc 型。

### 不一致分类（27 条）

| 类别 | 数量 | 说明 |
|------|------|------|
| 摘要信息不足（pred=E，gold 有方向） | 13 | 摘要单句不足以推出方向，gold 存在"过度自信"标注 |
| 方向真冲突（pred 与 gold 都是 A/B/D 且不同） | 14 | 其中至少 3-4 条摘要明确支持 pred，**gold 疑似错误** |
| gold=E 但模型有方向 | 0 | — |

### 各 gold 类别一致率

| gold | 一致率 |
|------|--------|
| A (increases) | 20/34 (59%) |
| B (decreases) | 14/22 (64%) |
| D (modulates) | 2/7 (29%) |

D 类题目一致率最低：模型倾向猜具体方向，说明 **D 类（方向不明）题目的区分度差**，需在题面设计上强化。

### 确认/疑似 gold 错误的条目（摘要证据核实）

| id | 因子 | 元件 | 标 gold | 摘要证据 | 正确方向 |
|----|------|------|--------|---------|---------|
| L3_001_8 | RepID | globin promoter | A | "euchromatinization ... in RepID-depleted cells" → RepID 缺失增强活性 | **B** |
| L3_001_15 | DNA methylation | promoter (K562) | A | "reduction of ABO expression is strongly associated with DNA methylation" → 甲基化与表达降低相关 | **B** |
| L3_002_55 | enhancer RNAs | enhancer | B | "depletion of eRNAs silences BCL11A" → eRNA 缺失抑制 → eRNA 增强 | **A** |
| L3_001_33 | H3K4me1 | enhancer | B | "loss of H3K4me1 reduced H3K27ac" → 丢失降低活性 → H3K4me1 增强 | **A** |
| L3_002_2 | BCL11A | enhancer | A | 摘要描述的是"BCL11A erythroid enhancer"被破坏损害转录——BCL11A 是 enhancer 的靶基因而非作用因子，**题目语义可能反向** | 需重查 |

**注**：上表 4 条为高置信 gold 错误；其余冲突多为 A↔D 或 B↔D（模型与标注对"方向是否明确"判断不同），属可讨论项。

---

## 修复建议

1. **立即修复**（5 条 high-confidence 错误）：
   - L3_001_8：gold A→B，effect/evidence.direction 同步改为 decreases
   - L3_001_15：gold A→B；并同步修正 KB 中 DNA methylation→promoter (K562) 的 increases 条目
   - L3_002_55：gold B→A
   - L3_001_33：gold B→A
   - L3_002_2：核查题目语义，若 factor/element 关系反向则删除或重写

2. **字段统一**：L3_002_8 / L3_002_26 / L3_002_37 / L3_003_29 的 evidence.direction 与 effect 对齐（优先采用明确方向值）。

3. **KB 冲突处理**：
   - MYCN→enhancer（3 种方向）：逐条回原文核实，保留有摘要支持的
   - DNA methylation→promoter：修正 K562 increases（见上）

4. **因子名清洗**：23 条描述性因子名改为可考实体名，或标记为不可用并删除对应题目。

5. **非法 effect**：6 条映射到规范词表或删除记录。

6. **D 类题目**：gold=D 一致率仅 29%，建议将"方向不明"类题目的摘要证据补强，或减少 D 类占比。

---

## 复现方式

```bash
# T1-T4（无模型）
cd /workspace/zzc/BioDesign-Agent
PYTHONPATH=src python scripts/verify_gold.py --out tmp/gold_report.json --audit tmp/gold_audit.md

# T5（需 GPU）
PYTHONPATH=src python scripts/verify_gold.py --model models/Qwen2.5-32B-AWQ --device cuda:0 \
  --t5-limit 63 --out tmp/gold_t5.json
```

产物：`tmp/gold_audit_report.json`、`tmp/gold_audit_t5.json`、`tmp/gold_audit_list.md`（30 条人工抽检清单）
