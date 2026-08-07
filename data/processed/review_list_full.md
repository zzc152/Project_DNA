# 科学断言知识库 —— 验证复核完整报告（自包含版）

> **本文件自包含**：背景、流程、数据统计、验证方法、验证结果、34 条问题记录的前因后果与处置方案，全部在此，无需再查其他文件。
> 生成日期：2026-08-06（**v2**：吸收人工审核意见，升级 relation 粒度 + error pattern 标签 + 目标 claim schema）｜ 数据：`data/processed/knowledge_base_clean.jsonl`（631 条）
> 配套表格：`review_list.csv`（同 34 条，Excel 可直接打开）

---

## 一、背景：这条知识库是怎么来的（前因）

**模块管线**：Module 1（文献抽取原始断言）→ 本流程（去重 → 清洗 → 验证）→ Module 2（模板填充）

1. **Module 1 产出**：从 PubMed 文献摘要抽取"因子 → 对 promoter/enhancer 活性的作用"类科学断言，写入 `knowledge_base_raw.jsonl`。
2. **去重**（本流程第 1 步）：精确去重 + 实体指纹去重 + 语义去重，共合并 2 条语义重复记录；另识别**互补关系** 6 条（同 PMID 内互补）+ 8 条（跨 PMID 互补）。
3. **清洗**（本流程第 2 步）：统一 factor/effect/cell_line 字段格式、规范方向术语、补齐来源 PMID，产出 `knowledge_base_clean.jsonl`。
4. **验证**（本流程第 3 步）：三层验证（见下节），产出 `validation_report.jsonl` + `validation_summary.json` + 本复核清单。

**最终数据规模**：631 条 = 150 条统计来源（大规模 QTL/eQTL 统计证据）+ 481 条文献来源；平均置信度 **0.8929**。

---

## 二、验证方法（三层验证）

| 层 | 检查项 | 工具 | 结果 |
|---|---|---|---|
| A | **统计可靠性** check_stat_record：样本量、效应量、P 值 | 规则脚本 | 150/150 通过，0 失败 |
| B1 | **文献一致性** check_lit_consistency：claim 与 evidence 语句主体/方向是否一致 | 规则脚本 | 481/481 通过，0 失败 |
| B2+C | **Qwen 大模型判定**：对每条 claim 给 supported / unsupported / unclear；并对 design_rule 检查推理链自洽 chain_consistent | Qwen2.5-7B-Instruct（CUDA） | 见下节 |

**Qwen V2 prompt 关键设计**：按 claim_type 区分判定标准——
- `mechanistic`：要求证据逐字支持（严格）；
- `design_rule`：只验证前提（底层关系 + 方向正确即 supported）。

---

## 三、验证结果（后果 · 现状）

| 指标 | 数值 |
|---|---|
| 输入总数 | 631 |
| 统计检查失败 | 0 |
| 文献一致性失败 | 0 |
| Qwen 已判定 | 481 |
| supported | 448（93.1%） |
| unsupported | 31（6.4%） |
| unclear | 2（0.4%） |
| 推理链自洽 | 471/481（97.9%） |
| 问题记录（需人工复核） | **34 条**（unsupported 31 + unclear 2 + chain 不自洽额外命中 1） |

**报告完整性修复**：`validation_report.jsonl` 曾因程序化阶段与 Qwen 阶段重复 append 产生 962 行；已修复为按 idx 合并（保留 stat_checks/checks_passed），重跑后为 **631 行**（481 文献 + 150 统计，无重复）。

---

## 四、34 条问题记录 · 逐条前因后果

### 分类汇总（v2：error pattern 标签）

| 类别 | 数量 | error pattern | 最终处置 |
|---|---|---|---|
| A. Qwen 误判-可放行 | 6 | 模型逻辑局限（非抽取错误） | ✅ 保留 + 补 chain 反推句 |
| B. 方向/主客体反 | 5 | **incorrect extraction**（causal direction） | ⚠️ 修正（schema 化 + relation 粒度） |
| C. 实验体系错 | 3 | **incorrect extraction**（context grounding） | ⚠️ 修正（system/relation/强度） |
| D. 过度泛化 | 6 | **unsupported specificity / overclaim** | 🔍 降级/限定范围 |
| E. 证据间接 | 8 | **evidence gap**（relation 过强） | 🔧 改变 relation 类型 |
| F. 证据不相关 | 2 | **hallucinated biological relation** | ❌ 删除 |
| G. 推理链缺陷 | 4 | **reasoning chain risk** | ⚠️ 修正 chain + 降 design 强度（#447/#449 或删） |

### 目标 Claim Schema（v2 升级）

**核心教训**：不要只修 direction，还要检查 **relation 粒度**——"expression → activity"最容易二次过推。修正时统一按以下 schema 建模：

| 字段 | 说明 | 示例（#258） |
|---|---|---|
| factor | 因子/实体 | HDAC2 |
| relation | 关系类型（非仅方向） | negative_regulator |
| target | 靶标 | hSVCT2 expression |
| target_level | 靶标层级 | expression（≠ promoter activity） |
| evidence_type | 证据类型 | knockdown |
| direction | 方向 | decrease |

**修正铁律**：证据只到 expression 就写 expression；只到 binding/interaction 就写 binding；只到 association 就写 correlation。**不跨层级推断**。

---

### A 类：Qwen 误判 → 保留（6 条）

> **共同前因**：Qwen 对"敲低/耗竭/缺失 → 活性上升"这类 **loss-of-function 证据**，未反推"因子存在时 → 活性下降"，机械套用"证据须逐字支持"而误判 unsupported。经人工对照原始证据，claim 逻辑均成立。

#### #86 ｜ H3K27me3 → decreases（cis-regulatory element）
- 来源：PMID 40186025 ｜ chain 自洽 ✓
- **原始 claim**：文献证据表明，在未指定细胞系的实验中，因子 H3K27me3 降低（decreases）cis-regulatory element 活性；机制上，H3K27 methylation is selectively depleted from primed GBP genes, facilitating hyperactivation.
- **证据原文**：H3K27 methylation is selectively depleted from primed GBP genes, facilitating hyperactivation（即"甲基化被耗竭 → 基因超激活"）。
- **前因（为何判 unsupported）**：Qwen 认为证据描述的是"耗竭后的效应"，未直接写"存在时抑制"。
- **后果（处置）**：**保留**。反推逻辑：甲基化耗竭→超激活 ⇒ 存在 H3K27me3 时→抑制，claim 成立。

#### #109 ｜ BCL11A → no_effect（promoter）
- 来源：PMID 39393056 ｜ chain 自洽 ✓
- **原始 claim**：因子 BCL11A 对 promoter 活性无显著影响（no_effect）；机制上，BCL11A does not bind to the γ-globin promoters at statistically significant frequency, while TR4 does.
- **证据原文**：BCL11A does not bind to the γ-globin promoters at statistically significant frequency（不结合 = 无效应）。
- **前因**：Qwen 将"does not bind"理解成"证据未提及"，而非"结论即 no_effect"。
- **后果（处置）**：**保留**。证据本身即支持 no_effect。

#### #137 ｜ ETO2 → decreases（enhancer）
- 来源：PMID 32960220 ｜ chain 自洽 ✓
- **原始 claim**：因子 ETO2 降低（decreases）enhancer 活性；机制上，Loss of ETO2 leads to increased enhancer activity and expression of target genes.
- **证据原文**：Loss of ETO2 → increased enhancer activity（敲除→活性↑）。
- **前因**：同 A 类共性——loss-of-function 反推未转过弯。
- **后果（处置）**：**保留**。ETO2 存在→降低，成立。

#### #286 ｜ ERRα → decreases（STING 转录）
- 来源：PMID 41708589 ｜ chain 自洽 ✗（仅该条 chain 标记问题，建议顺带补 chain 一句"ERRα 存在→抑制 STING 转录"）
- **原始 claim**：因子 ERRα 降低（decreases）sting gene transcription 活性；机制上，ERRα depletion leads to increased STING gene transcription.
- **证据原文**：ERRα depletion → increased STING gene transcription（耗竭→转录↑）。
- **前因**：同 loss-of-function 反推问题；chain 中未补反推句导致不自洽。
- **后果（处置）**：**保留**；顺带在 reasoning_chain 补充"ERRα 存在→抑制 STING 转录"以消除 chain 不自洽。

#### #272 ｜ SP3 → required_for（promoter）
- 来源：PMID 41456103 ｜ chain 自洽 ✗
- **原始 claim**：因子 Sp3 对 promoter 活性是必需的（required_for）；机制上，Sp3 transcription factor binds to the GC box 4 most involved in coronin-1 expression.
- **证据原文**：Sp3 binds to GC box 4 most involved in coronin-1 expression（结合关键 GC box 参与表达）。
- **前因**：Qwen 认为"结合"≠"必需"（未做敲除验证）；但该 GC box 是表达必需元件。
- **后果（处置）**：**保留**；chain 补"Sp3 结合必需 GC box → 支持 required_for"。

#### #418 ｜ DNA methylation → decreases（enhancer）
- 来源：PMID 33637692（K562） ｜ chain 自洽 ✗
- **原始 claim**：因子 DNA methylation 降低（decreases）enhancer 活性；机制上，DNA methylation inhibited MYB expression.
- **证据原文**：DNA methylation inhibited MYB expression（甲基化→抑制 MYB 表达）。
- **前因**：Qwen 认为证据只到"抑制基因表达"，未直接到"降低 enhancer 活性"。
- **后果（处置）**：**保留**。链条可接：甲基化→抑制 MYB 表达（MYB 为红系关键 TF）→降低其靶 enhancer 活性；chain 补全该桥接句。

---

### B 类：方向/主客体反 → 修正（5 条）｜ error pattern：**incorrect extraction（causal direction）**

> **共同前因**：Module 1 抽取时把文献中"被动/比较"句式的主客体或方向写反；这类错误是事实性错误，**必须修正**。
> **v2 修正要点**：不止修 direction——**必须按目标 schema 重建 relation 粒度**。证据只到 expression 就写 expression，不跨到 promoter activity。
> **附加价值**：这 5 条是高价值 error pattern，应收集为 **Module 1 抽取模型的 hard negative / evaluation case**。

#### #11 ｜ enhancer ↔ promoter 主客体颠倒（required_for）
- 来源：PMID 41617689 ｜ chain 自洽 ✓
- **原始 claim（错）**：因子 enhancer 对 promoter 活性是必需的（required_for）；机制上，Promoter activity may be necessary but not sufficient for enhancer function.
- **证据原文**：Promoter activity may be necessary but not sufficient for enhancer function（**promoter 对 enhancer** 是必要的）。
- **前因**：抽取时将被动句主语（promoter）与宾语（enhancer）颠倒。
- **后果（处置）**：**修正**主客体对调 + schema 化：
  ```
  factor: promoter
  relation: required_for
  target: enhancer function
  target_level: enhancer activity
  evidence_type: necessity analysis
  direction: positive
  ```

#### #236 ｜ DNMT1 方向反（increases → decreases）
- 来源：PMID 36810730（SK-N-SH） ｜ chain 自洽 ✓
- **原始 claim（错）**：因子 DNMT1 增强（increases）promoter 活性；机制上，Activated ERK enhanced DNMT1-mediated hypermethylation of FKBP52 promoter, reducing GR-mediated...
- **证据原文**：hypermethylation of FKBP52 promoter, **reducing**...（甲基化→**抑制**）。
- **前因**：抽取把"reducing"（抑制）误写成 increases。
- **后果（处置）**：**修正** effect → decreases，且 **target 限定 FKBP52 promoter**（证据直连），不泛化所有 promoter：
  ```
  factor: DNMT1
  relation: negative_regulator（经启动子高甲基化）
  target: FKBP52 promoter activity
  target_level: promoter activity
  evidence_type: hypermethylation assay
  direction: decrease
  ```

#### #250 ｜ NFIB 方向反 + 细胞系错（required_for → decreases；SK-N-SH → SH-SY5Y）
- 来源：PMID 35151899 ｜ chain 自洽 ✗
- **原始 claim（错）**：因子 NFIB 对 cdon 活性是必需的（required_for）；机制上，Silencing NFIB leads to upregulation of CDON in SH-SY5Y cells.
- **证据原文**：Silencing NFIB → upregulation of CDON（敲低→CDON **上调**），且证据细胞为 SH-SY5Y。
- **前因**：方向反（敲低→上调 ⇒ NFIB 存在→抑制，非必需）；细胞系抽取写错。
- **后果（处置）**：**修正**方向 + 细胞系 + **relation 粒度（证据在 expression 层，不写 promoter activity）**：
  ```
  factor: NFIB
  relation: negative_regulator
  target: CDON expression
  target_level: expression
  evidence_type: knockdown
  direction: decrease
  cell_line: SH-SY5Y
  ```

#### #258 ｜ HDAC2 方向反（required_for → decreases）
- 来源：PMID 36096242（SK-N-SH） ｜ chain 自洽 ✓
- **原始 claim（错）**：因子 HDAC2 对 promoter 活性是必需的（required_for）；机制上，Knockdown of HDAC2 significantly increased hSVCT2 functional expression.
- **证据原文**：Knockdown of HDAC2 → increased hSVCT2 expression（敲低→表达↑）。
- **前因**：方向反——证据表明 HDAC2 是**抑制**因子。
- **后果（处置）**：**修正**——**不写 "HDAC2 decreases promoter activity"**，而是按 schema：
  ```
  factor: HDAC2
  relation: negative_regulator
  target: hSVCT2 expression
  target_level: expression
  evidence_type: knockdown
  direction: decrease
  ```

#### #399 ｜ Ghd7/Ehd1 主客体反（modulates）
- 来源：PMID 37722508 ｜ chain 自洽 ✓
- **原始 claim（错）**：因子 Ghd7 and Ehd1 调节（modulates）promoter 活性；机制上，OsDDE9 regulates the expression of Ghd7 and Ehd1.
- **证据原文**：OsDDE9 regulates the expression of Ghd7 and Ehd1（**OsDDE9 调节 Ghd7/Ehd1**，主客体写反）。
- **前因**：把"被调节的基因"当成"调节因子"。
- **后果（处置）**：**修正**主客体对调，且 **relation 落在 expression 层**（OsDDE9 不是直接作用于 promoter）：
  ```
  factor: OsDDE9
  relation: regulates expression（上游）
  target: Ghd7 / Ehd1 expression
  target_level: expression
  evidence_type: expression analysis
  direction: modulate
  ```

**B 类总结**：5 条全部进入人工修正集；修正后入库，同时复制到 Module 1 的 hard-negative 评测集。

---

### C 类：实验体系错 → 修正（3 条）｜ error pattern：**incorrect extraction（context grounding）**

> **共同前因**：Module 1 把证据来源组织/细胞系张冠李戴。v2 修正不只改 cell_line，还要区分 **experimental system mismatch** 与 relation 强度问题。

#### #224 ｜ OGT（SH-SY5Y → SK-N-SH）★★★★
- 来源：PMID 38192280 ｜ chain 自洽 ✓
- **原始 claim（错）**：在 SK-N-SH（神经母细胞瘤细胞）实验中，OGT 降低 promoter 活性；机制上，Knock-down of OGT in **SH-SY5Y** increases ATF4 protein and mRNA expression...
- **证据原文**：实验在 SH-SY5Y 中进行；且证据只到 **ATF4 expression**（无 luciferase reporter / promoter assay / ChIP）。
- **前因**：cell_line 错 + relation 过强（expression 证据写成 promoter activity）。
- **后果（处置）**：**修正** cell_line → SH-SY5Y，relation 从 promoter activity **降为 expression regulation**：
  ```
  factor: OGT
  relation: represses
  target: ATF4 expression
  target_level: expression
  evidence_type: knockdown
  direction: decrease
  cell_line: SH-SY5Y
  ```

#### #226 ｜ ATF4（小鼠脑组织 → SK-N-SH）★★★★★
- 来源：PMID 38192280 ｜ chain 自洽 ✓
- **原始 claim（错）**：在 SK-N-SH 实验中，ATF4 对 promoter 活性是必需的；机制上，ATF4 occupancy increases at the ATF5 promoter site in **brains isolated from TMG treated mice**...
- **证据原文**：小鼠脑组织（体内）实验；且 **occupancy ↑ 只能说明 ATF4 结合增加**，不能推出 required_for promoter activity（除非有 ATF4 KO / reporter / perturbation）。
- **前因**：两个问题——① in vitro ≠ in vivo（system 错，非单纯 cell_line 错）；② causal strength 过强。
- **后果（处置）**：**修正** system 字段 + 降低因果强度：
  ```
  system_type: in vivo tissue
  organism: mouse
  tissue: brain
  factor: ATF4
  relation: binds（occupancy 增加）
  target: ATF5 promoter
  target_level: binding occupancy
  evidence_type: ChIP（occupancy）
  direction: increase（occupancy）
  ```

#### #254 ｜ YY1（SH-SY5Y → SK-N-SH）★★
- 来源：PMID 36096242 ｜ chain 自洽 ✓
- **原始 claim（错）**：在 SK-N-SH 实验中，YY1 对 promoter 活性是必需的；机制上，YY1 drives the SLC23A2 promoter activity... in VPA-treated **SH-SY5Y** cells.
- **证据原文**：实验在 SH-SY5Y 中进行；"drives SLC23A2 promoter activity" 有实验支持（relation 合理）。
- **前因**：纯 cell_line 写错，relation 无需改。
- **后果（处置）**：**修正** cell_line → SH-SY5Y，其余保留。

**C 类总结**

| ID | 问题 | 严重程度 |
|---|---|---|
| #224 | cell line 错 + expression/activity 混淆 | ★★★★ |
| #226 | cell line 错 + in vivo/in vitro 混淆 + causal 过强 | ★★★★★ |
| #254 | 纯 cell line 错 | ★★ |

---

### D 类：过度泛化 → 降级/限定（6 条）｜ error pattern：**unsupported specificity / overclaim**

> **共同前因**：Module 1 把"单基因/单位点/特定 context"的证据泛化为普遍规律。**v2 区分两类问题**：①关系错误（relation 过强，如 association→regulation）；②范围限定问题（scope 过大）。处理是降低 claim 强度，不是删除。

#### #14 ｜ EGR-1 binding motif（限定为 CSF1R 基因座）
- 来源：PMID 41538316 ｜ chain 自洽 ✓
- **原始 claim**：EGR-1 binding motif 对 enhancer-promoter interaction 活性是必需的；机制上，A single EGR-1 binding motif dictates activation of CSF1R.
- **证据原文**：单 motif 决定 **CSF1R** 的激活（specific locus）。
- **前因**：单基因证据泛化为普遍 enhancer-promoter interaction 必需。
- **后果（处置）**：**限定**为 CSF1R 基因座；且"enhancer-promoter interaction"可能过强——**若文献无 Hi-C/3C 证据，不写 interaction**，更稳写 "CSF1R enhancer activity"：
  ```
  EGR-1 binding motif is required for CSF1R enhancer activation
  ```

#### #28 ｜ targeted activation of this regulatory element（increases）
- 来源：PMID 40991835 ｜ chain 自洽 ✓
- **原始 claim**：targeted activation of this regulatory element 增强 cis-regulatory element 活性；机制上，The targeted activation of this regulatory element promotes differentiation of these aggressive...
- **证据原文**：靶向激活该调控元件→促 AML 分化。
- **前因**：把"activation（干预）"当成"measured enhancer activity"——**activation intervention ≠ activity measurement**（若实验是 CRISPRa，activation 是干预手段，不是活性读数）。
- **后果（处置）**：**改写**为干预描述，不写成"增强 CRE 活性"：
  ```
  Activation of this regulatory element promotes AML differentiation
  ```

#### #34 ｜ FOXP4（increases）
- 来源：PMID 40228145 ｜ chain 自洽 ✓
- **原始 claim**：FOXP4 增强（increases）promoter 活性；机制上，FOXP4 enhanced ESCC susceptibility and tumor growth by transcriptionally activating CYP26B1 and MYC.
- **证据原文**：FOXP4 转录激活 **CYP26B1 和 MYC** 基因。
- **前因**：transcriptionally activating gene ≠ 实验测 promoter activity。
- **后果（处置）**：**降级**为转录激活描述：
  ```
  FOXP4 activates CYP26B1 and MYC transcription
  ```

#### #43 ｜ SPI1/GABPB/STAT3（modulates）
- 来源：PMID 39870618 ｜ chain 自洽 ✓
- **原始 claim**：SPI1, GABPB and STAT3 调节（modulates）enhancer 活性；机制上，Many detected QTLs influence the predicted binding of myeloid transcription factors, including SPI1, GABPB and STAT3.
- **证据原文**：QTL 影响这几个 TF 的**结合预测**（遗传关联，非直接功能实验）。
- **前因**：**association ≠ mechanism**——variant→predicted binding change 不能推出 TF→enhancer activity。
- **后果（处置）**：**降级**改写为 variant 层描述（甚至不应放进 factor→enhancer，而应是 variant→TF binding motif）：
  ```
  Genetic variants may affect predicted binding sites of SPI1/GABPB/STAT3
  ```

#### #470 ｜ transcription factor binding motifs（modulates）
- 来源：PMID 34038741 ｜ chain 自洽 ✓
- **原始 claim**：TF binding motifs 调节（modulates）promoter 活性；机制上，caQTL variants frequently d[isrupt]...
- **证据原文**：caQTL 变异频繁破坏 TF 结合位点（variant-level 证据）。
- **前因**：位点破坏→泛化为"motif 普遍调节 promoter"（variant-level ≠ motif-level universal rule）。
- **后果（处置）**：**限定**为 variant-level：
  ```
  caQTL variants disrupting TF binding motifs can alter promoter accessibility/activity
  ```

#### #479 ｜ NF-Y（required_for，MAFA → 所有 promoter）⚠️ 需额外确认
- 来源：PMID 39189982 ｜ chain 自洽 ✓
- **原始 claim**：NF-Y 对 promoter 活性是必需的；机制上，Two E-box elements and a CCAAT motif, which bind NeuroD1 and [NF-Y]...
- **证据原文**：NF-Y 结合 **MAFA 启动子**的 CCAAT motif（单启动子证据）。
- **前因**：**binding ≠ required_for**——若原文只有 binding，不能写 required_for；需要 motif mutation / NF-Y knockdown / reporter 佐证。
- **后果（处置）**：按实验类型二选一：
  - 仅 binding：`NF-Y binds MAFA promoter`
  - 有 depletion→活性下降：`NF-Y required_for MAFA promoter activity`（限定 MAFA）

**D 类总结**

| ID | 判断 | 处理 |
|---|---|---|
| #14 | 正确 | 限定 CSF1R，避免 interaction 过强 |
| #28 | 正确但需调整 | activation intervention ≠ activity measurement |
| #34 | 正确 | expression activation，不等同 promoter activity |
| #43 | 正确 | association 不能升级为 regulation |
| #470 | 正确 | variant-level 限定 |
| #479 | 需要额外确认 | binding ≠ required_for |

---

### E 类：证据间接 → 改变 relation 类型（8 条）｜ error pattern：**evidence gap（relation 过强）**

> **共同前因**：文献支持 A 影响 B（phenotype/expression），Module 1 抽成 A 调节 promoter/enhancer activity——典型 **phenotype/expression ≠ promoter activity**。
> **v2 处理原则**：不是简单"降级"，而是**改变 relation 类型**，使 relation 与证据层级对齐：`required_for ↓ regulates expression ↓ associated with`。

#### #54 ｜ LEF1（required_for → gene regulation）
- 来源：PMID 40918098 ｜ chain 自洽 ✓
- **原始 claim**：LEF1 对 promoter 活性是必需的；机制上，LEF1 knockdown enhanced inflammatory responses and ROS production in vitro...
- **证据原文**：敲低 LEF1→炎症/ROS 上升（LEF1↓→inflammation↑，支持 LEF1 抑制炎症）。
- **前因**：把"抑制炎症反应"当"对 promoter 必需"（promoter 层面无证据）。
- **后果（处置）**：**改为 gene regulation**——**不要加"经 β-catenin/TCF 通路"**（原文无直接证据不补）：
  ```
  factor: LEF1
  relation: negatively_regulates
  target: inflammatory response-related genes
  target_level: gene expression
  evidence_type: knockdown
  direction: decrease
  ```

#### #160 ｜ DNMT3B（decreases → DNMT3B expression regulation）
- 来源：PMID 28220825（K562） ｜ chain 自洽 ✓
- **原始 claim**：DNMT3B 降低（decreases）promoter 活性；机制上，N1IC and Ets-1 suppressed the DNMT3B level in K562 cells.
- **证据原文**：N1IC/Ets-1 抑制 DNMT3B **水平**（证据只到 DNMT3B expression）。
- **前因**：原证据是 N1IC/Ets1→DNMT3B expression，不是 DNMT3B→promoter。**不要再加一层"DNMT3B→methylation→promoter repression"**（除非文献有 methylation assay）。
- **后果（处置）**：**改为 expression regulation**（保守）：
  ```
  factor: N1IC / Ets-1
  relation: regulates
  target: DNMT3B expression
  target_level: expression
  evidence_type: expression assay
  direction: decrease
  ```
  若确需保留 DNMT3B 功能：`DNMT3B may regulate promoter methylation`（标 indirect）。

#### #170 ｜ OGA（decreases → promoter interaction）
- 来源：PMID 27231347 ｜ chain 自洽 ✓
- **原始 claim**：OGA 降低（decreases）promoter 活性；机制上，Mutation of the GATA site to GAGA significantly reduces OGT and OGA promoter interactions in β-YAC...
- **证据原文**：GATA 位点突变→OGT/OGA 与 promoter 互作减少（测的是 **interaction**）。
- **前因**：interaction ≠ activity（binding → interaction → activity 三级不能跨越）。
- **后果（处置）**：**改为 promoter interaction**：
  ```
  factor: GATA motif
  relation: affects
  target: OGT/OGA promoter interaction
  target_level: promoter interaction（非 activity）
  evidence_type: interaction assay
  direction: modulate
  ```

#### #199 ｜ ALKBH5（decreases → transcription regulation）
- 来源：PMID 38227578（HepG2） ｜ chain 自洽 ✓
- **原始 claim**：ALKBH5 降低（decreases）promoter 活性；机制上，Silencing the demethylase reduced the level of BCP-RNAs and host gene transcripts.
- **证据原文**：沉默去甲基酶→BCP-RNA/宿主基因转录本下降（测转录本）。
- **前因**：转录本下降→跳跃到"promoter 活性"；且"去甲基化"来自 ALKBH5 身份，不一定来自该实验。
- **后果（处置）**：**改为 transcription regulation**（不提"去甲基化"）：
  ```
  factor: ALKBH5
  relation: positively_regulates
  target: BCP-RNA / host gene transcription
  target_level: transcription
  evidence_type: silencing
  direction: increase
  ```

#### #292 ｜ DNA methylation（decreases → 去 direction）⚠️ 最特殊
- 来源：PMID 41141393/28669403 ｜ chain 自洽 ✓
- **原始 claim**：DNA methylation 降低（decreases）promoter 活性；机制上，Hypermethylation at promoters of up-regulated genes is observed in C/G fusion-positive acute mega...
- **证据原文**：**上调基因**的启动子反而高甲基化（方向与经典知识相反）。
- **前因**：经典知识 methylation→↓transcription，但此处观察 methylation↑ 与 expression↑ 共存——可能原因：context-specific methylation / enhancer methylation / non-CpG methylation / 时间动态 / 相关而非因果。
- **后果（处置）**：**去掉 direction，改写为共现描述**（不强行给方向）：
  ```
  promoter hypermethylation co-occurs with gene upregulation
  ```

#### #323 ｜ chromatin accessibility（modulates → correlation）
- 来源：PMID 38281519 ｜ chain 自洽 ✓
- **原始 claim**：chromatin accessibility 调节（modulates）promoter 活性；机制上，ATAC-seq peaks with high chromatin accessibility located in both prom[oters]...
- **证据原文**：ATAC-seq 观察到高可及性位点（观察性，非因果）。
- **前因**：correlation ≠ causation。
- **后果（处置）**：**改为 correlation**（若保留 modulates 则 confidence 必须下降）：
  ```
  Chromatin accessibility correlates with promoter activity
  ```

#### #375 ｜ TSKU（modulates，unclear）→ 建议删除 enhancer claim
- 来源：PMID 40301309（HepG2） ｜ chain 自洽 ✓
- **原始 claim**：TSKU 调节（modulates）enhancer 活性；机制上，TSKU influences expression of vitamin D-metabolizing enzymes.
- **证据原文**：TSKU 影响维生素 D 代谢酶表达——**enhancer 完全未出现**。
- **前因**：证据只到 gene expression；无 enhancer 信息。
- **后果（处置）**：**删除 enhancer claim**（若知识库目标是 promoter/enhancer regulation，保留"TSKU regulates vitamin D enzyme expression"无 enhancer 信息，不进入知识库）。

#### #381 ｜ ODC1（increases → 删除 promoter claim）
- 来源：PMID 41359169 ｜ chain 自洽 ✗
- **原始 claim**：ODC1 增强（increases）promoter 活性；机制上，ODC1 was consistently upregulated in primary and metastatic tumors.
- **证据原文**：ODC1 在肿瘤中上调（关联，未测 promoter；甚至不能说明 ODC1 是 regulator）。
- **前因**：表达变化→跳跃到"ODC1→promoter activity"。
- **后果（处置）**：**删除 promoter claim**，降为关联描述：
  ```
  ODC1 expression is elevated in tumors
  ```

**E 类总结**

| ID | 建议 |
|---|---|
| #54 | 改为 gene regulation |
| #160 | 改为 DNMT3B expression regulation |
| #170 | 改 promoter interaction |
| #199 | 改 transcription regulation |
| #292 | 去掉 direction |
| #323 | correlation |
| #375 | 建议删除 enhancer claim |
| #381 | 删除 promoter claim |

---

### F 类：证据不相关 → 删除（2 条）｜ error pattern：**hallucinated biological relation**

> **共同前因**：Module 1 把与 promoter/enhancer 完全无关的文献结论硬套成增强子断言（variant phenotype / 其他通路）。**直接删除**，保留会污染知识库。

#### #193 ｜ rs2294510（increases）
- 来源：PMID 38365720（HepG2） ｜ chain 自洽 ✗
- **原始 claim**：rs2294510 增强（increases）promoter 活性；机制上，Mutant HepG2 cells exhibited increased lipid accumulation, a hallmark of NAFLD.
- **证据原文**：突变细胞脂质积累↑（NAFLD 表型）——variant phenotype ≠ promoter activity，完全不同层级。
- **后果（处置）**：**删除**。

#### #367 ｜ NCOA4（increases）
- 来源：PMID 40758996 ｜ chain 自洽 ✓
- **原始 claim**：NCOA4 增强（increases）promoter 活性；机制上，NCOA4 expression modulated sensitivity to ferroptosis-inducing agents.
- **证据原文**：NCOA4 调节铁死亡敏感性——无 promoter / enhancer / transcription regulation 任何关联。
- **后果（处置）**：**删除**。

---

### G 类：推理链缺陷 → 修正 chain + 降低 design 强度（4 条）｜ error pattern：**reasoning chain risk**

> **⚠️ G 类风险最高**：design_rule 是二阶推理 `observation → mechanism → design principle`，比一阶 mechanistic claim 更容易过推。
> **v2 处理原则**：chain 修正 + **降低 design 强度**；并为 design_rule 增加元数据：
> ```json
> { "rule_confidence": "low", "requires_validation": true }
> ```
> 因为 scientific discovery：correlation ≠ causation ≠ design principle。

#### #57 ｜ DNA methylation（chain 方向自相矛盾 + design 过强）
- 来源：PMID 40918098 ｜ chain 自洽 ✗ ｜ unsupported
- **原始 claim**：设计建议（文献增强）：应**避免**引入 DNA methylation 相关特征，因其与高活性**负相关**（decreases）。
- **chain 问题**：chain 中"原始发现"写 increases，但建议是"避免（负相关）"——方向自相矛盾；且"avoid"已是工程推断（文献只支持负相关）。
- **后果（处置）**：**修正 chain** 统一方向 + **降低 design 强度**（不绝对 avoid）：
  ```
  DNA methylation-associated features should be considered carefully
  due to negative association with activity
  ```

#### #425 ｜ H3 acetylation K4/K9（chain 未明确 + avoid 过强）
- 来源：PMID 28216155（K562） ｜ chain 自洽 ✗ ｜ unsupported
- **原始 claim**：应**避免**引入 H3 acetylation at lysine 4 and 9 相关特征，因其与高活性负相关（decreases）。
- **chain 问题**：chain 未明确"应避免引入"，且末步"需功能实验"削弱结论；从 mark correlation 到 design avoid 跨度大。
- **后果（处置）**：**修正 chain**——保留 "associated with activity" 表述，**避免绝对 avoid**；功能实验建议移入备注。

#### #447 ｜ acetylating and de-methylating agents（chain 缺步）⚠️ 高度可疑
- 来源：PMID 27625068（K562） ｜ chain 自洽 ✗ ｜ unclear
- **原始 claim**：应**包含** acetylating and de-methylating agents 相关特征，因其与高活性正相关（increases）。
- **chain 问题**：chain 第二步缺乏"处理与活性关系"的直接支持。
- **更深层问题**：**agent treatment（药物干预）不是 DNA sequence design feature**——若系统产出是 DNA 序列设计规则，这条原则上不成立。
- **后果（处置）**：**评估删除**（若设计对象允许 epigenetic intervention 才保留）；否则删除并记录为 design-rule 反例。

#### #449 ｜ acetylating and de-methylating agents（chain 缺步）⚠️ 高度可疑
- 来源：PMID 27625068 ｜ chain 自洽 ✗ ｜ unclear
- **原始 claim**：应包含 acetylating and de-methylating agents 相关特征，因其与高活性正相关（increases）。
- **chain 问题**：第二、三步均缺乏直接证据支持（同 #447，且来源组织未指定）。
- **后果（处置）**：**同 #447**：评估删除；或与 #447 合并去重（同一 PMID 同一结论）。

**G 类总结**

| ID | 判断 |
|---|---|
| #57 | chain 修，同时降低 design 强度 |
| #425 | chain 修，避免绝对 avoid |
| #447 | 高度可疑，可能删除（agent treatment ≠ DNA sequence feature） |
| #449 | 同 #447 |

---

## 五、处置执行计划（v2：按 error pattern 分流）

| 步骤 | 操作 | 涉及 | error pattern | 预期效果 |
|---|---|---|---|---|
| 1 | B 类 5 条 **schema 化修正**（方向/主客体 + relation 粒度）→ 同时收集为 Module 1 hard negative | #11 #236 #250 #258 #399 | incorrect extraction | 消除事实性错误 |
| 2 | C 类 3 条修正（实验体系 + relation 强度） | #224 #226 #254 | incorrect extraction | 消除 context grounding 错 |
| 3 | D 类 6 条降级/限定范围 | #14 #28 #34 #43 #470 #479 | unsupported specificity / overclaim | 消除过泛化 |
| 4 | E 类 8 条**改变 relation 类型**（与证据层级对齐） | #54 #160 #170 #199 #292 #323 #375 #381 | evidence gap | relation 不再过强 |
| 5 | F 类 2 条删除 | #193 #367 | hallucinated biological relation | 知识库 -2 |
| 6 | G 类：chain 修正 + 降低 design 强度；#447/#449 评估删除；design_rule 增加 `rule_confidence`/`requires_validation` | #57 #425 #447 #449 | reasoning chain risk | design rule 风险可控 |
| 7 | 重跑 `validate_knowledge_base.py` | 全部 | — | supported 率提升 + chain 自洽 |

**标签体系（供数据集标注复用）**：
- B/C → `incorrect extraction`
- D → `unsupported specificity / overclaim`
- E → `evidence gap`
- F → `hallucinated biological relation`
- G → `reasoning chain risk`

**最终知识库预期**：删除 F 2 条 + 可能删除 E #375、G #447/#449（3 条）→ 约 **626 条**；平均置信度保持 ≥0.8；全部 chain 自洽（或 design_rule 标注 requires_validation）；B/C/D/E 修正后的 schema 全部与证据层级对齐，可直接供 Module 2 模板填充。
