# 知识整理模块（module_knowledge_curation）构建详细报告

> 覆盖范围：claim 版本知识库从**数据来源 → 构建方法 → 数据结构 → 清洗筛选 → 验证修正**的全流程。
> 数据基准：`data/processed/knowledge_base_clean.jsonl`（最终 **617 条**）
> 生成日期：2026-08-06

---

## 〇、模块定位

本模块是 BioDesign-Agent 的知识底座：把"**什么因素/机制影响调控元件（enhancer/promoter）活性**"的原始证据（大规模统计 + 文献）整理为**统一 schema 的科学结论（claim）**，经过清洗、筛选、三层验证与两轮人工修正，产出可直接指导序列设计的高置信知识库。

```
Module 1 原始断言 ──► 去重/清洗 ──► 验证 ──► 人工复核 ──► Module 2 模板填充
   (两路数据源)      (curate)      (validate)   (v2/v3 修正)   (序列设计)
```

---

## 一、Claim 构建来源（两个数据源）

### Source 1：统计知识源 —— MPRA 高通量实验

| 环节 | 脚本 | 产出 |
|---|---|---|
| 数据本体 | Gosai et al. 2022 MPRA（K562/HepG2/SK-N-SH 三系，200bp 序列，高/低活性） | `data/raw/` MPRA h5 |
| 高低活性分组 | `scripts/analysis/split_mpra_groups.py` | `high/low_<CELL>_indices.npy` |
| 基序富集 | `scripts/analysis/run_motif_enrichment.py` + `src/knowledge/motif_enrichment.py`：JASPAR2024 白名单基序（按细胞系分组）FIMO 扫描（p<1e-4）→ Fisher 精确检验 → BH-FDR，筛选 enriched/depleted（FDR<0.01, OR>1.5） | `motif_enrichment_<CELL>.csv` |
| GC 含量 | `scripts/analysis/run_gc_analysis.py`：200bp 序列 GC 分箱 vs 活性 | `gc_summary.json` |
| DNA 形状 | `scripts/analysis/run_dna_shape_analysis.py` + `src/knowledge/dna_shape.py`：DNAshapeR 预测 MGW/ProT/Roll/HelT → Welch t-test + Cohen's d | `shape_summary.json` |

**规模**：约 **150 条**统计 claim。

### Source 2：文献知识源 —— PubMed 摘要（LLM 抽取）

| 环节 | 脚本 | 产出 |
|---|---|---|
| 定向检索 | `scripts/download/download_pubmed_activity.py`：按 6 组主题查询（增强子活性决定因素 / K562 红系 / HepG2 肝 / SK-N-SH 神经 / 染色质表观 / 序列特征），relevance 排序取 top + pmid 去重 | `data/raw/abstracts_activity.jsonl` |
| 知识抽取 | `scripts/extract/extract_activity_knowledge.py` + `src/extraction/extractor.py`（BioExtractor 封装 Qwen2.5-7B-Instruct 本地推理，`local_files_only`，left-padding，OOM 自动降批重试）| `activity_extractions.jsonl`（findings） |
| 提示词 | `src/extraction/prompts_activity.py`：严格约束 + few-shot 正反例（禁止幻觉：cell_line 未指明用 `not_specified`；effect 用受控词表；evidence_sentence 逐字引用原文） | — |

**规模**：约 **481 条**文献 claim。

---

## 二、构建方法

### 2.1 统计 claims（`scripts/build/build_stat_claims.py`）

- 从三路分析结果逐条生成中文科学结论句子，携带数值证据：
  - `evidence.p_value` / `effect_size`（OR / Cohen's d / 倍数）/ `direction`（enriched_in_high / depleted_in_high）
- **置信度由统计显著度映射**：`confidence_from_fdr`（FDR≤1e-50→0.97；≤1e-10→0.95；≤1e-3→0.85；否则 0.75）、`confidence_from_p` 同理。
- **推理链 = 方法链**：数据来源 → 分析步骤 → 显著度 → 方向解释；基序类还注入 JASPAR TF 功能注释（文献增强）。
- 每条 claim 是**完整可独立理解的句子**，非字段拼接。

### 2.2 文献 claims（`scripts/build/build_literature_claims.py`）

1. **effect 同义词归一化**到受控词表：`increases / decreases / required_for / modulates / no_effect`（enhances→increases、represses→decreases 等约 30 个同义词）。
2. **cell_line 别名归一化**：`SH-SY5Y→SK-N-SH`（亚系）、`HepG2-NTCP→HepG2`、`neuroblastoma→SK-N-SH`；无法识别回退 `not_specified`。
3. **跨文献证据合并**：同 `(cell_line, regulatory_element, factor, effect)` 跨多 PMID 时合并 `pmids + evidence_sentences + mechanisms`（不重复）。
4. **claim 类型分派**：
   - 有 mechanism → `mechanistic`；仅关联 → `association`；
   - 方向明确（increases/required_for/decreases）→ **额外生成 1 条 `design_rule`**（设计建议：正向→"应包含 X 特征，因其与高活性正相关"；负向→"应避免引入"）。
5. **因子类型启发式推断** `infer_factor_type`：正则分类 TF / motif / sequence_feature / epigenetic / other。
6. **置信度按 effect 预设**：required_for=0.92、increases/decreases=0.90、modulates=0.85、no_effect=0.80。

---

## 三、Claim 结构（统一 schema）

### 3.1 8 键 schema（唯一权威格式）

```json
{
  "claim": "完整、可独立理解的中文科学结论句子",
  "claim_type": "association | mechanistic | design_rule",
  "entities": {
    "factor": "因子（TF/基序/序列特征/表观修饰/变异）",
    "factor_type": "TF | motif | sequence_feature | epigenetic | variant | other",
    "cell_line": "K562 | HepG2 | SK-N-SH | not_specified",
    "regulatory_element": "enhancer | promoter | cis-regulatory element | ...",
    "effect": "increases | decreases | required_for | modulates | no_effect | ...",
    "pmids": ["..."],
    "context": "上下文限定（调控层级/靶基因/表型）"
  },
  "evidence": {
    "source": "PubMed 数据源 / MPRA 统计源",
    "direction": "increases_activity | decreases_activity | ...",
    "p_value": null,
    "effect_size": null
  },
  "reasoning_chain": ["来源", "证据原文", "机制描述", "因子类型判定/主客体修正"],
  "experimental_context": {
    "assay": "MPRA | literature",
    "cell_type": "...",
    "sequence_length": 200,
    "source_type": "PubMed"
  },
  "confidence": 0.0
}
```

### 3.2 relation 粒度升级（v2 教训）

初版只有 5 个 effect，修正中暴露"**expression→activity 二次过推**"问题。验证器 `EFFECT_DIRECTION_MAP` 扩展至 **16 个枚举**，强制**证据层级 ↔ relation 粒度对齐**：

| 证据层级 | 允许的 relation | 示例 |
|---|---|---|
| 活性/表达变化 | increases / decreases / required_for / modulates / no_effect | #236 DNMT1 decreases FKBP52 promoter |
| 结合/占据 | **binds** (binds_target) | #212 ASCL1 binds neuronal enhancers |
| 相关/共现 | **correlates** / **co_occurs_with** | #323 accessibility |
| 变异破坏 | **disrupts** | #470 caQTL |
| 表型促进 | **promotes** | #28 AML 分化 |
| 表达调控（方向未明） | **regulates_expression** / **positively_regulates** / **negatively_regulates** | #399 OsDDE9 |
| 影响互作/结合位点 | **affects_interaction** / **affects** (affects_binding) | #170 GATA / #43 QTLs |

**铁律**：证据只到 expression 就写 expression；只到 binding/interaction 就写 binding；只到 association 就写 correlation。**不跨层级推断**。

---

## 四、构建后的清洗和筛选

### 4.1 Curate 清洗（`scripts/curate/curate_knowledge_base.py`，633→631）

1. **去重（3 层）**：
   - 精确去重（claim 完全一致，保留最高 confidence）
   - 实体指纹去重（按记录类型分派：文献/统计 TF 基序/DNA shape/GC/其他；指纹 + PMID + 方向一致才去重，不同 PMID 是独立证据不删）
   - 命名变体去重（同 PMID + 同方向 + factor 大小写/连字符变体，如 PGC1A-En1 vs Pgc1a-En1）
2. **字段清洗**：必填字段检查（claim/claim_type/entities/evidence.source/confidence）、基因名大小写规范化（GENE_NORMALIZATION 白名单，如 Spi1→SPI1）、claim_type 一致性校验、p_value∈[0,1]。
3. **方向语义修复**（人工确认的抽取错误）：#206/207 DNA methylation（represses 却标 increases→decreases）、#474 Osa（"required to constrain"=负调控却标 required_for→decreases）。
4. **上下文标注增强**：对"方向看似矛盾实为不同调控层级"的记录补 `context`（Osa 两条、MYCN 分靶基因组 361/390/394 与 371/393/410 等）。
5. **互补关系识别（context-dependent 原则）**：同 PMID/跨 PMID 同 factor 不同效应 → 先挖底层逻辑（关系对象/测量层级不同）→ 白名单可解释则标注 context 区分 biological axis（互补非矛盾）；未覆盖则输出 `complementary_unclassified` 供人工复核。
6. **一致性验证**：同指纹方向冲突检测（排除已标注 context 的互补组）、MPRA 记录 p_value/effect_size 有效性。

### 4.2 Validate 三层验证（`scripts/curate/validate_knowledge_base.py`）

| 层 | 检查项 | 工具 | 结果 |
|---|---|---|---|
| A | **统计可靠性**：p_value∈[0,1] 非 NaN、effect_size 有效、direction 与显著性一致（富集/缺失必须 p<0.05）、推理链完整 | 规则脚本 | 150 条全过 |
| B1 | **文献一致性**：claim 与 evidence 主体/方向一致 | 规则脚本 | 481 条全过 |
| B2 | **推理链自洽性**：chain 各步骤逻辑自洽且导向 claim | Qwen2.5-7B | 逐条判定 |
| C | **外部事实核对**：supported / unsupported / unclear 三分类（mechanistic 要求证据逐字支持；design_rule 只验证前提关系+方向） | Qwen2.5-7B | 见 4.3 |

输出：`validation_report.jsonl`（逐条复核清单）+ `validation_summary.json`（汇总）。

> ⚠️ 运行要点：验证支持断点续跑，但**重跑前必须删除旧 `validation_report.jsonl`**，否则旧判定被误用（started_at==finished_at 假完成）。

### 4.3 人工复核与两轮批量修正（34 条问题记录）

初次验证（631 条）：supported 93.1%（448/481）、unsupported 31、unclear 2、chain 不自洽 1 → **34 条进入人工复核**（`data/processed/review_list_full.md` 自包含报告，逐条前因后果）。

| 类别 | 数量 | error pattern | 处置 |
|---|---|---|---|
| A. Qwen 误判-可放行 | 6 | 模型逻辑局限（loss-of-function 反推） | 保留 + 补 chain 反推句 |
| B. 方向/主客体反 | 5 | incorrect extraction (causal direction) | 修正（schema 化重建）→ **hard_negative_module1.jsonl** |
| C. 实验体系错 | 3 | incorrect extraction (context grounding) | 修正（SH-SY5Y 归一等） |
| D. 过度泛化 | 6 | unsupported specificity / overclaim | 降级/限定范围 |
| E. 证据间接 | 8 | evidence gap（relation 过强） | 改变 relation 类型（#375/#381 删除） |
| F. 证据不相关 | 2 | hallucinated biological relation | 删除 |
| G. 推理链缺陷 | 4 | reasoning chain risk | 修正 chain + 降 design 强度（#447/#449 删除→design-rule 反例） |

**fix_v2.py（631→625）**：修 25 条 + 删 6 条（#193/#367/#375/#381/#447/#449）；输出 `hard_negative_module1.jsonl`（5 条 B 类原始错误 claim + 修正对照，供 Module 1 评测）+ `design_rule_negatives.jsonl`（#447/#449）。

**fix_v3.py（625→620，用户决策：只处理可优化的，其余删除）**：
- 修 3 条：#43（QTLs modulates→affects 降级）、#212（ASCL1 required_for→binds）、#225（ATF4 chain 因果强度改写）
- 删 5 条：#86/#109/#137/#271/#285（A 类 Qwen 反推局限，人工确认成立但验证无法通过）

**最终收尾（620→617）**：
- 删 #43（Qwen 仍否）、#249（VPA 是药物处理非序列特征 + 细胞系 SK-N-SH 应为 SH-SY5Y）、#255（HDAC2 方向错误，与 #236 矛盾）
- #249 追加 `design_rule_negatives.jsonl`（→ 3 条）

---

## 五、最终结果与指标

| 指标 | 修正前 (631) | fix_v2 (625) | fix_v3 (620) | **最终 (617)** |
|---|---|---|---|---|
| 统计检查失败 | 0 | 0 | 0 | **0** |
| 文献一致性失败 | 0 | 0 | 0 | **0** |
| Qwen supported | 93.1% | 98.53% | 99.36% | **100% (467/467)** |
| unsupported | 31 | 7 | 1 | **0** |
| unclear | 2 | 0 | 2 | **0** |
| chain 一致率 | 97.9% | 99.37% | 99.57% | **100%** |

**累计处置**：删除 14 条（fix_v2 6 + fix_v3 5 + 最终 3），修正 28 条。

**输出文件清单**：

| 文件 | 内容 |
|---|---|
| `knowledge_base_clean.jsonl` | **最终知识库 617 条**（150 统计 + 467 文献） |
| `validation_report.jsonl` / `validation_summary.json` | 三层验证逐条结果与汇总 |
| `hard_negative_module1.jsonl` | 5 条 B 类因果方向错误反例（供 Module 1 评测） |
| `design_rule_negatives.jsonl` | 3 条 design-rule 过度外推反例（#447/#449/#249，agent treatment ≠ DNA feature） |
| `fix_v2_report.json` | v2 逐条处置明细 |
| `review_list_full.md` / `review_list.csv` | 34 条问题记录自包含复核报告 |

---

## 六、关键经验（工程教训）

1. **重跑验证前必须删旧报告**：`validate_knowledge_base.py` 断点续跑按 `qwen_verdict` 存在跳过，旧判定会导致假完成。
2. **relation 粒度是二次过推重灾区**：先修方向不够，必须让证据层级决定 relation 类型（expression→expression、binding→binds）。
3. **Qwen loss-of-function 反推是系统性局限**：敲低→活性上升 ⇒ 存在→抑制 这类反推，7B 模型常判 unsupported；人工确认成立后可保留或按用户决策删除。
4. **design_rule 是独立反例类别**：药物处理（agent）、非序列特征、细胞系错误等"设计建议越权"应单独收集（`design_rule_negatives.jsonl`），供后续模块规避。
5. **互补关系 ≠ 矛盾**：同 factor 不同方向多为 context-dependent（调控层级/靶基因/表型不同），先挖底层逻辑再决定标注，不轻易合并删除。
