# BioDesign-Agent 模块一 + 模块二 总结

> 覆盖范围：**模块一（知识整理 module_knowledge_curation）** + **模块二（数据合成与评测 data_synthesis / benchmark）**
> 本文档为两个模块的最小实现总结，逐文件标明作用（临时文件、模型权重、运行日志不列入）。
> 生成日期：2026-08-07

---

## 〇、项目总览

BioDesign-Agent 的目标是：从 **PubMed 文献 + MPRA 高通量实验** 出发，构建 **基因调控序列设计 Agent** 的知识底座与可消费数据。

```
模块一：知识整理（Module 1）
  PubMed/MPRA 原始数据 ──► 分析/抽取 ──► Claim 构建 ──► 清洗 ──► 三层验证 ──► 人工修正
                                   └──────────► 知识库 knowledge_base_clean.jsonl (617 条)

模块二：数据合成 + 评测（Module 2）
  知识库 ──► 模板填充 (drafts) ──► LLM 增强 (enhanced) ──► 质量过滤 (filtered)
        ──► L3 封闭式 Benchmark (l3_benchmark_v2) ──► 7B/32B 能力评测 + 诚实性评估
```

分层定位（详见 `data_architecture.md`）：**知识库 = L0 原料层**，**训练数据/benchmark = L2 产品层**，二者由不同改造管线（管线 A / 管线 B）从同一知识库产出。

---

## 一、模块一：知识整理（module_knowledge_curation）

### 1.1 定位

把"什么因素/机制影响调控元件（enhancer/promoter）活性"的原始证据整理为**统一 schema 的科学结论（claim）**，经清洗、三层验证与两轮人工修正，产出高置信知识库。**最终产物：`knowledge_base_clean.jsonl`（617 条）**。

### 1.2 流程与文件作用

| 环节 | 文件 | 作用 |
|------|------|------|
| **① 数据下载** | `scripts/download/download_pubmed_activity.py` | 按 6 组主题（增强子活性决定因素/K562 红系/HepG2 肝/SK-N-SH 神经/染色质表观/序列特征）定向检索 PubMed，relevance 排序取 top + pmid 去重 → `abstracts_activity.jsonl` |
| | `scripts/download/download_pubmed_mpra.py` | 下载 MPRA 相关文献摘要 |
| | `scripts/download/dedupe_abstracts.py` | 摘要去重（按 pmid/title） |
| **② 序列分析（统计源）** | `scripts/analysis/split_mpra_groups.py` | MPRA 200bp 序列按高/低活性分组（三细胞系 K562/HepG2/SK-N-SH）→ `high/low_<CELL>_indices.npy` |
| | `scripts/analysis/run_motif_enrichment.py` | 基序富集主流程：JASPAR2024 白名单 FIMO 扫描 → Fisher 精确检验 → BH-FDR → `motif_enrichment_<CELL>.csv` |
| | `scripts/analysis/run_gc_analysis.py` | 序列 GC 含量分箱 vs 活性关联 → `gc_summary.json` |
| | `scripts/analysis/run_dna_shape_analysis.py` | DNAshapeR 预测 MGW/ProT/Roll/HelT → Welch t-test + Cohen's d → `shape_summary.json` |
| | `scripts/analysis/run_gc_motif_enrichment.py` | GC-motif 联合富集分析（基序 GC 偏好 × 活性）→ `gc_motif_claims_<CELL>.jsonl` |
| | `scripts/analysis/run_motif_gc_features.py` | 基序 GC 特征表 → `motif_gc_features.csv/json` |
| **③ 知识抽取（文献源）** | `scripts/extract/extract_activity_knowledge.py` | 调用本地 Qwen2.5-7B-Instruct 从摘要批量抽取结构化 findings（cell_line/element/factor/effect/mechanism/evidence_sentence）→ `activity_extractions.jsonl` |
| | `src/extraction/extractor.py` | `BioExtractor`：封装本地 LLM 推理（local_files_only、left-padding、OOM 自动降批重试、批量推理） |
| | `src/extraction/prompts_activity.py` | 抽取提示词：严格约束 + few-shot 正反例（禁止幻觉：cell_line 未指明用 `not_specified`；effect 受控词表；evidence_sentence 逐字引用原文） |
| | `src/extraction/prompts.py` | 通用抽取提示词（旧版文献抽取） |
| | `scripts/extract/extract_knowledge.py` | 通用知识抽取入口（旧版） |
| | `scripts/extract/extract_dna_shape_table.py` | DNA 形状特征表抽取 |
| | `scripts/extract/review_extractions.py` | 抽取结果人工复核工具 |
| **④ Claim 构建** | `scripts/build/build_stat_claims.py` | 统计源 claims：三路分析结果 → 中文科学结论句（携带 p_value/effect_size/direction），置信度由 FDR/p 映射，推理链=方法链 → `knowledge_base_stat.jsonl` |
| | `scripts/build/build_literature_claims.py` | 文献源 claims：effect 同义词归一化（enhances→increases 等 ~30 词）、cell_line 别名归一化（SH-SY5Y→SK-N-SH）、跨文献证据合并、claim_type 分派（mechanistic/association/design_rule）、因子类型启发式推断 → `knowledge_base_clean.jsonl` |
| **⑤ 清洗筛选** | `scripts/curate/curate_knowledge_base.py` | 三层去重（精确/实体指纹/命名变体）+ 字段清洗 + 方向语义修复 + 上下文标注 + 互补关系识别 → 633→631 |
| **⑥ 三层验证** | `scripts/curate/validate_knowledge_base.py` | A 统计可靠性（p 值/direction 一致性）+ B1 文献一致性 + B2 推理链自洽性（Qwen 判定）+ C 外部事实核对（supported/unsupported/unclear）→ `validation_report.jsonl` / `validation_summary.json` |
| **⑦ 人工修正** | `scripts/curate/fix_v2.py` | 第一轮批量修正：34 条问题记录 → 修 25 + 删 6 → `hard_negative_module1.jsonl`（5 条因果方向错误反例）、`design_rule_negatives.jsonl`、`fix_v2_report.json` |
| | `scripts/curate/fix_v3.py` | 第二轮：只处理可优化条目，修 3 删 5 |
| **⑧ 知识组件** | `src/knowledge/motif_enrichment.py` | 基序富集核心实现（FIMO + Fisher + FDR） |
| | `src/knowledge/gc_analysis.py` | GC 分箱分析实现 |
| | `src/knowledge/dna_shape.py` | DNA 形状分析实现 |
| | `src/knowledge/mpra_loader.py` | MPRA 数据加载与预处理 |
| | `src/knowledge/triplet_builder.py` | 旧版三元组构建（head/relation/tail，已废弃，保留兼容） |

### 1.3 产物清单

| 文件 | 内容 |
|------|------|
| `knowledge_base_clean.jsonl` | **最终知识库 617 条**（150 统计 + 467 文献），8 键 schema |
| `knowledge_base_stat.jsonl` | 统计源 claims（633 条） |
| `knowledge_base_gc_motif.jsonl` | GC-motif 层 claims（33 条，`gc_association`） |
| `validation_report.jsonl` / `validation_summary.json` | 三层验证逐条结果与汇总 |
| `hard_negative_module1.jsonl` | 5 条因果方向错误反例（供模块一评测） |
| `design_rule_negatives.jsonl` | 3 条 design-rule 过度外推反例 |
| `fix_v2_report.json` | v2 逐条处置明细 |
| `review_list_full.md` / `review_list.csv` | 34 条问题记录自包含复核报告 |
| `curation_report.jsonl` / `curation_summary.json` | 清洗过程明细 |

### 1.4 关键指标

| 指标 | 最终值 |
|------|--------|
| 统计检查失败 | 0 |
| 文献一致性失败 | 0 |
| Qwen supported 比例 | 100% (467/467) |
| chain 一致率 | 100% |
| 累计处置 | 删 14 条、修 28 条 |

---

## 二、模块二：数据合成 + 评测（data_synthesis / benchmark）

### 2.1 定位

从模块一知识库出发，完成两件事：
1. **管线 B（训练数据）**：模板填充 + LLM 增强 + 质量过滤 → 可微调数据集（`enhanced_full2.jsonl` 908 条，L1-L5 五级分层）。
2. **管线 A（Benchmark）**：L3 封闭式评测（v2.1），修复 v1"题面泄露答案"与 v2"蒙 A 得分"问题，支持 7B/32B 能力对比与**诚实性评估**。

### 2.2 流程与文件作用

| 环节 | 文件 | 作用 |
|------|------|------|
| **① 模板定义** | `templates/L1_recall.yaml` | L1 记忆题模板（5 个：事实回忆） |
| | `templates/L2_association.yaml` | L2 关联题模板（4 个：实体-效应关联） |
| | `templates/L3_mechanistic.yaml` | L3 机制推理模板（4 个：多跳因果链解释） |
| | `templates/L4_design.yaml` | L4 设计题模板（4 个：序列设计建议） |
| | `templates/L5_optimization.yaml` | L5 优化题模板（3 个：多目标优化/失败分析） |
| **② 模板填充** | `src/data_synthesis/template_engine.py` | `TemplateFiller`：解析 yaml 模板 → 从知识库按占位符/类型约束采样填充 → 轮询均衡 → `drafts_full2.jsonl`（192 条 L3 等） |
| **③ LLM 增强** | `src/data_synthesis/llm_enhancer.py` | `LLM Enhancer`：L1/L2 规则直出（direct 答案）+ L3/L4/L5 调用 Qwen2.5-32B-AWQ 三卡并行生成，严格禁止编造文献外事实（system prompt 约束）→ `enhanced_full2.jsonl`（908 条） |
| **④ 质量过滤** | `src/data_synthesis/quality_filter.py` | 规则校验（结构完整性：L4 须 3+ 条准则、L5 须失败原因/改进/预期效果、实体引用、PMID 一致性）+ 可选 LLM judge（判同义不同形）→ `filtered_full2.jsonl` |
| **⑤ 辅助工具** | `src/data_synthesis/gene_names.py` | 内置常用基因名表：判定 regulatory_element 前缀（如 "FKBP52 promoter"）是否基因特异 |
| | `scripts/regen_l12_l4.py` | 局部重生成：L1/L2 新直出逻辑（补核心实体）、L3/L5 保留、L4 断点续跑 → 统一 9 键 schema |
| | `scripts/sub_drafts.py` | 用更新模板 + 原 KB 记录重新填充问题记录 instruction，写临时 drafts 供小批量验证 |
| | `scripts/run_full_pipeline.sh` | 一键全量流水线：重建 drafts → 三卡并行增强 → 规则版 + judge 版质量校验 |
| **⑥ Benchmark** | `src/data_synthesis/l3_benchmark_v2.py` | **L3 封闭式评测 v2.1**（当前主线）：三题型（rc 给摘要/know 不给/data 统计）+ A-E 五选项 + 正负方向平衡（`_balance_directions`）+ 诚实指令 + 不诚实惩罚分 + 三卡并行 + 断点续跑 + score-only 评分 |
| | `src/data_synthesis/l3_benchmark.py` | L3 benchmark v1（旧版：题面泄露答案，已废弃） |
| | `scripts/verify_l3_v2.py` | 验证脚本：题面构造/解析/评分/平衡/诚实性指标 7 组测试 |

### 2.3 训练数据产物 schema（统一 9 键）

`enhanced_full2.jsonl` 每条强制统一 9 键：`id` / `level` / `template_id` / `instruction` / `input` / `output` / `raw_output_placeholder` / `direct`（true=规则直出 / false=LLM 生成）/ `metadata`（8 键：source_claim_ids / claim_type / entities_used / evidence / reasoning_chain / gene_specific / llm_model / quality_score）。

**规模**：908 条 = L1(192) + L2(200) + L3(196) + L4(192) + L5(128)。

### 2.4 L3 Benchmark v2.1 关键设计

| 特性 | 说明 |
|------|------|
| 三题型 | rc（给摘要原文判方向/判结论）/ know（不给原文纯知识）/ data（统计判读） |
| 五选项 | A 增强 / B 减弱 / C 无影响 / D 方向不明 / E 无法判断 |
| 平衡 | `_balance_directions`：gold A/B 各半，杜绝"蒙 A 得分"（v2 教训：7B 蒙 A 拿 0.596 acc → v2.1 平衡后 know 型 acc 暴跌 0.048 铁证） |
| 诚实指令 | "知之为知之，不知为不知"，信息不足选"无法判断"，禁猜测 |
| 评分 | 正确性 60（know 型答 E 给诚实分 30）+ 忠实性 40（实体引用 20 + 无编造 20） |
| 诚实性指标 | claimed / confident_acc（声称知道中答对率）/ overclaim_rate / dishonest（confident_acc<0.5 扣分 `(0.5-confident_acc)*40`，上限 20） |
| 并行 | 多 worker 分片（三卡 GPU0/2/3）+ 断点续跑 + score-only 一致性 |

**评测结果（v2.1 清洗后 130 条）**：32B avg 64.7 > 7B avg 61.7；rc 型 acc 0.556 双双持平；答 A/B 时 32B 88.6% / 7B 85.1% 准确（有真知识）；答 C/D 时均 ~0%（无"方向不明"概念）；confident_acc 0.573 / 0.545（未显著高于 50% 随机）。

### 2.5 数据质量修复（方案 1：误标清洗）

发现：22 条 modulates/binds 记录中至少 5 条摘要句含明确方向词被误标（如 "enhanced NRF2 binding"、"increased GLI3 binding"、"required to activate enhancers"、"inactive enhancers/devoid of H3K27ac"）→ 清洗脚本改判为 increases/decreases，产出 `drafts_full2_clean.jsonl`（130 条：A:59/B:59/D:10/C:2），修复后 32B/7B 评分均上升且排序稳定。

---

## 三、关键数据流（端到端）

```
data/raw/abstracts_activity.jsonl (文献)
        │  extract_activity_knowledge.py (Qwen2.5-7B)
        ▼
activity_extractions.jsonl (findings, 259 条)
        │  build_literature_claims.py
        ▼
knowledge_base_clean.jsonl (617 条, 8 键 schema)  ◄── 模块一最终产物
        │  template_engine.py (TemplateFiller)
        ▼
drafts_full2.jsonl (模板填充草稿)
        │  llm_enhancer.py (32B-AWQ 三卡)
        ▼
enhanced_full2.jsonl (908 条, 9 键 schema)  ◄── 训练数据（管线 B）
        │  quality_filter.py
        ▼
filtered_full2.jsonl (质量过滤后)
        │  l3_benchmark_v2.py (封闭式题面 + 评分)
        ▼
7B/32B 评测分数 + 诚实性报告  ◄── Benchmark（管线 A）
```

---

## 四、运行方式

```bash
# 模块一：知识库构建（含抽取/分析）
python scripts/extract/extract_activity_knowledge.py --input data/raw/abstracts_activity.jsonl
python scripts/build/build_literature_claims.py --output data/processed/knowledge_base_clean.jsonl
python scripts/curate/curate_knowledge_base.py
python scripts/curate/validate_knowledge_base.py
python scripts/curate/fix_v2.py && python scripts/curate/fix_v3.py

# 模块二：数据合成（一键流水线）
bash scripts/run_full_pipeline.sh 64 full2

# 模块二：L3 评测（三卡并行 + 评分）
python src/data_synthesis/l3_benchmark_v2.py --drafts data/synthetic/drafts_full2_clean.jsonl \
    --model models/Qwen2.5-32B-AWQ --device cuda:0 --worker 3 --worker-id 0
python src/data_synthesis/l3_benchmark_v2.py --drafts ... --score-only tmp/xxx.jsonl --out tmp/scores.jsonl

# 验证
python scripts/verify_l3_v2.py
```

---

*文档状态：模块一 + 模块二最小实现总结（2026-08-07）。临时文件（tmp/、日志）、模型权重（models/、Qwen2.5-7B-Instruct/）不列入。*
