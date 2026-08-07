# Module 2 模板库 · 占位符词汇表与字段映射规范

本目录存放 Module 2（Scientific Reasoning Data Factory）各层级指令模板（YAML）。
模板由人工设计（2-3 个示例起步），Agent 负责解析占位符并填充断言。

## 模板 YAML 格式

```yaml
- id: L1_001
  level: L1
  description: "..."
  instruction_template: |
    ...{placeholder}...
  expected_output_type: "short factual answer"
```

## 占位符词汇表（TemplateFiller 归一化依据）

### 通用占位符（所有记录类型可用）

| 占位符 | 来源字段 | 说明 |
|--------|----------|------|
| `{claim_text}` | `claim` | 完整科学结论句子 |
| `{claim_type}` | `claim_type` | mechanistic / association / design_rule |
| `{reasoning_chain}` | `reasoning_chain` | 数组 → 编号列表文本（1. ... 2. ...） |
| `{evidence_direction}` | `evidence.direction` | 16 枚举方向（increases_activity 等） |
| `{p_value}` | `evidence.p_value` | 统计记录可用 |
| `{effect_size}` | `evidence.effect_size` | 统计记录可用 |
| `{confidence}` | `confidence` | 0-1 |

### 文献类占位符（claim_type = mechanistic / design_rule，共 546 条）

| 占位符 | 来源字段 |
|--------|----------|
| `{factor}` | `entities.factor` |
| `{factor_type}` | `entities.factor_type`（TF/motif/sequence_feature/epigenetic/variant/other） |
| `{cell_line}` | `entities.cell_line`（K562/HepG2/SK-N-SH/not_specified） |
| `{regulatory_element}` | `entities.regulatory_element`（enhancer/promoter/...） |
| `{effect}` | `entities.effect`（increases/decreases/required_for/...） |
| `{pmids}` | `entities.pmids`（列表 → "PMID 1、PMID 2"） |
| `{context}` | `entities.context`（可选，有则填） |

### 统计类占位符（claim_type = association，共 71 条）

| 占位符 | 来源字段 | 记录子类 |
|--------|----------|----------|
| `{tf}` | `entities.tf` | TF 基序富集 |
| `{motif}` | `entities.motif` | TF 基序富集（JASPAR ID，如 MA0475.3） |
| `{gc_range}` | `entities.gc_range` | GC 含量 |
| `{gc_lo}` / `{gc_hi}` | `entities.gc_lo` / `entities.gc_hi` | GC 含量 |
| `{shape}` | `entities.shape` | DNA 形状（MGW/ProT/Roll/HelT） |
| `{feature}` | `entities.feature` | DNA 形状（如 MGW_mean） |
| `{direction}` | `evidence.direction` | enriched_in_high / depleted_in_high / optimal / suppressive / high / low |

## 模板→断言匹配规则（TemplateFiller 实现要点）

1. 模板中出现的**所有**占位符必须能从同一条断言解析，否则跳过该断言。
2. 含 `{motif}`/`{tf}`/`{gc_range}`/`{shape}` 的模板 → **只采样 association（统计类）记录**。
3. 含 `{factor}`/`{regulatory_element}`/`{effect}` 的模板 → 只采样文献类记录（mechanistic/design_rule）。
4. `{cell_line}` 两类记录都有；`not_specified` 的记录在需要具体细胞系的模板中应过滤。
5. 模板可声明 `require_claim_type`（如 L3 → mechanistic、L4 → design_rule），填充前过滤。
6. 多实体占位符（如 `{tf}`+`{motif}`）必须来自同一条断言（同源约束）。

## 层级-类型映射（Module 1 衔接，§7）

- L1/L2：混合类型（单/双跳事实）
- L3：优先 `mechanistic`（reasoning_chain 直接填充推理步骤）
- L4：优先 `design_rule`（设计建议 → 序列设计任务）
- L5：优先 `mechanistic` + design_rule（反馈迭代）

## 模板文件清单（当前种子，逐步扩充至 ~50 个）

| 文件 | 层级 | 模板数 | 状态 |
|------|------|--------|------|
| L1_recall.yaml | L1 知识回忆 | 3 | 种子 |
| L2_association.yaml | L2 知识关联 | 3 | 种子 |
| L3_mechanistic.yaml | L3 机制推理 | 3 | 种子 |
| L4_design.yaml | L4 序列设计 | 3 | 种子 |
| L5_optimization.yaml | L5 闭环优化 | 2 | 种子 |

> 种子模板为骨架，验证填充流程跑通后按 §6 规模规划扩充（8-12 模板/层级）。
