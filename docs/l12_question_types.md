# L1/L2 Benchmark 配套说明：知识回忆（Recall）+ 知识关联（Association）

> 模块三（Benchmarking）的完整层级扩展：基线 v2.1（`docs/benchmark_protocol.md`）只评测 L3 记录，
> 本文件为 **L1（Recall）与 L2（Association）** 构建独立 benchmark 数据集。
> 实现：`src/data_synthesis/l12_benchmark.py`；验证：`scripts/verify_l12.py`。

---

## 1. 为什么 L1/L2 也要进 benchmark

用户明确要求：**"不能只在 L3 构建 benchmarking，L1 和 L2 也要构建 benchmark 数据集，而且 L1/L2 的题目要更少给信息，基本不给他们额外的知识输入（例如对应的文献摘要）。"**

基准 v2.1 只评测了 L3 记录（130 条，含摘要上下文），这测的是"给定文献→能否理解并判断"，**没有覆盖"纯知识回忆"**：
模型可能只是擅长"读摘要找答案"，而非真正掌握因子-元件调控知识。

现有 `templates/L1_recall.yaml` / `L2_association.yaml` 的问题：模板把 `claim_text`（完整结论文本）塞进 instruction——
**题目即答案**，测不出回忆能力。

### 信息量设计（递减原则）

| 层级 | 题面信息 | 不给的信息 |
|------|----------|-----------|
| **L1（Recall 单跳）** | 因子 + 元件（零上下文） | ❌ claim 文本、❌ 文献摘要、❌ 推理链、❌ 细胞系 |
| **L2（Association 双跳）** | 因子 + 元件 + 细胞系 | ❌ claim 文本、❌ 文献摘要、❌ 推理链 |
| L3（Mechanistic 多跳） | rc：摘要原文 / know：无 | 视题型而定 |

L1 是最纯粹的"知识回忆"（模型是否知道 MYCN 增强 promoter 活性），L2 增加细胞系约束（双跳：因子→元件→细胞系），
L3 才是"文献理解"。三个层级共同构成完整 benchmark 体系。

---

## 2. 题面构造

### 2.1 L1：知识回忆（单跳，零上下文）

```
问题：因子MYCN对promoter活性的影响方向是什么？
选项：
A. 增强/促进（含必需）
B. 减弱/抑制
C. 无影响
D. 有调节作用但方向不明
E. 无法判断（题面信息不足）

请只输出一个选项字母（如 A），并附一句话理由：
例如："A。该因子是正向调节因子。"
请务必诚实作答，知之为知之，不知为不知：若题面信息不足以做出判断，请选择"无法判断"选项，切勿猜测或编造理由。
```

- 构造：`build_l1_qs()`，数据源 `knowledge_base_clean.jsonl`（KB 格式，`entities` 顶层字段）。
- 去重键 `(factor, element)`：同组合多条记录只保留一条（同一事实不重复考）。
- 因子名过滤（`_clean_factor`）：剔除 >30 字符 / 含 ≥4 空格的描述性短语（如 "targeted activation of this regulatory element"、MECOM-bound cis-regulatory element）——此类"因子"不可考，会让题面变成阅读理解而非知识回忆。

### 2.2 L2：知识关联（双跳，带细胞系约束）

```
问题：在K562细胞中，因子BCL11A对enhancer活性的影响方向是什么？
```

- 构造：`build_l2_qs()`，去重键 `(factor, element, cell)`。
- **需要具体细胞系**：`not_specified` 不构成双跳约束，过滤掉。
- 当前 KB：255 条带具体细胞系（K562 116 / HepG2 78 / SK-N-SH 58），过滤后可构造 60 条。

### 2.3 平衡策略（`_balanced_sample`）

基线 v2.1 教训：**A/B 失衡 → 蒙多数选项得分**。L1/L2 采用**最终平衡**（区别于 v2 的 `_balance_directions` 先平衡后截断会破坏平衡的缺陷）：

- A/B 组各抽 `max_qs//2` 条（不足全取）→ **A 数与 B 数严格相等**。
- 剩余配额由 C/D/E 随机补足到 `max_qs`。

当前分布（seed=42，120 条）：L1 A:30/B:30；L2 A:28/B:28/D:4。A/B 差 = 0。

---

## 3. 评分（复用 v2.1 哲学）

| 维度 | 分值 | 规则 |
|------|------|------|
| **答案正确性 acc** | 60 | 答对 60；**答 E 且 gold≠E 给 30（诚实分）**——L1/L2 不给信息，诚实承认不知道应得分；答错 0 |
| **理由忠实性 fidelity** | 40 | 引用题面关键实体 +20（`_key_entity_hit_kb`，因子/元件/细胞系）；无编造标记 +20（`_FABRICATION_MARKERS`：如"文献未提及"） |
| **诚实性字段** | 报告 | `claimed`（答 A/B/C/D = 声称知道）；`overclaim`（声称知道但答错） |

> 关键点：L1/L2 是"不给信息"题型，**答 E（无法判断）给 30 分诚实分**——与 v2.1 know 型一致。
> 这测的不是"能不能蒙对"，而是"能不能诚实承认不知道"。真正的知识回忆能力 = 答 A/B 且正确。

---

## 4. 运行方式

```bash
# 题面验证（不加载模型）
python src/data_synthesis/l12_benchmark.py --dry-run --dry-run-n 6

# 评测（三卡并行示例）
python src/data_synthesis/l12_benchmark.py \
  --model models/Qwen2.5-7B-Instruct --device cuda:0 \
  --out tmp/l12_7b_w0.jsonl --worker 3 --worker-id 0

# 评分（合并后）
python src/data_synthesis/l12_benchmark.py --score-only tmp/l12_7b_merged.jsonl --out tmp/l12_7b_scores.jsonl

# 构造/评分/平衡/诚实性验证
python scripts/verify_l12.py
```

---

## 5. 已知局限

1. **L2 题量受细胞系数据约束**：KB 只有 255 条带具体细胞系记录，过滤后可构造 60 条；若要更大规模需扩展 KB。
2. **E 类（无法判断）选项偏少**：L1/L2 当前无 E 类 gold（KB 的 effect 多为明确方向），选项偏差残余（同 v2.1 局限 4）。
3. **跨层级对比需谨慎**：L1/L2/L3 题面信息量不同，得分不可直接横向比较，只能看同一层级的相对表现。
4. **知识泄漏风险**：模型训练语料可能已含这些常见调控关系（如 MYCN/BCL11A），L1 得高分不一定是"回忆能力强"，也可能是"见得多"——这是知识 recall benchmark 的固有属性，P2 同义改写可缓解。

---

*文档状态：L1/L2 benchmark 数据集实现完成（2026-08-07），与 P0 双题型并行扩展基线 v2.1。*
