# P0 新题型配套说明：矛盾检测 + Evidence-level 推理

> 模块三（Benchmarking）P0 阶段新增的两个评测题型，扩展基线 v2.1（见 `docs/benchmark_protocol.md`）。
> 实现：`src/data_synthesis/l3_benchmark_p0.py`（复用 v2.1 的模型加载/生成/评分/报告框架）。

---

## 1. 为什么加这两个题型

基线 v2.1 暴露的核心问题：

- **know 型 acc ≈ 0.05-0.06**（不给原文时模型无法判断方向，只会蒙 A/B 或 E）——但这只测了"事实回忆"。
- **confident_acc ≈ 0.55（不显著高于随机）**——模型"声称知道"但接近瞎猜。
- 基线只测了 **claim 内部**的方向（"X 对 enhancer 是增强还是减弱"），**没有测跨 claim 的逻辑一致性与证据层级判断**。

两个 P0 题型针对这两个空缺：

| 题型 | 测什么 | 针对基线的哪个局限 |
|------|--------|-------------------|
| **矛盾检测** | 跨 claim 逻辑一致性（两条结论是否矛盾/一致/无关） | 局限 5：claim-level 推理 → 上升到 reasoning 层面 |
| **evidence-level** | 证据层级判断（表达证据 ≠ 元件活性证据） | D/E 类错误的直接测量：模型常把"表达上调"误推为"元件活性增强" |

---

## 2. 题型 A：矛盾检测（contradiction）

### 构造逻辑

数据源：`knowledge_base_clean.jsonl`（617 条），按 `(factor, regulatory_element, cell_line)` 分组，组内按 effect 语义分三类：

- `pos`：increases / required_for（增强/必需 → 正向）
- `neg`：decreases（减弱 → 负向）
- `unknown`：modulates / binds / no_effect（方向依赖上下文/无方向）

**三类题目（分层采样保证平衡，各取约 1/3）：**

| gold | 构造方式 | 例子（真实数据） |
|------|----------|------------------|
| A 相互矛盾 | 同组内 `pos × neg` 笛卡尔积（同因子同元件方向相反） | MYCN: "增强 promoter 活性"(increases) vs "降低 promoter 活性"(decreases) |
| B 相互支持 | 同组内 pos 内两两 / neg 内两两（同向） | required_for vs increases（都是正向） |
| C 无关或无法判断 | 同组 pos/neg vs unknown + 跨组随机配对 | BRG1 required_for vs modulates（方向不明） |

### 平衡策略

- 基线 v2.1 的教训：**选项失衡 → 模型蒙多数选项得分**（v1 的 A:111 vs B:51）。
- 因此本题型**分层采样**：A/B/C 各取约 `max_pairs/3`，不足的类别全保留，合并后随机打乱。
- 真实矛盾对（A）天然稀少（同组正负方向冲突在真实文献中很少），当前知识库可构造 12 对（MYCN 组为主）——这本身就是有价值的知识库特征（**真实科学矛盾罕见**）。

### 题面示例

```
结论 1："文献证据表明……因子MYCN增强（increases）promoter活性……"
结论 2："文献证据表明……因子MYCN降低（decreases）promoter活性……"
问题：这两条结论之间是什么关系？
A. 相互矛盾  B. 相互支持（一致）  C. 无关或无法判断
```

---

## 3. 题型 B：evidence-level 推理（evidence）

### 构造逻辑

数据源：每条记录的 `reasoning_chain` 中提取"摘要原文证据"句，按证据类型分类：

| 证据类型 | 关键词 | gold |
|----------|--------|------|
| 活性类（act） | reporter / luciferase / MPRA / enhancer activity / promoter activity / 活性实验 / 报告基因 | **A 能推出** |
| 表达类（expr） | RNA-seq / mRNA / expression / transcript / 表达 / 转录 | **B 不能推出（证据层级不足）** |
| 其他（other） | 无以上关键词（机制描述等） | C 无法判断 |

结论模板：`因子{factor}调节{element}活性`（如"因子MYCN调节promoter活性"）。

### 测什么

核心是**证据层级判断**：
- 表达证据（RNA-seq 显示上调）**不能**推出"元件活性增强"——表达 ≠ 活性。
- 活性证据（luciferase/reporter/MPRA）**能**推出。

这正是 D/E 类错误（模型把表达证据误推为活性结论）的最直接测量。

---

## 4. 评分（满分 100，与 v2.1 哲学一致）

| 分项 | 分值 | 规则 |
|------|------|------|
| 答案正确性 | 60 | 答对 60；gold≠C 时答 C（无法判断）给 30 诚实分；其余 0 |
| 理由忠实性 | 40 | 引用题面关键实体 20 + 无自补编造 20 |

诚实性字段（复用于 v2.1 的 report）：`claimed`（答 A/B 视为声称知道）/ `overclaim` / 全局 honest 报告。

---

## 5. 运行方式

```bash
# 只查看题面（dry-run，不加载模型）
PYTHONPATH=src python src/data_synthesis/l3_benchmark_p0.py --dry-run

# 全量评测（并行示例，三卡）
PYTHONPATH=src python src/data_synthesis/l3_benchmark_p0.py \
  --model models/Qwen2.5-32B-AWQ --device cuda:0 \
  --out tmp/p0_eval_32b_w0.jsonl --worker 3 --worker-id 0

# 只评分
PYTHONPATH=src python src/data_synthesis/l3_benchmark_p0.py \
  --score-only tmp/p0_eval_32b_merged.jsonl --out tmp/p0_scores.jsonl

# 验证（构造 + 评分 + 诚实性）
PYTHONPATH=src python scripts/verify_l3_p0.py
```

默认参数：`--n-contradict 60 --n-evidence 60` → 共 ~112 条题面（矛盾检测 A:12/B:20/C:20 + evidence A:30/B:30）。

---

## 6. 已知局限（P0 阶段）

1. **矛盾对样本少**：gold=A 仅 12 条（真实知识库中同组正负冲突罕见），统计功效有限；后续可扩大知识库后补充。
2. **evidence 结论模板化**：`因子X调节Y活性` 是统一模板，未区分"活性增强/减弱"的细粒度；后续可扩展为方向性结论。
3. **证据提取**：取 reasoning_chain 中第一条含"摘要原文证据"的步骤，个别记录证据句可能不完整。
4. **无人类基线**：P1 计划补充 Human baseline（5 人 × 100 题）。

---

*P0 阶段完成（2026-08-07）。下一阶段：P1 Human baseline。*
