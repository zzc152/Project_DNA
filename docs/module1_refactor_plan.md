# 模块一（知识整理）重构蓝图

> 日期：2026-08-08
> 背景：全链路问题盘点后决定**不逐点修补，从模块 1 重构**，从源头消除问题类别。
> 本文档 = 问题清单 → 根因分析 → 模块 1 重构设计 → 实施顺序。

---

## 一、问题全景（是什么）

按数据流水线归类，全项目共 **7 类问题**：

| # | 问题 | 量化证据 | 所在环节 |
|---|------|----------|----------|
| P1 | **gold 方向错误 / 过度自信** | T1 5 条矛盾；T5 摘要↔gold 一致率仅 57.1%，13/63 条摘要信息不足仍标方向 | KB 抽取 → 标注 |
| P2 | **相关/共现被当因果** | 3 条 KB 记录（co_occurs_with/correlates 类）已删；L3_003_64 相关→required_for | KB 抽取 |
| P3 | **因子名噪声（描述性短语当因子）** | 21/189 L3 条 factor_valid=false（11.1%）；KB 29 条 factor>30 字符 | KB 抽取 → schema |
| P4 | **修饰标记当因子（类型错误）** | 24/189 L3 条（12.7%）：H3K4me3、DNA methylation 非可作用因子 | KB 抽取 → schema |
| P5 | **评测不可增量对比（题面漂移）** | final 删 3 条 → 21 条 rc↔know 翻转 + 12 条 acc 翻转 | 评测脚本 |
| P6 | **验证体系只查字段、不查语义** | T1-T5 漏掉全部 P3/P4；T5 仅模型辅助 | 验证 |
| P7 | **统计功效低 / 选项偏差 / 知识泄漏** | L3 每模板 40-60 条；D 档 5%；L1/L2 测常见调控 | 评测设计 |

---

## 二、根因分析（为什么发生）

### 根因 R1：抽取是"单次自由文本输出"，无受控约束

`prompts.py / prompts_activity.py` 让 LLM 输出自由文本 JSON：
- 因子名**无词表约束** → P3（"acetylating and de-methylating agents"）
- **无实体类型字段校验** → P4（修饰标记进 factor 槽位）
- **无因果强度维度** → P2（共现/相关被无差别提取成 required_for/increases）

### 根因 R2：schema 缺"证据→结论"的强制绑定

KB 8 键 schema 中 evidence 与 entities 是平级字段，**没有机制保证"结论必须由证据推出"**：
- 摘要只说"binds/primed"（L3_002_37）→ gold 标 required_for
- 摘要只说"correlates with"（L3_003_64）→ gold 标 required_for
- → P1、P2 的温床

### 根因 R3：出题/评测依赖数据集合顺序

`qtype = idx % 2`、`rng.sample` 都以**全量列表索引/集合**为随机源：
- 数据增删 1 条 → 题面连锁漂移 → P5
- 这不是数据问题，是**管线非确定性**问题

### 根因 R4：验证与构造脱节

T1-T5 是事后补的验证，**验证维度 = 构造时没约束的维度**：
- 构造时不限制 factor 类型 → 验证时也无从校验（T1-T5 盲区）→ P6
- 验证只做字段级（矛盾/冲突/非法值），不做语义级（实体类型/因果链）→ P3/P4 漏网

### 根因 R5：一次性构建，无版本化

数据是"跑一遍得到"，没有：
- 数据版本号 → 无法回溯"哪个版本的评测对应哪个版本的数据"
- 变更日志 → 改了什么、为什么改不可追踪
- → P7 的评测对比困难 + 信任崩塌的放大器

---

## 三、模块 1 重构设计（怎么解决）

### 3.0 核心原则（贯穿所有设计）

1. **约束前置**：能在构造时硬约束的，绝不留到验证时靠抽查
2. **证据绑定**：任何结论字段必须能从 evidence 推出，否则拒绝入库
3. **确定性**：一切随机性由稳定 ID 派生（`hash(id)`），与集合顺序无关
4. **版本化**：数据 + 代码 + 评测结果三者绑定版本，可回溯
5. **留痕**：每条记录保留全链路 provenance

### 3.1 新 KB Schema v3（解决 R1/R2）

```jsonc
{
  "id": "KB_0001",
  "claim": "BCL11A required_for enhancer activity in K562",
  "relation_kind": "causation",          // causation | co_occurrence | correlation | association
  "causal_strength": "direct_evidence",  // direct_evidence | indirect | inferential | unknown
  "entities": {
    "factor": { "name": "BCL11A", "id": "HGNC:13218", "type": "TF" },  // 受控实体
    "element": { "name": "enhancer", "id": "SO:0000165", "type": "regulatory_element" },
    "cell_line": "K562",
    "effect": "required_for"             // 枚举：increases/decreases/required_for/modulates/no_effect/unknown
  },
  "evidence": {
    "pmids": ["41191525"],
    "quotes": ["..."],                    // 原文引用（支持结论的唯一依据）
    "direction_support": "explicit"       // explicit | implicit | absent —— 摘要是否明确给了方向
  },
  "provenance": { "extractor": "qwen2.5-7b-v2", "date": "...", "review": "auto|manual", "status": "released" }
}
```

**关键变化**：
- `relation_kind` 把"因果"与"共现/相关"**在 schema 层分离** → P2 结构性消失
- `factor` 变为受控对象 `{name, id, type}` → P3/P4 结构性消失
- `evidence.direction_support` 强制标注摘要方向明确度 → 摘要不足的题 **gold 必须是 unknown**，P1 从源头消除

### 3.2 受控词表 + 实体链接（解决 P3/P4）

建立最小 ontology（优先复用现有库）：

| 实体类型 | 词表来源 | 举例 |
|----------|----------|------|
| TF / 基因 | HGNC 符号 | BCL11A、MYCN |
| 组蛋白修饰 | 规范名（H3K4me3、H3K27ac…） | 仅允许 `修饰` 类型，禁入 factor |
| 元件 | SO 术语 | enhancer、promoter、silencer |
| 细胞系 | Cellosaurus | K562、HepG2 |

**写入规则**：
- LLM 输出 → 实体链接到词表 → **解析失败则拒绝入库**（fail-closed）
- 修饰标记只能作为 `element.modification` 或独立实体，**禁止出现在 factor 槽位**

### 3.3 抽取 Prompt v2（解决 P1/P2）

`prompts_activity.py` 重构要点：
1. **强制输出受控枚举**：effect、relation_kind、direction_support 全部单选枚举，不许自由文本
2. **证据三段论**：要求输出 `claim → quote → direction_support` 三元组，quote 缺失则 direction_support=absent → effect 强制 unknown
3. **区分共现与因果**：显式指令"观察性共现（correlates/co-occurs）不得输出为因果 effect"，直接输出 relation_kind=correlation

### 3.4 流水线重构（pipeline v2）

```
extract(LLM) → validate(schema+词表+证据链) → curate(规则/人工) → review(抽检) → release(版本化)
     ↑              ↑                            ↑
  可重跑           失败即拒                  高疑虑优先
```

- 每步产物独立落盘（raw → validated → curated → reviewed → released）
- `validate` 层自动校验：枚举合法性、实体链接、证据完整性、方向支持度
- 新**增量流程**：数据变更只重新跑受影响部分，用 `hash(id)` 保证题面稳定

### 3.5 确定性评测（解决 P5）

- 所有 `rng` 改为 `random.Random(hash(record_id))` 派生
- `qtype` 改为 `hash(id) % 2`
- 评测脚本独立于集合顺序：**同一 id 永远同一题面**，数据增删不再漂移
- 评测结果落盘时带 `data_version` 字段

### 3.6 验证体系 v3（解决 P6）

分层校验，构造约束前置后验证只需兜底：

| 层 | 校验内容 | 方式 |
|----|----------|------|
| Schema 层 | 枚举/必填/类型 | 自动（构造时已做） |
| 语义层 | 实体类型合法、因果链完整、方向支持度 | 自动规则 + LLM 辅助 |
| 人工层 | 高疑虑优先抽检（随机 + 可疑样本加权） | 人工仲裁 |

新增 **T6 因子类型校验**（自动正则/词表：修饰标记、多实体列表、描述性短语 → 标记），**T7 因果链校验**（relation_kind=causation 必须有 direct_evidence 或 explicit 方向）。

### 3.7 版本化（解决 P7 的对比问题）

- 数据版本 = 内容哈希（或递增号），发布时生成 `release_manifest.json`（版本号 + 记录数 + 变更日志）
- 评测产物绑定数据版本：`eval_result_v3_on_data_v7.json`
- 对比评测只允许**同数据版本**，不同版本用配对（同 id）对比

---

## 四、重构实施顺序（最小可用优先）

| 步骤 | 动作 | 产出 | 依赖 |
|------|------|------|------|
| S1 | 抽取 Prompt v2 + schema v3 落地 | 新 KB 抽取器（约束枚举+证据绑定） | 无 |
| S2 | 受控词表 + 实体链接（fail-closed） | ontology.json + 链接模块 | S1 |
| S3 | 新 KB 全量重抽取（沿用 614 条 PMID 源） | KB_v2 | S1+S2 |
| S4 | 出题脚本确定性改造（hash(id)） | L1/L2/L3 题面稳定 | S3 |
| S5 | 验证 v3（T6/T7 + 分层） | 自动校验管线 | S3 |
| S6 | 三卡重跑评测 + 版本化产物 | 新基线（data_v2 可对比） | S4+S5 |

> S1-S2 是半天到一天工作量；S3 是重抽取（GPU 小时级）；S4-S6 是半天。**整个重构 3-5 天可完成。**

---

## 五、重构时注意事项（踩坑清单）

1. **不要重做全部内容**：旧 614 条 KB 的价值在于 PMID 来源与人工核查结论，重抽取应**复用 PMID 列表**而不是重新选文献
2. **MYCN 教训要带进 schema**：`(factor, element, cell_line)` 键仍缺靶基因维度 → schema v3 需预留 `target_gene` 或记录级 context 字段，避免生物学真实再被误判为冲突
3. **LLM 抽取要校准**：重抽取后必须用人工样本（≥50 条）测新抽取器的准确率，不达标不上线
4. **证据三段论是硬门槛**：摘要无明确方向 → gold 必须 unknown，宁可"少而准"不要"多而错"
5. **评测对比纪律**：重构后旧分数全部标注"旧版不可比"，基线以 v2 为准，禁止新旧混讲
6. **D 档（unknown/modulates）从 schema 层就保证存在**：词表/枚举里 effect=unknown 是合法值，让方向不明成为一等公民，避免后续补 D 档
7. **每个字段变化都要能解释**：release_manifest 记录变更原因，维护可审计性（这已是本项目的优势，保持）

---

## 六、重构后预期收益

| 问题 | 重构前 | 重构后 |
|------|--------|--------|
| P1 方向错误/过度自信 | 5-13 条，靠事后修 | schema 层强制：无证据即 unknown |
| P2 相关当因果 | 3+ 条已删 | relation_kind 分离，结构不可能 |
| P3 因子名噪声 | 21 条标记未清洗 | 词表 fail-closed，进不了库 |
| P4 修饰当因子 | 24 条未处理 | 类型约束，进不了 factor 槽 |
| P5 评测漂移 | 21 条题面翻转 | hash(id) 确定性，永不复现 |
| P6 验证盲区 | 事后补、字段级 | 分层前置校验 |
| P7 对比困难 | 无版本 | 数据/评测/代码三绑定 |

**一句话**：重构的核心不是"再多洗一遍数据"，而是**把约束从"事后验证"前移到"构造即强制"**——这 7 类问题大多会在写入 KB 之前就被 schema 和词表挡住。
