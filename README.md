# BioDesign-Agent

基因调控序列设计 Agent。核心工作流：PubMed 文献/MPRA 统计知识 → **claim 版本知识库** → 指导序列设计。

## 知识库（Claim 版本，唯一权威）

**所有下游操作统一基于 `data/processed/knowledge_base_stat.jsonl`**（8 键 schema）：

```json
{
  "claim": "...",
  "claim_type": "association | mechanistic | design_rule",
  "entities": { "factor": "...", "factor_type": "...", "cell_line": "...", "regulatory_element": "...", "effect": "...", "pmids": [...] },
  "evidence": { "source": "...", "direction": "...", "p_value": ..., "effect_size": ... },
  "reasoning_chain": ["..."],
  "experimental_context": { "assay": "...", "cell_type": "...", "source_type": "PubMed" },
  "confidence": 0.0
}
```

构建脚本（追加模式，重跑前需先截断）：

| 脚本 | 用途 |
|------|------|
| `scripts/build/build_stat_claims.py` | 统计知识源（MPRA 分析结果）claims |
| `scripts/build/build_literature_claims.py` | 文献知识源（activity 抽取 findings）claims |

> ⚠️ 旧版三元组格式（head/relation/tail）已废弃，文件移至 `data/backup/`，脚本移至 `scripts/legacy/`。

## 项目结构

```
BioDesign-Agent/
├── configs/           # 配置文件（default_config.yaml）
├── data/
│   ├── raw/           # 原始数据（abstracts_*.jsonl）
│   ├── samples/       # 测试样本
│   ├── processed/     # 处理产物（含知识库 knowledge_base_stat.jsonl）
│   ├── backup/        # 废弃知识库备份
│   └── logs/          # 运行日志
├── scripts/
│   ├── download/      # 数据下载（PubMed/MPRA）
│   ├── extract/       # 知识抽取（LLM 抽取 findings）
│   ├── analysis/      # 序列分析（MPRA 分组/基序/GC/DNA shape）
│   ├── build/         # ★ claim 知识库构建（主线）
│   ├── legacy/        # 废弃脚本（旧三元组流程）
│   ├── debug/         # 调试工具
│   └── shell/         # shell 流水线脚本
├── src/
│   ├── extraction/    # 抽取模块（prompts、extractor）
│   └── knowledge/     # 知识模块（triplet_builder、mpra_loader 等）
├── tests/             # 测试
├── docs/              # 文档
├── notebooks/         # Notebook
└── Qwen2.5-7B-Instruct/  # 本地模型（safetensors）
```

## 环境

- Python 3.10（conda 环境：`/workspace/zzc/envs/project_800`）
- 依赖：见 `requirements.txt`（torch、transformers、accelerate、biopython、tqdm）

## 快速开始

```bash
# 1. 下载文献
python scripts/download/download_pubmed_activity.py --output data/raw/abstracts_activity.jsonl

# 2. 抽取知识 findings
python scripts/extract/extract_activity_knowledge.py --input data/raw/abstracts_activity.jsonl

# 3. 构建 claim 知识库（追加）
python scripts/build/build_literature_claims.py --output data/processed/knowledge_base_stat.jsonl
```
