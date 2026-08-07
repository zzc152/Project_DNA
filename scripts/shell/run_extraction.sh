#!/bin/bash
# Module 1 知识抽取：在远程服务器 /workspace/zzc/BioDesign-Agent 下执行
set -e

# 激活 conda 环境 project_800（环境路径：/workspace/zzc/envs/project_800）
export PATH="/workspace/zzc/envs/project_800/bin:$PATH"
CONDA_BASE="$(conda info --base 2>/dev/null || echo /opt/conda)"
if [ -f "$CONDA_BASE/etc/profile.d/conda.sh" ]; then
    source "$CONDA_BASE/etc/profile.d/conda.sh"
    conda activate /workspace/zzc/envs/project_800
else
    echo "[WARN] 未找到 conda.sh，尝试直接使用 python"
fi

cd /workspace/zzc/BioDesign-Agent
echo "工作目录: $(pwd)"
echo "Python: $(which python) $(python --version 2>&1)"
echo "GPU: $(nvidia-smi --query-gpu=index,memory.total --format=csv,noheader 2>/dev/null | head -4)"

# 1. 创建测试集
echo ""
echo "===== 步骤1: 创建测试集 ====="
mkdir -p data/samples
head -10 data/raw/abstracts.jsonl > data/samples/test_abstracts.jsonl
wc -l data/samples/test_abstracts.jsonl

# 2. 小批量测试
echo ""
echo "===== 步骤2: 小批量测试（10 条）====="
python scripts/extract/extract_knowledge.py \
    --input data/samples/test_abstracts.jsonl \
    --output data/processed/raw_extractions_test.jsonl \
    --batch-size 4

echo ""
echo "测试集抽取结果样例:"
head -2 data/processed/raw_extractions_test.jsonl

# 3. 全量抽取
echo ""
echo "===== 步骤3: 全量抽取（500 条）====="
python scripts/extract/extract_knowledge.py \
    --input data/raw/abstracts.jsonl \
    --output data/processed/raw_extractions.jsonl \
    --batch-size 8

# 4. 构建知识库（旧版三元组流程，已废弃，仅保留路径）
echo ""
echo "===== 步骤4: 构建知识三元组（已废弃，请改用 build/ 下 claim 版本）====="
python scripts/legacy/build_knowledge_base.py \
    --input data/processed/raw_extractions.jsonl \
    --output data/processed/knowledge_base.jsonl

echo ""
echo "===== 全部完成 ====="
ls -lh data/processed/
