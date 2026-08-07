#!/usr/bin/env bash
# -*- coding: utf-8 -*-
# run_full_pipeline.sh — 全量重建+增强+质量校验（一条命令自动完成）
# 用法: bash scripts/run_full_pipeline.sh [per-template] [out_suffix]
#   默认 per-template=64（14 模板 → 理论 896 条；受 KB 匹配限制）
#   三卡并行增强（GPU0/2/3）→ wait 等全部 worker → 规则版校验 → judge 版后台
# 预估: per-template=64 约 4-6 小时（< 9 小时）
set -e

cd /workspace/zzc/BioDesign-Agent
PY=/workspace/zzc/envs/project_800/bin/python
PT=${1:-64}
SUF=${2:-full}
DRF=data/synthetic/drafts_${SUF}.jsonl
ENH=data/synthetic/enhanced_${SUF}.jsonl
FLT=data/synthetic/filtered_${SUF}.jsonl
mkdir -p tmp

echo "========== [1/3] 重建 drafts (per-template=$PT) =========="
rm -f "$DRF"
$PY src/data_synthesis/template_engine.py \
    --kb data/processed/knowledge_base_clean.jsonl \
    --templates templates/ \
    --out "$DRF" \
    --per-template "$PT" \
    --seed 42 2>&1 | grep -E "产出|层级|模板|匹配" || true

echo "========== [2/3] 三卡并行增强 =========="
rm -f "$ENH" tmp/enhance_${SUF}_*.log
for i in 0 1 2; do
    case $i in
        0) gpu=0;;
        1) gpu=2;;
        2) gpu=3;;
    esac
    CUDA_VISIBLE_DEVICES=$gpu nohup $PY src/data_synthesis/llm_enhancer.py \
        --drafts "$DRF" --out "$ENH" --worker 3 --worker-id $i \
        > tmp/enhance_${SUF}_$i.log 2>&1 &
done
wait || true   # 等三个 worker 完成（容忍单个 worker 失败）
echo "[增强完成] 产出: $(wc -l < "$ENH") 条"
grep -h "完成" tmp/enhance_${SUF}_*.log || true

echo "========== [3/3] 质量校验 =========="
# 规则版（秒级，立即出结果）
CUDA_VISIBLE_DEVICES=0 $PY src/data_synthesis/quality_filter.py \
    --in "$ENH" --out /dev/null 2>/dev/null \
    | grep -E "校验|✗" > tmp/qf_${SUF}_rule.log || true
echo "[规则版] $(cat tmp/qf_${SUF}_rule.log | grep 校验 || echo 无)"
# LLM judge 版（后台，判同义不同形）
CUDA_VISIBLE_DEVICES=0 nohup $PY src/data_synthesis/quality_filter.py \
    --in "$ENH" --out "$FLT" --llm-judge \
    > tmp/qf_${SUF}_judge.log 2>&1 &

echo "========== 主流程 DONE（judge 后台继续）=========="
echo "drafts:   $DRF"
echo "enhanced: $ENH"
echo "filtered: $FLT (judge 完成后)"
echo "日志: tmp/enhance_${SUF}_{0,1,2}.log / tmp/qf_${SUF}_{rule,judge}.log"
