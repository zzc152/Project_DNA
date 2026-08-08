#!/bin/bash
# 用修复后 final 数据重跑三层级评测（三卡并行：cuda:0/2/3）
# 每轮一个模型一个层级，worker 3 分片到三卡
set -u
cd /workspace/zzc/BioDesign-Agent
export PYTHONPATH=src
PY=/workspace/zzc/envs/project_800/bin/python
DRAFTS=data/synthetic/drafts_full2_clean_final.jsonl
KB=data/processed/knowledge_base_clean_final.jsonl
GPUS=(0 2 3)
log() { echo "[$(date '+%H:%M:%S')] $1"; }

# ---------- 1) L3 32B 三卡并行 ----------
log "=== L3 32B 三卡并行 ==="
for i in 0 1 2; do
  gpu=${GPUS[$i]}
  nohup $PY src/data_synthesis/l3_benchmark_v2.py \
    --drafts $DRAFTS --model models/Qwen2.5-32B-AWQ --device cuda:$gpu \
    --out tmp/l3_final_32b_w$i.jsonl --worker 3 --worker-id $i \
    > tmp/l3_final_32b_w$i.log 2>&1 &
done
wait
log "=== L3 32B 完成 ==="

# ---------- 2) L3 7B 三卡并行 ----------
log "=== L3 7B 三卡并行 ==="
for i in 0 1 2; do
  gpu=${GPUS[$i]}
  nohup $PY src/data_synthesis/l3_benchmark_v2.py \
    --drafts $DRAFTS --model Qwen2.5-7B-Instruct --device cuda:$gpu \
    --out tmp/l3_final_7b_w$i.jsonl --worker 3 --worker-id $i \
    > tmp/l3_final_7b_w$i.log 2>&1 &
done
wait
log "=== L3 7B 完成 ==="

# ---------- 3) L1/L2 32B 三卡并行 ----------
log "=== L1/L2 32B 三卡并行 ==="
for i in 0 1 2; do
  gpu=${GPUS[$i]}
  nohup $PY src/data_synthesis/l12_benchmark.py \
    --kb $KB --model models/Qwen2.5-32B-AWQ --device cuda:$gpu \
    --out tmp/l12_final_32b_w$i.jsonl --worker 3 --worker-id $i \
    > tmp/l12_final_32b_w$i.log 2>&1 &
done
wait
log "=== L1/L2 32B 完成 ==="

# ---------- 4) L1/L2 7B 三卡并行 ----------
log "=== L1/L2 7B 三卡并行 ==="
for i in 0 1 2; do
  gpu=${GPUS[$i]}
  nohup $PY src/data_synthesis/l12_benchmark.py \
    --kb $KB --model Qwen2.5-7B-Instruct --device cuda:$gpu \
    --out tmp/l12_final_7b_w$i.jsonl --worker 3 --worker-id $i \
    > tmp/l12_final_7b_w$i.log 2>&1 &
done
wait
log "=== L1/L2 7B 完成 ==="

# ---------- 5) 合并 + 评分 ----------
log "=== 合并 + 评分 ==="
cat tmp/l3_final_32b_w{0,1,2}.jsonl > tmp/l3_final_32b_merged.jsonl
cat tmp/l3_final_7b_w{0,1,2}.jsonl  > tmp/l3_final_7b_merged.jsonl
cat tmp/l12_final_32b_w{0,1,2}.jsonl > tmp/l12_final_32b_merged.jsonl
cat tmp/l12_final_7b_w{0,1,2}.jsonl  > tmp/l12_final_7b_merged.jsonl

$PY src/data_synthesis/l3_benchmark_v2.py \
  --drafts $DRAFTS --score-only tmp/l3_final_32b_merged.jsonl \
  --out tmp/l3_final_32b_scores.jsonl 2>&1 | tail -5
$PY src/data_synthesis/l3_benchmark_v2.py \
  --drafts $DRAFTS --score-only tmp/l3_final_7b_merged.jsonl \
  --out tmp/l3_final_7b_scores.jsonl 2>&1 | tail -5
$PY src/data_synthesis/l12_benchmark.py \
  --kb $KB --score-only tmp/l12_final_32b_merged.jsonl \
  --out tmp/l12_final_32b_scores.jsonl 2>&1 | tail -5
$PY src/data_synthesis/l12_benchmark.py \
  --kb $KB --score-only tmp/l12_final_7b_merged.jsonl \
  --out tmp/l12_final_7b_scores.jsonl 2>&1 | tail -5

log "=== 全部完成 ==="
