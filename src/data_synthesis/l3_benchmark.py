#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
l3_benchmark.py — L3 层级能力 benchmark（Module 2）

设计（2026-08-07）：
  L3 是当前数据中最好的 benchmark 层级——每道题带 ground truth
  （entities_used 方向/实体 + reasoning_chain 参照），测"忠实机制解释"。

评分维度（满分 100）：
  ① 方向一致性（30）：输出中 effect 方向词与 ground truth 匹配
  ② 实体覆盖（30）：factor/tf(10) + regulatory_element(10) + cell_line(10)
  ③ 结构完整性（20）：L3_001/002 分步+结论；L3_003 总体判断
  ④ 忠实性（20）：无自补越界句（链外机制/文献没提+但…）；出现即 0 分

用法：
  # 生成 + 评分（单进程，--worker 分片可并行）
  CUDA_VISIBLE_DEVICES=3 python src/data_synthesis/l3_benchmark.py \
      --drafts data/synthetic/drafts_full2.jsonl \
      --model Qwen2.5-7B-Instruct --out tmp/l3_eval.jsonl

  # 仅评分已有输出
  python src/data_synthesis/l3_benchmark.py --score-only tmp/l3_eval.jsonl
"""

import argparse
import json
import os
import re
import sys
from collections import Counter, defaultdict


# ---------------------------------------------------------------------------
# ① 方向一致性
# ---------------------------------------------------------------------------

# effect ground truth → 输出中应出现的方向词（中英）
_DIRECTION_WORDS = {
    "increases": ["增加", "增强", "促进", "激活", "上调", "正相关",
                  "increases", "increased", "positive", "upregulat"],
    "decreases": ["降低", "抑制", "减少", "下调", "负相关", "decreases",
                  "decreased", "negative", "downregulat", "repress"],
    "required_for": ["必需", "必要", "需要", "required", "essential",
                     "necessary", "依赖"],
    "modulates": ["调节", "调控", "modulat"],
    "no_effect": ["无影响", "不影响", "no effect", "no significant"],
    "cell_specific_enriched": ["特异", "富集", "cell-specific", "enriched",
                               "specific to"],
}
# L3_003 的 evidence.direction 可能是 cell_specific_enriched 等统计方向
_DIRECTION_ALIAS = {
    "increases_activity": "increases",
    "decreases_activity": "decreases",
    "required_for_activity": "required_for",
    "cell_specific_enriched": "cell_specific_enriched",
    "cell_specific_depleted": "cell_specific_enriched",  # 近似复用
}


def _expected_direction(rec: dict) -> str:
    """该记录期望的方向键（规范化）。"""
    ents = rec["metadata"].get("entities_used", {}) or {}
    ev = rec["metadata"].get("evidence", {}) or {}
    raw = (ents.get("effect") or ev.get("direction") or "").strip()
    raw = _DIRECTION_ALIAS.get(raw, raw)
    return raw.lower()


def _direction_hit(output: str, direction: str) -> bool:
    """输出是否出现与期望方向一致的方向词。"""
    words = _DIRECTION_WORDS.get(direction)
    if not words:
        return False
    out_l = output.lower()
    return any(w.lower() in out_l for w in words)


# ---------------------------------------------------------------------------
# ② 实体覆盖
# ---------------------------------------------------------------------------

def _norm_ent(v) -> str:
    if isinstance(v, (list, tuple)):
        v = v[0] if v else ""
    v = str(v or "").strip()
    if v.lower() in ("not_specified", "none", "other", "tf"):
        return ""
    return v


# 元件中英对照（regulatory_element 取值有限，中文回答时精确匹配失败）
_EL_ZH = {
    "enhancer": ["增强子", "增强元件"],
    "promoter": ["启动子", "启动元件"],
    "super-enhancer": ["超级增强子", "超增强子"],
    "cis-regulatory element": ["顺式调控元件", "顺式调节元件"],
    "tf binding motifs": ["转录因子结合基序", "结合基序"],
}


def _ent_hit(name: str, output: str) -> bool:
    """实体是否被提及：精确子串 → 大小写不敏感 → 去连字符/空格归一化
    → 英文 token 级（len>=3）→ 中英对照。"""
    if str(name) in output:
        return True
    out_l = output.lower()
    # 大小写不敏感
    if str(name).lower() in out_l:
        return True
    # 去连字符/空格归一化（ETS1 vs Ets-1、FOXA1 vs FoxA1）
    norm = re.sub(r"[\s\-]", "", str(name)).lower()
    if norm and norm in re.sub(r"[\s\-]", "", output).lower():
        return True
    # 英文 token 级（中文回答常保留英文核心词，如 macroH2A、MEF2C）
    toks = [t for t in re.split(r"[^A-Za-z0-9]+", str(name))
            if len(t) >= 3 and not t.isdigit()]
    if toks and toks[0].lower() in out_l:
        return True
    # 中英对照
    for zh in _EL_ZH.get(str(name).lower(), []):
        if zh in output:
            return True
    return False


def _entity_checks(rec: dict) -> dict:
    """返回 {实体名: 得分权重}。"""
    ents = rec["metadata"].get("entities_used", {}) or {}
    checks = {}
    # factor 型（L3_001/002）
    f = _norm_ent(ents.get("factor"))
    if f:
        checks[f] = 10
    tf = _norm_ent(ents.get("tf"))
    if tf:
        checks[tf] = 10
    # 元件 + 细胞系各 10 分
    el = _norm_ent(ents.get("regulatory_element"))
    if el:
        checks[el] = 10
    cl = _norm_ent(ents.get("cell_line"))
    if cl:
        checks[cl] = 10
    return checks


# ---------------------------------------------------------------------------
# ③ 结构完整性
# ---------------------------------------------------------------------------

_STEP_RE = re.compile(r"(?:^|\n)\s*(?:\d+[\.、\)）]|step\s*\d+|\*\*\d+\.|第一步|第二步|第三步)")
_CONCLUSION_RE = re.compile(r"因此|总之|综上|综上所述|结论|conclusion|in summary|to summarize")
_VERDICT_RE = re.compile(r"支持|不支持|一致|冲突|总体判断|总体结论|verdict|support|consistent|contradict")


def _structure_ok(rec: dict, output: str) -> tuple:
    """返回 (是否通过, 明细描述)。"""
    tid = rec["template_id"]
    if tid == "L3_003":
        ok = bool(_VERDICT_RE.search(output))
        return ok, ("有总体判断" if ok else "缺总体判断")
    # L3_001/002：分步 + 结论
    has_step = bool(_STEP_RE.search(output))
    has_concl = bool(_CONCLUSION_RE.search(output))
    ok = has_step and has_concl
    detail = []
    if not has_step:
        detail.append("无分步编号")
    if not has_concl:
        detail.append("无结论句")
    return ok, ("+".join(detail) if detail else "分步+结论完整")


# ---------------------------------------------------------------------------
# ④ 忠实性（自补越界复用 quality_filter 标记）
# ---------------------------------------------------------------------------

_FABRICATION_MARKERS = [
    "文献未提及", "文献未直接提及", "文献中没有", "原文未提及",
    "虽然文献未", "尽管文献未",
]


# ---------------------------------------------------------------------------
# 评分主函数
# ---------------------------------------------------------------------------

def score_l3(rec: dict, output: str) -> dict:
    """对单条 L3 生成输出评分。返回维度明细 + 总分。"""
    out = output or ""
    dims = {}

    # ① 方向一致性 30
    exp_dir = _expected_direction(rec)
    dir_ok = _direction_hit(out, exp_dir) if exp_dir else True
    dims["方向一致性(30)"] = (30 if dir_ok else 0)

    # ② 实体覆盖 30
    checks = _entity_checks(rec)
    hit = sum(w for e, w in checks.items() if _ent_hit(e, out))
    dims["实体覆盖(30)"] = hit
    missed = [e for e, w in checks.items() if not _ent_hit(e, out)]

    # ③ 结构完整性 20
    struct_ok, struct_detail = _structure_ok(rec, out)
    dims["结构完整性(20)"] = (20 if struct_ok else 0)

    # ④ 忠实性 20（出现任一自补越界句 → 0 分）
    fab = [m for m in _FABRICATION_MARKERS if m in out]
    dims["忠实性(20)"] = (0 if fab else 20)

    total = sum(dims.values())
    return {
        "id": rec["id"],
        "template_id": rec["template_id"],
        "dims": dims,
        "total": total,
        "direction": exp_dir,
        "missed_entities": missed,
        "fabrication": fab,
        "structure": struct_detail,
    }


# ---------------------------------------------------------------------------
# 评测主流程：加载模型 → 生成 → 评分
# ---------------------------------------------------------------------------

def load_model(model_path: str, device: str):
    from transformers import AutoModelForCausalLM, AutoTokenizer
    tok = AutoTokenizer.from_pretrained(model_path, local_files_only=True)
    model = AutoModelForCausalLM.from_pretrained(
        model_path, device_map={"": device},
        torch_dtype="auto", local_files_only=True)
    return model, tok


def generate(model, tok, instruction: str, max_new_tokens: int) -> str:
    msgs = [{"role": "user", "content": instruction}]
    prompt = tok.apply_chat_template(
        msgs, tokenize=False, add_generation_prompt=True)
    inputs = tok(prompt, return_tensors="pt")
    inputs = {k: v.to(next(model.parameters()).device) for k, v in inputs.items()}
    out = model.generate(
        **inputs, max_new_tokens=max_new_tokens,
        do_sample=False,  # greedy，基准测试不引入采样噪声
    )
    new = out[0][inputs["input_ids"].shape[1]:]
    return tok.decode(new, skip_special_tokens=True).strip()


def report(results: list) -> dict:
    """汇总报告：按模板 + 整体。"""
    by_tpl = defaultdict(list)
    for r in results:
        by_tpl[r["template_id"]].append(r)
    summary = {"by_template": {}, "overall": {}}
    all_totals = [r["total"] for r in results]
    summary["overall"] = {
        "n": len(results),
        "avg_score": round(sum(all_totals) / len(all_totals), 1) if all_totals else 0,
        "pass80": sum(1 for t in all_totals if t >= 80),
    }
    for tid, rs in sorted(by_tpl.items()):
        tot = [r["total"] for r in rs]
        dim_avg = {}
        for k in rs[0]["dims"]:
            dim_avg[k] = round(sum(r["dims"][k] for r in rs) / len(rs), 1)
        summary["by_template"][tid] = {
            "n": len(rs),
            "avg_score": round(sum(tot) / len(tot), 1),
            "dims": dim_avg,
        }
    return summary


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--drafts", default="data/synthetic/drafts_full2.jsonl")
    ap.add_argument("--model", default="Qwen2.5-7B-Instruct")
    ap.add_argument("--device", default="cuda:3")
    ap.add_argument("--out", default="tmp/l3_eval.jsonl")
    ap.add_argument("--max-new-tokens", type=int, default=512)
    ap.add_argument("--limit", type=int, default=0, help="仅前 N 条 L3（0=全部）")
    ap.add_argument("--worker", type=int, default=1)
    ap.add_argument("--worker-id", type=int, default=0)
    ap.add_argument("--score-only", default="",
                    help="仅对已有输出文件评分，不生成（传入路径）")
    args = ap.parse_args()

    # 评分已有输出
    if args.score_only:
        results = []
        for l in open(args.score_only, encoding="utf-8"):
            if not l.strip():
                continue
            rec = json.loads(l)
            results.append(score_l3(rec, rec.get("output", "")))
        print(json.dumps(report(results), ensure_ascii=False, indent=1))
        # 逐条明细落盘（--out 直接作为明细输出路径）
        with open(args.out, "w", encoding="utf-8") as f:
            for r in results:
                f.write(json.dumps(r, ensure_ascii=False) + "\n")
        return

    # 生成 + 评分
    drafts = [json.loads(l) for l in open(args.drafts, encoding="utf-8")]
    l3 = [d for d in drafts if d["level"] == "L3"]
    if args.limit > 0:
        l3 = l3[:args.limit]
    if args.worker > 1:
        l3 = [d for i, d in enumerate(l3) if i % args.worker == args.worker_id]
    existing = {}
    if os.path.exists(args.out):
        for l in open(args.out, encoding="utf-8"):
            if l.strip():
                r = json.loads(l)
                existing[r["id"]] = r
    pending = [d for d in l3 if d["id"] not in existing]
    print(f"[L3 评测] 待测 {len(pending)} / {len(l3)} | "
          f"[worker {args.worker_id}/{args.worker}]", flush=True)

    model, tok = None, None
    if pending:
        print(f"[加载模型] {args.model} -> {args.device}", flush=True)
        model, tok = load_model(args.model, args.device)
        print("[模型就绪]", flush=True)

    os.makedirs(os.path.dirname(args.out) or ".", exist_ok=True)
    n = 0
    for d in pending:
        try:
            output = generate(model, tok, d["instruction"], args.max_new_tokens)
            rec = dict(d)
            rec["output"] = output
            with open(args.out, "a", encoding="utf-8") as f:
                f.write(json.dumps(rec, ensure_ascii=False) + "\n")
            s = score_l3(rec, output)
            n += 1
            print(f"  {s['id']} total={s['total']} "
                  f"[dir_hit={s['dims']['方向一致性(30)']>0} "
                  f"ents={s['dims']['实体覆盖(30)']}/30 miss={s['missed_entities']}]",
                  flush=True)
        except Exception as e:  # 单条失败不中断
            print(f"  ✗ {d['id']} 生成失败: {e}", flush=True)
    print(f"[完成] 新增 {n} 条", flush=True)

    # 汇总评分（含断点已生成的）
    all_recs = []
    for l in open(args.out, encoding="utf-8"):
        if l.strip():
            all_recs.append(json.loads(l))
    results = [score_l3(r, r.get("output", "")) for r in all_recs]
    print(json.dumps(report(results), ensure_ascii=False, indent=1))
    with open(args.out.replace(".jsonl", "_scores.jsonl"), "w",
              encoding="utf-8") as f:
        for r in results:
            f.write(json.dumps(r, ensure_ascii=False) + "\n")


if __name__ == "__main__":
    main()
