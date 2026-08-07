# -*- coding: utf-8 -*-
"""模块三：L1/L2 benchmark（知识回忆 + 知识关联）。

基线 v2.1（l3_benchmark_v2.py）只评测了 L3 记录（130 条）。
本文件为 L1/L2 层级构建独立 benchmark 数据集，核心原则：

**信息量递减**（用户明确要求）：
  - L1（Recall，单跳）：零上下文。只给 (factor, element)，问影响方向。
    不给 claim 文本、不给文献摘要、不给推理链 —— 纯知识回忆。
  - L2（Association，双跳）：带细胞系约束。给 (factor, element, cell)，问影响方向。
    仍不给 claim/摘要/推理链 —— 需要因子→元件→细胞系的关联知识。
  - L3（Mechanistic，v2.1）：rc 给摘要原文 / know 不给 —— 多跳推理。
  对比现有 L1/L2 模板（instruction 把 claim_text 完整塞进题面）：
  本 benchmark 删除题面中的 claim 文本与摘要，避免"照抄答案"。

数据源：knowledge_base_clean.jsonl（KB 格式，entities 顶层字段）。
评分：复用 v2.1 哲学（答案 60 + 忠实性 40 + 诚实性），答 E 给诚实分。

用法：
  python src/data_synthesis/l12_benchmark.py --dry-run
  python src/data_synthesis/l12_benchmark.py --model <M> --device cuda:N \
      --out tmp/l12_eval.jsonl --worker 3 --worker-id 0
  python src/data_synthesis/l12_benchmark.py --score-only tmp/l12_eval.jsonl --out tmp/l12_scores.jsonl
"""

import argparse
import json
import os
import random
import sys
from collections import defaultdict

try:
    from l3_benchmark_v2 import (  # noqa: E402
        _norm_ent, _gold_direction, _OPTS_DIR, _HONESTY_PROMPT, _FABRICATION_MARKERS,
        load_model, generate, report, parse_answer, _balance_directions,
    )
except ModuleNotFoundError:
    sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", ".."))
    from data_synthesis.l3_benchmark_v2 import (  # noqa: E402
        _norm_ent, _gold_direction, _OPTS_DIR, _HONESTY_PROMPT, _FABRICATION_MARKERS,
        load_model, generate, report, parse_answer, _balance_directions,
    )


# ---------------------------------------------------------------------------
# KB 字段读取（KB 格式：entities 是顶层字段）
# ---------------------------------------------------------------------------

def _kb_factor(r):
    return _norm_ent((r.get("entities") or {}).get("factor"))


def _kb_element(r):
    return _norm_ent((r.get("entities") or {}).get("regulatory_element"))


def _kb_cell(r):
    return _norm_ent((r.get("entities") or {}).get("cell_line"))


def _kb_effect(r):
    return (r.get("entities") or {}).get("effect")


def _has_effect(r) -> bool:
    e = _kb_effect(r)
    return bool(e) and str(e).lower() not in ("none", "null", "")


def _clean_factor(r) -> str:
    """过滤不可考的因子名（描述性短语、过长名），返回合法因子名或 None。"""
    f = _kb_factor(r)
    if not f:
        return None
    # 排除明显是句子/短语的噪声名（如 "targeted activation of this regulatory element"）
    if len(f) > 30 or f.count(" ") >= 4:
        return None
    return f

# ---------------------------------------------------------------------------
# 题面构造
# ---------------------------------------------------------------------------

def _build_q(qid: str, level: str, factor: str, element: str, cell: str,
             effect: str, n_ctx: int = 0) -> dict:
    """构造一道方向判断题。

    n_ctx 控制上下文信息量：
      0 = L1 零上下文（只给因子+元件）
      1 = L2 带细胞系（双跳）
    """
    gold = _gold_direction(effect)
    if n_ctx >= 1 and cell:
        q = (f"问题：在{cell}细胞中，因子{factor}对{element or '靶元件'}活性的影响方向是什么？\n"
             f"选项：\n{_OPTS_DIR}")
    else:
        q = (f"问题：因子{factor}对{element or '靶元件'}活性的影响方向是什么？\n"
             f"选项：\n{_OPTS_DIR}")
    q += ("\n\n请只输出一个选项字母（如 A），并附一句话理由：\n"
          "例如：“A。该因子是正向调节因子。”\n"
          f"{_HONESTY_PROMPT}")
    return {"qtype": level, "question": q, "gold": gold, "opts": _OPTS_DIR,
            "factor": factor, "element": element, "cell": cell, "effect": effect}


def _balanced_sample(qs: dict, rng: random.Random, max_qs: int) -> dict:
    """A/B 严格各半 + CDE 补位到 max_qs。

    与 v2 的 _balance_directions 区别：v2 只削多余 A/B，不处理超量截断
    （先平衡后截断会破坏平衡）。这里 A 组/B 组各抽 max_qs//2 条（不足则全取），
    剩余配额由 C/D/E 随机补足，保证最终 A 数与 B 数严格相等。
    """
    a_ids = [i for i, q in qs.items() if q["gold"] == "A"]
    b_ids = [i for i, q in qs.items() if q["gold"] == "B"]
    other = [(i, q) for i, q in qs.items() if q["gold"] not in ("A", "B")]
    keep = min(len(a_ids), len(b_ids), max_qs // 2)
    sel = set(rng.sample(a_ids, keep)) | set(rng.sample(b_ids, keep))
    space = max_qs - 2 * keep
    if space > 0 and other:
        sel |= set(i for i, _ in rng.sample(other, min(space, len(other))))
    return {i: q for i, q in qs.items() if i in sel}


def build_l1_qs(recs: list, rng: random.Random, max_qs: int = 60) -> dict:
    """L1 知识回忆（单跳，零上下文）。

    去重键 (factor, element)：同组合多条记录只保留一条（同一事实不重复考）。
    """
    pool = [r for r in recs if _has_effect(r) and _clean_factor(r) and _kb_element(r)]
    seen, qs = set(), {}
    for r in pool:
        key = (_clean_factor(r), _kb_element(r))
        if key in seen:
            continue
        seen.add(key)
        qid = f"L1_{len(qs) + 1:04d}"
        qs[qid] = _build_q(qid, "L1", _clean_factor(r), _kb_element(r),
                           _kb_cell(r), _kb_effect(r), n_ctx=0)
        if len(qs) >= max_qs * 3:  # 预采样足够候选供平衡
            break
    return _balanced_sample(qs, rng, max_qs)


def build_l2_qs(recs: list, rng: random.Random, max_qs: int = 60) -> dict:
    """L2 知识关联（双跳，带细胞系约束）。

    去重键 (factor, element, cell)：需要具体细胞系（not_specified 不构成双跳约束）。
    """
    pool = [r for r in recs if _has_effect(r) and _clean_factor(r) and _kb_element(r)
            and _kb_cell(r)]
    seen, qs = set(), {}
    for r in pool:
        key = (_clean_factor(r), _kb_element(r), _kb_cell(r))
        if key in seen:
            continue
        seen.add(key)
        qid = f"L2_{len(qs) + 1:04d}"
        qs[qid] = _build_q(qid, "L2", _clean_factor(r), _kb_element(r),
                           _kb_cell(r), _kb_effect(r), n_ctx=1)
        if len(qs) >= max_qs * 3:
            break
    return _balanced_sample(qs, rng, max_qs)


def build_all(recs: list, rng: random.Random,
              n_l1: int = 60, n_l2: int = 60) -> dict:
    qs = {}
    qs.update(build_l1_qs(recs, rng, n_l1))
    qs.update(build_l2_qs(recs, rng, n_l2))
    return qs


# ---------------------------------------------------------------------------
# 评分（复用 v2 哲学，实体从 KB 读）
# ---------------------------------------------------------------------------

def _key_entity_hit_kb(q: dict, output: str) -> bool:
    """理由中是否引用题面关键实体（因子/元件/细胞系）。"""
    out = re_sub_nospace(output or "")
    for name in (q.get("factor"), q.get("element"), q.get("cell")):
        if not name:
            continue
        n = re_sub_nospace(str(name))
        if n and n in out:
            return True
    return False


def re_sub_nospace(s: str) -> str:
    import re
    return re.sub(r"[\s\-]", "", s).lower()


def score_l12(qid: str, q: dict, output: str) -> dict:
    """评分：答案正确性 60 + 理由忠实性 40 + 诚实性字段。

    L1/L2 均为"不给信息"题型 → 答 E（无法判断）给 30 诚实分（同 v2 know 型）。
    """
    out = output or ""
    ans = parse_answer(out)
    gold = q["gold"]

    if ans == gold:
        acc = 60
    elif ans == "":
        acc = 0
    elif ans == "E" and gold != "E":
        acc = 30  # 诚实分：题面确实没给信息
    else:
        acc = 0

    ent_hit = _key_entity_hit_kb(q, out)
    fab = [m for m in _FABRICATION_MARKERS if m in out]
    fid = (20 if ent_hit else 0) + (0 if fab else 20)

    claimed = ans in ("A", "B", "C", "D")
    overclaim = claimed and ans != gold

    total = acc + fid
    return {
        "id": qid, "template_id": q["qtype"], "qtype": q["qtype"],
        "gold": gold, "answer": ans, "acc": acc, "fidelity": fid,
        "ent_hit": ent_hit, "fabrication": fab,
        "claimed": claimed, "overclaim": overclaim, "total": total,
    }


# ---------------------------------------------------------------------------
# 主流程
# ---------------------------------------------------------------------------

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--kb", default="data/processed/knowledge_base_clean.jsonl")
    ap.add_argument("--model", default="Qwen2.5-7B-Instruct")
    ap.add_argument("--device", default="cuda:3")
    ap.add_argument("--out", default="tmp/l12_eval.jsonl")
    ap.add_argument("--max-new-tokens", type=int, default=256)
    ap.add_argument("--seed", type=int, default=42)
    ap.add_argument("--worker", type=int, default=1)
    ap.add_argument("--worker-id", type=int, default=0)
    ap.add_argument("--n-l1", type=int, default=60)
    ap.add_argument("--n-l2", type=int, default=60)
    ap.add_argument("--score-only", default="")
    ap.add_argument("--dry-run", action="store_true")
    ap.add_argument("--dry-run-n", type=int, default=6)
    args = ap.parse_args()

    recs = [json.loads(l) for l in open(args.kb, encoding="utf-8")]
    rng = random.Random(args.seed)
    questions = build_all(recs, rng, args.n_l1, args.n_l2)

    if args.dry_run:
        for i, (qid, q) in enumerate(questions.items()):
            if i >= args.dry_run_n:
                break
            print("=" * 70)
            print(f"id={qid} | qtype={q['qtype']} | gold={q['gold']}")
            print(q["question"])
        from collections import Counter
        c = Counter((q["qtype"], q["gold"]) for q in questions.values())
        print("=" * 70)
        print(f"[dry-run] 共构造 {len(questions)} 条题面")
        for k, v in sorted(c.items()):
            print(f"  {k[0]} gold={k[1]}: {v}")
        return

    ids = list(questions.keys())
    if args.worker > 1:
        ids = [i for idx, i in enumerate(ids) if idx % args.worker == args.worker_id]

    if args.score_only:
        results = []
        for l in open(args.score_only, encoding="utf-8"):
            if not l.strip():
                continue
            rec = json.loads(l)
            q = questions.get(rec.get("id", ""))
            if not q:
                continue
            results.append(score_l12(rec.get("id"), q, rec.get("output", "")))
        print(json.dumps(report(results), ensure_ascii=False, indent=1))
        with open(args.out, "w", encoding="utf-8") as f:
            for r in results:
                f.write(json.dumps(r, ensure_ascii=False) + "\n")
        return

    existing = {}
    if os.path.exists(args.out):
        for l in open(args.out, encoding="utf-8"):
            if l.strip():
                r = json.loads(l)
                existing[r["id"]] = r
    pending = [i for i in ids if i not in existing]
    print(f"[L1/L2 评测] 待测 {len(pending)} / {len(ids)} | "
          f"[worker {args.worker_id}/{args.worker}]", flush=True)

    model, tok = None, None
    if pending:
        print(f"[加载模型] {args.model} -> {args.device}", flush=True)
        model, tok = load_model(args.model, args.device)
        print("[模型就绪]", flush=True)

    os.makedirs(os.path.dirname(args.out) or ".", exist_ok=True)
    n = 0
    for qid in pending:
        try:
            q = questions[qid]
            output = generate(model, tok, q["question"], args.max_new_tokens)
            rec = {"id": qid, "template_id": q["qtype"], "qtype": q["qtype"],
                   "question": q["question"], "gold": q["gold"],
                   "output": output}
            with open(args.out, "a", encoding="utf-8") as f:
                f.write(json.dumps(rec, ensure_ascii=False) + "\n")
            s = score_l12(qid, q, output)
            n += 1
            print(f"  {s['id']} total={s['total']} gold={s['gold']} "
                  f"ans={s['answer']} [{s['qtype']}]", flush=True)
        except Exception as e:
            print(f"  ✗ {qid} 生成失败: {e}", flush=True)
    print(f"[完成] 新增 {n} 条", flush=True)

    all_recs = []
    for l in open(args.out, encoding="utf-8"):
        if l.strip():
            r = json.loads(l)
            q = questions.get(r.get("id", ""))
            if q:
                all_recs.append(score_l12(r["id"], q, r.get("output", "")))
    print(json.dumps(report(all_recs), ensure_ascii=False, indent=1))
    with open(args.out + "_scores.jsonl", "w", encoding="utf-8") as f:
        for r in all_recs:
            f.write(json.dumps(r, ensure_ascii=False) + "\n")


if __name__ == "__main__":
    main()
