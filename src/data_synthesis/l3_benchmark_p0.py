# -*- coding: utf-8 -*-
"""模块三 P0：两类新增评测题型（扩展 L3 封闭式评测）。

基线（l3_benchmark_v2.py）评测的是"方向判断"（rc/know/data 三题型）。
本文件在基线上新增两个题型，针对更深的科学推理能力：

题型 A：矛盾检测（contradiction）
  - 输入：同一 (factor, element, cell) 的两条真实知识库 claim（或跨组构造对）。
  - 任务：判断两条 claim 是「相互矛盾 / 相互支持 / 无关或无法判断」。
  - 检验能力：跨 claim 的逻辑一致性推理（scientific reasoning），而非单条事实回忆。
  - 构造：
      矛盾对：同组 effect ∈ {increases, decreases} 方向相反（如 MYCN: increases vs decreases）
      支持对：同组 effect 同向（如 required_for vs increases）
      无关对：不同 factor 或不同 element（跨组随机配对）

题型 B：evidence-level 推理（evidence）
  - 输入：一条证据句（摘要原文）+ 一个结论。
  - 任务：判断该证据能否推出该结论（能推出 / 不能推出 / 无法判断）。
  - 检验能力：证据层级判断——表达证据（RNA-seq/mRNA/表达水平）≠ 元件活性证据
    （reporter/MPRA/活性实验）。这正是 D/E 类错误的核心：模型常把"表达上调"
    误推为"元件活性增强"。
  - 构造：
      能推出：活性类证据（reporter/luciferase/MPRA/enhancer activity）+ 对应活性结论
      不能推出：表达类证据（RNA-seq/expression/mRNA/transcription）+ 活性结论

评分（满分 100，与 v2 一致的哲学）：
  答案正确性 60 + 理由忠实性 40（引用题面关键实体 20 + 无自补编造 20）
  诚实性：声称知道但答错 = overclaim；E/无法判断 视为"不知道"。

用法：
  python src/data_synthesis/l3_benchmark_p0.py --dry-run          # 只看题面
  python src/data_synthesis/l3_benchmark_p0.py --model <M> --device cuda:N \
      --out tmp/p0_eval.jsonl --worker 3 --worker-id 0            # 评测
  python src/data_synthesis/l3_benchmark_p0.py --score-only tmp/p0_eval.jsonl \
      --out tmp/p0_scores.jsonl                                   # 只评分
"""

import argparse
import json
import os
import random
import re
import sys
from collections import defaultdict

# 兼容两种运行方式：脚本直跑（python src/data_synthesis/l3_benchmark_p0.py）
# 或包导入（PYTHONPATH=src 后 from data_synthesis.l3_benchmark_p0 import ...）
try:
    from l3_benchmark_v2 import (  # noqa: E402
        load_model, generate, report,
        _FABRICATION_MARKERS, _HONESTY_PROMPT,
    )
except ModuleNotFoundError:
    sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", ".."))
    from data_synthesis.l3_benchmark_v2 import (  # noqa: E402
        load_model, generate, report,
        _FABRICATION_MARKERS, _HONESTY_PROMPT,
    )

# ---------------------------------------------------------------------------
# 方向语义（判断 claim 对是"矛盾/一致"）
# ---------------------------------------------------------------------------

_POS_EFFECTS = {"increases", "required_for", "increases_activity",
                "required_for_activity", "promotes", "promotes_activity"}
_NEG_EFFECTS = {"decreases", "decreases_activity", "disrupts", "represses"}
_UNKNOWN_EFFECTS = {"modulates", "binds", "no_effect", "affects_interaction",
                    "co_occurs_with", "correlates", "regulates_expression", None}

_ABS_KEYWORDS = ["摘要原文证据", "原始发现", "文献机制描述"]
_EXPR_KEYWORDS = ["rna-seq", "mrna", "expression", "transcript", "express",
                  "表达", "转录", "mrna 水平"]
_ACT_KEYWORDS = ["reporter", "luciferase", "mpra", "enhancer activity",
                 "promoter activity", "活性实验", "报告基因", "增强子活性",
                 "启动子活性", "活性检测"]


def _dir_sem(effect: str) -> str:
    """effect → 语义：pos / neg / unknown。"""
    e = (effect or "").lower()
    if e in _POS_EFFECTS:
        return "pos"
    if e in _NEG_EFFECTS:
        return "neg"
    return "unknown"


def _short_claim(rec: dict, max_len: int = 200) -> str:
    """截断 claim 文本，避免题面过长。"""
    c = rec.get("claim", "")
    return c if len(c) <= max_len else c[:max_len] + "…"


def _claim_entities(rec: dict) -> tuple:
    ent = rec.get("entities", {})
    return (ent.get("factor"), ent.get("regulatory_element"),
            ent.get("cell_line"))


# ---------------------------------------------------------------------------
# 题型 A：矛盾检测题面构造
# ---------------------------------------------------------------------------

_OPTS_CONTRADICT = ("A. 相互矛盾\nB. 相互支持（一致）\nC. 无关或无法判断")


def build_contradiction_pairs(recs: list, rng: random.Random,
                              max_pairs: int = 60) -> list:
    """构造矛盾检测题（真实同组 claim 对 + 跨组无关对）。

    返回 list of dict: {id, qtype, question, gold, opts, claim_a, claim_b}

    平衡策略：分别构造 A/B/C 三类，再分层采样合并，保证三类都进入最终集，
    避免"支持对数量淹没矛盾对"导致的类别缺失。
    """
    # 按 (factor, element, cell) 分组
    groups = defaultdict(list)
    for r in recs:
        groups[_claim_entities(r)].append(r)

    cat = {"A": [], "B": [], "C": []}
    pid = 0

    # 1) 同组内构造矛盾对（A：pos×neg 笛卡尔积，尽量多配）、
    #    支持对（B：同向两两）、无法判断对（C：pos/neg vs 方向不明）
    for key, members in groups.items():
        pos = [r for r in members if _dir_sem(r["entities"].get("effect")) == "pos"]
        neg = [r for r in members if _dir_sem(r["entities"].get("effect")) == "neg"]
        unk = [r for r in members if _dir_sem(r["entities"].get("effect")) == "unknown"]
        # 矛盾对：所有 pos × 所有 neg（同组同因子同元件，方向相反 → 真矛盾）
        for a in pos:
            for b in neg:
                pid += 1
                cat["A"].append({
                    "id": f"P0_CT_{pid:04d}", "qtype": "contradiction",
                    "question": _contradict_q(_short_claim(a), _short_claim(b)),
                    "gold": "A", "opts": _OPTS_CONTRADICT,
                    "claim_a": _short_claim(a), "claim_b": _short_claim(b),
                })
        # 支持对：pos 内两两 + neg 内两两（同向 → 一致）
        for i in range(len(pos)):
            for j in range(i + 1, len(pos)):
                pid += 1
                cat["B"].append({
                    "id": f"P0_CT_{pid:04d}", "qtype": "contradiction",
                    "question": _contradict_q(_short_claim(pos[i]), _short_claim(pos[j])),
                    "gold": "B", "opts": _OPTS_CONTRADICT,
                    "claim_a": _short_claim(pos[i]), "claim_b": _short_claim(pos[j]),
                })
        for i in range(len(neg)):
            for j in range(i + 1, len(neg)):
                pid += 1
                cat["B"].append({
                    "id": f"P0_CT_{pid:04d}", "qtype": "contradiction",
                    "question": _contradict_q(_short_claim(neg[i]), _short_claim(neg[j])),
                    "gold": "B", "opts": _OPTS_CONTRADICT,
                    "claim_a": _short_claim(neg[i]), "claim_b": _short_claim(neg[j]),
                })
        # 无法判断对：pos/neg vs 方向不明（modulates/binds/no_effect → 无法判定）
        for a in (pos + neg):
            for b in unk:
                pid += 1
                cat["C"].append({
                    "id": f"P0_CT_{pid:04d}", "qtype": "contradiction",
                    "question": _contradict_q(_short_claim(a), _short_claim(b)),
                    "gold": "C", "opts": _OPTS_CONTRADICT,
                    "claim_a": _short_claim(a), "claim_b": _short_claim(b),
                })

    # 2) 跨组"无关对"（不同 factor 或不同 element）→ C
    keys = list(groups.keys())
    rng.shuffle(keys)
    for i in range(0, len(keys) - 1, 2):
        k1, k2 = keys[i], keys[i + 1]
        if k1[0] == k2[0] and k1[1] == k2[1]:
            continue
        r1 = rng.choice(groups[k1])
        r2 = rng.choice(groups[k2])
        pid += 1
        cat["C"].append({
            "id": f"P0_CT_{pid:04d}", "qtype": "contradiction",
            "question": _contradict_q(_short_claim(r1), _short_claim(r2)),
            "gold": "C", "opts": _OPTS_CONTRADICT,
            "claim_a": _short_claim(r1), "claim_b": _short_claim(r2),
        })

    # 3) 分层采样：每类各取约 1/3，保证平衡；数量不足的类全保留
    per = max(1, max_pairs // 3)
    rng.shuffle(cat["A"]); rng.shuffle(cat["B"]); rng.shuffle(cat["C"])
    pairs = (cat["A"][:per] + cat["B"][:per] + cat["C"][:per])
    rng.shuffle(pairs)
    return pairs[:max_pairs]


def _contradict_q(claim_a: str, claim_b: str) -> str:
    return (f"以下是两条关于基因调控的科学结论：\n\n"
            f"结论 1：\n“{claim_a}”\n\n"
            f"结论 2：\n“{claim_b}”\n\n"
            f"问题：这两条结论之间是什么关系？\n"
            f"选项：\n{_OPTS_CONTRADICT}\n\n"
            f"请只输出一个选项字母（如 A），并附一句话理由。\n"
            f"{_HONESTY_PROMPT}")


# ---------------------------------------------------------------------------
# 题型 B：evidence-level 推理题面构造
# ---------------------------------------------------------------------------

_OPTS_EVIDENCE = ("A. 能推出\nB. 不能推出（证据层级不足）\nC. 无法判断")


def _extract_evidence_sent(rec: dict) -> str:
    """从 reasoning_chain 中提取"摘要原文证据"句。"""
    rc = rec.get("reasoning_chain", [])
    if isinstance(rc, list):
        for step in rc:
            s = str(step)
            if any(k in s for k in _ABS_KEYWORDS):
                return s
    return ""


def _evid_type(sent: str) -> str:
    """证据类型：expr（表达类）/ act（活性类）/ other。"""
    low = sent.lower()
    has_expr = any(k in low for k in _EXPR_KEYWORDS)
    has_act = any(k in low for k in _ACT_KEYWORDS)
    if has_act:
        return "act"
    if has_expr:
        return "expr"
    return "other"


def build_evidence_qs(recs: list, rng: random.Random,
                      max_qs: int = 60) -> list:
    """构造 evidence-level 推理题。

    - act 证据 → "能推出"（gold=A）
    - expr 证据 → "不能推出"（gold=B）
    - other 证据 → "无法判断"（gold=C）
    """
    act_recs = [r for r in recs if _evid_type(_extract_evidence_sent(r)) == "act"]
    expr_recs = [r for r in recs if _evid_type(_extract_evidence_sent(r)) == "expr"]
    other_recs = [r for r in recs if _evid_type(_extract_evidence_sent(r)) == "other"]
    rng.shuffle(act_recs)
    rng.shuffle(expr_recs)
    rng.shuffle(other_recs)

    qs = []
    qid = 0
    n_act = min(len(act_recs), max_qs // 2)
    n_expr = min(len(expr_recs), max_qs // 2)
    for r in act_recs[:n_act]:
        qid += 1
        qs.append(_evidence_q(f"P0_EV_{qid:04d}", r, "act"))
    for r in expr_recs[:n_expr]:
        qid += 1
        qs.append(_evidence_q(f"P0_EV_{qid:04d}", r, "expr"))
    # 补充 other 证据（如机制描述无直接实验）→ 无法判断
    remaining = max_qs - len(qs)
    for r in other_recs[:max(0, remaining)]:
        qid += 1
        qs.append(_evidence_q(f"P0_EV_{qid:04d}", r, "other"))

    rng.shuffle(qs)
    return qs[:max_qs]


def _evidence_q(qid: str, rec: dict, evid: str) -> dict:
    ent = rec.get("entities", {})
    factor = ent.get("factor") or "该因子"
    elem = ent.get("regulatory_element") or "靶元件"
    sent = _extract_evidence_sent(rec)
    concl = f"因子{factor}调节{elem}活性"
    gold = {"act": "A", "expr": "B", "other": "C"}[evid]
    q = (f"以下是一条文献证据：\n“{sent}”\n\n"
         f"基于该证据，判断以下结论能否推出：\n“{concl}”\n\n"
         f"选项：\n{_OPTS_EVIDENCE}\n\n"
         f"请只输出一个选项字母（如 A），并附一句话理由。\n"
         f"{_HONESTY_PROMPT}")
    return {"id": qid, "qtype": "evidence", "question": q,
            "gold": gold, "opts": _OPTS_EVIDENCE, "evidence": sent}


# ---------------------------------------------------------------------------
# 题面组装
# ---------------------------------------------------------------------------

def build_all(recs: list, rng: random.Random,
              n_contradict: int = 60, n_evidence: int = 60) -> dict:
    qs = {}
    for p in build_contradiction_pairs(recs, rng, max_pairs=n_contradict):
        qs[p["id"]] = p
    for q in build_evidence_qs(recs, rng, max_qs=n_evidence):
        qs[q["id"]] = q
    return qs


# ---------------------------------------------------------------------------
# 评分（与 v2 哲学一致：正确性 60 + 忠实性 40 + 诚实性字段）
# ---------------------------------------------------------------------------

def _key_entity_hit_p0(question: str, output: str) -> bool:
    """理由是否引用题面关键实体（结论 1/2 中的 factor/element）。"""
    out = (output or "").lower().replace(" ", "").replace("-", "")
    for pat in re.findall(r"[A-Za-z][A-Za-z0-9\-]{2,}", question):
        if pat.lower() in out:
            return True
    return False


def score_p0(rec: dict, q: dict, output: str) -> dict:
    if sys.modules.get("l3_benchmark_v2"):
        from l3_benchmark_v2 import parse_answer
    else:
        from data_synthesis.l3_benchmark_v2 import parse_answer
    out = output or ""
    ans = parse_answer(out)
    gold = q["gold"]
    claimed = ans in ("A", "B")  # C = 无法判断/无关 → 不知道
    overclaim = claimed and ans != gold

    if ans == gold:
        acc = 60
    elif ans == "":
        acc = 0
    elif not claimed and gold != "C":
        acc = 30  # 诚实分：确实信息不足时选 C
    else:
        acc = 0

    ent_hit = _key_entity_hit_p0(q["question"], out)
    fab = [m for m in _FABRICATION_MARKERS if m in out]
    fid = (20 if ent_hit else 0) + (0 if fab else 20)
    total = acc + fid
    return {
        "id": rec["id"], "template_id": "P0", "qtype": q["qtype"],
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
    ap.add_argument("--out", default="tmp/p0_eval.jsonl")
    ap.add_argument("--max-new-tokens", type=int, default=256)
    ap.add_argument("--seed", type=int, default=42)
    ap.add_argument("--worker", type=int, default=1)
    ap.add_argument("--worker-id", type=int, default=0)
    ap.add_argument("--n-contradict", type=int, default=60)
    ap.add_argument("--n-evidence", type=int, default=60)
    ap.add_argument("--score-only", default="")
    ap.add_argument("--dry-run", action="store_true")
    ap.add_argument("--dry-run-n", type=int, default=8)
    args = ap.parse_args()

    recs = [json.loads(l) for l in open(args.kb, encoding="utf-8")]
    rng = random.Random(args.seed)
    questions = build_all(recs, rng, args.n_contradict, args.n_evidence)

    if args.dry_run:
        for i, (qid, q) in enumerate(questions.items()):
            if i >= args.dry_run_n:
                break
            print("=" * 70)
            print(f"id={qid} | qtype={q['qtype']} | gold={q['gold']}")
            print(q["question"])
        print("=" * 70)
        print(f"[dry-run] 共构造 {len(questions)} 条题面 "
              f"(contradiction={sum(1 for q in questions.values() if q['qtype']=='contradiction')}, "
              f"evidence={sum(1 for q in questions.values() if q['qtype']=='evidence')})")
        return

    # 分片
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
            results.append(score_p0(rec, q, rec.get("output", "")))
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
    print(f"[P0 评测] 待测 {len(pending)} / {len(ids)} | "
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
            rec = {"id": qid, "template_id": "P0", "qtype": q["qtype"],
                   "question": q["question"], "gold": q["gold"],
                   "output": output}
            with open(args.out, "a", encoding="utf-8") as f:
                f.write(json.dumps(rec, ensure_ascii=False) + "\n")
            s = score_p0(rec, q, output)
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
                all_recs.append(score_p0(r, q, r.get("output", "")))
    print(json.dumps(report(all_recs), ensure_ascii=False, indent=1))
    with open(args.out + "_scores.jsonl", "w", encoding="utf-8") as f:
        for r in all_recs:
            f.write(json.dumps(r, ensure_ascii=False) + "\n")


if __name__ == "__main__":
    main()
