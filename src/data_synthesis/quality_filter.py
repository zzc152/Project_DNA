#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
quality_filter.py — Module 2 第四步：质量过滤器（v1 规则版）

目的：对 llm_enhancer 产出的增强记录做规则级质量校验，标记可疑条目。
设计（用户已定方向，2026-08-06）：
  ① 规则校验：结构完整性（L4 须有 3+ 条准则、L5 须有失败原因/改进方案/预期效果）
  ② 知识一致性：output 中引用的 PMID 必须在来源知识库 pmids 内（防幻觉引用新文献）；
                 实体必须落在 instruction/metadata 实体集合内（防跑题/替换）
  ③ 幻觉标记：output 中出现"文献未提及"类自注 → 标记为模型自补越界（违反 D1 不编造）

校验通过 → passed=True；有 issue → passed=False 并附 issues 列表。
评分：规则版 0-10（结构完整性为主），LLM 自检评分留作后续可选步骤（--llm-score）。

用法：
    python src/data_synthesis/quality_filter.py \
        --in data/synthetic/enhanced.jsonl \
        --out data/synthetic/filtered.jsonl
    # --report 输出校验报告（默认打印摘要）
"""

import argparse
import json
import re

# LLM judge 用：延迟 import（llm_enhancer 会拉 torch/awq，非 judge 模式不加载）
_generate_fn = None


def _get_generate():
    """惰性导入 llm_enhancer.generate（避免普通校验加载重依赖）。"""
    global _generate_fn
    if _generate_fn is None:
        from llm_enhancer import generate
        _generate_fn = generate
    return _generate_fn


# ---------------------------------------------------------------------------
# 🔧 boilerplate：IO
# ---------------------------------------------------------------------------

def load_records(path: str) -> list:
    with open(path, encoding="utf-8") as f:
        return [json.loads(l) for l in f if l.strip()]


def dump_records(recs: list, path: str) -> None:
    with open(path, "w", encoding="utf-8") as f:
        for r in recs:
            f.write(json.dumps(r, ensure_ascii=False) + "\n")


# ---------------------------------------------------------------------------
# 🧠 design decision：知识一致性（防幻觉引用 + 防跑题）
# ---------------------------------------------------------------------------

_PMID_RE = re.compile(r"\b\d{7,8}\b")


def _known_pmids(rec: dict) -> set:
    """该记录来源知识库的合法 PMID 集合（entities_used.pmids + reasoning_chain 内提取）。"""
    known = set()
    ents = rec.get("metadata", {}).get("entities_used", {}) or {}
    for p in ents.get("pmids", []) or []:
        if p:
            known.add(str(p))
    rc = rec.get("metadata", {}).get("reasoning_chain")
    if isinstance(rc, list):
        for s in rc:
            for m in _PMID_RE.findall(str(s)):
                known.add(m)
    return known


def _cited_pmids(text: str) -> set:
    return set(_PMID_RE.findall(text))


def _entity_names(rec: dict) -> set:
    """该记录涉及的**内容实体**名集合（用于防跑题检查）。
    只检查内容实体字段：factor/tf/regulatory_element/cell_line/gene_specific。
    排除结构化控制词：effect（increases/decreases）、factor_type（TF/other）、
    motif（MAxxxx 技术标识）、not_specified 占位值。"""
    names = set()
    ents = rec.get("metadata", {}).get("entities_used", {}) or {}
    # 内容实体字段（白名单）
    for k in ("factor", "tf", "regulatory_element", "cell_line"):
        v = ents.get(k)
        if isinstance(v, str) and v and len(v) <= 40:
            v = v.strip()
            # 排除控制词/占位值
            if v.lower() in ("not_specified", "none", "other", "tf"):
                continue
            names.add(v)
    gene = rec.get("metadata", {}).get("gene_specific")
    if gene:
        names.add(gene)
    return names


# ---------------------------------------------------------------------------
# 🧠 design decision：实体一致性判断（两层，用户设计 2026-08-06）
# 第一层：精确子串匹配（output 中出现实体原名 → 提及）。
# 第二层：未精确匹配的实体不直接判错，引入 LLM 判断——
#   output 中的相关表述是否与原始知识表达一致（同一意思、无精度损失、
#   无指代不清、无指代其他事物）。LLM 判"一致"则通过，判"不一致"才是问题。
# ---------------------------------------------------------------------------

JUDGE_SYSTEM = (
    "你是严格的知识一致性审查员。判断模型回答中是否出现并一致地表达了"
    "给定的原始知识实体。"
    "要求：1) 回答中出现该实体（同义/同指）且表达一致、无精度损失 → 一致；"
    "2) 回答中出现了相关内容但指代不清、精度下降、或指代了另一个事物 → 不一致；"
    "3) 回答完全没有提到该实体 → 未提及。"
    "只输出：一致 或 不一致 或 未提及。"
)


def judge_entities_consistency(model, tok, output: str,
                               entities: list) -> dict:
    """用 LLM 判断 entities 中每个实体在 output 中如何被表达。
    返回 {entity: True/False/None}：
      True=一致表达（不算问题）
      False=提到但表达不一致/指代漂移（问题）
      None=完全未提及（由 check_record 按层级完整性规则处理）"""
    result = {}
    if not entities:
        return result
    prompt = (
        "原始知识实体：\n"
        + "\n".join(f"- {e}" for e in entities)
        + f"\n\n模型回答：\n{output[:1500]}\n\n"
        "请逐个判断每个实体在回答中的出现情况（一致 / 不一致 / 未提及）。"
        "对每个实体输出一行：\n"
        "<实体名>：一致 或 不一致 或 未提及"
    )
    msgs = [
        {"role": "system", "content": JUDGE_SYSTEM},
        {"role": "user", "content": prompt},
    ]
    try:
        gen = _get_generate()
        # temperature 必须 >0（do_sample=True 时）；0.01 ≈ greedy
        judge_out = gen(model, tok, msgs, 256, 0.01)
    except Exception as e:  # LLM 失败则保守视为不一致
        print(f"  [judge 异常] {e}", flush=True)
        for e_ in entities:
            result[e_] = False
        return result
    import re
    for ent in entities:
        # 找 "<ent>：一致/不一致/未提及" 行
        m = re.search(
            rf"{re.escape(ent)}[：:]\s*(一致|不一致|未提及)", judge_out
        )
        if not m:
            result[ent] = None  # LLM 未明确判定 → 交给完整性规则
            continue
        tag = m.group(1)
        if tag == "一致":
            result[ent] = True
        elif tag == "不一致":
            result[ent] = False
        else:
            result[ent] = None
    return result


# 模型自补越界信号（违反 D1：不得编造文献之外事实）
_FABRICATION_MARKERS = [
    "文献未提及", "文献未直接提及", "文献中没有", "原文未提及",
    "虽然文献未", "尽管文献未",
]
# 合理但需验证的表述（弱信号，仅记录不否决）
_VERIFY_MARKERS = ["需通过功能实验验证", "需要功能实验验证", "需实验验证", "尚未验证"]


# ---------------------------------------------------------------------------
# 🧠 design decision：层级结构规则（L4/L5 新设计校验）
# ---------------------------------------------------------------------------

def _count_list_items(text: str) -> int:
    """统计输出中的编号条目数（1. 2. 3. 或 **1. 或 准则1 等）。"""
    n = len(re.findall(r"(?:^|\n)\s*\*{0,2}\s*\d+[\.、）)]", text))
    n += len(re.findall(r"准则\s*\d", text))
    return n


def check_record(rec: dict, judge=None) -> dict:
    """对单条记录做规则校验。
    judge: 可选函数 judge(output, entities) -> {entity: True/False}，
           用于对未精确匹配的实体做 LLM 一致性复核。
    """
    rec_id = rec["id"]
    level = rec["level"]
    output = rec.get("output") or ""
    issues = []
    warnings = []

    # --- ① 通用结构 ---
    if not output.strip():
        issues.append("output 为空")
        return {"id": rec_id, "level": level, "passed": False,
                "score": 0, "issues": issues, "warnings": warnings}
    if len(output) < 30:
        issues.append(f"output 过短（{len(output)} 字符）")
    ph = re.findall(r"\{[a-z_]+\}", output)
    if ph:
        issues.append(f"占位符残留: {ph[:3]}")

    # --- ② 知识一致性：PMID 引用必须落在来源知识库内 ---
    cited = _cited_pmids(output)
    known = _known_pmids(rec)
    fake_pmids = cited - known
    if fake_pmids:
        issues.append(f"引用了知识库之外的 PMID: {sorted(fake_pmids)}")

    # --- ③ 幻觉标记：模型自补越界声明 ---
    for mk in _FABRICATION_MARKERS:
        if mk in output:
            idx = output.find(mk)
            snippet = output[max(0, idx - 20): idx + 25].replace("\n", " ")
            issues.append(f"模型自补越界（{mk}）: ...{snippet}...")
            break
    for vk in _VERIFY_MARKERS:
        if vk in output:
            warnings.append(f"含待验证表述（{vk}）")

    # --- ④ 实体一致性（精确子串 + 可选 LLM 复核；分层核心实体） ---
    # 核心实体（按层级）：答案应涵盖的对象；缺失/表达漂移 → issue
    #   L1/L2 短答案：核心 = regulatory_element（effect/cell_line 由模板语义决定，
    #                 factor 在问题不在答案，非核心）
    #   L3/L4/L5 长答案：核心 = factor + cell_line
    ents = _entity_names(rec)
    if ents:
        # 拆出核心实体（按层级）与非核心实体（仅查表达一致性，不查缺失）
        if level in ("L1", "L2"):
            core_ents = {e for e in ents if e in ("enhancer", "promoter",
                                                  "cis-regulatory element",
                                                  "super-enhancer",
                                                  "TF binding motifs")}
        else:
            core_ents = set(ents)  # L3+ 全部内容实体都是核心
        optional_ents = ents - core_ents

        def _check_missing(ent: str) -> bool:
            """未精确命中时，该实体是否构成'缺失/漂移'问题。
            精确命中 → False（无问题）。"""
            if str(ent) in output:
                return False
            return True

        # 第一层：精确子串
        missing_core = [e for e in core_ents if _check_missing(e)]
        missing_opt = [e for e in optional_ents if _check_missing(e)]

        if judge is not None:
            # 第二层：LLM 复核所有未精确命中的实体（核心+可选）
            all_unmatched = missing_core + missing_opt
            if all_unmatched:
                verdict = judge(output, all_unmatched)
                # True=一致（OK）；False=提到但漂移（issue）；
                # None=未提及 → 核心实体算 issue，可选实体忽略
                bad_core = [e for e in missing_core
                            if verdict.get(e, None) is False]
                bad_opt = [e for e in missing_opt
                           if verdict.get(e, None) is False]
                miss_core = [e for e in missing_core
                             if verdict.get(e, None) is None]
                still_bad = bad_core + bad_opt + miss_core
            else:
                still_bad = []
        else:
            # 无 judge：核心实体缺失 → issue；可选实体缺失 → 忽略（防误伤）
            still_bad = list(missing_core)
        if still_bad:
            issues.append(f"output 未一致提及实体: {still_bad[:5]}")

    # --- ⑤ 层级结构规则 ---
    if level == "L4":
        n_items = _count_list_items(output)
        # 有 ≥3 个编号条目即视为准则形式（不再强制"准则/条件"字样，避免误报）
        if n_items < 3:
            if "准则" not in output and "条件" not in output:
                issues.append("L4 缺少'准则'/'条件'输出形式")
            issues.append(f"L4 设计准则条目过少（{n_items} < 3）")
    elif level == "L5":
        # 按模板区分要素要求（L5_001 失败诊断闭环三要素；
        # L5_002 迭代优化两要素：改进方案+验证实验）
        tpl = rec.get("template_id", "")
        if "L5_002" in tpl:
            for need in ("改进方案", "验证实验"):
                if need not in output:
                    issues.append(f"L5 缺少要素: {need}")
        else:
            for need in ("失败原因", "改进方案", "预期效果"):
                if need not in output:
                    issues.append(f"L5 缺少要素: {need}")
    elif level == "L3":
        if "PMID" not in output and "来源" not in output:
            warnings.append("L3 未标注来源/PMID")

    # --- 评分（0-10，规则版） ---
    score = 10
    score -= min(5, 3 * len(issues))          # 每个 issue 最多扣 3 分
    score -= min(2, len(warnings))            # 每个 warning 扣 1 分
    score = max(0, score)

    return {"id": rec_id, "level": level, "passed": len(issues) == 0,
            "score": score, "issues": issues, "warnings": warnings}


# ---------------------------------------------------------------------------
# 主流程
# ---------------------------------------------------------------------------

def main():
    ap = argparse.ArgumentParser(description="Quality Filter 质量过滤器（规则版）")
    ap.add_argument("--in", dest="inp", default="data/synthetic/enhanced.jsonl")
    ap.add_argument("--out", default="data/synthetic/filtered.jsonl")
    ap.add_argument("--threshold", type=float, default=6.0,
                    help="评分阈值，低于此分的记录被标记丢弃（默认 6.0）")
    ap.add_argument("--ids", default="", help="仅校验指定 id（逗号分隔），空=全部")
    ap.add_argument("--llm-judge", action="store_true",
                    help="对未精确匹配的实体启用 LLM 一致性复核（较慢）")
    ap.add_argument("--model", default="models/Qwen2.5-32B-AWQ",
                    help="LLM 判断用的模型（--llm-judge 时）")
    ap.add_argument("--device", default="cuda:0",
                    help="LLM 判断用设备（--llm-judge 时）")
    args = ap.parse_args()

    recs = load_records(args.inp)
    if args.ids:
        want = set(x.strip() for x in args.ids.split(",") if x.strip())
        recs = [r for r in recs if r["id"] in want]

    # 可选：LLM 一致性判断器
    judge = None
    if args.llm_judge:
        from llm_enhancer import load_model, generate
        print(f"[加载判断模型] {args.model} -> {args.device}", flush=True)
        model, tok = load_model(args.model, args.device)
        print("[判断模型就绪]", flush=True)

        def judge(output: str, entities: list) -> dict:
            return judge_entities_consistency(model, tok, output, entities)

    results = []
    kept = []
    dropped = []
    for r in recs:
        chk = check_record(r, judge=judge)
        results.append(chk)
        r["metadata"]["quality_score"] = chk["score"]
        r["metadata"]["quality_passed"] = chk["passed"]
        r["metadata"]["quality_issues"] = chk["issues"]
        (kept if chk["passed"] else dropped).append(r)

    # 输出报告
    print(f"[校验] {len(results)} 条 | 通过 {len(kept)} | 标记 {len(dropped)}")
    for chk in sorted(results, key=lambda c: c["score"]):
        flag = "✓" if chk["passed"] else "✗"
        print(f"  {flag} {chk['id']:12} score={chk['score']} "
              f"issues={len(chk['issues'])} warnings={len(chk['warnings'])}")
        for i in chk["issues"]:
            print(f"      - {i}")
        for w in chk["warnings"]:
            print(f"      ~ {w}")

    if args.out:
        dump_records(kept + dropped, args.out)
        print(f"[写出] {args.out} （{len(kept) + len(dropped)} 条，含标记）")


if __name__ == "__main__":
    main()
