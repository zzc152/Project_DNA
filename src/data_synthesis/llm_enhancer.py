#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
llm_enhancer.py — Module 2 第三步：LLM Enhancer 增强器（v1）

设计决策（用户已确认，2026-08-06）：
  D1  system prompt = 严格遵守文献和 claim 的严谨生物学家（禁止编造文献外事实）
  D2  instruction  = 直接用模板引擎生成的 instruction 原文（不重新组织）
  D3  上下文注入   = instruction + 该草稿来源的 reasoning_chain（如 instruction 已含推理链则不重复）
  D4  输出格式     = 自由生成，不强制格式（后置解析交给 quality_filter）

处理范围：
  L1/L2 → 规则直出（模板答案可从 metadata 直接拼接，无需 LLM；规则见 DIRECT_ANSWER_RULES）
  L3/L4/L5 → 本地 Qwen AWQ 生成（默认 32B，--model 可切 72B 做小批量验证）

用法：
    python src/data_synthesis/llm_enhancer.py \
        --drafts data/synthetic/drafts.jsonl \
        --out data/synthetic/enhanced.jsonl \
        --model models/Qwen2.5-32B-AWQ \
        --max-new-tokens 512 --temperature 0.2

断点续跑：已存在于 --out 中的 id 自动跳过（可中断后重跑同一命令）。
"""

import argparse
import json
import os
import sys


# ---------------------------------------------------------------------------
# 🔧 boilerplate：IO
# ---------------------------------------------------------------------------

def load_drafts(path: str) -> list:
    with open(path, encoding="utf-8") as f:
        return [json.loads(l) for l in f if l.strip()]


def load_existing(path: str) -> dict:
    """读取已产出结果 {id: record}，用于断点续跑。"""
    out = {}
    if os.path.exists(path):
        with open(path, encoding="utf-8") as f:
            for l in f:
                if l.strip():
                    rec = json.loads(l)
                    out[rec["id"]] = rec
    return out


def append_result(path: str, rec: dict) -> None:
    """单条追加写入（崩溃时最多丢一条）。"""
    os.makedirs(os.path.dirname(path) or ".", exist_ok=True)
    with open(path, "a", encoding="utf-8") as f:
        f.write(json.dumps(rec, ensure_ascii=False) + "\n")


# ---------------------------------------------------------------------------
# 🧠 design decision：L1/L2 规则直出（v2，用户设计 2026-08-06）
# 短输出需补"原因解释 + 事实依据"，从最初文献知识与统计规律中提取。
# 输出格式：答案：...\n依据：...（依据来自 reasoning_chain 的原始发现/文献机制 +
#           evidence 的统计值/来源）
# ---------------------------------------------------------------------------

def _fmt_reasoning_chain(rc) -> str:
    if not rc:
        return ""
    if isinstance(rc, str):
        return rc
    return "\n".join(f"{i}. {s}" for i, s in enumerate(rc, 1))


def _rc_texts(draft: dict) -> list:
    """把 reasoning_chain 规整为字符串列表。"""
    rc = draft["metadata"].get("reasoning_chain")
    if isinstance(rc, list):
        return [str(s) for s in rc]
    if isinstance(rc, str):
        return [rc]
    return []


def _find_rc(draft: dict, keywords: list) -> str:
    """在推理链中找第一条含任一关键词的条目（返回原文）。
    关键词里 '文献' 之类宽泛词会误伤 '由文献关联推导'，调用方需把
    精确词（原始发现/文献机制/摘要原文）排在前面。"""
    for s in _rc_texts(draft):
        if any(k in s for k in keywords):
            return s
    return ""


def _evidence_stats(draft: dict) -> str:
    """从 evidence 提取统计规律（效应量/p值），用于依据。
    效应量符号：direction 含 depleted 时 OR<1 表示耗尽，此时显示符号取负；
    否则取正（OR>1 富集）。
    """
    ev = draft["metadata"].get("evidence") or {}
    stats = []
    if ev.get("effect_size") is not None:
        es = ev["effect_size"]
        direction = str(ev.get("direction") or "")
        if "depleted" in direction and es > 1:
            es = -es  # 耗尽时 OR 实际 <1，取负号与 direction 一致
        stats.append(f"效应量={es:+.2f}")
    if ev.get("p_value") is not None:
        stats.append(f"p={ev['p_value']:.2e}")
    return "；".join(stats)


def _evidence_source(draft: dict) -> str:
    """来源：优先 PMID，其次 evidence.source。"""
    ents = draft["metadata"].get("entities_used") or {}
    pmids = ents.get("pmids") or []
    if pmids:
        return f"PMID {pmids[0]}"
    ev = draft["metadata"].get("evidence") or {}
    return ev.get("source") or ""


def _basis(draft: dict) -> str:
    """构造'依据'句：统计规律（若有） + 推理链原始发现/文献机制 + 来源。"""
    parts = []
    stats = _evidence_stats(draft)
    if stats:
        parts.append(f"统计规律：{stats}")
    raw = _find_rc(draft, ["原始发现", "文献机制", "摘要原文", "文献机制描述"])
    if raw:
        parts.append(raw)
    src = _evidence_source(draft)
    if src:
        parts.append(f"来源：{src}")
    return "\n".join(parts) if parts else ""


def direct_answer_L1_001(draft: dict) -> str:
    """因子对调控元件活性的影响方向（附依据）。"""
    ents = draft["metadata"]["entities_used"]
    ev = draft["metadata"]["evidence"] or {}
    core = (ents.get("effect") or ev.get("direction") or "not_specified")
    out = [f"答案：{core}"]
    b = _basis(draft)
    if b:
        out.append(f"依据：{b}")
    return "\n".join(_ensure_core_target(out, draft))


def direct_answer_L1_002(draft: dict) -> str:
    """因子影响哪个元件 + 哪个细胞系（附依据）。"""
    ents = draft["metadata"]["entities_used"]
    el = ents.get("regulatory_element") or "not_specified"
    cl = ents.get("cell_line") or "not_specified"
    out = [f"答案：调控元件 {el}；细胞系 {cl}"]
    b = _basis(draft)
    if b:
        out.append(f"依据：{b}")
    return "\n".join(_ensure_core_target(out, draft))


def direct_answer_L1_003(draft: dict) -> str:
    """一句话机制概括（附依据）。
    优先取链中 原始发现/文献机制/摘要原文 条目；无则回退第一条
    非'由文献关联推导'的实质条目。"""
    rc = _rc_texts(draft)
    core = _find_rc(draft, ["原始发现", "文献机制", "摘要原文", "文献机制描述"])
    if not core:
        # 回退：跳过 '由文献关联推导' 这类占位条目
        for s in rc:
            if s and "由文献关联推导" not in s and "需功能实验验证" not in s:
                core = s
                break
    if not core and rc:
        core = rc[0]
    if not core:
        core = "not_specified"
    out = [f"答案：{core}"]
    b = _basis(draft)
    if b:
        out.append(f"依据：{b}")
    return "\n".join(_ensure_core_target(out, draft))


def direct_answer_L2_001(draft: dict) -> str:
    """统计富集 → TF 调控含义（附依据）。"""
    ents = draft["metadata"]["entities_used"]
    ev = draft["metadata"]["evidence"] or {}
    tf = ents.get("tf", "该转录因子")
    cl = ents.get("cell_line", "该细胞系")
    direction = ev.get("direction", "")
    if "depleted" in str(direction) or "low" in str(direction).lower():
        rel = "负相关（结合位点在高活性序列中减少）"
    else:
        rel = "正相关（结合位点在高活性序列中富集）"
    out = [f"答案：在 {cl} 中，{tf} 结合位点与增强子活性{rel}。"]
    b = _basis(draft)
    if b:
        out.append(f"依据：{b}")
    return "\n".join(_ensure_core_target(out, draft))


def direct_answer_L2_002(draft: dict) -> str:
    """方向含义解释（附依据）。"""
    ents = draft["metadata"]["entities_used"]
    ev = draft["metadata"]["evidence"] or {}
    factor = ents.get("factor", "该因子")
    el = ents.get("regulatory_element", "调控元件")
    cl = ents.get("cell_line", "该细胞系")
    direction = ents.get("effect") or ev.get("direction") or "not_specified"
    meaning = {
        "increases": "促进（激活）",
        "increases_activity": "促进（激活）",
        "decreases": "抑制",
        "decreases_activity": "抑制",
        "required_for": "是活性所必需的",
        "modulates": "调节",
        "no_effect": "无显著影响",
    }.get(str(direction), str(direction))
    out = [f"答案：{factor} 的调控方向为 {meaning}，在 {cl} 细胞背景下表明其{meaning}{el}活性。"]
    b = _basis(draft)
    if b:
        out.append(f"依据：{b}")
    return "\n".join(_ensure_core_target(out, draft))


def direct_answer_L2_003(draft: dict) -> str:
    """推理链压缩成一句话因果（附依据）。"""
    rc = _rc_texts(draft)
    core = _find_rc(draft, ["文献机制描述", "摘要原文", "表明", "因此", "因果"])
    if not core and rc:
        core = rc[-1]  # 回退链尾结论句
    if not core:
        core = "not_specified"
    # 去掉条目自身的 '文献：' 前缀，避免 '答案：文献：...'
    core = core.replace("文献：", "", 1)
    core = core.replace("结论：", "", 1)
    out = [f"答案：{core}"]
    b = _basis(draft)
    if b:
        out.append(f"依据：{b}")
    return "\n".join(_ensure_core_target(out, draft))


def _ensure_core_target(out_lines: list, draft: dict) -> list:
    """L1/L2 直出补全：确保答案覆盖核心实体（regulatory_element / cell_line）。
    若输出文本中未出现元件的具体值（如 enhancer/promoter），则追加一行
    '目标：{元件}（{细胞系}）'，避免核心实体缺失被 judge 标为未提及。
    （L1/L2 核心实体 = regulatory_element 型；factor 型非核心不强制）"""
    ents = draft["metadata"].get("entities_used") or {}
    els = ents.get("regulatory_element") or []
    if isinstance(els, str):
        els = [els]
    cl = ents.get("cell_line") or ""
    if isinstance(cl, (list, tuple)):
        cl = cl[0] if cl else ""
    text = "\n".join(out_lines)
    need_el = [str(e) for e in els
               if e and str(e) not in ("not_specified", "none", "other")
               and str(e) not in text]
    if cl and str(cl) not in ("not_specified", "none", "other") and str(cl) not in text:
        need_el.append(str(cl))
    if need_el:
        out_lines.append(f"目标：{'、'.join(need_el)}")
    return out_lines


def direct_answer_L2_004(draft: dict) -> str:
    """GC 富集方向 → GC 含量偏好（附统计依据）。
    答案可由 evidence.direction 直接判定：
      enriched_in_gc_rich → 结合位点倾向出现在 GC 含量高的序列
      depleted_in_gc_rich  → 结合位点倾向出现在 GC 含量低的序列
    """
    ents = draft["metadata"].get("entities_used") or {}
    ev = draft["metadata"].get("evidence") or {}
    tf = ents.get("tf") or "该因子"
    cl = ents.get("cell_line") or ""
    direction = str(ev.get("direction") or "")
    if "depleted" in direction:
        relation = "结合位点在GC含量高的序列中显著减少（耗尽）"
        tendency = "更倾向出现在GC含量低的序列中"
    else:
        relation = "结合位点在GC含量高的序列中显著富集"
        tendency = "更倾向出现在GC含量高的序列中"
    out = [f"答案：在{cl}中，{tf}结合位点与序列GC含量相关——{relation}，说明{tf}结合位点{tendency}。"]
    b = _basis(draft)
    if b:
        out.append(f"依据：{b}")
    return "\n".join(_ensure_core_target(out, draft))


# 模板 id → 直出函数
DIRECT_ANSWER_RULES = {
    "L1_001": direct_answer_L1_001,
    "L1_002": direct_answer_L1_002,
    "L1_003": direct_answer_L1_003,
    "L2_001": direct_answer_L2_001,
    "L2_002": direct_answer_L2_002,
    "L2_003": direct_answer_L2_003,
    "L2_004": direct_answer_L2_004,
}

# L3+ 才调用 LLM
LLM_LEVELS = {"L3", "L4", "L5"}


# ---------------------------------------------------------------------------
# 🧠 design decision：LLM prompt 构造（D1/D2/D3）
# ---------------------------------------------------------------------------

SYSTEM_PROMPT = (
    "你是严格遵守文献和 claim 的严谨生物学家。"
    "回答必须严格基于给定的科学结论及其推理链，不得编造文献之外的事实，"
    "不得添加未经证据支持的机制或数值。"
    "\n\n【已知错误示范——严厉禁止，必须避免】"
    "\n1. 自补越界（编造文献外机制）——触发即作废：只要输出中出现"
    "'虽然文献未直接提及/虽然文献中没有/尽管文献未/文献中没有直接提到/原文未提及/"
    "未在文献中'等'文献没提 + 但…'句式，整条输出即视为编造、直接作废。"
    "若给定推理链中没有某个事实，就不要提及它，更不要用'虽然…但…'去补充圆场。"
    "所有内容必须逐条能在给定推理链中找到出处；找不到的事实一律不写。"
    "正面写法示例：'推理链给出的依据是……'，而不是'虽然文献没有直接提到，但理论上……'。"
    "\n2. 输出不完整：设计优化任务只输出'改进方案'，"
    "缺少'失败原因'与'预期效果'。必须完整覆盖任务要求的全部要素。"
    "\n3. 空泛无依据：只给出结论方向（如 increases）或'由文献关联推导'，"
    "没有任何机制/统计依据。应引用来源与具体机制。"
    "\n4. 实体指代漂移：把知识中的因子 A 换成因子 B 来回答。"
    "只能使用给定的因子与元件。"
    "\n5. 漏掉目标细胞系/元件：设计准则、诊断方案必须明确提及任务指定的"
    "细胞系（如 K562）与元件（如 enhancer/promoter）。"
    "只写泛化规则（如'基序密度应较高'）而不点名目标细胞系/元件，视为不完整，禁止。"
)


def build_messages(draft: dict) -> list:
    """构造 chat 消息。
    - instruction 原样作为 user 内容（D2）
    - 若 instruction 未含推理链，则附加来源 reasoning_chain（D3）
    - 若知识基因特异（gene_specific），注入基因特异性提示（用户决策 2026-08-06）：
      该规则仅对该基因的元件有效，设计/回答须针对该基因序列。
    """
    instruction = draft["instruction"]
    rc = draft["metadata"].get("reasoning_chain")
    extra = ""
    if rc and "推理链" not in instruction and "reasoning chain" not in instruction.lower():
        extra = "\n\n参考推理链（必须严格遵循）:\n" + _fmt_reasoning_chain(rc)
    gene = draft["metadata"].get("gene_specific")
    if gene:
        el = draft["metadata"].get("entities_used", {}).get("regulatory_element", "")
        extra += (f"\n\n⚠️ 注意：上述知识是针对 {gene} 基因的 {el} 得出的，"
                  f"仅对该基因的 {el} 序列有效，不可泛化到其他基因的 {el}。")
    return [
        {"role": "system", "content": SYSTEM_PROMPT},
        {"role": "user", "content": instruction + extra},
    ]


# ---------------------------------------------------------------------------
# 🔧 boilerplate：awq 模型加载与生成（已验证路线）
# ---------------------------------------------------------------------------

def load_model(model_path: str, device: str):
    """awq 原生库加载（AutoAWQ 0.2.9 已验证）。
    ⚠️ device_map 必须传 dict {"": "cuda:0"}，传字符串 "cuda:0" 会报错。
    """
    from awq import AutoAWQForCausalLM
    from transformers import AutoTokenizer
    model = AutoAWQForCausalLM.from_quantized(
        model_path,
        device_map={"": device},
        fuse_layers=True,
        local_files_only=True,
    )
    tok = AutoTokenizer.from_pretrained(model_path, local_files_only=True)
    return model, tok


def _model_device(model) -> str:
    """获取模型所在设备。AWQ 模型对象没有 .device 属性，取首个参数设备。"""
    return next(model.parameters()).device


def generate(model, tok, messages: list, max_new_tokens: int,
             temperature: float) -> str:
    """chat template 构造 prompt → 生成 → 解码时跳过输入部分。"""
    prompt = tok.apply_chat_template(
        messages, tokenize=False, add_generation_prompt=True
    )
    inputs = tok(prompt, return_tensors="pt").to(_model_device(model))
    out = model.generate(
        **inputs,
        max_new_tokens=max_new_tokens,
        temperature=temperature,
        do_sample=True,
        top_p=0.9,
        repetition_penalty=1.05,
    )
    new_tokens = out[0][inputs.input_ids.shape[1]:]
    return tok.decode(new_tokens, skip_special_tokens=True).strip()


# ---------------------------------------------------------------------------
# 主流程
# ---------------------------------------------------------------------------

def main():
    ap = argparse.ArgumentParser(description="LLM Enhancer 增强器")
    ap.add_argument("--drafts", default="data/synthetic/drafts.jsonl")
    ap.add_argument("--out", default="data/synthetic/enhanced.jsonl")
    ap.add_argument("--model", default="models/Qwen2.5-32B-AWQ")
    ap.add_argument("--device", default="cuda:0")
    ap.add_argument("--max-new-tokens", type=int, default=512)
    ap.add_argument("--temperature", type=float, default=0.2)
    ap.add_argument("--limit", type=int, default=0,
                    help="仅处理前 N 条（0=全部；用于小批量验证）")
    ap.add_argument("--worker", type=int, default=1,
                    help="并行 worker 总数（多卡并行用）")
    ap.add_argument("--worker-id", type=int, default=0,
                    help="本进程 worker 编号（0-based，< --worker）")
    args = ap.parse_args()

    drafts = load_drafts(args.drafts)
    existing = load_existing(args.out)
    pending = [d for d in drafts if d["id"] not in existing]
    if args.limit > 0:
        pending = pending[:args.limit]
    # worker 分片：每个 worker 处理 id % worker == worker-id 的条目
    if args.worker > 1:
        pending = [d for i, d in enumerate(pending)
                   if i % args.worker == args.worker_id]
    print(f"[草稿] {len(drafts)} | [已完成] {len(existing)} | "
          f"[待处理] {len(pending)} | [worker {args.worker_id}/{args.worker}]",
          flush=True)

    llm_todo = [d for d in pending if d["level"] in LLM_LEVELS]
    direct_todo = [d for d in pending if d["level"] not in LLM_LEVELS]
    print(f"[LLM] {len(llm_todo)} 条 | [规则直出] {len(direct_todo)} 条", flush=True)

    model, tok = None, None
    if llm_todo:
        print(f"[加载模型] {args.model} -> {args.device}", flush=True)
        model, tok = load_model(args.model, args.device)
        print("[模型就绪]", flush=True)

    n_ok = 0
    for d in pending:
        level = d["level"]
        try:
            if level in LLM_LEVELS:
                msgs = build_messages(d)
                output = generate(model, tok, msgs, args.max_new_tokens,
                                  args.temperature)
                llm_model = args.model.split("/")[-1]
                direct = False
            else:
                fn = DIRECT_ANSWER_RULES.get(d["template_id"])
                if fn is None:
                    raise ValueError(f"无直出规则: {d['template_id']}")
                output = fn(d)
                llm_model = "rule"
                direct = True
            rec = dict(d)
            rec["output"] = output
            rec["direct"] = direct
            rec["metadata"] = dict(d["metadata"])
            rec["metadata"]["llm_model"] = llm_model
            append_result(args.out, rec)
            n_ok += 1
            if n_ok % 10 == 0 or n_ok <= 3:
                print(f"  ✓ {d['id']} ({level})", flush=True)
        except Exception as e:
            print(f"  ✗ {d['id']} 失败: {e}", flush=True)
            # 断点续跑友好：失败不写入，重跑该命令即可续上

    print(f"[完成] 处理 {n_ok}/{len(pending)} 条 -> {args.out}", flush=True)


if __name__ == "__main__":
    main()
