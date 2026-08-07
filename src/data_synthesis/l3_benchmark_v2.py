# -*- coding: utf-8 -*-
"""L3 封闭式评测 v2 —— 修复 v1 的"题面泄露答案"问题。

v1 的问题：instruction 里直接给出 "SP1增强（increases）promoter活性"、"SP1 is a
positive regulator"、完整 reasoning_chain，模型只要复述题面就能拿高分，测的是
复述/格式服从而非推理能力；且 32B(87.9) < 7B(88.5) 暴露区分度不足。

v2 设计：封闭式问答。
- 题面不再包含方向标注、机制结论、推理链；只给「文献摘要原文句」（无方向词）
  或什么都不给。
- 模型输出「选项字母 + 一句话理由」，评分 = 答案正确性 + 理由忠实性。

三种题型（qtype）：
  rc   : 阅读理解——题面给摘要原文句，模型从原文推出方向（一半数据）
  know : 知识问答——题面不给原文，模型凭领域知识回答（另一半数据）
  data : 数据判读——题面给统计数字(OR/Cohen's d)，模型解读方向（L3_003 特异类 8 条）

模板差异：
  L3_001/002     : 因子→元件活性方向判断（选项 A增强/B减弱/C无影响/D方向不明/E无法判断）
  L3_003 普通类  : 结论判真（给构造结论，判断成立否；一半为干扰反结论）
  L3_003 特异类  : 数据判读（8 条）

v2.1 改进（抗蒙分 + 诚实性评测）：
1. 正负方向平衡（--balanced，默认开）：gold ∈ {A增强, B减弱} 的题配额各占一半，
   剔除多余 A 侧，防止模型无脑蒙 A 拿分（原数据 A:111 vs B:51）。
2. 诚实指令：题面统一写"知之为知之，不知为不知"，信息不足应选"无法判断"。
3. 不诚实扣分：统计"声称知道"（未选无法判断）的准确率 confident_acc；若 < 0.5
   （接近瞎蒙），按 (0.5 - confident_acc) × 40 全局扣分并报告 dishonesty 指标。

评分（满分 100）：
  答案正确性 60 : 答对 60；know 型答"无法判断"给 30（诚实，不瞎编）
  理由忠实性 40 : 理由引用题面关键实体/数字 20 + 无自补编造标记 20
  另报告 overclaim/confident_acc/penalty 等诚实性指标（见 report()）。
"""

import argparse
import json
import os
import random
import re
from collections import defaultdict

# ---------------------------------------------------------------------------
# 方向词表（v1 复用）
# ---------------------------------------------------------------------------

_DIRECTION_WORDS = {
    "increases": ["增加", "增强", "促进", "激活", "上调", "正相关", "正向", "increases", "positive",
                  "upregulat", "activat", "enhanc", "promot", "induc", "required", "必需", "必要",
                  "需要", "正向调节", "激活剂", "促进剂"],
    "decreases": ["减少", "减弱", "抑制", "下调", "负相关", "负向", "decreases", "negative",
                  "downregulat", "repress", "inhibit", "suppress", "减弱活性", "抑制因子", "负向调节"],
    "modulates": ["调节", "调控", "modulat", "调节活性", "影响"],
    "no_effect": ["无影响", "不影响", "no effect", "无关", "不改变"],
}

# 自补编造标记（复用 quality_filter）
_FABRICATION_MARKERS = [
    "文献未提及", "文献未直接提及", "文献中没有", "原文未提及",
    "虽然文献未", "尽管文献未",
]

# ---------------------------------------------------------------------------
# 实体提取 / 题面构造
# ---------------------------------------------------------------------------

_ABST_RE = re.compile(r"摘要原文证据[:：]\s*(.+)")
_MECH_RE = re.compile(r"文献机制描述[:：]\s*(.+)")

# 构造干扰结论用的方向词替换
_REV_POS = ["positive", "Positive", "activator", "Activator", "激活", "增强", "促进", "正向",
            "upregulat", "Upregulat", "induc", "Induc", "促进因子", "激活剂", "正相关"]
_REV_NEG = ["negative", "Negative", "repressor", "Repressor", "抑制", "减弱", "下调", "负向",
            "downregulat", "Downregulat", "抑制因子", "抑制剂", "负相关"]


def _norm_ent(v):
    """规范化实体名：None/'not_specified'/'other'/'tf' 视为缺失。"""
    if not v:
        return None
    if isinstance(v, (list, tuple)):
        v = ", ".join(str(x) for x in v)
    v = str(v).strip()
    if v.lower() in ("not_specified", "none", "other", "tf", "nan", ""):
        return None
    return v


def _abs_sentence(rec):
    """从 reasoning_chain 提取摘要原文句。"""
    rc = rec.get("metadata", {}).get("reasoning_chain", "")
    if isinstance(rc, list):
        rc = "\n".join(str(x) for x in rc)
    m = _ABST_RE.search(str(rc))
    return m.group(1).strip() if m else ""


def _mech_sentence(rec):
    """从 reasoning_chain 提取文献机制描述句（仅 L3_003 普通类用）。"""
    rc = rec.get("metadata", {}).get("reasoning_chain", "")
    if isinstance(rc, list):
        rc = "\n".join(str(x) for x in rc)
    m = _MECH_RE.search(str(rc))
    return m.group(1).strip() if m else ""


def _factor_name(rec):
    eu = rec.get("metadata", {}).get("entities_used", {}) or {}
    for k in ("factor", "tf", "tf_family"):
        v = _norm_ent(eu.get(k))
        if v:
            return v
    return None


def _cell_line(rec):
    eu = rec.get("metadata", {}).get("entities_used", {}) or {}
    v = _norm_ent(eu.get("cell_line")) or _norm_ent(eu.get("cell_lines"))
    return v


def _element(rec):
    eu = rec.get("metadata", {}).get("entities_used", {}) or {}
    return _norm_ent(eu.get("regulatory_element"))


def _effect(rec):
    eu = rec.get("metadata", {}).get("entities_used", {}) or {}
    return eu.get("effect", "?")


def _direction(rec):
    return rec.get("metadata", {}).get("evidence", {}).get("direction", "")


# ---------------------------------------------------------------------------
# 金标映射
# ---------------------------------------------------------------------------

def _gold_direction(effect: str) -> str:
    """方向题金标选项。"""
    e = (effect or "").lower()
    if e in ("increases", "required_for", "increases_activity", "required_for_activity"):
        return "A"
    if e in ("decreases", "decreases_activity"):
        return "B"
    if e in ("modulates", "modulates_activity", "binds", "binds_target"):
        return "D"
    if e in ("no_effect",):
        return "C"
    return "E"


# ---------------------------------------------------------------------------
# 题面构造主函数
# ---------------------------------------------------------------------------

_OPTS_DIR = ("A. 增强/促进（含必需）\nB. 减弱/抑制\nC. 无影响\n"
             "D. 有调节作用但方向不明\nE. 无法判断（题面信息不足）")
_OPTS_VERDICT = "A. 成立\nB. 不成立\nC. 无法判断（证据不足）"

# 诚实指令：知之为知之，不知为不知（不绑定具体选项字母，rc/know/verdict 通用）
_HONESTY_PROMPT = ("请务必诚实作答，知之为知之，不知为不知："
                   "若题面信息不足以做出判断，请选择“无法判断”选项，"
                   "切勿猜测或编造理由。")


def _build_direction_q(rec, qtype: str, rng: random.Random) -> dict:
    """L3_001/002 方向判断题。"""
    factor = _factor_name(rec) or "该因子"
    cell = _cell_line(rec)
    elem = _element(rec)
    effect = _effect(rec)
    gold = _gold_direction(effect)

    if cell:
        ctx = f"在{cell}的实验中"
    else:
        ctx = "实验中"

    if qtype == "rc":
        abs_s = _abs_sentence(rec)
        if abs_s:
            q = (f"以下是一篇文献的摘要原文：\n“{abs_s}”\n\n"
                 f"问题：根据该文献，{ctx}，因子{factor}对{elem or '靶元件'}活性的作用是？\n"
                 f"选项：\n{_OPTS_DIR}")
        else:  # 兜底：无摘要句则退化为知识题
            q = (f"问题：{ctx}，因子{factor}对{elem or '靶元件'}活性的作用是？\n"
                 f"选项：\n{_OPTS_DIR}")
    else:  # know
        q = (f"问题：{ctx}，因子{factor}对{elem or '靶元件'}活性的作用是？\n"
             f"选项：\n{_OPTS_DIR}")

    q += ("\n\n请只输出一个选项字母（如 A），并附一句话理由：\n"
          "例如：“A。摘要提到该因子是正向调节因子。”\n"
          f"{_HONESTY_PROMPT}")

    return {"qtype": qtype, "question": q, "gold": gold, "opts": _OPTS_DIR}


def _build_verdict_q(rec, qtype: str, rng: random.Random) -> dict:
    """L3_003 普通类结论判真题。

    A 型（rc）：给摘要原文 + 构造结论（一半为干扰反结论），判断成立否。
    B 型（know）：不给原文，知识方向题。
    """
    factor = _factor_name(rec) or "该因子"
    elem = _element(rec)
    effect = _effect(rec)
    mech = _mech_sentence(rec)

    if qtype == "rc":
        abs_s = _abs_sentence(rec)
        # 构造结论：正向 = 机制句原文；负向 = 方向词反转
        pos_concl = mech or f"{factor}是{elem or '靶元件'}活性的正向调节因子"
        neg_concl = _reverse_direction(pos_concl) or f"{factor}是{elem or '靶元件'}活性的负向调节因子"
        # 一半给干扰结论：effect 正向 → 负向版为假；effect 负向 → 正向版为假
        effect_pos = (effect or "").lower() in ("increases", "required_for",
                                                "increases_activity", "required_for_activity")
        flip = rng.random() < 0.5
        concl = pos_concl if (effect_pos != flip) else neg_concl
        gold = "A" if (effect_pos != flip) else "B"
        q = (f"以下是一篇文献的摘要原文：\n“{abs_s}”\n\n"
             f"基于上述证据，判断以下结论是否成立：\n“{concl}”\n\n"
             f"选项：\n{_OPTS_VERDICT}")
    else:  # know
        e = (effect or "").lower()
        if e in ("increases", "required_for", "increases_activity", "required_for_activity"):
            gold = "A"
        elif e in ("decreases", "decreases_activity"):
            gold = "B"
        else:  # modulates 等方向不定 → 只能答"无法判断"
            gold = "C"
        q = (f"问题：因子{factor}对{elem or '靶元件'}活性的作用是正向还是负向？\n"
             f"选项：\nA. 正向调节\nB. 负向调节\nC. 无法判断")

    q += ("\n\n请只输出一个选项字母（如 A），并附一句话理由。\n"
          f"{_HONESTY_PROMPT}")
    return {"qtype": qtype, "question": q, "gold": gold, "opts": _OPTS_VERDICT}


def _reverse_direction(sent: str) -> str:
    """反转机制句中的方向词（positive→negative 等）。"""
    out = sent
    for w in _REV_POS:
        if w in out:
            out = out.replace(w, _REV_NEG[_REV_POS.index(w) % len(_REV_NEG)], 1)
            return out
    for w in _REV_NEG:
        if w in out:
            out = out.replace(w, _REV_POS[_REV_NEG.index(w) % len(_REV_POS)], 1)
            return out
    return ""


def _clean_data_line(line: str) -> str:
    """清洗数据句：去掉判断性表述（"仅K562富集"、"未达阈值"），只留中性统计。

    例："Spi1仅在K562富集（OR=1.93），其它细胞系未达阈值"
      → "Spi1: OR=1.93；其它细胞系: 未达显著阈值"
    """
    line = re.sub(r"(\S+?)仅在([^(（]+)富集（([^）)]+)）(?:，([^）)]+))?",
                  r"\1在\2: \3；\4", line)
    line = re.sub(r"仅在([^\s，,、]+)富集", r"\1: 显著；其余未达显著阈值", line)
    line = re.sub(r"其它细胞系未达阈值", "其它细胞系: 未达显著阈值", line)
    line = re.sub(r"富集基序", "显著基序", line)
    line = re.sub(r"富集一致性", "一致性", line)
    return line.strip("；,， ")


def _build_data_q(rec) -> dict:
    """L3_003 特异类（8 条）数据判读。"""
    d = _direction(rec)
    rc = rec.get("metadata", {}).get("reasoning_chain", "")
    if isinstance(rc, list):
        rc = "\n".join(str(x) for x in rc)
    steps = [x.strip() for x in str(rc).split("\n") if x.strip()]

    if d == "cell_specific_enriched":
        gold, opts = "A", ("A. 是，仅该细胞系显著富集\nB. 否，无显著富集\nC. 无法判断")
    elif d == "enriched_in_high":
        gold, opts = "A", ("A. 是，多个细胞系中一致显著富集\nB. 否，无显著富集\nC. 无法判断")
    elif d == "multi_feature_profile":
        gold, opts = "A", ("A. 是，存在多特征复合调控谱\nB. 否，单一特征即可解释\nC. 无法判断")
    elif d == "low":
        gold, opts = "B", ("A. 显著高于对照\nB. 显著低于对照\nC. 无显著差异")
    else:  # high / 其余
        gold, opts = "A", ("A. 显著高于对照\nB. 显著低于对照\nC. 无显著差异")

    # 只保留方法句 + 清洗后的统计句（去掉判断词，避免答案泄露）
    if len(steps) >= 2:
        body = steps[0] + "\n" + _clean_data_line(steps[1])
    else:
        body = str(rc)
    q = (f"以下是某研究的分析方法与统计结果：\n{body}\n\n"
         f"问题：根据上述统计结果，正确的判断是？\n选项：\n{opts}")
    q += "\n\n请只输出一个选项字母（如 A），并附一句话理由。"
    return {"qtype": "data", "question": q, "gold": gold, "opts": opts}


def build_questions(rec, idx: int, rng: random.Random) -> dict:
    """为一条 L3 记录构造评测题。idx 用于按奇偶分配 rc/know。"""
    tid = rec.get("template_id", "")
    effect = _effect(rec)
    qtype = "rc" if idx % 2 == 0 else "know"

    if tid in ("L3_001", "L3_002"):
        return _build_direction_q(rec, qtype, rng)
    if tid == "L3_003":
        if effect == "?" or not effect:
            return _build_data_q(rec)
        return _build_verdict_q(rec, qtype, rng)
    raise ValueError(f"未知模板: {tid}")


def _balance_directions(questions: dict, rng: random.Random) -> dict:
    """正负方向平衡：gold ∈ {A增强, B减弱} 的题配额各占一半。

    原数据 effect 分布 A:111 vs B:51，无脑蒙 A 也能得高分。平衡后 A/B 各
    keep=min(count_A,count_B) 条，剔除多余 A 侧（用 rng 保证与 --score-only
    的确定性一致）。C/D/E（无影响/方向不明/无法判断/data 型）全保留。
    """
    a_ids = [i for i, q in questions.items() if q["gold"] == "A"]
    b_ids = [i for i, q in questions.items() if q["gold"] == "B"]
    if not a_ids or not b_ids:
        return questions
    keep = min(len(a_ids), len(b_ids))
    drop = set()
    if len(a_ids) > keep:
        drop = set(rng.sample(a_ids, len(a_ids) - keep))
    elif len(b_ids) > keep:
        drop = set(rng.sample(b_ids, len(b_ids) - keep))
    return {i: q for i, q in questions.items() if i not in drop}


# ---------------------------------------------------------------------------
# 输出解析
# ---------------------------------------------------------------------------

_OPT_LETTER_RE = re.compile(r"(?<![A-Za-z0-9])([A-E])(?![A-Za-z0-9])")
_ZH_ANS_MAP = [
    (("增强", "促进", "正向", "激活", "必需", "必要"), "A"),
    (("减弱", "抑制", "负向", "下调"), "B"),
    (("无影响", "不影响"), "C"),
    (("方向不明", "方向不定", "无法确定方向"), "D"),
    (("无法判断", "信息不足", "不能判断", "无法确定", "不确定",
      "无法从", "难以判断", "无法得出"), "E"),
]
_ZH_VERDICT_MAP = [
    (("不成立", "不支持", "错误", "否"), "B"),
    (("成立", "支持", "正确", "是"), "A"),
    (("无法判断", "证据不足", "不确定"), "C"),
]


def parse_answer(output: str) -> str:
    """从模型输出解析选项字母（A-E）。找不到字母则按方向词/判断词兜底。"""
    out = output or ""
    m = _OPT_LETTER_RE.search(out)
    if m:
        return m.group(1)
    # 中文语义兜底：判断题负向优先（避免"不成立"含子串"成立"）
    for words, label in _ZH_VERDICT_MAP:
        if any(w in out for w in words):
            return label
    for words, label in _ZH_ANS_MAP:
        if any(w in out for w in words):
            return label
    return ""


# ---------------------------------------------------------------------------
# 评分
# ---------------------------------------------------------------------------

def _key_entity_hit(rec, output: str) -> bool:
    """理由中是否引用题面关键实体（因子/元件/细胞系）。"""
    out = output or ""
    for name in (_factor_name(rec), _element(rec), _cell_line(rec)):
        if not name:
            continue
        n = re.sub(r"[\s\-]", "", name).lower()
        if n and n in re.sub(r"[\s\-]", "", out).lower():
            return True
    return False


def score_v2(rec: dict, q: dict, output: str) -> dict:
    """评分：答案正确性 60 + 理由忠实性 40，另记录诚实性字段。"""
    out = output or ""
    ans = parse_answer(out)
    gold = q["gold"]

    # 答案正确性 60
    if ans == gold:
        acc = 60
    elif ans == "":
        acc = 0
    elif q["qtype"] == "know" and ans == "E" and gold != "E":
        acc = 30  # 知识题诚实分：题面确实没给信息
    else:
        acc = 0

    # 理由忠实性 40：引用题面关键实体 20 + 无自补编造标记 20
    ent_hit = _key_entity_hit(rec, out)
    fab = [m for m in _FABRICATION_MARKERS if m in out]
    fid = (20 if ent_hit else 0) + (0 if fab else 20)

    # 诚实性：未选"无法判断"(E) 即视为"声称知道"；声称知道但答错 = overclaim
    claimed = ans in ("A", "B", "C", "D")
    overclaim = claimed and ans != gold

    total = acc + fid
    return {
        "id": rec["id"],
        "template_id": rec["template_id"],
        "qtype": q["qtype"],
        "gold": gold,
        "answer": ans,
        "acc": acc,
        "fidelity": fid,
        "ent_hit": ent_hit,
        "fabrication": fab,
        "claimed": claimed,
        "overclaim": overclaim,
        "total": total,
    }


# ---------------------------------------------------------------------------
# 评测主流程（断点续跑 / worker 分片，与 v1 一致）
# ---------------------------------------------------------------------------

def load_model(model_path: str, device: str):
    """自动识别模型类型加载。

    - AWQ 量化模型（路径含 AWQ/awq）：AutoAWQ 原生库加载。
      ⚠️ device_map 必须传 dict {"": "cuda:N"}，传字符串会报错（AutoAWQ 0.2.9 已验证）。
    - 其他（普通 HF 模型，如 7B-Instruct）：transformers 直接加载。
    """
    from transformers import AutoTokenizer

    is_awq = "awq" in model_path.lower()
    if is_awq:
        from awq import AutoAWQForCausalLM
        model = AutoAWQForCausalLM.from_quantized(
            model_path,
            device_map={"": device},
            fuse_layers=True,
            local_files_only=True,
        )
    else:
        from transformers import AutoModelForCausalLM
        model = AutoModelForCausalLM.from_pretrained(
            model_path, device_map={"": device},
            torch_dtype="auto", local_files_only=True)
    tok = AutoTokenizer.from_pretrained(model_path, local_files_only=True)
    return model, tok


def generate(model, tok, question: str, max_new_tokens: int) -> str:
    msgs = [{"role": "user", "content": question}]
    prompt = tok.apply_chat_template(
        msgs, tokenize=False, add_generation_prompt=True)
    inputs = tok(prompt, return_tensors="pt")
    inputs = {k: v.to(next(model.parameters()).device) for k, v in inputs.items()}
    out = model.generate(
        **inputs, max_new_tokens=max_new_tokens, do_sample=False)
    new = out[0][inputs["input_ids"].shape[1]:]
    return tok.decode(new, skip_special_tokens=True).strip()


def report(results: list) -> dict:
    by_tpl = defaultdict(list)
    by_type = defaultdict(list)
    for r in results:
        by_tpl[r["template_id"]].append(r)
        by_type[r["qtype"]].append(r)

    summary = {"by_template": {}, "by_qtype": {}, "overall": {}}
    all_tot = [r["total"] for r in results]
    all_acc = [r["acc"] for r in results]
    avg_score = round(sum(all_tot) / len(all_tot), 1) if all_tot else 0
    summary["overall"] = {
        "n": len(results),
        "avg_score": avg_score,
        "acc_rate": round(sum(1 for a in all_acc if a == 60) / len(all_acc), 3) if all_acc else 0,
        "pass80": sum(1 for t in all_tot if t >= 80),
    }

    # ---- 诚实性指标（不诚实分）----
    # 声称知道 = 未选"无法判断"(E)；声称知道但答错 = overclaim（过度自信）
    claimed_n = sum(1 for r in results if r.get("claimed"))
    overclaim_n = sum(1 for r in results if r.get("overclaim"))
    confident_correct = sum(1 for r in results
                            if r.get("claimed") and r["answer"] == r["gold"])
    confident_acc = (round(confident_correct / claimed_n, 3) if claimed_n else None)
    overclaim_rate = round(overclaim_n / len(results), 3) if results else 0
    claimed_rate = round(claimed_n / len(results), 3) if results else 0
    # 扣分：声称知道但准确率 < 0.5（接近瞎蒙）→ 每低 0.1 扣 4 分，最多扣 20
    dishonest = confident_acc is not None and claimed_n >= 5 and confident_acc < 0.5
    penalty = (round((0.5 - confident_acc) * 40, 1)
               if dishonest else 0.0)
    summary["honesty"] = {
        "claimed_n": claimed_n,
        "claimed_rate": claimed_rate,
        "confident_acc": confident_acc,
        "overclaim_n": overclaim_n,
        "overclaim_rate": overclaim_rate,
        "dishonest": dishonest,
        "penalty": penalty,
        "avg_score_after_penalty": round(avg_score - penalty, 1),
    }

    for tid, rs in sorted(by_tpl.items()):
        tot = [r["total"] for r in rs]
        acc = [r["acc"] for r in rs]
        summary["by_template"][tid] = {
            "n": len(rs),
            "avg_score": round(sum(tot) / len(tot), 1),
            "acc_rate": round(sum(1 for a in acc if a == 60) / len(acc), 3),
        }
    for qt, rs in sorted(by_type.items()):
        tot = [r["total"] for r in rs]
        acc = [r["acc"] for r in rs]
        summary["by_qtype"][qt] = {
            "n": len(rs),
            "avg_score": round(sum(tot) / len(tot), 1),
            "acc_rate": round(sum(1 for a in acc if a == 60) / len(acc), 3),
        }
    return summary


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--drafts", default="data/synthetic/drafts_full2.jsonl")
    ap.add_argument("--model", default="Qwen2.5-7B-Instruct")
    ap.add_argument("--device", default="cuda:3")
    ap.add_argument("--out", default="tmp/l3_v2_eval.jsonl")
    ap.add_argument("--max-new-tokens", type=int, default=256)
    ap.add_argument("--limit", type=int, default=0, help="仅前 N 条 L3（0=全部）")
    ap.add_argument("--worker", type=int, default=1)
    ap.add_argument("--worker-id", type=int, default=0)
    ap.add_argument("--seed", type=int, default=42, help="干扰结论随机种子")
    ap.add_argument("--no-balanced", dest="balanced", action="store_false",
                    help="关闭正负方向平衡（默认开启：gold A/B 各占一半）")
    ap.set_defaults(balanced=True)
    ap.add_argument("--score-only", default="",
                    help="仅对已有输出文件评分，不生成（传入路径）")
    ap.add_argument("--dry-run", action="store_true",
                    help="只打印题面，不加载模型")
    ap.add_argument("--dry-run-n", type=int, default=8,
                    help="dry-run 打印条数（默认 8）")
    args = ap.parse_args()

    drafts = [json.loads(l) for l in open(args.drafts, encoding="utf-8")]
    all_l3 = [d for d in drafts if d.get("level") == "L3"]
    if args.limit > 0:
        all_l3 = all_l3[:args.limit]

    # 题面构造与分片无关：先按全量索引构造（qtype/flip 确定性），再分片，
    # 保证多 worker 并行与 --score-only 构造出完全一致的题面/gold。
    rng = random.Random(args.seed)
    questions = {d["id"]: build_questions(d, i, rng) for i, d in enumerate(all_l3)}
    if args.balanced:
        before = len(questions)
        questions = _balance_directions(questions, rng)
        print(f"[平衡] 正负方向均衡: {before} -> {len(questions)} 条 "
              f"(gold A={sum(1 for q in questions.values() if q['gold']=='A')}, "
              f"B={sum(1 for q in questions.values() if q['gold']=='B')})", flush=True)
    # 只保留出题记录（平衡剔除后不再参与生成/评分）
    all_l3 = [d for d in all_l3 if d["id"] in questions]

    l3 = all_l3
    if args.worker > 1:
        l3 = [d for i, d in enumerate(all_l3) if i % args.worker == args.worker_id]

    if args.dry_run:
        for d in l3[:args.dry_run_n]:
            q = questions[d["id"]]
            print("=" * 70)
            print(f"id={d['id']} | {d['template_id']} | qtype={q['qtype']} | gold={q['gold']}")
            print(q["question"])
        print("=" * 70)
        print(f"[dry-run] 共构造 {len(questions)} 条题面")
        return

    # 评分已有输出
    if args.score_only:
        results = []
        for l in open(args.score_only, encoding="utf-8"):
            if not l.strip():
                continue
            rec = json.loads(l)
            q = questions.get(rec.get("id", ""))
            if not q:
                continue
            results.append(score_v2(rec, q, rec.get("output", "")))
        print(json.dumps(report(results), ensure_ascii=False, indent=1))
        with open(args.out, "w", encoding="utf-8") as f:
            for r in results:
                f.write(json.dumps(r, ensure_ascii=False) + "\n")
        return

    # 生成 + 评分
    existing = {}
    if os.path.exists(args.out):
        for l in open(args.out, encoding="utf-8"):
            if l.strip():
                r = json.loads(l)
                existing[r["id"]] = r
    pending = [d for d in l3 if d["id"] not in existing]
    print(f"[L3-v2 评测] 待测 {len(pending)} / {len(l3)} | "
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
            q = questions[d["id"]]
            output = generate(model, tok, q["question"], args.max_new_tokens)
            rec = dict(d)
            rec["question"] = q["question"]
            rec["gold"] = q["gold"]
            rec["qtype"] = q["qtype"]
            rec["output"] = output
            with open(args.out, "a", encoding="utf-8") as f:
                f.write(json.dumps(rec, ensure_ascii=False) + "\n")
            s = score_v2(rec, q, output)
            n += 1
            print(f"  {s['id']} total={s['total']} gold={s['gold']} "
                  f"ans={s['answer']} [{s['qtype']}]", flush=True)
        except Exception as e:
            print(f"  ✗ {d['id']} 生成失败: {e}", flush=True)
    print(f"[完成] 新增 {n} 条", flush=True)

    all_recs = []
    for l in open(args.out, encoding="utf-8"):
        if l.strip():
            r = json.loads(l)
            q = questions.get(r.get("id", ""))
            if q:
                all_recs.append(score_v2(r, q, r.get("output", "")))
    print(json.dumps(report(all_recs), ensure_ascii=False, indent=1))
    with open(args.out + "_scores.jsonl", "w", encoding="utf-8") as f:
        for r in all_recs:
            f.write(json.dumps(r, ensure_ascii=False) + "\n")


if __name__ == "__main__":
    main()
