#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
template_engine.py — Module 2 第一步：TemplateFiller 填充引擎

功能：
1. 解析 templates/*.yaml（种子 14 个模板）
2. 从 knowledge_base_clean.jsonl 加载 617 条断言
3. 占位符 {field} → 断言字段归一化填充（支持通用/文献类/统计类三类映射）
4. 断言类型筛选：require_claim_type / require_effect / 占位符→记录类型推断
5. 同源约束：模板内多个实体占位符必须来自同一条断言
6. 轮询采样：每模板 3-5 个实体组合，优先覆盖不同断言（均匀性）
7. 输出草稿 JSONL：含 template_id/level/instruction/raw_output_placeholder/metadata

用法：
    python src/data_synthesis/template_engine.py \
        --kb data/processed/knowledge_base_clean.jsonl \
        --templates templates/ \
        --out data/synthetic/drafts.jsonl \
        --per-template 4
"""

import argparse
import glob
import json
import os
import random
import re
import sys
from collections import Counter, defaultdict

# 保证从 src/data_synthesis 目录运行时能 import（无第三方依赖，仅标准库）

PLACEHOLDER_RE = re.compile(r"\{(\w+)\}")

# ---------------------------------------------------------------------------
# 占位符 → 字段归一化映射
# ---------------------------------------------------------------------------

# 通用占位符：任何记录类型均可解析
GENERIC_FIELDS = {
    "claim_text": "claim",
    "claim_type": "claim_type",
    "reasoning_chain": "reasoning_chain",
    "evidence_direction": ("evidence", "direction"),
    "p_value": ("evidence", "p_value"),
    "effect_size": ("evidence", "effect_size"),
    "confidence": "confidence",
}

# 文献类占位符（claim_type = mechanistic / design_rule）
LITERATURE_FIELDS = {
    "factor": ("entities", "factor"),
    "factor_type": ("entities", "factor_type"),
    "cell_line": ("entities", "cell_line"),
    "regulatory_element": ("entities", "regulatory_element"),
    "effect": ("entities", "effect"),
    "pmids": ("entities", "pmids"),
    "context": ("entities", "context"),
}

# 统计类占位符（claim_type = association）
STAT_FIELDS = {
    "tf": ("entities", "tf"),
    "motif": ("entities", "motif"),
    "gc_range": ("entities", "gc_range"),
    "gc_lo": ("entities", "gc_lo"),
    "gc_hi": ("entities", "gc_hi"),
    "gc_mean": ("entities", "gc_mean"),
    "gc_std": ("entities", "gc_std"),
    "shape": ("entities", "shape"),
    "feature": ("entities", "feature"),
    "direction": ("evidence", "direction"),
}

# GC-motif 关联层占位符（claim_type = gc_association，来自 knowledge_base_gc_motif.jsonl）
# 独立类型：与活性类 association 互斥，避免"高GC序列"模板与"高活性序列"模板串扰
GC_ASSOC_FIELDS = {
    "gc_class": ("entities", "gc_class"),
    "motif_gc": ("entities", "motif_gc"),
    "gc_rich_q": ("metadata", "gc_rich_q"),
    "gc_poor_q": ("metadata", "gc_poor_q"),
}

# 归一化映射：所有占位符 → (路径或字段名, 记录类型约束)
PLACEHOLDER_MAP = {}
for _ph, _spec in GENERIC_FIELDS.items():
    PLACEHOLDER_MAP[_ph] = (_spec, None)
for _ph, _spec in LITERATURE_FIELDS.items():
    PLACEHOLDER_MAP[_ph] = (_spec, ("mechanistic", "design_rule"))
for _ph, _spec in STAT_FIELDS.items():
    PLACEHOLDER_MAP[_ph] = (_spec, ("association",))
for _ph, _spec in GC_ASSOC_FIELDS.items():
    PLACEHOLDER_MAP[_ph] = (_spec, ("gc_association",))

# 单元格 line 需要具体细胞系：过滤 not_specified
CONCRETE_CELL_LINE_TEMPLATE_MARKERS = {
    "在 {cell_line} 细胞背景下", "在 {cell_line} 中", "的 {cell_line} 高活性",
}

# 特殊占位符需要格式化处理
LIST_PLACEHOLDERS = {"reasoning_chain", "pmids"}

# ---------------------------------------------------------------------------
# regulatory_element 白名单（用户决策 2026-08-06）
# - 过滤非元件取值（过程/表型/表达词等），只保留真正的调控元件
# - 带基因名的元件（如 foxp4 promoter / il-6 promoter）保留细粒度，
#   并提取基因名 → 下游注入"基因特异性"知识
# ---------------------------------------------------------------------------

# 合法元件词根：命中任一即视为调控元件
VALID_ELEMENT_KEYWORDS = (
    "promoter", "enhancer", "insulator", "element", "region",
    "intronic", "gene body", "motif", "binding site", "domain",
    "locus control", "utr", "exon", "gene",
)

# 排除词：过程/表型/表达/通路等，出现任一即视为非元件
INVALID_ELEMENT_KEYWORDS = (
    "expression", "migration", "invasion", "proliferation",
    "outgrowth", "differentiation", "pathway", "transcription",
    "upregulation", "activity", "function", "methylation level",
)

# 基因名 + 元件词模式：如 "foxp4 promoter" / "il-6 promoter" / "β-globin gene"
GENE_ELEMENT_RE = re.compile(
    r"^([a-z0-9][a-z0-9\-_\.]*)\s+(promoter|enhancer|gene)$", re.IGNORECASE
)

# 基因名判定：内置常用基因名表 + 知识库 factor 自洽（用户决策 2026-08-06）
from gene_names import is_gene_name, register_kb_genes  # noqa: E402


def parse_regulatory_element(value) -> tuple:
    """校验 regulatory_element 合法性并提取基因名。
    返回 (is_valid: bool, gene: str|None)。gene 非 None 表示知识基因特异。
    判定：前缀命中内置基因名表（归一化后）→ 基因特异。"""
    if not value:
        return False, None
    v = str(value).strip().lower()
    # 排除非元件词
    if any(k in v for k in INVALID_ELEMENT_KEYWORDS):
        return False, None
    # 必须命中合法元件词根
    if not any(k in v for k in VALID_ELEMENT_KEYWORDS):
        return False, None
    # 提取基因名（细粒度保留）："foxp4 promoter" → gene="foxp4"
    m = GENE_ELEMENT_RE.match(v)
    token = m.group(1) if m else None
    if token and is_gene_name(token):
        return True, token  # 前缀是基因名 → 基因特异
    return True, None  # 合法元件，但非基因特异


def deep_get(record: dict, path) -> object:
    """按路径取字段：str 为直接键，tuple 为嵌套键序列。取不到返回 None。"""
    if isinstance(path, str):
        return record.get(path)
    cur = record
    for key in path:
        if not isinstance(cur, dict):
            return None
        cur = cur.get(key)
    return cur


def fmt_value(ph: str, value) -> str:
    """占位符值 → 文本。列表转编号列表/顿号连接，其余 str()。"""
    if value is None:
        return None
    if ph == "reasoning_chain" and isinstance(value, list):
        return "\n".join(f"{i}. {s}" for i, s in enumerate(value, 1))
    if ph == "pmids" and isinstance(value, list):
        return "、".join(f"PMID {p}" for p in value)
    if ph == "gc_class":
        # GC 类别中文化：gc_rich -> GC-rich 基序 / at_rich -> AT-rich 基序 / neutral -> 中性GC偏好基序
        return {"gc_rich": "GC-rich基序", "at_rich": "AT-rich基序",
                "neutral": "中性GC偏好基序"}.get(str(value), str(value))
    if isinstance(value, float):
        # 科学计数法简短化，避免 1.27e-257 被截断
        return f"{value:.2e}" if abs(value) < 1e-4 or abs(value) > 1e4 else str(round(value, 4))
    return str(value)


def extract_placeholders(template_text: str) -> list:
    """提取模板文本中的占位符列表（保序、去重）。"""
    seen = set()
    out = []
    for m in PLACEHOLDER_RE.finditer(template_text):
        ph = m.group(1)
        if ph not in seen:
            seen.add(ph)
            out.append(ph)
    return out


def infer_record_types(placeholders: list) -> set:
    """由占位符推断可用的记录类型集合（None 表示不限制）。"""
    types = set()
    unconstrained = False
    for ph in placeholders:
        spec, constraint = PLACEHOLDER_MAP.get(ph, (None, None))
        if spec is None:
            continue
        if constraint is None:
            unconstrained = True
        else:
            types.update(constraint)
    # 若既有通用占位符又有类型限制占位符 → 取交集约束
    if unconstrained and types:
        return types
    if types:
        return types
    return None  # 全通用 → 不限制


def record_matches_placeholders(record: dict, placeholders: list) -> bool:
    """检查记录能否解析模板所有占位符（同源约束 + 字段存在性）。"""
    for ph in placeholders:
        spec, _ = PLACEHOLDER_MAP.get(ph, (None, None))
        if spec is None:
            return False  # 未知占位符 → 无法填充
        val = deep_get(record, spec)
        if val is None:
            return False
        if isinstance(val, list) and len(val) == 0:
            return False
    return True


def fill_template(record: dict, instruction: str) -> str:
    """用记录字段填充模板文本，返回填充后的指令。"""
    out = instruction
    for ph in extract_placeholders(instruction):
        spec, _ = PLACEHOLDER_MAP.get(ph, (None, None))
        val = deep_get(record, spec)
        if val is None:
            raise ValueError(f"占位符 {{{ph}}} 在记录 {record.get('claim_type')} 中无法解析")
        out = out.replace("{" + ph + "}", fmt_value(ph, val))
    return out


# ---------------------------------------------------------------------------
# 断言筛选
# ---------------------------------------------------------------------------

def filter_records(records: list, template: dict) -> list:
    """按模板约束筛选可用断言：
    1. require_claim_type
    2. require_effect（entities.effect 匹配）
    3. 占位符→记录类型推断
    4. 全占位符可解析（同源约束）
    5. 具体细胞系要求（叙述型上下文中的 {cell_line} → 排除 not_specified；
       问答型如"哪个细胞系"豁免，not_specified 是合法答案）
    """
    plist = extract_placeholders(template["instruction_template"])
    req_types = template.get("require_claim_type")
    req_effects = template.get("require_effect")

    # 具体细胞系判定：
    # - 叙述型信号模式（出现则要求具体细胞系）
    # - 问答型豁免模式（如"哪个细胞系"→ not_specified 可作答案）
    text = template["instruction_template"]
    need_concrete_cell = False
    if "{cell_line}" in text:
        concrete_signals = (
            "in {cell_line} cells", "in {cell_line}",
            "在 {cell_line} 中", "在 {cell_line} 细胞背景下",
            "的 {cell_line} 高活性", "{cell_line} 细胞背景",
            "用 {cell_line} 的实验证据", "{cell_line} 中",
        )
        exempt_signals = ("哪个细胞系", "哪一种细胞系")
        need_concrete_cell = (
            any(s in text for s in concrete_signals)
            and not any(s in text for s in exempt_signals)
        )

    # regulatory_element 白名单：模板含 {regulatory_element} 时启用
    need_valid_element = "{regulatory_element}" in text

    candidates = []
    for r in records:
        ct = r.get("claim_type")
        # 1. require_claim_type（支持 str 或 list）
        if req_types is not None:
            wanted = [req_types] if isinstance(req_types, str) else req_types
            if ct not in wanted:
                continue
        # 2. require_effect
        if req_effects is not None:
            eff = r.get("entities", {}).get("effect")
            wanted = [req_effects] if isinstance(req_effects, str) else req_effects
            if eff not in wanted:
                continue
        # 3. 占位符推断类型
        inferred = infer_record_types(plist)
        if inferred is not None and ct not in inferred:
            continue
        # 4. 同源：所有占位符可解析
        if not record_matches_placeholders(r, plist):
            continue
        # 5. 具体细胞系
        if need_concrete_cell and r.get("entities", {}).get("cell_line") == "not_specified":
            continue
        # 6. regulatory_element 白名单（过滤过程/表型/表达词等非元件）
        if need_valid_element:
            el = r.get("entities", {}).get("regulatory_element")
            ok, _gene = parse_regulatory_element(el)
            if not ok:
                continue
        candidates.append(r)
    return candidates


# ---------------------------------------------------------------------------
# 轮询采样
# ---------------------------------------------------------------------------

def round_robin_sample(candidates: list, k: int, rng: random.Random,
                       key_func=None) -> list:
    """轮询采样 k 个不同断言；不足 k 个时循环复用（保证 k 组合）。
    key_func 用于多样性（如按 claim 去重）。"""
    if not candidates:
        return []
    key_func = key_func or (lambda r: r.get("claim"))
    # 按 key 分组
    groups = defaultdict(list)
    for c in candidates:
        groups[key_func(c)].append(c)
    keys = list(groups.keys())
    rng.shuffle(keys)
    picked = []
    idx = 0
    while len(picked) < k and keys:
        key = keys[idx % len(keys)]
        picked.append(groups[key][0])
        idx += 1
        if idx % len(keys) == 0:
            # 一轮过后打乱顺序再进入下一轮（提高均匀性）
            rng.shuffle(keys)
    return picked


# ---------------------------------------------------------------------------
# 主流程
# ---------------------------------------------------------------------------

def load_templates(templates_dir: str) -> list:
    """加载 templates/*.yaml（跳过 README），校验 id 唯一。"""
    import yaml
    templates = []
    ids = set()
    for f in sorted(glob.glob(os.path.join(templates_dir, "L*.yaml"))):
        with open(f, encoding="utf-8") as fh:
            data = yaml.safe_load(fh)
        if not data:
            continue
        for t in data:
            t["_source_file"] = os.path.basename(f)
            tid = t.get("id")
            if tid in ids:
                raise ValueError(f"重复模板 id: {tid}")
            ids.add(tid)
            templates.append(t)
    return templates


def load_kb(path: str) -> list:
    records = []
    factors = []
    with open(path, encoding="utf-8") as f:
        for line in f:
            line = line.strip()
            if line:
                rec = json.loads(line)
                records.append(rec)
                f_ = rec.get("entities", {}).get("factor")
                if f_:
                    factors.append(f_)
    # 知识库 factor 并入基因名表（内部自洽：本项目语境下的基因）
    register_kb_genes(factors)
    return records


def build_draft(template: dict, record: dict, instruction: str, idx: int) -> dict:
    """构造草稿样本（含溯源元数据）。"""
    # 基因特异性：regulatory_element 带基因名时（如 foxp4 promoter），
    # 该知识仅对该基因的元件有效 → 供 llm_enhancer 注入特异性知识
    gene = None
    el = record.get("entities", {}).get("regulatory_element")
    if el:
        _ok, gene = parse_regulatory_element(el)
    return {
        "id": f"{template['id']}_{idx}",
        "level": template.get("level"),
        "template_id": template["id"],
        "instruction": instruction,
        "input": "",
        "raw_output_placeholder": "",  # L1/L2 可直接由模板答案规则生成；L3+ 待 LLM 增强
        "metadata": {
            "source_claim_ids": [record.get("claim_id", record.get("_idx"))],
            "claim_type": record.get("claim_type"),
            "entities_used": record.get("entities", {}),
            "reasoning_chain": record.get("reasoning_chain"),
            "evidence": record.get("evidence", {}),
            "gene_specific": gene,  # str|None：知识是否基因特异
            "llm_model": "",  # 由 llm_enhancer 填写
            "quality_score": None,
        },
    }


def main():
    ap = argparse.ArgumentParser(description="TemplateFiller 填充引擎")
    ap.add_argument("--kb", default="data/processed/knowledge_base_clean.jsonl")
    ap.add_argument("--templates", default="templates/")
    ap.add_argument("--out", default="data/synthetic/drafts.jsonl")
    ap.add_argument("--per-template", type=int, default=4,
                    help="每个模板采样组合数（3-5）")
    ap.add_argument("--seed", type=int, default=42)
    args = ap.parse_args()

    rng = random.Random(args.seed)

    templates = load_templates(args.templates)
    records = load_kb(args.kb)
    # 记录行号作为 fallback id
    for i, r in enumerate(records):
        r.setdefault("_idx", i)

    print(f"[templates] {len(templates)} 个 | [kb] {len(records)} 条", flush=True)

    drafts = []
    stats = Counter()
    unmatched = []
    for t in templates:
        candidates = filter_records(records, t)
        picked = round_robin_sample(candidates, args.per_template, rng)
        n_pick = len(picked)
        stats[f"{t['id']}"] = n_pick
        if n_pick == 0:
            unmatched.append(t["id"])
            continue
        for i, rec in enumerate(picked):
            instruction = fill_template(rec, t["instruction_template"])
            draft = build_draft(t, rec, instruction, i + 1)
            drafts.append(draft)

    # 输出
    os.makedirs(os.path.dirname(args.out) or ".", exist_ok=True)
    with open(args.out, "w", encoding="utf-8") as f:
        for d in drafts:
            f.write(json.dumps(d, ensure_ascii=False) + "\n")

    print(f"[产出] {len(drafts)} 条草稿 -> {args.out}", flush=True)
    by_level = Counter(d["level"] for d in drafts)
    print("[层级分布]", dict(by_level), flush=True)
    if unmatched:
        print("[⚠ 无匹配断言模板]", unmatched, flush=True)
    else:
        print("[全部模板均有匹配断言]", flush=True)

    # 抽样展示 2 条
    for d in drafts[:2]:
        print("--- 样例 ---", flush=True)
        print(f"[{d['id']}] level={d['level']}", flush=True)
        print("instruction:", d["instruction"][:200], flush=True)


if __name__ == "__main__":
    main()
