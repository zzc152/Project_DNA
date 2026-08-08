"""知识库受控词表校验 + 细胞系严格过滤 + 打回队列生成器。

决策依据（用户判决 2026-08-08）:
1. **JASPAR motif 受控**：因子若为 motif/TF，必须能映射到 JASPAR2024 词表
   （data/vocab/jaspar2024_vertebrates_core.json，2059 个脊椎动物 CORE）。
   词表外 → 拒绝入库，标记 factor_status=unverified，进打回队列待二次确认。
2. **细胞系严格模式**：只保留 K562 / HepG2 / SK-N-SH 三种细胞系的记录。
   not_specified 打回（待人工裁决）；其他细胞系/缺失 → 标记 out_of_scope。
3. **因子类型约束**：组蛋白修饰（H3K4me3 等）、DNA methylation、染色质
   accessibility 等"标记/状态"禁止作为 factor（路径 A：禁止入 factor 槽）。
   这类记录 factor_status=wrong_type，进打回队列。

输出:
- <out_dir>/kb_curated.jsonl      —— 通过全部校验的记录
- <out_dir>/kb_rejected.jsonl     —— 被拒记录（含结构化拒绝理由）
- <out_dir>/kb_curation_report.json —— 汇总报告（分桶统计 + 明细）

用法:
    python scripts/curate_kb_vocab.py \
        --kb data/processed/knowledge_base_clean_final.jsonl \
        --vocab data/vocab/jaspar2024_vertebrates_core.json \
        --out tmp/curation
"""

import argparse
import json
import re
from collections import Counter
from pathlib import Path

# ---- 用户判决：严格模式只保留三种细胞系 ----
ALLOWED_CELL_LINES = {"K562", "HepG2", "SK-N-SH"}

# ---- 禁止作为 factor 的"标记/状态"模式（路径 A） ----
# 命中这些模式 → factor 类型错误（wrong_type），拒绝入 factor 槽
FORBIDDEN_FACTOR_PATTERNS = [
    re.compile(r"H[34][A-Z]\d", re.I),              # H3K4me3, H4K16ac 等组蛋白修饰
    re.compile(r"\bH[0-9]\s", re.I),                # "H3 acetylation at lysine..."
    re.compile(r"histone\s+\w*", re.I),             # histone acetylation/methylation
    re.compile(r"\bme[123]\b", re.I),               # me1/me2/me3 后缀
    re.compile(r"\bac\b|\bph\b|\bub\b", re.I),      # ac/ph/ub 修饰后缀（严格单独词）
    re.compile(r"DNA\s*methylation", re.I),         # DNA methylation
    re.compile(r"chromatin\s+accessibility", re.I), # chromatin accessibility
    re.compile(r"nucleosome\s+occupancy", re.I),    # nucleosome occupancy
    re.compile(r"epigenetic\s+status", re.I),       # 描述性短语
    re.compile(r"acetylating|de-methylating|deacetyl", re.I),  # 动作短语
]

# ---- 描述性短语/多实体当因子的启发式 ----
# 命中 → factor 名噪声（unverified）
PHRASE_PATTERNS = [
    re.compile(r"\b(tf|tfs|transcription factors?|the machinery|regulators?)\b", re.I),
    re.compile(r"\b(downstream|upstream)\b", re.I),   # "TLR4 downstream NF-κB..."
    re.compile(r"\b(dependent|mediated|induced)\b", re.I),  # "RET-dependent..."
    re.compile(r",\s*(and|&)?\s*", re.I),             # 多实体列表 "GATA1, TAL1"
    re.compile(r"\s+and\s+", re.I),                   # "X and Y"
]
MAX_FACTOR_LEN = 30  # >30 字符视为描述性短语（与 _clean_factor 一致）

# ---- 双 schema 检测：MPRA 记录用 entities.tf + entities.motif(MAxxxx.x) ----
MPRA_MOTIF_RE = re.compile(r"^MA\d+\.\d+$", re.I)

# ---- 染色质修饰酶/表观酶（JASPAR 不收，因为不直接结合 DNA） ----
# 这些是合法 factor（prompt 明确允许"染色质修饰酶"），但 factor_type 常被误标为 TF。
# 命中 → 修正 factor_type=enzyme，豁免 JASPAR 词表强制校验。
ENZYME_HINTS = [
    "HDAC", "DNMT", "KMT", "KDM", "PRMT", "SETD", "SET", "EP300", "CREBBP",
    "WDR5", "ALKBH", "BRD", "SMYD", "KAT", "TET", "UTX", "EZH", "SUV39",
    "EED", "SUZ12", "MLL", "COMPASS", "PCAF", "CBP", "P300", "MOF", "TIP60",
    "SIRT", "LSD", "JMJD", "NSD", "DOT1", "EHMT", "G9A", "GLP", "SIN3",
    "NuRD", "SWI/SNF", "BAF", "PBAF", "INO80", "ISWI", "CHD", "ACF",
]
# 精确名单：短名/别名（避免前缀误伤，如 SET 会误伤 SETX）
ENZYME_EXACT = {"WDR5", "BRD4", "EP300", "CBP", "P300", "G9a", "G9A", "SUV39H1"}

# ---- 氨基酸突变因子模式（S87N / R89W / R85Q —— 不是实体名，是突变描述） ----
MUTATION_PATTERN = re.compile(r"^[A-Z][a-z]?\d+[A-Z][a-z]?$")

# ---- 明确非 TF 实体的噪声名单（人工核实过的打回） ----
NOISE_FACTORS = {
    "VPA",           # 丙戊酸（药物）
    "N1IC",          # Notch 胞内结构域（蛋白片段）
    "NUP153",        # 核孔蛋白
    "5'HS2",         # DNase 高敏位点（序列区域）
    "SEMA3C_L1ME4a", # 重复元件注释
    "OsDDE9",        # 非人类基因
    "PGC1A-hEn1",    # 增强子克隆名
    "PGC1A-En1",
    "MEF2C enhancer mutation",  # 描述性
}

# ---- 基因别名 → JASPAR 词表名（人工核实） ----
# 同一基因在文献用别名，JASPAR 用 HGNC 正式名
GENE_ALIASES = {
    "SREBP1": "SREBF1",  # 脂代谢转录因子，JASPAR 正式名 SREBF1
    "SREBP-1": "SREBF1",
    "FOG1": None,        # ZFPM1 转录辅因子，不直接结合 DNA → 词表无，打回
    "BRG1": "SMARCA4",   # 染色质重塑酶（SWI/SNF ATPase）
    "TAL1": None,        # JASPAR 只有 GATA1::TAL1 / TAL1::TCF3 复合体 → 走复杂匹配
}


def classify_factor(factor: str) -> tuple[str, list[str]]:
    """对 factor 分类。返回 (status, reasons)。

    status ∈ {ok, wrong_type, unverified}:
    - ok:        可作为因子（但 motif/TF 仍需词表命中，见 check_vocab）
    - wrong_type: 组蛋白修饰/甲基化/状态标记 —— 禁止入 factor 槽（路径 A）
    - unverified: 描述性短语/多实体/过长 —— 打回待人工
    """
    reasons = []
    if not factor or not factor.strip():
        return "unverified", ["empty_factor"]

    # 1) 类型错误：标记/状态当因子（硬规则，优先级最高）
    for pat in FORBIDDEN_FACTOR_PATTERNS:
        if pat.search(factor):
            reasons.append(f"wrong_type:matches_pattern:{pat.pattern}")
    if reasons:
        return "wrong_type", reasons

    # 2) 描述性短语/多实体/过长（启发式，打回人工）
    for pat in PHRASE_PATTERNS:
        if pat.search(factor):
            reasons.append(f"unverified:matches_pattern:{pat.pattern}")
    if len(factor) > MAX_FACTOR_LEN:
        reasons.append(f"unverified:too_long:{len(factor)}chars")
    if reasons:
        return "unverified", reasons

    return "ok", []


def check_vocab(factor: str, factor_type: str, vocab: dict, motif_names: set) -> tuple[str, list[str]]:
    """对 ok 类因子做词表校验。返回 (vocab_status, reasons)。

    规则：
    - factor_type in {TF, motif} → 必须在 JASPAR 词表中（matrix_id 或 name 命中）
    - 复合体命名：词表名 'GATA1::TAL1' 拆 '::' 匹配单因子
    - 家族前缀：factor 'CREB' 前缀命中 'CREB1'（len>=4 防误伤）
    - 词表未命中 → unverified（打回，待二次确认）
    - 非 TF/motif（如 enzyme、sequence_feature）→ 不强制词表（标记 type 供人工看）
    """
    if factor_type not in ("TF", "motif"):
        return "ok", [f"non_lexical_type:{factor_type}"]

    # 别名映射（SREBP1 → SREBF1）
    canonical = GENE_ALIASES.get(factor)
    if canonical:
        factor = canonical

    # 规范化匹配：去掉常见噪声字符
    norm = re.sub(r"[_\-\.]", "", factor).lower()
    if norm in motif_names:
        return "ok", ["vocab_hit:exact"]

    # 复合体拆分：词表名 'gata1::tal1' → {gata1, tal1}
    parts = {re.sub(r"[_\-\.]", "", p).lower() for p in norm.split("::")}
    if any(p in motif_names for p in parts):
        return "ok", ["vocab_hit:complex"]

    # 家族前缀匹配：factor 是某词表名的前缀（CREB ⊂ CREB1）
    if len(norm) >= 4:
        if any(n.startswith(norm) for n in motif_names):
            return "ok", ["vocab_hit:family_prefix"]

    return "unverified", [f"vocab_miss:{factor}"]


def is_chromatin_enzyme(factor: str) -> bool:
    """判断因子是否为染色质修饰酶（JASPAR 不收，豁免词表）。"""
    if factor in ENZYME_EXACT:
        return True
    up = factor.upper()
    return any(h in up for h in ENZYME_HINTS)


def curate_record(rec: dict, vocab: dict, motif_names: set) -> dict:
    """对单条 KB 记录做全量校验。返回 {record, status, reasons, bucket}。

    双 schema 支持：
    - MPRA 型：entities.motif 形如 MAxxxx.x（自带 JASPAR ID）→ 专用校验
    - 标准型：entities.factor + factor_type + cell_line
    """
    ent = rec.get("entities", {}) or {}
    cell_line = ent.get("cell_line", "")
    motif_id = ent.get("motif", "")
    factor = ent.get("factor", "")
    factor_type = ent.get("factor_type", "")
    reasons = []

    # ---- Schema C：GC 含量/序列特征分析（无 factor/motif，entities 含序列特征键） ----
    # 覆盖：gc_range / gc_optimal / gc_mean / feature / shape 等变体
    SEQ_FEATURE_KEYS = {"gc_range", "gc_optimal", "gc_suppressive", "gc_mean", "gc_std",
                        "feature", "shape", "n_enriched_motifs"}
    if not factor and not motif_id and (SEQ_FEATURE_KEYS & set(ent.keys())):
        if cell_line not in ALLOWED_CELL_LINES:
            if cell_line == "not_specified":
                return {"record": rec, "status": "rejected", "bucket": "gc_cell_line_not_specified",
                        "reasons": ["strict_cell_line:not_specified（需人工裁决）"]}
            return {"record": rec, "status": "rejected", "bucket": "gc_cell_line_out_of_scope",
                    "reasons": [f"strict_cell_line:{cell_line}（非三系）"]}
        return {"record": rec, "status": "kept", "bucket": "kept_sequence_feature",
                "reasons": [f"sequence_feature:keys:{sorted(SEQ_FEATURE_KEYS & set(ent.keys()))}"]}

    # ---- Schema A：MPRA 记录（tf + motif + cell_line） ----
    if MPRA_MOTIF_RE.match(str(motif_id)):
        # 1) 细胞系严格模式
        if cell_line not in ALLOWED_CELL_LINES:
            if cell_line == "not_specified":
                return {"record": rec, "status": "rejected", "bucket": "mpra_cell_line_not_specified",
                        "reasons": ["strict_cell_line:not_specified（需人工裁决）"]}
            return {"record": rec, "status": "rejected", "bucket": "mpra_cell_line_out_of_scope",
                    "reasons": [f"strict_cell_line:{cell_line}（非三系）"]}
        # 2) motif 必须在 JASPAR 词表
        if motif_id in vocab:
            return {"record": rec, "status": "kept", "bucket": "kept_mpra",
                    "reasons": [f"mpra:motif_vocab_hit:{motif_id}"]}
        return {"record": rec, "status": "rejected", "bucket": "mpra_motif_missing",
                "reasons": [f"mpra:motif_not_in_vocab:{motif_id}"]}

    # ---- Schema B：标准记录（factor + factor_type + cell_line） ----
    # 1) 细胞系严格模式
    if cell_line not in ALLOWED_CELL_LINES:
        if cell_line == "not_specified":
            status, bucket = "rejected", "cell_line_not_specified"
            reasons.append("strict_cell_line:not_specified（需人工裁决）")
        else:
            status, bucket = "rejected", "cell_line_out_of_scope"
            reasons.append(f"strict_cell_line:{cell_line}（非三系）")
        return {"record": rec, "status": status, "bucket": bucket, "reasons": reasons}

    # 2) 噪声名单（明确非 TF 实体）
    if factor in NOISE_FACTORS:
        return {"record": rec, "status": "rejected", "bucket": "factor_noise",
                "reasons": [f"noise:{factor}"]}

    # 3) 因子类型校验（路径 A：标记/状态禁止入 factor 槽）
    fstatus, freasons = classify_factor(factor)
    reasons += freasons
    if fstatus == "wrong_type":
        return {"record": rec, "status": "rejected", "bucket": "factor_wrong_type", "reasons": reasons}

    # 4) 染色质修饰酶：豁免词表，修正 factor_type
    if is_chromatin_enzyme(factor):
        ent["factor_type"] = "enzyme"
        return {"record": rec, "status": "kept", "bucket": "kept_enzyme",
                "reasons": ["enzyme:chromatin_modifier（豁免 JASPAR）"]}

    # 5) 词表校验（TF/motif 必须命中 JASPAR）
    if fstatus == "ok":
        vstatus, vreasons = check_vocab(factor, factor_type, vocab, motif_names)
        reasons += vreasons
        if vstatus == "unverified":
            # 氨基酸突变（S87N/R89W）→ 归为噪声，不是单纯词表未命中
            if MUTATION_PATTERN.match(factor):
                return {"record": rec, "status": "rejected", "bucket": "factor_mutation",
                        "reasons": [f"mutation_entity:{factor}（非 TF 名）"]}
            return {"record": rec, "status": "rejected", "bucket": "factor_vocab_miss", "reasons": reasons}

    # 6) 通过（unverified 启发式但无硬伤 → 仍打回人工，因为需二次确认）
    if fstatus == "unverified":
        return {"record": rec, "status": "rejected", "bucket": "factor_phrase", "reasons": reasons}

    return {"record": rec, "status": "kept", "bucket": "kept", "reasons": reasons}


def main():
    ap = argparse.ArgumentParser(description="知识库受控词表校验 + 细胞系严格过滤 + 打回队列")
    ap.add_argument("--kb", default="data/processed/knowledge_base_clean_final.jsonl")
    ap.add_argument("--vocab", default="data/vocab/jaspar2024_vertebrates_core.json")
    ap.add_argument("--out", default="tmp/curation")
    args = ap.parse_args()

    out_dir = Path(args.out)
    out_dir.mkdir(parents=True, exist_ok=True)

    # 加载词表
    vocab_data = json.load(open(args.vocab))
    vocab = {mid: info["name"] for mid, info in vocab_data["motifs"].items()}
    motif_names = {re.sub(r"[_\-\.]", "", n).lower() for n in vocab.values()}

    # 加载 KB
    recs = [json.loads(l) for l in open(args.kb) if l.strip()]
    print(f"输入 KB: {len(recs)} 条")

    results = [curate_record(r, vocab, motif_names) for r in recs]
    kept = [r for r in results if r["status"] == "kept"]
    rejected = [r for r in results if r["status"] == "rejected"]

    # 写输出
    with open(out_dir / "kb_curated.jsonl", "w") as f:
        for r in kept:
            f.write(json.dumps(r["record"], ensure_ascii=False) + "\n")
    with open(out_dir / "kb_rejected.jsonl", "w") as f:
        for r in rejected:
            f.write(json.dumps({
                "record": r["record"],
                "reject_bucket": r["bucket"],
                "reasons": r["reasons"],
            }, ensure_ascii=False) + "\n")

    # 汇总报告
    buckets = Counter(r["bucket"] for r in results)
    report = {
        "input_count": len(recs),
        "kept_count": len(kept),
        "rejected_count": len(rejected),
        "buckets": dict(buckets),
        "bucket_details": {},
    }
    # 每个 bucket 列几个代表性原因
    for b in buckets:
        sample_reasons = Counter()
        for r in results:
            if r["bucket"] == b:
                for why in r["reasons"]:
                    sample_reasons[why] += 1
        report["bucket_details"][b] = dict(sample_reasons.most_common(8))

    with open(out_dir / "kb_curation_report.json", "w") as f:
        json.dump(report, f, ensure_ascii=False, indent=1)

    print(f"通过: {len(kept)} | 打回: {len(rejected)}")
    print("分桶:")
    for b, c in buckets.most_common():
        print(f"  {b}: {c}")
    print(f"报告: {out_dir}/kb_curation_report.json")


if __name__ == "__main__":
    main()
