# -*- coding: utf-8 -*-
"""
知识库验证脚本（validate step）

输入: data/processed/knowledge_base_clean.jsonl (Module 2 清洗后的知识库, 631 条)
输出:
  - data/processed/validation_report.jsonl   逐条验证结果（复核清单，含 Qwen 裁决）
  - data/processed/validation_summary.json   汇总统计

验证内容（三层，对应 Module 1 → Module 2 之间的质量关卡）:
  [A] 统计可靠性回查（程序化，150 条统计记录，无需模型）
      - p_value ∈ [0,1] 且非 NaN
      - effect_size 有效性（数值、>0）
      - direction 与 p_value 显著性一致性（富集/缺失结论必须 p<0.05）
      - 推理链完整性
  [B] 知识一致性（程序化 + Qwen）
      - B1 程序化: 文献记录 effect ↔ evidence.direction 映射一致性
      - B2 Qwen:  推理链自洽性（reasoning_chain 各步骤是否逻辑自洽并导向 claim）
  [C] 外部事实核对（Qwen 三分类，481 条文献记录）
      - 对每条文献 claim，结合来源证据与推理链，
        判定 supported（证据直接支持）/ unsupported（证据矛盾或无关）/ unclear（证据不足）

验证结论字段（写入 report）:
  - stat_checks: [通过项列表]（程序化）
  - qwen_verdict: "supported" | "unsupported" | "unclear"
  - chain_consistent: true | false | null
  - qwen_reason: Qwen 给出的理由

用法:
  python scripts/curate/validate_knowledge_base.py                              # 全量验证
  python scripts/curate/validate_knowledge_base.py --limit 20                   # 仅前 20 条（含 LLM）
  python scripts/curate/validate_knowledge_base.py --no-llm                     # 跳过 Qwen（仅程序化）
  python scripts/curate/validate_knowledge_base.py --batch-size 8 --max-memory '{"0":"8GiB","1":"8GiB"}'
"""
from __future__ import annotations

import argparse
import json
import logging
import math
import sys
import time
from pathlib import Path

import torch

ROOT = Path(__file__).resolve().parent.parent.parent
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from src.extraction.extractor import BioExtractor, parse_json_response  # noqa: E402

logging.basicConfig(level=logging.INFO, format="%(levelname)s %(message)s")
logger = logging.getLogger("validate")

# ---------------- 配置 ----------------
# 文献记录 effect → evidence.direction 合法映射
# v2: 新增 relation 粒度枚举（与证据层级对齐：expression/binding/interaction/correlation）
EFFECT_DIRECTION_MAP = {
    "increases": "increases_activity",
    "decreases": "decreases_activity",
    "required_for": "required_for_activity",
    "modulates": "modulates_activity",
    "no_effect": "no_effect",
    # ---- v2 扩展：relation 粒度与证据层级对齐 ----
    "binds": "binds_target",                        # occupancy/binding 证据（#226 ATF4）
    "correlates": "correlates_with",                # 观察性相关（#323 accessibility）
    "co_occurs_with": "co_occurs_with",             # 共现、无方向（#292 methylation）
    "disrupts": "disrupts",                         # 变异破坏 motif（#470 caQTL）
    "promotes": "promotes_phenotype",               # 干预促表型（#28 AML 分化）
    "regulates_expression": "regulates_expression",  # 调节表达、方向未明（#399 OsDDE9）
    "positively_regulates": "positively_regulates_expression",  # 正调控表达（#199 ALKBH5）
    "negatively_regulates": "negatively_regulates_expression",  # 负调控表达（#54/#160/#224）
    "affects_interaction": "affects_interaction",   # 影响 promoter 互作（#170 GATA）
    "affects": "affects_binding",                   # 影响预测结合位点（#43 QTLs association）
}
# 需要 p 值显著的统计方向（富集/缺失类结论）
# 注: effect_size 负值合法（Cohen's d / logFC / 均值差）
SIGNIFICANT_DIRECTIONS = {"enriched_in_high", "depleted_in_high", "cell_specific_enriched"}

VALIDATION_SYSTEM_PROMPT = (
    "You are a rigorous biomedical knowledge verification expert. "
    "Given a scientific claim, its type, its evidence source, and a "
    "reasoning chain, you must judge: (1) whether the claim is supported "
    "by the provided evidence (verdict), and (2) whether the reasoning "
    "chain is logically self-consistent and leads to the claim "
    "(chain_consistent).\n"
    "STRICT RULES:\n"
    "- claim_type tells you what KIND of claim this is:\n"
    "    mechanistic : a factual statement about how a factor affects a "
    "                  regulatory element. verdict='supported' ONLY if the "
    "                  evidence literally states this relationship.\n"
    "    design_rule : a DESIGN RECOMMENDATION derived by inference from "
    "                  the evidence (e.g. 'when designing sequences, include "
    "                  motif X because it correlates with high activity'). "
    "                  For design_rule, you verify the PREMISE: does the "
    "                  evidence support the underlying relationship (e.g. "
    "                  motif X correlates with high activity) with the "
    "                  correct direction? The recommendation itself is an "
    "                  expected inference, NOT a statement the evidence must "
    "                  literally contain. If the evidence supports the "
    "                  premise with the same direction, use 'supported'.\n"
    "- verdict must be 'supported' ONLY if the evidence supports the claim "
    "  (literally for mechanistic; at the premise level with correct "
    "  direction for design_rule). If the evidence contradicts the claim "
    "  or is completely unrelated, use 'unsupported'. If the evidence is "
    "  insufficient to decide, use 'unclear'.\n"
    "- chain_consistent must be true if each step follows logically from "
    "  the previous one and the final step supports the claim. For "
    "  design_rule, a final inference step (evidence -> recommendation) is "
    "  legitimate and does NOT make the chain inconsistent. Only mark false "
    "  if a step truly contradicts the claim or is illogical.\n"
    "- reason: one sentence in Chinese explaining your judgment.\n"
    "Output ONLY a JSON object: {\"verdict\": \"supported|unsupported|unclear\", "
    "\"chain_consistent\": true|false, \"reason\": \"...\"}"
)

VALIDATION_USER_TEMPLATE = """【科学断言类型】
{claim_type}

【科学断言】
{claim}

【证据来源】
{source}

【推理链】
{chain}

请验证该断言：
1) verdict: 断言是否被提供的证据直接支持（mechanistic）或前提是否被证据以正确方向支持（design_rule）？只能是 "supported" / "unsupported" / "unclear"
2) chain_consistent: 推理链各步骤是否逻辑自洽并导向该断言？只能是 true / false
3) reason: 一句话中文理由

仅输出 JSON：{{"verdict": "...", "chain_consistent": true/false, "reason": "..."}}"""


def build_validation_messages(payload: str) -> list[dict]:
    """构造验证 Chat 消息（payload 为已经格式化好的用户内容）。"""
    return [
        {"role": "system", "content": VALIDATION_SYSTEM_PROMPT},
        {"role": "user", "content": payload},
    ]


# ---------------- [A] 统计可靠性回查（程序化） ----------------
def check_stat_record(r: dict) -> list[str]:
    """对统计记录做可靠性回查，返回通过项描述列表；发现问题抛异常由上层收集。"""
    ev = r.get("evidence", {})
    checks = []
    problems = []

    p = ev.get("p_value")
    if p is None:
        problems.append("p_value 缺失")
    else:
        try:
            p = float(p)
            if math.isnan(p):
                problems.append("p_value 为 NaN")
            elif not (0 <= p <= 1):
                problems.append(f"p_value 越界: {p}")
            else:
                checks.append(f"p_value={p:.3g} 在 [0,1] 内")
                # 方向与显著性一致性
                direction = ev.get("direction", "")
                if direction in SIGNIFICANT_DIRECTIONS and p >= 0.05:
                    problems.append(f"direction={direction} 但 p_value={p:.3g} ≥ 0.05, 富集结论不显著")
                elif direction in SIGNIFICANT_DIRECTIONS:
                    checks.append(f"direction={direction} 且 p<0.05, 显著性一致")
        except (TypeError, ValueError):
            problems.append(f"p_value 非数值: {p}")

    es = ev.get("effect_size")
    if es is None:
        checks.append("effect_size 缺失（允许，如非效应量类记录）")
    else:
        try:
            es = float(es)
            if math.isnan(es):
                problems.append("effect_size 为 NaN")
            else:
                # 负值合法: Cohen's d / logFC / 均值差 等可为负（负向效应）
                checks.append(f"effect_size={es:.4g} 数值有效（{'负值' if es < 0 else '非负'}）")
        except (TypeError, ValueError):
            problems.append(f"effect_size 非数值: {es}")

    # 推理链完整性（≥2 步: design_rule 类精简链也合法）
    chain = r.get("reasoning_chain") or []
    if len(chain) < 2:
        problems.append(f"推理链过短: {len(chain)} 步")
    else:
        checks.append(f"推理链完整: {len(chain)} 步")

    if problems:
        raise ValueError("; ".join(problems))
    return checks


# ---------------- [B1] 知识一致性（程序化） ----------------
def check_lit_consistency(r: dict) -> list[str]:
    """文献记录 effect ↔ direction 一致性检查，返回通过项；问题抛异常。"""
    e = r.get("entities", {})
    ev = r.get("evidence", {})
    effect = e.get("effect")
    direction = ev.get("direction")
    checks, problems = [], []

    if effect is None:
        problems.append("entities.effect 缺失")
    if direction is None:
        problems.append("evidence.direction 缺失")
    if effect is not None and direction is not None:
        expected = EFFECT_DIRECTION_MAP.get(effect)
        if expected is None:
            problems.append(f"effect 值非法: {effect}")
        elif direction != expected:
            problems.append(f"effect={effect} 与 direction={direction} 不一致(期望 {expected})")
        else:
            checks.append(f"effect↔direction 一致 ({effect}→{direction})")

    # 必填实体
    for key in ("factor", "regulatory_element"):
        if not e.get(key):
            problems.append(f"entities.{key} 缺失")

    chain = r.get("reasoning_chain") or []
    if len(chain) < 2:
        problems.append(f"推理链过短: {len(chain)} 步")
    else:
        checks.append(f"推理链完整: {len(chain)} 步")

    if problems:
        raise ValueError("; ".join(problems))
    return checks


# ---------------- 输出 ----------------
def load_rows(path: Path) -> list[dict]:
    rows = []
    with open(path, "r", encoding="utf-8") as f:
        for line in f:
            line = line.strip()
            if line:
                rows.append(json.loads(line))
    return rows


def format_payload(r: dict) -> str:
    """格式化单条记录的验证载荷（含类型 / claim / 来源 / 推理链）。"""
    chain = r.get("reasoning_chain") or []
    chain_text = "\n".join(f"{i+1}. {s}" for i, s in enumerate(chain))
    ev = r.get("evidence", {})
    source = ev.get("source", "未知")
    return VALIDATION_USER_TEMPLATE.format(
        claim_type=r.get("claim_type", "unknown"),
        claim=r.get("claim", ""),
        source=source,
        chain=chain_text,
    )


def main():
    ap = argparse.ArgumentParser(description="知识库验证（统计回查 + Qwen 外部事实核对）")
    ap.add_argument("--input", default=str(ROOT / "data/processed/knowledge_base_clean.jsonl"))
    ap.add_argument("--output", default=str(ROOT / "data/processed/validation_report.jsonl"))
    ap.add_argument("--summary", default=str(ROOT / "data/processed/validation_summary.json"))
    ap.add_argument("--model", default=str(ROOT / "Qwen2.5-7B-Instruct"))
    ap.add_argument("--batch-size", type=int, default=8)
    ap.add_argument("--max-new-tokens", type=int, default=256)
    ap.add_argument("--temperature", type=float, default=0.1)
    ap.add_argument("--limit", type=int, default=0, help="只处理前 N 条（0=全部）")
    ap.add_argument("--no-llm", action="store_true", help="跳过 Qwen 验证（仅程序化检查）")
    ap.add_argument("--max-memory", default=None, help="GPU 显存上限 JSON，如 '{\"0\":\"8GiB\",\"1\":\"8GiB\"}'")
    args = ap.parse_args()

    max_memory = None
    if args.max_memory:
        try:
            raw = json.loads(args.max_memory)
            max_memory = {int(k): v for k, v in raw.items()}
        except (json.JSONDecodeError, ValueError):
            logger.error("--max-memory 不是合法 JSON: %s", args.max_memory)
            sys.exit(1)

    inp = Path(args.input)
    rows = load_rows(inp)
    if args.limit > 0:
        rows = rows[: args.limit]
    logger.info("输入: %s (%d 条, limit=%d)", inp, len(rows), args.limit or len(rows))

    out_path = Path(args.output)
    out_path.parent.mkdir(parents=True, exist_ok=True)

    stats = {"input": len(rows), "no_llm": args.no_llm, "started_at": time.strftime("%Y-%m-%d %H:%M:%S")}
    report_lines: list[dict] = []
    stat_fail = lit_consistency_fail = 0

    # ---- [A] + [B1] 程序化验证 ----
    llm_items = []  # 待 Qwen 验证的载荷
    for i, r in enumerate(rows):
        is_stat = "factor" not in r.get("entities", {})
        entry = {"idx": i, "record_type": "stat" if is_stat else "lit",
                 "claim": r.get("claim", "")[:120]}
        try:
            if is_stat:
                entry["stat_checks"] = check_stat_record(r)
            else:
                entry["stat_checks"] = check_lit_consistency(r)
            entry["checks_passed"] = True
        except ValueError as e:
            entry["checks_passed"] = False
            entry["stat_errors"] = str(e)
            if is_stat:
                stat_fail += 1
            else:
                lit_consistency_fail += 1
        if not is_stat:
            llm_items.append((i, r))
        report_lines.append(entry)

    logger.info("程序化检查: 统计记录失败 %d 条, 文献一致性失败 %d 条", stat_fail, lit_consistency_fail)
    stats["stat_checks_failed"] = stat_fail
    stats["lit_consistency_failed"] = lit_consistency_fail

    # ---- [B2] + [C] Qwen 验证（文献记录）----
    if not args.no_llm and llm_items:
        # 断点续跑：已输出文件中带 qwen_verdict 的 idx 跳过
        done = set()
        if out_path.exists():
            for line in out_path.read_text(encoding="utf-8").splitlines():
                line = line.strip()
                if not line:
                    continue
                try:
                    rec = json.loads(line)
                except json.JSONDecodeError:
                    continue
                if rec.get("qwen_verdict"):
                    done.add(rec["idx"])
        pending = [(i, r) for i, r in llm_items if i not in done]
        logger.info("Qwen 验证: 共 %d 条文献记录, 待处理 %d 条", len(llm_items), len(pending))

        if pending:
            logger.info("加载模型: %s ...", args.model)
            extractor = BioExtractor(model_path=args.model, max_memory=max_memory)
            start = time.time()
            batch_size = args.batch_size

            for k in range(0, len(pending), batch_size):
                batch = pending[k : k + batch_size]
                items = [{"pmid": str(i), "abstract": format_payload(r)} for i, r in batch]

                current_bs = batch_size
                batch_results = None
                while batch_results is None:
                    try:
                        batch_results = extractor._extract_one_batch(
                            items[:current_bs],
                            max_new_tokens=args.max_new_tokens,
                            temperature=args.temperature,
                            use_chat_template=True,
                            prompt_builder=build_validation_messages,
                            fields=("verdict", "chain_consistent", "reason"),
                        )
                    except torch.cuda.OutOfMemoryError:
                        torch.cuda.empty_cache()
                        current_bs = max(1, current_bs // 2)
                        logger.warning("OOM! 降批至 %d 重试", current_bs)
                        if current_bs == 1:
                            logger.warning("单条仍 OOM, 跳过 idx=%s", batch[0][0])
                            batch_results = []

                # 解析并写盘（立即追加）
                with open(out_path, "a", encoding="utf-8") as f:
                    for (idx, _r), res in zip(batch, batch_results):
                        verdict, chain_ok, reason = None, None, None
                        parsed = parse_json_response(res.get("raw_output", ""))
                        if parsed:
                            verdict = parsed.get("verdict")
                            chain_ok = parsed.get("chain_consistent")
                            reason = parsed.get("reason")
                        entry = {
                            "idx": idx,
                            "record_type": "lit",
                            "claim": _r.get("claim", "")[:120],
                            "qwen_verdict": verdict,
                            "chain_consistent": chain_ok,
                            "qwen_reason": reason,
                            "raw_output": res.get("raw_output", "")[:200],
                        }
                        f.write(json.dumps(entry, ensure_ascii=False) + "\n")
                        # 更新内存中已有的程序化条目（保留 stat_checks/checks_passed），避免重复写入
                        for _e in report_lines:
                            if _e.get("idx") == idx and _e.get("record_type") == "lit":
                                _e.update({k: v for k, v in entry.items() if k not in ("idx", "record_type")})
                                break
                        else:
                            report_lines.append(entry)

                done_pct = min(k + batch_size, len(pending))
                if (k // batch_size) % 5 == 0 or done_pct == len(pending):
                    elapsed = time.time() - start
                    logger.info("Qwen 进度: %d/%d (%.1f 秒, 累计 %.1f 分钟)",
                                done_pct, len(pending), elapsed, elapsed / 60)

            elapsed = time.time() - start
            logger.info("Qwen 验证完成: %d 条, 耗时 %.1f 分钟 (%.2f 秒/条)",
                        len(pending), elapsed / 60, elapsed / max(len(pending), 1))
            del extractor
            if torch.cuda.is_available():
                torch.cuda.empty_cache()

    # ---- 汇总 ----
    if not args.no_llm:
        # 重新读全量输出（含断点续跑历史）
        final_entries = {}
        if out_path.exists():
            for line in out_path.read_text(encoding="utf-8").splitlines():
                line = line.strip()
                if not line:
                    continue
                try:
                    rec = json.loads(line)
                except json.JSONDecodeError:
                    continue
                final_entries[rec["idx"]] = rec
        lit_rows = [r for r in rows if "factor" in r.get("entities", {})]
        verdicts = [final_entries[i]["qwen_verdict"] for i, _ in llm_items if i in final_entries and final_entries[i].get("qwen_verdict")]
        chains = [final_entries[i]["chain_consistent"] for i, _ in llm_items if i in final_entries and final_entries[i].get("chain_consistent") is not None]
        stats["lit_total"] = len(lit_rows)
        stats["qwen_verified"] = len(verdicts)
        stats["qwen_verdict_dist"] = {v: verdicts.count(v) for v in ("supported", "unsupported", "unclear")}
        stats["qwen_chain_consistent"] = sum(1 for c in chains if c is True)
        stats["qwen_chain_total"] = len(chains)
        stats["qwen_support_rate"] = round(verdicts.count("supported") / max(len(verdicts), 1), 4)
        stats["qwen_chain_consistent_rate"] = round(sum(1 for c in chains if c is True) / max(len(chains), 1), 4)

    stats["finished_at"] = time.strftime("%Y-%m-%d %H:%M:%S")
    stats["principles"] = [
        "验证遵循三层标准：①统计可靠性（p_value 有效性+方向显著性一致性+effect_size 合理性）；"
        "②知识一致性（effect↔direction 映射+推理链完整性）；"
        "③外部事实核对（Qwen2.5-7B 本地模型对文献 claim 做 supported/unsupported/unclear 三分类，"
        "并结合推理链自洽性判断）。verdict=unsupported 或 chain_consistent=false 的记录"
        "进入人工复核清单。"
    ]

    # 断点续跑：合并历史文件中已完成的 Qwen 结果（避免本次跳过时丢失）
    if not args.no_llm and out_path.exists():
        _hist = {}
        for _line in out_path.read_text(encoding="utf-8").splitlines():
            _line = _line.strip()
            if not _line:
                continue
            try:
                _rec = json.loads(_line)
            except json.JSONDecodeError:
                continue
            if _rec.get("qwen_verdict"):
                _hist[_rec["idx"]] = _rec
        for _e in report_lines:
            if _e.get("record_type") == "lit" and not _e.get("qwen_verdict") and _e.get("idx") in _hist:
                _h = _hist[_e["idx"]]
                _e.update({k: v for k, v in _h.items() if k not in ("idx", "record_type")})

    # 输出 report（程序化 + Qwen 合并结果，每 idx 仅一条）
    with open(out_path, "w", encoding="utf-8") as f:
        for e in report_lines:
            f.write(json.dumps(e, ensure_ascii=False) + "\n")
    with open(Path(args.summary), "w", encoding="utf-8") as f:
        json.dump(stats, f, ensure_ascii=False, indent=2)

    logger.info("验证报告: %s", out_path)
    logger.info("汇总: %s", json.dumps(stats, ensure_ascii=False, default=str))
    if stats.get("qwen_verdict_dist", {}).get("unsupported"):
        logger.warning("发现 %d 条 unsupported（见 report, 人工复核）",
                       stats["qwen_verdict_dist"]["unsupported"])
    if stats.get("qwen_chain_consistent", 0) < stats.get("qwen_chain_total", 0):
        logger.warning("发现 %d 条推理链不自洽（见 report, 人工复核）",
                       stats.get("qwen_chain_total", 0) - stats.get("qwen_chain_consistent", 0))


if __name__ == "__main__":
    main()
