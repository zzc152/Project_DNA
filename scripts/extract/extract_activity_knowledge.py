"""文献知识源（Source 2）抽取入口：三系调控元件活性影响原因。

使用本地 Qwen2.5-7B-Instruct 从 abstracts_activity.jsonl 中抽取
"什么因素/机制影响调控元件活性"的结构化知识（findings 列表）。

用法（项目根目录 /workspace/zzc/BioDesign-Agent 下运行）:
    python scripts/extract/extract_activity_knowledge.py                          # 全量抽取
    python scripts/extract/extract_activity_knowledge.py --limit 10               # 仅前 10 条
    python scripts/extract/extract_activity_knowledge.py --input data/raw/abstracts_activity.jsonl
"""

import argparse
import json
import logging
import sys
import time
from pathlib import Path

import torch

# 确保项目根目录可被 import
ROOT = Path(__file__).resolve().parent.parent
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from src.extraction.extractor import BioExtractor  # noqa: E402
from src.extraction.prompts_activity import build_activity_messages  # noqa: E402

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s | %(levelname)s | %(message)s",
    datefmt="%H:%M:%S",
)
logger = logging.getLogger("extract_activity_knowledge")

# 输出字段（BioExtractor 动态字段支持）
ACTIVITY_FIELDS = ("findings",)


def load_abstracts(path: Path) -> list[dict]:
    """读取 jsonl 摘要文件，每行一个 JSON 对象。"""
    items: list[dict] = []
    with open(path, "r", encoding="utf-8") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            try:
                items.append(json.loads(line))
            except json.JSONDecodeError:
                logger.warning("跳过无法解析的行: %s", line[:100])
    return items


def main() -> None:
    parser = argparse.ArgumentParser(description="从 PubMed 摘要抽取三系调控元件活性影响原因知识")
    parser.add_argument("--input", default="data/raw/abstracts_activity.jsonl", help="输入摘要 jsonl")
    parser.add_argument("--output", default="data/processed/activity_extractions.jsonl", help="输出 jsonl")
    parser.add_argument("--model", default="./Qwen2.5-7B-Instruct", help="本地模型路径")
    parser.add_argument("--batch-size", type=int, default=8, help="批量推理大小")
    parser.add_argument("--max-new-tokens", type=int, default=512)
    parser.add_argument("--temperature", type=float, default=0.1)
    parser.add_argument("--limit", type=int, default=0, help="只处理前 N 条（0=全部）")
    parser.add_argument("--no-chat-template", action="store_true", help="使用字符串拼接提示词")
    parser.add_argument(
        "--max-memory",
        default=None,
        help="每张 GPU 显存上限（JSON，如 '{\"0\": \"8GiB\", \"1\": \"4GiB\"}'）；超出部分自动放 CPU",
    )
    args = parser.parse_args()

    # 解析 max_memory（键自动转整数）
    max_memory = None
    if args.max_memory:
        try:
            raw = json.loads(args.max_memory)
            max_memory = {int(k): v for k, v in raw.items()}
        except (json.JSONDecodeError, ValueError):
            logger.error("--max-memory 不是合法 JSON: %s", args.max_memory)
            sys.exit(1)

    input_path = Path(args.input)
    if not input_path.exists():
        logger.error("输入文件不存在: %s（请在项目根目录运行）", input_path)
        sys.exit(1)
    output_path = Path(args.output)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    items = load_abstracts(input_path)
    if args.limit > 0:
        items = items[: args.limit]
    logger.info("共加载 %d 条摘要", len(items))
    if not items:
        logger.error("没有可处理的摘要，退出")
        sys.exit(1)

    extractor = BioExtractor(model_path=args.model, max_memory=max_memory)

    # 断点续跑：只跳过解析成功的条目（parsed=true）；解析失败的允许下次重试
    done_pmids: set[str] = set()
    n_failed_prev = 0
    if output_path.exists():
        for line in output_path.read_text(encoding="utf-8").splitlines():
            line = line.strip()
            if not line:
                continue
            try:
                rec = json.loads(line)
            except json.JSONDecodeError:
                continue
            if rec.get("parsed"):
                done_pmids.add(str(rec.get("pmid")))
            else:
                n_failed_prev += 1
        logger.info("断点续跑：已有 %d 条解析成功（跳过），%d 条解析失败将重试",
                    len(done_pmids), n_failed_prev)

    pending = [it for it in items if str(it.get("pmid")) not in done_pmids]
    logger.info("待处理 %d 条（batch_size=%d, temperature=%.1f）...",
                len(pending), args.batch_size, args.temperature)
    if not pending:
        logger.info("全部已完成，无需抽取")
        sys.exit(0)

    results: list[dict] = []
    start = time.time()
    batch_size = args.batch_size
    for i in range(0, len(pending), batch_size):
        batch = pending[i : i + batch_size]

        # OOM 自动降批重试（不中断整个进程；若降到 1 仍失败则跳过该条）
        current_bs = batch_size
        batch_results = None
        while batch_results is None:
            try:
                batch_results = extractor._extract_one_batch(
                    batch[:current_bs],
                    max_new_tokens=args.max_new_tokens,
                    temperature=args.temperature,
                    use_chat_template=not args.no_chat_template,
                    prompt_builder=build_activity_messages,
                    fields=ACTIVITY_FIELDS,
                )
            except torch.cuda.OutOfMemoryError:
                torch.cuda.empty_cache()
                current_bs = max(1, current_bs // 2)
                logger.warning("OOM! 降批至 %d 重试 (pmid=%s...)", current_bs, batch[0].get("pmid"))
                if current_bs == 1:
                    logger.warning("单条仍 OOM，跳过 pmid=%s", batch[0].get("pmid"))
                    batch_results = []

        # 立即追加写盘（断点续跑基础）
        with open(output_path, "a", encoding="utf-8") as f:
            for r in batch_results:
                f.write(json.dumps(r, ensure_ascii=False) + "\n")

        results.extend(batch_results)

    # 统计
    n_ok = sum(1 for r in results if r.get("parsed"))
    n_findings = sum(len(r.get("findings") or []) for r in results)
    n_empty = sum(1 for r in results if not (r.get("findings") or []))
    elapsed = time.time() - start
    logger.info("=" * 70)
    logger.info("抽取完成: %d 条（解析成功 %d，含 findings %d 条，空结果 %d 条）",
                len(results), n_ok, n_findings, n_empty)
    logger.info("耗时 %.1f 分钟（%.1f 秒/条）", elapsed / 60, elapsed / max(len(results), 1))
    logger.info("输出文件: %s", output_path)
    logger.info("=" * 70)


if __name__ == "__main__":
    main()
