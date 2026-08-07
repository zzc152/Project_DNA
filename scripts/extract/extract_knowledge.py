"""Module 1 知识抽取入口脚本。

用法（在项目根目录 /workspace/zzc/BioDesign-Agent 下运行）:
    python scripts/extract/extract_knowledge.py                          # 全量抽取
    python scripts/extract/extract_knowledge.py --limit 10               # 仅前 10 条（小批量测试）
    python scripts/extract/extract_knowledge.py --input data/samples/test_abstracts.jsonl
"""

import argparse
import json
import logging
import sys
import time
from collections import Counter
from pathlib import Path

import torch

# 确保项目根目录可被 import
ROOT = Path(__file__).resolve().parent.parent
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from src.extraction.extractor import BioExtractor  # noqa: E402

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s | %(levelname)s | %(message)s",
    datefmt="%H:%M:%S",
)
logger = logging.getLogger("extract_knowledge")


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
    parser = argparse.ArgumentParser(description="从 PubMed 摘要抽取 DNA 调控知识")
    parser.add_argument("--input", default="data/raw/abstracts.jsonl", help="输入摘要 jsonl")
    parser.add_argument("--output", default="data/processed/raw_extractions.jsonl", help="输出 jsonl")
    parser.add_argument("--model", default="./Qwen2.5-7B-Instruct", help="本地模型路径")
    parser.add_argument("--batch-size", type=int, default=8, help="批量推理大小")
    parser.add_argument("--max-new-tokens", type=int, default=256)
    parser.add_argument("--temperature", type=float, default=0.1)
    parser.add_argument("--limit", type=int, default=0, help="只处理前 N 条（0=全部）")
    parser.add_argument("--no-chat-template", action="store_true", help="使用字符串拼接提示词")
    parser.add_argument(
        "--max-memory",
        default=None,
        help="每张 GPU 显存上限（JSON，如 '{\"0\": \"8GiB\", \"1\": \"4GiB\"}'）；超出部分自动放 CPU",
    )
    args = parser.parse_args()

    # 解析 max_memory（键自动转整数，如 '{"0": "8GiB"}' → {0: "8GiB"})
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

    # 断点续跑：输出文件里已有的 pmid 直接跳过
    done_pmids: set[str] = set()
    if output_path.exists():
        for line in output_path.read_text(encoding="utf-8").splitlines():
            line = line.strip()
            if not line:
                continue
            try:
                done_pmids.add(str(json.loads(line).get("pmid")))
            except json.JSONDecodeError:
                pass
        logger.info("断点续跑：已有 %d 条完成，跳过", len(done_pmids))

    pending = [it for it in items if str(it.get("pmid")) not in done_pmids]
    logger.info("待处理 %d 条（batch_size=%d, temperature=%.1f）...", len(pending), args.batch_size, args.temperature)
    if not pending:
        logger.info("全部已完成，无需抽取")
        sys.exit(0)

    results: list[dict] = []
    start = time.time()
    batch_size = args.batch_size
    for i in range(0, len(pending), batch_size):
        batch = pending[i : i + batch_size]

        # OOM 自动降批重试（不中断整个进程；若降到 1 仍失败则跳过该条，绝不 crash）
        current_bs = batch_size
        batch_results = None
        while batch_results is None:
            try:
                batch_results = extractor._extract_one_batch(
                    batch[:current_bs],
                    max_new_tokens=args.max_new_tokens,
                    temperature=args.temperature,
                    use_chat_template=not args.no_chat_template,
                )
            except torch.cuda.OutOfMemoryError:
                torch.cuda.empty_cache()
                current_bs = max(1, current_bs // 2)
                logger.warning("OOM! 降批至 %d 重试 (pmid=%s...)", current_bs, batch[0].get("pmid"))
                if current_bs == 1:
                    # 单条也 OOM：跳过这一条，不阻塞整体
                    logger.warning("单条仍 OOM，跳过 pmid=%s", batch[0].get("pmid"))
                    batch_results = []

        # 立即追加写盘（断点续跑基础）
        with open(output_path, "a", encoding="utf-8") as f:
            for r in batch_results:
                f.write(json.dumps(r, ensure_ascii=False) + "\n")
        results.extend(batch_results)

        # 当前批次可能因 OOM 只处理了一部分，若 batch 未耗尽需继续下一小批
        if current_bs < batch_size and batch_results:
            remaining = batch[current_bs:]
            while remaining:
                try:
                    sub = extractor._extract_one_batch(
                        remaining[:1],
                        max_new_tokens=args.max_new_tokens,
                        temperature=args.temperature,
                        use_chat_template=not args.no_chat_template,
                    )
                    with open(output_path, "a", encoding="utf-8") as f:
                        for r in sub:
                            f.write(json.dumps(r, ensure_ascii=False) + "\n")
                    results.extend(sub)
                    remaining = remaining[1:]
                except torch.cuda.OutOfMemoryError:
                    torch.cuda.empty_cache()
                    logger.warning("OOM! 单条跳过 pmid=%s", remaining[0].get("pmid"))
                    remaining = remaining[1:]

        logger.info("进度: %d/%d", min(i + batch_size, len(pending)), len(pending))

    elapsed = time.time() - start

    # 统计信息
    total = len(results)
    success = sum(1 for r in results if r.get("parsed"))
    rate = success / total * 100 if total else 0
    tf_cnt = sum(len(r.get("tf", [])) for r in results)
    gene_cnt = sum(len(r.get("gene", [])) for r in results)
    motif_cnt = sum(len(r.get("motif", [])) for r in results)
    disease_cnt = sum(len(r.get("disease", [])) for r in results)
    relation_cnt = sum(1 for r in results if r.get("relation"))

    logger.info("=" * 60)
    logger.info("抽取统计:")
    logger.info("  总摘要数        : %d", total)
    logger.info("  成功抽取数      : %d (%.1f%%)", success, rate)
    logger.info("  转录因子实体数  : %d", tf_cnt)
    logger.info("  靶基因实体数    : %d", gene_cnt)
    logger.info("  结合基序实体数  : %d", motif_cnt)
    logger.info("  疾病实体数      : %d", disease_cnt)
    logger.info("  含调控关系数    : %d", relation_cnt)
    logger.info("  总耗时          : %.1f s (平均 %.2f s/条)", elapsed, elapsed / total if total else 0)
    logger.info("=" * 60)

    # 实体频率 TOP10
    all_entities = Counter()
    for r in results:
        for key in ("tf", "gene", "motif", "disease"):
            for ent in r.get(key, []):
                if ent:
                    all_entities[str(ent)] += 1
    logger.info("实体频率 TOP10: %s", all_entities.most_common(10))


if __name__ == "__main__":
    main()
