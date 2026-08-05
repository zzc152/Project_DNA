"""合并多个 PubMed 摘要 jsonl 并按 pmid 去重。

用法（在项目根目录运行）:
    python scripts/dedupe_abstracts.py \
        --input data/raw/abstracts.jsonl data/raw/abstracts_targeted.jsonl \
        --output data/raw/abstracts_merged.jsonl
"""

import argparse
import json
import logging
from pathlib import Path

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s | %(levelname)s | %(message)s",
    datefmt="%H:%M:%S",
)
logger = logging.getLogger("dedupe_abstracts")


def load_jsonl(path: Path) -> list[dict]:
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
    parser = argparse.ArgumentParser(description="合并并去重 PubMed 摘要")
    parser.add_argument("--input", nargs="+", required=True, help="输入 jsonl 文件（可多个）")
    parser.add_argument("--output", default="data/raw/abstracts_merged.jsonl", help="输出文件")
    args = parser.parse_args()

    seen: dict[str, dict] = {}
    for p in args.input:
        items = load_jsonl(Path(p))
        logger.info("读取 %s: %d 条", p, len(items))
        for it in items:
            pmid = str(it.get("pmid"))
            if pmid and pmid not in seen:
                seen[pmid] = it

    logger.info("合并去重后共 %d 条", len(seen))
    out = Path(args.output)
    out.parent.mkdir(parents=True, exist_ok=True)
    with open(out, "w", encoding="utf-8") as f:
        for it in seen.values():
            f.write(json.dumps(it, ensure_ascii=False) + "\n")
    logger.info("已写入 %s", out)


if __name__ == "__main__":
    main()
