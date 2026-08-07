"""Module 1 · DNA 形状查询表提取。

从 Bioconductor DNAshapeR 包源码（src/TableCompiler.cpp）中提取内置的
标准 DNAshape 五核苷酸查询表（987 条，覆盖全部 1024 个 ACGT 5-mer：
非回文 5-mer 与其反向互补共享一个条目）。

每行 91 字段: [5mer] + 90 数值，关键字段（按 load_data_from_vector 布局）:
    dv[0]    MGW   (minor groove width)  平均值
    dv[6]    ProT  (propeller twist)     平均值
    dv[15]   roll1  (5-mer 第 2 个 step 的 roll)
    dv[18]   roll2  (5-mer 第 1 个 step 的 roll)
    dv[21]   twist1 (5-mer 第 2 个 step 的 HelT)
    dv[24]   twist2 (5-mer 第 1 个 step 的 HelT)
后续字段为 MC/XRC/MD 模拟来源的补充参数（本分析不使用）。

用法（项目根目录）:
    python scripts/extract/extract_dna_shape_table.py [--source <TableCompiler.cpp>] [--out <tsv>]
"""

from __future__ import annotations

import argparse
import itertools
import logging
import re
import sys
from pathlib import Path

logging.basicConfig(level=logging.INFO, format="%(asctime)s | %(levelname)s | %(message)s", datefmt="%H:%M:%S")
logger = logging.getLogger("extract_dna_shape_table")

ROOT = Path(__file__).resolve().parent.parent

# 服务器上 DNAshapeR 包源码路径（可复现源）
DEFAULT_SOURCE = "/tmp/d140/DNAshapeR/src/TableCompiler.cpp"


def extract_query_table(source: Path) -> list[list[str]]:
    """从 C++ 源码解析 QueryTable[]，返回每行 91 个字段的列表。"""
    text = source.read_text(encoding="utf-8", errors="replace")
    m = re.search(r"const std::string QueryTable\[\] ={(.*?)};", text, re.S)
    if not m:
        raise ValueError(f"未找到 QueryTable[]: {source}")
    raw_rows = re.findall(r'"([^"]+)"', m.group(1))
    rows: list[list[str]] = []
    for r in raw_rows:
        parts = r.split()
        if len(parts) != 91:
            logger.warning("跳过异常行 (%d 字段): %s...", len(parts), parts[0] if parts else "?")
            continue
        rows.append(parts)
    return rows


def validate(rows: list[list[str]]) -> None:
    """校验: 987 条、5-mer 唯一、1024 标准 5-mer 全覆盖（含反向互补）。"""
    pents = [r[0] for r in rows]
    if len(pents) != len(set(pents)):
        raise ValueError("存在重复 5-mer")
    all5 = set("".join(p) for p in itertools.product("ACGT", repeat=5))
    std = set(p for p in pents if set(p) <= set("ACGT"))
    comp = str.maketrans("ACGT", "TGCA")

    def rc(s: str) -> str:
        return s.translate(comp)[::-1]

    missing = all5 - std
    bad = [m for m in missing if rc(m) not in std]
    if bad:
        raise ValueError(f"以下 5-mer 及其反向互补均不在表中: {bad[:10]}")
    logger.info("校验通过: %d 条, %d 个标准 5-mer, %d 个简并(M/Q)条目",
                len(rows), len(std), len(pents) - len(std))


def main() -> None:
    parser = argparse.ArgumentParser(description="提取 DNAshape 五核苷酸查询表")
    parser.add_argument("--source", default=DEFAULT_SOURCE, help="DNAshapeR 的 TableCompiler.cpp 路径")
    parser.add_argument("--out", default="data/processed/dna_shape_querytable.tsv")
    args = parser.parse_args()

    src = Path(args.source)
    if not src.exists():
        logger.error("源码不存在: %s\n请先解压 DNAshapeR 包（本脚本从 TableCompiler.cpp 提取内置表）。", src)
        sys.exit(1)

    rows = extract_query_table(src)
    logger.info("解析到 %d 条查询表条目", len(rows))
    validate(rows)

    out_path = ROOT / args.out
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with open(out_path, "w", encoding="utf-8") as f:
        for parts in rows:
            f.write("\t".join(parts) + "\n")
    logger.info("已保存: %s (%d 行)", out_path, len(rows))


if __name__ == "__main__":
    main()
