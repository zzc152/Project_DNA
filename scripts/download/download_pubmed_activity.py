"""从 PubMed 定向下载"三个细胞系调控元件活性影响原因"相关文献摘要。

与 download_pubmed_mpra.py 的差异：本脚本聚焦 *活性影响原因*（determinants /
mechanisms / factors that affect enhancer-promoter activity），按细胞系分主题：
  - K562（红系）: GATA/KLF/NFE2/FLI1 等红系因子如何决定增强子活性
  - HepG2（肝） : HNF/CEBP/FOXA 等肝富集因子如何决定增强子活性
  - SK-N-SH（神经）: NEUROG/ASCL/NEUROD/SOX 等神经因子如何决定增强子活性
  - 通用活性机制: 染色质状态 / 表观修饰 / 序列特征如何影响增强子活性

每个主题按 relevance 排序取前 per_query 条，合并按 pmid 去重后截断到 --total。

用法（项目根目录）:
    python scripts/download/download_pubmed_activity.py --output data/raw/abstracts_activity.jsonl
    python scripts/download/download_pubmed_activity.py --max-results 60 --total 120 --min-year 2015
"""

import argparse
import json
import logging
import re
import sys
import time
import xml.etree.ElementTree as ET
from pathlib import Path

import requests

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s | %(levelname)s | %(message)s",
    datefmt="%H:%M:%S",
)
logger = logging.getLogger("download_pubmed_activity")

EMAIL = "1007037391@qq.com"  # NCBI 要求提供联系方式
BASE = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"

# 主题优先级：数值越小优先级越高（合并后超量时先保留优先级高的）。
# 每个主题的查询词围绕"细胞系 + 调控元件 + 活性 + 影响原因/机制"。
QUERIES = [
    # 1. 通用：增强子/启动子活性决定因素与机制（核心背景）
    (0, "enhancer_activity_determinants",
     '("enhancer" OR "promoter" OR "cis-regulatory element" OR "regulatory element") '
     'AND ("activity" OR "expression") '
     'AND ("determinant" OR "required" OR "necessary" OR "mechanism" OR "underlie" '
     'OR "govern" OR "control" OR "modulate" OR "influence") '
     'AND ("transcription factor" OR "chromatin" OR "epigenetic" OR "DNA sequence") '
     'NOT review[Publication Type]'),
    # 2. K562 红系：决定增强子活性的红系因子
    (1, "K562_erythroid_determinants",
     '("K562") AND ("erythroid" OR "globin" OR GATA1 OR GATA2 OR KLF1 OR NFE2 OR FLI1) '
     'AND ("enhancer" OR "promoter" OR "regulatory element" OR "cis-regulatory") '
     'AND ("activity" OR "expression" OR "transcription") '
     'AND ("determinant" OR "required" OR "necessary" OR "regulate" OR "control" '
     'OR "activate" OR "repress" OR "mechanism") '
     'NOT review[Publication Type]'),
    # 3. HepG2 肝：决定增强子活性的肝富集因子
    (2, "HepG2_liver_determinants",
     '("HepG2") AND ("liver" OR "hepatocyte" OR HNF4A OR HNF1A OR CEBPA OR CEBPB '
     'OR FOXA1 OR FOXA2 OR ONECUT1) '
     'AND ("enhancer" OR "promoter" OR "regulatory element" OR "cis-regulatory") '
     'AND ("activity" OR "expression" OR "transcription") '
     'AND ("determinant" OR "required" OR "necessary" OR "regulate" OR "control" '
     'OR "activate" OR "repress" OR "mechanism") '
     'NOT review[Publication Type]'),
    # 4. SK-N-SH 神经母细胞瘤：决定增强子活性的神经因子
    (3, "SKNSH_neural_determinants",
     '("SK-N-SH" OR "neuroblastoma") AND (NEUROG2 OR ASCL1 OR NEUROD1 OR NEUROD2 '
     'OR SOX10 OR PHOX2B OR "neuronal" OR "neural" OR "neurogenesis") '
     'AND ("enhancer" OR "promoter" OR "regulatory element" OR "cis-regulatory") '
     'AND ("activity" OR "expression" OR "transcription") '
     'AND ("determinant" OR "required" OR "necessary" OR "regulate" OR "control" '
     'OR "activate" OR "repress" OR "mechanism") '
     'NOT review[Publication Type]'),
    # 5. 染色质/表观遗传对增强子活性的影响（机制背景）
    (4, "chromatin_epigenetic",
     '("enhancer" OR "promoter" OR "regulatory element") '
     'AND ("activity" OR "expression") '
     'AND ("chromatin state" OR "histone modification" OR "H3K27ac" OR "H3K4me1" '
     'OR "nucleosome" OR "DNA methylation" OR "accessibility" OR "ATAC-seq") '
     'AND ("determinant" OR "required" OR "mechanism" OR "influence" OR "modulate") '
     'NOT review[Publication Type]'),
    # 6. 序列特征对增强子活性的影响（GC/形状/基序位置，MPRA 机制背景）
    (5, "sequence_features",
     '("enhancer" OR "promoter" OR "regulatory element" OR "cis-regulatory") '
     'AND ("activity" OR "expression") '
     'AND ("GC content" OR "DNA shape" OR "motif position" OR "spacing" OR '
     '"sequence composition" OR "nucleotide" OR "single nucleotide") '
     'AND ("determinant" OR "influence" OR "affect" OR "determine" OR "mechanism") '
     'NOT review[Publication Type]'),
]

FETCH_BATCH = 50
SLEEP_BETWEEN_REQUESTS = 0.5


def esearch(query: str, max_results: int, min_year: int | None) -> list[str]:
    """relevance 排序检索，返回 PMID 列表。"""
    term = query
    if min_year:
        term = f"({query}) AND {min_year}/01/01:3000/12/31[dp]"
    params = {
        "db": "pubmed", "term": term, "retmax": max_results,
        "sort": "relevance", "retmode": "json", "tool": "BioDesign-Agent",
        "email": EMAIL,
    }
    for attempt in range(3):
        try:
            resp = requests.get(f"{BASE}/esearch.fcgi", params=params, timeout=60)
            resp.raise_for_status()
            data = resp.json()
            return data.get("esearchresult", {}).get("idlist", [])
        except Exception as e:  # noqa: BLE001
            logger.warning("esearch 第 %d 次失败 (%s)，重试...", attempt + 1, e)
            time.sleep(2 * (attempt + 1))
    logger.error("esearch 最终失败: %s", query[:80])
    return []


def _elem_text(parent: ET.Element, tag: str) -> str:
    """取 parent 下第一个 tag 子元素的文本（可能为 None）。"""
    node = parent.find(tag)
    if node is None or node.text is None:
        return ""
    return node.text.strip()


def _abstract_text(abstract: ET.Element) -> str:
    """拼接 Abstract/AbstractText 各段的文本。"""
    parts = []
    for at in abstract.findall("AbstractText"):
        label = at.get("Label")
        text = (at.text or "").strip()
        if text:
            parts.append(f"{label}: {text}" if label else text)
    return " ".join(parts)


def fetch_batch(pmids: list[str]) -> list[dict]:
    """批量 efetch 摘要（XML），提取元数据。"""
    out: list[dict] = []
    for attempt in range(3):
        try:
            resp = requests.get(
                f"{BASE}/efetch.fcgi",
                params={
                    "db": "pubmed", "id": ",".join(pmids),
                    "rettype": "abstract", "retmode": "xml",
                    "tool": "BioDesign-Agent", "email": EMAIL,
                },
                timeout=120,
            )
            resp.raise_for_status()
            root = ET.fromstring(resp.text)
            for article in root.findall(".//PubmedArticle"):
                try:
                    cit = article.find("MedlineCitation")
                    if cit is None:
                        continue
                    pmid = _elem_text(cit, "PMID")
                    art = cit.find("Article")
                    if art is None:
                        continue

                    title = ""
                    t_node = art.find("ArticleTitle")
                    if t_node is not None:
                        title = "".join(t_node.itertext()).strip()

                    abstract = ""
                    abs_node = art.find("Abstract")
                    if abs_node is not None:
                        abstract = _abstract_text(abs_node)

                    journal = ""
                    j_node = art.find("Journal")
                    if j_node is not None:
                        journal = _elem_text(j_node, "Title")

                    year = ""
                    if j_node is not None:
                        issue = j_node.find("JournalIssue")
                        if issue is not None:
                            pub_date = issue.find("PubDate")
                            if pub_date is not None:
                                year = _elem_text(pub_date, "Year")
                                if not year:
                                    md = pub_date.find("MedlineDate")
                                    if md is not None and md.text:
                                        m = re.search(r"\d{4}", md.text)
                                        year = m.group(0) if m else ""

                    doi = ""
                    aid_list = art.find("ArticleIdList")
                    if aid_list is not None:
                        for aid in aid_list.findall("ArticleId"):
                            if aid.get("IdType") == "doi":
                                doi = (aid.text or "").strip()
                                break

                    if abstract:
                        out.append({
                            "pmid": pmid,
                            "title": title,
                            "abstract": abstract,
                            "journal": journal,
                            "year": year,
                            "doi": doi,
                        })
                except Exception:  # noqa: BLE001
                    continue
            return out
        except Exception as e:  # noqa: BLE001
            logger.warning("efetch 第 %d 次失败 (%s)，重试...", attempt + 1, e)
            time.sleep(2 * (attempt + 1))
    logger.error("efetch 最终失败: 共 %d 个 PMID", len(pmids))
    return []


def main() -> None:
    parser = argparse.ArgumentParser(description="下载三系调控元件活性影响原因相关 PubMed 摘要")
    parser.add_argument("--output", default="data/raw/abstracts_activity.jsonl")
    parser.add_argument("--max-results", type=int, default=60,
                        help="每个主题最多取多少条（relevance 排序）")
    parser.add_argument("--total", type=int, default=240, help="去重后目标总条数")
    parser.add_argument("--min-year", type=int, default=2014,
                        help="仅保留该年份之后的文献（默认 2014，可设 0 关闭）")
    args = parser.parse_args()

    # 阶段 1: 各主题检索 ID，记录 (priority, query_name)
    id_meta: dict[str, tuple[int, str]] = {}
    for priority, name, query in QUERIES:
        ids = esearch(query, args.max_results, args.min_year)
        logger.info("主题 [%s] 命中 %d 条", name, len(ids))
        for pmid in ids:
            if pmid not in id_meta:
                id_meta[pmid] = (priority, name)
        time.sleep(SLEEP_BETWEEN_REQUESTS)

    all_pmids = list(id_meta.keys())
    logger.info("合并去重前共 %d 个 PMID，开始下载摘要...", len(all_pmids))

    # 阶段 2: 分批下载
    records: dict[str, dict] = {}
    for i in range(0, len(all_pmids), FETCH_BATCH):
        batch = all_pmids[i:i + FETCH_BATCH]
        items = fetch_batch(batch)
        for it in items:
            records[str(it["pmid"])] = it
        logger.info("下载进度 %d/%d (有效 %d)", i + len(batch), len(all_pmids), len(records))
        time.sleep(SLEEP_BETWEEN_REQUESTS)

    # 阶段 3: 配额制截断——每个主题保底 per_topic 条，剩余名额按优先级填充
    per_topic = max(1, args.total // len(QUERIES))
    by_topic: dict[str, list[tuple[str, dict]]] = {name: [] for _, name, _ in QUERIES}
    for pmid, it in records.items():
        by_topic[id_meta[pmid][1]].append((pmid, it))
    # 主题内部按年份倒序
    for name in by_topic:
        by_topic[name].sort(key=lambda kv: -int(kv[1].get("year") or 0))

    ranked: list[tuple[str, dict]] = []
    # 第一轮：每主题保底 per_topic 条
    for _, name, _ in QUERIES:
        ranked.extend(by_topic[name][:per_topic])
    # 第二轮：剩余名额按主题优先级继续填充
    remaining = args.total - len(ranked)
    if remaining > 0:
        for _, name, _ in QUERIES:
            if remaining <= 0:
                break
            extra = by_topic[name][per_topic:per_topic + remaining]
            ranked.extend(extra)
            remaining -= len(extra)
    ranked = ranked[: args.total]

    # 阶段 4: 写文件（带主题标注）
    out = Path(args.output)
    out.parent.mkdir(parents=True, exist_ok=True)
    with open(out, "w", encoding="utf-8") as f:
        for pmid, it in ranked:
            priority, query_name = id_meta[pmid]
            it["topic"] = query_name
            it["priority"] = priority
            f.write(json.dumps(it, ensure_ascii=False) + "\n")

    # 阶段 5: 统计
    from collections import Counter
    n_by_topic = Counter(it["topic"] for _, it in ranked)
    n_by_year = Counter(str(it.get("year") or "?") for _, it in ranked)
    logger.info("=" * 70)
    logger.info("完成: 共下载 %d 条摘要 -> %s", len(ranked), out)
    logger.info("主题分布:")
    for name, n in sorted(n_by_topic.items(), key=lambda x: -x[1]):
        logger.info("  %-28s %d", name, n)
    logger.info("年份分布(前5): %s",
                ", ".join(f"{y}:{n}" for y, n in n_by_year.most_common(5)))
    logger.info("=" * 70)


if __name__ == "__main__":
    main()
