"""从 PubMed 定向下载与 MPRA 数据集（Gosai 2024）相关的高质量论文摘要。

覆盖主题（与 200nt 调控序列表达分析直接相关）：
  1. MPRA / 高通量报告实验中的调控序列与转录因子基序
  2. DNA 形状 (DNA shape) 对转录因子结合与调控活性的影响
  3. 基序出现位置 (motif position) 对增强子/启动子活性的影响
  4. 细胞类型特异性调控：K562 (红系) / HepG2 (肝) / SK-N-SH (神经)
  5. 增强子序列特征 (GC/结构/架构) 与表达的关系

每个主题按 relevance 排序取前 per_query 条，合并按 pmid 去重后，
按主题优先级排序截断到 --total 条（默认 300）。

用法（在项目根目录运行）:
    python scripts/download/download_pubmed_mpra.py --output data/raw/abstracts_mpra.jsonl
    python scripts/download/download_pubmed_mpra.py --max-results 80 --total 300 --min-year 2013
"""

import argparse
import json
import logging
import re
import sys
import time
from pathlib import Path

from Bio import Entrez

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s | %(levelname)s | %(message)s",
    datefmt="%H:%M:%S",
)
logger = logging.getLogger("download_pubmed_mpra")

Entrez.email = "1007037391@qq.com"  # 必填

# 主题优先级：数值越小优先级越高（合并后超量时先保留优先级高的）
QUERIES = [
    # 1. MPRA 高通量报告实验：调控序列 × 转录因子基序（核心）
    (0, "MPRA_enhancer_motif",
     '("MPRA" OR "massively parallel reporter assay") '
     'AND ("enhancer" OR "promoter" OR "regulatory sequence" OR "cis-regulatory") '
     'AND ("transcription factor" OR motif OR "binding site") '
     'NOT review[Publication Type]'),
    # 2. DNA 形状对转录因子结合与调控活性的影响
    (1, "DNA_shape",
     '("DNA shape" OR "DNAshape" OR "minor groove width" OR "helix twist" '
     'OR "DNA structural") '
     'AND ("transcription factor binding" OR "TF binding" OR enhancer OR regulatory '
     'OR "binding affinity") '
     'NOT review[Publication Type]'),
    # 3. 基序/结合位点位置对增强子活性的影响（位置效应）
    (2, "motif_position",
     '("motif position" OR "binding site position" OR "position-dependent" '
     'OR "position effect" OR "spacing") '
     'AND (enhancer OR promoter OR "regulatory") '
     'AND (activity OR expression OR "transcription") '
     'NOT review[Publication Type]'),
    # 4. MPRA 细胞类型特异性 / 序列-表达定量
    (3, "MPRA_cellspecific",
     '("MPRA" OR "massively parallel reporter assay" OR "STARR-seq" OR "lentiMPRA") '
     'AND ("cell type-specific" OR "cell-type-specific" OR "cell line" OR specific) '
     'AND (expression OR activity OR determinant) '
     'NOT review[Publication Type]'),
    # 5. K562 红系细胞：红系因子（GATA/KLF/NFE2）与增强子
    (4, "K562_erythroid",
     '("K562") AND (GATA1 OR GATA2 OR KLF1 OR "NFE2" OR "erythroid" OR "globin") '
     'AND (enhancer OR promoter OR regulatory OR "gene expression") '
     'NOT review[Publication Type]'),
    # 6. HepG2 肝细胞：肝富集因子（HNF/CEBP/FOXA）与增强子
    (5, "HepG2_liver",
     '("HepG2") AND (HNF4A OR HNF1A OR CEBPA OR CEBPB OR FOXA1 OR FOXA2 '
     'OR "hepatocyte" OR "liver") '
     'AND (enhancer OR promoter OR regulatory OR "gene expression") '
     'NOT review[Publication Type]'),
    # 7. SK-N-SH / 神经母细胞瘤：神经因子（NEUROG/ASCL/NEUROD）与增强子
    (6, "SKNSH_neural",
     '("SK-N-SH" OR "neuroblastoma") AND (NEUROG2 OR ASCL1 OR NEUROD1 OR NEUROD2 '
     'OR SOX10 OR "neuronal" OR "neural") '
     'AND (enhancer OR promoter OR regulatory OR "gene expression") '
     'NOT review[Publication Type]'),
    # 8. 增强子序列特征（GC 含量/架构/组合）与表达
    (7, "enhancer_determinants",
     '("enhancer activity" OR "enhancer function" OR "cis-regulatory") '
     'AND ("sequence features" OR "GC content" OR architecture OR grammar OR "sequence determinants") '
     'AND (expression OR activity OR prediction) '
     'NOT review[Publication Type]'),
]

FETCH_BATCH = 50
SLEEP_BETWEEN_REQUESTS = 0.5


def esearch(query: str, max_results: int, min_year: int | None) -> list[str]:
    """relevance 排序检索，返回 PMID 列表。"""
    term = query
    if min_year:
        term = f"({query}) AND {min_year}/01/01:3000/12/31[dp]"
    for attempt in range(3):
        try:
            handle = Entrez.esearch(
                db="pubmed", term=term, retmax=max_results, sort="relevance"
            )
            record = Entrez.read(handle)
            handle.close()
            return record["IdList"]
        except Exception as e:  # noqa: BLE001
            logger.warning("esearch 第 %d 次失败 (%s)，重试...", attempt + 1, e)
            time.sleep(2 * (attempt + 1))
    logger.error("esearch 最终失败: %s", query[:80])
    return []


def fetch_batch(pmids: list[str]) -> list[dict]:
    """批量 efetch 摘要，提取元数据。"""
    out: list[dict] = []
    for attempt in range(3):
        try:
            handle = Entrez.efetch(
                db="pubmed", id=pmids, rettype="abstract", retmode="xml"
            )
            records = Entrez.read(handle)
            handle.close()
            for article in records["PubmedArticle"]:
                try:
                    cit = article["MedlineCitation"]
                    pmid = cit["PMID"]
                    art = cit["Article"]
                    title = art.get("ArticleTitle", "")

                    abs_parts = art.get("Abstract", {}).get("AbstractText", [])
                    abstract = " ".join(str(t) for t in abs_parts) if abs_parts else ""

                    journal = ""
                    if "Journal" in art:
                        journal = art["Journal"].get("Title", "")

                    year = ""
                    pub_date = art.get("Journal", {}).get("JournalIssue", {}).get(
                        "PubDate", {})
                    year = str(pub_date.get("Year", ""))
                    if not year and "MedlineDate" in pub_date:
                        m = re.search(r"\d{4}", str(pub_date["MedlineDate"]))
                        year = m.group(0) if m else ""

                    doi = ""
                    for aid in art.get("ArticleIdList", []):
                        if aid.attributes.get("IdType") == "doi":
                            doi = str(aid)
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
    parser = argparse.ArgumentParser(description="下载 MPRA 相关定向 PubMed 摘要")
    parser.add_argument("--output", default="data/raw/abstracts_mpra.jsonl")
    parser.add_argument("--max-results", type=int, default=80,
                        help="每个主题最多取多少条（relevance 排序）")
    parser.add_argument("--total", type=int, default=300, help="去重后目标总条数")
    parser.add_argument("--min-year", type=int, default=2013,
                        help="仅保留该年份之后的文献（默认 2013，可设 0 关闭）")
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

    # 阶段 3: 按主题优先级排序并截断
    ranked = sorted(
        records.items(),
        key=lambda kv: (id_meta[kv[0]][0], -int(kv[1].get("year") or 0)),
    )
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
        logger.info("  %-24s %d", name, n)
    logger.info("年份分布(前5): %s",
                ", ".join(f"{y}:{n}" for y, n in n_by_year.most_common(5)))
    logger.info("=" * 70)


if __name__ == "__main__":
    main()
