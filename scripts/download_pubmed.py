import argparse
import json
import time
from Bio import Entrez

Entrez.email = "1007037391@qq.com"  # 必填，换成你的邮箱

# 定向查询：只检索包含"具体调控事件"的文章（有 TF→靶基因 调控动作），
# 按相关性排序取前 N 条，避免检索到方法学/资源/综述类文章。
# 查询词经实测调优：过严（多短语 AND）只命中 52 条；本版本命中约 5000 条，
# relevance 排序后取前 200，实体密度最高。
QUERY = (
    '("transcription factor" OR "TF") '
    'AND ("target gene" OR "activates transcription" OR "represses transcription") '
    'AND human'
)
MAX_RESULTS = 200

def fetch_ids(query, max_results):
    handle = Entrez.esearch(db="pubmed", term=query, retmax=max_results, sort="relevance")
    record = Entrez.read(handle)
    handle.close()
    return record["IdList"]
def fetch_abstracts(id_list, batch_size=50):
    abstracts = []
    for i in range(0, len(id_list), batch_size):
        batch = id_list[i:i+batch_size]
        handle = Entrez.efetch(db="pubmed", id=batch, rettype="abstract", retmode="xml")
        records = Entrez.read(handle)
        handle.close()
        for article in records["PubmedArticle"]:
            try:
                pmid = article["MedlineCitation"]["PMID"]
                title = article["MedlineCitation"]["Article"]["ArticleTitle"]
                abs_list = article["MedlineCitation"]["Article"].get("Abstract", {}).get("AbstractText", [])
                abstract = " ".join(abs_list) if abs_list else ""
                if abstract:
                    abstracts.append({"pmid": pmid, "title": title, "abstract": abstract})
            except:
                continue
        time.sleep(0.5)
    return abstracts

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="下载 PubMed 摘要（定向查询）")
    parser.add_argument("--output", default="data/raw/abstracts_targeted.jsonl",
                        help="输出文件（默认 data/raw/abstracts_targeted.jsonl，避免覆盖原始 500 条）")
    parser.add_argument("--max-results", type=int, default=MAX_RESULTS,
                        help="下载条数（默认 %d）" % MAX_RESULTS)
    args = parser.parse_args()

    ids = fetch_ids(QUERY, args.max_results)
    print(f"Found {len(ids)} IDs. Downloading abstracts...")
    abstracts = fetch_abstracts(ids)
    with open(args.output, "w") as f:
        for item in abstracts:
            f.write(json.dumps(item) + "\n")
    print(f"Saved {len(abstracts)} abstracts to {args.output}")