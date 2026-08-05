import json
import time
from Bio import Entrez

Entrez.email = "1007037391@qq.com"  # 必填，换成你的邮箱
QUERY = "(transcription factor binding site) OR (regulatory DNA motif) AND human"
MAX_RESULTS = 500

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
    ids = fetch_ids(QUERY, MAX_RESULTS)
    print(f"Found {len(ids)} IDs. Downloading abstracts...")
    abstracts = fetch_abstracts(ids)
    with open("data/raw/abstracts.jsonl", "w") as f:
        for item in abstracts:
            f.write(json.dumps(item) + "\n")
    print(f"Saved {len(abstracts)} abstracts to data/raw/abstracts.jsonl")