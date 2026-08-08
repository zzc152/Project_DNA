"""本地下载 JASPAR2024 脊椎动物 CORE 词表（用于 git 入库）。

与远程 data/vocab/jaspar2024_vertebrates_core.json 同源同法：
    https://jaspar.elixir.no/api/v1/matrix/?format=json&tax_group=vertebrates&collection=CORE&page_size=100&page={N}
"""
import json
import urllib.request

BASE = "https://jaspar.elixir.no/api/v1/matrix/?format=json&tax_group=vertebrates&collection=CORE&page_size=100&page={}"
HEADERS = {"User-Agent": "Mozilla/5.0"}
OUT = "D:/code_home/project/data/vocab/jaspar2024_vertebrates_core.json"

all_motifs = {}
page = 1
while True:
    req = urllib.request.Request(BASE.format(page), headers=HEADERS)
    data = json.load(urllib.request.urlopen(req, timeout=60))
    results = data.get("results", [])
    if not results:
        break
    for m in results:
        all_motifs[m["matrix_id"]] = {
            "name": m.get("name"),
            "base_id": m.get("base_id"),
            "version": m.get("version"),
            "collection": m.get("collection"),
        }
    print(f"page {page}: {len(results)} (累计 {len(all_motifs)})", flush=True)
    page += 1

payload = {"source": BASE, "api": "jaspar2024", "count": len(all_motifs), "motifs": all_motifs}
with open(OUT, "w", encoding="utf-8") as f:
    json.dump(payload, f, ensure_ascii=False, indent=1)
print(f"完成: {len(all_motifs)} motifs -> {OUT}")
