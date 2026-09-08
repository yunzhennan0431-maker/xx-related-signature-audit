"""
扩充主题词表:第3.2节里指出继续深挖同一批候选池边际收益已很低,更有效的路径是扩大主题词表
本身。原23个种子主题(见01号脚本)已覆盖铁死亡/焦亡/铜死亡/二硫死亡/坏死性凋亡/自噬/乏氧/
糖酵解/EMT/血管生成/免疫/炎症/细胞衰老/N6-甲基腺苷/线粒体/脂代谢/胆固醇代谢/氧化应激/
谷氨酰胺代谢/凋亡/外泌体/肿瘤干细胞。本脚本新增一批近几年肿瘤转录组学"XX-related signature"
文献里同样高频但尚未纳入检索的过程/通路关键词,只查真实命中数(esearch的count字段不受retmax
影响),用于评估是否值得投入后续下载与人工复核。

合规访问:仅用NCBI官方E-utilities(esearch),遵守默认限速。
"""
import json
import time
from pathlib import Path

import requests

OUT_DIR = Path(__file__).resolve().parent.parent / "data"
OUT_DIR.mkdir(parents=True, exist_ok=True)

EUTILS_BASE = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi"
TOOL_NAME = "SignatureDL_audit_pilot"

# 新增候选主题:细胞死亡新模式(近年trendy)、RNA修饰(m6A之外)、翻译后修饰、代谢子类、
# 应激反应通路、微环境细胞互作相关关键词。避免和原23个种子重复。
NEW_THEMES = [
    "anoikis", "PANoptosis", "lactylation", "cuproptosis-related lncRNA",
    "N7-methylguanosine", "N1-methyladenosine", "5-methylcytosine",
    "mitophagy", "ferritinophagy", "cGAS-STING", "ubiquitination",
    "endoplasmic reticulum stress", "immunogenic cell death",
    "one-carbon metabolism", "amino acid metabolism", "purine metabolism",
    "neutrophil extracellular trap", "telomere", "SUMOylation",
    "alternative splicing", "cancer-associated fibroblast", "circular RNA",
    "chromatin regulator", "histone modification", "DNA damage repair",
]

QUERY_TEMPLATE = (
    '("{theme}-related" OR "{theme} related" OR "{theme}-associated" OR "{theme} associated") '
    'AND ("prognostic signature" OR "risk signature" OR "prognostic model" OR "risk model") '
    'AND (Cox OR LASSO)'
)


def esearch(term: str, retmax: int = 20) -> dict:
    params = {
        "db": "pubmed",
        "term": term,
        "retmode": "json",
        "retmax": retmax,
        "tool": TOOL_NAME,
    }
    resp = requests.get(EUTILS_BASE, params=params, timeout=30)
    resp.raise_for_status()
    return resp.json()


def main():
    results = []
    for theme in NEW_THEMES:
        query = QUERY_TEMPLATE.format(theme=theme)
        try:
            data = esearch(query)
            count = int(data["esearchresult"]["count"])
            pmids = data["esearchresult"]["idlist"]
        except Exception as e:
            print(f"[ERROR] {theme}: {e}")
            count, pmids = 0, []
        print(f"{theme}: {count} hits")
        results.append({"theme": theme, "count": count, "sample_pmids": pmids})
        time.sleep(0.4)

    with open(OUT_DIR / "new_theme_search_results.json", "w") as f:
        json.dump(results, f, ensure_ascii=False, indent=2)

    total = sum(r["count"] for r in results)
    print(f"\n新增{len(NEW_THEMES)}个主题合计命中(未去重): {total}")
    print(f"写入 {OUT_DIR}/new_theme_search_results.json")


if __name__ == "__main__":
    main()
