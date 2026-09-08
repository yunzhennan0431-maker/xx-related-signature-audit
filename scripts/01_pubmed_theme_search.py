"""
阶段1第一步:对预设的主题词表逐个查询PubMed候选论文池规模。
不下载全文,只统计每个主题的命中数量+抽取PMID列表,用于评估跨主题审计的实际规模,
并为下一步(全文/补充材料可得性抽样核查)提供样本来源。

合规访问:仅用NCBI官方E-utilities(esearch),遵守默认限速(无API key时<=3请求/秒)。
"""
import json
import time
from pathlib import Path

import requests

OUT_DIR = Path(__file__).resolve().parent.parent / "data"
OUT_DIR.mkdir(parents=True, exist_ok=True)

EUTILS_BASE = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi"
TOOL_NAME = "SignatureDL_audit_pilot"

# 主题种子词表:肿瘤转录组学"XX-related prognostic signature"文献里高频出现的过程/通路名词。
# 来源:研究者领域经验(常见于该类"建模模板"论文标题/关键词) + MSigDB过程类基因集命名习惯。
# 这是初版种子词表,阶段1执行中可能根据检索结果补充遗漏主题或剔除命中过低的噪声词。
THEMES = [
    "ferroptosis", "cuproptosis", "pyroptosis", "disulfidptosis", "necroptosis",
    "autophagy", "hypoxia", "glycolysis", "epithelial-mesenchymal transition",
    "angiogenesis", "immune", "inflammation", "cellular senescence",
    "N6-methyladenosine", "mitochondrial", "lipid metabolism",
    "cholesterol metabolism", "oxidative stress", "glutamine metabolism",
    "apoptosis", "exosome", "cancer stem cell", "cuproptosis-related ferroptosis",
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
    for theme in THEMES:
        query = QUERY_TEMPLATE.format(theme=theme)
        try:
            data = esearch(query)
            count = int(data["esearchresult"]["count"])
            pmids = data["esearchresult"]["idlist"]
        except Exception as e:
            count = None
            pmids = []
            print(f"[ERROR] theme={theme}: {e}")
        results.append({"theme": theme, "query": query, "count": count, "sample_pmids": pmids})
        print(f"{theme:35s} count={count}")
        time.sleep(0.4)  # 无API key,限速<=3请求/秒,留出余量

    out_path = OUT_DIR / "pubmed_theme_search_results.json"
    with open(out_path, "w", encoding="utf-8") as f:
        json.dump(results, f, ensure_ascii=False, indent=2)
    print(f"\nSaved to {out_path}")

    total = sum(r["count"] for r in results if r["count"] is not None)
    print(f"\n主题数: {len(results)}, 命中总数(未去重,主题间可能重叠): {total}")


if __name__ == "__main__":
    main()
