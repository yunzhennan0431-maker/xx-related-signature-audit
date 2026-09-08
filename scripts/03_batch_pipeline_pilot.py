"""
阶段1批量试点:串起"主题检索→PMC全文链接筛选→OA批量分发子集(S3桶)可得性筛选→下载文章包"
整条机械化流水线,在一个中等规模样本(而非全量几千篇)上跑出真实的两层转化率统计,
并为下一步LLM结构化抽取准备好本地文章包。

三层输出,写入 analysis_output/data/:
  - candidate_pool.csv        主题检索去重后的候选PMID池
  - pipeline_conversion.csv   每个PMID的PMC链接/OA桶可得性核查结果(逐行落盘,支持中断续跑)
  - pmc_packages/PMCxxxx/     实际下载的文章包(仅"OA桶里能找到"的论文)

合规性:所有网络请求走NCBI E-utilities(esearch/elink,官方API)和PMC Cloud Service公开S3桶
(pmc-oa-opendata,匿名HTTPS,NLM 2026年官方文档确认的现行合规渠道),不抓取pmc.ncbi.nlm.nih.gov
人类可读页面(会被reCAPTCHA拦截,见pilot_extraction_notes.md)。
"""
import csv
import re
import time
from pathlib import Path
from typing import Optional

import requests

DATA_DIR = Path(__file__).resolve().parent.parent / "data"
PKG_DIR = DATA_DIR / "pmc_packages"
EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
BUCKET = "https://pmc-oa-opendata.s3.amazonaws.com/"
TOOL = "SignatureDL_audit_pilot"

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
PER_THEME_CAP = 30  # 每个主题最多纳入的候选PMID数,控制本轮试点规模


def esearch(term: str, retmax: int) -> list[str]:
    r = requests.get(f"{EUTILS}/esearch.fcgi", params={
        "db": "pubmed", "term": term, "retmode": "json", "retmax": retmax, "tool": TOOL,
    }, timeout=30)
    r.raise_for_status()
    return r.json()["esearchresult"]["idlist"]


def elink_pmc(pmid: str) -> Optional[str]:
    r = requests.get(f"{EUTILS}/elink.fcgi", params={
        "dbfrom": "pubmed", "db": "pmc", "id": pmid, "retmode": "json", "tool": TOOL,
    }, timeout=30)
    r.raise_for_status()
    ls = r.json()["linksets"][0]
    for ldb in ls.get("linksetdbs", []):
        if ldb.get("linkname") == "pubmed_pmc":
            return ldb["links"][0]
    return None


def bucket_files(pmcid_num: str) -> list[str]:
    r = requests.get(BUCKET, params={"list-type": "2", "prefix": f"PMC{pmcid_num}"}, timeout=30)
    r.raise_for_status()
    return re.findall(r"<Key>(.*?)</Key>", r.text)


def download_package(keys: list[str], pmcid_num: str) -> None:
    out_dir = PKG_DIR / f"PMC{pmcid_num}"
    out_dir.mkdir(parents=True, exist_ok=True)
    for key in keys:
        fname = key.split("/")[-1]
        dest = out_dir / fname
        if dest.exists():
            continue
        # 补充材料/表格优先,大PDF/大图跳过不下载(不是抽取需要的内容,且是超时主因)
        if fname.lower().endswith((".pdf", ".tif", ".tiff")):
            continue
        for attempt in range(3):
            try:
                r = requests.get(BUCKET + key, timeout=(10, 120))
                r.raise_for_status()
                dest.write_bytes(r.content)
                break
            except requests.exceptions.RequestException as e:
                if attempt == 2:
                    print(f"    [WARN] 下载失败,跳过 {key}: {e}")
                else:
                    time.sleep(2)


def step1_build_pool():
    pool = {}  # pmid -> theme (先到先得,记录首次命中的主题)
    for theme in THEMES:
        query = QUERY_TEMPLATE.format(theme=theme)
        pmids = esearch(query, retmax=PER_THEME_CAP)
        new = 0
        for pmid in pmids:
            if pmid not in pool:
                pool[pmid] = theme
                new += 1
        print(f"  {theme:35s} +{new:3d} new (query returned {len(pmids)})")
        time.sleep(0.4)

    path = DATA_DIR / "candidate_pool.csv"
    with open(path, "w", newline="", encoding="utf-8") as f:
        w = csv.writer(f)
        w.writerow(["pmid", "theme"])
        for pmid, theme in pool.items():
            w.writerow([pmid, theme])
    print(f"\n候选池: {len(pool)} 篇(去重后),写入 {path}")
    return pool


def step2_conversion_pipeline(pool: dict):
    path = DATA_DIR / "pipeline_conversion.csv"
    # 支持中断续跑:已处理过的pmid跳过
    done = set()
    if path.exists():
        with open(path, encoding="utf-8") as f:
            for row in csv.DictReader(f):
                done.add(row["pmid"])

    mode = "a" if path.exists() else "w"
    with open(path, mode, newline="", encoding="utf-8") as f:
        w = csv.writer(f)
        if mode == "w":
            w.writerow(["pmid", "theme", "pmcid", "has_pmc_link", "in_oa_bucket", "n_files", "downloaded"])
        for i, (pmid, theme) in enumerate(pool.items()):
            if pmid in done:
                continue
            try:
                pmcid = elink_pmc(pmid)
                time.sleep(0.35)
                in_bucket, n_files, downloaded = False, 0, False
                if pmcid:
                    keys = bucket_files(pmcid)
                    if keys:
                        in_bucket, n_files = True, len(keys)
                        download_package(keys, pmcid)
                        downloaded = True
                    time.sleep(0.1)
                w.writerow([pmid, theme, pmcid or "", bool(pmcid), in_bucket, n_files, downloaded])
            except requests.exceptions.RequestException as e:
                print(f"  [WARN] pmid={pmid} 请求失败,记为待重试并跳过: {e}")
                w.writerow([pmid, theme, "", "", "", "", ""])
            f.flush()
            if (i + 1) % 20 == 0:
                print(f"  processed {i + 1}/{len(pool)}")
    print(f"\n转化结果写入 {path}")


def main():
    print("=== Step 1: 构建候选池 ===")
    pool = step1_build_pool()
    print("\n=== Step 2: PMC链接 + OA桶可得性 + 下载 ===")
    step2_conversion_pipeline(pool)


if __name__ == "__main__":
    main()
