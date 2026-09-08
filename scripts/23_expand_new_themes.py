"""
在22号脚本摸底的基础上,正式把新主题词表纳入候选池:去重(和已有1043篇候选池比对)、
idconv转PMCID、OA桶可得性核查、下载完整文章包。复用17号脚本已验证稳定的curl子进程方案
(避免requests在这台机器上对S3桶请求偶发卡死的问题)。

排除"cuproptosis-related lncRNA"——22号脚本摸底显示这个query和已有的cuproptosis主题
候选池20/20样本完全重合,是同一批论文的子集查询,不是真正的新主题,纳入只会浪费idconv/
下载配额,不产生新增候选。
"""
import csv
import re
import subprocess
import time
from pathlib import Path
from typing import Dict, Optional

import requests

DATA_DIR = Path(__file__).resolve().parent.parent / "data"
PKG_DIR = DATA_DIR / "pmc_packages"
EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
IDCONV = "https://www.ncbi.nlm.nih.gov/pmc/utils/idconv/v1.0/"
BUCKET = "https://pmc-oa-opendata.s3.amazonaws.com/"
TOOL = "SignatureDL_audit_newthemes"
SKIP_EXTENSIONS = {".pdf", ".jpg", ".jpeg", ".png", ".tif", ".tiff"}

NEW_THEMES = [
    "anoikis", "PANoptosis", "lactylation", "N7-methylguanosine",
    "N1-methyladenosine", "5-methylcytosine", "mitophagy", "ferritinophagy",
    "cGAS-STING", "ubiquitination", "endoplasmic reticulum stress",
    "immunogenic cell death", "one-carbon metabolism", "amino acid metabolism",
    "purine metabolism", "neutrophil extracellular trap", "telomere",
    "SUMOylation", "alternative splicing", "cancer-associated fibroblast",
    "circular RNA", "chromatin regulator", "histone modification",
    "DNA damage repair",
]
RETMAX = 250
QUERY_TEMPLATE = (
    '("{theme}-related" OR "{theme} related" OR "{theme}-associated" OR "{theme} associated") '
    'AND ("prognostic signature" OR "risk signature" OR "prognostic model" OR "risk model") '
    'AND (Cox OR LASSO)'
)


def get_with_retry(url, max_retries=6, **kwargs):
    kwargs.setdefault("timeout", 10)
    last_exc = None
    for attempt in range(max_retries):
        try:
            r = requests.get(url, **kwargs)
            r.raise_for_status()
            return r
        except requests.exceptions.RequestException as e:
            last_exc = e
            wait = 2 ** attempt
            print(f"  [retry {attempt+1}/{max_retries}] {url} failed: {e}; sleep {wait}s")
            time.sleep(wait)
    raise last_exc


def esearch(term: str, retmax: int) -> list:
    r = get_with_retry(f"{EUTILS}/esearch.fcgi", params={
        "db": "pubmed", "term": term, "retmode": "json", "retmax": retmax, "tool": TOOL,
    })
    return r.json()["esearchresult"]["idlist"]


def batch_idconv(pmids: list) -> Dict[str, Optional[str]]:
    result = {}
    batch = 200
    for i in range(0, len(pmids), batch):
        chunk = pmids[i:i + batch]
        r = get_with_retry(IDCONV, params={"ids": ",".join(chunk), "format": "json", "tool": TOOL}, timeout=60)
        data = r.json()
        for rec in data.get("records", []):
            result[str(rec["requested-id"])] = rec.get("pmcid")
        time.sleep(0.4)
    return result


def curl_get(url: str, max_retries=4, timeout=15) -> Optional[bytes]:
    for attempt in range(max_retries):
        try:
            result = subprocess.run(
                ["curl", "-sS", "--max-time", str(timeout), "--retry", "0", url],
                capture_output=True, timeout=timeout + 5,
            )
            if result.returncode == 0:
                return result.stdout
            print(f"  [curl retry {attempt+1}/{max_retries}] {url} curl exit={result.returncode}: "
                  f"{result.stderr.decode(errors='replace')[:200]}")
        except subprocess.TimeoutExpired:
            print(f"  [curl retry {attempt+1}/{max_retries}] {url} subprocess超时,强制杀掉重试")
        time.sleep(2 ** attempt)
    return None


def bucket_files(pmcid_num: str) -> list:
    body = curl_get(f"{BUCKET}?list-type=2&prefix=PMC{pmcid_num}", timeout=15)
    if body is None:
        return []
    return re.findall(r"<Key>(.*?)</Key>", body.decode(errors="replace"))


def download_package(keys: list, pmcid_num: str) -> None:
    out_dir = PKG_DIR / f"PMC{pmcid_num}"
    out_dir.mkdir(parents=True, exist_ok=True)
    for key in keys:
        fname = key.split("/")[-1]
        ext = Path(fname).suffix.lower()
        if ext in SKIP_EXTENSIONS:
            continue
        dest = out_dir / fname
        if dest.exists():
            continue
        content = curl_get(BUCKET + key, timeout=30)
        if content is not None:
            dest.write_bytes(content)


def main():
    existing_pmids = set()
    for fn in ["candidate_pool.csv", "candidate_pool_expand.csv"]:
        fp = DATA_DIR / fn
        if fp.exists():
            with open(fp, encoding="utf-8") as f:
                for row in csv.DictReader(f):
                    existing_pmids.add(row["pmid"])
    print(f"已有候选池(两轮合计): {len(existing_pmids)}")

    new_pool = []
    for theme in NEW_THEMES:
        query = QUERY_TEMPLATE.format(theme=theme)
        pmids = esearch(query, retmax=RETMAX)
        new_pmids = [p for p in pmids if p not in existing_pmids]
        print(f"  {theme:28s} 命中{len(pmids)}, 新增{len(new_pmids)}")
        for p in new_pmids:
            new_pool.append((p, theme))
        time.sleep(0.4)

    seen = set()
    dedup_pool = []
    for p, t in new_pool:
        if p in seen:
            continue
        seen.add(p)
        dedup_pool.append((p, t))
    print(f"\n新增主题候选池(去重后): {len(dedup_pool)}")

    with open(DATA_DIR / "candidate_pool_newthemes.csv", "w", newline="", encoding="utf-8") as f:
        w = csv.writer(f)
        w.writerow(["pmid", "theme"])
        for p, t in dedup_pool:
            w.writerow([p, t])

    out_path = DATA_DIR / "pipeline_conversion_newthemes.csv"
    done = set()
    if out_path.exists():
        with open(out_path, encoding="utf-8") as f:
            for row in csv.DictReader(f):
                done.add(row["pmid"])
        print(f"续跑: 已完成{len(done)}条, 跳过")

    pmids_only = [p for p, t in dedup_pool]
    print("\n批量idconv查询...")
    pmcid_map = batch_idconv(pmids_only)
    n_found = sum(1 for v in pmcid_map.values() if v)
    print(f"idconv结果: {n_found}/{len(pmids_only)} 有PMCID ({n_found/len(pmids_only):.1%})")

    mode = "a" if out_path.exists() else "w"
    with open(out_path, mode, newline="", encoding="utf-8") as f:
        w = csv.writer(f)
        if mode == "w":
            w.writerow(["pmid", "theme", "pmcid", "has_pmc_link", "in_oa_bucket", "n_files", "downloaded"])
        for i, (pmid, theme) in enumerate(dedup_pool):
            if pmid in done:
                continue
            pmcid_full = pmcid_map.get(pmid)
            pmcid_num = pmcid_full.replace("PMC", "") if pmcid_full else None
            in_bucket, n_files, downloaded = False, 0, False
            if pmcid_num:
                keys = bucket_files(pmcid_num)
                if keys:
                    in_bucket, n_files = True, len(keys)
                    download_package(keys, pmcid_num)
                    downloaded = True
                time.sleep(0.1)
            w.writerow([pmid, theme, pmcid_full or "", bool(pmcid_num), in_bucket, n_files, downloaded])
            f.flush()
            if (i + 1) % 50 == 0:
                print(f"  processed {i + 1}/{len(dedup_pool)}")

    print(f"\n写入 {out_path}")


if __name__ == "__main__":
    main()
