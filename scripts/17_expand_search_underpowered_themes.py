"""
阶段1重启:针对性扩大检索,不是加新主题,而是给已经命中大量论文、但阶段1第一轮
只取了前30篇的9个主题,把retmax提高、把之前没拉到的论文补进候选池。

关键发现(2026-09-03):阶段1第一轮`01_pubmed_theme_search.py`跑23个主题时,
每个主题固定只取retmax=30,但很多主题实际命中远不止30篇——cuproptosis 349条、
autophagy/pyroptosis各252条、disulfidptosis 153条、hypoxia 131条、glycolysis 94条,
之前"候选池已榨干"这个结论只是"这30篇的信息含量榨干了",不是PubMed里真的没有更多候选。

这一步:对这9个主题重新检索(retmax提到100),和已有candidate_pool.csv去重,
只保留新PMID,过一遍idconv+OA桶可得性筛查+下载(复用04脚本逻辑),
再用05/06的候选表格扫描器打分,为下一步人工核实做准备。
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
TOOL = "SignatureDL_audit_expand"
SKIP_EXTENSIONS = {".pdf", ".jpg", ".jpeg", ".png", ".tif", ".tiff"}

EXPAND_THEMES = [
    "cuproptosis", "autophagy", "pyroptosis", "disulfidptosis", "hypoxia",
    "glycolysis", "cellular senescence", "mitochondrial", "exosome",
]
NEW_RETMAX = 100
QUERY_TEMPLATE = (
    '("{theme}-related" OR "{theme} related" OR "{theme}-associated" OR "{theme} associated") '
    'AND ("prognostic signature" OR "risk signature" OR "prognostic model" OR "risk model") '
    'AND (Cox OR LASSO)'
)


def get_with_retry(url, max_retries=6, **kwargs):
    # 这台机器的代理偶尔会"连接建立"但卡住不返回数据,requests的timeout有时候盖不住这种
    # 半死不活的隧道连接(踩过坑:进程CPU时间几乎是0但挂着不退出)。缩短单次timeout,
    # 配合更多重试次数,让它宁可多失败几次快速重试,也不要卡在一次连接上不动。
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
    """S3桶相关请求踩过坑:requests在这台机器上偶尔会"连接建立"但卡住不返回数据,
    进程CPU时间几乎是0但挂着不退出,timeout参数盖不住这种半死不活的隧道连接
    (重现过2次,curl对同一个URL反而立刻正常返回)。改成shell out到curl,
    用--max-time做进程级硬超时,卡住了curl子进程会被subprocess.run自己的timeout强制杀掉,
    不会一直挂着。"""
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
    with open(DATA_DIR / "candidate_pool.csv", encoding="utf-8") as f:
        for row in csv.DictReader(f):
            existing_pmids.add(row["pmid"])
    print(f"已有候选池: {len(existing_pmids)}")

    new_pool = []  # (pmid, theme)
    for theme in EXPAND_THEMES:
        query = QUERY_TEMPLATE.format(theme=theme)
        pmids = esearch(query, retmax=NEW_RETMAX)
        new_pmids = [p for p in pmids if p not in existing_pmids]
        print(f"  {theme:25s} 命中{len(pmids)}, 新增{len(new_pmids)}")
        for p in new_pmids:
            new_pool.append((p, theme))
        time.sleep(0.4)

    # 去重(同一pmid可能出现在多个主题里,先到先得)
    seen = set()
    dedup_pool = []
    for p, t in new_pool:
        if p in seen:
            continue
        seen.add(p)
        dedup_pool.append((p, t))
    print(f"\n新候选池(去重后): {len(dedup_pool)}")

    with open(DATA_DIR / "candidate_pool_expand.csv", "w", newline="", encoding="utf-8") as f:
        w = csv.writer(f)
        w.writerow(["pmid", "theme"])
        for p, t in dedup_pool:
            w.writerow([p, t])

    # idconv + OA桶 + 下载
    out_path = DATA_DIR / "pipeline_conversion_expand.csv"
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
