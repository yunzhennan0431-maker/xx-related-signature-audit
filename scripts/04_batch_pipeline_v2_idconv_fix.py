"""
修正版流水线:用官方ID Converter批量接口(idconv,一次最多200个id)替换逐篇`elink`查询。

踩坑记录(详见pilot_extraction_notes.md/日志.md):03版脚本用elink逐篇查询PMC链接,
在长时间批量运行中部分主题出现"疑似0%转化率"的异常结果,独立用idconv复核后证实是
elink在长时间高频请求下出现静默失败(返回200但结果为空,不报错),而非真实的PMC覆盖率差异。
idconv是批量接口,请求数从617次降到4次,既更快也避免了这个问题。

同时保留S3桶可得性检查与包下载逻辑,复用02/03脚本已验证过的部分。
"""
import csv
import re
import time
from pathlib import Path
from typing import Dict, Optional

import requests

DATA_DIR = Path(__file__).resolve().parent.parent / "data"
PKG_DIR = DATA_DIR / "pmc_packages"
IDCONV = "https://www.ncbi.nlm.nih.gov/pmc/utils/idconv/v1.0/"
BUCKET = "https://pmc-oa-opendata.s3.amazonaws.com/"
TOOL = "SignatureDL_audit_pilot"


def get_with_retry(url, max_retries=4, **kwargs):
    """带指数退避重试的GET,踩过一次代理/连接超时把整个批量任务打断的坑之后加的。"""
    kwargs.setdefault("timeout", 30)
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


def batch_idconv(pmids: list) -> Dict[str, Optional[str]]:
    """批量PMID->PMCID映射,官方权威接口,一次最多200个id。"""
    result = {}
    batch = 200
    for i in range(0, len(pmids), batch):
        chunk = pmids[i:i + batch]
        r = get_with_retry(IDCONV, params={
            "ids": ",".join(chunk), "format": "json", "tool": TOOL,
        }, timeout=60)
        data = r.json()
        for rec in data.get("records", []):
            result[str(rec["requested-id"])] = rec.get("pmcid")
        time.sleep(0.4)
    return result


def bucket_files(pmcid_num: str) -> list:
    r = get_with_retry(BUCKET, params={"list-type": "2", "prefix": f"PMC{pmcid_num}"}, timeout=30)
    return re.findall(r"<Key>(.*?)</Key>", r.text)


# 结构化抽取用不到整篇PDF和图片(体积大、内容用不上),只拉正文XML/元数据/纯文本和
# 真正含数据的补充材料格式;PDF/图片大文件多次触发本机代理的连接不稳定,跳过后
# 既省带宽又减少失败面。
SKIP_EXTENSIONS = {".pdf", ".jpg", ".jpeg", ".png", ".tif", ".tiff", ".gif"}


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
        r = get_with_retry(BUCKET + key, timeout=60)
        dest.write_bytes(r.content)


def main():
    pool = []
    with open(DATA_DIR / "candidate_pool.csv", encoding="utf-8") as f:
        for row in csv.DictReader(f):
            pool.append((row["pmid"], row["theme"]))
    pmids = [p for p, _ in pool]

    print(f"批量idconv查询 {len(pmids)} 个PMID...")
    pmcid_map = batch_idconv(pmids)
    n_found = sum(1 for v in pmcid_map.values() if v)
    print(f"idconv结果: {n_found}/{len(pmids)} 有PMCID ({n_found/len(pmids):.1%})")

    out_path = DATA_DIR / "pipeline_conversion_v2.csv"
    done = set()
    if out_path.exists():
        with open(out_path, encoding="utf-8") as f:
            for row in csv.DictReader(f):
                done.add(row["pmid"])
        print(f"续跑:已完成 {len(done)} 条,跳过")

    mode = "a" if out_path.exists() else "w"
    with open(out_path, mode, newline="", encoding="utf-8") as f:
        w = csv.writer(f)
        if mode == "w":
            w.writerow(["pmid", "theme", "pmcid", "has_pmc_link", "in_oa_bucket", "n_files", "downloaded"])
        for i, (pmid, theme) in enumerate(pool):
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
            if (i + 1) % 100 == 0:
                print(f"  processed {i + 1}/{len(pool)}")

    print(f"\n修正后结果写入 {out_path}")


if __name__ == "__main__":
    main()
