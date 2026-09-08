"""
通过PMC Cloud Service(公开S3桶pmc-oa-opendata,匿名HTTPS访问,无需AWS账号)
拉取指定PMCID的完整文章包(正文XML/PDF/元数据JSON/图片/补充材料原始文件)。

替代已于2026-08-24下线的旧FTP/oa.fcgi服务。桶内路径规则:
  https://pmc-oa-opendata.s3.amazonaws.com/PMC{id}.{version}/PMC{id}.{version}.{ext}
  补充材料保留原始文件名,如 .../PMC{id}.{version}/Table2.xlsx

用法示例:
  python3 02_pmc_s3_fetch.py PMC13481506          # 列出该文章在桶里的所有文件
  python3 02_pmc_s3_fetch.py PMC13481506 --download  # 下载全部文件到 data/pmc_packages/PMC13481506/
"""
import argparse
import re
import sys
import time
from pathlib import Path

import requests

BUCKET_BASE = "https://pmc-oa-opendata.s3.amazonaws.com/"
OUT_ROOT = Path(__file__).resolve().parent.parent / "data" / "pmc_packages"


def list_package_files(pmcid: str) -> list[str]:
    """返回该PMCID在桶里的完整key列表(含版本号前缀),找不到则返回空列表。"""
    resp = requests.get(BUCKET_BASE, params={"list-type": "2", "prefix": pmcid}, timeout=30)
    resp.raise_for_status()
    keys = re.findall(r"<Key>(.*?)</Key>", resp.text)
    return keys


def download_package(pmcid: str, keys: list[str]) -> None:
    out_dir = OUT_ROOT / pmcid
    out_dir.mkdir(parents=True, exist_ok=True)
    for key in keys:
        fname = key.split("/")[-1]
        dest = out_dir / fname
        if dest.exists():
            continue
        r = requests.get(BUCKET_BASE + key, timeout=60)
        r.raise_for_status()
        dest.write_bytes(r.content)
        print(f"  saved {fname} ({len(r.content)} bytes)")
        time.sleep(0.1)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("pmcid", help="e.g. PMC13481506 (不带版本号后缀)")
    ap.add_argument("--download", action="store_true", help="下载全部文件,默认只列出")
    args = ap.parse_args()

    keys = list_package_files(args.pmcid)
    if not keys:
        print(f"{args.pmcid}: 未在pmc-oa-opendata桶中找到(不属于OA批量分发子集,或PMCID有误)")
        sys.exit(1)

    print(f"{args.pmcid}: 找到 {len(keys)} 个文件")
    for k in keys:
        print(" ", k)

    if args.download:
        download_package(args.pmcid, keys)


if __name__ == "__main__":
    main()
