"""
通用批量dump工具:给定一批PMCID,把它们排名第1(或指定)候选的完整内容dump成文本,
供人工/LLM逐篇判读+抽取,替代08脚本的纯规则列识别。
用法: python3 09_dump_batch.py PMC1 PMC2 PMC3 ...
"""
import json
import sys
from pathlib import Path

import openpyxl

DATA_DIR = Path(__file__).resolve().parent.parent / "data"
PKG_DIR = DATA_DIR / "pmc_packages"
OUT_DIR = DATA_DIR / "batch_dump"
OUT_DIR.mkdir(exist_ok=True)

candidates = {r["pmcid"]: r for r in json.loads((DATA_DIR / "extraction_candidates_v2.json").read_text(encoding="utf-8"))}


def dump_xlsx(path: Path, max_rows=60):
    lines = []
    try:
        wb = openpyxl.load_workbook(path, read_only=True, data_only=True)
    except Exception as e:
        return [f"[ERROR loading: {e}]"]
    for ws in wb.worksheets[:2]:
        lines.append(f"--- sheet: {ws.title} (max_row={ws.max_row}, max_col={ws.max_column}) ---")
        for i, row in enumerate(ws.iter_rows(min_row=1, max_row=max_rows, values_only=True)):
            cells = [str(c) if c is not None else "" for c in row]
            if any(c.strip() for c in cells):
                lines.append(" | ".join(cells))
            if i >= max_rows:
                break
        lines.append("")
    return lines


def main(pmcids):
    for pmcid in pmcids:
        rec = candidates.get(pmcid)
        if not rec or not rec.get("candidates"):
            print(f"{pmcid}: NO CANDIDATES")
            continue
        top = rec["candidates"][0]
        out_txt = OUT_DIR / f"{pmcid}.txt"
        header = f"=== {pmcid} :: source={top['source']} score={top.get('total_score', top.get('score'))} ===\ncaption: {top.get('caption','')}\n\n"
        if top["source"] == "xml_table_wrap":
            lines = [" | ".join(row) for row in (top.get("rows") or [])]
        else:
            fpath = PKG_DIR / pmcid / top["file"]
            if fpath.suffix.lower() == ".docx":
                lines = ["[docx,需另外转换]"]
            elif fpath.suffix.lower() == ".csv":
                try:
                    lines = fpath.read_text(encoding="utf-8", errors="ignore").splitlines()[:60]
                except Exception as e:
                    lines = [f"[ERROR: {e}]"]
            elif fpath.exists():
                lines = dump_xlsx(fpath)
            else:
                lines = [f"[FILE NOT FOUND: {fpath}]"]
        out_txt.write_text(header + "\n".join(lines), encoding="utf-8")
        print(f"{pmcid}: {len(lines)} lines -> {out_txt}")


if __name__ == "__main__":
    main(sys.argv[1:])
