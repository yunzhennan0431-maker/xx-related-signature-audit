"""
05的补丁:很多期刊(Elsevier"Multimedia component"、部分"Supplementary Information"等)
在XML里完全不给补充材料的描述性标题,只有"Click here for additional data file"这类占位符,
单靠caption关键词打分覆盖不到。这一步改成直接打开候选Excel/CSV文件本身,读表头列名做打分,
表头信息(如"coefficient"/"HR"/"gene"/"symbol")比publisher给不给caption更可靠。

同时把没被05当作候选(caption没提到)但实际就是唯一/主要表格文件的情况也捞回来——
如果一个文章包补充材料里只有1-2个表格文件,不管caption写了什么,都值得直接打开看一眼。
"""
import json
import re
import zipfile
from pathlib import Path

try:
    import openpyxl
except ImportError:
    openpyxl = None

PKG_DIR = Path(__file__).resolve().parent.parent / "data" / "pmc_packages"
IN_PATH = Path(__file__).resolve().parent.parent / "data" / "extraction_candidates.json"
OUT_PATH = Path(__file__).resolve().parent.parent / "data" / "extraction_candidates_v2.json"

HEADER_POS = [
    (r"coefficient", 6), (r"\bcoef\b", 5), (r"\bhr\b", 3), (r"\bhazard", 4),
    (r"beta", 3), (r"\bgene ?symbol\b", 4), (r"^symbol$", 4), (r"\bgene\b", 2),
    (r"p ?value", 1), (r"risk ?score", 3),
]

# 和05脚本里的clinical_covariate_penalty同一个模式:"risk score作为临床协变量"表是
# 目前批量抽取里最大宗的假阳性来源,xlsx来源的候选也要查(top_rows只有前3行,能抓到的有限,
# 但聊胜于无,真正可靠的判断还是靠人工/LLM读)。
CLINICAL_VAR_TOKENS = re.compile(
    r"^(age|gender|sex|stage|grade|smoking|metastasis|race|tumor ?size|t ?stage|n ?stage|m ?stage|"
    r"differentiation|bmi|histolog\w*|clinical ?stage|residual tumor|venous invasion|radiation|"
    r"chemotherapy|alcohol)\b",
    re.I,
)


def peek_headers(path: Path):
    """返回(sheet名, 表头列表, 数据行数)的列表,读取失败返回None。"""
    if openpyxl is None or path.suffix.lower() not in (".xlsx", ".xls"):
        return None
    try:
        wb = openpyxl.load_workbook(path, read_only=True, data_only=True)
        out = []
        for ws in wb.worksheets:
            # 很多期刊的补充表第1行是纯文字标题(如"Supplementary Table S1. xxx"),
            # 真正的列名在第2-3行,所以前3行都读,分别保留,打分时全部纳入。
            top_rows = []
            for row in ws.iter_rows(min_row=1, max_row=3, values_only=True):
                cells = [str(c) for c in row if c is not None]
                if cells:
                    top_rows.append(cells)
            header = top_rows[0] if top_rows else []
            out.append({
                "sheet": ws.title, "header": header, "top_rows": top_rows, "max_row": ws.max_row,
            })
        return out
    except Exception:
        return None


def score_headers(sheets) -> int:
    if not sheets:
        return 0
    best = 0
    for s in sheets:
        text = " ".join(" ".join(row) for row in s.get("top_rows", [s["header"]])).lower()
        score = sum(w for pat, w in HEADER_POS if re.search(pat, text))
        clinical_hits = sum(
            1 for row in s.get("top_rows", []) if row and CLINICAL_VAR_TOKENS.match(row[0].strip())
        )
        if clinical_hits >= 2:
            score -= 8
        # 行数太少(<5)大概率不是完整signature表,行数太多(>500)大概率是全基因组DEG表而非最终模型
        if s["max_row"] and s["max_row"] < 5:
            score -= 2
        best = max(best, score)
    return best


def main():
    records = json.loads(IN_PATH.read_text(encoding="utf-8"))
    for rec in records:
        if "candidates" not in rec:
            continue
        pkg_dir = PKG_DIR / rec["pmcid"]
        for c in rec["candidates"]:
            if c.get("source") == "xml_table_wrap":
                # 05已经用caption+实际行内容打过总分,不需要再打开文件(也没有文件)
                continue
            fpath = pkg_dir / c["file"]
            if not fpath.exists():
                c["headers"] = None
                continue
            sheets = peek_headers(fpath)
            c["headers"] = sheets
            content_score = score_headers(sheets)
            c["content_score"] = content_score
            c["total_score"] = c["score"] + content_score
        rec["candidates"].sort(key=lambda c: -c.get("total_score", c["score"]))

    OUT_PATH.write_text(json.dumps(records, ensure_ascii=False, indent=2), encoding="utf-8")

    n_strong = sum(
        1 for r in records
        if r.get("candidates") and r["candidates"][0].get("total_score", 0) >= 5
    )
    print(f"总文章包: {len(records)}")
    print(f"最高分候选文件total_score>=5(强候选): {n_strong}")
    print(f"写入 {OUT_PATH}")


if __name__ == "__main__":
    main()
