"""
阶段1第五步:对全部502篇文章包,取排名第1的候选表格,尝试自动识别"基因symbol列"和
"系数/HR列"并抽取成结构化行,产出第一版(未经人工核查)的抽取结果表。

这是机械化的列识别,不是语义判断"这个表是不是真的最终模型"——那一步已经在21篇试点里
论证过,打分第1名不保证一定对(pilot显示52.4%严格命中、71.4%至少能拿到基因列表)。
这里的产出必须标注为"未核查",后续要在全量结果里抽一批新样本(不能复用调过参的21篇)
做独立人工核查,才能报出可信的最终准确率。

输出: analysis_output/data/extracted_signatures_raw.csv (逐基因一行) +
      analysis_output/data/extraction_summary.csv (逐论文一行,含是否成功抽取)
"""
import csv
import re
from pathlib import Path

import openpyxl

DATA_DIR = Path(__file__).resolve().parent.parent / "data"
PKG_DIR = DATA_DIR / "pmc_packages"
CANDIDATES_PATH = DATA_DIR / "extraction_candidates_v2.json"
CONVERSION_PATH = DATA_DIR / "pipeline_conversion_v2.csv"
OUT_ROWS = DATA_DIR / "extracted_signatures_raw.csv"
OUT_SUMMARY = DATA_DIR / "extraction_summary.csv"

GENE_PATTERN = re.compile(r"^[A-Z][A-Z0-9\-\.]{2,14}$")  # 粗略gene symbol形状:大写开头,无空格,长度>=3,不太长
MIN_SCORE_TO_ATTEMPT = 0  # 负分候选(如已知是引物表/动物实验表)不再尝试抽取,直接跳过
GENE_HEADER_PATTERN = re.compile(r"gene ?symbol|^gene$|^id$|^symbol$", re.I)
COEF_HEADER_PATTERNS = [
    re.compile(r"coefficient|^coef$", re.I),
    re.compile(r"^beta$|weight", re.I),
    re.compile(r"^hr$|hazard ?ratio", re.I),
]


def load_json(path):
    import json
    return json.loads(path.read_text(encoding="utf-8"))


def load_full_rows(pmcid: str, candidate: dict, max_rows: int = 200):
    """拿到候选表格的完整行(不只是打分用的前3行)。table_wrap已经在05里存了rows(上限30/60行内);
    xlsx/csv需要重新打开文件读全部。"""
    if candidate["source"] == "xml_table_wrap":
        return candidate.get("rows") or []
    fpath = PKG_DIR / pmcid / candidate["file"]
    if not fpath.exists():
        return []
    try:
        if fpath.suffix.lower() == ".csv":
            import csv as csvmod
            with open(fpath, encoding="utf-8", errors="ignore") as f:
                return [row for row in csvmod.reader(f)][:max_rows]
        elif fpath.suffix.lower() in (".xlsx", ".xls"):
            wb = openpyxl.load_workbook(fpath, read_only=True, data_only=True)
            ws = wb.worksheets[0]
            rows = []
            for row in ws.iter_rows(min_row=1, max_row=max_rows, values_only=True):
                cells = [str(c) if c is not None else "" for c in row]
                if any(c.strip() for c in cells):
                    rows.append(cells)
            return rows
    except Exception:
        return []
    return []


def find_header_row(rows):
    """很多表第1行是纯标题句(如'Supplementary Table S1. xxx'),真正表头在第2-3行。
    找第一行'列数>=2且大多数格子是短词'的行当表头。"""
    for i, row in enumerate(rows[:5]):
        non_empty = [c for c in row if c.strip()]
        if len(non_empty) >= 2 and all(len(c) < 40 for c in non_empty):
            return i
    return 0


def extract_gene_coef(rows):
    if len(rows) < 2:
        return None, None, []
    header_idx = find_header_row(rows)
    header = rows[header_idx]
    data_rows = rows[header_idx + 1:]

    gene_col = None
    for i, h in enumerate(header):
        if GENE_HEADER_PATTERN.search(h):
            gene_col = i
            break
    if gene_col is None:
        # 没有明确表头命中,退化成:哪一列的值看起来最像基因symbol
        best_col, best_frac = None, 0
        for i in range(len(header)):
            vals = [r[i] for r in data_rows[:20] if i < len(r) and r[i]]
            if not vals:
                continue
            frac = sum(1 for v in vals if GENE_PATTERN.match(v)) / len(vals)
            if frac > best_frac:
                best_col, best_frac = i, frac
        if best_frac >= 0.6:
            gene_col = best_col

    coef_col = None
    for pat in COEF_HEADER_PATTERNS:
        for i, h in enumerate(header):
            if pat.search(h):
                coef_col = i
                break
        if coef_col is not None:
            break

    if gene_col is None:
        return None, None, []

    extracted = []
    for r in data_rows:
        if gene_col >= len(r):
            continue
        gene = r[gene_col].strip()
        if not gene or not GENE_PATTERN.match(gene):
            continue
        coef = None
        if coef_col is not None and coef_col < len(r):
            try:
                coef = float(r[coef_col])
            except (ValueError, TypeError):
                coef = r[coef_col].strip() or None
        extracted.append((gene, coef))
    return gene_col, coef_col, extracted


def main():
    candidates = {r["pmcid"]: r for r in load_json(CANDIDATES_PATH)}
    theme_map = {}
    with open(CONVERSION_PATH, encoding="utf-8") as f:
        for row in csv.DictReader(f):
            if row.get("pmcid"):
                theme_map[row["pmcid"]] = {"theme": row["theme"], "pmid": row["pmid"]}

    n_no_candidate = 0
    n_extracted = 0
    n_failed_extraction = 0

    with open(OUT_ROWS, "w", newline="", encoding="utf-8") as fr, \
         open(OUT_SUMMARY, "w", newline="", encoding="utf-8") as fs:
        rw = csv.writer(fr)
        rw.writerow(["pmcid", "pmid", "theme", "gene", "coefficient", "candidate_score", "source", "caption_or_file"])
        sw = csv.writer(fs)
        sw.writerow(["pmcid", "pmid", "theme", "status", "n_genes_extracted", "candidate_score", "source", "caption_or_file"])

        for pmcid, rec in candidates.items():
            meta = theme_map.get(pmcid, {})
            theme = meta.get("theme", "")
            pmid = meta.get("pmid", "")
            cands = rec.get("candidates") or []
            if not cands:
                n_no_candidate += 1
                sw.writerow([pmcid, pmid, theme, "no_candidate", 0, "", "", ""])
                continue

            top = cands[0]
            score = top.get("total_score", top.get("score"))
            label = top.get("caption") or top.get("file") or ""
            if score is not None and score < MIN_SCORE_TO_ATTEMPT:
                n_no_candidate += 1
                sw.writerow([pmcid, pmid, theme, "top_candidate_score_too_low", 0, score, top["source"], label])
                continue
            rows = load_full_rows(pmcid, top)
            gene_col, coef_col, extracted = extract_gene_coef(rows)

            if not extracted:
                n_failed_extraction += 1
                sw.writerow([pmcid, pmid, theme, "extraction_failed", 0, score, top["source"], label])
                continue

            n_extracted += 1
            sw.writerow([pmcid, pmid, theme, "extracted", len(extracted), score, top["source"], label])
            for gene, coef in extracted:
                rw.writerow([pmcid, pmid, theme, gene, coef, score, top["source"], label])

    total = len(candidates)
    print(f"总文章包(有候选记录): {total}")
    print(f"无任何候选文件: {n_no_candidate}")
    print(f"有候选但列识别/抽取失败: {n_failed_extraction}")
    print(f"成功抽取到>=1个基因: {n_extracted}")
    print(f"写入 {OUT_ROWS} 和 {OUT_SUMMARY}")


if __name__ == "__main__":
    main()
