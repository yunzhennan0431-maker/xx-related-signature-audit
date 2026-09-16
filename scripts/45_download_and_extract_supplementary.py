"""
对134篇正文.txt里找不到最终模型系数、但候选池同义重复/miss池复查判断需要进一步核实的论文,
通过PMC OA S3桶下载完整文章包(含补充材料原始文件Excel/docx),尝试从补充表格里提取
最终模型的基因-系数对应关系,补充44号脚本已跑出的167模型归因结果。

复用02号脚本的list_package_files/download_package,不做修改。
"""
import json
import os
import re
import sys
import importlib.util

spec = importlib.util.spec_from_file_location("s02", "analysis_output/scripts/02_pmc_s3_fetch.py")
s02 = importlib.util.module_from_spec(spec)
spec.loader.exec_module(s02)

import pandas as pd

PKG_DIR = "analysis_output/data/pmc_packages"

CLINICAL = {'age', 'gender', 'sex', 'race', 'stage', 'grade', 'characteristic', 'characteristics',
            'covariates', 'variables', 'gene symbol', 'gene', 'genes', 'id', 'ensembl', 'lncrna',
            'symbol', 'description', 'relevance score', 'p value', 'pvalue', 'hr', 'coefficient',
            'coef', 'beta', 'univariate', 'multivariate', '95% ci', 'ci', 'name'}


def gene_like(token):
    t = str(token).strip()
    if not t or len(t) < 2 or len(t) > 25:
        return False
    if not re.match(r'^[A-Za-z][A-Za-z0-9\-\._]{1,24}$', t):
        return False
    return t.lower() not in CLINICAL


def find_coef_column(df):
    """在表格列名里找系数/coefficient/beta/weight这类列,返回列名或None。"""
    for col in df.columns:
        cl = str(col).lower()
        if any(k in cl for k in ('coefficient', 'coef', 'beta', 'weight', 'regression coef')):
            return col
    return None


def extract_from_excel(path):
    try:
        xl = pd.ExcelFile(path)
    except Exception:
        return []
    best_genes = []
    for sheet in xl.sheet_names:
        try:
            df = xl.parse(sheet, header=0)
        except Exception:
            continue
        if df.empty or df.shape[0] < 2:
            continue
        coef_col = find_coef_column(df)
        if coef_col is None:
            continue
        gene_col = df.columns[0]
        genes = [g for g in df[gene_col].astype(str).tolist() if gene_like(g)]
        # 系数列必须是数值且行数与基因数匹配(排除表头误判)
        try:
            numeric_coefs = pd.to_numeric(df[coef_col], errors="coerce")
            n_valid = numeric_coefs.notna().sum()
        except Exception:
            n_valid = 0
        if len(genes) >= 2 and n_valid >= len(genes) * 0.7:
            if len(genes) > len(best_genes):
                best_genes = genes
    return best_genes


def extract_from_docx(path):
    try:
        import docx
        d = docx.Document(path)
    except Exception:
        return []
    best_genes = []
    for t in d.tables:
        if not t.rows:
            continue
        header = [c.text.strip().lower() for c in t.rows[0].cells]
        coef_idx = None
        for i, h in enumerate(header):
            if any(k in h for k in ('coefficient', 'coef', 'beta', 'weight')):
                coef_idx = i
                break
        if coef_idx is None:
            continue
        genes, coefs = [], []
        for row in t.rows[1:]:
            cells = [c.text.strip() for c in row.cells]
            if len(cells) <= coef_idx:
                continue
            if gene_like(cells[0]):
                genes.append(cells[0])
                coefs.append(cells[coef_idx])
        try:
            n_numeric = sum(1 for c in coefs if re.match(r'^[\-−]?\d+\.?\d*$', c.replace(' ', '')))
        except Exception:
            n_numeric = 0
        if len(genes) >= 2 and n_numeric >= len(genes) * 0.7 and len(genes) > len(best_genes):
            best_genes = genes
    return best_genes


def process_paper(pid, keys):
    out_dir = os.path.join(PKG_DIR, pid)
    try:
        s02.download_package(pid, keys)
    except Exception as e:
        return {"pmcid": pid, "status": f"download_failed: {e}", "genes": []}

    candidates = []
    for k in keys:
        fname = k.split("/")[-1]
        ext = fname.lower().rsplit(".", 1)[-1] if "." in fname else ""
        fpath = os.path.join(out_dir, fname)
        if not os.path.exists(fpath):
            continue
        if ext in ("xlsx", "xls"):
            genes = extract_from_excel(fpath)
            if genes:
                candidates.append((fname, genes))
        elif ext == "docx":
            genes = extract_from_docx(fpath)
            if genes:
                candidates.append((fname, genes))

    if not candidates:
        return {"pmcid": pid, "status": "no_coef_table_in_supplementary", "genes": []}
    best = max(candidates, key=lambda x: len(x[1]))
    return {"pmcid": pid, "status": f"found_in:{best[0]}", "genes": best[1]}


def main():
    oa_listing = json.load(open("/tmp/oa_bucket_listing.json"))
    results = []
    for i, (pid, keys) in enumerate(oa_listing.items()):
        r = process_paper(pid, keys)
        results.append(r)
        print(f"[{i+1}/{len(oa_listing)}] {pid}: {r['status']} ({len(r['genes'])} genes)")

    json.dump(results, open("analysis_output/data/supplementary_extraction_results.json", "w"),
              ensure_ascii=False, indent=2)


if __name__ == "__main__":
    main()
