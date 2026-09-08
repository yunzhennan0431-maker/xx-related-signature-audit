"""
阶段1第二步(机械化部分):对502个已下载文章包逐个扫描,
定位"最可能包含基因signature列表+系数"的补充材料文件,为后续LLM语义抽取圈定候选范围。

不做语义判断(哪个文件真的是最终风险模型系数表,需要读懂论文方法学,是LLM该做的事),
这一步只做:
  1. 解析正文XML,拿到标题、supplementary-material的标题(caption)与对应文件名。
  2. 关键词打分,给每个候选文件一个粗略优先级(高分不代表一定对,只是缩小人工/LLM要看的范围)。
  3. zip包(常见于MDPI等期刊)解压后列出内部文件名,同样打分。
  4. **论文正文自带的内嵌表格(<table-wrap>)**——21篇试点样本复查时发现,有些论文的最终
     风险系数表根本不在补充材料里,而是直接放在正文的Table 1/2里(如PMC11836339)。早期版本
     只扫了<supplementary-material>,完全漏掉这一类,是真实的覆盖盲区,这版补上。内嵌表格
     不需要另外下载文件,直接把表格内容解析出来存进候选记录里(rows字段),打分时连caption
     和实际表头/数据一起看,不需要再经06脚本内容打分那一步。

输出: analysis_output/data/extraction_candidates.json,每篇论文一条记录。
"""
import json
import re
import zipfile
from pathlib import Path

PKG_DIR = Path(__file__).resolve().parent.parent / "data" / "pmc_packages"
OUT_PATH = Path(__file__).resolve().parent.parent / "data" / "extraction_candidates.json"

# 命中这些词加分(越靠前权重越高),命中这些词减分(通常是单基因分析/无关材料)。
# "coefficient"单独出现容易误中动物实验里的"lung coefficient"(肺系数)这类无关用法
# (试点复查时实测踩过这个坑,见PMC11836339案例),所以拆成"coefficient跟gene同句出现"
# 给高分、"coefficient单独出现"给低分,并对动物实验相关词给负分。
POS_KEYWORDS = [
    (r"coefficients?.{0,60}genes?|genes?.{0,60}coefficients?", 6),
    (r"\bcoefficient", 1),  # 单独出现的coefficient弱信号,可能是别的意思
    (r"\brisk score", 5), (r"\bLASSO", 4), (r"\bmulti-?variate cox", 4),
    (r"\bprognostic (model|signature)", 4), (r"\brisk (model|signature)", 4),
    (r"\bcandidate genes?\b", 3), (r"\bhub genes?\b", 2), (r"\bsignature genes?\b", 4),
    (r"\bunivariate cox", 2), (r"\bgene list\b", 3),
]
NEG_KEYWORDS = [
    (r"\bprimer", -5), (r"\bantibody", -5), (r"\bGSEA enrichment results of the candidate gene\b", -4),
    (r"\bclinical characteristics\b", -2), (r"\bIC50\b", -2), (r"\bquality control\b", -3),
    (r"\bsub-?clustering\b", -3),
    (r"\bbody weight\b", -6), (r"\blung (weight|coefficient)\b", -6),
    (r"\brats?\b", -4), (r"\bmice\b", -4), (r"\bmouse\b", -4), (r"\bxenograft\b", -4),
]


# 交互式批量抽取过程中发现的第三类、也是目前最大宗的假阳性来源:"risk score作为临床协变量,
# 和Age/Gender/Stage/T/N/M等一起做单因素/多因素Cox回归"这类验证表——caption同样含
# "multivariate cox regression""risk score"这些正向关键词会被打高分,但完全不含基因层面系数。
# 区分特征:表格的行标签(通常是第一列)里出现2个以上这种临床变量名,而不是基因symbol。
CLINICAL_VAR_TOKENS = re.compile(
    r"^(age|gender|sex|stage|grade|smoking|metastasis|race|tumor ?size|t ?stage|n ?stage|m ?stage|"
    r"differentiation|bmi|histolog\w*|clinical ?stage|residual tumor|venous invasion|radiation|"
    r"chemotherapy|alcohol)\b",
    re.I,
)


def clinical_covariate_penalty(rows) -> int:
    """行标签(每行第1个非空格子)里命中>=2个临床变量名,判定为'risk score作为协变量'表,扣重分。"""
    hits = 0
    for row in rows[:20]:
        if not row:
            continue
        first_cell = row[0].strip()
        if CLINICAL_VAR_TOKENS.match(first_cell):
            hits += 1
    return -8 if hits >= 2 else 0


def score_text(text: str) -> int:
    score = 0
    low = text.lower()
    for pat, w in POS_KEYWORDS:
        if re.search(pat, low):
            score += w
    for pat, w in NEG_KEYWORDS:
        if re.search(pat, low):
            score += w
    return score


# 兼容旧名字,05内部其他地方还在用score_caption
score_caption = score_text


def parse_table_wrap_rows(table_html: str, max_rows: int = 30):
    """把<table-wrap>里<table>的<tr>/<td|th>解析成行列表,粗暴够用,不追求完美还原rowspan/colspan。"""
    rows = []
    for tr_m in re.finditer(r"<tr[^>]*>(.*?)</tr>", table_html, re.S):
        cells = re.findall(r"<t[dh][^>]*>(.*?)</t[dh]>", tr_m.group(1), re.S)
        cells = [re.sub(r"<[^>]+>", "", c).strip() for c in cells]
        cells = [c for c in cells if c]
        if cells:
            rows.append(cells)
        if len(rows) >= max_rows:
            break
    return rows


def parse_package(pkg_dir: Path) -> dict:
    pmcid = pkg_dir.name
    json_files = list(pkg_dir.glob("*.json"))
    meta = json.loads(json_files[0].read_text(encoding="utf-8")) if json_files else {}
    xml_files = list(pkg_dir.glob("*.xml"))
    xml_text = xml_files[0].read_text(encoding="utf-8", errors="ignore") if xml_files else ""

    candidates = []
    # 1. 命名清晰的独立补充文件(Frontiers这类,caption和文件名分开出现)
    for m in re.finditer(r"<supplementary-material[^>]*>(.*?)</supplementary-material>", xml_text, re.S):
        block = m.group(1)
        cap_m = re.search(r"<caption>(.*?)</caption>", block, re.S)
        caption = re.sub(r"<[^>]+>", " ", cap_m.group(1)).strip() if cap_m else ""
        href_m = re.search(r'xlink:href="([^"]+)"', block)
        href = href_m.group(1) if href_m else None
        if href and re.search(r"\.(xlsx|xls|csv|docx?)$", href, re.I):
            candidates.append({"file": href, "caption": caption, "score": score_caption(caption), "source": "xml_supp"})

    # 1b. 论文正文自带的内嵌表格(<table-wrap>),不是补充材料,是Table 1/2这种。
    #     没有单独的文件,内容直接解析进rows字段;打分同时看caption和实际表头/前几行内容,
    #     行数太多(>50)大概率是DEG/GSEA这类全量结果而非精炼过的最终模型,扣分。
    for m in re.finditer(r"<table-wrap[^>]*>(.*?)</table-wrap>", xml_text, re.S):
        block = m.group(1)
        cap_m = re.search(r"<caption>(.*?)</caption>", block, re.S)
        caption = re.sub(r"<[^>]+>", " ", cap_m.group(1)).strip() if cap_m else ""
        table_m = re.search(r"<table[^>]*>(.*?)</table>", block, re.S)
        rows = parse_table_wrap_rows(table_m.group(1)) if table_m else []
        content_text = " ".join(" ".join(r) for r in rows[:3])
        score = score_text(caption) + score_text(content_text)
        score += clinical_covariate_penalty(rows)
        if len(rows) > 50:
            score -= 3
        candidates.append({
            "file": None, "caption": caption, "rows": rows,
            "score": score, "total_score": score, "content_score": None,
            "source": "xml_table_wrap",
        })

    # 2. zip包(MDPI等):解压列出内部文件,caption信息有限,只能靠文件名/sheet名弱打分
    for zpath in pkg_dir.glob("*.zip"):
        try:
            with zipfile.ZipFile(zpath) as zf:
                extract_dir = pkg_dir / (zpath.stem + "_extracted")
                extract_dir.mkdir(exist_ok=True)
                zf.extractall(extract_dir)
                for name in zf.namelist():
                    if re.search(r"\.(xlsx|xls|csv|docx?)$", name, re.I):
                        candidates.append({
                            "file": str((extract_dir / name).relative_to(pkg_dir)),
                            "caption": f"(来自zip: {zpath.name}, 内部文件名: {name})",
                            "score": score_caption(name),
                            "source": "zip",
                        })
        except zipfile.BadZipFile:
            pass

    candidates.sort(key=lambda c: -c["score"])
    return {
        "pmcid": pmcid,
        "pmid": meta.get("pmid"),
        "title": meta.get("title"),
        "doi": meta.get("doi"),
        "n_candidates": len(candidates),
        "candidates": candidates,
    }


def main():
    records = []
    for pkg_dir in sorted(PKG_DIR.iterdir()):
        if not pkg_dir.is_dir():
            continue
        try:
            rec = parse_package(pkg_dir)
        except Exception as e:
            rec = {"pmcid": pkg_dir.name, "error": str(e)}
        records.append(rec)

    OUT_PATH.write_text(json.dumps(records, ensure_ascii=False, indent=2), encoding="utf-8")

    n_with_candidates = sum(1 for r in records if r.get("n_candidates", 0) > 0)
    n_with_positive = sum(1 for r in records if r.get("candidates") and r["candidates"][0]["score"] > 0)
    print(f"总文章包: {len(records)}")
    print(f"至少有1个候选表格文件: {n_with_candidates}")
    print(f"最高分候选文件分数>0(较有把握): {n_with_positive}")
    print(f"写入 {OUT_PATH}")


if __name__ == "__main__":
    main()
