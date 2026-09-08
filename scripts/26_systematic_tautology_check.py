"""
系统性核查全部76个模型是否存在"候选基因池本身被单细胞分析预先限定为特定细胞类型marker基因"
这一同义重复陷阱(此前只核查了3个免疫/TME显式命名模型,发现1例`PMC10791883`)。

方法:对每个模型的PMC全文,用一组红旗短语(marker genes of X / single-cell...marker /
cell-type-specific marker等)做正则扫描,命中的候选逐一人工核查上下文,判断该单细胞marker
分析是(a)用于限定/筛选最终建模的候选基因池本身,还是(b)独立于候选池构建的下游验证/
细胞组成表征分析(此类不构成同义重复,是这类文献里的常规做法)。
"""
import re
from pathlib import Path

import pandas as pd

main = pd.read_csv("analysis_output/data/phase3_meta_analysis_main.csv")
models = main["model"].tolist()
PKG_DIR = Path("analysis_output/data/pmc_packages")

RED_FLAGS = [
    r"marker genes? of\b",
    r"single-cell.{0,80}marker",
    r"scRNA-?seq.{0,80}marker",
    r"cell-type-specific marker",
    r"markers? (?:genes? )?(?:of|for|identif\w+).{0,60}(?:population|cluster|subpopulation|cell type|CAF|macrophage|fibroblast|T cell|endothelial)",
    r"identif\w+.{0,50}marker genes",
]
PATTERN = re.compile("|".join(RED_FLAGS), re.IGNORECASE)

# 人工核查结论(2026-09-06逐一读取红旗命中上下文后判定,见`日志.md`阶段39)。
# 2026-09-07补充:`PMC10140328`当时因原始复核误判为MISS,不在76模型分析集里,没有被这次
# 系统性核查扫描到;后来在一致性抽样检验扩大样本时被独立复核会话发现是原始遗漏(见阶段40),
# 补录为HIT后核实其候选池构建方法(scRNA-seq CAF marker基因∩WGCNA模块基因),确认是
# 与`PMC10791883`完全同构的第2例同义重复。补充验证:用本脚本的红旗正则重新扫描该论文原文,
# 确认能正常命中(4处匹配,含"scRNA-seq marker genes"等),说明当初若这篇论文在76模型分析集
# 里,本脚本本可以正常捕获它——问题出在它当时还没进入分析集,不是扫描规则本身有漏洞。
MANUAL_VERDICTS = {
    "PMC10140328_CAF_BLCA": "确认同义重复(2026-09-07补充)——候选池=scRNA-seq鉴定的CAF marker基因∩bulk WGCNA模块基因交集,与PMC10791883完全同构的候选池构建模板",
    "PMC10046686_pyroptosis_lncRNA_STAD": "假阳性(无关语境,'prognostic marker'泛指)",
    "PMC10294479_autophagy_HCC": "假阳性(无关语境)",
    "PMC10414028_hypoxia_LUAD": "非同义重复——候选池214个HRG来自MSigDB数据库,单细胞TAM marker分析是候选池确定后的下游巨噬细胞亚型表征,不影响候选池构成",
    "PMC10791883_CAF_BRCA": "确认同义重复——候选池=bulk WGCNA模块基因∩单细胞CAF marker基因交集",
    "PMC11382408_cuproptosis_HCC": "非同义重复——候选池是CRG数据库基因∩单细胞'高/低铜死亡活性细胞'差异表达基因∩bulk DEG三方交集,筛选依据是通路活性打分而非特定细胞类型身份,单细胞聚类注释是独立的下游细胞类型表征分析",
    "PMC11760997_senescence_STAD": "假阳性(验证单个基因在成纤维细胞中表达,非候选池构建)",
    "PMC12084325_disulfidptosis_lncRNA_AML": "非同义重复——候选9个DRG来自二硫死亡基因数据库,单细胞marker分析用于AddModuleScore下游打分展示,不影响候选池",
    "PMC12219082_necroptosis_SKCM": "假阳性(标准单细胞聚类细胞类型注释工作流,非候选池构建)",
    "PMC12647489_cholesterol_HCC": "假阳性(单细胞细胞类型比例分析用于下游免疫浸润关联,候选池来自凋亡+胆固醇代谢基因数据库交集)",
    "PMC12891303_CAF_BLCA": "非同义重复——候选池=DEG∩转移相关基因∩基底膜基因三方交集,单细胞分析只验证2个基因",
    "PMC12891315_CAF_UCEC": "候选池确实限定为Wnt通路相关CAF基因,但归因结果落在Tprolif而非CAF/Fibroblast——候选池标签未被证实,是错配证据而非同义重复,按常规计入H1a",
    "PMC13490468_ATIC_HCC": "假阳性(单个基因ATIC的下游单细胞表达验证,非候选池构建)",
    "PMC10287976_chromatin_UVM": "假阳性(红旗命中在参考文献列表条目里,与本研究方法无关)",
    "PMC12971589_CAF_OSCC": "确认同义重复(第4例)——标题即'Single-cell RNA sequencing identifies cancer-associated fibroblast marker genes...',候选池=FindAllMarkers对成纤维细胞聚类簇(cluster 7)直接提取的CAF marker基因全集,无bulk WGCNA交集稀释,比PMC10791883/PMC10140328更彻底的同义重复构造",
    "PMC10900655_lactylation_CRC": "假阳性——红旗命中在标准单细胞聚类细胞类型注释工作流描述里;候选池23个LRG基因'基于既往研究选取'(引用外部文献),单细胞重分析是候选池确定后对最终模型基因的下游TME分布验证(附带报告部分基因预测的细胞类型),不影响候选池构成。附注:作者自述的这一下游单细胞验证结果(11个基因预测于上皮/内皮细胞、2个预测于T/NK、5个预测于单核/巨噬细胞)与本研究z-score归因结果部分吻合(约7/16个可比基因方向一致),可作为第5条任务(作者自述单细胞观察)的第8个候选案例,未来补充整理。",
}

TAUTOLOGICAL_FINAL = ["PMC10791883_CAF_BRCA", "PMC10140328_CAF_BLCA", "PMC12971589_CAF_OSCC"]
# 2026-09-07 batch10(92模型集)红旗扫描:15个新模型均未命中红旗正则,与批次自身候选池构建
# 方式(多来自现成基因集数据库,如MSigDB/CRG/DRG库)一致,未发现新的同义重复。另注:同批次
# 复核中人工阅读全文发现的PMC12575298(HNSCC,候选池=CellChat+CSOmap细胞通讯分析预先限定
# 的APOE+TAM亚群特征基因)是第3例同义重复,但因该癌种本地无参考单细胞数据无法跑归因,不在
# `phase3_meta_analysis_main.csv`的92模型集内,故不会被本脚本的正则扫描覆盖到,只在
# `15_cross_theme_meta_analysis.py`的TAUTOLOGICAL_MODELS里单独留痕,不计入TAUTOLOGICAL_FINAL。


def main_check():
    flagged = []
    for model in models:
        pmcid = model.split("_")[0]
        pkg = PKG_DIR / pmcid
        txt_files = list(pkg.glob(f"{pmcid}*.txt")) if pkg.exists() else []
        if not txt_files:
            continue
        text = txt_files[0].read_text(encoding="utf-8", errors="replace")
        if PATTERN.search(text):
            flagged.append(model)

    print(f"全部模型数: {len(models)}")
    print(f"红旗短语命中数(需人工核查): {len(flagged)}")
    print(f"人工核查后确认同义重复数: {len(TAUTOLOGICAL_FINAL)}")
    print()
    for model in flagged:
        verdict = MANUAL_VERDICTS.get(model, "未记录人工核查结论")
        print(f"{model}: {verdict}")

    unflagged_but_recorded = set(TAUTOLOGICAL_FINAL) - set(flagged)
    if unflagged_but_recorded:
        print("\n警告:以下已确认同义重复的模型未被红旗短语扫描到,需要检查扫描规则覆盖是否完整:")
        print(unflagged_but_recorded)


if __name__ == "__main__":
    main_check()
