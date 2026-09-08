"""
阶段3:跨主题meta分析。设计方案见`跨主题meta分析_设计方案.md`,这里是执行代码。

核心统计:每个模型做精确二项检验(H0: 基因归因在该数据集K个候选细胞类型间等概率随机),
配合归一化熵集中度作为连续型补充指标,19个模型的p值按同一自然家族做BH-FDR校正。
另做排除<20细胞候选类型的敏感性分析,以及众数细胞类型的三大类归类统计。
"""
from collections import Counter

import numpy as np
import pandas as pd
from scipy import stats

DATA_DIR = "analysis_output/data/"
COMBINED_PATH = DATA_DIR + "phase2_pilot_all_combined.csv"

# 每个模型对应的数据集,数据集的候选细胞类型及各自参考细胞数(从各阶段2脚本运行输出里摘录)
DATASET_CELLTYPE_COUNTS = {
    "GSE81861": {  # CRC/结肠癌共用
        "Epithelial": 432, "Tcell": 45, "Bcell": 35, "Macrophage": 29,
        "Fibroblast": 26, "Endothelial": 6, "MastCell": 4,
    },
    "GSE176078": {  # BRCA
        "T-cells": 35214, "Cancer Epithelial": 24489, "Myeloid": 9675, "Endothelial": 7605,
        "CAFs": 6573, "PVL": 5423, "Normal Epithelial": 4355, "Plasmablasts": 3524, "B-cells": 3206,
    },
    "GSE151530": {  # HCC
        "T cells": 19587, "Malignant cells": 17164, "TAMs": 5171, "TECs": 2646,
        "B cells": 2291, "CAFs": 1459,
    },
    "GSE167297": {  # STAD
        "CD8T": 8749, "B": 4628, "Epithelial": 2797, "Plasma": 2490, "DC": 1119,
        "Endothelial": 854, "Mono/Macro": 854, "Fibroblasts": 655, "Mast": 318,
    },
    "SKCM_GSE120575": {  # 黑色素瘤
        "CD8Tex": 3922, "CD4Tconv": 3341, "CD8T": 2843, "B": 1467, "Mono/Macro": 1398,
        "NK": 1037, "Tprolif": 850, "Treg": 839, "Plasma": 310, "DC": 284,
    },
    "Glioma_GSE131928": {  # 胶质瘤
        "MES-like Malignant": 5203, "Mono/Macro": 3969, "OPC-like Malignant": 1798,
        "NPC-like Malignant": 1368, "Oligodendrocyte": 394, "Malignant": 367,
        "AC-like Malignant": 262, "CD8Tex": 197,
    },
    "AML_GSE116256": {  # 急性髓系白血病
        "Malignant": 12489, "CD4Tconv": 6323, "Progenitor": 3447, "Mono/Macro": 2893,
        "HSC": 2523, "Promonocyte": 2029, "EryPro": 2019, "NK": 1738, "CD8T": 1524,
        "GMP": 1314, "Plasma": 1227, "B": 422, "Tprolif": 400,
    },
    "KIRC_GSE159115": {  # 肾透明细胞癌
        "Malignant": 9027, "Epithelial": 5477, "Mono/Macro": 4730, "Endothelial": 3798,
        "Pericytes": 3059, "CD8T": 1247, "Plasma": 249, "Erythroblasts": 82,
    },
    "OS_GSE162454": {  # 骨肉瘤
        "Mono/Macro": 16682, "Malignant": 9443, "Fibroblasts": 8006, "Osteoblasts": 4131,
        "CD4Tconv": 4085, "CD8Tex": 2010, "Plasmocytes": 1208, "Endothelial": 979,
    },
    "UCEC_GSE139555": {  # 子宫内膜癌
        "CD4Tconv": 5319, "CD8T": 3892, "Treg": 1748, "CD8Tex": 1449, "Tprolif": 284, "Fibroblasts": 66,
    },
    "PAAD_GSE154778": {  # 胰腺癌
        "Malignant": 9519, "Fibroblasts": 2218, "Mono/Macro": 1578, "Epithelial": 1069,
        "CD8T": 421, "Plasma": 148,
    },
    "NSCLC_GSE131907": {  # 肺腺癌——注:该数据集匹配到的细胞子集只标注了3种免疫细胞大类,
        # 不含恶性/上皮/基质细胞,candidate типes覆盖极不完整,结果解读需要格外谨慎(见文档第5节)。
        "CD4Tconv": 42763, "B": 20287, "CD8T": 4174,
    },
    "CESC_GSE168652": {  # 宫颈癌
        "Malignant": 12034, "SMC": 4761, "Endometrial stromal cells": 2737, "Endothelial": 1781,
        "Fibroblasts": 1153, "CD8T": 365, "Mono/Macro": 167,
    },
    "BLCA_GSE130001": {  # 膀胱癌——注:该数据集meta.tsv里全部细胞的malignancy字段都标"Stromal cells"
        # (即使major-lineage标Epithelial),推测这是TISCH2原始收录时只保留了基质compartment的
        # 子集,不含免疫细胞和明确的恶性细胞,K=4且候选类型覆盖极不完整,与NSCLC_GSE131907
        # (只有免疫细胞、缺基质/恶性)构成相反方向的覆盖偏差,解读需格外谨慎(见文档第5节)。
        "Epithelial": 3772, "Fibroblasts": 144, "Endothelial": 143, "Myofibroblasts": 70,
    },
    # 弥补"15个新癌种候选一直悬着"这一局限,2026-09新下载8个TISCH2数据集覆盖其中11个模型
    # (另有2个:肾上腺皮质癌基因数太少且无数据集覆盖、Wilms瘤TISCH2无对应数据集,继续搁置)。
    "HNSC_GSE139324": {  # 头颈鳞癌;与PMC12575298(候选池同义重复第3例)原文使用的是同一数据集
        "CD4Tconv": 36357, "Mono/Macro": 27599, "CD8Tex": 16478, "B": 16439, "CD8T": 14799,
        "NK": 7235, "Treg": 4579, "Tprolif": 3938, "Plasma": 1581, "DC": 1218, "Mast": 498,
    },
    "OSCC_GSE172577": {  # 口腔鳞癌(头颈鳞癌亚型,TISCH2有专门数据集,不与HNSC_GSE139324合并)
        "Keratinocytes": 19314, "Ductal": 8529, "Malignant": 7914, "Gland": 7556,
        "Mono/Macro": 5422, "NK": 2674, "CD8T": 2416, "Endothelial": 1652, "Mast": 1140,
        "Epithelial": 465, "Treg": 421,
    },
    "UVM_GSE139829": {  # 葡萄膜黑色素瘤(与皮肤黑色素瘤SKCM_GSE120575生物学不同,不复用)
        "Malignant": 79105, "CD4Tconv": 7835, "CD8Tex": 7252, "Mono/Macro": 5663,
        "Plasma": 1986, "Endothelial": 909, "CD8T": 481, "B": 472,
    },
    "ALL_GSE132509": {  # 急性淋巴细胞白血病
        "Malignant": 21370, "Erythrocytes": 4281, "CD8T": 4145, "CD4Tconv": 3692,
        "B": 1874, "Mono/Macro": 1412, "Tprolif": 1162,
    },
    "OV_GSE154600": {  # 卵巢癌
        "Epithelial": 8669, "Mono/Macro": 8333, "Fibroblasts": 7622, "CD8T": 6471,
        "CD8Tex": 4355, "Malignant": 2320, "Plasma": 1298, "Endothelial": 1186,
        "Myofibroblasts": 1067, "B": 763, "Tprolif": 499,
    },
    "PRAD_GSE172301": {  # 前列腺癌
        "Fibroblasts": 27270, "Epithelial": 26135, "Endothelial": 17595, "Mono/Macro": 13052,
        "CD8T": 11551, "B": 3930, "Mast": 3416, "SMC": 733, "Plasma": 339,
    },
    "MB_GSE155446": {  # 髓母细胞瘤——注:候选类型K=4明显偏窄(只标注恶性/单核巨噬/CD8T/中性粒细胞
        # 四类,不含B/T CD4/内皮/基质),解读需谨慎,与NSCLC_GSE131907、BLCA_GSE130001同类局限。
        "Malignant": 32926, "Mono/Macro": 2665, "CD8T": 1590, "Neutrophils": 264,
    },
}
# OS_GSE162454/NSCLC_GSE131907的候选类型清单复用上面已有条目,不重复定义

MODEL_TO_DATASET = {
    "PMC13190656_ferroptosis_CRC": "GSE81861",
    "PMC11656000_pyroptosis_CRC": "GSE81861",
    "PMC11384320_m6A_ferroptosis_COAD": "GSE81861",
    "PMC11297001_oxidative_stress_COAD": "GSE81861",
    "PMC13083393_autophagy_BRCA": "GSE176078",
    "PMC13406258_PCDgenes_BRCA": "GSE176078",
    "PMC11855622_mitophagyE3_BRCA": "GSE176078",
    "PMC12647489_cholesterol_HCC": "GSE151530",
    "PMC13490468_ATIC_HCC": "GSE151530",
    "PMC9906058_ferroptosis_cuproptosis_HCC": "GSE151530",
    "PMC11577587_glutamine_HCC": "GSE151530",
    "PMC11760997_senescence_STAD": "GSE167297",
    "PMC12399450_mitochondrial_STAD": "GSE167297",
    "PMC12219082_necroptosis_SKCM": "SKCM_GSE120575",
    "PMC12905907_disulfidptosis_SKCM": "SKCM_GSE120575",
    "PMC9420977_cuproptosis_Glioma": "Glioma_GSE131928",
    "PMC12239816_hypoxia_Glioma": "Glioma_GSE131928",
    "PMC12378281_disulfidptosis_AML": "AML_GSE116256",
    "PMC12402452_exosome_AML": "AML_GSE116256",
    "PMC10184413_cuproptosis_ferroptosis_KIRC": "KIRC_GSE159115",
    "PMC12047218_glycolysis_OS": "OS_GSE162454",
    "PMC10469729_glutamine_UCEC": "UCEC_GSE139555",
    "PMC11504279_hypoxia_PAAD": "PAAD_GSE154778",
    "PMC13333213_pyroptosis_NSCLC": "NSCLC_GSE131907",
    "PMC10564413_senescence_CESC": "CESC_GSE168652",
    "PMC12906242_m6A_OS": "OS_GSE162454",
    "PMC10475834_cuproptosis_lncRNA_LUAD": "NSCLC_GSE131907",
    "PMC10085621_hypoxia_NSCLC": "NSCLC_GSE131907",
    "PMC11564170_cuproptosis_AML": "AML_GSE116256",
    "PMC12084325_disulfidptosis_lncRNA_AML": "AML_GSE116256",
    "PMC9086515_ferroptosis_Glioma": "Glioma_GSE131928",
    "PMC9048552_autophagy_lncRNA_SKCM": "SKCM_GSE120575",
    "PMC9524961_hypoxia_HCC": "GSE151530",
    "PMC10294479_autophagy_HCC": "GSE151530",
    "PMC10164321_pyroptosis_HCC": "GSE151530",
    "PMC11491388_disulfidptosis_lncRNA_HCC": "GSE151530",
    "PMC11914417_disulfidptosis_CRC": "GSE81861",
    # batch6(第三轮扩大检索复核)追加,11模型,零新增数据成本(另有3个lncRNA模型因参考
    # 面板匹配基因数=1,判定为不可分析而排除,见日志)
    "PMC9767988_pyroptosis_lncRNA_STAD": "GSE167297",
    "PMC9614251_glycolysis_lncRNA_STAD": "GSE167297",
    "PMC11178724_cuproptosis_methylation_SKCM": "SKCM_GSE120575",
    "PMC9342864_autophagy_Glioma": "Glioma_GSE131928",
    "PMC9561419_glycolysis_immune_LUSC": "NSCLC_GSE131907",
    "PMC10414028_hypoxia_LUAD": "NSCLC_GSE131907",
    "PMC9550247_autophagy_LUAD": "NSCLC_GSE131907",
    "PMC11382408_cuproptosis_HCC": "GSE151530",
    "PMC8605142_glycolysis_HCC": "GSE151530",
    "PMC8377503_glycolysis_colon": "GSE81861",
    "PMC9341065_immune_autophagy_BLCA": "BLCA_GSE130001",
    # batch7(第四轮扩大检索复核,score降到3-4区间)追加,12模型,零新增数据成本(另有2个模型
    # 因参考面板匹配基因数=0/1不可分析而排除,见日志)
    "PMC9805482_autophagy_OS": "OS_GSE162454",
    "PMC9662042_autophagy_lncRNA_KIRC": "KIRC_GSE159115",
    "PMC9124146_autophagy_lncRNA_PAAD": "PAAD_GSE154778",
    "PMC9097333_autophagy_GBM": "Glioma_GSE131928",
    # 注:该模型实际是ESCA+STAD+LIHC+CRC四种消化系统癌种联合建模(原文18个候选HRlncRNA取
    # 四癌种交集后LASSO筛选),不对应单一癌种,这里借用STAD单细胞参考做归因只是权宜之计,
    # 结果解读需要注明这个方法学错配,讨论部分单独说明,不能和其他单癌种模型同等看待。
    "PMC9092832_hypoxia_lncRNA_digestive_panCancer": "GSE167297",
    "PMC10046686_pyroptosis_lncRNA_STAD": "GSE167297",
    "PMC12210501_ferroptosis_LUAD": "NSCLC_GSE131907",
    "PMC7847640_glycolysis_lncRNA_UCEC": "UCEC_GSE139555",
    "PMC8812245_glycolysis_READ": "GSE81861",
    "PMC11483459_cuproptosis_CRC": "GSE81861",
    "PMC11920466_COAD": "GSE81861",
    "PMC9035888_autophagy_BRCA": "GSE176078",
    # batch8(第五轮扩大检索复核,score降到1-2区间)追加,仅1模型零新增数据成本
    # (35篇复核命中率骤降到8.6%,另有2个候选因参考面板匹配基因数=0/1不可分析而排除)
    "PMC8994655_autophagy_lncRNA_STAD": "GSE167297",
    # batch9(主题词表扩充后新一轮检索,命中率62.9%的富矿)追加,15模型零新增数据成本
    # (另有1个模型因NSCLC_GSE131907候选面板极窄+lncRNA匹配基因数=1不可分析而排除)
    "PMC7880321_immune_BLCA": "BLCA_GSE130001",
    "PMC12891303_CAF_BLCA": "BLCA_GSE130001",
    "PMC11634757_lactylation_STAD": "GSE167297",
    "PMC9578220_chromatin_immune_UCEC": "UCEC_GSE139555",
    "PMC12891315_CAF_UCEC": "UCEC_GSE139555",
    "PMC9614380_aminoacid_KIRC": "KIRC_GSE159115",
    "PMC11310873_lncRNA_LUAD": "NSCLC_GSE131907",
    "PMC9843944_autophagy_LUAD": "NSCLC_GSE131907",
    "PMC8379743_m6A_LUAD": "NSCLC_GSE131907",
    "PMC9618960_ICD_LGG": "Glioma_GSE131928",
    "PMC10932653_chromatin_HCC": "GSE151530",
    "PMC10662873_anoikis_HCC": "GSE151530",
    "PMC10531437_lipid_ERstress_CRC": "GSE81861",
    "PMC11085026_ESR_BRCA": "GSE176078",
    "PMC10791883_CAF_BRCA": "GSE176078",
    # 2026-09-07:扩大批次9抽取一致性抽样检验时(样本20→35)发现的1处真实抽取遗漏——
    # 原始批次9复核只读了Sheet1(候选池,983行),没往下翻到Sheet2(9行,7基因最终系数),
    # 独立复核会话读到了完整文件,揪出这个假阴性。补录为命中,零新增数据成本复用BLCA_GSE130001。
    "PMC10140328_CAF_BLCA": "BLCA_GSE130001",
    # batch10(主题词表二次扩充后检索,命中率48.6%)追加,15模型零新增数据成本(另有2个命中因
    # 是全新癌种[HNSCC头颈鳞癌、ALL急性淋巴细胞白血病]本地无参考单细胞数据,暂不纳入本轮分析,
    # 其中HNSCC模型[PMC12575298]全文核查发现候选池同样是单细胞marker预先限定的第3例同义重复
    # 候选池,见下方TAUTOLOGICAL_MODELS注释,虽不能跑归因但记入同义重复模板的系统性证据)
    "PMC12354743_NETs_Glioma": "Glioma_GSE131928",
    "PMC10980369_ERstress_SKCM": "SKCM_GSE120575",
    "PMC9361349_DNArepair_SKCM": "SKCM_GSE120575",
    "PMC12179157_PANoptosis_SKCM": "SKCM_GSE120575",
    "PMC10076677_anoikis_lncRNA_OS": "OS_GSE162454",
    "PMC13355270_mitophagy_BLCA": "BLCA_GSE130001",
    "PMC12055729_lactylation_UCEC": "UCEC_GSE139555",
    "PMC9870742_anoikis_LUAD": "NSCLC_GSE131907",
    "PMC8291813_cGASSTING_STAD": "GSE167297",
    "PMC9468367_ICD_lncRNA_HCC": "GSE151530",
    "PMC9720162_ICD_lncRNA_HCC": "GSE151530",
    "PMC10616163_chromatin_HCC": "GSE151530",
    "PMC11866598_telomere_HCC": "GSE151530",
    "PMC9636133_mitophagy_CRC": "GSE81861",
    # 2026-09 batch11复核中发现:PMC13114778(乳腺癌,批次10计入)的23基因系数表与
    # PMC10900655(结直肠癌,2024年发表于J Transl Med,早于PMC13114778的2026年发表)逐基因
    # 逐系数完全一致(精确到小数点后15位),几乎可以确定是同一模型输出的重复使用/误标癌种,
    # 而非两个独立推导的模型。按"每个模型必须独立衍生"这一基本假设,不能计为两个样本点。
    # PMC10900655发表更早,判定为原始/合法模型并按其实际声称的结直肠癌重新归因(见下方
    # PMC10900655_lactylation_CRC条目、`31_pmc10900655_crc_attribution.py`);PMC13114778
    # 判定为疑似重复/误标癌种,已从模型集中移除,不再作为独立乳腺癌模型统计(原BRCA归因结果
    # 未曾达到FDR显著,移除不影响此前报告的显著模型计数,但会影响总N、功效不足计数等描述性统计)。
    "PMC10900655_lactylation_CRC": "GSE81861",
    # batch11(继续扩大候选池复核,命中率25.7%)追加,6模型零新增数据成本(另有1个命中
    # PMC10932661因8个lncRNA只匹配上1个基因、n=1不可分析而排除;1个命中PMC12973827因
    # 神经母细胞瘤本地无参考数据暂未纳入)
    "PMC10807923_anoikis_OV": "OV_GSE154600",
    "PMC10469923_anoikis_PRAD": "PRAD_GSE172301",
    "PMC10683111_anoikis_lncRNA_HNSCC": "HNSC_GSE139324",
    "PMC9465161_m7G_lncRNA_COAD": "GSE81861",
    "PMC11411761_anoikis_BRCA": "GSE176078",
    "PMC10339815_DDR_lncRNA_STAD": "GSE167297",
    # 2026-09-07 batch10核查表格式候选(多sheet)时发现的第2处遗漏(与批次9的PMC10140328同类
    # 问题,原始复核只看了第一个sheet的候选基因数据库列表,没往后翻到第5列"Multivariate Cox
    # significant genes"最终10基因系数):补录为命中,零新增数据成本复用KIRC_GSE159115。
    "PMC9523360_telomere_KIRC": "KIRC_GSE159115",
    # 2026-09 新下载8个数据集覆盖的11个新癌种模型中,9个可纳入定量分析(2个不可纳入:
    # PMC12575298因候选池同义重复单独处理见TAUTOLOGICAL_MODELS注释;PMC9945660因10个lncRNA
    # 只匹配上1个基因,n=1对二项检验无意义,按历史惯例判定不可分析而排除,不计入此处)。
    "PMC8791745_hypoxia_lncRNA_HNSCC": "HNSC_GSE139324",
    "PMC12971589_CAF_OSCC": "OSCC_GSE172577",
    "PMC10599730_ERstress_OSCC": "OSCC_GSE172577",
    "PMC10287976_chromatin_UVM": "UVM_GSE139829",
    "PMC9709208_ICD_UVM": "UVM_GSE139829",
    "PMC12078231_ubiquitination_ALL": "ALL_GSE132509",
    "PMC8464158_glycolysis_lncRNA_OV": "OV_GSE154600",
    "PMC11550951_cuproptosis_lncRNA_PRAD": "PRAD_GSE172301",
    "PMC12504172_nucleotide_MB": "MB_GSE155446",
}

# 众数细胞类型的三大类归类(人工判断,留痕以便复核)。
# 注:GSE81861的"Epithelial"是正常+肿瘤上皮混合(前序课题原始数据结构如此),不是纯恶性细胞,
#     单独归为Epithelial_mixed而非Malignant;AML的HSC/GMP/Promonocyte/EryPro/Progenitor这类
#     造血分化阶段类型,在白血病里克隆本身可能跨越多个分化阶段,不等同于典型实体瘤"正常免疫细胞"
#     概念,但该数据集未标注它们属于恶性克隆,故仍按"Immune_Hematopoietic"处理,讨论部分需要说明
#     这个AML特有的歧义。
CELLTYPE_CATEGORY = {
    "Epithelial": "Epithelial_mixed", "Cancer Epithelial": "Malignant", "Normal Epithelial": "Epithelial_mixed",
    "Malignant cells": "Malignant", "Malignant": "Malignant",
    "AC-like Malignant": "Malignant", "NPC-like Malignant": "Malignant",
    "MES-like Malignant": "Malignant", "OPC-like Malignant": "Malignant",
    "Endothelial": "Stromal_Vascular", "Fibroblast": "Stromal_Vascular", "Fibroblasts": "Stromal_Vascular",
    "CAFs": "Stromal_Vascular", "PVL": "Stromal_Vascular", "TECs": "Stromal_Vascular",
    "Oligodendrocyte": "Stromal_Vascular", "Pericytes": "Stromal_Vascular", "SMC": "Stromal_Vascular",
    "Endometrial stromal cells": "Stromal_Vascular",
    # Osteoblasts:该数据集(骨肉瘤OS_GSE162454)把"Malignant"和"Osteoblasts"分开标注,
    # 后者视为非恶性的成骨谱系基质细胞,归入Stromal_Vascular而非Malignant,但要承认这个
    # 归类本身有争议——骨肉瘤本来就是成骨细胞谱系来源的肿瘤,"正常成骨细胞"和"恶性细胞"
    # 的转录组边界不像其他癌种的"肿瘤vs基质"那样清晰,讨论部分需要单独说明。
    "Osteoblasts": "Stromal_Vascular",
    "Myofibroblasts": "Stromal_Vascular",
    "Erythroblasts": "Immune_Hematopoietic", "Plasmocytes": "Immune_Hematopoietic",
    "Tcell": "Immune_Hematopoietic", "T cells": "Immune_Hematopoietic", "T-cells": "Immune_Hematopoietic",
    "CD8T": "Immune_Hematopoietic", "CD8Tex": "Immune_Hematopoietic", "CD4Tconv": "Immune_Hematopoietic",
    "Treg": "Immune_Hematopoietic", "Tprolif": "Immune_Hematopoietic",
    "Bcell": "Immune_Hematopoietic", "B cells": "Immune_Hematopoietic", "B-cells": "Immune_Hematopoietic", "B": "Immune_Hematopoietic",
    "Plasma": "Immune_Hematopoietic", "Plasmablasts": "Immune_Hematopoietic",
    "Macrophage": "Immune_Hematopoietic", "TAMs": "Immune_Hematopoietic", "Mono/Macro": "Immune_Hematopoietic",
    "Myeloid": "Immune_Hematopoietic", "DC": "Immune_Hematopoietic",
    "Mast": "Immune_Hematopoietic", "MastCell": "Immune_Hematopoietic", "NK": "Immune_Hematopoietic",
    "HSC": "Immune_Hematopoietic", "GMP": "Immune_Hematopoietic", "Promonocyte": "Immune_Hematopoietic",
    "EryPro": "Immune_Hematopoietic", "Progenitor": "Immune_Hematopoietic",
    # 2026-09新增数据集(HNSC/OSCC/UVM/ALL/OV/PRAD/MB)引入的细胞类型标注:
    "Erythrocytes": "Immune_Hematopoietic", "Neutrophils": "Immune_Hematopoietic",
    # Keratinocytes/Gland/Ductal是OSCC_GSE172577里非恶性的正常上皮谱系分化状态(角质形成细胞、
    # 腺体/导管上皮),与该数据集里单独标注的"Malignant"区分,按"Epithelial"同等逻辑归入
    # Epithelial_mixed而非Malignant。
    "Keratinocytes": "Epithelial_mixed", "Gland": "Epithelial_mixed", "Ductal": "Epithelial_mixed",
}

MIN_CELLS_SENSITIVITY = 20

# 2026-09-06发现:异质性检验发现"ImmuneTME"主题组显著比例偏高后,核查发现其中至少1个模型
# 的候选基因池本身就是靠单细胞分析预先限定为某种TME细胞marker基因(不是广义通路/过程基因池),
# 这种情况下归因结果落在该细胞类型属于同义重复,不构成"命名-机制脱节"证据,必须和其余"广义
# 候选池→意外集中于某TME细胞"的模型分开统计,不能混在H1a里一起计数。
# 已逐一核实候选池构建方法(见`日志.md`阶段38):
#  - PMC10791883(CAF_BRCA):候选池="bulk WGCNA模块基因 ∩ 单细胞CAF marker基因",
#    归因结果恰好落在CAFs——同义重复,标记为Tautological。
#  - PMC12891315(CAF_UCEC):候选池同样限定为"Wnt-related CAF基因",但归因结果落在Tprolif
#    (增殖性T细胞)而非CAF/Fibroblast——候选池的"CAF"标签本身未被单细胞归因证实,不是同义
#    重复,反而是更强的错配证据,不做特殊处理,按常规H1a/Immune_Hematopoietic计入。
#  - PMC12891303(CAF_BLCA):候选池="DEG ∩ 转移相关基因 ∩ 基底膜基因"三方交集,不限定于CAF
#    marker,单细胞分析只是后续对2个基因的机制验证,不影响候选池构成——不做特殊处理。
TAUTOLOGICAL_MODELS = {
    "PMC10791883_CAF_BRCA": "候选基因池经单细胞分析预先限定为CAF marker基因,归因于CAFs属同义重复,不计入H1a的脱节证据",
    # 2026-09-07确认第2例:方法完全同构(scRNA-seq鉴定CAF marker基因 ∩ bulk WGCNA模块基因),
    # 归因于Fibroblasts同样属同义重复。两例均是"WGCNA∩单细胞CAF marker"这一特定候选池构建
    # 模板,提示这一具体模板本身系统性容易产生同义重复,不是随机个案,论文讨论部分需要指出。
    "PMC10140328_CAF_BLCA": "候选基因池同构造方法(单细胞CAF marker∩bulk WGCNA模块基因),归因于Fibroblasts属同义重复,不计入H1a的脱节证据",
    # 2026-09-07 batch10全文核查发现第3例,构造模板不同(不是WGCNA∩marker交集,而是"细胞通讯+
    # 空间邻近分析→直接选定单一TME亚群的特征基因作为候选池"):PMC12575298(HNSCC头颈鳞癌)
    # 原文明确写道"我们从与IL32+NK亚群有显著相互作用且空间可及的APOE+TAM亚群中选取620个
    # 特征基因作为候选基因"(见原文第258行),候选池已被单细胞分析预先限定为APOE+TAM(髓系/
    # 巨噬细胞谱系)marker基因,若归因结果落在髓系/巨噬细胞类别同样属同义重复。该模型因HNSCC
    # 本地无参考单细胞数据集,暂无法跑z-score归因验证是否真落在该类别,故未加入下方
    # MODEL_TO_DATASET/DATASET_CELLTYPE_COUNTS,只作为"候选池同义重复"证据独立记录:第3个
    # 案例出现了第2种不同的构造模板(细胞通讯/空间邻近分析限定候选池),提示这类同义重复陷阱
    # 不局限于"WGCNA∩单细胞marker"这一种具体方法,而是"用单细胞分析预先限定候选池到特定细胞
    # 类型"这一更一般的做法本身就容易产生同义重复,论文讨论部分需要更新为"至少两类不同构造
    # 模板"而非"同一模板的两个案例"。
    "PMC12575298_NKTAM_HNSCC": "候选基因池经细胞通讯(CellChat)+空间邻近(CSOmap)分析预先限定为APOE+TAM亚群特征基因,候选池构建方法本身已构成同义重复的方法学缺陷,与构建方式无关,故仍排除在H1a/H1b定量计数之外。"
        "2026-09补充:下载该论文原文实际使用的同一数据集(HNSC_GSE139324)后用本研究z-score方法独立跑了归因(见`30_new_celltype_datasets_attribution.py`),**结果并未清晰复现论文'APOE+TAM髓系细胞富集'的框架性结论**——23个模型基因分布较为分散(Treg 5、Mast 5、Tprolif 4、Mono/Macro仅3、DC 3、Plasma 2、NK 1),Mono/Macro不是众数类别。这一经验结果不改变'候选池构建方法本身存在同义重复缺陷'这一判断(缺陷在于候选池限定步骤本身,不取决于LASSO最终选中的子集是否真的表现出该细胞类型特征),但如实说明:即使候选池被人为限定,最终入选的具体基因子集也未必在独立复核时表现出与限定标签一致的富集模式,提示'候选池同义重复'和'归因结果与候选池标签一致'是两个可能脱钩的问题,不能假设两者总是同步发生。该模型仍不纳入93→102模型的定量分析集,该经验结果只作为方法学讨论的补充证据。",
    # 2026-09弥补"15个新癌种候选"局限、新下载8个数据集跑归因时,红旗正则扫描发现第4例:
    # PMC12971589(口腔鳞癌CAF模型)标题本身即"Single-cell RNA sequencing identifies
    # cancer-associated fibroblast marker genes..."——用FindAllMarkers对成纤维细胞聚类簇
    # (cluster 7)直接提取CAF marker基因作为候选基因池全集,不经过bulk WGCNA交集这一步,是比
    # PMC10791883/PMC10140328(需要bulk∩单细胞两步交集)更彻底的同义重复构造(候选池100%
    # 来自单细胞marker,无bulk环节稀释)。已下载对应数据集OSCC_GSE172577独立跑归因(见
    # `30_new_celltype_datasets_attribution.py`),结果同样未复现"落在Fibroblasts"这一框架性
    # 结论(8个基因分布于Endothelial 2、Malignant 2、Treg 1、Mono/Macro 1、Epithelial 1、
    # Mast 1,无一落在Fibroblasts)——与PMC12575298一致,再次印证"候选池同义重复"与"归因结果
    # 匹配候选池标签"是两个可能脱钩的问题。纳入MODEL_TO_DATASET计入102模型描述性统计,但排除
    # 在H1a/H1b脱节证据计数之外。
    "PMC12971589_CAF_OSCC": "候选池=单细胞FindAllMarkers对成纤维细胞聚类簇直接提取的CAF marker基因全集(无bulk WGCNA交集稀释,比PMC10791883/PMC10140328更彻底的同义重复构造),不计入H1a/H1b;经验归因结果未落在Fibroblasts(落在Endothelial/Malignant/Treg/Mono-Macro/Epithelial/Mast),与PMC12575298一样呈现'候选池同义重复'和'归因匹配候选池标签'脱钩的模式",
}


def entropy_concentration(counts, K):
    p = np.array(counts) / np.sum(counts)
    p = p[p > 0]
    H = -np.sum(p * np.log(p))
    Hmax = np.log(K)
    return 1 - H / Hmax if Hmax > 0 else np.nan


def analyze(df, exclude_rare=False):
    rows = []
    for model, g in df.groupby("model"):
        dataset = MODEL_TO_DATASET[model]
        celltype_counts = DATASET_CELLTYPE_COUNTS[dataset]
        candidate_types = set(celltype_counts.keys())
        if exclude_rare:
            candidate_types = {ct for ct, n in celltype_counts.items() if n >= MIN_CELLS_SENSITIVITY}
            g = g[g["top_celltype"].isin(candidate_types)]
        K = len(candidate_types)
        n_genes = len(g)
        if n_genes == 0 or K == 0:
            continue

        vc = g["top_celltype"].value_counts()
        top_type = vc.index[0]
        top_count = vc.iloc[0]
        max_prop = top_count / n_genes

        # 精确二项检验: H0 p=1/K
        pval = stats.binomtest(top_count, n_genes, 1 / K, alternative="greater").pvalue

        ent_conc = entropy_concentration(vc.values, K)
        category = CELLTYPE_CATEGORY.get(top_type, "Unclassified")

        rows.append({
            "model": model, "dataset": dataset, "n_genes": n_genes, "K_celltypes": K,
            "top_celltype": top_type, "top_count": top_count, "max_prop": round(max_prop, 3),
            "entropy_concentration": round(ent_conc, 3), "binom_pvalue": pval,
            "top_category": category, "tautological_pool": model in TAUTOLOGICAL_MODELS,
        })
    res = pd.DataFrame(rows)
    if len(res):
        res["fdr_qvalue"] = stats.false_discovery_control(res["binom_pvalue"], method="bh")
        res["significant_fdr05"] = res["fdr_qvalue"] < 0.05
        res["underpowered_n_lt6"] = res["n_genes"] < 6  # 精确二项检验在n<6时几乎不可能显著,单独标注
    return res


def main():
    df = pd.read_csv(COMBINED_PATH)

    print("=" * 70)
    print("主分析(全部候选细胞类型)")
    print("=" * 70)
    main_res = analyze(df, exclude_rare=False)
    print(main_res[["model", "n_genes", "K_celltypes", "top_celltype", "max_prop",
                     "entropy_concentration", "binom_pvalue", "fdr_qvalue", "significant_fdr05",
                     "underpowered_n_lt6", "top_category"]].to_string(index=False))
    main_res.to_csv(DATA_DIR + "phase3_meta_analysis_main.csv", index=False)

    print(f"\n显著模型数(FDR<0.05): {main_res['significant_fdr05'].sum()}/{len(main_res)}")
    print(f"检验功效不足(n_genes<6): {main_res['underpowered_n_lt6'].sum()}/{len(main_res)}")
    print("\n众数细胞类型三大类分布(全部模型,不论是否显著):")
    print(main_res["top_category"].value_counts())
    print("\n众数细胞类型三大类分布(仅FDR显著模型):")
    print(main_res[main_res["significant_fdr05"]]["top_category"].value_counts())

    sig = main_res[main_res["significant_fdr05"]]
    sig_taut = sig[sig["tautological_pool"]]
    sig_genuine = sig[~sig["tautological_pool"]]
    if len(sig_taut):
        print(f"\n其中{len(sig_taut)}个模型候选池预先限定为特定TME细胞marker基因,归因结果同义重复"
              f"(不是命名-机制脱节证据),已从下面H1a/H1b计数中剔除,单独列出:")
        print(sig_taut[["model", "top_celltype"]].to_string(index=False))
    h1a = sig_genuine[sig_genuine["top_category"].isin(["Stromal_Vascular", "Immune_Hematopoietic"])]
    h1b = sig_genuine[sig_genuine["top_category"] == "Malignant"]
    print(f"\n剔除同义重复模型后,H1a(TME细胞集中,真实脱节证据)={len(h1a)},"
          f"H1b(肿瘤细胞自身集中)={len(h1b)},其余={len(sig_genuine)-len(h1a)-len(h1b)}")

    print("\n" + "=" * 70)
    print("敏感性分析(排除参考细胞数<20的候选类型)")
    print("=" * 70)
    sens_res = analyze(df, exclude_rare=True)
    print(sens_res[["model", "n_genes", "K_celltypes", "top_celltype", "max_prop",
                     "binom_pvalue", "fdr_qvalue", "significant_fdr05"]].to_string(index=False))
    sens_res.to_csv(DATA_DIR + "phase3_meta_analysis_sensitivity.csv", index=False)

    merged = main_res[["model", "top_celltype", "significant_fdr05"]].merge(
        sens_res[["model", "top_celltype", "significant_fdr05"]], on="model",
        suffixes=("_main", "_sensitivity"))
    flipped = merged[
        (merged["top_celltype_main"] != merged["top_celltype_sensitivity"])
        | (merged["significant_fdr05_main"] != merged["significant_fdr05_sensitivity"])
    ]
    print(f"\n主分析与敏感性分析结果不一致的模型数: {len(flipped)}/{len(merged)}")
    if len(flipped):
        print(flipped.to_string(index=False))

    print(f"\n写入 {DATA_DIR}phase3_meta_analysis_main.csv 和 phase3_meta_analysis_sensitivity.csv")

    # 新增敏感性分析:排除候选细胞类型覆盖已知不完整的数据集(NSCLC_GSE131907只标注3种
    # 免疫细胞大类、BLCA_GSE130001推测只收录基质compartment),检验主结论是否依赖这两个
    # 覆盖有缺陷的数据集。
    print("\n" + "=" * 70)
    print("敏感性分析2(排除候选细胞类型覆盖已知不完整的数据集:NSCLC_GSE131907、BLCA_GSE130001)")
    print("=" * 70)
    INCOMPLETE_DATASETS = {"NSCLC_GSE131907", "BLCA_GSE130001"}
    excluded_models = {m for m, d in MODEL_TO_DATASET.items() if d in INCOMPLETE_DATASETS}
    df_cov = df[~df["model"].isin(excluded_models)]
    cov_res = analyze(df_cov, exclude_rare=False)
    cov_res.to_csv(DATA_DIR + "phase3_meta_analysis_coverage_sensitivity.csv", index=False)
    print(f"排除{len(excluded_models)}个模型(来自{INCOMPLETE_DATASETS}),剩余{len(cov_res)}个模型")
    print(f"显著模型数(FDR<0.05): {cov_res['significant_fdr05'].sum()}/{len(cov_res)}"
          f"  (主分析为 {main_res['significant_fdr05'].sum()}/{len(main_res)})")
    cov_sig = cov_res[cov_res["significant_fdr05"]]
    cov_sig_genuine = cov_sig[~cov_sig["tautological_pool"]]
    cov_h1a = cov_sig_genuine[cov_sig_genuine["top_category"].isin(["Stromal_Vascular", "Immune_Hematopoietic"])]
    cov_h1b = cov_sig_genuine[cov_sig_genuine["top_category"] == "Malignant"]
    print(f"剔除同义重复后,H1a={len(cov_h1a)}, H1b={len(cov_h1b)}  (主分析为 H1a={len(h1a)}, H1b={len(h1b)})")
    print(f"写入 {DATA_DIR}phase3_meta_analysis_coverage_sensitivity.csv")


if __name__ == "__main__":
    main()
