# -*- coding: utf-8 -*-
"""
v2 批量处置脚本 —— 对应 review_list_full.md「五、处置执行计划」

对 data/processed/knowledge_base_clean.jsonl（631 条）执行人工审核后的批量修正：
  - B/C 类 8 条：schema 化重建（方向/主客体 + relation 粒度）→ 同时收集为 Module 1 hard negative
  - D 类 6 条：降级 / 限定范围（消除 overclaim）
  - E 类 8 条：改变 relation 类型，与证据层级对齐（#375/#381 删除）
  - F 类 2 条：删除（hallucinated biological relation）
  - G 类 4 条：chain 修正 + 降低 design 强度（#57/#425 保留，#447/#449 删除并记为 design-rule 反例）
  - A 类 3 条：补 chain 反推/桥接句（#286/#272/#418），其余保留不动

输出：
  - data/processed/knowledge_base_clean.jsonl   修正后知识库（删除 6 条 → 625 条）
  - data/processed/hard_negative_module1.jsonl  B 类 5 条（原始错误 claim + 修正对照，供 Module 1 评测）
  - data/processed/design_rule_negatives.jsonl  G 类 2 条 design-rule 反例（agent treatment ≠ DNA feature）
  - data/processed/fix_v2_report.json           逐条处置明细

用法: python scripts/curate/fix_v2.py
"""
from __future__ import annotations

import json
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent.parent
KB_PATH = ROOT / "data/processed/knowledge_base_clean.jsonl"

# ---------------------------------------------------------------
# 待删除 idx（6 条）
#   F: #193 rs2294510(variant phenotype)、#367 NCOA4(ferroptosis)
#   E: #375 TSKU(无 enhancer 证据)、#381 ODC1(仅表达关联)
#   G: #447/#449 acetylating agents(agent treatment ≠ DNA sequence feature)
# ---------------------------------------------------------------
DELETES = {193, 367, 375, 381, 447, 449}

# ---------------------------------------------------------------
# 修正映射: idx -> {字段: 新值}
#   claim / entities(部分字段) / evidence / reasoning_chain(整链) /
#   experimental_context / confidence
# ---------------------------------------------------------------
FIXES: dict[int, dict] = {
    # ================= B 类：方向/主客体反 → schema 化重建（5 条）=================
    11: {
        "claim": "文献证据表明，在未指定细胞系的实验中，因子promoter对enhancer功能是必需的（required_for）；机制上，Promoter activity may be necessary but not sufficient for enhancer function.（来源：PMID 41617689）。",
        "entities": {
            "factor": "promoter",
            "factor_type": "other",
            "regulatory_element": "enhancer function",
            "context": "context=element_combination; promoter 活性对 enhancer 功能是必需但非充分条件 (necessary but not sufficient)",
        },
        "reasoning_chain": [
            "来源: PMID 41617689",
            "摘要原文证据: that promoter activity may be necessary but not sufficient for enhancer function.",
            "文献机制描述: Promoter activity may be necessary but not sufficient for enhancer function.",
            "因子类型判定: other；细胞系: not_specified；元件: enhancer function",
            "主客体修正: 证据主语是 promoter（对 enhancer 必需），修正原抽取的 enhancer→promoter 颠倒（error pattern: incorrect extraction / causal direction）",
        ],
    },
    236: {
        "claim": "文献证据表明，在SK-N-SH（神经母细胞瘤细胞）的实验中，因子DNMT1降低（decreases）FKBP52 promoter活性；机制上，Activated ERK enhanced DNA methyltransferase 1 (DNMT1)-mediated hypermethylation of FKBP52 promoter, reducing GR-mediated mitochondrial dysfunction and cell apoptosis.（来源：PMID 36810730）。",
        "entities": {
            "effect": "decreases",
            "regulatory_element": "FKBP52 promoter",
            "context": "context=hypermethylation; DNMT1 经 FKBP52 promoter 高甲基化抑制其活性 (PMID 36810730)",
        },
        "evidence": {"direction": "decreases_activity"},
        "reasoning_chain": [
            "来源: PMID 36810730",
            "摘要原文证据: The activated ERK then enhanced DNA methyltransferase 1 (DNMT1)-mediated hypermethylation of FKBP52 promoter, reducing GR-mediated mitochondrial dysfunction and cell apoptosis, the effects of which were reversed by knocking down DNMT1.",
            "文献机制描述: Activated ERK enhanced DNA methyltransferase 1 (DNMT1)-mediated hypermethylation of FKBP52 promoter, reducing GR-mediated mitochondrial dysfunction and cell apoptosis.",
            "因子类型判定: TF；细胞系: SK-N-SH；元件: FKBP52 promoter",
            "方向修正: 高甲基化→抑制（reducing），原抽取误写 increases，改为 decreases（error pattern: incorrect extraction / causal direction）",
        ],
    },
    250: {
        "claim": "文献证据表明，在SH-SY5Y（神经母细胞瘤细胞）的实验中，因子NFIB降低（decreases）CDON表达；机制上，Silencing NFIB leads to upregulation of CDON in SH-SY5Y cells.（来源：PMID 35151899）。",
        "entities": {
            "factor": "NFIB",
            "factor_type": "TF",
            "cell_line": "SH-SY5Y",
            "regulatory_element": "CDON expression",
            "effect": "decreases",
        },
        "evidence": {"direction": "decreases_activity"},
        "experimental_context": {"cell_type": "SH-SY5Y"},
        "reasoning_chain": [
            "来源: PMID 35151899",
            "摘要原文证据: silencing NFIB leads to upregulation of CDON in SH-SY5Y cells",
            "文献机制描述: Silencing NFIB leads to upregulation of CDON in SH-SY5Y cells.",
            "因子类型判定: TF；细胞系: SH-SY5Y；元件: CDON expression",
            "方向修正: 敲低→CDON 上调 ⇒ NFIB 存在→抑制，原抽取误写 required_for，改为 decreases（error pattern: incorrect extraction / causal direction）",
            "关系粒度修正: 证据在 expression 层，不写 promoter activity",
        ],
    },
    258: {
        "claim": "文献证据表明，在SK-N-SH（神经母细胞瘤细胞）的实验中，因子HDAC2降低（decreases）hSVCT2表达；机制上，Knockdown of HDAC2 significantly increased hSVCT2 functional expression.（来源：PMID 36096242）。",
        "entities": {
            "factor": "HDAC2",
            "factor_type": "TF",
            "regulatory_element": "hSVCT2 expression",
            "effect": "decreases",
        },
        "evidence": {"direction": "decreases_activity"},
        "reasoning_chain": [
            "来源: PMID 36096242",
            "摘要原文证据: Knockdown of HDAC2, a predominant isoform in neuronal systems, significantly increased hSVCT2 functional expression.",
            "文献机制描述: Knockdown of HDAC2 significantly increased hSVCT2 functional expression.",
            "因子类型判定: TF；细胞系: SK-N-SH；元件: hSVCT2 expression",
            "方向修正: 敲低→表达↑ ⇒ HDAC2 存在→抑制，原抽取误写 required_for，改为 decreases（error pattern: incorrect extraction / causal direction）",
            "关系粒度修正: 证据在 expression 层，不写 promoter activity",
        ],
    },
    399: {
        "claim": "文献证据表明，在未指定细胞系的实验中，因子OsDDE9调节（regulates）Ghd7 and Ehd1表达；机制上，OsDDE9 regulates the expression of Ghd7 and Ehd1 under both short-day and long-day conditions, which promotes heading date.（来源：PMID 37722508）。",
        "entities": {
            "factor": "OsDDE9",
            "factor_type": "other",
            "regulatory_element": "Ghd7 and Ehd1 expression",
            "effect": "regulates_expression",
        },
        "evidence": {"direction": "regulates_expression"},
        "reasoning_chain": [
            "来源: PMID 37722508",
            "摘要原文证据: OsDDE9 is a nuclear-localized protein expressed ubiquitously, which promotes heading date by regulating the expression of Ghd7 and Ehd1 under both short-day and long-day conditions.",
            "文献机制描述: OsDDE9 regulates the expression of Ghd7 and Ehd1 under both short-day and long-day conditions, which promotes heading date.",
            "因子类型判定: other；细胞系: not_specified；元件: Ghd7 and Ehd1 expression",
            "主客体修正: OsDDE9 是调节者，Ghd7/Ehd1 是被调节基因；原抽取把被调节基因当因子（error pattern: incorrect extraction / causal direction）",
            "关系粒度修正: 证据在 expression 层，不写 promoter 活性",
        ],
    },

    # ================= C 类：实验体系错 → 修正（3 条）=================
    224: {
        "claim": "文献证据表明，在SH-SY5Y（神经母细胞瘤细胞）的实验中，因子O-GlcNAc transferase (OGT)降低（decreases）ATF4表达；机制上，Knock-down of O-GlcNAc transferase (OGT) in SH-SY5Y increases ATF4 protein and mRNA expression, suggesting that OGT negatively regulates ATF4 expression.（来源：PMID 38192280）。",
        "entities": {
            "cell_line": "SH-SY5Y",
            "regulatory_element": "ATF4 expression",
            "effect": "decreases",
        },
        "evidence": {"direction": "decreases_activity"},
        "experimental_context": {"cell_type": "SH-SY5Y"},
        "reasoning_chain": [
            "来源: PMID 38192280",
            "摘要原文证据: knock-down of O-GlcNAc transferase (OGT), the enzyme that adds O-GlcNAc, in SH-SY5Y increases ATF4 protein and mRNA expression.",
            "文献机制描述: Knock-down of O-GlcNAc transferase (OGT) in SH-SY5Y increases ATF4 protein and mRNA expression, suggesting that OGT negatively regulates ATF4 expression.",
            "因子类型判定: other；细胞系: SH-SY5Y；元件: ATF4 expression",
            "实验体系修正: 证据细胞系为 SH-SY5Y，原抽取写 SK-N-SH（error pattern: incorrect extraction / context grounding）",
            "关系粒度修正: 证据只到 ATF4 expression（无 reporter/ChIP），不写 promoter activity",
        ],
    },
    226: {
        "claim": "文献证据表明，在小鼠脑组织（in vivo）的实验中，因子ATF4结合（binds）ATF5 promoter；机制上，ATF4 occupancy increases at the ATF5 promoter site in brains isolated from TMG treated mice suggesting that O-GlcNAc is regulating ATF4 targeted gene expression.（来源：PMID 38192280）。",
        "entities": {
            "cell_line": "not_specified",
            "system_type": "in vivo tissue",
            "organism": "mouse",
            "tissue": "brain",
            "regulatory_element": "ATF5 promoter",
            "effect": "binds",
        },
        "evidence": {"direction": "binds_target"},
        "experimental_context": {"cell_type": "brain tissue"},
        "reasoning_chain": [
            "来源: PMID 38192280",
            "摘要原文证据: ATF4 occupancy increases at the ATF5 promoter site in brains isolated from TMG treated mice suggesting that O-GlcNAc is regulating ATF4 targeted gene expression.",
            "文献机制描述: ATF4 occupancy increases at the ATF5 promoter site in brains isolated from TMG treated mice suggesting that O-GlcNAc is regulating ATF4 targeted gene expression.",
            "因子类型判定: TF；体系: in vivo tissue（mouse brain）；元件: ATF5 promoter",
            "实验体系修正: 证据为小鼠脑组织体内实验，原抽取写 SK-N-SH 细胞系（error pattern: incorrect extraction / context grounding）",
            "因果强度修正: occupancy 增加只能说明 ATF4 结合增加（binds），不能推出 required_for promoter activity（无 KO/reporter/perturbation）",
        ],
    },
    254: {
        "claim": "文献证据表明，在SH-SY5Y（神经母细胞瘤细胞）的实验中，因子YY1对promoter活性是必需的（required_for）；机制上，YY1 drives the SLC23A2 promoter activity, protein and mRNA expression levels were markedly upregulated in VPA-treated SH-SY5Y cells.（来源：PMID 36096242）。",
        "entities": {"cell_line": "SH-SY5Y"},
        "experimental_context": {"cell_type": "SH-SY5Y"},
        "reasoning_chain": [
            "来源: PMID 36096242",
            "摘要原文证据: Yin Yang-1 (YY1), a transcription factor that drives the SLC23A2 promoter activity, protein and mRNA expression levels were markedly upregulated in VPA-treated SH-SY5Y cells and mice brain.",
            "文献机制描述: YY1 drives the SLC23A2 promoter activity, protein and mRNA expression levels were markedly upregulated in VPA-treated SH-SY5Y cells.",
            "因子类型判定: TF；细胞系: SH-SY5Y；元件: promoter",
            "实验体系修正: 证据细胞系为 SH-SY5Y，原抽取写 SK-N-SH（error pattern: incorrect extraction / context grounding）",
        ],
    },

    # ================= D 类：过度泛化 → 降级/限定（6 条）=================
    14: {
        "claim": "文献证据表明，在未指定细胞系的实验中，因子EGR-1 binding motif对CSF1R enhancer激活是必需的（required_for）；机制上，A single EGR-1 binding motif dictates activation of CSF1R.（来源：PMID 41538316）。",
        "entities": {"regulatory_element": "CSF1R enhancer", "effect": "required_for"},
        "evidence": {"direction": "required_for_activity"},
        "reasoning_chain": [
            "来源: PMID 41538316",
            "摘要原文证据: we find that a single EGR-1 binding motif dictates activation of CSF1R.",
            "文献机制描述: A single EGR-1 binding motif dictates activation of CSF1R.",
            "因子类型判定: sequence_feature；细胞系: not_specified；元件: CSF1R enhancer",
            "范围限定: 单基因座证据，限定为 CSF1R（不泛化到所有 enhancer-promoter interaction；error pattern: unsupported specificity / overclaim）",
        ],
    },
    28: {
        "claim": "文献证据表明，在未指定细胞系的实验中，因子targeted activation of this regulatory element促进（promotes）AML分化；机制上，The targeted activation of this regulatory element promotes differentiation of these aggressive AMLs and reduces leukemia burden in vivo.（来源：PMID 40991835）。",
        "entities": {
            "regulatory_element": "AML differentiation",
            "effect": "promotes",
        },
        "evidence": {"direction": "promotes_phenotype"},
        "reasoning_chain": [
            "来源: PMID 40991835",
            "摘要原文证据: the targeted activation of this regulatory element promotes differentiation of these aggressive AMLs and reduces leukemia burden in vivo.",
            "文献机制描述: The targeted activation of this regulatory element promotes differentiation of these aggressive AMLs and reduces leukemia burden in vivo.",
            "因子类型判定: other；细胞系: not_specified；元件: AML differentiation",
            "干预≠活性读数: activation 是干预手段（如 CRISPRa），不是 measured enhancer activity；改写为表型描述（error pattern: unsupported specificity / overclaim）",
        ],
    },
    34: {
        "claim": "文献证据表明，在未指定细胞系的实验中，因子FOXP4增强（increases）CYP26B1和MYC转录；机制上，FOXP4 enhanced ESCC susceptibility and tumor growth by transcriptionally activating CYP26B1 and MYC.（来源：PMID 40228145）。",
        "entities": {"regulatory_element": "CYP26B1 and MYC transcription", "effect": "increases"},
        "evidence": {"direction": "increases_activity"},
        "reasoning_chain": [
            "来源: PMID 40228145",
            "摘要原文证据: FOXP4 enhanced ESCC susceptibility and tumor growth by transcriptionally activating CYP26B1 and MYC.",
            "文献机制描述: FOXP4 enhanced ESCC susceptibility and tumor growth by transcriptionally activating CYP26B1 and MYC.",
            "因子类型判定: TF；细胞系: not_specified；元件: CYP26B1 and MYC transcription",
            "关系粒度修正: 转录激活基因 ≠ 实验测 promoter activity；降级为 transcription 描述（error pattern: unsupported specificity / overclaim）",
        ],
    },
    43: {
        "claim": "文献证据表明，在未指定细胞系的实验中，遗传变异（QTLs）调节（modulates）SPI1、GABPB和STAT3的预测结合位点；机制上，Many detected QTLs influence the predicted binding of myeloid transcription factors, including SPI1, GABPB and STAT3.（来源：PMID 39870618）。",
        "entities": {
            "factor": "QTLs",
            "factor_type": "variant",
            "regulatory_element": "predicted binding sites of SPI1/GABPB/STAT3",
            "effect": "modulates",
        },
        "evidence": {"direction": "modulates_activity"},
        "reasoning_chain": [
            "来源: PMID 39870618",
            "摘要原文证据: many detected QTLs overlap disease susceptibility loci and influence the predicted binding of myeloid transcription factors, including SPI1, GABPB and STAT3.",
            "文献机制描述: Many detected QTLs influence the predicted binding of myeloid transcription factors, including SPI1, GABPB and STAT3.",
            "因子类型判定: variant；细胞系: not_specified；元件: SPI1/GABPB/STAT3 预测结合位点",
            "association≠mechanism: variant→predicted binding change，不能推出 TF→enhancer activity；降级为 variant 层描述（error pattern: unsupported specificity / overclaim）",
        ],
    },
    470: {
        "claim": "文献证据表明，在未指定细胞系的实验中，caQTL变异破坏（disrupts）肝脏表达转录因子的结合基序；机制上，caQTL variants frequently disrupt binding motifs of transcription factors expressed in liver.（来源：PMID 34038741）。",
        "entities": {
            "factor": "caQTL variants",
            "factor_type": "variant",
            "regulatory_element": "TF binding motifs",
            "effect": "disrupts",
        },
        "evidence": {"direction": "disrupts"},
        "reasoning_chain": [
            "来源: PMID 34038741",
            "摘要原文证据: The caQTL variants are enriched in liver tissue promoter and enhancer states and frequently disrupt binding motifs of transcription factors expressed in liver.",
            "文献机制描述: caQTL variants frequently disrupt binding motifs of transcription factors expressed in liver.",
            "因子类型判定: variant；细胞系: not_specified；元件: TF binding motifs",
            "范围限定: variant-level 证据（位点破坏），不能泛化为 motif 普遍调节 promoter；限定为 caQTL 变异破坏基序（error pattern: unsupported specificity / overclaim）",
        ],
    },
    479: {
        "claim": "文献证据表明，在未指定细胞系的实验中，因子NF-Y对MAFA promoter活性是必需的（required_for）；机制上，Two E-box elements and a CCAAT motif, which bind NeuroD1 and ubiquitous NF-Y transcription factors, respectively, were necessary for transcriptional activation of the MAFA promoter by CREB.（来源：PMID 39189982）。",
        "entities": {"regulatory_element": "MAFA promoter", "effect": "required_for"},
        "evidence": {"direction": "required_for_activity"},
        "reasoning_chain": [
            "来源: PMID 39189982",
            "摘要原文证据: Two E-box elements and a CCAAT motif, which bind NeuroD1 and ubiquitous NF-Y transcription factors, respectively, were necessary for transcriptional activation of the MAFA promoter by CREB.",
            "文献机制描述: Two E-box elements and a CCAAT motif, which bind NeuroD1 and ubiquitous NF-Y transcription factors, respectively, were necessary for transcriptional activation of the MAFA promoter by CREB.",
            "因子类型判定: other；细胞系: not_specified；元件: MAFA promoter",
            "范围限定: 单启动子证据（MAFA），不泛化到所有 promoter；motif 为 necessary（有 mutation 支持），required_for 限定于 MAFA（error pattern: unsupported specificity / overclaim）",
        ],
    },

    # ================= E 类：证据间接 → 改变 relation 类型（8 条，#375/#381 删除）=================
    54: {
        "claim": "文献证据表明，在未指定细胞系的实验中，因子LEF1降低（decreases）炎症反应相关基因表达；机制上，LEF1 knockdown enhanced inflammatory responses and ROS production in vitro, indicating that LEF1 is required for the repression of inflammatory responses and ROS levels.（来源：PMID 40918098）。",
        "entities": {
            "regulatory_element": "inflammatory response-related gene expression",
            "effect": "decreases",
        },
        "evidence": {"direction": "decreases_activity"},
        "reasoning_chain": [
            "来源: PMID 40918098",
            "摘要原文证据: LEF1 knockdown enhanced inflammatory responses and ROS production in vitro.",
            "文献机制描述: LEF1 knockdown enhanced inflammatory responses and ROS production in vitro, indicating that LEF1 is required for the repression of inflammatory responses and ROS levels.",
            "因子类型判定: TF；细胞系: not_specified；元件: inflammatory response-related gene expression",
            "relation 类型修正: 证据在基因表达/表型层（敲低→炎症↑），不是 promoter 必需；改为 gene regulation（不补 β-catenin/TCF 通路，原文无直接证据）（error pattern: evidence gap / relation 过强）",
        ],
    },
    160: {
        "claim": "文献证据表明，在K562（红系前体细胞）的实验中，因子N1IC and Ets-1降低（decreases）DNMT3B表达；机制上，N1IC and Ets-1 suppressed the DNA methyltransferase 3B (DNMT3B) level in K562 cells.（来源：PMID 28220825）。",
        "entities": {
            "factor": "N1IC and Ets-1",
            "factor_type": "TF",
            "regulatory_element": "DNMT3B expression",
            "effect": "decreases",
        },
        "evidence": {"direction": "decreases_activity"},
        "reasoning_chain": [
            "来源: PMID 28220825",
            "摘要原文证据: N1IC and Ets-1 suppressed the DNA methyltransferase 3B (DNMT3B) level in K562 cells.",
            "文献机制描述: N1IC and Ets-1 suppressed the DNA methyltransferase 3B (DNMT3B) level in K562 cells.",
            "因子类型判定: TF；细胞系: K562；元件: DNMT3B expression",
            "relation 类型修正: 证据只到 DNMT3B expression 水平；原 claim 的 DNMT3B→promoter 无证据，且不加 DNMT3B→methylation→promoter 层级（error pattern: evidence gap / relation 过强）",
        ],
    },
    170: {
        "claim": "文献证据表明，在β-YAC骨髓细胞的实验中，GATA位点突变降低（decreases）OGT和OGA的promoter相互作用；机制上，Mutation of the GATA site to GAGA significantly reduces OGT and OGA promoter interactions in β-globin locus yeast artificial chromosome (β-YAC) bone marrow cells.（来源：PMID 27231347）。",
        "entities": {
            "factor": "GATA motif",
            "factor_type": "sequence_feature",
            "regulatory_element": "OGT/OGA promoter interaction",
            "effect": "affects_interaction",
        },
        "evidence": {"direction": "affects_interaction"},
        "reasoning_chain": [
            "来源: PMID 27231347",
            "摘要原文证据: mutation of the GATA site to GAGA significantly reduces OGT and OGA promoter interactions in β-globin locus yeast artificial chromosome (β-YAC) bone marrow cells.",
            "文献机制描述: Mutation of the GATA site to GAGA significantly reduces OGT and OGA promoter interactions in β-globin locus yeast artificial chromosome (β-YAC) bone marrow cells.",
            "因子类型判定: sequence_feature；体系: β-YAC bone marrow cells；元件: OGT/OGA promoter interaction",
            "relation 类型修正: 证据测的是 promoter interaction（非 activity）；binding→interaction→activity 三级不跨越（error pattern: evidence gap / relation 过强）",
        ],
    },
    199: {
        "claim": "文献证据表明，在HepG2（肝癌细胞）的实验中，因子ALKBH5增强（increases）BCP-RNA和宿主基因转录；机制上，Silencing the demethylase reduced the level of BCP-RNAs and host gene transcripts.（来源：PMID 38227578）。",
        "entities": {
            "regulatory_element": "BCP-RNA / host gene transcription",
            "effect": "increases",
        },
        "evidence": {"direction": "increases_activity"},
        "reasoning_chain": [
            "来源: PMID 38227578",
            "摘要原文证据: silencing the demethylase reduced the level of BCP-RNAs and host gene (CA9, NDRG1, VEGFA, BNIP3, FUT11, GAP and P4HA1) transcripts",
            "文献机制描述: Silencing the demethylase reduced the level of BCP-RNAs and host gene transcripts.",
            "因子类型判定: TF；细胞系: HepG2；元件: BCP-RNA / host gene transcription",
            "relation 类型修正: 证据测转录本水平，改为 transcription regulation；不提去甲基化机制（来自 ALKBH5 身份而非该实验）（error pattern: evidence gap / relation 过强）",
        ],
    },
    292: {
        "claim": "文献证据表明，在未指定细胞系的实验中，promoter高甲基化与基因上调共现（co_occurs_with）；机制上，Hypermethylation at promoters of up-regulated genes is observed in C/G fusion-positive acute megakaryoblastic leukemia.（来源：PMID 41141393、PMID 28669403）。",
        "entities": {
            "factor": "DNA methylation",
            "factor_type": "epigenetic",
            "regulatory_element": "gene upregulation",
            "effect": "co_occurs_with",
        },
        "evidence": {"direction": "co_occurs_with"},
        "reasoning_chain": [
            "来源: PMID 41141393、PMID 28669403",
            "摘要原文证据: This multi-omics analysis reveals a distinct hypermethylation pattern at promoters of up-regulated genes in C/G",
            "文献机制描述: Hypermethylation at promoters of up-regulated genes is observed in C/G fusion-positive acute megakaryoblastic leukemia.",
            "因子类型判定: epigenetic；细胞系: not_specified；元件: gene upregulation",
            "relation 类型修正: 观察为 methylation↑ 与 expression↑ 共现（与经典负相关相反），去掉 direction 改共现描述，不强行给方向（error pattern: evidence gap / relation 过强）",
        ],
    },
    323: {
        "claim": "文献证据表明，在未指定细胞系的实验中，chromatin accessibility与promoter区域甲基化水平相关（correlates）；机制上，ATAC-seq peaks with high chromatin accessibility located in both promoter (≤ 2 kb from TSS) and distal (> 2 kb from TSS) regions corresponded to low methylation levels.（来源：PMID 38281519）。",
        "entities": {
            "regulatory_element": "promoter methylation level",
            "effect": "correlates",
        },
        "evidence": {"direction": "correlates_with"},
        "reasoning_chain": [
            "来源: PMID 38281519",
            "摘要原文证据: ATAC-seq peaks with high chromatin accessibility located in both promoter (≤ 2 kb from TSS) and distal (> 2 kb from TSS) regions corresponded to low methylation levels.",
            "文献机制描述: ATAC-seq peaks with high chromatin accessibility located in both promoter (≤ 2 kb from TSS) and distal (> 2 kb from TSS) regions corresponded to low methylation levels.",
            "因子类型判定: epigenetic；细胞系: not_specified；元件: promoter methylation level",
            "relation 类型修正: ATAC-seq 观察性证据（可及性↔低甲基化），correlation≠causation；改为 correlates（error pattern: evidence gap / relation 过强）",
        ],
    },

    # ================= G 类：推理链缺陷 → chain 修正 + 降 design 强度（4 条，#447/#449 删除）=================
    57: {
        "claim": "设计建议（文献增强）：在设计未指定细胞系的promoter序列时，应审慎考虑DNA methylation相关特征，因其与高活性呈负相关，需实验验证后再定取舍（来源：PMID 40918098）。",
        "entities": {
            "effect": "decreases",
            "rule_confidence": "low",
            "requires_validation": True,
        },
        "evidence": {"direction": "decreases_activity"},
        "reasoning_chain": [
            "由文献关联推导: PMID 40918098",
            "原始发现: DNA methylation 与 promoter 高活性呈负相关（decreases）",
            "设计建议: 审慎考虑相关特征（低置信，需实验验证）",
        ],
        "confidence": 0.7,
    },
    425: {
        "claim": "设计建议（文献增强）：在设计K562（红系前体细胞）的promoter序列时，应审慎考虑H3 acetylation at lysine 4 and 9 (K4 and K9)相关特征，因其与高活性呈负相关，需实验验证后再定取舍（来源：PMID 28216155）。",
        "entities": {
            "effect": "decreases",
            "rule_confidence": "low",
            "requires_validation": True,
        },
        "evidence": {"direction": "decreases_activity"},
        "reasoning_chain": [
            "由文献关联推导: PMID 28216155",
            "原始发现: H3 acetylation (K4/K9) 与 promoter 高活性呈负相关（decreases）",
            "设计建议: 审慎考虑相关特征（低置信，需实验验证）",
        ],
        "confidence": 0.7,
    },

    # ================= A 类：补 chain 反推/桥接句（3 条，其余保留）=================
    272: {
        "reasoning_chain": [
            "来源: PMID 41456103",
            "摘要原文证据: Sp3 binding regulation likely depended on chromatin accessibility.",
            "文献机制描述: Sp3 transcription factor binds to the GC box 4 most involved in coronin-1 expression.",
            "因子类型判定: other；细胞系: not_specified；元件: promoter",
            "反推: Sp3 结合对 coronin-1 表达必需的 GC box 4 → 支持 required_for（loss-of-function 反推）",
        ],
    },
    286: {
        "reasoning_chain": [
            "来源: PMID 41708589",
            "摘要原文证据: ERRα depletion leads to further activation of STING gene transcription and TBK1-IRF3 pathway, accompanied by increased type I interferon (IFN) and IFN-stimulated gene (ISG) expression.",
            "文献机制描述: ERRα depletion leads to increased STING gene transcription.",
            "因子类型判定: other；细胞系: not_specified；元件: sting gene transcription",
            "反推: ERRα 耗竭→STING 转录激活 ⇒ ERRα 存在→抑制 STING 转录（loss-of-function 反推）",
        ],
    },
    418: {
        "reasoning_chain": [
            "来源: PMID 33637692",
            "摘要原文证据: DNA methylation inhibited MYB expression.",
            "文献机制描述: DNA methylation inhibited MYB expression.",
            "因子类型判定: epigenetic；细胞系: K562；元件: enhancer",
            "桥接: 甲基化→抑制 MYB 表达；MYB 为红系关键 TF → 降低其靶 enhancer 活性",
        ],
    },
}


def deep_merge(base: dict, patch: dict) -> dict:
    """dict 递归合并（patch 覆盖 base 的同名键；嵌套 dict 递归）。"""
    out = dict(base)
    for k, v in patch.items():
        if isinstance(v, dict) and isinstance(out.get(k), dict):
            out[k] = deep_merge(out[k], v)
        else:
            out[k] = v
    return out


def main() -> None:
    rows = []
    with open(KB_PATH, encoding="utf-8") as f:
        for line in f:
            line = line.strip()
            if line:
                rows.append(json.loads(line))

    n_before = len(rows)
    report: list[dict] = []
    hard_negatives: list[dict] = []
    design_negatives: list[dict] = []

    # ---------- 第一步：按原 idx 应用修正（删除前定位） ----------
    for i, r in enumerate(rows):
        # G 类 #447/#449 → design-rule 反例（先收集，即使随后删除）
        if i in (447, 449):
            design_negatives.append({
                "idx": i,
                "error_pattern": "reasoning chain risk (design-rule overreach)",
                "claim": r.get("claim"),
                "reason": "agent treatment（药物干预）不是 DNA sequence design feature；若产出为序列设计规则，此条不成立",
                "pmids": r.get("entities", {}).get("pmids"),
            })
        if i in DELETES:
            continue
        fix = FIXES.get(i)
        if fix is None:
            continue
        merged = deep_merge(r, fix)
        rows[i] = merged
        report.append({
            "idx": i,
            "claim_type": merged.get("claim_type"),
            "new_claim": merged.get("claim"),
            "changed_fields": list(fix.keys()),
        })
        # B 类 5 条 → Module 1 hard negative
        if i in (11, 236, 250, 258, 399):
            hard_negatives.append({
                "idx": i,
                "error_pattern": "incorrect extraction (causal direction)",
                "original_claim": r.get("claim"),
                "corrected_claim": merged.get("claim"),
                "corrected_entities": merged.get("entities"),
                "pmids": merged.get("entities", {}).get("pmids"),
            })

    # ---------- 第二步：删除（从大到小，避免索引漂移） ----------
    deleted_records = []
    for i in sorted(DELETES, reverse=True):
        deleted_records.append({"idx": i, "claim": rows[i].get("claim")})
        del rows[i]

    # ---------- 写回 ----------
    with open(KB_PATH, "w", encoding="utf-8") as f:
        for r in rows:
            f.write(json.dumps(r, ensure_ascii=False) + "\n")

    hn_path = ROOT / "data/processed/hard_negative_module1.jsonl"
    with open(hn_path, "w", encoding="utf-8") as f:
        for r in hard_negatives:
            f.write(json.dumps(r, ensure_ascii=False) + "\n")

    dn_path = ROOT / "data/processed/design_rule_negatives.jsonl"
    with open(dn_path, "w", encoding="utf-8") as f:
        for r in design_negatives:
            f.write(json.dumps(r, ensure_ascii=False) + "\n")

    rpt_path = ROOT / "data/processed/fix_v2_report.json"
    with open(rpt_path, "w", encoding="utf-8") as f:
        json.dump({
            "input_rows": n_before,
            "output_rows": len(rows),
            "deleted": len(DELETES),
            "fixed": len(report),
            "hard_negatives": len(hard_negatives),
            "design_rule_negatives": len(design_negatives),
            "deleted_records": deleted_records,
            "fix_details": report,
        }, f, ensure_ascii=False, indent=1)

    print(f"[fix_v2] 输入 {n_before} 条 → 输出 {len(rows)} 条（删除 {len(DELETES)}，修正 {len(report)}）")
    print(f"[fix_v2] hard_negative_module1.jsonl: {len(hard_negatives)} 条")
    print(f"[fix_v2] design_rule_negatives.jsonl: {len(design_negatives)} 条")
    print(f"[fix_v2] 报告: {rpt_path}")


if __name__ == "__main__":
    main()
