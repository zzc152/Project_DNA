"""[DEPRECATED] 方案A: 分层过滤知识库，产出"干净"的行动型知识库。

⚠️ 已废弃：本脚本产出旧版三元组格式（knowledge_base_actionable.jsonl）。
项目现统一使用 claim 版本知识库（knowledge_base_stat.jsonl，8 键 schema），
由 build_stat_claims.py + build_literature_claims.py 构建。
旧文件已移至 data/backup/。请勿再使用本脚本。

问题诊断（knowledge_base_actionable.jsonl, 3946 条）:
  1. 自指噪音 head==tail (140 条)
  2. 疾病关联 associated_with（~448 条，对序列设计无用）
  3. 泛化/描述性短语（85 条）
  4. binds_motif 但 tail 不是基序（475 条，如 AHR->FOXA1 是 TF-TF 互作）
  5. regulates 的 head 是信号蛋白/非 TF（4EBP1、ADRBK1 等被模型误放入 tf 字段）

关键做法:
  knowledge_base_full.jsonl 只保留 head/relation/tail，丢失了抽取时的
  tf/gene/motif/disease 分类，无法做精确过滤。因此本脚本直接从
  raw_extractions*.jsonl（保留分类）重建文献知识库，并结合:
    - flagged v4 幻觉清单（排除不在摘要中的实体）
    - 黑名单泛化词（复用 review_extractions.BLACKLIST_WORDS）
    - JASPAR TF 白名单（regulates 的 head 必须是真 TF）
    - 基序特征词（binds_motif 的 tail 必须是真基序）
    - 自指/疾病关联直接丢弃

输出分层（knowledge_type 字段）:
  - motif_activity   : MPRA 统计基序-活性效应（34 条，直接保留）
  - gc_activity      : MPRA 统计 GC-活性偏好（6 条，直接保留）
  - sequence_actionable : 文献 binds_motif（TF->真基序），可直接转序列约束
  - gene_regulation  : 文献 regulates（真TF->基因），背景调控知识
  其余全部丢弃。

用法（项目根目录）:
    python scripts/filter_knowledge_base.py
"""

from __future__ import annotations

import argparse
import json
import logging
import re
import sys
from collections import Counter
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s | %(levelname)s | %(message)s",
    datefmt="%H:%M:%S",
)
logger = logging.getLogger("filter_knowledge_base")

SOURCE_LIT = "PubMed literature extraction (Qwen2.5-7B-Instruct)"

# ---------------------------------------------------------------------------
# 黑名单泛化词（与 review_extractions.BLACKLIST_WORDS 一致）
# ---------------------------------------------------------------------------
BLACKLIST_WORDS = [
    "master regulators", "transcription factors", "tfs", "regulators",
    "activators", "silencers", "bhlh", "homeodomain", "hd",
    "transfer rnas", "dna sequences", "binding sites", "coordinator",
    "composite motifs", "motifs", "disease", "autoimmune diseases",
    "disease risk alleles", "tn5", "basic helix-loop-helix",
]
_SPECIAL_EXACT = {"DISEASE", "DISEASES"}


def is_blacklisted(entity: str) -> bool:
    """实体（已大写）是否命中泛化词黑名单（"disease" 仅精确匹配）。"""
    for w in BLACKLIST_WORDS:
        wu = w.upper()
        if wu in _SPECIAL_EXACT:
            if entity == wu:
                return True
        elif wu in entity:
            return True
    return False


# ---------------------------------------------------------------------------
# JASPAR TF 白名单: 用于过滤 regulates 的 head（必须是真 TF）
# ---------------------------------------------------------------------------
_JASPAR_MOTIF_TFS = {
    # 从 motifs_*.meme 提取的 19 个基序对应 TF
    "RUNX1", "GATA1", "GATA2", "GABPA", "SP1", "SPI1", "TAL1", "TCF3",
    "YY1", "ETS1", "MYB", "CTCF", "NFE2", "NFE2L2", "ELF1", "FLI1",
    "KLF1", "MAF", "NRF1",
}
_COMMON_TFS = {
    # 常见/高频转录因子（大写，含连字符写法）
    "TP53", "P53", "MYC", "FOXA1", "FOXA2", "FOXO1", "FOXO3", "FOXP3",
    "NFY", "NF-Y", "NFYA", "NFYB", "HNF4A", "HNF4", "HNF1", "HNF1A",
    "HNF1B", "HNF3", "CEBPA", "CEBPB", "CEBP", "TWIST1", "NGN2",
    "STAT1", "STAT3", "STAT5", "STAT6", "AP-1", "AP1", "FOS", "C-FOS",
    "JUN", "C-JUN", "FOSB", "JUNB", "GATA-1", "GATA-3", "GATA-4",
    "PU.1", "SP3", "SP2", "ONECUT1", "E2F", "E2F1", "E2F-1", "E2F4",
    "OCT-1", "OCT1", "OCT-2", "OCT2", "POU2F1", "POU2F2", "NF-KB",
    "NF-KAPPA B", "NFKB1", "NFKB2", "RELA", "NF-E2", "COUP-TF",
    "PIT-1", "POU1F1", "IRF-1", "IRF1", "IRF4", "ATF-3", "ATF3",
    "ATF4", "ATF", "MITF", "BATF", "SOX2", "SOX9", "SOX10", "SRY",
    "HOXA9", "HOXB4", "MEIS1", "PBX1", "ERG", "FLI-1", "TAL1::TCF3",
    "GATA1::TAL1", "MAF::NFE2", "HIF1A", "ARNT", "USF1", "USF2",
    "MAX", "MAD", "MXI1", "TFAP2A", "AP-2", "AP2", "TFAP2",
    "PITX1", "NKX2-1", "NKX2.1", "TTF-1", "PAX5", "PAX6", "PAX2",
    "RUNX2", "RUNX3", "AML-1", "CBFB", "ETV4", "ETV1", "PEA3",
    "SPIB", "ELF4", "FOSL1", "FRA-1", "JUNB", "JUND", "CREB1", "CREB",
    "ATF2", "XBP1", "NFATC1", "NFATC2", "GATA6", "SNAI1", "SNAIL",
    "SNAI2", "SLUG", "TWIST2", "ZEB1", "ZEB2", "TCF7", "TCF7L2",
    "LEF1", "TCF4", "TCF12", "HES1", "HEY1", "NOTCH1", "RBPJ",
    "CBF1", "SMAD3", "SMAD4", "SMAD1", "SMAD2", "SMAD5",
    "GLI1", "GLI2", "GLI3", "ZIC1", "ZIC2", "NR3C1", "GR",
    "GLUCOCORTICOID RECEPTOR", "ESR1", "AR", "ANDROGEN RECEPTOR",
    "RARA", "RAR", "RXRA", "VDR", "PPARG", "PPARA", "NR1H4",
    "SREBP1", "SREBF1", "SREBP2", "CHOP", "DDIT3", "C/EBP",
    "C/EBP ALPHA", "C/EBP BETA", "NRF2", "BACH1", "BACH2",
    "KAPPA B", "OCT", "OCTAMER", "POU", "E-BOX FACTORS",
    "BASIC HELIX-LOOP-HELIX PROTEINS", "HLH", "MEF2", "MEF2A",
    "MEF2C", "MEF2D", "NFAT", "ETS", "GABP", "GABPB1",
    "ELF2", "ETS2", "ERF", "ETV3", "ETV5", "SPDEF",
    "KLF2", "KLF4", "KLF5", "KLF6", "KLF9", "EGR1", "EGR2",
    "WT1", "ZNF263", "ZNF143", "ZBTB33", "ZFP", "YY1AP",
    "REST", "NRSF", "MECP2", "CTCFL", "BORIS", "THAP1",
    "NRF1", "GABP", "GABPA", "SFPI1", "ETS1", "FLI1", "ERG",
    "FELV1", "PLAG1", "PLAGL1", "TFIID", "TBP", "TFIIB",
    "TFIIA", "TFIIH", "RNA POLYMERASE II",
    # 高频真 TF 补充（从 tf 字段频率探测获得）
    "GATA3", "OCT4", "POU5F1", "CRX", "MYOD", "MYOD1", "MYF5",
    "NFE2L1", "OTF-1", "OTF1", "HNF-4", "HNF4A", "NR1I3", "SP4",
    "RFX1", "NANOG", "THAP11", "NF-1", "NF1", "NFKB", "NFKB1",
    "E47", "TCF3", "NFI", "NFIA", "NFIB", "NFIC", "NFIX", "USF",
    "USF1", "USF2", "EBF", "EBF1", "SRF", "E2A", "TCF3", "FOXM1",
    "STAT4", "TWIST", "TWIST1", "P63", "TP63", "KLF", "FOXA",
    "GATA", "GATA2", "SNAIL", "SLUG", "SNAI1", "SNAI2", "ZEB1",
    "ZEB2", "SOX", "SOX2", "SOX9", "SOX10", "HOX", "HOXA9",
    "HOXB4", "PBX1", "MEIS1", "MEIS", "RBPJ", "CBF1", "HES1",
    "HEY1", "HEY2", "NRF1", "NRF2", "NFE2", "NF-E2", "MAF",
    "MAFB", "MAFK", "MAFG", "FOS", "FOSB", "FOSL1", "FOSL2",
    "JUN", "JUNB", "JUND", "ATF", "ATF1", "ATF2", "ATF3", "ATF4",
    "ATF6", "CREB", "CREB1", "CREM", "CEBP", "CEBPA", "CEBPB",
    "CEBPD", "CEBPG", "SP1", "SP2", "SP3", "SP4", "SP5",
    "E2F", "E2F1", "E2F2", "E2F3", "E2F4", "E2F5", "DP1",
    "TFAP2", "TFAP2A", "TFAP2B", "TFAP2C", "AP-2", "AP2",
    "IRF", "IRF1", "IRF2", "IRF3", "IRF4", "IRF7", "IRF8",
    "STAT", "STAT1", "STAT2", "STAT3", "STAT4", "STAT5",
    "STAT5A", "STAT5B", "STAT6", "SMAD", "SMAD1", "SMAD2",
    "SMAD3", "SMAD4", "SMAD5", "SMAD7", "RUNX", "RUNX1",
    "RUNX2", "RUNX3", "AML-1", "AML1", "CBFB", "CBF",
    "GABPA", "GABPB", "GABPB1", "GABPB2", "ELF", "ELF1",
    "ELF2", "ELF3", "ELF4", "ETS", "ETS1", "ETS2", "ERG",
    "FLI1", "FLI-1", "ERF", "ETV1", "ETV3", "ETV4", "ETV5",
    "ETV6", "SPIB", "SPIC", "PU.1", "SPI1", "SFPI1", "FBI-1",
    "ZBTB7A", "P300", "EP300", "CBP", "CREBBP", "RB", "RB1",
    "RBL1", "RBL2", "P107", "P130", "CAR", "NR1I3", "CONSTITUTIVE "
    "ANDROSTANE RECEPTOR", "PPARA", "PPARD", "PPARG", "NR1H3",
    "NR1H4", "NR1D1", "REV-ERBA", "RORA", "RORC", "ESR2",
    "NR3C2", "MR", "GATA4", "GATA5", "GATA6", "GATA-2",
    "GATA-3", "GATA-4", "GATA-6", "PIT-1", "PITX1", "PITX2",
    "NKX2-1", "NKX2-2", "NKX2-5", "TTF-1", "PAX", "PAX2",
    "PAX3", "PAX5", "PAX6", "PAX7", "PAX8", "FOX", "FOXA1",
    "FOXA2", "FOXA3", "FOXC1", "FOXD3", "FOXE1", "FOXF1",
    "FOXH1", "FOXK1", "FOXL2", "FOXM1", "FOXO1", "FOXO3",
    "FOXO4", "FOXP1", "FOXP2", "FOXP3", "HNF", "HNF1",
    "HNF1A", "HNF1B", "HNF3", "HNF3A", "HNF3B", "HNF4",
    "HNF4A", "HNF6", "ONECUT1", "OC-2", "OC-3", "TALE",
    "HOXA", "HOXB", "HOXC", "HOXD", "MSX1", "MSX2", "DLX1",
    "DLX2", "DLX3", "DLX4", "DLX5", "LHX1", "LHX2", "LHX3",
    "ISL1", "ISL2", "POU", "POU2F1", "POU2F2", "POU3F2",
    "POU4F1", "POU5F1", "BRN2", "BRN3", "BRN4", "OCT-1",
    "OCT-2", "OCT-4", "OCT1", "OCT2", "OCT4", "NANOG",
    "SOX2", "KLF4", "MYC", "MYCN", "MYCL", "MAX", "MAD",
    "MXI1", "MNT", "USF", "TCF3", "E12", "E47", "E2A",
    "HEB", "TCF12", "TCF4", "TCF7", "TCF7L1", "TCF7L2",
    "LEF1", "HIC1", "ZIC1", "ZIC2", "ZIC3", "GLI1", "GLI2",
    "GLI3", "SAL1", "SALL1", "SALL4", "WT1", "EGR1", "EGR2",
    "EGR3", "EGR4", "NAB1", "NAB2", "KLF1", "KLF2", "KLF3",
    "KLF4", "KLF5", "KLF6", "KLF7", "KLF8", "KLF9", "KLF10",
    "KLF11", "KLF12", "KLF13", "KLF15", "SP", "SP1", "SP2",
    "SP3", "SP4", "SP5", "SP6", "SP7", "SP8", "ZBTB",
    "ZBTB7A", "ZBTB16", "PLZF", "ZFP36", "TTP", "GATA",
    "TAL1", "TAL2", "SCL", "LYL1", "LMO1", "LMO2", "LDB1",
    "NFE2", "BACH1", "BACH2", "MAF", "MAFA", "MAFB", "MAFG",
    "MAFK", "NRF1", "NRF2", "NFE2L1", "NFE2L2", "NRF3",
    "THAP11", "REST", "NRSF", "RCOR1", "COREST", "MECP2",
    "MBD1", "MBD2", "MBD3", "MBD4", "CTCF", "CTCFL", "BORIS",
    "YY1", "YY2", "ZNF143", "ZNF263", "ZNF281", "ZNF148",
    "ZNF281", "THAP1", "SMARCA4", "BRG1", "SMARCB1", "SNF5",
    "ARID1A", "BAF", "SWI/SNF", "CHD1", "CHD2", "CHD4",
    "MI-2", "HMG", "HMGA1", "HMGA2", "HMGB1", "HMGB2", "SOX9",
    "SRY", "TBP", "TFIID", "TFIIB", "TFIIA", "TFIIE",
    "TFIIF", "TFIIH", "SPT", "TAF1", "TAF", "BRF", "SL1",
    "TFIIIB", "TFIIIC", "RNA POLYMERASE I", "RNA POLYMERASE III",    # 高频真 TF 补充（第二轮探测）
    "C/EBPB", "C/EBPEPSILON", "IKAROS", "IKZF1", "MAZ", "HSF1",
    "CDX1", "CDX2", "WRKY", "BCL6", "BCL6B", "CREL", "P50",
    "NFKB1", "NOTCH", "NOTCH1", "ICN1", "ICN", "BRD4", "P73",
    "TP73", "TAP73", "E2A-PBX1", "PBX1", "HNF-1", "HNF1",
    "TFII-I", "GTF2I", "SIX1", "SIX2", "SIX3", "LEF-1", "LEF1",
    "ETS-1", "ETS1", "VP16", "ERA", "ESR1", "ZTA", "BZLF1",
    "AHR", "SOX11", "SOX4", "SOX5", "SOX6", "SOX17",
    "E2A-PBX1", "PAX3-FOXO1", "EWS-FLI1", "TMPRSS2-ERG",
    "HOXD13", "PITX3", "PROX1", "NR2E1", "TLX", "NR2F1",
    "COUP-TFI", "NR2F2", "COUP-TFII", "NR2F6", "EAR2",
    "HNF4G", "HNF4A", "HNF4B", "PPARGC1A", "PGC1A",
    "TFCP2", "LBP", "CP2", "YBX1", "YB-1", "NF-YA",
    "NF-YB", "NF-YC", "DUX4", "DUX4-FL", "ZNF750", "GRHL2",
    "OVOL1", "OVOL2", "TEAD1", "TEAD2", "TEAD3", "TEAD4",
    "YAP1", "TAZ", "WWTR1", "MAML1", "MAML2", "NCID",
    "RBPJK", "SUH", "DLL1", "CBF1/RBPJK", "MEF2B", "TCF7L1",
    "TCF7L2", "TCF3", "MESP1", "MESP2", "TBX1", "TBX2",
    "TBX3", "TBX4", "TBX5", "TBX6", "TBX18", "TBX21", "T-BET",
    "EOMES", "TBR1", "TBR2", "EOMES", "NEUROD1", "NEUROD2",
    "NEUROG1", "NEUROG2", "NGN2", "NGN3", "ASCL1", "MASH1",
    "ASCL2", "MATH1", "ATOH1", "PTF1A", "OLIG1", "OLIG2",
    "NEUROG3", "BHLHE40", "DEC1", "BHLHE41", "DEC2",
    "CLOCK", "BMAL1", "ARNTL", "PER1", "PER2", "CRY1", "CRY2",
    "RORA", "NR1D1", "REV-ERBA", "DBP", "HLF", "TEF",
    # 高频真 TF 补充（第三轮探测）
    "PDX1", "USF-1", "RFX", "RFX1", "RFX2", "RFX3", "RFX4",
    "GRHL3", "GRHL2", "GRHL1", "HIF2A", "EPAS1", "TEL", "ETV6",
    "RELA", "P65", "GSX", "GSX1", "GSX2", "E2F-3", "C/EBPA",
    "NF KAPPA B", "RELA/P65", "P65/RELA", "PRB",
    "TCF3::PBX1", "TEL/AML1", "AML1-ETO", "ETV6-RUNX1",
    "RUNX1-ETO", "C-EBP", "C/EBP DELTA", "C/EBP GAMMA",
    "C/EBP EPSILON", "C/EBP ZETA", "GADD153",
    "HNF4ALPHA", "HNF4BETA", "C/EBPALPHA", "C/EBPBETA",
    "C/EBPGAMMA", "C/EBPDELTA", "NF-YA", "NF-YB", "NF-YC",
    "FOXA1A", "FOXA1B", "GATA-1A", "GATA1A",
    # 全称/酵母 TF（文献中出现的真 TF）
    "ACTIVATOR PROTEIN-1", "ACTIVATOR PROTEIN 1", "ABF1", "ACE1",
    "ACE2", "GAL4", "GCN4", "ADR1", "ADR-1", "MIG1", "PHO4",
    "STE12", "SWI5", "MBF", "SCB", "MCB", "CBF1", "CPF1",
    "HAP1", "HAP2", "HAP3", "HAP4", "HAP5", "BAS1", "BAS2",
    "PHO2", "RAP1", "GRF1", "GCR1", "GCR2", "REB1", "TUF",
}
TF_WHITELIST = _JASPAR_MOTIF_TFS | _COMMON_TFS

# 明确不是 TF 的信号蛋白/方法名（即使出现在 tf 字段也剔除）
_NON_TF_PAT = re.compile(
    r"(KINASE|RECEPTOR|LIGASE|PHOSPHATASE|DEMETHYLASE|METHYLTRANSFERASE|"
    r"HISTONE ACETYLTRANSFERASE|SELEX|CUT&RUN|CUT & RUN|CHIP-SEQ|CHIP SEQ|"
    r"MPRA|CRISPR|RNA-SEQ|ATAC-SEQ|MOTIF SCAN|SCREEN|ASSAY|SUMO\b|"
    r"\bTF$\b|\bTFS\b)", re.I)
# 精确匹配的泛化/非 TF 实体（避免用宽泛词误伤 TF 全称）
_NON_TF_EXACT = {
    "TF", "TFS", "E-BOX", "E BOX", "ETS FAMILY", "GENE", "PROTEIN",
    "BINDING SITES", "BINDING SITE", "DNA SEQUENCES", "MOTIFS",
    "TRANSCRIPTION FACTOR", "TRANSCRIPTION FACTORS", "REGULATORS",
    "ACTIVATORS", "SILENCERS", "HOMEODOMAIN", "MASTER REGULATORS",
    "MASTER REGULATORS OF IMMUNE DIFFERENTIATION",
}

# 希腊字母 -> 拉丁字母（模型输出常带 Unicode 希腊字母，如 NF-ΚB/HIF1Α）
_GREEK_MAP = str.maketrans({
    "Α": "A", "Β": "B", "Γ": "G", "Δ": "D", "Ε": "E", "Ζ": "Z",
    "Η": "E", "Θ": "TH", "Ι": "I", "Κ": "K", "Λ": "L", "Μ": "M",
    "Ν": "N", "Ξ": "X", "Ο": "O", "Π": "P", "Ρ": "R", "Σ": "S",
    "Τ": "T", "Υ": "Y", "Φ": "PH", "Χ": "CH", "Ψ": "PS", "Ω": "O",
    "α": "a", "β": "b", "γ": "g", "δ": "d", "ε": "e", "ζ": "z",
    "η": "e", "θ": "th", "ι": "i", "κ": "k", "λ": "l", "μ": "m",
    "ν": "n", "ξ": "x", "ο": "o", "π": "p", "ρ": "r", "σ": "s",
    "τ": "t", "υ": "y", "φ": "ph", "χ": "ch", "ψ": "ps", "ω": "o",
})


# ---------------------------------------------------------------------------
# 基序特征: binds_motif 的 tail 必须是真基序
# ---------------------------------------------------------------------------
_MOTIF_KW = re.compile(
    r"\b(MOTIF|SITE|BOX|ELEMENT|CONSENSUS|TATA|CCAAT|GC-RICH|AT-RICH|"
    r"RESPONSE ELEMENT|ENHANCER|PROMOTER|BINDING SEQUENCE|RECOGNITION "
    r"SEQUENCE|RE)\b", re.I)
# IUPAC 核酸序列样式（如 CACACA、ACAGATGG、(A/T)GATA(A/G)、5'-GTGCGTG-3'）
_IUPAC_SEQ = re.compile(r"^[ACGTURYSWKMBDHVN()/-]{4,30}$", re.I)
_SEQ_5P = re.compile(r"5'?-?[ACGTN]{4,}")

# 明确不是基序的 tail（TF 名/基因名误放 motif 字段）
_NON_MOTIF_KW = re.compile(
    r"(KINASE|RECEPTOR|GENE|PROTEIN\b|COMPLEX\b|PATHWAY|SIGNALING|"
    r"EXPRESSION|REGULATION|CELL LINE|TRANSFECTION)", re.I)

# 常见基序/响应元件缩写白名单（无特征词但确为真基序）
_MOTIF_NAMES = {
    "AP-1", "AP1", "AP-2", "AP2", "OCTAMER", "CRE", "CRE2", "ARE",
    "ERE", "GRE", "PRE", "SRE", "TRE", "AHRE", "IIARE", "DRE", "XRE",
    "MRE", "BRE", "INR", "DPE", "ISRE", "GAS", "TATA", "CCAAT",
    "GC BOX", "GC BOXES", "G-BOX", "E-BOX", "E BOX", "E2F SITE",
    "SP1 SITE", "P63RE", "ETS", "ETS1", "YY1", "SP1", "CACACA",
    "G-QUADRUPLEX", "G4", "GA-RICH", "TG-RICH", "CG-RICH", "AT-RICH",
    "PALINDROME", "TFBS", "PWM",
}
# 以 RE 结尾的响应元件缩写（AHRE/IIARE/P63RE...）
_RE_ABBR = re.compile(r"^[A-Z0-9]{2,6}RE$")


def looks_like_motif(tail: str) -> bool:
    """tail 是否是"真基序"（含基序特征词 / 缩写白名单 / IUPAC 序列）。"""
    t = tail.strip()
    if _NON_MOTIF_KW.search(t) and not _MOTIF_KW.search(t):
        return False
    if _MOTIF_KW.search(t):
        return True
    if t in _MOTIF_NAMES or _RE_ABBR.match(t):
        return True
    if _IUPAC_SEQ.match(t) or _SEQ_5P.search(t):
        return True
    return False


# ---------------------------------------------------------------------------
# 实体归一化（与 review_extractions.normalize 一致）
# ---------------------------------------------------------------------------
_TAG_RE = re.compile(r"<[^>]+>")


def normalize(s) -> str:
    s = _TAG_RE.sub("", str(s))
    s = s.upper()
    s = s.translate(_GREEK_MAP)
    s = re.sub(r"\s+", " ", s).strip()
    return s.strip(" .,;:!?\"'()[]{}<>-—–")


def load_flagged(flagged_path: Path) -> dict[str, set[str]]:
    """加载 flagged 幻觉清单（reason=not_in_abstract），返回 pmid -> 排除实体。"""
    result: dict[str, set[str]] = {}
    if not flagged_path.exists():
        return result
    with open(flagged_path, "r", encoding="utf-8") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            try:
                item = json.loads(line)
            except json.JSONDecodeError:
                continue
            if item.get("reason") != "not_in_abstract":
                continue
            pmid = str(item.get("pmid"))
            ent = str(item.get("entity", "")).strip()
            if pmid and ent:
                result.setdefault(pmid, set()).add(normalize(ent))
    return result


def build_triples(record: dict) -> list[tuple[str, str, str]]:
    """从抽取记录构造三元组（与 build_knowledge_base.build_triples 一致）。"""
    triples: set[tuple[str, str, str]] = set()
    tfs = {normalize(e) for e in record.get("tf", []) if e}
    genes = {normalize(e) for e in record.get("gene", []) if e}
    motifs = {normalize(e) for e in record.get("motif", []) if e}
    diseases = {normalize(e) for e in record.get("disease", []) if e}

    for tf in tfs:
        for g in genes:
            triples.add((tf, "regulates", g))
        for m in motifs:
            triples.add((tf, "binds_motif", m))
        for d in diseases:
            triples.add((tf, "associated_with", d))
    for g in genes:
        for d in diseases:
            triples.add((g, "associated_with", d))
    return list(triples)


def load_literature_records(extraction_files: list[Path],
                            flagged_files: list[Path]) -> list[dict]:
    """读取所有抽取记录，应用幻觉/黑名单过滤，返回清洗后的记录。"""
    records: list[dict] = []
    for path in extraction_files:
        if not path.exists():
            logger.warning("抽取文件不存在: %s", path)
            continue
        with open(path, "r", encoding="utf-8") as f:
            for line in f:
                line = line.strip()
                if not line:
                    continue
                try:
                    rec = json.loads(line)
                except json.JSONDecodeError:
                    continue
                if not rec.get("parsed"):
                    continue
                records.append(rec)

    flagged: dict[str, set[str]] = {}
    for fp in flagged_files:
        f = load_flagged(fp)
        for pmid, ents in f.items():
            flagged.setdefault(pmid, set()).update(ents)

    cleaned: list[dict] = []
    n_hal = n_bl = 0
    for rec in records:
        pmid = str(rec.get("pmid", ""))
        exclude = flagged.get(pmid, set())
        for field in ("tf", "gene", "motif", "disease"):
            keep = []
            for e in rec.get(field, []):
                ne = normalize(e)
                if not ne:
                    continue
                if ne in exclude:
                    n_hal += 1
                    continue
                if is_blacklisted(ne):
                    n_bl += 1
                    continue
                keep.append(e)
            rec[field] = keep
        cleaned.append(rec)
    logger.info("  幻觉排除实体: %d, 黑名单排除实体: %d", n_hal, n_bl)
    return cleaned


def literature_triplets(records: list[dict]) -> list[tuple[str, str, str]]:
    """从清洗记录生成全部文献三元组并去重。"""
    all_triples: set[tuple[str, str, str]] = set()
    for rec in records:
        all_triples.update(build_triples(rec))
    return sorted(all_triples)


def make_actionable(head: str, rel: str, tail: str, ktype: str) -> dict:
    """把过滤后的三元组包装为行动型记录。"""
    if rel == "binds_motif":
        recommendation = (
            f"{head} binds motif {tail} (literature evidence). Including this "
            f"motif in a candidate sequence creates a binding opportunity for "
            f"{head}; use it as a sequence constraint when designing for {head} "
            f"activity."
        )
    elif rel == "regulates":
        recommendation = (
            f"{head} regulates {tail} (literature evidence). To drive {tail} "
            f"expression, include binding sites for the TF {head} in the "
            f"candidate sequence."
        )
    else:  # 兜底
        recommendation = f"{head} {rel} {tail} (literature evidence)."
    return {
        "head": head,
        "relation": rel,
        "tail": tail,
        "effect_size": None,
        "confidence": {"evidence_type": "literature", "source": "PubMed"},
        "recommendation": recommendation,
        "source": SOURCE_LIT,
        "knowledge_type": ktype,
    }


def main() -> None:
    parser = argparse.ArgumentParser(description="分层过滤知识库（方案A）")
    parser.add_argument("--outdir", default="data/processed")
    parser.add_argument("--output", default="data/processed/knowledge_base_actionable.jsonl")
    parser.add_argument("--old-actionable", default="data/processed/knowledge_base_actionable.jsonl",
                        help="现有 actionable 库（提取其中的 MPRA 统计类型1/2）")
    parser.add_argument("--extractions", nargs="+", default=[
        "data/processed/raw_extractions.jsonl",
        "data/processed/raw_extractions_mpra.jsonl",
    ])
    parser.add_argument("--flagged", nargs="+", default=[
        "data/processed/flagged_extractions_v4.jsonl",
        "data/processed/flagged_extractions_mpra_v4.jsonl",
    ])
    args = parser.parse_args()

    outdir = ROOT / args.outdir
    output_path = ROOT / args.output

    # ---------- 1) MPRA 统计层（类型1/2）：从现有 actionable 库提取 ----------
    stat_triplets: list[dict] = []
    old_path = ROOT / args.old_actionable
    if old_path.exists():
        with open(old_path, "r", encoding="utf-8") as f:
            for line in f:
                line = line.strip()
                if not line:
                    continue
                try:
                    t = json.loads(line)
                except json.JSONDecodeError:
                    continue
                if t.get("knowledge_type") in ("motif_activity", "gc_activity"):
                    stat_triplets.append(t)
    logger.info("MPRA 统计层（类型1/2）: %d 条", len(stat_triplets))

    # ---------- 2) 文献层重建 ----------
    extraction_files = [ROOT / p for p in args.extractions]
    flagged_files = [ROOT / p for p in args.flagged]
    records = load_literature_records(extraction_files, flagged_files)
    logger.info("清洗后抽取记录: %d", len(records))

    triples = literature_triplets(records)
    logger.info("文献三元组（去重后）: %d", len(triples))

    seq_actionable: list[dict] = []
    gene_reg: list[dict] = []
    dropped = Counter()
    dropped_examples: dict[str, list[str]] = {"self_ref": [], "disease": [],
                                              "bad_binds": [], "bad_regulates": []}

    for head, rel, tail in triples:
        # 自指删除
        if head == tail:
            dropped["self_ref"] += 1
            if len(dropped_examples["self_ref"]) < 5:
                dropped_examples["self_ref"].append(f"{head}->{tail}")
            continue
        if rel == "associated_with":
            dropped["disease_assoc"] += 1
            if len(dropped_examples["disease"]) < 5:
                dropped_examples["disease"].append(f"{head}->{tail}")
            continue
        if rel == "binds_motif":
            if head in TF_WHITELIST and head not in _NON_TF_EXACT \
                    and not _NON_TF_PAT.search(head) \
                    and looks_like_motif(tail):
                seq_actionable.append(make_actionable(head, rel, tail,
                                                      "sequence_actionable"))
            else:
                dropped["bad_binds_motif"] += 1
                if len(dropped_examples["bad_binds"]) < 5:
                    dropped_examples["bad_binds"].append(f"{head}->{tail}")
            continue
        if rel == "regulates":
            # head 必须是真 TF（白名单），且不是信号蛋白/方法名
            if head in TF_WHITELIST and head not in _NON_TF_EXACT \
                    and not _NON_TF_PAT.search(head):
                gene_reg.append(make_actionable(head, rel, tail, "gene_regulation"))
            else:
                dropped["bad_regulates_head"] += 1
                if len(dropped_examples["bad_regulates"]) < 5:
                    dropped_examples["bad_regulates"].append(f"{head}->{tail}")
            continue
        # 其他关系（理论上没有）丢弃
        dropped[f"other:{rel}"] += 1

    # ---------- 3) 写输出 ----------
    all_out = stat_triplets + seq_actionable + gene_reg
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with open(output_path, "w", encoding="utf-8") as f:
        for t in all_out:
            f.write(json.dumps(t, ensure_ascii=False) + "\n")

    logger.info("=" * 60)
    logger.info("方案A 过滤统计:")
    logger.info("  序列行动层  motif_activity     : %d", sum(1 for t in stat_triplets if t["knowledge_type"]=="motif_activity"))
    logger.info("  序列行动层  gc_activity        : %d", sum(1 for t in stat_triplets if t["knowledge_type"]=="gc_activity"))
    logger.info("  序列行动层  binds_motif(真基序): %d", len(seq_actionable))
    logger.info("  基因调控层  regulates(真TF)    : %d", len(gene_reg))
    logger.info("  保留总数                       : %d", len(all_out))
    logger.info("  丢弃明细                       : %s", dict(dropped))
    for k, ex in dropped_examples.items():
        if ex:
            logger.info("  丢弃样例 %s: %s", k, " | ".join(ex))
    logger.info("  输出文件                       : %s", output_path)
    logger.info("=" * 60)


if __name__ == "__main__":
    main()
