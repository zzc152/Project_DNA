#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
gene_names.py — 内置常用基因名表（用户决策 2026-08-06）

用途：判定 regulatory_element 前缀（如 "FKBP52 promoter" 的 FKBP52）是否为基因名。
判定方式：前缀命中原样/归一化后的基因名表 → 是基因特异；否则非基因。

设计说明：
- 不追求全 HGNC，只内置本项目语境下常见基因 + 知识库高频 factor
- 判定用归一化形式：去连字符/空格、全大写（GATA-1 → GATA1, gata2 → GATA2）
- 这是"常用基因名表"，内部知识库用它自洽（factor 交叉验证）
"""

# 归一化：去连字符/点/空格 → 全大写
def _norm(name: str) -> str:
    return name.replace("-", "").replace(".", "").replace(" ", "").upper()


# β 等希腊字母归一化为拉丁字母（β-globin → BGLOBIN，与 HBB 别名一致）
_GREEK_MAP = {"β": "B", "α": "A", "γ": "G", "δ": "D", "ε": "E", "κ": "K"}


def _norm(name: str) -> str:
    s = name.replace("-", "").replace(".", "").replace(" ", "")
    out = []
    for ch in s:
        if ch in _GREEK_MAP:
            out.append(_GREEK_MAP[ch])
        else:
            out.append(ch)
    return "".join(out).upper()


# 内置常用基因名（原始写法，运行时归一化）。覆盖本项目出现的高频基因。
_RAW_GENE_NAMES = [
    # --- 转录因子 / 调控因子（知识库高频） ---
    "MYCN", "CREB", "ETS1", "WDR5", "SMARCB1", "BRG1", "FOXP4", "KLF4",
    "BCLAF1", "E2F1", "NFIB", "YY1", "VAV1", "ETO2", "GATA2", "GATA-1",
    "GATA-2", "GATA3", "GATA6", "BCL11A", "MECOM", "ZNF410", "CTCF",
    "TFAP2A", "TFAP2C", "LEF1", "BRD4", "SMARCA5", "PBX1", "PBX3", "KAT7",
    "SOX2", "SOX4", "SOX11", "CBP", "TAL1", "TR4", "THRA", "RUNX1",
    "RUNX2", "HDAC1", "HDAC2", "OGT", "SMAD3", "NUP153", "ZNF692",
    "ALKBH5", "FOXO1", "FUBP3", "ASCL1", "HAND2", "ATF4", "ATF5",
    "BACH1", "KDM4C", "DNMT1", "DNMT3A", "DNMT3B", "SREBP1", "PHOX2B",
    "MBD2", "NEAT1", "DDX27", "HOXB13", "HOXA10", "KMT2B", "RYBP",
    "EP300", "SMYD2", "CRX", "OTX2", "GLI3", "CEBPB", "CEBPA", "MYB",
    "FOG1", "LYAR", "PRMT5", "SP1", "SP3", "FOXA1", "CaMKII", "NeuroD1",
    "SNRPA", "RepID", "OGA", "MEF2A", "NCOA4", "ODC1", "Osa", "KLF6",
    # --- 常见通用基因 ---
    "TP53", "MYC", "RAS", "AKT", "ERK", "STAT3", "STAT5", "JAK2",
    "EGFR", "VEGFA", "HIF1A", "HIF-1A", "AR", "ER", "ESR1", "PTEN",
    "RB1", "CDK2", "CDK4", "CCND1", "BCL2", "BAX", "CASP3", "MAPK1",
    "MAPK3", "PIK3CA", "KRAS", "NRAS", "BRAF", "RAF1", "MEK1", "MEK2",
    "IL6", "IL-6", "TNF", "TGFB1", "TGF-β1", "TGFB1", "CSF1R", "FKBP52",
    "MAFA", "BCL2", "FOXP4", "NFKB1", "RELA", "P65", "MTOR", "AMPK",
    "SIRT1", "HDAC3", "HDAC6", "DNMT1", "TET2", "TET1", "EZH2", "SUZ12",
    "BMI1", "RING1B", "PRC2", "KDM5A", "KDM6A", "KDM6B", "UTX", "JMJD3",
    "BRCA1", "BRCA2", "ATM", "ATR", "CHEK1", "CHEK2", "WEE1", "PLK1",
    "AURKA", "AURKB", "MCM2", "MCM7", "PCNA", "TOP1", "TOP2A",
    # --- 珠蛋白/红细胞系 ---
    "HBB", "HBA1", "HBA2", "GATA1", "KLF1", "NFE2", "TAL1", "LMO2",
    "ERG", "FLI1", "PU.1", "SPI1", "IRF4", "IRF8", "CEBPA",
    "β-globin", "β globin",  # HBB 家族别称（归一化 → BGLOBIN）
    # --- 本项目已出现但易遗漏的常用基因 ---
    "CYP26B1", "Ghd7", "SVCT2", "hSVCT2", "STING", "TMEM173",
    "EHD1", "BDNF", "GHD7",
]

# 归一化后的基因名集合（用于判定）
GENE_NAMES = {_norm(g) for g in _RAW_GENE_NAMES}

# 补充：从知识库 factor 自动并入（模块加载时由调用方传入，见 register_kb_genes）
_KB_GENES = set()


def register_kb_genes(factors) -> None:
    """将知识库 factor 中的'纯基因符号'并入基因名表（内部自洽）。"""
    import re
    pat = re.compile(r"^[A-Za-z][A-Za-z0-9_.]{1,14}$")
    for f in factors:
        s = str(f).strip()
        if pat.match(s) and s.lower() not in (
            "enhancer", "promoter", "hypoxia", "erna", "eqtls",
            "venlafaxine", "dnam",
        ):
            _KB_GENES.add(_norm(s))
    # 刷新 GENE_NAMES（引用同一 set，避免重建）
    GENE_NAMES.update(_KB_GENES)


def is_gene_name(token: str) -> bool:
    """判定 token 是否为基因名（归一化后查表）。"""
    if not token:
        return False
    return _norm(token) in GENE_NAMES


if __name__ == "__main__":
    # 自检：本项目带前缀元件判定
    tests = ["FKBP52", "core", "cis-regulatory", "foxp4", "bcl2", "globin",
             "CSF1R", "ATF5", "MAFA", "il-6", "β-globin", "hypoxia",
             "enhancer", "neuronal"]
    for t in tests:
        print(f"  {t!r:20} -> {'基因' if is_gene_name(t) else '非基因'}")
