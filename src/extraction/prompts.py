"""提示词模板：用于指导 Qwen2.5-7B-Instruct 从 PubMed 摘要中抽取 DNA 调控知识。

输出 JSON 格式约定：
{
    "tf": [转录因子列表],
    "gene": [靶基因列表],
    "motif": [DNA 结合基序列表],
    "disease": [相关疾病列表],
    "relation": "调控关系描述（最多两句话）"
}

注意：
- 不要使用 "return empty dict" 之类的表述——实测会让模型倾向于输出空 JSON
  而不去提取摘要中实际存在的实体。改为明确要求列出所有实体，未命中的类型
  返回空列表 []。
- 必须给出精确的实体定义 + 排除规则，否则模型会把 "transcription factors"、
  "transfer RNAs"、"DNA sequences" 这类描述性词汇误当成实体名。
- Few-shot 示例（FEW_SHOT_EXAMPLES）对模型约束力最强：正例展示正确抽取，
  反例直接展示此前实测出现的幻觉（25363779 中凭空输出 TP53/MYC/E-box 等），
  并给出该摘要的正确空列表输出。
"""

EXTRACTION_SYSTEM_PROMPT = (
    "You are a biomedical information extraction assistant specializing in "
    "human gene regulation and transcription factors. "
    "Extract ONLY entities that literally appear as specific names in the "
    "abstract. "
    "STRICT RULES: (1) Never invent or hallucinate entities that are not "
    "explicitly written in the text. (2) If the abstract mentions no specific "
    "TF, gene, motif, or disease names, return an empty list for that field. "
    "(3) Do not generalize: generic descriptions such as 'transcription "
    "factors', 'regulators', 'binding sites', or 'disease' are NOT entities. "
    "(4) Only include a name if you saw that exact name in the abstract."
)

# 实体定义与排除规则（注入到 user 提示词中，帮助模型区分真正的命名实体
# 与描述性/类别性词汇）
ENTITY_DEFINITIONS = (
    'Definitions of the entity types:\n'
    '- "tf": specific transcription factor proteins by their gene symbol or '
    'protein name, e.g. TP53, MYC, TWIST1, NGN2, FOXA1. '
    'Do NOT include generic or category words such as "transcription '
    'factors", "TFs", "master regulators", "activators", "silencers", '
    '"transfer RNAs", protein family names such as "bHLH" or "homeodomain '
    '(HD)", non-TF proteins/enzymes (e.g. Tn5), or any name that does not '
    'appear verbatim in the abstract.\n'
    '- "gene": specific target genes that appear verbatim in the abstract, '
    'especially genes described as regulated by transcription factors, '
    'e.g. "regulates MYOD1 expression".\n'
    '- "motif": specific DNA binding motifs or sequence patterns, e.g. E-box, '
    'CACACA, GC box, CpG island. Do NOT include generic words such as "DNA '
    'sequences", "binding sites", "motifs", protein-family names, or custom '
    'names invented by the paper (e.g. "Coordinator", "composite motifs").\n'
    '- "disease": specific disease names, e.g. breast cancer, type 2 '
    'diabetes. Do NOT include phrases like "disease risk alleles", '
    '"disease associations", or generic "autoimmune diseases" unless a '
    'specific disease name is given.\n'
    '\n'
    'HARD BLACKLIST: never output any of these terms in any field, even if '
    'they appear in the abstract: "master regulators", "Tn5", "transcription '
    'factors", "TFs", "regulators", "activators", "silencers", "bHLH", '
    '"homeodomain", "HD", "transfer RNAs", "DNA sequences", "binding sites", '
    '"Coordinator", "composite motifs", "motifs", "disease", "autoimmune '
    'diseases", "disease risk alleles". These are generic descriptions, '
    'protein families, or tools, not named entities.\n'
)

# Few-shot 示例：模型对示例的服从度远高于文字描述。
# 正例（EXAMPLE 1/2）取自真实摘要的正确抽取结果；
# 反例（EXAMPLE 3）展示模型此前犯过的幻觉错误（把原文没有的
# TP53/MYC/E-box 等编造出来），并给出该摘要的正确空列表输出。
FEW_SHOT_EXAMPLES = """Follow these examples exactly:

=== EXAMPLE 1 (correct) ===
Abstract: "We experimentally profile DNA methylation and combine this with published occupancy profiles of five distinct TFs (CTCF, CEBPA, HNF4A, ONECUT1, FOXA1) in the liver of five mammalian species."
JSON: {"tf": ["CTCF", "CEBPA", "HNF4A", "ONECUT1", "FOXA1"], "gene": [], "motif": [], "disease": [], "relation": "TF binding occupancy of CTCF, CEBPA, HNF4A, ONECUT1 and FOXA1 is profiled in liver."}

=== EXAMPLE 2 (correct) ===
Abstract: "Many TFs, including canonical activators such as NRF1, NFY and Sp1, activate or repress transcription initiation depending on their precise position relative to the TSS."
JSON: {"tf": ["NRF1", "NFY", "Sp1"], "gene": [], "motif": [], "disease": [], "relation": "NRF1, NFY and Sp1 activate or repress transcription initiation depending on their position relative to the TSS."}

=== EXAMPLE 3 (WRONG - hallucination, NEVER do this) ===
Abstract: "Causal variants tend to occur near binding sites for master regulators of immune differentiation."  (no specific TF, gene, motif or disease name appears in the text)
WRONG JSON (do NOT output): {"tf": ["TP53", "MYC", "TWIST1", "NGN2", "FOXA1"], "motif": ["E-box", "CACACA", "GC box"], "disease": ["autoimmune diseases"]}
CORRECT JSON: {"tf": [], "gene": [], "motif": [], "disease": [], "relation": "No specific entities are named; only generic terms such as binding sites and master regulators appear."}
"""

EXTRACTION_USER_TEMPLATE = (
    "Read the abstract below and extract every specific transcription factor "
    "(TF), target gene, DNA binding motif, and disease entity mentioned.\n\n"
    f"{ENTITY_DEFINITIONS}\n"
    f"{FEW_SHOT_EXAMPLES}\n"
    'Respond ONLY in this exact JSON format: '
    '{"tf": [...], "gene": [...], "motif": [...], "disease": [...], "relation": "..."}. '
    "Use an empty list [] for entity types with no matches. "
    "relation should be a short description of the regulatory relationship "
    "(max 2 sentences).\n\n"
    "Abstract: {abstract}\n\nJSON:"
)


def build_prompt(abstract_text: str) -> str:
    """构造提示词（字符串拼接形式，供非 chat 模板使用）。"""
    return (
        f"{EXTRACTION_SYSTEM_PROMPT}\n"
        f"{EXTRACTION_USER_TEMPLATE.replace('{abstract}', abstract_text)}"
    )


def build_messages(abstract_text: str) -> list[dict]:
    """构造 Chat 消息格式（推荐用于 Qwen2.5-Instruct 的 apply_chat_template）。"""
    return [
        {"role": "system", "content": EXTRACTION_SYSTEM_PROMPT},
        {
            "role": "user",
            "content": EXTRACTION_USER_TEMPLATE.replace("{abstract}", abstract_text),
        },
    ]
