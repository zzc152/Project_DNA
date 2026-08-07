"""提示词模板 · 文献知识源（Source 2）：三系调控元件活性影响原因抽取。

与 prompts.py（三元组抽取）不同，本模板聚焦"因果知识"：
    什么因素/机制 影响 哪些调控元件 在哪个细胞系中的活性。

输出 JSON schema（每条摘要可包含多条 finding）:
{
    "findings": [
        {
            "cell_line": "K562" | "HepG2" | "SK-N-SH" | "not_specified",
            "regulatory_element": "enhancer" | "promoter" | "cis-regulatory element" | ...,
            "factor": "影响活性的具体因素（TF / 基序 / 序列特征 / 表观修饰）",
            "effect": "increases" | "decreases" | "required_for" | "modulates" | "no_effect",
            "mechanism": "机制描述（1-2 句，原文支持的因果解释）",
            "evidence_sentence": "摘要中支持该结论的最直接证据句（原文）"
        }
    ]
}

设计要点（吸取上一版抽取教训）:
- 严禁幻觉：cell_line / factor / regulatory_element 必须出自原文；摘要未
  指明细胞系时 cell_line 用 "not_specified"，绝不猜测。
- effect 用受控词汇表，避免自由文本带来的下游解析困难。
- mechanism 必须是原文可支持的因果描述，不是模型编造的通识。
- evidence_sentence 必须逐字来自摘要（允许小幅截断）。
- Few-shot 正反例对约束力最强，反例直接展示此前常见的幻觉模式。
"""

ACTIVITY_SYSTEM_PROMPT = (
    "You are a biomedical knowledge extraction expert specializing in "
    "gene regulation and cis-regulatory elements (enhancers, promoters). "
    "Your task: extract KNOWLEDGE ABOUT WHAT DETERMINES OR AFFECTS THE "
    "ACTIVITY of regulatory elements (enhancers/promoters) — i.e. causal "
    "or mechanistic findings — from each abstract. "
    "STRICT RULES: (1) Never invent facts. Every extracted factor, cell "
    "line and regulatory element must literally appear in the abstract. "
    "(2) Only extract statements about factors that AFFECT regulatory "
    "element activity/expression. Ignore background, methods, and "
    "findings about other topics. (3) If the abstract does not mention a "
    "specific cell line, set cell_line to \"not_specified\" — never guess "
    "one. (4) The mechanism and effect must be supported by the abstract, "
    "not by general biology knowledge. (5) Use only the controlled "
    "vocabulary for effect. (6) Quote the evidence sentence verbatim from "
    "the abstract (you may trim it, but do not paraphrase into new words)."
)

ACTIVITY_FIELD_DEFINITIONS = (
    'Definitions of the JSON fields:\n'
    '- "cell_line": one of "K562", "HepG2", "SK-N-SH", or "not_specified". '
    'Use the exact cell line name if the abstract states experiments in '
    'K562, HepG2, or SK-N-SH cells. If the abstract studies other cells '
    '(e.g. HEK293T, mouse cells) or no cells at all, use "not_specified". '
    'NEVER infer a cell line that is not named.\n'
    '- "regulatory_element": the type of regulatory element whose activity '
    'is discussed, e.g. "enhancer", "promoter", "cis-regulatory element", '
    '"enhancer-promoter interaction", "silencer", "insulator", "locus '
    'control region". Use a concise lowercase term that matches the paper. '
    'Do not invent elements not mentioned.\n'
    '- "factor": the specific entity that affects the element activity — '
    'a transcription factor (e.g. GATA1, KLF1, HNF4A, NEUROG2), a DNA '
    'motif (e.g. E-box, GATA motif), a sequence feature (e.g. GC content, '
    'DNA shape, motif position), or an epigenetic mark (e.g. H3K27ac, '
    'DNA methylation, chromatin accessibility). Use the exact name from '
    'the abstract. Do NOT use generic terms like "transcription factors", '
    '"chromatin", "the motif" unless that is all the abstract says.\n'
    '- "effect": the direction of influence, from this controlled set:\n'
    '    "increases"     — factor boosts element activity/expression\n'
    '    "decreases"     — factor reduces element activity/expression\n'
    '    "required_for"  — factor is necessary for element activity\n'
    '    "modulates"     — factor changes activity, direction context-dependent\n'
    '    "no_effect"     — factor has no effect on activity\n'
    '- "mechanism": a 1-2 sentence mechanistic explanation of HOW the '
    'factor affects activity, strictly based on what the abstract states '
    '(e.g. "GATA1 binds the enhancer and recruits co-activators to '
    'activate transcription"). If the abstract gives no mechanism, '
    'describe only the observed relationship.\n'
    '- "evidence_sentence": the single most direct sentence from the '
    'abstract that supports this finding, quoted verbatim.\n'
)

# Few-shot：正例展示正确抽取（含"not_specified"与受控词汇用法），
# 反例展示幻觉/越权推断（编造细胞系、编造机制、把背景当因果）。
ACTIVITY_FEW_SHOT_EXAMPLES = """Follow these examples exactly:

=== EXAMPLE 1 (correct) ===
Abstract: "We show that the erythroid transcription factor GATA1 is required for the activity of the beta-globin locus control region. Deletion of the GATA1 binding sites in the HS2 enhancer abolished enhancer-driven reporter expression in K562 cells, indicating that GATA1 occupancy is necessary for enhancer function."
JSON: {"findings": [{"cell_line": "K562", "regulatory_element": "enhancer", "factor": "GATA1", "effect": "required_for", "mechanism": "GATA1 binding sites in the HS2 enhancer are necessary for enhancer-driven reporter expression; GATA1 occupancy is required for enhancer function.", "evidence_sentence": "Deletion of the GATA1 binding sites in the HS2 enhancer abolished enhancer-driven reporter expression in K562 cells."}]}

=== EXAMPLE 2 (correct, no cell line) ===
Abstract: "Nucleosome positioning at promoters is a major determinant of transcriptional output. Here we show that high nucleosome occupancy near the transcription start site represses promoter activity independent of the cell type."
JSON: {"findings": [{"cell_line": "not_specified", "regulatory_element": "promoter", "factor": "nucleosome occupancy", "effect": "decreases", "mechanism": "High nucleosome occupancy near the transcription start site represses promoter activity.", "evidence_sentence": "high nucleosome occupancy near the transcription start site represses promoter activity independent of the cell type."}]}

=== EXAMPLE 3 (WRONG - hallucination, NEVER do this) ===
Abstract: "Enhancers are important for gene regulation. We studied chromatin features associated with active enhancers."  (no cell line, no specific factor, no causal claim)
WRONG JSON (do NOT output): {"findings": [{"cell_line": "K562", "regulatory_element": "enhancer", "factor": "H3K27ac", "effect": "increases", "mechanism": "H3K27ac is a well-known mark of active enhancers that promotes transcription.", "evidence_sentence": "We studied chromatin features associated with active enhancers."}]}
CORRECT JSON: {"findings": []}

=== EXAMPLE 4 (WRONG - cell line invented, NEVER do this) ===
Abstract: "HNF4A binding is required for hepatocyte-specific enhancer activity."  (no cell line is named in the abstract)
WRONG JSON (do NOT output): {"findings": [{"cell_line": "HepG2", "regulatory_element": "enhancer", "factor": "HNF4A", "effect": "required_for", "mechanism": "HNF4A is a liver-enriched factor that activates hepatocyte enhancers.", "evidence_sentence": "HNF4A binding is required for hepatocyte-specific enhancer activity."}]}
CORRECT JSON: {"findings": [{"cell_line": "not_specified", "regulatory_element": "enhancer", "factor": "HNF4A", "effect": "required_for", "mechanism": "HNF4A binding is required for hepatocyte-specific enhancer activity.", "evidence_sentence": "HNF4A binding is required for hepatocyte-specific enhancer activity."}]}
"""

ACTIVITY_USER_TEMPLATE = (
    "Read the abstract below and extract every finding about what "
    "determines or affects the activity of regulatory elements "
    "(enhancers/promoters).\n\n"
    f"{ACTIVITY_FIELD_DEFINITIONS}\n"
    f"{ACTIVITY_FEW_SHOT_EXAMPLES}\n"
    'Respond ONLY with a JSON object in this exact format: '
    '{"findings": [{"cell_line": "...", "regulatory_element": "...", '
    '"factor": "...", "effect": "...", "mechanism": "...", '
    '"evidence_sentence": "..."}]}. '
    "Include ALL distinct findings from the abstract. If the abstract "
    "contains no finding about what affects regulatory element activity, "
    'return {"findings": []}.\n\n'
    "Abstract: {abstract}\n\nJSON:"
)


def build_activity_messages(abstract_text: str) -> list[dict]:
    """构造 Chat 消息格式（供 Qwen2.5-Instruct 的 apply_chat_template）。"""
    return [
        {"role": "system", "content": ACTIVITY_SYSTEM_PROMPT},
        {
            "role": "user",
            "content": ACTIVITY_USER_TEMPLATE.replace("{abstract}", abstract_text),
        },
    ]


def build_activity_prompt(abstract_text: str) -> str:
    """构造字符串拼接提示词（非 chat 模板时使用）。"""
    return (
        f"{ACTIVITY_SYSTEM_PROMPT}\n"
        f"{ACTIVITY_USER_TEMPLATE.replace('{abstract}', abstract_text)}"
    )
