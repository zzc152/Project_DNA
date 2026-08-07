"""深入调试：测试不同提示词变体，找出模型输出空 JSON 的原因。

用法:
    python scripts/debug/debug_prompts.py [max_memory_json]
"""
import json
import sys

import torch
from transformers import AutoModelForCausalLM, AutoTokenizer

MODEL_PATH = "./Qwen2.5-7B-Instruct"

ABSTRACT = (
    "Transcription factors (TFs) recognizing DNA motifs within regulatory regions drive cell identity. "
    "Despite recent advances, their specificity remains incompletely understood. Here, we address this "
    "by contrasting two TFs, Neurogenin-2 (NGN2) and MyoD1, which recognize ubiquitous E-box motifs yet "
    "drive distinct cell fates toward neurons and muscles, respectively. Upon induction in mouse "
    "embryonic stem cells, we mapped genome-wide binding of NGN2 and MyoD1 over time."
)

# 各种提示词变体
VARIANTS = {
    "v1_current": (
        "You are a biomedical information extraction assistant. "
        "Extract relationships about human gene regulation from the given abstract. "
        "Output ONLY a JSON object with keys: tf (list), gene (list), "
        "motif (list), disease (list), relation (string, max 2 sentences). "
        "If no relevant info, return empty dict {}.",
        'Example:\nAbstract: "GATA1 binds to WGATAR motif and activates EPOR. Mutations cause anemia."\n'
        'Output:\n{"tf": ["GATA1"], "gene": ["EPOR"], "motif": ["WGATAR"], "disease": ["anemia"], '
        '"relation": "GATA1 activates EPOR via WGATAR binding; mutations cause anemia."}\n\n'
        "Now extract from the following abstract and output ONLY the JSON:\nAbstract: " + ABSTRACT,
    ),
    "v2_no_empty": (
        "You are a biomedical information extraction assistant. "
        "Extract transcription factor (TF), target gene, DNA motif, and disease entities from the abstract. "
        "Return ONLY a JSON object with keys: tf (list), gene (list), motif (list), disease (list), "
        "relation (string, max 2 sentences). Always extract as much as possible.",
        "Abstract: " + ABSTRACT + "\n\nOutput JSON:",
    ),
    "v3_chinese": (
        "你是一名生物医学信息抽取助手。请从给定的摘要中抽取转录因子(TF)、靶基因、DNA结合基序和疾病实体。"
        "只输出一个 JSON 对象，包含以下键：tf (列表), gene (列表), motif (列表), disease (列表), relation (字符串，最多两句话)。"
        "请尽可能多地抽取实体，不要遗漏。",
        "摘要：" + ABSTRACT + "\n\n请输出 JSON：",
    ),
    "v4_direct_question": (
        "You are a helpful biology assistant.",
        'Read the abstract below and list every transcription factor, gene, DNA motif, and disease mentioned. '
        'Respond ONLY in this JSON format: {"tf": [...], "gene": [...], "motif": [...], "disease": [...], '
        '"relation": "..."}. If an entity type has no matches, use an empty list [].\n\n'
        "Abstract: " + ABSTRACT + "\n\nJSON:",
    ),
}


def load_model(max_memory):
    tokenizer = AutoTokenizer.from_pretrained(MODEL_PATH, local_files_only=True, trust_remote_code=True)
    if tokenizer.pad_token is None:
        tokenizer.pad_token = tokenizer.eos_token
    model = AutoModelForCausalLM.from_pretrained(
        MODEL_PATH,
        dtype=torch.float16,
        device_map="auto",
        local_files_only=True,
        trust_remote_code=True,
        max_memory=max_memory,
    )
    model.eval()
    return tokenizer, model


def main():
    max_memory = None
    if len(sys.argv) > 1:
        try:
            raw = json.loads(sys.argv[1])
            max_memory = {int(k): v for k, v in raw.items()}
        except (json.JSONDecodeError, ValueError):
            print(f"无效 max_memory: {sys.argv[1]}")
            sys.exit(1)

    tokenizer, model = load_model(max_memory)
    print("模型加载完成\n")

    for name, (system, user) in VARIANTS.items():
        messages = [
            {"role": "system", "content": system},
            {"role": "user", "content": user},
        ]
        prompt = tokenizer.apply_chat_template(messages, tokenize=False, add_generation_prompt=True)
        inputs = tokenizer(prompt, return_tensors="pt").to(model.device)
        with torch.no_grad():
            outputs = model.generate(**inputs, max_new_tokens=200, do_sample=False)
        new_tokens = outputs[0][inputs["input_ids"].shape[1]:]
        text = tokenizer.decode(new_tokens, skip_special_tokens=True)
        print("=" * 60)
        print(f"[{name}] 生成 {len(new_tokens)} tokens")
        print("OUTPUT:")
        print(text)
        print()


if __name__ == "__main__":
    main()
