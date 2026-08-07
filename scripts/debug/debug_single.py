"""单条调试：对比不同生成配置下模型的抽取输出质量。

用法:
    python scripts/debug/debug_single.py [max_memory_json]

定位批量抽取输出退化（复读指令/空 JSON）的原因。
"""
import json
import sys
import time

import torch
from transformers import AutoModelForCausalLM, AutoTokenizer

MODEL_PATH = "./Qwen2.5-7B-Instruct"

# 用 pmid=40780181（含 NGN2/MyoD1/E-box）作为调试样例
ABSTRACT = (
    "Transcription factors (TFs) recognizing DNA motifs within regulatory regions drive cell identity. "
    "Despite recent advances, their specificity remains incompletely understood. Here, we address this "
    "by contrasting two TFs, Neurogenin-2 (NGN2) and MyoD1, which recognize ubiquitous E-box motifs yet "
    "drive distinct cell fates toward neurons and muscles, respectively. Upon induction in mouse "
    "embryonic stem cells, we mapped genome-wide binding of NGN2 and MyoD1 over time."
)

SYSTEM_PROMPT = (
    "You are a biomedical information extraction assistant. "
    "Extract relationships about human gene regulation from the given abstract. "
    "Output ONLY a JSON object with keys: tf (list), gene (list), "
    "motif (list), disease (list), relation (string, max 2 sentences). "
    "If no relevant info, return empty dict {}."
)

FEW_SHOT = (
    'Example:\nAbstract: "GATA1 binds to WGATAR motif and activates EPOR. Mutations cause anemia."\n'
    'Output:\n{"tf": ["GATA1"], "gene": ["EPOR"], "motif": ["WGATAR"], "disease": ["anemia"], '
    '"relation": "GATA1 activates EPOR via WGATAR binding; mutations cause anemia."}'
)


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


def run(tokenizer, model, messages, label, **gen_kwargs):
    prompt = tokenizer.apply_chat_template(messages, tokenize=False, add_generation_prompt=True)
    inputs = tokenizer(prompt, return_tensors="pt").to(model.device)
    start = time.time()
    with torch.no_grad():
        outputs = model.generate(**inputs, **gen_kwargs)
    elapsed = time.time() - start
    new_tokens = outputs[0][inputs["input_ids"].shape[1]:]
    text = tokenizer.decode(new_tokens, skip_special_tokens=True)
    print("=" * 60)
    print(f"[{label}] 耗时 {elapsed:.1f}s, 生成 {len(new_tokens)} tokens")
    print("RAW OUTPUT:")
    print(text)
    print()


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
    print("模型加载完成")

    user_content = f"{FEW_SHOT}\nNow extract from the following abstract and output ONLY the JSON:\nAbstract: {ABSTRACT}"
    messages = [
        {"role": "system", "content": SYSTEM_PROMPT},
        {"role": "user", "content": user_content},
    ]

    # 配置1: 贪心解码（do_sample=False）
    run(tokenizer, model, messages, "greedy", max_new_tokens=200, do_sample=False)

    # 配置2: 采样 temperature=0.1 top_p=0.9（当前批量抽取用的）
    run(tokenizer, model, messages, "sample t=0.1", max_new_tokens=200, do_sample=True, temperature=0.1, top_p=0.9)

    # 配置3: 采样 temperature=0.7（标准）
    run(tokenizer, model, messages, "sample t=0.7", max_new_tokens=200, do_sample=True, temperature=0.7, top_p=0.9)


if __name__ == "__main__":
    main()
