import sys
import json
import time
import torch
from transformers import AutoModelForCausalLM, AutoTokenizer

# 模型路径，相对于项目根目录（或直接用绝对路径）
MODEL_PATH = "./Qwen2.5-7B-Instruct"  # 如果脚本在根目录运行，这个相对路径有效

def main():
    # 可选的 max_memory 参数（JSON 字符串，键为 GPU 编号）
    max_memory = None
    if len(sys.argv) > 1:
        try:
            raw = json.loads(sys.argv[1])
            max_memory = {int(k): v for k, v in raw.items()}
        except (json.JSONDecodeError, ValueError):
            print(f"无效的 max_memory 参数: {sys.argv[1]}")
            sys.exit(1)
    print(f"max_memory: {max_memory}")

    print("=" * 50)
    print("Loading tokenizer...")
    tokenizer = AutoTokenizer.from_pretrained(
        MODEL_PATH,
        local_files_only=True,
        trust_remote_code=True
    )
    print("Tokenizer loaded.")

    print("Loading model (this may take 10-20 seconds)...")
    start = time.time()
    model = AutoModelForCausalLM.from_pretrained(
        MODEL_PATH,
        dtype=torch.float16,
        device_map="auto",          # 自动分配到多张 GPU
        local_files_only=True,
        trust_remote_code=True,
        max_memory=max_memory
    )
    load_time = time.time() - start
    print(f"Model loaded in {load_time:.2f}s")
    print(f"Model device: {model.device}")
    print(f"Model dtype: {model.dtype}")

    # 测试简单 prompt（用 chat 模板）
    messages = [
        {"role": "system", "content": "You are a helpful assistant."},
        {"role": "user", "content": "列举K562细胞中高活性序列富集的3个基序"},
    ]
    prompt = tokenizer.apply_chat_template(messages, tokenize=False, add_generation_prompt=True)
    print(f"\nInput: {prompt[:200]}")
    inputs = tokenizer(prompt, return_tensors="pt").to(model.device)

    print("Generating...")
    start = time.time()
    with torch.no_grad():
        outputs = model.generate(
            **inputs,
            max_new_tokens=100,
            do_sample=False,
            temperature=1.0
        )
    gen_time = time.time() - start
    response = tokenizer.decode(outputs[0], skip_special_tokens=True)
    print(f"Generation took {gen_time:.2f}s")
    print(f"Output:\n{response}")

    # 显存使用情况
    if torch.cuda.is_available():
        for i in range(torch.cuda.device_count()):
            alloc = torch.cuda.memory_allocated(i) / 1024**3
            reserved = torch.cuda.memory_reserved(i) / 1024**3
            print(f"GPU {i}: allocated {alloc:.2f} GB, reserved {reserved:.2f} GB")

    print("=" * 50)
    print("Test completed successfully!")

if __name__ == "__main__":
    main()