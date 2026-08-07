import torch
from transformers import AutoTokenizer
from awq import AutoAWQForCausalLM

model_path = "/workspace/zzc/BioDesign-Agent/models/Qwen2.5-72B-AWQ"  # 确认你的实际路径

tokenizer = AutoTokenizer.from_pretrained(model_path, local_files_only=True)

# 方法1：自动分配（推荐）
model = AutoAWQForCausalLM.from_quantized(
    model_path,
    device_map="auto",          # 让 accelerate 自动分配到单卡或多卡
    fuse_layers=True,
    local_files_only=True
)

# 方法2：强制指定单卡（如果你确定只用 GPU 0）
# model = AutoAWQForCausalLM.from_quantized(
#     model_path,
#     device_map={'': 'cuda:0'},   # 这样写
#     fuse_layers=True,
#     local_files_only=True
# )

prompt = "用中文解释什么是转录因子"
inputs = tokenizer(prompt, return_tensors="pt").to("cuda")
outputs = model.generate(**inputs, max_new_tokens=100)
print(tokenizer.decode(outputs[0], skip_special_tokens=True))