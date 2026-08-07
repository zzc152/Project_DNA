"""BioExtractor：基于本地 Qwen2.5-7B-Instruct 的 DNA 调控知识抽取器。

封装模型加载、批量推理与 JSON 输出解析。所有模型加载均使用
``local_files_only=True``，不访问网络。
"""

import json
import logging
import re
from typing import Optional

import torch
from transformers import AutoModelForCausalLM, AutoTokenizer

from .prompts import build_messages, build_prompt

logger = logging.getLogger(__name__)

# 期望输出的字段
EXPECTED_FIELDS = ("tf", "gene", "motif", "disease", "relation")


def parse_json_response(text: str) -> Optional[dict]:
    """从模型输出中提取并解析 JSON。

    兼容以下几种常见情况：
    - 纯 JSON 对象
    - 被 ```json ... ``` 或 ``` ... ``` 代码块包裹
    - 前后有其他文字，但包含完整的 {...} 片段
    解析失败返回 None。
    """
    if not text:
        return None
    text = text.strip()

    # 1) 直接解析
    try:
        return json.loads(text)
    except json.JSONDecodeError:
        pass

    # 2) 提取 ```json ... ``` / ``` ... ``` 代码块
    m = re.search(r"```(?:json)?\s*(.*?)\s*```", text, re.DOTALL)
    if m:
        try:
            return json.loads(m.group(1))
        except json.JSONDecodeError:
            pass

    # 3) 提取第一个 { 到最后一个 } 之间的内容
    start, end = text.find("{"), text.rfind("}")
    if start != -1 and end > start:
        try:
            return json.loads(text[start : end + 1])
        except json.JSONDecodeError:
            pass

    return None


class BioExtractor:
    """从摘要批量抽取调控知识的抽取器。"""

    def __init__(
        self,
        model_path: str,
        device_map: str = "auto",
        dtype=torch.float16,
        local_files_only: bool = True,
        max_memory: Optional[dict] = None,
    ):
        logger.info("加载 tokenizer: %s", model_path)
        self.tokenizer = AutoTokenizer.from_pretrained(
            model_path,
            local_files_only=local_files_only,
            trust_remote_code=True,
        )
        if self.tokenizer.pad_token is None:
            self.tokenizer.pad_token = self.tokenizer.eos_token
        # decoder-only 模型批量生成必须 left-padding，
        # 否则模型会从右侧 pad token 处开始乱续写（复读提示词片段）
        self.tokenizer.padding_side = "left"

        logger.info(
            "加载模型: %s (dtype=%s, device_map=%s, max_memory=%s)",
            model_path,
            dtype,
            device_map,
            max_memory,
        )
        self.model = AutoModelForCausalLM.from_pretrained(
            model_path,
            dtype=dtype,
            device_map=device_map,
            local_files_only=local_files_only,
            trust_remote_code=True,
            max_memory=max_memory,
        )
        self.model.eval()

        n_gpus = torch.cuda.device_count() if torch.cuda.is_available() else 0
        logger.info("模型就绪，可用 GPU 数: %d", n_gpus)

    @torch.no_grad()
    def extract_batch(
        self,
        items: list[dict],
        batch_size: int = 8,
        max_new_tokens: int = 256,
        temperature: float = 0.1,
        use_chat_template: bool = True,
        prompt_builder=None,
        fields: tuple = EXPECTED_FIELDS,
    ) -> list[dict]:
        """批量抽取。

        Args:
            items: 输入记录列表，每个元素为 {"pmid": ..., "abstract": ...}
            batch_size: 推理批大小
            max_new_tokens: 最大生成 token 数
            temperature: 采样温度（任务要求 0.1）
            use_chat_template: True 使用 apply_chat_template，否则用字符串拼接
            prompt_builder: 可选。自定义 prompt 构造器
                （callable: abstract -> messages 或 prompt 字符串）。
                默认使用 prompts.build_messages / build_prompt。
            fields: 可选。输出字段元组（默认 EXPECTED_FIELDS）。
                对每个字段，若解析结果含该 key 则保留原值，
                否则给默认值（relation 给空串，其余给空列表）。

        Returns:
            结果列表，每个元素包含 pmid、抽取字段以及解析状态与原始输出。
            单条 JSON 解析失败不中断，记录 raw_output 并标记 parsed=False。
        """
        results: list[dict] = []
        total = len(items)
        for i in range(0, total, batch_size):
            batch = items[i : i + batch_size]
            results.extend(
                self._extract_one_batch(
                    batch,
                    max_new_tokens,
                    temperature,
                    use_chat_template,
                    prompt_builder=prompt_builder,
                    fields=fields,
                )
            )
            logger.info("进度: %d/%d", min(i + batch_size, total), total)

        return results

    @torch.no_grad()
    def _extract_one_batch(
        self,
        batch: list[dict],
        max_new_tokens: int = 256,
        temperature: float = 0.1,
        use_chat_template: bool = True,
        prompt_builder=None,
        fields: tuple = EXPECTED_FIELDS,
    ) -> list[dict]:
        """抽取单个批次，返回该批结果列表。供 OOM 降批重试时逐批调用。"""
        abstracts = [it.get("abstract", "") for it in batch]

        if prompt_builder is not None:
            prompts = []
            for a in abstracts:
                built = prompt_builder(a)
                # 支持两种返回：字符串（直接使用）或消息列表（走 chat template）
                if isinstance(built, str):
                    prompts.append(built)
                else:
                    prompts.append(
                        self.tokenizer.apply_chat_template(
                            built, tokenize=False, add_generation_prompt=True
                        )
                    )
        elif use_chat_template:
            prompts = [
                self.tokenizer.apply_chat_template(
                    build_messages(a),
                    tokenize=False,
                    add_generation_prompt=True,
                )
                for a in abstracts
            ]
        else:
            prompts = [build_prompt(a) for a in abstracts]

        inputs = self.tokenizer(
            prompts,
            return_tensors="pt",
            padding=True,
            truncation=True,
            max_length=2048,
        ).to(self.model.device)

        outputs = self.model.generate(
            **inputs,
            max_new_tokens=max_new_tokens,
            do_sample=True,
            temperature=temperature,
            top_p=0.9,
            pad_token_id=self.tokenizer.pad_token_id,
            eos_token_id=self.tokenizer.eos_token_id,
        )

        # 只保留新生成的 token
        new_tokens = outputs[:, inputs["input_ids"].shape[1] :]
        texts = self.tokenizer.batch_decode(new_tokens, skip_special_tokens=True)

        results: list[dict] = []
        for item, text in zip(batch, texts):
            parsed = parse_json_response(text)
            result: dict = {
                "pmid": item.get("pmid"),
                "parsed": parsed is not None,
                "raw_output": text,
            }
            # 字段默认值：relation/mechanism 等字符串字段给空串，其余给空列表
            for key in fields:
                result[key] = "" if key in ("relation", "mechanism") else []
            if parsed is not None:
                for key in fields:
                    val = parsed.get(key)
                    if val is None:
                        continue
                    if isinstance(val, list):
                        result[key] = val
                    elif isinstance(val, str):
                        result[key] = val
                    elif isinstance(val, dict):
                        result[key] = val
                    else:
                        result[key] = [val] if val else result[key]
            else:
                logger.warning(
                    "JSON 解析失败 (pmid=%s): %s...",
                    item.get("pmid"),
                    text[:120].replace("\n", " "),
                )
            results.append(result)

        return results
