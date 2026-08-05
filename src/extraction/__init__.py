"""知识抽取模块：从 PubMed 摘要中提取 DNA 调控知识。"""
from .extractor import BioExtractor, parse_json_response

__all__ = ["BioExtractor", "parse_json_response"]
