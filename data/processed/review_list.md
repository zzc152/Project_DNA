# 知识库验证复核清单（V2）

> 生成时间：2026-08-06 | 数据源：`data/processed/knowledge_base_clean.jsonl`（631 条）
> 验证：Qwen2.5-7B 三分类 + 推理链自洽性 | 问题记录：**33 条**（unsupported 31 + unclear 2，其中不自洽 10）
> 详细 CSV：远程 `data/processed/review_list.csv`（UTF-8 BOM，Excel 可直接打开）

## 汇总

| 建议操作 | 数量 | 说明 |
|---|---|---|
| 保留（Qwen 误判） | 6 | loss-of-function 逻辑模型未转过弯，claim 实际成立 |
| 修正 | 8 | 方向/主客体反 5 条 + 细胞系标错 3 条 |
| 修正 chain | 4 | design_rule 推理链矛盾/缺陷 |
| 复核降级/改写 | 13 | 过度泛化/证据间接，Qwen 严格但可辩护 |
| 剔除 | 2 | 证据完全不相关 |

---

## A. 保留（Qwen 误判，6 条）✅

| idx | factor | effect | 说明 |
|---|---|---|---|
| #86 | H3K27me3 | decreases | 证据"甲基化耗竭→超激活"⇒存在时抑制 |
| #109 | BCL11A | no_effect | "不结合 γ-globin 启动子"恰好支持 no_effect |
| #137 | ETO2 | decreases | "ETO2 loss→活性↑"⇒ETO2 存在→降低 |
| #286 | ERRα | decreases | "depletion→STING 转录↑"⇒存在→抑制 |
| #272 | SP3 | required_for | 结合 GC box 参与 coronin-1 表达，支持必需 |
| #418 | DNA methylation | decreases | 抑制 MYB 表达→降低 enhancer 活性 |

## B. 修正（真问题，8 条）⚠️

| idx | factor | 类别 | 问题 |
|---|---|---|---|
| #11 | enhancer | 方向/主客体反 | 证据"promoter 对 enhancer 必需"，claim 写成反了 |
| #236 | DNMT1 | 方向反 | 证据 hypermethylation...reducing（抑制），claim 写 increases |
| #250 | NFIB | 方向反+细胞系错 | 敲低→CDON 上调⇒NFIB 抑制；且 SH-SY5Y 写成 SK-N-SH |
| #258 | HDAC2 | 方向反 | 敲低→hSVCT2 表达↑⇒HDAC2 抑制 |
| #399 | Ghd7 and Ehd1 | 方向/主客体反 | OsDDE9 调节 Ghd7/Ehd1，被写成反了 |
| #224 | O-GlcNAc transferase (OGT) | 细胞系标错 | SH-SY5Y 写成 SK-N-SH |
| #226 | ATF4 | 细胞系标错 | 小鼠 brain 组织写成 SK-N-SH |
| #254 | YY1 | 细胞系标错 | SH-SY5Y 写成 SK-N-SH |

## C. 修正 chain（design_rule，4 条）⚠️

| idx | factor | 问题 |
|---|---|---|
| #57 | DNA methylation | chain 写"原始发现: increases"，但建议是"避免(负相关)"，方向自相矛盾 |
| #425 | H3 acetylation K4/K9 | chain 未明确"应避免引入"，末步"需功能实验"削弱结论 |
| #447 | acetylating and de-methylating agents | chain 第二步缺活性关系直接支持 |
| #449 | acetylating and de-methylating agents | chain 第二三步缺直接证据支持 |

## D. 复核降级/改写（13 条）🔍

| idx | factor | effect | 问题 |
|---|---|---|---|
| #14 | EGR-1 binding motif | required_for | CSF1R 激活→泛化 enhancer-promoter interaction |
| #28 | targeted activation of this regulatory element | increases | AML 分化→泛化 CRE 活性 |
| #34 | FOXP4 | increases | 转录激活基因→泛化 promoter 活性 |
| #43 | SPI1, GABPB and STAT3 | modulates | QTL 影响结合→泛化 TF 调节 enhancer |
| #54 | LEF1 | required_for | 敲低增强炎症/ROS→未直接说明 promoter 必需 |
| #160 | DNMT3B | decreases | 抑制 DNMT3B 水平→未直接说明降 promoter 活性 |
| #170 | OGA | decreases | GATA 突变减互作→未直接说明 OGA 降活性 |
| #199 | ALKBH5 | decreases | 去甲基酶沉默降转录本→未直接说明 |
| #292 | DNA methylation | decreases | 上调基因高甲基化，方向需人工核定 |
| #323 | chromatin accessibility | modulates | 高可及性位点→未直接说明调节 |
| #375 | TSKU | modulates | 维生素 D 代谢酶→跳跃到 enhancer |
| #381 | ODC1 | increases | 肿瘤上调→未直接说明增强 promoter |
| #470 | transcription factor binding motifs | modulates | caQTL 破坏位点→泛化 |
| #479 | NF-Y | required_for | MAFA 启动子→泛化所有 promoter |

> ⚠️ 注：D 类实际 14 条（#14/28/34/43/54/160/170/199/292/323/375/381/470/479），与汇总表中 13 条有出入，以 CSV 为准。

## E. 剔除（证据不相关，2 条）❌

| idx | factor | 问题 |
|---|---|---|
| #193 | rs2294510 | 证据是脂质积累(NAFLD)，claim 说增强 promoter，无关 |
| #367 | NCOA4 | 证据是铁死亡敏感性，claim 说增强 promoter，无关 |

---

## 后续操作建议

1. **自动修正 B、C 类（12 条）**：在 `curate_knowledge_base.py` 或直接改 knowledge_base_clean.jsonl，按上述问题说明修正方向/细胞系/chain
2. **剔除 E 类（2 条）**：删除 #193、#367
3. **D 类人工复核**：逐条对照原文决定降级（modulates）或改写 claim 措辞
4. **重跑验证**：修正后重跑 `validate_knowledge_base.py` 确认 supported 率提升
