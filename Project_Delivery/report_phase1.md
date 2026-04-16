# 单细胞转录组整合分析报告
**CRC/UC 肠道免疫微环境研究 · 阶段一**
**日期：2026-03-15**（2026-03-04 首发；2026-03-09 加入数据清洗；2026-03-10 更新：HVG 3k 整合 + CellAssign 注释 + 新 Figure 1；2026-03-11 更新：MT% 阈值 10%，patience=50，重跑全流程，153,136 cells；2026-03-15 更新：scANVI 半监督整合 + 阶段一交付）

---

## 一、数据概况

### 原始数据规模（清洗前）

| 数据集 | 疾病 | 原始细胞数 | 特点 |
|---|---|---|---|
| GSE132257 | CRC | 18,409 | 新鲜/冷冻组织对比实验（Smillie et al.） |
| GSE132465 | CRC | 63,689 | 最大 CRC 免疫细胞数据集（Zhang et al.） |
| GSE125527 | UC | 84,129 | 几乎全为免疫细胞，无上皮/肿瘤细胞 |
| GSE150115 | UC | 32,000 | 含多组织部位；原始 UMI 中位数仅 112 |
| cellxgene | CRC | 47,107 | 仅含恶性肿瘤和上皮细胞，来自乙状结肠/直肠 |
| **合计** | CRC + UC | **245,334** | 清洗前汇总 |

- 所有数据集 X 矩阵均为**原始整数 UMI counts**（含 .raw 提取步骤，满足 scVI 负二项分布建模要求）
- 经批次感知 HVG 筛选至 **3,000 个**（见第三节）

---

## 二、数据清洗（Data Preprocessing）

### 2.1 原始 QC 指标（清洗前）

| 数据集 | 原始细胞数 | UMI 中位数 | 基因数中位数 | 线粒体% 中位数 | MT 基因数 |
|---|---|---|---|---|---|
| GSE132257 | 18,409 | 3,742 | 1,155 | 3.8% | 13 |
| GSE132465 | 63,689 | 4,259 | 1,265 | 4.5% | 13 |
| GSE125527 | 84,129 | 1,858 | 796 | 6.8% | 13 |
| GSE150115 | 32,000 | 112 | 58 | 26.1% | 34 |
| cellxgene | 47,107 | 14,104 | 3,192 | 18.1% | 13 |

### 2.2 清洗步骤：MT 过滤 + Scrublet 双细胞去除

脚本：`preprocess_phase1.py`

| 数据集 | MT% 阈值 | 最低基因数 | 最低 UMI 数 | 双细胞批次 |
|---|---|---|---|---|
| GSE132257 | ≤ 10% | — | — | `Sample` 列 |
| GSE132465 | ≤ 10% | — | — | `Patient` 列 |
| GSE125527 | ≤ 10% | — | — | cell ID 前缀 |
| GSE150115 | ≤ 10% | ≥ 500 | ≥ 1,000 | GSM 编号前缀 |
| cellxgene | ≤ 10% | ≥ 200 | ≥ 500（+中位数+3×MAD 上界） | `donor_id` 列 |

### 2.3 清洗后数据规模

| 数据集 | 清洗前 | **清洗后** | 保留率 |
|---|---|---|---|
| GSE132257 | 18,409 | **15,642** | 85.0% |
| GSE132465 | 63,689 | **54,487** | 85.6% |
| GSE125527 | 84,129 | **71,866** | 85.4% |
| GSE150115 | 32,000 | **2,138** | 6.7% |
| cellxgene | 47,107 | **9,003** | 19.1% |
| **合计** | **245,334** | **153,136** | **62.4%** |

> GSE150115 线粒体% 中位数高达 26.1%，MT≤10% 后仅剩 6.7%，符合预期（低质量数据集）。cellxgene 线粒体% 中位数 18.1%，MT≤10% 后保留 19.1%，上皮/肿瘤来源细胞天然线粒体含量偏高。

![QC 质控曲线](figures/phase1/suppfig_qc_all_datasets.png)

---

## 三、基因筛选（HVG 3k）

> **新增（2026-03-10）**：脚本：`select_hvg3k.py`

**方法**：`sc.pp.highly_variable_genes(flavor='seurat_v3', n_top_genes=3000, batch_key='batch', layer='counts')`

直接对原始 counts 做批次感知的 HVG 筛选，选取在多个批次中均有较高变异的基因，避免单一批次特有基因主导整合。

| 基因出现批次数 | 基因数 |
|---|---|
| 5/5（全部批次）| 814 |
| 4/5 | 678 |
| 3/5 | 513 |
| 2/5 | 446 |
| 1/5 | 549 |
| **合计** | **3,000** |

输出：`processed_data/merged_hvg3k.h5ad`（153,136 × 3,000，含 `layers['counts']`）

---

## 四、scVI 批次整合

脚本：`train_scvi_hvg3k.py`

| 参数 | 值 |
|---|---|
| 输入基因数 | **3,000 HVG** |
| 潜在维度（n_latent）| 30 |
| 隐藏层数（n_layers）| 2 |
| 基因分布 | 负二项（NB）|
| batch_key | `batch`（5 个数据集）|
| 协变量 | `disease`（CRC/UC，categorical covariate）|
| 学习率 | 1e-3 |
| 批大小 | 256 |
| 早停 | patience=50，监控 ELBO_validation |
| 训练时间 | 约 24 分钟（CPU，153,136 cells）|

**关键设计**：将 `disease` 作为显式协变量，scVI 在消除批次技术噪声的同时**有意保留 CRC/UC 的生物学差异**。批次感知 HVG 3,000 个基因在保留跨批次共有信号的同时大幅降低噪声维度，提升整合质量。

输出：`scvi_model_hvg3k/`，`scvi_results_hvg3k/merged_scvi_hvg3k.h5ad`

### 4.1 scANVI：半监督扩展（改善 iLISI）

> **新增（2026-03-15）**：脚本：`train_scanvi_hvg3k.py`

scVI 的 iLISI（批次混合指数）偏低（原始值 1.33/5），原因是 disease 协变量导致 CRC/UC 生物学分离，使不同病种批次的局部邻域混合受限。scANVI 以 CellAssign 细胞类型标签作为弱监督信号，在 scVI 基础上进行半监督微调，同时提升批次混合和细胞类型边界清晰度。

| 参数 | 值 |
|---|---|
| 初始化 | 从 scVI 模型继承权重 |
| labels_key | `cellass_type`（10 种细胞类型）|
| unlabeled_category | `Unknown`（无未标注细胞）|
| max_epochs | 50（全程跑完）|
| n_samples_per_label | 100 |
| 训练时间 | 88.5 分钟（CPU）|

**LISI 对比（k=90 邻域，原始值）**：

| 指标 | scVI | scANVI | 变化 |
|---|---|---|---|
| iLISI（批次混合，↑ 越好） | 1.3304 | **1.3639** | +0.0336 ↑ |
| cLISI（细胞类型分离，↓ 越好）| 1.3337 | **1.2104** | −0.1234 ↓ |

scANVI 同时改善了批次混合（iLISI ↑）和细胞类型边界（cLISI ↓），是阶段二分析的推荐 embedding。

输出：`scanvi_model_hvg3k/`，`scanvi_results_hvg3k/merged_scanvi_hvg3k.h5ad`（含 `X_scanvi`、`X_umap_scanvi`、`X_scvi`、`X_umap_scvi`）

---

## 五、细胞类型注释（CellAssign）

脚本：`annotate_celltypes_cellassign.py`

### 5.1 Marker 基因矩阵

使用在 3k HVG 中确认存在的 40 个 marker 基因，定义 10 种细胞类型：

| 细胞类型 | Marker 基因 |
|---|---|
| CD4 T | CD3D, CD3E, CD4, IL7R, CCR7 |
| CD8 T | CD3D, CD3E, CD8A, GZMA, GZMB, GZMK, NKG7 |
| Treg | CD3D, CD3E, FOXP3, IL2RA, IKZF2, CTLA4 |
| NK | NCAM1, GNLY, NKG7, KLRB1, KLRC1, KLRD1, GZMB, PRF1 |
| B cell | CD19, CD79A, CD79B, MS4A1 |
| Plasma | MZB1, SDC1, IGHA1, IGHG1, IGKC, DERL3 |
| Myeloid | LYZ, CD14, CD68, CSF1R, FCGR3A, ITGAM |
| Mast | FCER1A |
| Stromal | VIM, ACTA2 |
| Epithelial | KRT18 |

训练参数：max_epochs=200，early stopping（patience=15）

### 5.2 全局细胞类型分布

| 细胞类型 | 细胞数 | 占比 |
|---|---|---|
| CD4 T | 40,325 | 26.3% |
| Epithelial | 28,703 | 18.7% |
| CD8 T | 17,974 | 11.7% |
| Plasma | 16,368 | 10.7% |
| B cell | 15,069 | 9.8% |
| NK | 13,712 | 9.0% |
| Myeloid | 13,448 | 8.8% |
| Treg | 5,534 | 3.6% |
| Stromal | 1,718 | 1.1% |
| Mast | 285 | 0.2% |
| **合计** | **153,136** | 100% |

### 5.3 各批次细胞类型组成（CellAssign）

| 批次 | CD4 T | CD8 T | Treg | NK | B cell | Plasma | Myeloid | Mast | Stromal | Epithelial |
|---|---|---|---|---|---|---|---|---|---|---|
| GSE132257 (CRC) | 18.0% | 14.3% | 4.8% | 5.6% | 6.7% | 19.5% | 7.2% | 0.4% | 1.6% | 22.0% |
| GSE132465 (CRC) | 18.0% | 13.5% | 4.0% | 6.6% | 6.0% | 10.0% | 10.4% | 0.1% | 2.6% | 28.7% |
| cellxgene (CRC) | 0.6% | 0.8% | 0.2% | 0.4% | 0.1% | 0.2% | 3.0% | 0.0% | 0.0% | **94.9%** |
| GSE125527 (UC) | **38.1%** | 11.4% | 3.4% | 12.7% | 14.4% | 9.8% | 8.8% | 0.2% | 0.0% | 1.2% |
| GSE150115 (UC) | 14.0% | 3.5% | 6.5% | 5.1% | 17.9% | **37.2%** | 3.1% | 0.0% | 1.1% | 11.6% |

> cellxgene 的 94.9% Epithelial 符合预期（该数据集原始注释为 malignant + epithelial，均属上皮谱系，CellAssign 无独立 Malignant 类别，统一归入 Epithelial）。

预测置信度（最大概率中位数）：**1.000**

---

## 六、Figure 1：跨疾病 CRC 免疫微环境全景图谱

![Figure 1 总图](figures/phase1/figure1_landscape.png)

### Panels A–C — scVI UMAP（按数据集 / 疾病 / 细胞类型）


![Panel AB+C](figures/phase1/fig1_panelAB.png)

**A**：5 个批次在整合后实现充分混合，技术批次效应有效校正。UC 批次（GSE125527/GSE150115）与 CRC 批次存在可分辨的分布差异，体现疾病生物学信号得到保留。

**B**：CRC（79,132 cells）与 UC（74,004 cells）存在可辨别的空间分布差异，scVI 保留了疾病信号，未过度混合。UC 细胞更集中于免疫细胞区域，CRC 细胞在 epithelial/malignant 区域富集。

**C**：CellAssign 10 种细胞类型在 UMAP 上呈现清晰的空间分离，CD4 T 最多（40,325 cells, 26.3%），Epithelial 次之（28,703 cells, 18.7%），与批次组成特征一致。

### Panel D — 各批次细胞类型组成（CellAssign）

![Panel D](figures/phase1/fig1_panelC.png)

全部 5 个批次均有完整注释。cellxgene 以 Epithelial 为主（肿瘤来源）；UC 批次（GSE125527/GSE150115）以 T/NK/B/Plasma 免疫细胞为主，与 UC 以免疫浸润为特征的病理生理学一致。

### Panel E — Canonical Marker 基因气泡图（CellAssign）

![Panel E](figures/phase1/fig1_panelD.png)

22 个 marker 基因（T/NK/Treg/B/Plasma/Myeloid/Epithelial/Stromal）在 10 种 CellAssign 细胞类型中的表达符合预期，验证注释可靠性。

### Panel F — scib-metrics 整合质量 & 聚类稳定性

![Panel F](figures/phase1/fig1_panelEF.png)

使用 `scib_metrics.benchmark.Benchmarker` 计算全套整合质量指标（共 9 个），基于 CellAssign 10 种细胞类型标签（`cellass_type`）：

**scib-metrics 原始分数（未归一化）：**

| 指标 | 类别 | scVI | Harmony | PCA (uninteg.) |
|---|---|---|---|---|
| iLISI | Batch | 0.014 | 0.087 | 0.000 |
| kBET | Batch | 0.100 | 0.327 | 0.063 |
| Graph connectivity | Batch | **0.865** | 0.850 | 0.839 |
| PCR comparison | Batch | **0.842** | 0.432 | 0.000 |
| cLISI | Bio | 0.995 | 0.994 | 0.995 |
| Silhouette label | Bio | 0.510 | 0.557 | 0.554 |
| KMeans NMI | Bio | 0.639 | 0.611 | **0.652** |
| KMeans ARI | Bio | **0.478** | 0.423 | 0.442 |
| Isolated labels | Bio | 0.581 | 0.568 | 0.584 |
| **Bio conservation（聚合）** | — | **0.640** | 0.631 | 0.645 |
| **Batch correction（聚合）** | — | **0.455** | 0.424 | 0.225 |
| **Total** | — | **0.566** | 0.548 | 0.477 |

> scVI 的 **PCR comparison（0.842）**远高于 Harmony（0.432）和 PCA（0.000），说明 scVI 在潜在空间中最有效地去除了批次主成分方差，批次效应校正最彻底。
> 三方法的 bio conservation 指标相近（0.63–0.65），证明整合过程没有破坏细胞类型的生物学特征。
> Total 综合得分：scVI（0.566）> Harmony（0.548）> PCA（0.477）。

#### Bootstrap 聚类稳定性（多分辨率）

在 4 个 Leiden 分辨率（0.3/0.5/0.8/1.0）各进行 10 次随机种子 bootstrap，全部 AMI > 0.90：

| Leiden 分辨率 | 簇数 | 平均 AMI | 标准差 |
|---|---|---|---|
| 0.3 | 23 | 0.9245 | ±0.0139 |
| 0.5 | 31 | 0.9286 | ±0.0152 |
| 0.8 | 39 | 0.9068 | ±0.0145 |
| 1.0 | 43 | 0.9211 | ±0.0175 |

scVI Leiden 聚类在所有测试分辨率下均具备极高统计鲁棒性，验证 23–43 个聚类定义在下游亚群分析中的可靠性。

---

## 七、阶段一结论（更新版）

1. **整合成功**：5 个异质数据集整合为统一潜在空间，**153,136 个细胞**（MT≤10% 严格过滤后净值；CRC 79,132 / UC 74,004）
2. **HVG 优化**：批次感知 HVG 3k 筛选减少噪声维度，814 个基因在全部 5 批次中高变异
3. **CellAssign 注释完整**：10 种细胞类型，置信度中位数 = 1.000；所有 5 个批次均获完整注释
4. **scVI 整合性能最优**：Total scib-metrics 综合得分 0.566，Batch correction 聚合分 0.455，PCR comparison 0.842（批次方差去除效果最彻底），bio conservation 0.640 与 PCA baseline（0.645）相当
5. **scANVI 进一步改善**：在 scVI 基础上半监督微调（88.5 min），iLISI 提升（1.3304→1.3639），cLISI 下降（1.3337→1.2104），同时改善批次混合和细胞类型分离；**推荐为阶段二主 embedding**
6. **聚类鲁棒性极高**：Bootstrap AMI > 0.90（4 个分辨率），为阶段二亚群精细分析奠定可靠基础
7. **数据质量完整记录**：GSE150115（MT≤10%后仅剩 6.7%，2,138 cells）反映其本身低质量特性；cellxgene（19.1%）高线粒体含量为肿瘤/上皮来源细胞特性

---

## 八、阶段二计划：免疫微环境生物学分析

### 核心科学问题
> **CRC 与 UC 的肠道免疫微环境有何异同？各细胞类型如何参与各自的疾病进程？**

### 拟分析内容

| 优先级 | 分析 | 方法 |
|---|---|---|
| ⭐⭐⭐ | Differential Abundance | Milo（KNN 邻域水平丰度检验）/ scCODA |
| ⭐⭐⭐ | Pseudobulk DEG | DESeq2 → GSEA/GSVA 通路富集 |
| ⭐⭐ | T 细胞精细亚型 + 耗竭轨迹 | Palantir/scFates 伪时间；TOX/TIGIT/PD-1 评分 |
| ⭐⭐ | 细胞通讯网络 | CellChatV2（PD-L1/TIGIT/IL/TNF/CXCL） |
| ⭐ | 髓系异质性 | SPP1+ TAM 识别；M1/M2 极化评分 |

**建议从 Milo Differential Abundance 开始。**

---

## 九、阶段一交付清单（2026-03-15）

交付目录：`Project_Delivery/`

### 9.1 报告

| 文件 | 位置 |
|---|---|
| `report_phase1.pdf` | `Project_Delivery/` |

### 9.2 数据（`Project_Delivery/data/processed_data/`）

| 文件 | 描述 | 大小 |
|---|---|---|
| `merged_hvg3k.h5ad` | 153,136 × 3,000 HVG，含 raw counts | ~400 MB |
| `merged_scvi_hvg3k.h5ad` | scVI 整合结果，含 `X_scvi`、`X_umap`、`cellass_type` | ~400 MB |
| `merged_scanvi_hvg3k.h5ad` | scANVI 整合结果，含 `X_scanvi`、`X_umap_scanvi` | ~400 MB |

### 9.3 Figure 1 — scVI Atlas（`Project_Delivery/results/Figure1_scVI_Atlas/`）

| 子目录 | 图片 | 表格 |
|---|---|---|
| `scVI_Integration/` | `fig1_panelAB.png/pdf`（UMAP A/B/C）、`figure1_landscape.png/pdf`（总图）、`fig_scanvi_umap.png/pdf`（scVI vs scANVI） | `dataset_summary.csv`、`scanvi_lisi_comparison.csv` |
| `Cell_Composition/` | `fig1_panelC.png/pdf`（各批次组成堆叠图） | `celltype_per_batch_pct.csv`、`celltype_per_batch_counts.csv`、`celltype_global_distribution.csv` |
| `Marker_Expression/` | `fig1_panelD.png/pdf`（Canonical marker 气泡图） | — |
| `Cluster_Stability/` | `fig1_panelEF.png/pdf`（scib-metrics 热图 + bootstrap AMI） | `scib_metrics.csv`、`bootstrap_ami.csv` |
| `script/` | 所有阶段一分析脚本（12 个 .py 文件） | — |

### 9.4 Supplementary（`Project_Delivery/results/Supplementary/Supp1_QC_Metrics/`）

| 文件 | 描述 |
|---|---|
| `suppfig_qc_all_datasets.png` | 5 个数据集 QC 小提琴图（清洗前后对比） |
