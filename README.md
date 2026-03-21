# CRC/UC 单细胞 RNA-seq 整合分析

**项目**：跨疾病阶段 CRC 免疫微环境全景图谱（scVI 整合）  
**数据**：5 个训练数据集 + 1 个外部验证集，共 ~240k 细胞  
**疾病**：结直肠癌（CRC）/ 溃疡性结肠炎（UC）

---

## 数据集

| 数据集 | 疾病 | 细胞数 | 来源 | 角色 |
|--------|------|--------|------|------|
| GSE132257 | CRC | 18,409 | GEO | 训练 |
| GSE132465 | CRC | 63,689 | GEO | 训练 |
| cellxgene | CRC | 26,945 | CellxGene | 训练（仅含恶性细胞+上皮亚型，含转移灶，**无免疫细胞**）|
| GSE125527 | UC  | 84,129 | GEO | 训练 |
| GSE150115 | UC  | 3,677  | GEO | 训练 |
| GSE146771 | CRC | 43,817 | GEO | **外部验证**（免疫细胞，TPM格式）|

训练集合并：**196,849 细胞 × 9,009 基因**

---

## 分析流程

```
raw data (data/)
    │
    ├─ convert_*.py          # Step 1: 各数据集转为 h5ad 格式
    │
    ├─ qc_cellxgene.py       # Step 2: 质量控制（各数据集 QC）
    ├─ filter_gse150115.py   #         GSE150115 特殊过滤
    │
    ├─ remerge_with_unified_annotations.py  # Step 3: 合并5数据集 → merged_5datasets_unified.h5ad（已弃用）
    ├─ rebuild_annotations.py               #         重建细胞类型注释
    │
    ├─ train_scvi_optimized.py  # Step 4: scVI 批次校正训练
    │                           #   → scvi_model_optimized/  (trained model)
    │                           #   → merged_annotated.h5ad  (主数据集)
    │
    ├─ annotate_subclusters.py  # Step 5: KNN 反向注释免疫亚群
    │                           #   以 GSE146771 Sub_Cluster 为参考
    │                           #   → annotation_results/
    │
    ├─ validate_external.py     # Step 6: 外部验证（GSE146771 投影）
    │
    ├─ retrain_scvi_6datasets.py  # Step 7: 6数据集 UMAP（GSE146771投影）
    │                              #   → scvi_results_6datasets/merged_6datasets_integrated.h5ad
    │
    └─ figure1_panel_*.py       # Step 8: Figure 1 各面板生成
                                #   → figures/phase1/
```

---

## 目录结构

```
changdao/
│
├── data/                          # 原始数据（只读）
│   ├── GSE132257_*.txt.gz
│   ├── GSE132465_*.txt.gz
│   ├── GSE125527_extracted/       # 解压的 UMI 表格
│   ├── GSE150115_extracted/
│   ├── GSE146771_*.txt.gz         # 外部验证（TPM格式）
│   └── cellxgene.h5ad
│
├── processed_data/                # 预处理后的数据
│   ├── GSE132257.h5ad             # 各数据集转换结果
│   ├── GSE132465.h5ad
│   ├── GSE125527.h5ad
│   ├── GSE150115.h5ad
│   ├── cellxgene_raw.h5ad
│   ├── GSE146771_tpm.h5ad         # 外部验证集（log1p TPM）
│   └── merged_annotated.h5ad      # ★ 主数据集（196,849 cells × 9,009 genes）
│                                  #   obsm: X_scvi (30-dim), X_umap (2-dim)
│                                  #   obs: batch, disease, unified_cell_type,
│                                  #        cell_type_grouped, unified_cell_subtype,
│                                  #        subcluster, subcluster_conf
│
├── scvi_model_optimized/          # ★ 训练好的 scVI 模型
│   └── model.pt                   #   架构: 9009 genes → 2 layers → 30 latent dims
│                                  #   训练参数: lr=1e-3, batch_size=256, early_stop
│
├── scvi_results_6datasets/        # 6数据集整合结果
│   └── merged_6datasets_integrated.h5ad  # 240,666 cells (5训练+GSE146771投影)
│                                          # obsm: X_scvi, X_umap
│
├── qc_results/                    # QC 报告图（各数据集）
├── annotation_results/            # 免疫亚群注释结果
│   ├── umap_subclusters.png
│   ├── subcluster_assignments.csv
│   └── consistency_heatmap.png
├── validation_results/            # 外部验证结果
│
├── figures/
│   └── phase1/                    # Figure 1 所有面板（PNG + PDF）
│       ├── figure1_panelA_batch.*         # UMAP 按数据集着色
│       ├── figure1_panelB_disease.*       # UMAP 按疾病着色（CRC/UC）
│       ├── figure1_panel_C_composition.*  # 细胞类型组成堆叠柱状图
│       ├── figure1_panel_D_dotplot.*      # Marker 基因 dotplot
│       └── figure1_panel_E_lisi.*         # LISI 整合质量评估
│
├── convert_gse132257.py           # 数据转换脚本（各数据集一个）
├── convert_gse132465.py
├── convert_gse125527_correct.py
├── convert_gse150115.py
├── convert_cellxgene.py
├── convert_gse146771.py
├── filter_gse150115.py            # GSE150115 质量过滤
├── qc_cellxgene.py                # CellxGene QC
├── remerge_with_unified_annotations.py  # 数据合并与注释统一
├── rebuild_annotations.py         # 注释重建
├── train_scvi_optimized.py        # scVI 训练（推荐入口）
├── retrain_scvi_6datasets.py      # 6数据集 UMAP 生成
├── annotate_subclusters.py        # 免疫亚群 KNN 注释
├── validate_external.py           # 外部验证
├── figure1_panel_C.py             # Figure 1C 生成脚本
├── figure1_panel_D.py             # Figure 1D 生成脚本
├── figure1_panel_E.py             # Figure 1E 生成脚本
└── requirements.txt               # Python 依赖
```

---

## 主数据集说明（merged_annotated.h5ad）

| 字段 | 说明 |
|------|------|
| `obs['batch']` | 数据集来源：GSE132257 / GSE132465 / GSE125527 / GSE150115 / cellxgene |
| `obs['disease']` | CRC 或 UC |
| `obs['unified_cell_type']` | 细粒度细胞类型（17类，含上皮亚型）|
| `obs['cell_type_grouped']` | 粗粒度分组（9类：T cells / NK / B / Myeloids / Mast / Epithelial / Stromal / Malignant / Unknown）|
| `obs['subcluster']` | 免疫细胞亚群（以 GSE146771 38个亚群为参考，KNN预测）|
| `obs['subcluster_conf']` | KNN 预测置信度（仅免疫细胞有值）|
| `obsm['X_scvi']` | scVI 30维 latent embedding（批次校正后）|
| `obsm['X_umap']` | UMAP 2维坐标 |

### 细胞类型组成

| 类型 | 细胞数 | 比例 |
|------|--------|------|
| T cells | 77,927 | 39.6% |
| B cells | 36,987 | 18.8% |
| Epithelial | 33,806 | 17.2% |
| Malignant | 16,664 | 8.5% |
| Myeloids | 15,972 | 8.1% |
| Stromal | 7,224 | 3.7% |
| NK cells | 7,224 | 3.7% |
| Unknown | 570 | 0.3% |
| Mast cells | 475 | 0.2% |

---

## scVI 模型参数

| 参数 | 值 |
|------|----|
| 基因数 | 9,009 |
| 隐层数 | 2 |
| 潜在维度 | 30 |
| 基因似然 | Negative Binomial |
| batch_key | `batch`（5个数据集）|
| covariate | `disease`（CRC/UC）|
| 学习率 | 1e-3 |
| Batch Size | 256 |
| Early stopping | patience=15 |
| 实际训练轮数 | ~35 轮（早停）|

> **注意**：GSE146771 不参与模型训练，仅通过 encoder 投影到 latent space 作为外部验证。

---

## Figure 1 面板

| 面板 | 文件 | 内容 | 数据来源 |
|------|------|------|----------|
| A | figure1_panelA_batch.* | UMAP（6数据集，含GSE146771） | merged_6datasets_integrated.h5ad |
| B | figure1_panelB_disease.* | UMAP（CRC vs UC） | merged_6datasets_integrated.h5ad |
| C | figure1_panel_C_composition.* | 细胞类型组成堆叠柱状图 | merged_annotated.h5ad |
| D | figure1_panel_D_dotplot.* | Marker 基因 dotplot（44基因×17类型）| merged_annotated.h5ad |
| E | figure1_panel_E_lisi.* | LISI 评分（scVI vs Harmony vs PCA）| merged_annotated.h5ad |

---

## 运行环境

```bash
pip install -r requirements.txt
```

主要依赖：`scanpy`, `scvi-tools`, `harmonypy`, `scib-metrics`, `matplotlib`, `seaborn`, `scikit-learn`

Python 3.12，建议使用 venv（项目已有 `venv/` 目录）。

---

## 运行顺序（重建流程）

```bash
# 1. 数据转换（各数据集）
python convert_gse132257.py
python convert_gse132465.py
python convert_gse125527_correct.py
python convert_gse150115.py && python filter_gse150115.py
python convert_cellxgene.py
python convert_gse146771.py

# 2. 合并与注释
python remerge_with_unified_annotations.py
python rebuild_annotations.py

# 3. scVI 训练（~20分钟）
python train_scvi_optimized.py

# 4. 免疫亚群注释
python annotate_subclusters.py

# 5. 6数据集 UMAP（含外部验证投影）
python retrain_scvi_6datasets.py

# 6. Figure 1 生成
python figure1_panel_C.py
python figure1_panel_D.py
python figure1_panel_E.py   # ~15分钟（含 PCA/Harmony 计算）
```
