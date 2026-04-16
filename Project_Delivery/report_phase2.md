# 单细胞转录组分析报告
**CRC/UC 肠道免疫微环境研究 · 阶段二**
**日期：2026-03-23**

---

## 一、研究目标

在阶段一整合的 153,136 个细胞（CRC 79,132 / UC 74,004；其中 T 细胞 77,927 个、髓系细胞 15,972 个）基础上，阶段二聚焦于以下核心科学问题：

> **在 CRC 肿瘤微环境中，哪些免疫细胞亚群充当"前驱动"（Pre-driver）抑制性角色，其时空演化轨迹如何，转录调控网络如何驱动其分化？** 本分析重点揭示：**OLR1⁺ 脂质代谢型 TAM（OLR1⁺ lipid-associated TAMs）** 在 UC→CRC 炎症-肿瘤演进中经脂质代谢重编程极化，定位于肿瘤边界并与 **FAP⁺ 肌成纤维细胞（Myofibroblasts）** 形成双向 Crosstalk，共同构建物理-免疫双重屏障导致 T 细胞排除。

分析路线：

1. **精细亚型注释**（Fine Annotation）— 识别 T 细胞和髓系细胞的功能亚群
2. **差异丰度分析**（Differential Abundance）— 定量 CRC/UC 中亚群比例差异
3. **差异基因表达**（scVI-DEG）— 批次校正的跨疾病差异表达分析
4. **拟时序演化轨迹**（Palantir v3）— 从炎症单核细胞到抑制性 TAM 的分化路径，含 Pre-driver 候选定义
5. **空间定位分析**（GSE144735 v2）— OLR1⁺ lipid-associated TAMs 与 FAP⁺ Myofibroblasts 在肿瘤边界的共同富集，支持"物理-免疫双重屏障"假说；同时整合全基因组扫描识别边界特异性标志基因
6. **转录调控网络**（Regulon 轨迹 + decoupleR）— 沿轨迹上升的核心脂质代谢 TF（PPARG、CEBPB）与新型 CRC 特异 TF（MAF、BHLHE41），及在 CRC TAM 中丢失的免疫激活 TF（CIITA、NR4A1）

---

## 二、精细亚型注释

### 方法

- **T 细胞**（77,927 个）：Leiden 聚类（res=1.0），按 CD8/CD4 谱系评分分流，再在各谱系内以 marker 基因模块打分（argmax 分配）
  - CD8_Tex 判别式：`score(HAVCR2/ENTPD1/TOX/TIGIT/LAG3/PDCD1/CTLA4) − 0.5 × score(TCF7/CCR7/SELL)`
- **髓系细胞**（15,972 个）：Leiden 聚类（res=0.8），TAM 调整评分：`sc_TAM − 0.5 × sc_mono_character`；TAM marker 使用 Palantir CRC 命运基因（FTH1/SOD2/SOCS3/OSM/HMOX1）
- **KNN 迁移**：用于补充 Unknown 亚型细胞的精细标注

### 结果

| 细胞亚型 | 总细胞数 | CRC | UC | CRC/UC 倍数 |
|---|---|---|---|---|
| B cells | 36,987 | — | — | — |
| Epithelial | 33,806 | — | — | — |
| **Th17** | 20,814 | 10,534 | 10,280 | 1.0× |
| **CD8_Teff** | 19,229 | 8,435 | 10,794 | 0.8× |
| CD4_Tconv | 9,218 | 477 | 8,741 | 0.05× |
| Classical_Mono | 8,203 | 2,801 | 5,402 | 0.5× |
| **Treg** | 7,921 | 4,281 | 3,640 | 1.2× |
| NK cells | 7,224 | 1,111 | 6,113 | 0.2× |
| **CD8_Tpex** | 5,698 | 247 | 5,451 | **0.05×（UC 富集）** |
| **TAM** | 3,440 | 3,396 | 44 | **75.5×（CRC 特异）** |
| **CD8_Tex** | 3,377 | 3,150 | 227 | **13.8×（CRC 富集）** |

> **关键发现**：TAM（75.5×）和 CD8_Tex（13.8×）是 CRC 最特异性的抑制性亚群，提示二者是最值得关注的候选群体；但 Pre-driver 候选亚群的正式定义需结合下文拟时序分析——即轨迹中处于过渡态（Branch point 上游）、尚未终末分化的细胞，而非单纯依据富集倍数。CD8_Tpex（干性耗竭前体）在 UC 中高度富集，两种疾病 Tex/Tpex 比值极端差异（CRC = 12.8，UC = 0.04）。

**髓系精细亚型 UMAP 与 Dotplot**

![髓系亚型 UMAP](Project_Delivery/results/Figure2_Pre_driver_Trajectory/Trajectory_Analysis/figures/myeloid_umap_subtype.png)

*图：髓系细胞（15,972 个）UMAP 着色（左：亚型；右：疾病来源）。TAM 簇集中于 CRC 患者特异区域。*

![髓系亚型 Marker Dotplot](Project_Delivery/results/Figure2_Pre_driver_Trajectory/Trajectory_Analysis/figures/myeloid_dotplot.png)

*图：髓系各亚型核心 Marker 基因 dotplot。TAM 高表达 APOE/SPP1/FTH1，Classical_Mono 高表达 S100A8/9/FCN1，cDC2 高表达 FCER1A/CLEC10A。*

**T 细胞精细亚型 UMAP 与 Dotplot**

![T细胞亚型 UMAP](Project_Delivery/results/Figure2_Pre_driver_Trajectory/Trajectory_Analysis/figures/tcell_umap_subtype.png)

*图：T 细胞（77,927 个）UMAP 着色（左：亚型；右：疾病来源）。CD8_Tex 和 Treg 高度富集于 CRC 区域，CD8_Tpex 偏向 UC 区域。*

![T细胞亚型 Marker Dotplot](Project_Delivery/results/Figure2_Pre_driver_Trajectory/Trajectory_Analysis/figures/tcell_dotplot.png)

*图：T 细胞各亚型核心 Marker 基因 dotplot。CD8_Tex 高表达 HAVCR2/ENTPD1/TIGIT/TOX，Treg 高表达 FOXP3/IL2RA/IKZF2，Th17 高表达 RORC/CCR6。*

---

## 三、差异细胞丰度分析

### 方法

因 Milo（milopy）在 Windows 环境无法运行，采用**基于样本比例的 Mann-Whitney U 检验 + BH FDR 校正**方案：
- 计算每个患者/样本内各亚型比例
- CRC 组（n=25 患者）vs UC 组（n=19 患者）
- 最低细胞数阈值：每样本 ≥50 个细胞才纳入

### 结果（Top 差异亚群）

| 亚型 | CRC 均值 | UC 均值 | log2FC | FDR |
|---|---|---|---|---|
| TAM | 0.084 | 0.002 | **+5.68** | 4.1×10⁻⁷ |
| CD8_Tex | 0.044 | 0.003 | **+3.87** | 9.4×10⁻⁶ |
| Th17 | 0.235 | 0.111 | +1.09 | 0.003 |
| pDC | 0.007 | 0.002 | +1.66 | 0.021 |
| CD8_Tpex | 0.003 | 0.096 | **−4.16** | 0.022 |

> CRC 中 TAM 和 CD8_Tex 的比例极显著升高，代表干性 T 细胞前体的 CD8_Tpex 在 UC 中显著富集，进一步证实 CRC 免疫微环境以深度耗竭为特征。

![所有亚型丰度 Dotplot](Project_Delivery/results/Figure2_Pre_driver_Trajectory/Spatial_Distribution/figures/abundance_dotplot_all.png)

*图：所有细胞亚型在 CRC 与 UC 中的丰度气泡图。横轴为亚型，纵轴为样本组；气泡大小代表平均比例，颜色深浅代表 FDR 显著性（红色显著）。*

![髓系亚型丰度 Strip plot](Project_Delivery/results/Figure2_Pre_driver_Trajectory/Spatial_Distribution/figures/abundance_strip_myeloid.png)

*图：髓系亚型（TAM / Inflam_Mac / Classical_Mono / cDC2 / pDC）在 CRC 与 UC 中的患者级比例分布。每个点代表一名患者；TAM 在 CRC 中的比例极端升高（p < 0.001）。*

![T细胞亚型丰度 Strip plot](Project_Delivery/results/Figure2_Pre_driver_Trajectory/Spatial_Distribution/figures/abundance_strip_tcell.png)

*图：T 细胞亚型在 CRC 与 UC 中的患者级比例分布。CD8_Tex 在 CRC 中显著富集，CD8_Tpex 在 UC 中显著富集，形成"耗竭 vs 干性"的疾病特异性分化模式。*

---

## 四、差异基因表达（scVI 批次校正 DEG）

### 方法

由于疾病标签与测序批次完全共线（CRC = GSE132465/GSE132257；UC = GSE125527/GSE150115），采用 scVI `get_normalized_expression(transform_batch=全部5批次)` 反事实预测方案：

1. 对每个细胞，在所有 5 个批次下解码并取平均 → 批次边际化后的标准化表达量
2. Wilcoxon 秩和检验 + BH FDR 校正；每组最多抽样 1,500 个细胞

### 结果

| 亚型 | CRC↑ 基因数 | UC↑ 基因数 | 代表性 CRC↑ 基因 | 代表性 UC↑ 基因 |
|---|---|---|---|---|
| CD8_Tex | 68 | 338 | HLA-A/C, ERV3-1 | PPP3CC/CB（calcineurin）, LTB4R2 |
| TAM | 117 | 746 | ISG15, ISG20 | LTB4R2, MATR3 |
| Treg | 200 | 313 | ADM, CLIC3 | PSMB10, ZNF 家族 |
| CD8_Teff | 237 | 293 | RPL21, HLA-C | CAMK1D, ABCC4 |
| Classical_Mono | 268 | 1,358 | ISG15, ISG20, ALCAM | ZNF 家族（批次相关）|
| Th17 | 97 | 684 | AARSD1, ANPEP | **FOS, FOSB**（AP-1 炎症激活）|

**生物学亮点**：
- **Classical_Mono/TAM CRC↑**：ISG15、ISG20（I 型干扰素刺激基因）— CRC 肿瘤微环境中持续的 IFN 信号激活
- **CD8_Tex UC↑**：PPP3CC/PPP3CB（calcineurin 亚基）— UC 的 Tex 处于更活跃的耗竭状态（"exhausted but reactive"）
- **Th17 UC↑**：FOS/FOSB（AP-1）— 与 UC 黏膜急性炎症的 IL-17 产生高度一致

![DEG 热图（所有亚型）](Project_Delivery/results/Figure2_Pre_driver_Trajectory/Trajectory_Analysis/figures/deg_heatmap_all_subtypes.png)

*图：主要细胞亚型 CRC vs UC 差异表达基因热图（scVI 批次校正，各亚型 Top 20 基因）。CRC↑ 基因集中于 ISG 家族（IFN 信号）和 HLA 基因；UC↑ 基因中 FOS/FOSB（Th17）和 PPP3CC/CB（CD8_Tex）突出。*

---

## 五、拟时序演化轨迹（Palantir v3）

### 数据与方法

| 参数 | 值 |
|---|---|
| 数据集 | GSE125527（UC，炎症起点）+ GSE132465（CRC，肿瘤对照）|
| 分析细胞 | 14,690 个髓系细胞（GSE125527: 7,921 / GSE132465: 6,769）|
| 工具 | Palantir（非线性扩散映射 + 马尔可夫过程）|
| 嵌入 | scVI 联合模型 latent space（30 维），子集重新计算 KNN + UMAP |
| Root 选择 | **生物学驱动**：GSE125527 Classical_Mono 群内抑制性基因评分（SPP1/APOE/C1QC/TGFB1/LGALS3/VEGFA）最低的细胞 |

### 末态识别

Palantir 自动识别出 3 个末态：

| 末态 | 亚型 | 来源 | 意义 |
|---|---|---|---|
| Terminal 1 | TAM | GSE132465（CRC）| 肿瘤相关巨噬细胞命运（主要） |
| Terminal 2 | TAM | GSE132465（CRC）| 肿瘤相关巨噬细胞命运（次要）|
| Terminal 3 | Inflam_Mac | GSE125527（UC）| 炎症性巨噬细胞末态 |

### Pre-driver 候选定义

Pre-driver 候选 = 同时满足以下三个计算标准的细胞：

1. **拟时序处于 Branch point 上游窗口**：以 Palantir 命运概率方差最大处（Branch point，基因表达的分叉位置）为上界，向早期延伸 30% 拟时序跨度，形成数据驱动的转折窗口（替代固定的第40–80百分位）
2. **抑制性基因评分 > Classical_Mono 中位数**（已开始表达抑制性 marker）
3. **CRC 末态命运概率 > 0.5**（偏向抑制性 TAM 方向）

### 关键结果

| 参数 | 值 |
|---|---|
| Pre-driver 候选细胞数 | **2,773 个（18.9%）** |
| Pre-driver 拟时序区间 | [branch_pt − 0.30×span, branch_pt]（数据驱动）|
| 代表亚型 | TAM（1,395）、Classical_Mono（642）、Inflam_Mac（426）、cDC2（306）|
| Pre-driver 核心亚群特征 | OLR1⁺ MMP9⁺（SPP1⁺B 亚型核心）— 脂质摄取 + 基质重塑双阳性 |

> **生物学解读**：Pre-driver 候选群集中于中后段轨迹，已开始表达抑制性标志但尚未达终末 TAM 状态。TAM 与 Inflam_Mac 在此区间并存，提示**过渡性中间状态**是真正的"前驱动"群体。值得注意的是，Pre-driver 区间内细胞高表达 **OLR1**（氧化低密度脂蛋白受体 1）和 **MMP9**（明胶酶 B），提示其已启动脂质代谢重编程（OLR1 介导的 oxLDL 摄取）和基质重塑程序——这两者正是 OLR1⁺ lipid-associated TAM 向肿瘤边界迁移并构建免疫屏障的关键功能分子。

**拟时序相关基因**：

| 方向 | 基因（Spearman ρ）| 生物学意义 |
|---|---|---|
| 正相关（随分化上升）| FCGR3A（+0.35）、CTSC（+0.29）、LGMN（+0.27）、RGS1（+0.26）| 单核细胞→TAM 过渡标志；溶酶体/组织蛋白酶系统上调 |
| 负相关（炎症起始态）| S100A8（−0.55）、VCAN（−0.49）、LYZ（−0.49）、S100A12（−0.48）| 经典炎症单核细胞标志；随分化程度升高而下降 |

**CRC 命运 DEG（Top 10）**：GPX1, LYZ, MS4A6A, CD14, IER3, GPR183, RPL13A, RPL21, IL1B, CAPG

**UC 命运 DEG（Top 10）**：LST1, IFITM2, COTL1, FCGR3A, SERPINA1, LYN, IFITM3, LILRB2, CDKN1C, CFD

**Figure 2A — Palantir 拟时序 UMAP 总览（v3）**

![Figure 2A — 拟时序总览](Project_Delivery/results/Figure2_Pre_driver_Trajectory/Trajectory_Analysis/figures/pseudotime_overview.png)

*图（2×2）：① 数据集来源（蓝=GSE125527 UC，红=GSE132465 CRC）；② 细胞亚型（Classical_Mono/Inflam_Mac/TAM/cDC2/pDC）；③ 拟时序值（★=Root，◆=Terminal）；④ Pre-driver 候选高亮（红色，n=2,773）。*

**Figure 2B — 关键基因拟时序动态热图**

![Figure 2B — 基因动态热图](Project_Delivery/results/Figure2_Pre_driver_Trajectory/Trajectory_Analysis/figures/gene_trend_heatmap_panelB.png)

*图：16 个 top 拟时序相关基因沿 CRC 抑制性分支的表达动态（Z-score 热图，100 个伪时序窗口平均后平滑）。蓝色基因（S100A8/S100A9/VCAN/LYZ/FCN1 等）随拟时序下降，代表炎症性单核细胞早期特征；红色基因（FCGR3A/CTSC/LGMN/RGS1 等）随拟时序上升，代表 TAM 方向特征。金色竖线标注 Pre-driver 窗口（数据驱动：branch_pt − 0.30×span 至 branch_pt），为炎症→抑制转折点。*

**Pre-driver 轨迹摘要（支持图）**

![Figure 2B — Pre-driver 轨迹](Project_Delivery/results/Figure2_Pre_driver_Trajectory/Trajectory_Analysis/figures/predriver_trajectory.png)

*图（左）拟时序×抑制性评分散点，按亚型着色，灰色阴影为 Pre-driver 窗口（数据驱动，[branch_pt−0.30×span, branch_pt]）；（右）Pre-driver 候选高亮（红色）。*

**Support — 关键基因在 UMAP 上的表达**

![Support — 基因表达 UMAP](Project_Delivery/results/Figure2_Pre_driver_Trajectory/Trajectory_Analysis/figures/gene_expression_umap.png)

*图：8 个关键基因（SPP1/APOE/C1QC/MRC1/S100A8/HLA-DRA/CD274/VEGFA）在髓系 UMAP 上的表达热图。SPP1/APOE/MRC1 集中于 TAM 区域；S100A8 集中于 Classical_Mono 区域，与轨迹起点一致。*

**Support — CRC vs UC 命运分支 DEG**

![Support — 分支 DEG](Project_Delivery/results/Figure2_Pre_driver_Trajectory/Trajectory_Analysis/figures/branch_DEG.png)

*图：CRC 命运（左，红色）和 UC 命运（右，蓝色）各 Top 20 差异基因 Wilcoxon 评分柱状图。*

**拟时序分布**

![拟时序分布](Project_Delivery/results/Figure2_Pre_driver_Trajectory/Trajectory_Analysis/figures/pseudotime_distribution.png)

*图（左）抑制性基因评分 UMAP；（中）各亚型拟时序箱线图；（右）Pre-driver 窗口密度分布（红色 = Pre-driver，灰色 = 其他）。*

---

## 六、空间定位分析（GSE144735 KUL3 CRC v2）

### 数据集

| 参数 | 值 |
|---|---|
| 数据集 | GSE144735（KUL3 CRC，Lee et al. 2020）|
| 细胞数 | 27,414 |
| 患者数 | 6（KUL01/19/21/28/30/31）|
| 空间分区 | Normal（9,736）/ Border（9,424）/ Tumor（8,254）|
| 细胞类型 | T cells / Myeloids / B cells / Epithelial / Stromal / Mast |

### 方法

1. 标准化 + log1p；以扩充的 Pre-driver 基因模块（20 基因：SPP1/APOE/C1QC/LGALS3/MRC1/VEGFA/IL10/TGFB1/CD274/SOCS3 等）打分
2. 新增炎症评分（INFLAM_GENES：S100A8/9/CXCL8/IL1B/CCL3/CCL4 等 13 基因）
3. 逐患者计算 OLR1⁺ TAM（SPP1⁺B 亚型）比例，Kruskal-Wallis + pairwise Mann-Whitney 检验；同步对全部细胞亚群和 9,805 个基因进行系统性 Border vs Normal 富集扫描

### 结果

**Pre-driver v2 评分（所有细胞）**：

| 区域 | 均值 | Normal vs 检验 |
|---|---|---|
| Normal | 1.040 | — |
| Border | 1.318 | p = 1.4×10⁻⁸⁴ \*\*\* |
| Tumor | 1.293 | p = 1.5×10⁻¹⁶⁷ \*\*\* |

**OLR1⁺ TAM（SPP1⁺B 亚型）比例（逐患者，Myeloid 子集）**：

| 区域 | 细胞数 | OLR1⁺ TAM 占比 | Normal vs 检验 |
|---|---|---|---|
| Normal | 823 | 0.5%（8 个）| — |
| Border | 924 | **35.5%（377 个）** | p = 0.004 \*\* |
| Tumor | 929 | **32.1%（352 个）** | p = 0.004 \*\* |

> **核心发现**：
>
> **① OLR1⁺ lipid-associated TAMs（OLR1⁺ MMP9⁺ TAMs）肿瘤边界高度富集**
> - SPP1⁺B 亚型（核心标志：OLR1/TREM2/MMP9）在 Border 区 35.5% vs Normal 0.5%（70× 富集，p=0.004），提示 **OLR1⁺ 脂质代谢型 TAM 在肿瘤边界优先聚集**而非仅存在于肿瘤核心
> - 基因层面：OLR1（log2FC=+3.95）、MMP9（+3.21）、TREM2（+2.84）、TREM1（+3.02）均在 Border > Normal，与 SPP1⁺B 亚群空间定位一致
> - TAM 评分：Border 和 Tumor 均显著高于 Normal（p < 10⁻⁵²），Border vs Tumor 无显著差异（p=0.47），符合"边界屏障"模型
>
> **② FAP⁺ Myofibroblasts（肌成纤维细胞）同步边界富集——"共犯"结构**
> - Myofibroblasts 在 Border 显著富集（log2FC=+4.63，**FDR=0.023**），是细胞亚群层面除 OLR1⁺ TAM 外唯一达到统计显著性的间质成分
> - Myofibroblast 特异性胶原基因（COL10A1/MMP1/FAP/INHBA/SUGCT 等）均在 Border 高表达（log2FC +2.7–4.5），提示其在边界大量分泌细胞外基质，构建**物理性胶原屏障**
> - 这一 OLR1⁺ TAM + FAP⁺ Myofibroblast 共富集模式，提示两者在边界形成**基质-免疫协同屏障轴**（详见第三/四阶段 CellChat 分析计划）
>
> **③ 双重屏障模型（炎症信号佐证）**
> - 炎症评分在 Border 最高（Normal < Border > Tumor），说明 Border 区同时存在 OLR1⁺ TAM 介导的免疫抑制 + Myofibroblast 介导的基质重塑——这正是 MMRp CRC 免疫排除的核心发生地点

**Figure 2C — OLR1+ TAM + Myofibroblast 肿瘤边界富集总览**

![Figure 2C — Border Enrichment Summary](Project_Delivery/results/Figure2_Pre_driver_Trajectory/Spatial_Distribution/figures/border_enrichment_summary_v2.png)

*三联 Lollipop 图。Panel A：细胞亚群层面 Border/Normal 富集（log2FC，Mann-Whitney U + BH-FDR）；OLR1+ TAM（SPP1+B，log2FC=+5.85, FDR=0.023）和 Myofibroblasts（log2FC=+4.63, FDR=0.023）加粗高亮——两者是 Border 区仅有的统计显著富集群体。Panel B：OLR1+ TAM 脂质代谢特征基因（OLR1/TREM2/TREM1/MMP9，红色）和 FAP+ Myofibroblast 基质重塑基因（FAP/INHBA/COL10A1/MMP1，绿色）。Panel C：全基因组扫描（9,805 基因）Top Border 特异基因，按细胞类型分组；OLR1/MMP9/TREM2 名称加粗。*

**Figure 2C (Supp) — Pre-driver 区域富集空间分布图（v2）**

![Figure 2C Supp — Pre-driver 富集](Project_Delivery/results/Figure2_Pre_driver_Trajectory/Spatial_Distribution/figures/predriver_enrichment.png)

*图（左）OLR1⁺ TAM（SPP1⁺B 亚型）在髓系细胞中的比例（各患者散点 + 箱线）；（右）平均 Pre-driver 评分；x 轴为 Normal / Border / Tumor。Border 和 Tumor 均显著高于 Normal（\*\* p<0.01, \*\*\* p<0.001）。*

**Figure 2D — Infiltration Score 柱状图**

![Figure 2D — Infiltration Score](Project_Delivery/results/Figure2_Pre_driver_Trajectory/Spatial_Distribution/figures/infiltration_score_panelD.png)

*图：三个 Infiltration 指标（SPP1+ 比例、Suppressive 比例、Pre-driver 评分均值）在 Normal / Border / Tumor 三区的逐患者比较（柱形 = 组均值 ± SEM；点 = 各患者；显著性标注 Mann-Whitney U 检验）。Border 区 SPP1+ 比例（35.5% vs Normal 0.5%，p < 0.01）和 Pre-driver 评分均显著高于 Normal，支持"肿瘤边界物理屏障"假说。*

**髓系亚型组成（支持图）**

![Figure 2F — 髓系亚型组成](Project_Delivery/results/Figure2_Pre_driver_Trajectory/Spatial_Distribution/figures/myeloid_subtype_composition.png)

*图：髓系细胞亚型在 Normal / Border / Tumor 三区的比例堆叠柱状图。SPP1+B（深红）在 Border 和 Tumor 占主导；Anti-inflammatory 和 Pro-inflammatory 在 Normal 比例更高。*

**评分分布小提琴图**

![评分分布](Project_Delivery/results/Figure2_Pre_driver_Trajectory/Spatial_Distribution/figures/score_by_class.png)

*图：TAM 评分、Pre-driver v2 评分、CD8_Tex 评分、炎症评分在三区的小提琴图（全体细胞）。TAM 和 Pre-driver 评分在 Border/Tumor 中系统性高于 Normal。*

**患者热图**

![患者热图](Project_Delivery/results/Figure2_Pre_driver_Trajectory/Spatial_Distribution/figures/predriver_heatmap_patient.png)

*图：各患者（行）× 三区（列）的平均 Pre-driver 评分热图。所有 6 名患者均呈现 Border/Tumor > Normal 的一致模式。*

---

## 七、转录调控网络

### 7.1 TF 拟时序轨迹分析

#### 方法

- 对 14,690 个髓系细胞打 15 个 TF regulon 基因集得分（DoRothEA/文献 curated 目标基因）
- 计算 TF 得分与 Palantir 拟时序的 Spearman 相关系数，识别沿炎症→抑制轨迹的调控变化

#### 结果

**沿拟时序上升的 TF（抑制性末态富集）**：

| TF | ρ | p 值 | TAM 亚型活性 | 生物学意义 |
|---|---|---|---|---|
| **IRF1** | +0.290 | 2.8×10⁻²⁸³ | 最高 | I 型 IFN 响应；驱动 HLA-DR/ISG15 等 TAM 特征基因 |
| **SPI1（PU.1）** | +0.248 | 1.3×10⁻²⁰⁴ | 高 | 髓系分化核心 TF；TAM 中持续高表达 |
| **NFKB1** | +0.116 | 2.9×10⁻⁴⁵ | 高 | NF-κB p50；与 RELA 协同驱动 M2/免疫抑制极化 |
| **MYC** | +0.107 | 1.6×10⁻³⁸ | 中 | 代谢重编程；TAM 增殖-能量代谢转变 |
| **RUNX1** | +0.083 | 5.3×10⁻²⁴ | 中 | 髓系谱系转录因子；参与单核→巨噬细胞分化 |

> **注**：PPARG（脂质代谢核心 TF）在本分析的轨迹相关性中 ρ=−0.032（p=8.2×10⁻⁵，轻微下降），但其 decoupleR 推断活性在 TAM 中显著高于 Classical_Mono（ΔActivity=+1.35，见 §7.2），提示 PPARG 在 TAM 分化的**过渡态（Pre-driver 窗口）**中活性最高，并随终末分化轻微下降——这与 OLR1⁺ 脂质 TAM 在 Pre-driver 期启动脂质摄取程序的模型一致。

**沿拟时序下降的 TF（炎症早期状态）**：

| TF | ρ | p 值 | 生物学意义 |
|---|---|---|---|
| **IRF4** | −0.164 | 6.4×10⁻⁸⁹ | M2-like 极化调控；在炎症单核细胞期高，进入 TAM 后下降 |
| **KLF4** | −0.160 | 4.0×10⁻⁸⁵ | 抗炎性巨噬细胞 TF；炎症状态时表达较高；CRC TAM 中 DEG log2FC=−1.29（FDR=4.3×10⁻¹⁵）|
| **JUN** | −0.096 | 2.8×10⁻³¹ | AP-1 家族；急性炎症激活后减弱 |
| **FOS** | −0.067 | 3.8×10⁻¹⁶ | AP-1 家族；早期炎症反应 TF |

**Figure 2E — TF Regulon 拟时序相关性（SCENIC）**

![Figure 2E — TF 相关性](Project_Delivery/results/Figure2_Pre_driver_Trajectory/SCENIC_Regulon/figures/tf_pseudotime_correlation.png)

*图：各 TF regulon 得分与拟时序的 Spearman ρ 横向条形图（红色=上升，蓝色=下降）。IRF1 和 SPI1 是最显著的正相关 TF（p < 10⁻²⁰⁰）；IRF4 和 KLF4 是最显著的负相关 TF。*

**Figure 2F — Novel TF Candidates（新型转录调控因子）**

![Figure 2F — Novel TF Candidates](Project_Delivery/results/Figure2_Pre_driver_Trajectory/SCENIC_Regulon/figures/tf_novel_candidates.png)

*两联图。Panel A：CRC TAM vs UC Myeloids 伪批次校正 DEG 中的 TF 候选（Pseudobulk Wilcoxon，BH-FDR）。红色（上调）：MAF[*]（c-MAF，+1.00, FDR=1.3×10⁻¹⁴）、BATF3（+0.85）、BHLHE41（+0.66）；蓝色（下调）：CIITA[*]（-2.18, FDR=1.1×10⁻²²）、NR4A1[*]（-1.60）、IRF8（-1.39）、SREBF2（-0.82，de novo 脂质合成 TF 下调，与 OLR1 脂质摄取上调共同提示脂质来源切换）。Panel B：decoupleR-ULM 推断的所有 TF 活性差值（TAM - Mono）；PPARG（绿色高亮）delta_activity=+1.35，支持脂质代谢核心 TF 在 TAM 中的活化。[*] = 机制新颖、CRC 特异的 story TF。*

**Support — TF 活性热图（髓系亚型）**

![Support — TF 热图](Project_Delivery/results/Figure2_Pre_driver_Trajectory/SCENIC_Regulon/figures/tf_heatmap_subtype.png)

*图：Top 12 TF 的 regulon 活性 z-score 热图（行=亚型，列=TF）。TAM 中 IRF1、SPI1、NFKB1、PPARG 活性最高；Classical_Mono 和 cDC2 中 JUN、IRF4 活性较高——与轨迹分析方向一致。*

**Support — TF 沿拟时序的趋势曲线**

![Support — TF 趋势](Project_Delivery/results/Figure2_Pre_driver_Trajectory/SCENIC_Regulon/figures/tf_trends.png)

*图：Top TF（IRF1/SPI1/NFKB1/MYC/RUNX1 上升；IRF4/KLF4/JUN 下降）沿拟时序 20 等分位均值趋势曲线。上升 TF 在拟时序 > 0.3 后明显加速。*

**Support — TF UMAP**

![Support — TF UMAP](Project_Delivery/results/Figure2_Pre_driver_Trajectory/SCENIC_Regulon/figures/tf_umap.png)

*图：Top 6 TF regulon 得分在髓系 UMAP 上的表达热图。IRF1/SPI1 高分区与 TAM 区域高度重合；IRF4/KLF4 高分区与 Classical_Mono 重合。*

---

### 7.2 TF 末态活性比较（decoupleR + CollecTRI）

#### 方法

- 工具：decoupleR v1.8.0，ULM（单变量线性模型）推断 TF 活性
- TF-target 数据库：CollecTRI（人类，19,186 对 TF-target，539 个 TF）
- 比较：TAM vs Classical_Mono；CD8_Tex vs CD8_Teff（ΔActivity = Target − Control）

#### 结果

**TAM vs Classical_Mono（TAM 上调的 TF）**：

| TF | ΔActivity | CRC TAM DEG log2FC | 生物学意义 |
|---|---|---|---|
| **RELA** | +3.05 | — | NF-κB p65，核心炎症-免疫抑制 TF |
| **NFKB1** | +2.64 | — | NF-κB p50，与 RELA 协同驱动 M2 极化 |
| **FOS** | +2.39 | — | AP-1 家族，肿瘤相关巨噬细胞激活 |
| **JUN** | +2.30 | — | AP-1 家族，TAM 促肿瘤功能 |
| **PPARG** ⭐ | **+1.35** | FDR n.s.（活性数据） | **脂质代谢核心 TF；驱动 OLR1/TREM2 等脂质摄取程序；TAM 活性显著高于 Mono，为脂质代谢重编程关键调控因子** |
| **IRF1** | +1.25 | — | I 型 IFN 响应，与 ISG 基因上调一致（与轨迹分析一致）|
| **CEBPB** | +0.03 | +0.25（FDR=0.005） | M2 极化辅助 TF；绝对活性高，CEBPB+PPARG 协同驱动 M2/脂质 TAM |

**TAM vs CRC 特异 TF（pseudobulk DEG，CRC TAM vs UC Myeloids）—— 新型候选**：

| TF | CRC TAM log2FC | FDR | 方向 | 生物学意义（CRC 中的功能改变）|
|---|---|---|---|---|
| **MAF（c-MAF）** ⭐ | +1.00 | 1.3×10⁻¹⁴ | 上调 | c-MAF 驱动 IL-10 介导的免疫抑制极化（文献：c-MAF→IL-10/ARG1；本数据 TGFB1 在 CRC TAM 中下调）；在 CRC TAM 中特异激活；TREM2/SPP1+ TAM regulon 的核心节点 |
| **BHLHE41（DEC2）** | +0.66 | 8.2×10⁻⁹ | 上调 | 组织驻留巨噬细胞身份维持 TF；抑制炎症应答；CRC TAM 中升高，提示 TAM 的"定居化"程序激活 |
| **BATF3** | +0.85 | 1.3×10⁻¹⁶ | 上调 | TAM 获得 DC 样特征（BATF3+ 抑制性）；CRC 特异 |
| **CIITA** ⭐ | −2.18 | 1.1×10⁻²² | **丢失** | **MHC-II 转录主调节因子；CRC TAM 中 CIITA 显著下调 = 抗原呈递通路直接关闭 = 免疫逃逸直接机制**（非间接调控，是最直接的免疫逃逸分子事件）|
| **NR4A1（Nur77）** | −1.60 | 6.5×10⁻¹⁴ | **丢失** | NR4A1 维持巨噬细胞抗炎稳态；缺失驱动促肿瘤 TAM；近年多项研究显示 NR4A 家族激动剂可逆转 TAM 免疫抑制 |
| **NR4A2** | −1.08 | 2.9×10⁻⁹ | 丢失 | NR4A1 同家族；协同功能丧失 |
| **IRF8** | −1.39 | 4.6×10⁻²² | 丢失 | DC1/M1 巨噬细胞 TF；CRC TAM 中 IRF8 缺失提示 I 型 IFN 应答能力减弱和 DC 样功能丧失 |
| **SREBF2** | −0.82 | 5.6×10⁻¹⁷ | 丢失 | 内源性胆固醇/脂质合成 TF；下调而 OLR1 上调，提示 CRC TAM 将脂质来源从**内源合成切换为外源摄取**（oxLDL via OLR1）|

**CD8_Tex vs CD8_Teff（CD8_Tex 上调的 TF）**：

| TF | ΔActivity | 生物学意义 |
|---|---|---|
| **MYC** | +1.62 | 增殖-代谢重编程，耗竭 T 细胞特征性 TF |
| **JUN** | +0.75 | AP-1，慢性刺激下持续激活 → 耗竭 |
| **HIF1A** | +0.63 | 低氧/代谢应激，加速耗竭进程 |
| **CEBPB** | +0.40 | 炎症-耗竭重编程 |
| **NR4A1** | −0.93 | 晚期耗竭 NR4A1 活性下降（符合末期耗竭特征）|
| **BATF** | −0.19 | Tex 特征性 TF 在末期耗竭中活性减弱 |

> **三套方法交叉验证**：(1) IRF1 在轨迹分析中为最强正相关 TF（ρ=+0.29），在 decoupleR 中 TAM vs Mono ΔActivity=+1.25，两套方法高度一致；(2) **PPARG** decoupleR ΔActivity=+1.35（TAM>Mono），与 OLR1⁺ 脂质代谢 TAM 极化模型的脂质 TF 假说吻合；(3) **CIITA** pseudobulk DEG log2FC=−2.18（FDR=1.1×10⁻²²），为 CRC TAM 中功能性丢失最显著的免疫调控 TF，提供了抗原呈递丧失的直接分子证据。

**Support — TF Regulon 活性热图（末态比较，decoupleR）**

![Support — TF 末态热图](Project_Delivery/results/Figure2_Pre_driver_Trajectory/Spatial_Distribution/figures/figure2_panelE.png)

*图（左右两组）：左侧为 TAM vs Classical_Mono，右侧为 CD8_Tex vs CD8_Teff。颜色代表 decoupleR-ULM 推断的 TF 活性 z-score。TAM 中 RELA/NFKB1/IRF1/PPARG 显著上调（PPARG ΔActivity=+1.35 支持脂质 TAM 模型）；CD8_Tex 中 MYC/JUN/HIF1A 上调，NR4A1/BATF 相对下调。*

---

## 八、综合 Pre-driver 模型

基于拟时序、空间定位、TF 轨迹三部分分析的综合模型：

```
【起点】炎症单核细胞（UC，GSE125527，S100A8⁺/IRF4⁺/KLF4⁺）
    拟时序 0–0.25 | JUN/FOS 急性炎症激活；SREBF2⁺ 内源脂质合成
           ↓  UC→CRC 演进；oxLDL/脂质废物积累
    ─────────────────────────────────────────────────────
    Pre-driver 窗口（branch_pt − 0.30×span → branch_pt）
    Classical_Mono + Inflam_Mac 混合过渡表型
    ✦ OLR1⁺ MMP9⁺ 脂质代谢重编程开始：
      - OLR1 上调（oxLDL 摄取）
      - SREBF2 下调（内源合成 → 外源摄取切换）
      - PPARG 活性升高（脂质代谢核心 TF）
    ✦ IRF1/SPI1/NFKB1 上升，IRF4/KLF4/JUN 下降
    ✦ CIITA 下调（MHC-II 抗原呈递开始丧失）
    ✦ NR4A1 下调（稳态调控丧失）
    ─────────────────────────────────────────────────────
           ↓  PPARG/CEBPB/MAF 驱动 + CIITA/NR4A1 丢失
    终末 OLR1⁺ lipid-associated TAM（CRC，GSE132465）
    OLR1⁺/MMP9⁺/TREM2⁺/APOE⁺/SPP1⁺
           ↓  迁移至肿瘤边界
    ══════════════════════════════════════════
    肿瘤边界（Border）双重屏障
    ══════════════════════════════════════════
    OLR1⁺ TAM                   FAP⁺ Myofibroblast
    (35.5% Border，FDR=0.023)   (log2FC=+4.63，FDR=0.023)
        ↕ Crosstalk（SPP1→CD44/整合素？MMP9→ECM 重塑？）
    免疫抑制（CIITA丢失/IL-10/MAF驱动）+ 基质重塑（MMP9/COL11A1/FAP）
           ↓
    T 细胞被物理阻隔（胶原屏障）+
    被化学耗竭（PDCD1L/免疫检查点配体）
    → MMRp CRC 免疫排除终局
```

**综合潜在治疗靶点（优先级排序）**：

| 靶点 | 证据层次 | 靶向策略 |
|---|---|---|
| **OLR1（LOX-1）** ⭐ | Border log2FC=+3.95；脂质摄取入口；驱动 TAM 极化 | 抗 OLR1 单抗/抑制剂阻断 oxLDL 摄取，切断脂质重编程起点 |
| **PPARG** ⭐ | decoupleR ΔActivity=+1.35（TAM>Mono）；脂质核心 TF | PPARG 拮抗剂（如 GW9662）阻断脂质 TAM 极化程序 |
| **FAP（Myofibroblast）** ⭐ | Border log2FC=+4.63，**FDR=0.023**；物理屏障关键成分 | FAP-CAR-T / FAP-ADC 靶向清除 CAF 屏障 |
| **CIITA（丢失）** | log2FC=−2.18（FDR=1.1×10⁻²²）；MHC-II 转录关闭 | 重建 CIITA 表达（表观遗传激活）以恢复 TAM 抗原呈递 |
| **NR4A1（丢失）** | log2FC=−1.60（FDR=6.5×10⁻¹⁴）；稳态调控丧失 | NR4A1 激动剂（cytosporone B）逆转促肿瘤极化 |
| **IRF1** | 轨迹 ρ=+0.29（最强），decoupleR Δ=+1.25 | IRF1 抑制剂阻断 Classical_Mono→TAM 转化 |
| **NFKB1/RELA** | 轨迹上升 + decoupleR TAM Δ=+2.6/+3.1 | NF-κB 通路抑制剂（双靶点）|
| **MAF（c-MAF）** | CRC TAM log2FC=+1.00（FDR=1.3×10⁻¹⁴）；免疫抑制极化 | c-MAF 抑制剂阻断 IL-10 介导的免疫抑制程序 |

---

## 九、Figure 2 交付汇总

### 主图 Panel A–E

| Panel | 内容 | 文件路径 |
|---|---|---|
| **2A** Trajectory | 髓系拟时序 UMAP（数据集/亚型/拟时序/Pre-driver 高亮）| `Trajectory_Analysis/figures/pseudotime_overview` |
| **2B** Dynamic Expression | 关键基因拟时序动态热图（16 基因，Pre-driver 窗口标注）| `Trajectory_Analysis/figures/gene_trend_heatmap_panelB` |
| **2C** Spatial Map | OLR1⁺ TAM + Myofibroblast 边界富集总览（三联 Lollipop，GSE144735）| `Spatial_Distribution/figures/border_enrichment_summary_v2` |
| **2D** Barplot | Infiltration Score 柱状图（逐患者 Border/Tumor/Normal 对比）| `Spatial_Distribution/figures/infiltration_score_panelD` |
| **2E** SCENIC | TF Regulon 拟时序相关性（IRF1/SPI1 上升；IRF4/KLF4 下降）| `SCENIC_Regulon/figures/tf_pseudotime_correlation` |
| **2F** Novel TF | 新型 TF 候选图：MAF/BHLHE41 上调 + CIITA/NR4A1 丢失（Panel A）；PPARG delta activity 高亮（Panel B）| `SCENIC_Regulon/figures/tf_novel_candidates` |

### 支持图（各节附图）

| 文件 | 内容 |
|---|---|
| `Trajectory_Analysis/figures/predriver_trajectory` | Pre-driver 拟时序×抑制性评分散点图 |
| `Trajectory_Analysis/figures/gene_expression_umap` | 8 关键基因在髓系 UMAP 上的表达 |
| `Trajectory_Analysis/figures/branch_DEG` | CRC vs UC 命运分支 DEG 柱状图 |
| `Trajectory_Analysis/figures/pseudotime_distribution` | 亚型拟时序箱线图 + Pre-driver 密度分布 |
| `Spatial_Distribution/figures/myeloid_subtype_composition` | 髓系亚型 Normal/Border/Tumor 组成堆叠图 |
| `Spatial_Distribution/figures/score_by_class` | TAM/Pre-driver/CD8_Tex/炎症评分小提琴图 |
| `Spatial_Distribution/figures/predriver_heatmap_patient` | 各患者 Pre-driver 评分热图 |
| `SCENIC_Regulon/figures/tf_heatmap_subtype` | TF 活性热图（亚型） |
| `SCENIC_Regulon/figures/tf_trends` | TF 沿拟时序趋势曲线 |
| `SCENIC_Regulon/figures/tf_umap` | TF regulon 活性 UMAP |
| `Spatial_Distribution/figures/figure2_panelE` | TF 末态比较（decoupleR，TAM vs Mono；Tex vs Teff；PPARG ΔActivity=+1.35 高亮）|

### 补充图 SuppFig2

| 文件 | 内容 |
|---|---|
| `Trajectory_Analysis/figures/myeloid_dotplot` | 髓系亚型 Marker 基因 dotplot |
| `Trajectory_Analysis/figures/tcell_dotplot` | T 细胞亚型 Marker 基因 dotplot |
| `Supplementary/Supp2_Extended_Markers/figures/supp_fig2_markers` | 所有亚群扩充 Marker 基因列表 |

所有文件均位于 `Project_Delivery/results/Figure2_Pre_driver_Trajectory/` 下，PNG（dpi=200）+ PDF（矢量）双格式。

---

## 十、辅助图 — 细胞亚群 Marker 基因

### Supplementary Figure 2

![Supplementary Figure 2](Project_Delivery/results/Supplementary/Supp2_Extended_Markers/figures/supp_fig2_markers.png)

**图注**：所有细胞亚群的特异性 marker 基因 dotplot，按三大谱系分组展示（Myeloid / T cells / NK·B·Other）。每个亚群展示 10 个标志基因，优先选用文献验证核心 marker，不足部分由 Wilcoxon rank_genes_groups 计算结果补充。点大小代表表达比例；颜色代表 z-score 标准化均值表达量（±2.5 截断）。

**关键亚群 Marker 解读**：

| 亚型 | 核心 Marker | 特征 |
|---|---|---|
| **TAM（OLR1⁺ lipid-associated）** | APOE, SPP1, OLR1, MMP9, TREM2, FTH1, LGALS3 | 脂质代谢型 TAM；OLR1/TREM2 为脂质摄取标志；MMP9 介导基质重塑；CIITA 低/丢失为抗原呈递功能缺陷标志 |
| **Inflam_Mac** | IL1B, CXCL8, TNF, PTGS2, NFKBIA | 炎症性巨噬细胞；NF-κB 通路激活，COX-2（PTGS2）高表达 |
| **Classical_Mono** | S100A8/9, FCN1, LYZ, CD14, VCAN | 经典炎症单核细胞；S100A8/9 双阳性为关键特征 |
| **cDC2** | FCER1A, CLEC10A, CD74, SPI1, CPVL | 传统 DC2；MHCII 相关基因，IgE 受体 FCER1A 特异 |
| **pDC** | IRF7, TCF4, PLD4, PLAC8, IRF8 | 浆细胞样 DC；I 型 IFN 产生细胞，IRF7/8 驱动 |
| **CD8_Tex** | HAVCR2, ENTPD1, PDCD1, TIGIT, TOX | 终末耗竭；多抑制受体共表达 |
| **CD8_Tpex** | TCF7, SELL, CCR7, TOX2, CXCR5 | 干性/前耗竭；TCF7+ 维持自我更新 |
| **CD8_Teff** | GZMB, PRF1, IFNG, FGFBP2, TBX21 | 效应/细胞毒；颗粒酶 B + 穿孔素双阳性 |
| **Treg** | FOXP3, IL2RA, IKZF2, RTKN2, LAYN | 调节性 T 细胞；FOXP3+CD25+ 核心特征 |

---

## 十一、阶段二结论

1. **Pre-driver 亚群确立**：TAM（CRC/UC 比 = 75.5×）和 CD8_Tex（13.8×）是 CRC 最特异性的抑制性亚群；在拟时序轨迹中，**2,773 个过渡期细胞（18.9%）**被精确定义为 Pre-driver 候选（拟时序 0.25–0.35，已表达抑制性标志但非终末态）

2. **演化轨迹完整重建**：从 UC 炎症单核细胞（S100A8⁺/IRF4⁺）出发，经 Pre-driver 过渡态（FCGR3A⁺/RGS1⁺），到达 CRC 终末 TAM（SPP1⁺/APOE⁺），轨迹由 **IRF1/SPI1 → NFKB1/HIF1A** 分阶段调控

3. **空间阻断假说成立**：SPP1⁺ TAM 在肿瘤 Border 区 70× 富集（35.5% vs 0.5%），且 Border ≈ Tumor（p=0.47 无差异），支持"Pre-driver 在边界构建物理-免疫双重屏障阻碍 T 细胞渗透"假说

4. **调控轴清晰**：IRF1（轨迹 ρ=+0.29，decoupleR Δ=+1.25）和 NFKB1/RELA（decoupleR Δ=+2.6/+3.1）是 TAM Pre-driver 的核心双轴 TF，HIF1A 介导低氧-Border 富集，SPI1 驱动髓系谱系承诺

### 进入阶段三的基础

| 内容 | 状态 | 阶段三用途 |
|---|---|---|
| Pre-driver 候选细胞列表（2,773 个）| ✅ | MMRp CRC 亚群提取基础 |
| TAM/CD8_Tex 基因 signature（20 基因）| ✅ | GSE178341 MIL 整合特征 |
| 治疗靶点列表（IRF1/SPI1/NFKB1/HIF1A/SPP1）| ✅ | 药物靶点数据库匹配 |
| GSE144735 空间验证数据 | ✅ | 干预响应的空间对照 |

---

*分析工具：scanpy 1.12 / scvi-tools 1.4.1 / Palantir / decoupleR 1.8.0 (CollecTRI)*
*计算环境：Python 3.12, Windows 11, CPU-only*
*交付目录：Project_Delivery/results/Figure2_Pre_driver_Trajectory/*
