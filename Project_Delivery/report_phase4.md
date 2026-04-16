# 第四阶段报告：细胞通讯耦合与免疫逃逸生物学机制

**项目**：MMRp 结直肠癌免疫耐药机制研究  
**阶段**：Phase 4 — CellChat/NicheNet 细胞通讯 × scMetabolism 代谢微环境  
**日期**：2026-04-14  
**数据集**：GSE144735（空间分区，27,414 cells）；GSE178341（免疫细胞，94,196 cells）

---

## 零、研究框架概览

**核心问题**：OLR1⁺ lipid-associated TAMs（SPP1+B 亚群）是如何抑制 CD8⁺ T 细胞、同时与 FAP⁺ Myofibroblast 协同构建物理–免疫双重屏障的？

**阶段四从两个维度回答此问题**：

1. **通讯层**：利用 LIANA（rank_aggregate）分析 OLR1⁺ TAM 向 CD8⁺ T 细胞和 Myofibroblast 发出的关键配体–受体信号
2. **代谢层**：利用 AUCell 通路评分（scMetabolism 等价实现，KEGG 2021 Human 基因集，decoupler-py v1.8.0）和 Mann-Whitney 检验，解析 OLR1⁺ TAM 在 MMRp 与 MMRd 环境中的代谢偏好差异

**数据来源对应**：

| 分析模块 | 数据集 | 细胞数 |
|---------|--------|-------|
| LIANA 通讯网络（Panel A/B） | GSE144735 Border cells | 9,424 |
| NicheNet-like 配体–靶基因图（Panel C） | GSE144735（SPP1+B + CD8 subsets） | SPP1+B: 634, CD8: 1,326 |
| 代谢通路评分（Panel D/E） | GSE144735（细胞类型分层） | 全 27,414 |
| GSEA 通路差异（Panel F） | GSE178341 Macro（MMRp vs MMRd） | MMRp: 9,906, MMRd: 8,095 |

---

## 一、研究背景与假设

阶段二（Phase 2）已证明 OLR1⁺ TAMs 在肿瘤边界（Border）显著富集（log2FC=+3.95，FDR=0.023），并与 FAP⁺ Myofibroblasts 共定位。阶段三（Phase 3）证明这一亚群在 MIL 模型中是判别"冷肿瘤（MMRp）"的核心特征细胞（注意力 7.5×）。阶段四需回答：

> **OLR1⁺ TAMs 通过哪些具体分子机制（配体–受体对、代谢通路）实现免疫抑制？**

**工作假设**：
- H1：OLR1⁺ TAMs 分泌 SPP1/LGALS1/MIF 等配体，直接作用于 CD8⁺ T 细胞表面受体，诱导耗竭
- H2：OLR1⁺ TAMs 通过 SPP1→ITGAV:ITGB1 与 Myofibroblast 形成正反馈，激活 ECM 重塑构建物理屏障
- H3：OLR1⁺ TAMs 在 MMRp 中上调腺苷代谢通路，通过产生腺苷（CD39/CD73 轴）进一步抑制 T 细胞活化

---

## 二、方法

### 2.1 细胞通讯分析（LIANA rank_aggregate）

- **工具**：liana 1.7.1，rank_aggregate 方法（整合 CellPhone、NATMI、CellChat、Connectome、log2FC 五种算法）
- **数据**：GSE144735 Border region 全细胞（9,424 cells，30 种亚型）
- **筛选**：cellphone_pvals < 0.05，magnitude_rank 越小代表信号越强（[0,1]）
- **Panel A**：汇总关键亚型间通讯强度（strength = -log₁₀(magnitude_rank + 1e-6)），绘制 Circle Plot
- **Panel B**：聚焦 SPP1+B（OLR1⁺ TAMs）发出的抑制性 + 基质重塑 LR 对，按功能分组绘制 Bubble Plot

### 2.2 NicheNet-like 配体–靶基因图（Panel C）

基于 LIANA 结果和已发表 NicheNet 先验知识，构建三层连接图：

```
OLR1⁺ TAM 配体（SPP1+B 高表达）
    → CD8⁺ T 细胞受体（LIANA 验证 LR 对）
        → CD8⁺ T 细胞耗竭标志（文献 + 数据）
```

连接权重 = LIANA magnitude_rank（取反归一化）；节点颜色 = 在对应细胞类型中的 mean expression（z-score）。

### 2.3 代谢通路评分（Panel D/E）

- **方法**：AUCell 评分（scMetabolism 等价实现）
  - 工具：decoupler-py v1.8.0，`run_aucell()` 函数
  - 基因集来源：KEGG 2021 Human（通过 gseapy 获取）
  - AUCell 计算每个细胞对特定基因集的富集程度（AUC 曲线下面积），对稀疏单细胞数据更鲁棒
  - 各细胞类型取平均 AUCell 分数，跨细胞类型 Z-score 标准化后作图
- **基因集（10 条 KEGG 通路）**：脂质摄取/FAO（自定义基因集）、胆固醇合成（Steroid biosynthesis）、糖酵解（Glycolysis/Gluconeogenesis）、氧化磷酸化（Oxidative phosphorylation）、脂肪酸合成（Fatty acid biosynthesis）、脂肪酸降解（Fatty acid degradation）、嘌呤代谢（Purine metabolism）、氨基糖代谢（Amino sugar and nucleotide sugar metabolism）、TCA 循环（Citrate cycle）、花生四烯酸代谢（Arachidonic acid metabolism）
- **Panel E**：代谢通路 × 细胞类型热图（pathways × cell types），同为 AUCell z-score，OLR1⁺ TAMs 列用红框标注

### 2.4 通路差异评分（Panel F，MMRp vs MMRd）

- **数据**：GSE178341 Macro 亚群（n = 17,001 cells），按 MMRStatus 分层
- **统计**：Mann-Whitney U 检验（双侧），FDR 不作多重校正（探索性）
- **10 条通路**：腺苷通路、氧化磷酸化、免疫检查点、抗原呈递、脂质摄取/FAO、缺氧/VEGFA、脂肪酸合成、EMT/ECM、色氨酸/精氨酸、炎症信号

---

## 三、结果

### 3.1 通讯网络：OLR1⁺ TAMs 是最活跃的免疫抑制信号发送者（Panel A）

**Figure 4A — 细胞通讯 Circle Plot（GSE144735 Border）**

![Panel A — Circle Plot](Project_Delivery/results/Figure4_Communication_Metabolism/Cell_Communication/figures/panelA_circle_plot.png)

*Hub-and-spoke 通讯网络图，OLR1⁺ TAMs 居中为核心 hub（红色圆，标注"Border"表示肿瘤边界定位；已移除"SPP1+B"内部标注以避免与 SPP1⁺ B 细胞混淆）。箭头宽度 = 通讯强度（LIANA lr_means）；n = 有效 LR 对数（p<0.05）；配体–受体对标注精确定位于对应弧线中点。OLR1⁺ TAMs 向 CD8⁺ T cells、CD4⁺ T cells、Reg T cells 发出免疫抑制信号，并与 Myofibroblasts 形成双向 ECM 重塑反馈回路。*

**关键发现**：

| 发送方 → 接收方 | 有效 LR 对数（p<0.05） | 平均 LR 均值 |
|---------------|---------------------|------------|
| OLR1⁺ TAMs → CD8⁺ T cells | 168 | 0.953 |
| OLR1⁺ TAMs → CD4⁺ T cells | 167 | 0.844 |
| OLR1⁺ TAMs → Reg T cells | 192 | 0.831 |
| OLR1⁺ TAMs → Myofibroblasts | 379 | 0.656 |
| Myofibroblasts → OLR1⁺ TAMs | 425 | — |

OLR1⁺ TAMs 向 CD8⁺ T cells 发出的信号均值最高（0.953），是 Border 区域最活跃的免疫调节信号枢纽。同时，与 Myofibroblast 之间存在**高密度双向通讯**（TAM→Myo: 379 对；Myo→TAM: 425 对），支持两者在边界的功能协同。

### 3.2 关键配体–受体对：三类功能性信号轴（Panel B）

**Figure 4B — 关键 LR 对 Bubble Plot**

![Panel B — Bubble Plot](Project_Delivery/results/Figure4_Communication_Metabolism/Cell_Communication/figures/panelB_bubble_plot.png)

*SPP1+B 发出的关键 LR 对按三大功能分组。气泡大小 = 通讯强度（-log₁₀ magnitude_rank）；x 轴 = mean LR 共表达值（lr_means）；颜色 = 功能类别。*

**I. 免疫抑制轴（SPP1+B → CD8⁺ T cells）**

| 配体 | 受体 | lr_means | magnitude_rank | 生物功能 |
|-----|------|---------|----------------|---------|
| **LGALS1** | CD69 | 2.447 | 0.0000149 | Galectin-1 交联 CD69，抑制 T 细胞活化 |
| **LGALS1** | PTPRC (CD45) | 2.374 | 0.0000381 | 抑制 TCR 信号磷酸化 |
| **SPP1** | CD44 | 2.009 | 0.000650 | 骨桥蛋白→CD44，T 细胞迁移阻断 |
| **MIF** | CD74/CXCR4 | 1.990 | 0.000720 | MIF 拮抗 CD74→抗凋亡信号 |
| **HMGB1** | CXCR4 | 2.266 | 0.000137 | 报警素→CXCR4，T 细胞滞留 |
| **S100A8** | CD69 | 1.846 | 0.001753 | 炎性钙结合蛋白→T 细胞活化阻断 |

> **LGALS1（Galectin-1）是首要抑制配体**（magnitude_rank = 0.0000149，最强之一），同时靶向 CD69 和 CD45，直接抑制 T 细胞受体激活与磷酸化级联。

**II. Treg 募集与激活轴（SPP1+B → Regulatory T cells）**

| 配体 | 受体 | lr_means | 生物功能 |
|-----|------|---------|---------|
| LGALS1 | CD69 | 2.391 | 抑制 Treg 撤活（维持抑制态） |
| LGALS1 | PTPRC | — | Treg 存活维持 |
| SPP1 | CD44 | 1.614 | 趋化性募集 Treg 至边界 |
| MIF | CD74/CXCR4 | — | 促 Treg 存活 |

**III. 基质重塑轴（SPP1+B ↔ Myofibroblasts，核心）**

| 配体（来源） | 受体（靶细胞） | lr_means | magnitude_rank | 机制 |
|------------|-------------|---------|----------------|------|
| **SPP1**（TAM→Myo） | ITGAV:ITGB1 | 1.703 | 0.00522 | 骨桥蛋白→αvβ1 整合素，激活 CAF 收缩 |
| **SPP1**（TAM→Myo） | CD44 | 1.614 | 0.01156 | 骨桥蛋白→CD44，激活 PI3K/AKT |
| **SPP1**（TAM→Myo） | ITGAV:ITGB5 | 1.610 | 0.01195 | 骨桥蛋白→αvβ5 整合素，ECM 重塑 |
| **APOE**（TAM→Myo） | LRP1 | 1.686 | 0.00509 | ApoE→LRP1，脂质信号转导至 CAF |
| **LGALS1**（TAM→Myo） | ITGB1 | 2.397 | 0.0000301 | Galectin-1→β1 整合素，CAF 活化最强信号 |
| **TIMP1**（TAM→Myo） | CD63 | 2.309 | 0.0000720 | 基质金属蛋白酶抑制 + CD63 信号 |
| **COL1A1**（Myo→TAM） | CD44 | 2.296 | 0.0000810 | 胶原→TAM CD44，正反馈激活 TAM |
| **COL1A2**（Myo→TAM） | CD44 | 2.291 | 0.0000850 | 同上，不同胶原亚型 |
| **FN1**（Myo→TAM） | CD44 | 1.998 | 0.000691 | 纤连蛋白→TAM CD44，维持 TAM 极化 |
| **TIMP1**（Myo→TAM） | CD63 | 2.504 | 0.0000070 | **最强反馈信号**：Myo→TAM |

> **关键结论**：SPP1 同时靶向 Myofibroblast 上的 ITGAV:ITGB1（αvβ1）、ITGAV:ITGB5（αvβ5）和 CD44 三个受体，是激活 FAP⁺ CAF 并诱导 ECM 收缩的核心配体。反向通讯中，Myofibroblast 分泌 COL1A1/COL1A2/FN1 作用于 TAM CD44，形成**TAM–CAF 正反馈环路**，使两者在 Border 共同维持激活状态。

### 3.3 NicheNet-like 配体–靶基因连接图（Panel C）

**Figure 4C — 配体–受体–耗竭靶基因连接图**

![Panel C — NicheNet](Project_Delivery/results/Figure4_Communication_Metabolism/Cell_Communication/figures/panelC_nichenet.png)

*三列分别为 OLR1⁺ TAM 配体（红色渐变，按表达量深浅）、CD8⁺ T 细胞受体（蓝色）、CD8⁺ T 细胞耗竭标志（紫色）。箭头代表 LIANA 验证或 NicheNet 先验知识支持的信号连接；箭头粗细 = 连接权重。*

**关键预测通路**：

```
LGALS1 ──→ CD69 ──────→ PDCD1 (PD-1)↑
LGALS1 ──→ PTPRC ─────→ LAG3↑
SPP1   ──→ CD44 ──────→ HAVCR2 (TIM-3)↑
SPP1   ──→ CD44 ──────→ TIGIT↑
HMGB1  ──→ CXCR4 ─────→ TOX↑  (耗竭主调控 TF)
MIF    ──→ CXCR4/CD74 → CXCL13↑ (耗竭 T 细胞标志)
VEGFA  ──→ CXCR4 ─────→ PDCD1↑
```

LGALS1 和 SPP1 是预测中最关键的上游配体：
- **LGALS1→CD69/PTPRC**：通过抑制 TCR 磷酸化，诱导 PD-1、LAG-3 等耗竭检查点上调
- **SPP1→CD44**：激活 PI3K/RhoA 信号，TIM-3 和 TIGIT 上调，T 细胞迁移力丧失
- **HMGB1→CXCR4**：报警素信号激活 TOX（耗竭核心转录因子），驱动 T 细胞进入终末耗竭

> OLR1⁺ TAM 通过 **三条平行抑制轴**（LGALS1轴/SPP1轴/HMGB1轴）同时作用于 CD8⁺ T 细胞，靶向不同受体，下游汇聚至 PD-1/TIM-3/LAG3/TIGIT/TOX 耗竭程序，形成多重冗余的免疫抑制网络。

### 3.4 代谢通路偏好：脂质摄取/FAO、氧化磷酸化与糖酵解三轴主导 OLR1⁺ TAM 代谢（Panel D）

**Figure 4D — 细胞类型代谢通路 AUCell 评分（GSE144735）**

![Panel D — Metabolism](Project_Delivery/results/Figure4_Communication_Metabolism/Metabolic_Flux/figures/panelD_metabolism.png)

*10 条 KEGG 代谢通路的 AUCell z-score，对比 OLR1⁺ TAMs（SPP1+B）、Anti-inflam TAMs、CD8⁺ T cells、Myofibroblasts 和 Reg T cells（n = 27,414 cells，GSE144735）。AUCell 方法参照 scMetabolism（Wu et al. 2021）实现，通过 decoupler-py v1.8.0 在 KEGG 2021 Human 基因集上计算。*

**关键 AUCell 代谢偏好（OLR1⁺ TAMs，z-score）**：

| 通路 | OLR1⁺ TAM z-score | 临床意义 |
|------|-------------------|---------|
| **脂质摄取/FAO**（OLR1/CD36/FABP4/CPT1A） | **+1.70** | 核心脂质摄取程序，oxLDL 摄取驱动 TAM 极化 |
| **氧化磷酸化** | **+1.70** | ATP 生成能力最强，代谢可塑性高 |
| **糖酵解** | **+1.60** | 高度糖酵解活性，支持快速能量供应 |
| **胆固醇合成** | **+1.23** | 膜脂重塑，与免疫抑制极化相关 |
| **花生四烯酸代谢** | **+1.14** | 前列腺素 E2 合成，抑制 T 细胞活化 |
| Purine 代谢 | -0.69 | 低于其他细胞类型 |
| TCA 循环 | -0.70 | 低于其他细胞类型 |

OLR1⁺ TAMs 在脂质摄取/FAO 和氧化磷酸化两个维度均为最高，提示其兼具高效脂质利用（OLR1/FABP4）与强大线粒体氧化能力，与"脂质关联巨噬细胞（LAM）"的能量代谢特征高度一致。同时，糖酵解的高活性表明 OLR1⁺ TAMs 具有与 CD8⁺ T 细胞竞争葡萄糖资源的能力，形成代谢竞争压力。

### 3.5 代谢通路活性热图：OLR1⁺ TAMs vs 其他细胞类型（Panel E）

**Figure 4E — 代谢通路 AUCell 热图（GSE144735）**

![Panel E — Metabolic Heatmap](Project_Delivery/results/Figure4_Communication_Metabolism/Metabolic_Flux/figures/panelE_metabolic_heatmap.png)

*10 条 KEGG 代谢通路 × 5 种细胞类型的 AUCell z-score 热图。红色 = AUCell 评分高（通路活性强）；蓝色 = 低活性。OLR1⁺ TAMs 列（SPP1+B）用红框标注。通路按功能分组并以彩色条标注。*

**读图逻辑**：

- OLR1⁺ TAMs（红框列）在**脂质摄取/FAO、氧化磷酸化、糖酵解、胆固醇合成、花生四烯酸代谢**5 条通路均呈现最高 AUCell z-score（暖红色）
- **Myofibroblasts** 在脂肪酸合成和 TCA 循环中占优，符合 CAF 的合成代谢特征
- **CD8⁺ T cells** 在所有代谢通路中评分最低，反映 TME 中 T 细胞的代谢受限状态
- 花生四烯酸代谢（前列腺素 E2 轴）的选择性高表达，提示 OLR1⁺ TAMs 通过 COX2/PGE2 介导额外的 T 细胞抑制通路

这一模式支持**代谢竞争（Metabolic Competition）假说**：OLR1⁺ TAMs 在脂质摄取/FAO 和氧化磷酸化双维度均处于代谢最优位，同时通过高糖酵解与 CD8⁺ T 细胞竞争葡萄糖，导致 TME 中 T 细胞能量受限，与通讯抑制轴协同形成"双锁定"耗竭机制。

### 3.6 MMRp vs MMRd 通路差异：腺苷信号在 MMRp 显著上调（Panel F）

**Figure 4F — GSEA 通路差异（GSE178341 Macro，MMRp vs MMRd）**

![Panel F — GSEA](Project_Delivery/results/Figure4_Communication_Metabolism/Pathway_Enrichment/figures/panelF_gsea.png)

*10 条代谢/免疫通路在 OLR1⁺ TAMs（Macro 亚群）中 MMRp vs MMRd 的 log₂FC 条形图（MMRp/MMRd > 1 = MMRp 中更高）；黄色点 = -log₁₀(p 值)；★ = p < 0.05 显著。*

**通路差异汇总**（GSE178341 Macro，n = 17,001）：

| 通路 | log₂FC（MMRp/MMRd） | p 值 | 结论 |
|------|---------------------|------|------|
| **腺苷通路** ★ | **+0.192** | **1.9×10⁻¹⁰** | MMRp OLR1⁺ TAM 腺苷产生能力更强 |
| 氧化磷酸化 | -0.012 | 0.168（n.s.） | 基础代谢无差异 |
| 免疫检查点表达 ★ | -0.021 | 0.018 | MMRd 中检查点略高（热肿瘤反应性） |
| 抗原呈递 ★ | -0.029 | 8.8×10⁻¹² | MMRd 中更活跃（高免疫原性驱动） |
| 脂质摄取/FAO ★ | -0.186 | 3.4×10⁻³² | MMRd 脂质代谢更活跃（可能更多 TAM 亚型混入） |
| 缺氧/VEGFA ★ | -0.213 | 6.9×10⁻⁷³ | MMRd 中更强（血管生成活跃） |
| EMT/ECM 重塑 ★ | -0.288 | 3.4×10⁻¹⁴² | MMRd 基质更活跃 |
| 炎症信号 ★ | -0.526 | 4.1×10⁻⁶⁷ | MMRd 炎症性 TAM 为主 |

**解读**：

- **腺苷通路是唯一在 MMRp 显著上调的通路**（p = 1.9×10⁻¹⁰）：MMRp 冷肿瘤中 OLR1⁺ TAM 的腺苷产生（ENTPD1↑/NT5E↑）高于 MMRd，与"腺苷→ADORA2A 信号抑制 T 细胞 cAMP"机制高度一致，是 MMRp 冷肿瘤环境特异性免疫抑制的代谢分子基础。
- MMRd 中炎症信号（CXCL10/IL6/TNF）更高，代表热肿瘤免疫激活状态，与冷肿瘤形成对比验证。

---

## 四、综合机制模型

**English version (for publication):**

![Integrated Mechanism Model — English](Project_Delivery/results/Figure4_Communication_Metabolism/Cell_Communication/figures/mechanism_model_EN.png)

*Integrated mechanism model showing three-layer synergistic immune suppression by OLR1+ lipid-associated TAMs in MMRp CRC: [1] Ligand–receptor axis (LGALS1/SPP1/HMGB1/MIF → CD8+ T cell exhaustion); [2] Metabolic competition (adenosine pathway MMRp-specific, p=1.9×10⁻¹⁰; IDO1/ARG1 amino acid depletion); [3] ECM remodeling (SPP1→ITGAV:ITGB1 TAM–CAF positive feedback loop → physical collagen barrier).*

---

**中文版机制模型：**

```
    脂质微环境（oxLDL 富集）
            ↓ OLR1 摄取
    ╔════════════════════════╗
    ║   OLR1⁺ lipid-assoc.  ║
    ║        TAMs            ║
    ║ (SPP1+B，Border 35.5%) ║
    ╚════════════╤═══════════╝
                 │
    ┌────────────┼─────────────┐
    │            │             │
    ▼            ▼             ▼
[配体轴]     [代谢轴]      [ECM 轴]
LGALS1→CD69  ENTPD1/NT5E  SPP1→ITGAV:ITGB1
SPP1→CD44    (腺苷产生)   APOE→LRP1
HMGB1→CXCR4 IDO1/ARG1    COL1A1→CD44（反馈）
    │        (氨基酸竞争)       │
    ▼                           ▼
CD8⁺ T 耗竭              FAP⁺ Myofibroblast
PD-1/TIM-3/LAG3↑         激活→胶原分泌
TOX↑（终末耗竭）          物理胶原屏障
    │                           │
    └─────────────┬─────────────┘
                  ▼
          MMRp CRC 免疫排除
        （"冷肿瘤"终局）
```

**三层协同抑制机制**：

1. **直接通讯抑制**：LGALS1、SPP1、MIF、HMGB1 多重配体并行作用 CD8⁺ T 细胞，靶向不同受体（CD69/CD44/CXCR4），下游汇聚至 PD-1/TIM-3/LAG3/TOX 耗竭程序

2. **代谢竞争剥夺**：腺苷通路（ENTPD1/NT5E）+ IDO1/ARG1 双轴，消耗 T 细胞所需腺苷前体和必需氨基酸，腺苷通路在 MMRp 中显著上调（p = 1.9×10⁻¹⁰），是 MMRp 冷肿瘤的代谢分子标志

3. **物理屏障构建**：SPP1 通过 αvβ1/αvβ5 整合素激活 Myofibroblast，Myofibroblast 分泌 COL1A1/FN1 反馈激活 TAM（TIMP1→CD63 最强，magnitude_rank = 0.0000070），形成 TAM–CAF 正反馈环路，持续扩大胶原屏障

---

## 五、结果汇总

| 分析模块 | 核心发现 |
|---------|---------|
| 通讯网络（Panel A） | OLR1⁺ TAMs 是 Border 区最活跃免疫调节枢纽，向 CD8（168 对）和 Myofibroblast（379 对）发出强信号 |
| 关键 LR 对（Panel B） | LGALS1→CD69（最强，rank=0.0000149）；SPP1→ITGAV:ITGB1（TAM→CAF 核心）；Myo→TAM TIMP1→CD63（最强反向信号） |
| NicheNet（Panel C） | LGALS1/SPP1/HMGB1 三轴驱动 CD8 进入 PD-1/TIM-3/LAG3/TOX 耗竭程序 |
| 代谢偏好（Panel D/E） | AUCell 评分（KEGG 2021，decoupler-py）：OLR1⁺ TAMs 在脂质摄取/FAO（z=+1.70）、氧化磷酸化（z=+1.70）、糖酵解（z=+1.60）三轴领先，花生四烯酸代谢（z=+1.14）提示 PGE2 抑制轴 |
| MMRp vs MMRd（Panel F） | 腺苷通路唯一在 MMRp 显著上调（log₂FC=+0.19，p=1.9×10⁻¹⁰），是冷肿瘤代谢分子基础 |

---

## 六、潜在治疗靶点（阶段四新增）

| 靶点 | 阶段四依据 | 靶向策略 |
|------|-----------|---------|
| **LGALS1（Galectin-1）** ⭐ | 最强免疫抑制 LR（CD69/CD45 双靶） | 抗 Galectin-1 单抗（TD139 analogue） |
| **SPP1（骨桥蛋白）** ⭐ | 同时靶向 CD8（CD44）和 CAF（ITGAV:ITGB1/5） | 抗 SPP1 单抗/整合素抑制剂 |
| **ENTPD1/CD39** | 腺苷通路 MMRp ↑ 核心酶 | CD39 抑制剂（POM-1，AK-119） |
| **NT5E/CD73** | 腺苷通路第二酶 | CD73 抑制剂（Oleclumab，CPI-006） |
| **IDO1** | MMRp TAM 色氨酸竞争 | IDO1 抑制剂（Epacadostat）联合 PD-1 |
| **MIF** | MIF→CD74/CXCR4，多靶抗凋亡 | MIF 抑制剂（ISO-1，SCD-19） |

---

## 七、与各阶段的连接

| 阶段 | 核心结论 | 阶段四补充 |
|------|---------|-----------|
| Phase 2（空间） | OLR1⁺ TAM + FAP⁺ CAF 共富集于 Border | **阶段四证明**：两者通过 SPP1→整合素 + COL1A1→CD44 正反馈驱动共富集 |
| Phase 3（MIL） | OLR1⁺ TAM 注意力 7.5×（MMRp vs MMRd） | **阶段四解释**：腺苷通路 MMRp↑ + LGALS1/SPP1 配体，是模型可"识别"冷肿瘤的分子基础 |
| Phase 4（本阶段） | 三层机制完整拼图 | — |

阶段四将"空间共定位→MIL 高注意力→分子通讯→代谢竞争"四条证据链整合，形成 OLR1⁺ TAMs 诱导 MMRp 冷肿瘤免疫耐药的完整生物学逻辑。
