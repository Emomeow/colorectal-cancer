#!/bin/bash

# 1. 创建顶层目录
mkdir -p Project_Delivery/{data/{raw_data,processed_data,metadata},script,results}

# 2. 批量生成 Figure 1 - 6 及 Supplementary 目录及其子结构
cd Project_Delivery/results

# Figure 1
mkdir -p Figure1_scVI_Atlas/{script,scVI_Integration/{figures,tables},Cell_Composition/{figures,tables},Marker_Expression/{figures,tables},Cluster_Stability/{figures,tables}}

# Figure 2
mkdir -p Figure2_Pre_driver_Trajectory/{script,Trajectory_Analysis/{figures,tables},Spatial_Distribution/{figures,tables},SCENIC_Regulon/{figures,tables}}

# Figure 3
mkdir -p Figure3_MIL_MMR_Prediction/{script,MIL_Model_Performance/{figures,tables},Attention_Analysis/{figures,tables},Augur_Predictability/{figures,tables},Transfer_Validation/{figures,tables}}

# Figure 4
mkdir -p Figure4_Communication_Metabolism/{script,Cell_Communication/{figures,tables},Metabolic_Flux/{figures,tables},Pathway_Enrichment/{figures,tables}}

# Figure 5
mkdir -p Figure5_Anatomical_Heterogeneity/{script,Anatomical_Distribution/{figures,tables},Transcriptional_Differences/{figures,tables},Clinical_Correlation/{figures,tables}}

# Figure 6
mkdir -p Figure6_Bulk_Clinical_Validation/{script,Deconvolution/{figures,tables},Survival_Analysis/{figures,tables}}

# Supplementary
mkdir -p Supplementary/{script,Supp1_QC_Metrics/{figures,tables},Supp2_Extended_Markers/{figures,tables},Supp3_MIL_Robustness/{figures,tables},Supp4_Normal_Control/{figures,tables},Supp5_Metabolism_Comm/{figures,tables}}

echo "Project_Delivery 目录结构创建完成！"
