##########################################

### health Brain sample analysis

### Begin 2025-4-2

### PMID: 39962241

### Author：honghui Zhao

##########################################

### load packages for analysis
library(Seurat)
library(tidyverse)
library(harmony)
library(scater)
library(reticulate)
library(clustree)
library(sctransform)
library(Cairo)
library(tidydr)
library(dplyr)
library(ggplot2)
library(cols4all)
options(bitmapType='cairo')

### setwd path
setwd('/sibcb1/douxiaoyanglab1/zhaohonghui/Glioma/health_barin')

pngpath <- '/sibcb1/douxiaoyanglab1/zhaohonghui/Glioma/health_barin/PDF_PNG'


sce <- readRDS('/sibcb1/douxiaoyanglab1/zhaohonghui/Glioma/health_barin/intergated_snRNA.rds')


#### basic plot

correct_order <- c("0-year-old stage", "4-year-old stage", "6-year-old stage",
                   "14-year-old stage", "20-year-old stage", "39-year-old stage",
                   "61-year-old stage", "62-year-old stage")


correct_order <- c(
  # 1号样本
  "postnatB_1_anterior cingulate cortex",
  "postnatB_1_caudate nucleus",
  "postnatB_1_dorsolateral prefrontal cortex",
  "postnatB_1_hippocampal formation",
  
  # 2号样本
  "postnatB_2_anterior cingulate cortex",
  "postnatB_2_caudate nucleus",
  "postnatB_2_dorsolateral prefrontal cortex",
  "postnatB_2_hippocampal formation",
  
  # 3号样本
  "postnatB_3_anterior cingulate cortex",
  "postnatB_3_dorsolateral prefrontal cortex",
  "postnatB_3_hippocampal formation",
  
  # 4号样本
  "postnatB_4_anterior cingulate cortex",
  "postnatB_4_caudate nucleus",
  "postnatB_4_dorsolateral prefrontal cortex",
  "postnatB_4_hippocampal formation",
  
  # 5号样本
  "postnatB_5_anterior cingulate cortex",
  "postnatB_5_caudate nucleus",
  "postnatB_5_dorsolateral prefrontal cortex",
  "postnatB_5_hippocampal formation",
  
  # 6号样本
  "postnatB_6_anterior cingulate cortex",
  "postnatB_6_caudate nucleus",
  "postnatB_6_dorsolateral prefrontal cortex",
  "postnatB_6_hippocampal formation",
  
  # 7号样本
  "postnatB_7_anterior cingulate cortex",
  "postnatB_7_caudate nucleus",
  "postnatB_7_dorsolateral prefrontal cortex",
  "postnatB_7_hippocampal formation",
  
  # 8号样本
  "postnatB_8_anterior cingulate cortex",
  "postnatB_8_caudate nucleus",
  "postnatB_8_dorsolateral prefrontal cortex",
  "postnatB_8_hippocampal formation",
  
  # 9号样本
  "postnatB_9_anterior cingulate cortex",
  "postnatB_9_caudate nucleus",
  "postnatB_9_dorsolateral prefrontal cortex",
  "postnatB_9_hippocampal formation",
  
  # 10号样本
  "postnatB_10_anterior cingulate cortex",
  "postnatB_10_caudate nucleus",
  "postnatB_10_dorsolateral prefrontal cortex",
  "postnatB_10_hippocampal formation"
)

qc_plots <- VlnPlot(
  object = sce,raster=FALSE,
  features = c("percent.rb"),# "nCount_RNA", "percent.mt","percent.rb"),
  group.by = "donor_tissue",          # 按分组列（需替换为您的列名）
  pt.size = 0,                 # 不显示点（避免过度密集）
  ncol = 1                   # 一行显示3个图
) +  guides(fill = "none")+
  theme(
    axis.title.x = element_blank(),  # 隐藏 x 轴标题
    axis.text.x = element_text(angle = 90, hjust = 1, size = 10),  # 旋转x轴标签
    strip.text = element_blank(),    # 移除分面标签（如 "nFeature_RNA"）
    strip.background = element_blank()  # 移除分面背景（可选）
  )

 for (i in seq_along(1)) {
   qc_plots[[i]][["data"]][["ident"]] <- factor(
     qc_plots[[i]][["data"]][["ident"]],
     levels = correct_order
   )
 }


# 保存为PDF或PNG（调整尺寸和分辨率）
ggsave("percent.rb_donor_tissue_QC.pdf", 
       plot = qc_plots, path = pngpath,
       width = 12, height = 6, dpi = 300)


### UMAP plot

umap_plots <- DimPlot(sce,raster=FALSE,reduction = 'umap_sct')
ggsave("umap_sct.pdf", 
       plot = umap_plots, path = pngpath,
       width = 6, height = 6, dpi = 300)

### clustree
clustree_plots <- clustree(sce@meta.data, prefix = "SCT_snn_res.")
ggsave("SCT_clustree.pdf", 
       plot = clustree_plots, path = pngpath,
       width = 10, height = 10, dpi = 300)


#### umap-sct
### UMAP plot

umapsct_plots <- DimPlot(sce,raster=FALSE,reduction = 'umap_sct',group.by = 'SCT_snn_res.1')
ggsave("umapsct_SCT_snn_res.1.pdf", 
       plot = umapsct_plots, path = pngpath,
       width = 6, height = 6, dpi = 300)


### cell percentage

# 生成计数矩阵并转换为按列（患者）的比例
cell_prop_df <- as.data.frame.table(
  prop.table(table(sce@meta.data$celltype, sce@meta.data$`IDH status`), margin = 2) * 100,
  responseName = "Proportion"
) 

colnames(cell_prop_df)[1] <- 'celltype'
colnames(cell_prop_df)[2] <- 'Mutation'

ggsave('cellprop_Mutation_bar.pdf',width = 8,height = 5)
ggplot(cell_prop_df, aes(x = Mutation, y = Proportion, fill = celltype)) +
  geom_col(color = "black", width = 0.8) +
  scale_fill_brewer(palette = "Set3") +
  labs(
    title = "Cell Type Composition by Mutation",
    x = "",  # 这里设置为空字符串即可移除X轴标签
    y = "Percentage (%)"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),  # 标题居中加粗
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
    legend.position = "right",
    legend.title = element_blank(),
    panel.grid.major.x = element_blank()
  )
dev.off()

ggsave('dimplot_split_patinet.pdf',width = 18,height = 6)
DimPlot(sce,label = T,split.by = 'patient_id') + NoLegend()
dev.off()
