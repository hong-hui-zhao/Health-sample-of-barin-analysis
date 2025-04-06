##########################################

### health Brain sample analysis

### Begin 2025-4-2

### PMID: 39962241

### Author：honghui Zhao

##########################################


library(Seurat)
library(Matrix)
library(tidyverse)
library(lisi)
library(readr)
library(ggplot2)
library(ggpubr)
options(bitmapType='cairo')
donor_tissue <- readRDS('/sibcb1/douxiaoyanglab1/zhaohonghui/Glioma/health_barin/harmony_donor_tissue_snRNA.rds')
donor <- readRDS('/sibcb1/douxiaoyanglab1/zhaohonghui/Glioma/health_barin/harmony_donor_snRNA.rds')
DimPlot(donor_tissue,raster=FALSE,reduction = 'umap_sct',group.by = 'SCT_snn_res.1')
DimPlot(donor,raster=FALSE,reduction = 'umap_sct',group.by = 'SCT_snn_res.1')


# 计算 donor_tissue 数据的 LISI 分数
# ILISI（批次效应评估）
ilisi_donor_tissue <- lisi::compute_lisi(
  X = Embeddings(donor_tissue, reduction = "harmony"),
  meta_data = donor_tissue@meta.data,
  label_colnames = "donor_tissue"
) %>% 
  dplyr::rename(ILISI = donor_tissue)

# CLISI（生物学效应评估）
clisi_donor_tissue <- lisi::compute_lisi(
  X = Embeddings(donor_tissue, reduction = "harmony_sct"),
  meta_data = donor_tissue@meta.data,
  label_colnames = "SCT_snn_res.1"
) %>% 
  dplyr::rename(CLISI = SCT_snn_res.1)

# 合并结果并添加数据集标签
donor_tissue_lisi <- dplyr::mutate(
  ilisi_donor_tissue,
  CLISI = clisi_donor_tissue$CLISI,
  dataset = "donor_tissue"
)

# 计算 donor 数据的 LISI 分数
# ILISI（批次效应评估）
ilisi_donor <- lisi::compute_lisi(
  X = Embeddings(donor, reduction = "harmony_sct"),
  meta_data = donor@meta.data,
  label_colnames = "donor_id"
) %>% 
  dplyr::rename(ILISI = donor_id)

# CLISI（生物学效应评估）
clisi_donor <- lisi::compute_lisi(
  X = Embeddings(donor, reduction = "harmony_sct"),
  meta_data = donor@meta.data,
  label_colnames = "SCT_snn_res.1"
) %>% 
  dplyr::rename(CLISI = SCT_snn_res.1)

# 合并结果并添加数据集标签
donor_lisi <- dplyr::mutate(
  ilisi_donor,
  CLISI = clisi_donor$CLISI,
  dataset = "donor"
)


ilisi_Nor <- lisi::compute_lisi(
  X = Embeddings(Nor, reduction = "harmony"),
  meta_data = Nor@meta.data,
  label_colnames = "donor_tissue"
) %>% 
  dplyr::rename(ILISI = donor_tissue)

# CLISI（生物学效应评估）
clisi_nor <- lisi::compute_lisi(
  X = Embeddings(Nor, reduction = "harmony"),
  meta_data = Nor@meta.data,
  label_colnames = "SCT_snn_res.1"
) %>% 
  dplyr::rename(CLISI = SCT_snn_res.1)

# 合并结果并添加数据集标签
donor_nor <- dplyr::mutate(
  ilisi_Nor,
  CLISI = clisi_nor$CLISI,
  dataset = "Nor"
)


# 合并两个数据集的结果
combined_lisi <- dplyr::bind_rows(donor_tissue_lisi, donor_lisi)

# 转换为长格式用于绘图
combined_lisi_long <- tidyr::pivot_longer(
  combined_lisi,
  cols = c(ILISI, CLISI),
  names_to = "Metric",
  values_to = "LISI_Score"
)

# 绘制分面箱型图
library(ggplot2)
ggplot(combined_lisi_long, aes(x = dataset, y = LISI_Score, fill = dataset)) +
  geom_boxplot(width = 0.6, outlier.size = 0.5) +
  facet_wrap(~Metric, scales = "free_y", ncol = 2) +
  scale_fill_manual(values = c("#E69F00", "#56B4E9")) +
  theme_bw(base_size = 12) +
  labs(x = "", y = "LISI Score", title = "LISI Score Comparison") +
  theme(legend.position = "top",
        strip.background = element_blank(),
        panel.grid.minor = element_blank())
ggsave('LISI_donor_tissue vs donor.pdf',width = 8,height = 4)


### 三个数据集 ------------------------
# 合并三个数据集的结果（添加Nor数据集）
combined_lisi <- dplyr::bind_rows(donor_tissue_lisi, donor_lisi, donor_nor)

# 转换为长格式用于绘图
combined_lisi_long <- tidyr::pivot_longer(
  combined_lisi,
  cols = c(ILISI, CLISI),
  names_to = "Metric",
  values_to = "LISI_Score"
)

# 绘制分面箱型图（适配三组比较）
ggplot(combined_lisi_long, 
       aes(x = dataset, y = LISI_Score, fill = dataset)) +
  geom_boxplot(width = 0.6, outlier.size = 0.5) +
  facet_wrap(~Metric, scales = "free_y", ncol = 2) +
  scale_fill_manual(
    values = c("#E69F00", "#56B4E9", "#009E73"), # 新增第三个颜色
    labels = c("Donor", "Donor_tissue", "Nor")    # 自定义图例标签
  ) +
  theme_bw(base_size = 12) +
  labs(x = "", y = "LISI Score", 
       title = "LISI Score Comparison (Three Datasets)") +
  theme(
    legend.position = "top",
    legend.title = element_blank(),  # 隐藏图例标题
    strip.background = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1) # 倾斜X轴标签
  ) +
  stat_compare_means(
    method = "wilcox.test",
    comparisons = list(
      c("donor", "donor_tissue"),
      c("donor", "Nor"),
      c("donor_tissue", "Nor")
    ),
    label = "p.signif"  # 显示显著性符号
  )

ggsave('LISI_three_datasets.pdf', width = 10, height = 6)
