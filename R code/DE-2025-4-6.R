##########################################

### health Brain sample analysis

### Begin 2025-4-2

### PMID: 39962241

### Author：honghui Zhao

##########################################


library(Seurat)
library(tidyverse)

### load data

# 设置工作路径
work_dir <- '/sibcb1/douxiaoyanglab1/zhaohonghui/Glioma/health_barin'
setwd(work_dir)

# 定义输出目录
result_dir <- '/sibcb1/douxiaoyanglab1/zhaohonghui/Glioma/health_barin/PDF_PNG'
if (!dir.exists(result_dir)) dir.create(result_dir, recursive = TRUE)

# 加载Seurat对象
sce <- readRDS('harmony_donor_tissue_snRNA.rds')
BiocManager::install("MAST")

####  DE gene

# 加载必要的包
library(Seurat)
library(ggplot2)
library(dplyr)
library(purrr)
library(MAST)
# 定义阶段顺序
stage_order <- c("infant", "Middle childhood", "Adolescence",
                 "Early adulthood", "Late adulthood")

Idents(sce) <- "life_stage"
table(sce@meta.data$donor_id,sce@meta.data$sex)
# 初始化一个空列表用于存储所有差异分析结果
all_degs <- list()
sce <- PrepSCTFindMarkers(sce, assay = "SCT", verbose = TRUE)
# 循环执行相邻阶段的差异分析
for (i in 1:(length(stage_order)-1)) {
  # 获取当前要比较的两个阶段
  group1 <- stage_order[i]
  group2 <- stage_order[i+1]
  comparison_name <- paste0(group2, "_vs_", group1)
  
  # 执行差异分析（使用默认参数，可根据需要调整）
  degs <- FindMarkers(sce,
                      ident.1 = group2,
                      ident.2 = group1,
                      logfc.threshold = 0.2, test.use = "MAST",
                      min.diff.pct = 0.1, 
                      )         
  
  # 添加基因名列（FindMarkers结果默认以基因名为行名）
  degs <- degs %>% mutate(gene = rownames(.))
  
  # 保存结果到列表
  all_degs[[comparison_name]] <- degs
  
  # 保存为CSV文件（替换空格为下划线）
  filename <- gsub(" ", "_", paste0(comparison_name, "_DEGs_MAST.csv"))
  write.csv(degs, file = filename, row.names = FALSE)
}

# 筛选持续变化的基因 -----------------------------------------------------------
### 筛选持续变化基因 ###
library(dplyr)
library(purrr)

logFC_list <- map(names(all_degs), ~ {
  comp_name <- .x
  all_degs[[comp_name]] %>% 
    filter(p_val_adj < 0.05) %>%   # 添加显著性过滤
    select(gene, avg_log2FC) %>% 
    setNames(c("gene", comp_name))
})

# 2. 合并所有比较数据（取基因交集）
combined_logFC <- purrr::reduce(logFC_list, inner_join, by = "gene")

# 3. 筛选持续变化基因
continuous_genes <- combined_logFC %>%
  rowwise() %>%
  mutate(
    all_positive = all(c_across(-gene) > 0), # 所有比较logFC>0
    all_negative = all(c_across(-gene) < 0)  # 所有比较logFC<0
  ) %>%
  filter(all_positive | all_negative) %>%
  ungroup()

# 4. 拆分上/下调基因
up_genes <- continuous_genes %>% filter(all_negative) %>% pull(gene)
down_genes <- continuous_genes %>% filter(all_positive) %>% pull(gene)

# 5. 保存结果
write.csv(data.frame(gene = up_genes), 
          "continuously_upregulated_genes.csv", row.names = FALSE)
write.csv(data.frame(gene = down_genes), 
          "continuously_downregulated_genes.csv", row.names = FALSE)

### 结果解读说明 ###
# 上調基因（up_genes）：在四个阶段过渡中持续上调（每个比较的logFC < 0）
# 即基因表达：infant < childhood < adolescence < adulthood
# 
# 下調基因（down_genes）：在四个阶段过渡中持续下调（每个比较的logFC > 0）
# 即基因表达：infant > childhood > adolescence > adulthood