library(Seurat)
library(dplyr)

# 定义m6A基因列表
m6a_genes <- c("METTL3", "METTL14", "METTL16", "ALKBH5", "FTO","NSUN2", "NSUN5",
               "MBD6", "MECP2",'IDH1', 'IDH2',"YTHDC1", "YTHDC2", "YTHDF1", 
               "YTHDF2", "YTHDF3","TET1", "TET2","NSUN7", "MBD5")

# 提取表达矩阵和元数据
expr_matrix <- GetAssayData(sce, assay = "SCT", slot = "data")[m6a_genes, ]
metadata <- sce@meta.data

# 转换年龄变量为数值型
metadata$dev_age <- as.numeric(gsub("[^0-9]", "", metadata$development_stage))
table(metadata$dev_age)

# 对life_stage_new进行数值化
stage_levels <- c("Early childhood", "Middle childhood", "Lately childhood",
                  "Early adulthood", "Middle adulthood", "Lately adulthood")
metadata$life_stage_num <- as.numeric(factor(metadata$life_stage_new, levels = stage_levels))

# 创建结果存储数据框
results <- data.frame(
  Gene = character(),
  CellType = character(),
  AgeType = character(),
  CorType = character(),
  Correlation = numeric(),
  Pvalue = numeric(),
  stringsAsFactors = FALSE
)

# 遍历每个细胞类型
celltypes <- unique(metadata$celltype)

for (ct in celltypes) {
  # 筛选当前细胞类型的细胞
  ct_cells <- rownames(metadata)[metadata$celltype == ct]
  
  # 提取对应的表达矩阵
  ct_expr <- expr_matrix[, ct_cells, drop = FALSE]
  
  # 提取元数据子集
  ct_meta <- metadata[ct_cells, c("life_stage_num", "dev_age")]
  
  # 遍历每个基因
  for (gene in m6a_genes) {
    # 获取当前基因的表达值
    gene_expr <- as.numeric(ct_expr[gene, ])
    
    # 计算相关性（Pearson 和 Spearman）
    calc_cor <- function(age_var, type = "pearson") {
      if (length(unique(age_var)) < 2) return(list(estimate = NA, p.value = NA))
      cor_test <- cor.test(gene_expr, age_var, method = type, exact = FALSE)
      list(estimate = cor_test$estimate, p.value = cor_test$p.value)
    }
    
    # 对life_stage_num计算相关性
    life_cor <- calc_cor(ct_meta$life_stage_num, "pearson")
    results <- rbind(results, data.frame(
      Gene = gene, CellType = ct, AgeType = "life_stage",
      CorType = "RCC", Correlation = life_cor$estimate, Pvalue = life_cor$p.value
    ))
    
    life_cor_spearman <- calc_cor(ct_meta$life_stage_num, "spearman")
    results <- rbind(results, data.frame(
      Gene = gene, CellType = ct, AgeType = "life_stage",
      CorType = "SCC", Correlation = life_cor_spearman$estimate, Pvalue = life_cor_spearman$p.value
    ))
    
    # 对dev_age计算相关性
    dev_cor <- calc_cor(ct_meta$dev_age, "pearson")
    results <- rbind(results, data.frame(
      Gene = gene, CellType = ct, AgeType = "development_age",
      CorType = "RCC", Correlation = dev_cor$estimate, Pvalue = dev_cor$p.value
    ))
    
    dev_cor_spearman <- calc_cor(ct_meta$dev_age, "spearman")
    results <- rbind(results, data.frame(
      Gene = gene, CellType = ct, AgeType = "development_age",
      CorType = "SCC", Correlation = dev_cor_spearman$estimate, Pvalue = dev_cor_spearman$p.value
    ))
  }
}

# 过滤无效结果
final_results <- results %>%
  filter(!is.na(Correlation)) %>%
  arrange(CellType, Gene, AgeType, CorType)

# 查看结果示例
head(final_results)

# 保存结果为CSV文件
write.csv(final_results, file = 'm6a_cor_stage.csv', row.names = FALSE)


##### --------------------------ALL

library(Seurat)
library(dplyr)

# 获取标准化后的表达矩阵
expr_matrix <- GetAssayData(sce, assay = "SCT", slot = "data")

# 创建空结果数据框
all_results <- data.frame(
  Gene = character(),
  CellType = character(),
  AgeType = character(),
  CorType = character(),
  Correlation = numeric(),
  Pvalue = numeric(),
  stringsAsFactors = FALSE
)

# 获取细胞类型
celltypes <- unique(metadata$celltype)

# 遍历每个细胞类型
for (ct in celltypes) {
  # 筛选当前细胞类型的细胞
  ct_cells <- rownames(metadata)[metadata$celltype == ct]
  
  # 提取对应表达矩阵
  ct_expr <- expr_matrix[, ct_cells, drop = FALSE]
  
  # 提取元数据子集
  ct_meta <- metadata[ct_cells, c("life_stage_num", "dev_age")]
  
  # 遍历每个基因
  for (gene in rownames(ct_expr)) {
    # 提取非零基因表达值
    gene_expr <- ct_expr[gene, ]
    non_zero_expr <- gene_expr[gene_expr != 0]  # 只考虑非零值
    
    # 如果该基因在该细胞类型中没有非零表达，跳过该基因
    if (length(non_zero_expr) < 2) {
      next
    }
    
    # 计算与life_stage_num和dev_age的相关性
    life_cor <- cor.test(non_zero_expr, ct_meta$life_stage_num, method = "pearson", exact = FALSE)
    dev_cor <- cor.test(non_zero_expr, ct_meta$dev_age, method = "pearson", exact = FALSE)
    
    # 将结果存入数据框
    all_results <- rbind(all_results, data.frame(
      Gene = gene, CellType = ct, AgeType = "life_stage",
      CorType = "RCC", Correlation = life_cor$estimate, Pvalue = life_cor$p.value
    ))
    
    all_results <- rbind(all_results, data.frame(
      Gene = gene, CellType = ct, AgeType = "development_age",
      CorType = "RCC", Correlation = dev_cor$estimate, Pvalue = dev_cor$p.value
    ))
  }
}

# 过滤无效结果
filtered_results <- all_results %>%
  filter(!is.na(Correlation)) %>%
  arrange(CellType, Gene, AgeType, CorType)

# 找到与发育阶段最相关的基因
top_genes <- filtered_results %>%
  filter(AgeType == "life_stage") %>%
  group_by(CellType) %>%
  top_n(10, abs(Correlation))  # 选择与life_stage相关性最强的前10个基因

# 查看最相关的基因
head(top_genes)
