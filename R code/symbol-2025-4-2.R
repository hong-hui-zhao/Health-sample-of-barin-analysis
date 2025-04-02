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
options(bitmapType='cairo')

### setwd path
setwd('/Glioma/health_barin')

rawsce <- readRDS('/Glioma/health_barin/intergated_snRNA.rds')
rawsce <- data
### Gene names and developmental names are modified
# 1. 检查列类型（确认是否因子）
class(rawsce@meta.data$development_stage)  # 应返回 "factor"

## 检查名称
table(rawsce@meta.data$development_stage)

# 2. 转换为字符型，因为之前遇见了“因子陷阱”（当尝试给因子变量赋一个不在其 Levels 中的值 时，R 会强制将其转换为 NA，并发出警告）
rawsce@meta.data$development_stage <- as.character(rawsce@meta.data$development_stage)

# 3. 执行替换操作
rawsce@meta.data$development_stage[rawsce@meta.data$development_stage == "infant stage"] <- "0-year-old stage"

table(rawsce@meta.data$development_stage)

rawsce@meta.data$development_stage <- factor(rawsce@meta.data$development_stage)

# 5. 验证结果
table(rawsce@meta.data$development_stage)

### 基因名称修改，rds中，只给出了 ensembl 名称 ，因此需要将其转化为symbol名称
### 原文中在 rawsce@assays[["RNA"]]@meta.features[["feature_name"]] 给出了部分symbol，我们根据这个进行symbol name 修改
ens <- as.matrix(rawsce@assays[["RNA"]]@data@Dimnames[[1]])

GRCh38_113 <- read.csv("~/Glioma/GRCh38.113.csv",header = TRUE,row.names =1)
GRCh38_113 <- na.omit(GRCh38_113)
# 创建映射关系
id_mapping <- setNames(GRCh38_113$gene_symbol, GRCh38_113$ensembl_gene_id)

# 替换Ensembl ID为基因Symbol，如果找不到对应关系，则保留原来的Ensembl ID
new_gene_names <- as.matrix(ifelse(ens %in% names(id_mapping), id_mapping[ens], ens))

ensembl_ids <- data.frame(ensembl_gene_id = rawsce@assays[["RNA"]]@data@Dimnames[[1]])
library(dplyr)
# 进行左连接，匹配基因Symbol
ensembl_ids <- ensembl_ids %>%
  left_join(GRCh38_113, by = "ensembl_gene_id") %>%
  mutate(gene_symbol = ifelse(is.na(gene_symbol), ensembl_gene_id, gene_symbol)) 

### 这里发现我们的结果与提取rds中的rawsce@assays[["RNA"]]@meta.features[["feature_name"]]的结果是一致的
symbol<- as.matrix(rawsce@assays[["RNA"]]@meta.features[["feature_name"]])
rawsce@assays[["RNA"]]@data@Dimnames[[1]] <- symbol
rawsce@assays[["RNA"]]@counts@Dimnames[[1]] <- symbol
new_genes <- gsub("_", "-", rownames(rawsce))
new_genes <- make.unique(new_genes)

### 手动计算，线粒体跟核糖体的比例
DefaultAssay(rawsce) <- 'RNA'
grep("^MT-", rownames(rawsce), value = TRUE)
mt_genes <- grep("^MT-", rownames(rawsce), value = TRUE)
rb_genes <- grep("^RPS|^RPL", rownames(rawsce), value = TRUE) 
rawsce[["percent.mt"]] <- colSums(rawsce[mt_genes, ]) / colSums(rawsce) * 100
rawsce[["percent.rb"]] <- colSums(rawsce[rb_genes, ]) / colSums(rawsce) * 100
# 创建新列 donor_tissue 
rawsce@meta.data$donor_tissue <- paste(
  rawsce@meta.data$donor_id, 
  rawsce@meta.data$tissue, 
  sep = "_"  # 用下划线连接两部分
)
# 检查前6行验证
head(rawsce@meta.data$donor_tissue)
#查看各组合的计数是否与原表一致
table(rawsce@meta.data$donor_tissue)
### 保存原始分析数据
saveRDS(rawsce,file = '/Glioma/health_barin/newsymbols.rds')
