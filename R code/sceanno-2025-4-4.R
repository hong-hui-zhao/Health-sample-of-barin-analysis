##########################################

### health Brain sample analysis

### Begin 2025-4-4

### PMID: 39962241
        
        
        
        

### Author：honghui Zhao

##########################################

### load packages for analysis
library(SingleR)
library(Seurat)
library(scCATCH)
setwd('/health_barin')

#### load data
sce <- readRDS('health_barin/intergated_snRNA.rds')
sce <- SetIdent(sce,value = "SCT_snn_res.0.4")

### load ref data

### singleR

BlueprintEncode <- readRDS('/Glioma/SingleR/BlueprintEncodeData.rds')
DatabaseImmune <- readRDS('/SingleR/DatabaseImmuneCellExpressionData.rds')
HPCA <- readRDS('/SingleR/HumanPrimaryCellAtlasData.rds')
Monaco <-readRDS('/SingleR/MonacoImmuneData.rds')
Novershtern <-readRDS('/SingleR/NovershternHematopoieticData.rds')


# 1. 加载所有人类参考数据集 ----------------------------------------------------
ref_list <- list(
  BlueprintEncode = BlueprintEncode,
  DatabaseImmune = DatabaseImmune,
  HPCA = HPCA,
  Monaco = Monaco,
  Novershtern = Novershtern
)
 
# 2. 数据预处理 --------------------------------------------------------------
# 假设您的Seurat对象名为scRNA
sce <- SetIdent(sce,value = "SCT_snn_res.0.4")
ens <- as.matrix(sce@assays[["SCT"]]@data@Dimnames[[1]])

GRCh38_113 <- read.csv("/Glioma/GRCh38.113.csv",header = TRUE,row.names =1)
GRCh38_113 <- na.omit(GRCh38_113)
# 创建映射关系
id_mapping <- setNames(GRCh38_113$gene_symbol, GRCh38_113$ensembl_gene_id)

# 替换Ensembl ID为基因Symbol，如果找不到对应关系，则保留原来的Ensembl ID
new_gene_names <- as.matrix(ifelse(ens %in% names(id_mapping), id_mapping[ens], ens))

ensembl_ids <- data.frame(ensembl_gene_id = sce@assays[["SCT"]]@data@Dimnames[[1]])
library(dplyr)
# 进行左连接，匹配基因Symbol
ensembl_ids <- ensembl_ids %>%
  left_join(GRCh38_113, by = "ensembl_gene_id") %>%
  mutate(gene_symbol = ifelse(is.na(gene_symbol), ensembl_gene_id, gene_symbol)) 

sce@assays[["SCT"]]@data@Dimnames[[1]] <- ensembl_ids$gene_symbol
sce@assays[["SCT"]]@counts@Dimnames[[1]] <- ensembl_ids$gene_symbol

testdata <- GetAssayData(sce, assay = 'SCT',slot = "count")  # 提取标准化表达矩阵
# 3. 执行多参考联合注释 -------------------------------------------------------
anno <- SingleR(
  test = testdata,
  ref = ref_list,  # 所有参考数据集组成的列表
  labels = list(   # 各数据集对应的标签列（与ref_list顺序一致）
    ref_list$BlueprintEncode$label.main,
    ref_list$DatabaseImmune$label.main,
    ref_list$HPCA$label.main,
    ref_list$Monaco$label.main, 
    ref_list$Novershtern$label.main
  ),
  clusters = sce@active.ident  # 按cluster注释（推荐）
)

# 4. 结果整合与可视化 ---------------------------------------------------------
# 创建注释对照表
celltype_df <- data.frame(
  ClusterID = rownames(anno),
  CombinedLabel = anno$pruned.labels,
  stringsAsFactors = FALSE
)
# 5. 注释结果导出 ------------------------------------------------------------
write.csv(celltype_df,file = '/Glioma/health_barin/singleR_anno.csv')

### scCATCH
# 创建scCATCH对象
obj <- createscCATCH(data = sce@assays$SCT$counts, cluster =  as.character(sce$SCT_snn_res.0.4))

# 查找每个簇的标记基因
obj <- findmarkergene(obj, species = "Human", marker = cellmatch,tissue = "Brain")
# 查找每个簇的细胞类型
obj <- findcelltype(obj)
celltype_scCATCH <- as.data.frame(obj@celltype)
celltype_scCATCH <- as.data.frame(obj@celltype$cell_type,obj@celltype$cluster)
write.csv(celltype_scCATCH,file = '/Glioma/health_barin/celltype_scCATCH.csv')


#### SCtype
lapply(c("dplyr","Seurat","HGNChelper","openxlsx"), library, character.only = T)
source("/sc-type-master/sc-type-master/R/gene_sets_prepare.R")
source("/sc-type-master/sc-type-master/R/sctype_score_.R")


db = "/sc-type-master/sc-type-master/ScTypeDB_full.xlsx"
tissue = "Brain"

# prepare gene sets
gs_list = gene_sets_prepare(db, tissue)

expr_mat <- as.matrix(GetAssayData(sce, assay = "SCT", slot = "data"))

expr <- as.matrix(sce@assays[["SCT"]]@data)
es.max = sctype_score(scRNAseqData = expr, scaled = F, gs = gs_list$gs_positive,
                      gs2 = gs_list$gs_negative)

cL_resutls = do.call("rbind", lapply(unique(sce@meta.data$SCT_snn_res.0.4), function(cl){
  es.max.cl = sort(rowSums(es.max[ ,rownames(sce@meta.data[sce@meta.data$SCT_snn_res.0.4==cl, ])]), decreasing = !0)
  head(data.frame(cluster = cl, type = names(es.max.cl), scores = es.max.cl, ncells = sum(sce@meta.data$SCT_snn_res.0.4==cl)), 10)
}))

sctype_scores = cL_resutls %>% group_by(cluster) %>% top_n(n = 1, wt = scores)  
write.csv(sctype_scores,file="/Glioma/health_barin/sctype_scores.csv")

saveRDS(sce,file = '/Glioma/health_barin/newsymbols.rds')


library(SCINA)
options(bitmapType='cairo')
# 从Seurat对象中获取标准化后的表达矩阵（确保使用正确的assay）
expression_matrix <- GetAssayData(sce, assay = "SCT", slot = "data") # 使用log-normalized数据

# 转置矩阵以满足SCINA的输入要求（基因为行，细胞为列）
expression_matrix <- as.matrix(expression_matrix)

# 运行SCINA
scina_results <- SCINA(
  exp = expression_matrix,
  signatures = celltype_markers,
  max_iter = 100,             # 最大迭代次数
  convergence_n = 10,         # 收敛判断窗口
  convergence_rate = 0.999,   # 收敛阈值
  sensitivity_cutoff = 0.9,   # 灵敏度阈值
  rm_overlap = TRUE,          # 移除重叠基因
  allow_unknown = TRUE,       # 允许未知细胞类型
  log_file = "SCINA.log"       # 日志文件
)


table(scina_results$cell_labels)

# 将结果添加到Seurat对象
sce$SCINA_celltype <- scina_results$cell_labels
table(sce$SCINA_celltype,sce$SCT_snn_res.0.4)
# 可视化结果
DimPlot(data, group.by = "SCINA_celltype",raster=FALSE)



# 按细胞类型分组标记基因（已调整部分可能拼写问题，如NEST改为NES）
celltype_markers <- list(
  # 兴奋性神经元（补充谷氨酸能神经元相关基因）
  Excitatory_neurons = c("SLC17A7", "SLC17A6", "SYNPR", "GRIN1", "CAMK2A", "GRIN2B", "GLS", "GRIN2A"),  
  # 抑制性神经元（补充GABA能神经元相关基因）
  Inhibitory_neurons = c("GAD1", "GAD2", "SLC32A1", "PVALB", "SST", "VIP","GABBR1", "GABBR2"),       
  
  # 星形胶质细胞（补充免疫调节相关基因）
  Astrocytes = c("GFAP", "AQP4", "SLC1A3", "ALDH1L1", "GLUL","S100B", "CD40", "CD80", "CD86", "C5AR1", "SLC1A2"),  
  
  # 少突胶质细胞（补充发育阶段特异性基因）
  Oligodendrocytes = c("MBP", "PLP1", "MOBP", "OLIG2", "CNP","OLIG1", "CLDN11", "MOG", "MAG", "GALC", "UGT8"), 
  
  # 小胶质细胞（补充激活状态标志物）
  Microglia = c("CX3CR1", "P2RY12", "TMEM119", "AIF1", "CD68","ITGAM", "C1QA", "ADGRE1", "TNF", "CCL4"), 
  
  # 内皮细胞（补充血管生成相关基因）
  Endothelial = c("CLDN5", "FLT1", "PECAM1", "CD34", "ICAM1","VWF", "A2M", "APOLD1", "ENG", "CDH5"),   
  
  # 少突胶质前体细胞（补充表面受体）
  OPC = c("PDGFRA", "CSPG4", "SOX10","LHFPL3", "MEGF11", "PCDH15"),                  
  
  # 神经前体细胞（补充增殖标志物）
  NPC = c("NES", "MKI67", "EOMES", "ASCL1","SOX2", "PAX6", "OTX2", "SMARCA4"),            
  
  # 新增特殊神经元类型（来自单细胞数据库）
  Cholinergic_neurons = c("CHAT", "SLC18A3", "ACHE"),          
  Dopaminergic_neurons = c("TH", "SLC6A3", "FOXA2", "KCNJ6", "NR4A2"),  
  Serotonergic_neurons = c("TPH1", "SLC6A4", "FEV", "HTR1A"),  
  
  # 新增胶质细胞亚型（参考最新脑图谱）
  Myelinating_Schwann = c("SOX10", "EGR2", "MBP", "MPZ"),      
  Non_myelinating_Schwann = c("SOX10", "GAP43", "NCAM1"),     
  Radial_glial = c("VIM", "NES", "PAX6", "HES1", "FABP7"),     
  
  # 新增肿瘤相关细胞（整合癌症数据库）
  Cancer_cells = c("CD44", "MBTPS2", "PARP1"),                
  Cancer_stem_cells = c("FUT4", "MSI1", "NES", "PROM1")        
)
