lapply(c("dplyr","Seurat","HGNChelper","openxlsx"), library, character.only = T)
setwd('/sibcb1/douxiaoyanglab1/zhaohonghui/Glioma/sc-type-master/sc-type-master')
data <- readRDS('/sibcb1/douxiaoyanglab1/zhaohonghui/Glioma/health_barin/harmony_integrated_SCT_snRNA.rds')
source("/sibcb1/douxiaoyanglab1/zhaohonghui/Glioma/sc-type-master/sc-type-master/R/gene_sets_prepare.R")
source("/sibcb1/douxiaoyanglab1/zhaohonghui/Glioma/sc-type-master/sc-type-master/R/sctype_score_.R")


db = "/sibcb1/douxiaoyanglab1/zhaohonghui/Glioma/sc-type-master/sc-type-master/ScTypeDB_full.xlsx"
tissue = "Brain"

# prepare gene sets
gs_list = gene_sets_prepare(db, tissue)
data <- SetIdent(data,value = "SCT_snn_res.0.3")

expr_mat <- as.matrix(GetAssayData(data, assay = "SCT", slot = "data"))
data <- SetIdent(data,value = "SCT_snn_res.0.3")

expr <- as.matrix(data@assays[["RNA"]]@data)
es.max = sctype_score(scRNAseqData = expr, scaled = F, gs = gs_list$gs_positive,
                      gs2 = gs_list$gs_negative)

save.image('/sibcb1/douxiaoyanglab1/zhaohonghui/Glioma/health_barin/sctype.RData')


cL_resutls = do.call("rbind", lapply(unique(data@meta.data$SCT_snn_res.0.3), function(cl){
  es.max.cl = sort(rowSums(es.max[ ,rownames(data@meta.data[data@meta.data$SCT_snn_res.0.3==cl, ])]), decreasing = !0)
  head(data.frame(cluster = cl, type = names(es.max.cl), scores = es.max.cl, ncells = sum(data@meta.data$SCT_snn_res.0.3==cl)), 10)
}))

sctype_scores = cL_resutls %>% group_by(cluster) %>% top_n(n = 1, wt = scores)  
write.csv(sctype_scores, "sctype_scores.csv")
