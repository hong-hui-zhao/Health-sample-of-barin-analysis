library(Seurat)
library(tidyverse)
library(harmony)
library(sctransform)
library(Cairo)

sce <- readRDS('/Glioma/health_barin/Panos.rds')

DefaultAssay(sce) <- 'RNA'
grep("^MT-", rownames(sce), value = TRUE)
mt_genes <- grep("^MT-", rownames(sce), value = TRUE)
rb_genes <- grep("^RPS|^RPL", rownames(sce), value = TRUE) 
sce[["percent.mt"]] <- colSums(sce[mt_genes, ]) / colSums(sce) * 100
sce[["percent.rb"]] <- colSums(sce[rb_genes, ]) / colSums(sce) * 100

# 检查是否有重复基因名
sce <- NormalizeData(sce, normalization.method = "LogNormalize") 
sce <- FindVariableFeatures(sce, selection.method = "vst", nfeatures = 3000)
sce <- ScaleData(sce, vars.to.regress = c("percent.mt", "percent.rb"))

sce <- SCTransform(
  sce,
  vars.to.regress = c("percent.mt", "percent.rb"),
  assay = "RNA",
  new.assay.name = "SCT"
)


# 传统流程分析
sce <- RunPCA(sce, assay = "RNA", npcs = 50, verbose = FALSE)
sce <- RunHarmony(sce, group.by.vars = "donor_tissue", assay.use = "RNA", reduction = "pca", dims.use = 1:20)
sce <- RunUMAP(sce, reduction = "harmony", dims = 1:20, reduction.name = "umap_rna")
sce <- FindNeighbors(sce, reduction = "harmony", dims = 1:20, assay = "RNA")
sce <- FindClusters(sce, resolution = seq(0.1, 1.0, by = 0.1), graph.name = "RNA_snn")

# SCTransform 流程分析
sce <- RunPCA(sce, assay = "SCT", npcs = 50, verbose = FALSE, reduction.name = "pca_sct")
sce <- RunHarmony(sce, group.by.vars = "donor_tissue", assay.use = "SCT", reduction = "pca_sct", dims.use = 1:15, reduction.name = "harmony_sct")
sce <- RunUMAP(sce, reduction = "harmony_sct", dims = 1:15, reduction.name = "umap_sct")
sce <- FindNeighbors(sce, reduction = "harmony_sct", dims = 1:15, assay = "SCT")
sce <- FindClusters(sce, resolution = seq(0.1, 1.0, by = 0.1), graph.name = "SCT_snn")


saveRDS(sce, file = "/intergated_snRNA.rds")


