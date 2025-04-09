library(Seurat)
library(dplyr)
library(ComplexHeatmap)
library(circlize)
library(RColorBrewer)

# 定义m6A相关基因
m6a_genes <- c("METTL3", "METTL14", "METTL16", "ALKBH5", "FTO","NSUN2", "NSUN5",
               "MBD6", "MECP2",'IDH1', 'IDH2',"YTHDC1", "YTHDC2", "YTHDF1", 
               "YTHDF2", "YTHDF3","TET1", "TET2","NSUN7", "MBD5")

# 设置默认assay
DefaultAssay(sce) <- "SCT"

# 创建组合标签（使用连字符连接）
sce$celltype_stage <- paste(sce$celltype, sce$life_stage_new, sep = "-")

# 获取有序的细胞类型和发育阶段信息
celltypes <- c("Astrocyte", "GABAergic neurons", "Microglia cell", 
               "Glutamatergic neurons", "Bergmann glial cell", "Lepotomeningeal cell",
               "Smooth muscle cell", "Endothelial cell", "Oligodendrocyte precursor cell",
               "Oligodendrocyte", "Ependymal cell")
stages <- c("Early childhood", "Middle childhood", "Lately childhood",
            "Early adulthood", "Middle adulthood", "Lately adulthood")

# 计算平均表达量
avg_expr <- AverageExpression(
  sce,
  group.by = "celltype_stage",
  features = m6a_genes,
  assays = "SCT"
)$SCT

# 生成正确的列顺序（celltypes × stages）
desired_order <- as.vector(outer(celltypes, stages, paste, sep = "-"))
valid_columns <- desired_order[desired_order %in% colnames(avg_expr)]
avg_expr <- avg_expr[, valid_columns, drop = FALSE]

# 数据标准化（行Z-Score）
mat_zscore <- t(scale(t(as.matrix(avg_expr))))

# 构建注释信息
celltype_anno <- sapply(strsplit(colnames(avg_expr), "-"), `[`, 1)
stage_anno <- sapply(strsplit(colnames(avg_expr), "-"), `[`, 2)

# 设置颜色映射
celltype_colors <- c(
  "Astrocyte" = "#FF0000", "GABAergic neurons" = "#FF8B00",
  "Microglia cell" = "#E8FF00", "Glutamatergic neurons" = "#5DFF00",
  "Bergmann glial cell" = "#00FF2E", "Lepotomeningeal cell" = "#00FFB9",
  "Smooth muscle cell" = "#00B9FF", "Endothelial cell" = "#002EFF",
  "Oligodendrocyte precursor cell" = "#5D00FF", "Oligodendrocyte" = "#E800FF",
  "Ependymal cell" = "#FF008B"
)

stage_colors <- colorRampPalette(brewer.pal(8, "Set2"))(length(stages))
names(stage_colors) <- stages

# 创建顶部注释
top_anno <- HeatmapAnnotation(
  CellType = factor(celltype_anno, levels = celltypes),
  Stage = factor(stage_anno, levels = stages),
  col = list(
    CellType = celltype_colors,
    Stage = stage_colors
  ),
  annotation_name_side = "left",
  gap = unit(1, "mm"),
  show_legend = c(TRUE, TRUE)
)

# 绘制热图
pdf('m6a_gene_heatmap.pdf',width = 12,height = 6)
Heatmap(
  mat_zscore,
  name = "Z-Score\nExpression",
  top_annotation = top_anno,
  cluster_columns = FALSE,
  cluster_rows = TRUE,
  show_column_names = FALSE,
  show_row_names = TRUE,
  column_split = factor(celltype_anno, levels = celltypes),
  column_title = NULL,
  row_names_gp = gpar(fontsize = 12),
  column_names_rot = 45,
  heatmap_legend_param = list(
    title_position = "leftcenter-rot",
    legend_height = unit(4, "cm"))
)
dev.off()
