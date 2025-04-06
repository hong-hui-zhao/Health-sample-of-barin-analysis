library(Seurat)
library(tidyverse)
library(patchwork)

# 定义阶段顺序
stage_order <- c("infant", "Middle childhood", "Adolescence",
                 "Early adulthood", "Late adulthood")

# 检查基因是否存在
m6a_genes <- c("METTL3", "METTL14", "METTL16", "ALKBH5", "FTO","NSUN2", "NSUN5",
               "MBD6", "MECP2",'IDH1', 'IDH2',"YTHDC1", "YTHDC2", "YTHDF1", 
               "YTHDF2", "YTHDF3","TET1", "TET2","NSUN7", "MBD5")
available_genes <- m6a_genes[m6a_genes %in% rownames(sce)]
missing_genes <- setdiff(m6a_genes, available_genes)
if(length(missing_genes) > 0) message("缺失基因：", paste(missing_genes, collapse = ", "))

# 提取数据并转换格式
exp_data <- FetchData(sce, vars = c(available_genes, "celltype", "life_stage")) %>% 
  pivot_longer(cols = all_of(available_genes),
               names_to = "gene",
               values_to = "expression") %>% 
  mutate(life_stage = factor(life_stage, levels = stage_order))

# 计算分组均值
summary_data <- exp_data %>%
  group_by(celltype, life_stage, gene) %>%
  summarise(mean_expression = mean(expression),
            .groups = "drop")

# 核心修改：嵌套式图形生成
plot_list <- summary_data %>%
  group_by(gene, celltype) %>%  # 按基因+细胞类型双重分组
  group_split() %>%            # 分割为独立数据集
  map(~ {
    ggplot(.x, aes(x = life_stage, y = mean_expression, fill = life_stage)) +
      geom_col(width = 0.6, show.legend = FALSE) +
      scale_fill_brewer(palette = "Pastel1") +
      labs(title = paste0(unique(.x$gene), " in ", unique(.x$celltype)),
           x = "Developmental Stage", y = "Expression Level") +
      theme_minimal(base_size = 12) +
      theme(
        axis.text.x = element_text(angle = 30, hjust = 1, color = "grey40"),
        plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
        panel.grid.major.x = element_blank(),
        panel.border = element_rect(fill = NA, color = "grey80", linewidth = 0.5)
      )
  })

# 智能分页保存（新增）
pdf("Gene_Celltype_Plots.pdf", width = 10, height = 8)
for(p in plot_list) {
  print(p + plot_annotation(tag_levels = "A"))  # 添加分页标签
}
dev.off()

write_csv(summary_data, "gene_expression_summary.csv")



library(Seurat)
library(ggplot2)
library(ggpubr)

# 检查基因列表中的基因是否存在于数据中
valid_genes <- m6a_genes[m6a_genes %in% rownames(sce)]

# 获取细胞类型信息
cell_types <- unique(sce$celltype) # 请确认meta.data中的列名是否为celltype

# 生成所有可能的细胞类型组合（两两比较）
comparisons <- combn(as.character(cell_types), 2, simplify = FALSE)

# 创建输出目录（如果不存在）
if (!dir.exists("vlnplots")) dir.create("vlnplots")

# 循环绘制每个基因
for (gene in valid_genes) {
  tryCatch({
    # 提取表达数据
    exp_data <- FetchData(sce, vars = c(gene, "celltype"))
    colnames(exp_data) <- c("expression", "celltype")
    
    # 创建基础小提琴图
    p <- ggplot(exp_data, aes(x = celltype, y = expression)) +
      geom_violin(aes(fill = celltype), scale = "width", trim = TRUE) +
      geom_boxplot(width = 0.1, fill = "white") +
      scale_fill_brewer(palette = "Set3") +
      theme_classic() +
      labs(title = gene, x = "Cell Type", y = "Expression Level") +
      theme(
        plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
        axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none"
      )
    
    # 保存图片
    ggsave(
      filename = paste0("vlnplots/", gene, "_vlnplot.png"),
      plot = p,
      width = max(6, length(cell_types)*1.5), # 根据细胞类型数量调整宽度
      height = 6,
      dpi = 300
    )
  }, error = function(e) {
    message(paste("Error processing", gene, ":", e$message))
  })
}



library(Seurat)
library(ggplot2)
library(ggpubr)

# 检查基因列表中的基因是否存在于数据中
valid_genes <- m6a_genes[m6a_genes %in% rownames(sce)]

# 生成所有可能的细胞类型组合（两两比较）
cell_types <- unique(sce$celltype)
comparisons <- combn(as.character(cell_types), 2, simplify = FALSE)

# 创建输出目录
if (!dir.exists("boxplots")) dir.create("boxplots")

# 循环绘制每个基因
for (gene in valid_genes) {
  tryCatch({
    # 提取表达数据
    exp_data <- FetchData(sce, vars = c(gene, "celltype"))
    colnames(exp_data) <- c("expression", "celltype")
    
    # 创建基础箱线图
    p <- ggplot(exp_data, aes(x = celltype, y = expression)) +
      geom_boxplot(
        aes(fill = celltype),
        width = 0.6,            # 调整箱型宽度
        outlier.size = 0.8,      # 异常点大小
        outlier.alpha = 0.5,     # 异常点透明度
        color = "black"          # 边框颜色
      ) +
      scale_fill_brewer(palette = "Set3") +
      theme_classic() +
      labs(
        title = paste(gene, "Expression"), 
        x = "Cell Type", 
        y = "Expression Level"
      ) +
      theme(
        plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
        axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
        axis.text.y = element_text(size = 10),
        legend.position = "none"
      )
    
    # 保存图片
    ggsave(
      filename = paste0("boxplots/", gene, "_boxplot.png"),
      plot = p,
      width = max(6, length(cell_types)*1.2), 
      height = 6,
      dpi = 300
    )
  }, error = function(e) {
    message(paste("Error processing", gene, ":", e$message))
  })
}



# 循环绘制每个基因
for (gene in valid_genes) {
  tryCatch({
    # 创建基础小提琴图
    p <- VlnPlot(sce,raster=FALSE, features = gene)
    # 保存图片
    ggsave(
      filename = paste0("vlnplots/", gene, "_vlnplot.png"),
      plot = p,
      width = max(6, length(cell_types)*1.5), # 根据细胞类型数量调整宽度
      height = 6,
      dpi = 300
    )
  }, error = function(e) {
    message(paste("Error processing", gene, ":", e$message))
  })
}
