# 加载包
devtools::install_github("BioSenior/ggVolcano")
library(ggplot2)
library(ggrepel)
library(readxl)
library(scales)
library(ggVolcano)
library(ggrepel)
library(ggfun)
library(grid)
library(tidyverse)
devtools::install_github("BioSenior/ggVolcano")
deg_data <- all_degs$`Middle childhood_vs_infant`
rownames(deg_data)<-deg_data$gene
res <- deg_data |> rownames_to_column("SYMBOL") |> arrange(desc(avg_log2FC))
df <- res |>  
  mutate(significant = case_when(avg_log2FC > 1 & p_val_adj < 0.05 ~ "Up",
                                 abs(avg_log2FC) < 1 | p_val_adj > 0.05 ~ "None",
                                 avg_log2FC < -0.5 & p_val_adj < 0.05 ~ "Down"))
df$significant |> as.factor()

head(df)
df_filtered <- df %>% 
  filter(!is.na(avg_log2FC), !is.na(p_val_adj), !is.infinite(avg_log2FC), !is.infinite(p_val_adj))

# Define the significance thresholds
log2FC_threshold <- 1
p_val_adj_threshold <- 0.05

# Add a column to categorize genes
df_filtered <- df_filtered %>%
  mutate(
    significant = case_when(
      abs(avg_log2FC) > log2FC_threshold & p_val_adj < p_val_adj_threshold ~ 
        ifelse(avg_log2FC > 0, "Up", "Down"),
      TRUE ~ "None"
    )
  )

# Select top 10 up-regulated and top 10 down-regulated genes
top_up <- df_filtered %>%
  tidyr::drop_na() %>%
  filter(significant == "Up") %>%
  arrange(desc(avg_log2FC), desc(-log10(p_val_adj))) %>%
  slice_head(n = 10)

top_down <- df_filtered %>%
  tidyr::drop_na() %>%
  filter(significant == "Down") %>%
  arrange(avg_log2FC, desc(-log10(p_val_adj))) %>%
  slice_head(n = 10)

top_genes <- bind_rows(top_up, top_down)


# 创建火山图
volcano_plot <- ggplot(df_filtered, aes(x = avg_log2FC, y = -log2(p_val_adj))) +
  geom_point(aes(color = avg_log2FC, size = -log10(p_val_adj))) +
  scale_color_gradientn(colours = c("#3288bd", "#66c2a5","#ffffbf", "#f46d43", "#9e0142"),
                        values = seq(0, 1, 0.2)) +
  scale_fill_gradientn(colours = c("#3288bd", "#66c2a5","#ffffbf", "#f46d43", "#9e0142"),
                       values = seq(0, 1, 0.2)) +
  geom_text_repel(
    data = top_genes,
    aes(label = SYMBOL),
    box.padding = 0.5,
    point.padding = 0.3,
    segment.color = "grey50",
    segment.size = 0.2,
    arrow = arrow(length = unit(0.01, "npc")),
    force = 2,
    max.overlaps = Inf,
    size = 4  # 增加基因标签的字体大小
  ) +
  geom_vline(xintercept = c(-log2FC_threshold, log2FC_threshold), linetype = 3) +
  geom_hline(yintercept = -log10(p_val_adj_threshold), linetype = 4) +
  scale_size(range = c(1, 7)) +
  labs(
    x = "log2(Fold Change)",
    y = "-log10(Adjusted p-value)"
  ) +
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    legend.background = element_rect(color = "#808080", linetype = 1),
    axis.text = element_text(size = 6, color = "#000000"),  # 增加轴刻度标签的字体大小
    axis.title = element_text(size = 8),  # 增加轴标题的字体大小
    legend.text = element_text(size = 4),  # 增加图例文本的字体大小
    legend.title = element_text(size = 6)  # 增加图例标题的字体大小
  ) +
  annotate(geom = "text", x = max(df_filtered$avg_log2FC, na.rm = TRUE) * 0.9, 
           y = 1, label = "p = 0.05", size = 8) +  # 增加p值标注的字体大小
  coord_cartesian(clip = "off") +
  annotation_custom(
    grob = grid::segmentsGrob(
      x0 = unit(0.1, "npc"), x1 = unit(0.3, "npc"),
      y0 = unit(0.95, "npc"), y1 = unit(0.95, "npc"),
      arrow = arrow(angle = 30, length = unit(0.1, "inches"), ends = "first", type = "closed"),
      gp = grid::gpar(lwd = 1.5, col = "#74add1")
    )
  ) +
  annotation_custom(
    grob = grid::textGrob(
      label = "Down",
      x = unit(0.2, "npc"), y = unit(0.98, "npc"),
      gp = grid::gpar(col = "#74add1", fontsize = 6)  # 增加"Down"标签的字体大小
    )
  ) +
  annotation_custom(
    grob = grid::segmentsGrob(
      x0 = unit(0.7, "npc"), x1 = unit(0.9, "npc"),
      y0 = unit(0.95, "npc"), y1 = unit(0.95, "npc"),
      arrow = arrow(angle = 30, length = unit(0.1, "inches"), ends = "last", type = "closed"),
      gp = grid::gpar(lwd = 1.5, col = "#d73027")
    )
  ) +
  annotation_custom(
    grob = grid::textGrob(
      label = "Up",
      x = unit(0.8, "npc"), y = unit(0.98, "npc"),
      gp = grid::gpar(col = "#d73027", fontsize = 6)  # 增加"Up"标签的字体大小
    )
  )

# 打印图表
print(volcano_plot)

# 保存图表
ggsave("volcano_plot.png", volcano_plot, width = 10, height = 8, dpi = 300)
ggsave("volcano_plot.pdf", volcano_plot, width = 10, height = 8, dpi = 300)
