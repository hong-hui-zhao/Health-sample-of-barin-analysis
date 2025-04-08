library(tidyverse)
# 生成计数矩阵并转换为按列（患者）的比例
cell_prop_df <- as.data.frame.table(
  prop.table(table(sce@meta.data$celltype, sce@meta.data$life_stage), margin = 2) * 100,
  responseName = "Proportion"
) 

colnames(cell_prop_df)[1] <- 'celltype'
colnames(cell_prop_df)[2] <- 'life_stage'

stage_order <- c("infant", "Middle childhood", "Adolescence",
                 "Early adulthood", "Late adulthood")

cell_prop_df$life_stage <- factor(
  cell_prop_df$life_stage, 
  levels = stage_order  # 核心设置顺序
)

ggplot(cell_prop_df, aes(x = life_stage, y = Proportion, fill = celltype)) +
  geom_col(color = "black", width = 0.8) +
  scale_x_discrete(limits = stage_order) + 
  scale_fill_brewer(palette = "Set3") +
  labs(
    title = "Cell Type Composition by life_stage",
    x = "",  # 这里设置为空字符串即可移除X轴标签
    y = "Percentage (%)"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),  # 标题居中加粗
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
    legend.position = "right",
    legend.title = element_blank(),
    panel.grid.major.x = element_blank()
  )
ggsave('cellprop_life_stage_bar.pdf',width = 8,height = 5)
write.csv(cell_prop_df,file = 'cell_prop_df.csv')


