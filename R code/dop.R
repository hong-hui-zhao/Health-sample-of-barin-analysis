###
DefaultAssay(sce) <- "SCT"
markers <- FindAllMarkers(sce, logfc.threshold = 0.25, min.pct = 0.25, only.pos = T,test.use = 'MAST')
top5 <- markers %>% group_by(cluster) %>% top_n(5, avg_log2FC)

BiocManager::install("muscat")

library(Seurat)
library(muscat)
# 载入所需包
library(Seurat)
library(ggplot2)

# 提取基因表达数据（使用 SCT assay 的 data slot）
gene_expression <- FetchData(sce, vars = c('ENSG00000286749.1', 'celltype'), slot = "data")
colnames(gene_expression) <- c("expression", "celltype")  # 重命名列
gene_expression <- gene_expression[gene_expression$expression >0,]
# 计算每个 celltype 的中位表达值
mean_expression <- gene_expression %>%
  group_by(celltype) %>%
  summarise(mean = mean(expression))
table(gene_expression$celltype)


# 绘制基础箱线图
p <- ggplot(gene_expression, aes(x = celltype, y = expression, fill = celltype)) +
  geom_boxplot(width = 0.6, outlier.size = 0.5) +  # 箱线图
  geom_text(
    data = mean_expression,
    aes(x = celltype, y = mean, label = round(mean, 2)),  # 标注中位值
    vjust = -0.5, size = 3, color = "black"
  ) +
  labs(
    title = "Expression of ENSG00000286749.1 (Median Values)",
    x = "celltype",
    y = "SCT Corrected Expression (log1p)"
  ) +
  theme_classic() +
  theme(legend.position = "none")  # 隐藏图例

# 显示图形
print(p)
VlnPlot(sce,feature = 'ENSG00000286749.1',raster=FALSE,slot = "data")

library(ggplot2)
p <- DotPlot(sce, features = top5$gene
             , group.by = "celltype")+coord_flip()
exp <- p$data
library(forcats)
exp$features.plot <- as.factor(exp$features.plot)
exp$features.plot <- fct_inorder(exp$features.plot)

exp$id <- as.factor(exp$id)
exp$id <- fct_inorder(exp$id)


ggplot(exp,aes(x=id,y= features.plot))+
  geom_point(aes(size=`pct.exp`,
                 color=`avg.exp.scaled`))+
  geom_point(aes(size=`pct.exp`,color=`avg.exp.scaled`),
             shape=21,color="black",stroke=1)+
  theme(panel.background =element_blank(),
        axis.line=element_line(colour="black", size = 1),
        panel.grid = element_blank(),
        axis.text.x=element_text(size=11,color="black",angle=90), 
        axis.text.y=element_text(size=11,color="black"))+
  scale_color_gradientn(colors = colorRampPalette(c("white", "#00C1D4", "#FFED99","#FF7600"))(10))+
  labs(x=NULL,y=NULL)


#c("#1E3163", "#00C1D4", "#FFED99","#FF7600")


p1 <- ggplot(exp,aes(x=id,y= features.plot))+
  geom_point(aes(size=`pct.exp`,
                 color=`avg.exp.scaled`))+
  geom_point(aes(size=`pct.exp`,color=`avg.exp.scaled`),
             shape=21,color="black",stroke=1)+
  theme(panel.background =element_blank(),
        axis.line=element_line(colour="black"),
        panel.grid = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.x=element_text(size=10,color="black"),
        axis.text.y=element_text(size=12,color="black"))+
  scale_color_gradientn(colors = colorRampPalette(c("white", "#00C1D4", "#FFED99","#FF7600"))(10))+
  labs(x=NULL,y=NULL)
ggsave('dop_marker.pdf',p1,width = 20,height = 20)






p2 <- p1+  
  geom_rect(aes(xmin = 1 - 0.6, xmax = 6 + 0.5, ymin = 1 - 0.6, ymax = 5 + 0.5),
            fill = "transparent", color = "black", size = 0.5)+
  geom_rect(aes(xmin = 7 - 0.5, xmax = 12 + 0.5, ymin = 6 - 0.5, ymax =10 + 0.5),
            fill = "transparent", color = "black", size = 0.5)+
  geom_rect(aes(xmin = 13 - 0.5, xmax = 18 + 0.5, ymin = 11 - 0.5, ymax =15 + 0.5),
            fill = "transparent", color = "black", size = 0.5)+
  geom_rect(aes(xmin = 19 - 0.5, xmax = 24 + 0.5, ymin = 16 - 0.5, ymax =20 + 0.5),
            fill = "transparent", color = "black", size = 0.5)+
  scale_x_discrete(labels = c('','','',"Macrophage",'','',
                              '','','',"T cell",'','',
                              '','','','mDC','','',
                              '','','','Neutrophil','',''))+
  theme(axis.text.x = element_text(angle=45,
                                   vjust=1, 
                                   size=11,
                                   hjust=1,
                                   color = 'black')) 







my36colors <-c('#E5D2DD', '#53A85F', '#F1BB72', '#F3B1A0', '#D6E7A3', '#57C3F3', '#476D87',
               '#E95C59', '#E59CC4', '#AB3282', '#23452F', '#BD956A', '#8C549C', '#585658',
               '#9FA3A8', '#E0D4CA', '#5F3D69', '#C5DEBA', '#58A4C3', '#E4C755', '#F7F398',
               '#AA9A59', '#E63863', '#E39A35', '#C1E6F3', '#6778AE', '#91D0BE', '#B53E2B',
               '#712820', '#DCC1DD', '#CCE0F5',  '#CCC9E6', '#625D9E', '#68A180', '#3A6963',
               '#968175')


library(dplyr)
library(ggplot2)

group <- colnames(exp) %>% as.data.frame() %>% 
  mutate(group=c(rep("BM1",1),rep("BM2",1),rep("BM3",1),rep("GM1",1),rep("GM2",1),rep('GM3',1))) %>%
  mutate(p="group") %>%
  ggplot(aes(.,y=p,fill=group))+
  geom_tile() + 
  theme_minimal()+xlab(NULL) + ylab(NULL) +
  theme(axis.text.x = element_blank(),
        axis.text.y=element_blank())+
  theme(legend.position = "none")



p2$layout$clip[grep("panel-2-\\d+", p2$layout$name)] <- "off" 
top <- ggplotGrob(group)
library(aplot)#拼图
p2 +annotation_custom(top,xmin=0,xmax=7,ymin=20,ymax = 22)+
  annotation_custom(top,xmin=6.1,xmax=13,ymin=20,ymax = 22)+
  annotation_custom(top,xmin=12.1,xmax=19,ymin=20,ymax = 22)+
  annotation_custom(top,xmin=18.1,xmax=25,ymin=20,ymax = 22)+
  theme(legend.direction = 'horizontal',
        legend.position = 'bottom',
        legend.justification=c(1,2),
        legend.key.width=unit(0.5,"cm"),
        legend.key.height = unit(0.3, "cm"),
        legend.title.position = 'top')
