# 加载所需包 ----------------------------------------------------------------
library(Seurat)       # 单细胞分析核心包
library(muscat)       # 多样本单细胞分析工具，用于伪批量差异分析
library(SingleCellExperiment) # 单细胞数据容器

# 数据预处理 ----------------------------------------------------------------
# 假设原始数据存储在sce对象中，筛选生命阶段为Young和Adult的细胞
sce <- data
sce <- sce[, sce$life_group %in% c("Adult","Older")] 

# 转换为SingleCellExperiment对象（muscat要求的数据格式）
sce_sce <- as.SingleCellExperiment(sce)

# 使用muscat的prepSCE函数标准化元数据格式 -------------------------------
# 参数说明：
# - kid: 定义细胞类型列的元数据名称（对应cluster_id）
# - gid: 定义分组列的元数据名称（对应group_id，此处为Young/Adult）
# - sid: 定义样本来源列的元数据名称（对应sample_id）
# - drop: 删除未使用的因子级别
sce_sce <- prepSCE(sce_sce,
                   kid = "celltype", 
                   gid = "life_group",
                   sid = "donor_id",
                   drop = TRUE)

# 确保组别顺序正确（Young作为对照组，Adult作为实验组） ------------------
sce_sce$group_id <- factor(sce_sce$group_id, levels = c("Adult","Older"))

# 检查分组样本量 ----------------------------------------------------------
table(sce_sce$group_id)  # 查看Young和Adult组的细胞数量分布

# 获取细胞类型(kid)和样本(sid)信息 -------------------------------------
nk <- length(kids <- levels(sce_sce$cluster_id))  # 提取所有细胞类型名称
ns <- length(sids <- levels(sce_sce$sample_id))   # 提取所有样本ID
names(kids) <- kids; names(sids) <- sids          # 命名便于后续引用

# 检查每个样本中细胞类型的分布 -------------------------------------------
t(table(sce_sce$cluster_id, sce_sce$sample_id))  # 输出样本×细胞类型矩阵

# 生成伪批量数据（Pseudobulk） --------------------------------------------
# 将每个样本(sid)中每个细胞类型(kid)的counts求和，生成伪批量表达矩阵
# 参数说明：
# - assay: 使用counts矩阵
# - fun: 聚合函数（sum表示求和）
# - by: 按细胞类型(cluster_id)和样本(sample_id)分组聚合
pb <- aggregateData(sce_sce,
                    assay = "counts", 
                    fun = "sum",
                    by = c("cluster_id", "sample_id"))

# 可视化样本相似性 -------------------------------------------------------
# 使用MDS方法展示伪批量样本在基因表达空间的分布
(pb_mds <- pbMDS(pb))  # 生成MDS图，检查样本是否按预期分组聚集

# 设置分组因子顺序（重要：影响差异分析的方向） --------------------------
pb$group_id <- factor(pb$group_id, levels = c("Adult","Older"))

# 执行伪批量差异分析 -----------------------------------------------------
# 使用DESeq2方法，针对每个细胞类型单独进行Young vs Adult比较
# 参数说明：
# - method: 指定差异分析工具（DESeq2/edgeR/limma-voom）
# - verbose: 关闭冗余输出
res <- pbDS(pb, method = 'DESeq2', verbose = FALSE)

# 转换结果格式 -----------------------------------------------------------
# 临时复制sce对象处理counts矩阵格式（muscat要求）
tmp <- sce_sce
counts(tmp) <- as.matrix(counts(tmp))  # 确保counts是标准矩阵格式

# 提取差异分析结果表 -----------------------------------------------------
# 参数说明：
# - bind: 按行合并所有细胞类型结果
# - frq: 不添加基因表达频率信息
# - cpm: 不添加CPM值
result_table <- resDS(tmp, res, bind = "row", frq = FALSE, cpm = FALSE)
rm(tmp)  # 删除临时对象

# 结果解读说明 -----------------------------------------------------------
# 最终结果result_table包含以下关键列：
# - gene: 基因名称
# - cluster_id: 对应的细胞类型（如OPC）
# - p_val: 原始p值
# - p_adj: 多重检验校正后的p值（如FDR）
# - logFC: Young组相对于Adult组的对数倍变化（log2 fold change）

# 加载所需包 ----------------------------------------------------------------
library(openxlsx) # 用于写入Excel多工作表

# 创建新的Excel工作簿 ------------------------------------------------------
wb <- createWorkbook() # 初始化Excel工作簿对象

# 获取所有细胞类型名称 -----------------------------------------------------
cell_types <- unique(result_table$cluster_id) # 从结果表中提取唯一的细胞类型
cell_types <- sort(cell_types) # 按字母顺序排序（可选）

# 检查细胞类型数量是否匹配预期 ----------------------------------------------
message("检测到细胞类型数量：", length(cell_types)) # 应输出11
stopifnot(length(cell_types) == 11) # 确保数量正确
table(result_table$cluster_id)
# 循环处理每个细胞类型 -----------------------------------------------------
for (ct in cell_types) {
  # 提取当前细胞类型的差异结果
  df_sub <- result_table[result_table$cluster_id == ct, ]
  
  # 删除cluster_id列（因为每个工作表已经按类型分开）
  df_sub$cluster_id <- NULL 
  
  # 清理行名
  rownames(df_sub) <- NULL
  
  # 创建工作表（注意：Excel工作表名称不能超过31字符/包含特殊符号）
  sheet_name <- gsub("[^a-zA-Z0-9]", "_", ct) # 替换特殊字符为下划线
  sheet_name <- substr(sheet_name, 1, 31)     # 截断到31字符
  
  # 添加工作表到工作簿
  addWorksheet(wb, sheetName = sheet_name)
  
  # 写入数据（带格式的表格）
  writeDataTable(wb, 
                 sheet = sheet_name, 
                 x = df_sub,
                 startCol = 1, 
                 startRow = 1, 
                 tableStyle = "TableStyleMedium2")
  
  # 设置列宽自适应（可选）
  setColWidths(wb, sheet_name, cols = 1:ncol(df_sub), widths = "auto")
}

# 保存Excel文件 -----------------------------------------------------------
output_file <- "adult_vs_older_CellType_DEG_Results.xlsx" # 定义输出文件名
saveWorkbook(wb, file = output_file, overwrite = TRUE) # 保存文件（覆盖已存在）

# 输出完成信息 -----------------------------------------------------------
message("\nExcel文件已生成：", output_file)
message("包含以下工作表：\n", paste(cell_types, collapse = "\n"))
