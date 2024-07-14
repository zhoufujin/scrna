# 读取之前保存的 Seurat 对象
seurat_obj <- readRDS("seurat_obj_step1.rds")

# 进行细胞簇的聚类
seurat_obj <- FindNeighbors(seurat_obj, dims = 1:30)
seurat_obj <- FindClusters(seurat_obj, resolution = 0.5)

# 运行 UMAP（或 t-SNE）进行降维可视化
seurat_obj <- RunUMAP(seurat_obj, dims = 1:30)

# 绘制 UMAP 图
umap_plot <- DimPlot(seurat_obj, reduction = 'umap', label = TRUE)

# 定义常见标志基因
markers <- c("CD3D", "MS4A1", "LYZ", "PPBP", "CD8A")

# 绘制标记基因的点图
dot_plot <- DotPlot(seurat_obj, features = markers) + RotatedAxis()

# 根据标志基因进行细胞类型注释
# 获取现有聚类的数量
num_clusters <- length(levels(seurat_obj))

# 定义新的聚类名称（示例名称，可以根据你的数据和分析需求调整）
new_cluster_names <- c("T cells", "B cells", "Monocytes", "Platelets", "CD8 T cells")
# 扩展或截断 new_cluster_names 以匹配 num_clusters 的长度
new_cluster_names <- rep(new_cluster_names, length.out = num_clusters)

# 设置新聚类的名字
names(new_cluster_names) <- levels(seurat_obj)
seurat_obj <- RenameIdents(seurat_obj, new_cluster_names)

# 保存 Seurat 对象到 .rds 文件
saveRDS(seurat_obj, file = "seurat_obj_step3.rds")

# 保存图形对象到 .rds 文件
saveRDS(list(umap_plot = umap_plot, dot_plot = dot_plot), file = "seurat_plots_step3.rds")
