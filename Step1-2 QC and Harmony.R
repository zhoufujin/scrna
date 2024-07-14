# 加载必要的包
library(Seurat)
library(harmony)
# 设置工作目录到数据文件夹
setwd('E:/其他文件/scrna/GSE189357_RAW/')

data_dir <- "GSM5699777_TD1"
print(data_dir)


sample_data <-  Read10X(data.dir = data_dir)
seurat_obj <- CreateSeuratObject(counts = sample_data)
# 检查数据
# 检查数据结构
str(sample_data)
print(seurat_obj)

# 数据预处理
seurat_obj <- NormalizeData(seurat_obj)
seurat_obj <- FindVariableFeatures(seurat_obj)
seurat_obj <- ScaleData(seurat_obj)

# 质控和标准化
seurat_obj <- RunPCA(seurat_obj,npcs = 30)

# 保存Seurat对象
saveRDS(seurat_obj, file = "seurat_obj_step1.rds")
# 打印Seurat对象的摘要
print(seurat_obj)


