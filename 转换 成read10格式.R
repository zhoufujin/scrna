# 1. 移动文件到各自的目录
move_files_to_dirs <- function() {
  files <- list.files(path = "GSE189357_RAW", pattern = "*.tsv.gz$|*.mtx.gz$", full.names = TRUE)

  for (file in files) {
    sample_name <- sub("_[^_]+$", "", basename(file))
    sample_dir <- file.path("GSE189357_RAW", sample_name)

    if (!dir.exists(sample_dir)) {
      dir.create(sample_dir, recursive = TRUE)
    }

    file.rename(file, file.path(sample_dir, basename(file)))
  }
}

# 运行文件移动函数
move_files_to_dirs()

# 2. 重命名每个目录中的文件以符合 Read10X 的要求
rename_files_in_dirs <- function() {
  sample_dirs <- list.dirs(path = "GSE189357_RAW", full.names = TRUE, recursive = FALSE)

  for (dir in sample_dirs) {
    files <- list.files(path = dir, full.names = TRUE)

    for (file in files) {
      if (grepl("barcodes.tsv.gz$", file)) {
        file.rename(file, file.path(dir, "barcodes.tsv.gz"))
      } else if (grepl("features.tsv.gz$", file)) {
        file.rename(file, file.path(dir, "features.tsv.gz"))
      } else if (grepl("matrix.mtx.gz$", file)) {
        file.rename(file, file.path(dir, "matrix.mtx.gz"))
      }
    }
  }
}

# 运行重命名函数
rename_files_in_dirs()

# 检查重命名后的文件
sample_dirs <- list.dirs(path = "GSE189357_RAW", full.names = TRUE, recursive = FALSE)
for (dir in sample_dirs) {
  cat("Checking directory:", dir, "\n")
  print(list.files(dir))
}

# 3. 读取10X数据
library(Seurat)

# 获取所有样本目录
sample_dirs <- list.dirs(path = "GSE189357_RAW", full.names = TRUE, recursive = FALSE)

# 读取每个样本目录中的数据
sce_list <- lapply(sample_dirs, function(dir) {
  Read10X(data.dir = dir)
})

# 查看结果
sce_list

#
