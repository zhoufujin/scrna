library(GEOquery)
library(curl)

# 定义显示下载进度的函数
curl_progress <- function() {
  function(down, dnow, up, unow) {
    if (down > 0) {
      pct <- dnow / down * 100
      cat(sprintf("\rDownload progress: %.2f%%, Downloaded %d of %d bytes", pct, dnow, down))
      flush.console()
    }
  }
}

download_with_progress <- function(geo_id) {
  # 获取 GEO 数据集的下载 URL
  urls <- getGEOSuppFiles(geo_id, baseDir = tempdir())
  url <- as.character(urls[1, "url"])

  # 定义下载文件的路径
  destfile <- file.path(tempdir(), paste0(geo_id, ".tar.gz"))

  # 检查文件是否已经部分下载
  resume_from <- 0
  if (file.exists(destfile)) {
    resume_from <- file.info(destfile)$size
    cat(sprintf("Resuming download from %d bytes\n", resume_from))
  } else {
    cat("Starting new download\n")
  }

  # 设置 curl 句柄
  h <- new_handle()
  handle_setopt(h, progressfunction = curl_progress())
  handle_setopt(h, noprogress = FALSE)
  if (resume_from > 0) {
    handle_setopt(h, resume_from_large = resume_from)
  }

  # 下载文件
  tryCatch({
    curl_download(url, destfile, handle = h)
    cat("\nDownload completed successfully\n")
  }, error = function(e) {
    cat("\nDownload failed: ", e$message, "\n")
  })

  # 使用 GEOquery 从下载的文件读取数据
  if (file.exists(destfile) && file.info(destfile)$size == 654202880) {
    gse <- getGEO(filename = destfile, GSEMatrix = TRUE)

    # 检查并处理数据
    if (length(gse) > 1) {
      gse <- gse[[1]]
    }

    return(gse)
  } else {
    stop("下载数据集失败，请检查网络连接或 GEOquery 包的版本。")
  }
}

# 下载示例数据集
gse <- download_with_progress("GSE189357")

# 检查对象是否成功创建
if (exists("gse")) {
    # 查看数据集的元数据
    meta_data <- pData(phenoData(gse))

    # 查看数据集的表达矩阵
    expr_data <- exprs(gse)

    # 查看元数据的前几行
    head(meta_data)

    # 查看表达矩阵的前几行和前几列
    head(expr_data[, 1:5])
} else {
    message("下载数据集失败，请检查网络连接或 GEOquery 包的版本。")
}

# 查看临时目录
tempdir()

# 列出临时目录中的文件
list.files(tempdir())
