BiocManager::install("GEOquery")
library(GEOquery)
gset <- getGEO("GSE10072", GSEMatrix = TRUE)
str(gset)
exprs(gset[[1]])

# View the first few rows of the expression data
head(exprs(gset[[1]]))

# Access the phenotype data
pData(gset[[1]])

# View the first few rows of the phenotype data
head(pData(gset[[1]]))