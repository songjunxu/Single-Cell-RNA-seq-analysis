#load packages
library(Seurat)
library(dplyr)
library(Matrix)
library(gdata)

#Read in NSCLC counts matrix
dirname <- "~/Documents/SingleCell"
counts_matrix_filename = paste0(dirname, "/filtered_gene_bc_matrices/GRCh38/")
counts <- Read10X(data.dir = counts_matrix_filename)  # Seurat function to read in 10x count data
# To minimize memory use on the docker - choose only the first 1000 cells
counts <- counts[,1:1000]
#Examine the sparse counts matrix
counts[1:10, 1:3]#most are zero=.
dim(counts) # report number of genes (rows) and number of cells (columns)
#object.size(counts) # size in bytes
#object.size(as.matrix(counts)) # size in bytes as dense matrix. Therefore choose sparse matrix

counts_per_cell <- Matrix::colSums(counts)#cell
counts_per_gene <- Matrix::rowSums(counts)
genes_per_cell <- Matrix::colSums(counts>0) # count gene only if it has non-zero reads mapped.
cells_per_gene <- Matrix::rowSums(counts>0) # only count cells where the gene is expressed

hist(log10(counts_per_cell+1),main='counts per cell',col='wheat')
hist(log10(genes_per_cell+1), main='genes per cell', col='wheat')
#plot counts gene distribution
plot(counts_per_cell, genes_per_cell, log='xy', col='wheat')
title('counts vs genes per cell')
#plot genes per cell, examine library complexity 
plot(sort(genes_per_cell), xlab='cell', log='y', main='genes per cell (ordered)')

#create a seurat object
#This object contains various “slots” (designated by seurat@slotname) that will store not only the raw count data, 
#but also the results from various computations below. 
seurat<-CreateSeuratObject(counts=counts, min.cells = 3, min.features  = 350, project = "10X_NSCLC")

##view the first 3 rows (genes) 3 columns (cells)
GetAssayData(object = seurat, slot = 'counts')[1:3, 1:3]
seurat@assays$RNA

# The number of genes and UMIs (nGene and nUMI) are automatically calculated for every object by Seurat.  For non-UMI
# data, nUMI represents the sum of the non-normalized values within a cell We calculate the percentage of mitochondrial
# genes here and store it in percent.mito using AddMetaData.  We use object@raw.data since this represents
# non-transformed and non-log-normalized counts The % of UMI mapping to MT-genes is a common scRNA-seq QC metric.
mito.genes <- grep(pattern = "^MT-", x = rownames(x = seurat@data), value = TRUE)
percent.mito <- Matrix::colSums(seurat@raw.data[mito.genes, ])/Matrix::colSums(seurat@raw.data)

# AddMetaData adds columns to object@meta.data, and is a great place to stash QC stats.  This also allows us to plot the
# metadata values using the Seurat's VlnPlot().
head(seurat@meta.data)  # Before adding