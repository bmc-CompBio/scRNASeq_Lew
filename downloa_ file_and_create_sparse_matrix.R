library(Matrix)

#download original digital expression matrix
download.file("https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE63472&format=file&file=GSE63472%5FP14Retina%5Fmerged%5Fdigital%5Fexpression%2Etxt%2Egz", 
              destfile = "GSE63472_P14Retina_merged_digital_expression.txt.gz", 
              mode = "wb") 

#cave! reading this table requires at least 32 GB of RAM! alternatively download ready-to-use sparse matrix (link in seurat script)
tsne_table<-read.table("GSE63472_P14Retina_merged_digital_expression.txt.gz", sep="\t", header=T, row.names=1)

tsne_table<-data.matrix(tsne_table)

sparse<-Matrix(tsne_table, sparse=T)

writeMM(obj = sparse, file="sparse_GSE63472_P14Retina_merged_digital_expression.mtx")

write(x = rownames(sparse), file = "genes.tsv")

write(x = colnames(sparse), file = "barcodes.tsv")

