library(Seurat)

#download and setup
download.file("https://syncandshare.lrz.de/dl/fi2nue8KDJ2b25qDe24KED35", 
              destfile = "sparse_GSE63472_P14Retina_merged_digital_expression.mtx")

counts<-readMM("sparse_GSE63472_P14Retina_merged_digital_expression.mtx")
colnames<-read.table("barcodes.tsv")
rownames<-read.table("genes.tsv")
colnames(counts)<-colnames[,1]
rownames(counts)<-rownames[,1]

#create object
nbt<-CreateSeuratObject(raw.data = counts, 
                        project = "P14Retina")
#log transform object to make more normally distributed
nbt=NormalizeData(nbt)

#find highly variable genes as they contribute more to cell heterogeneity; necessary for further analysis
nbt=FindVariableGenes(nbt,mean.function = ExpMean, 
                      dispersion.function = LogVMR, 
                      do.plot = FALSE)
var.genes<-nbt@var.genes

#each var.gene of a single cell is scaled relative to the expression in all other cells
nbt=ScaleData(nbt, do.par = T, num.cores=4, genes.use = var.genes)

#runPCA analysis
nbt=RunPCA(nbt, do.print=F, pcs.compute = 100)

#calculates how significant the association of a gene is with its respective PC, takes some time 
nbt=JackStraw(nbt,num.replicate = 200,do.par = T, num.cores=4, num.pc = 30)

#runTSNE analysis
nbt=RunTSNE(nbt, dims.use=1:25, 
                max_iter=1000, 
                num_threads=4, #for parallel computing install openmp branch from jkrijthe/Rtsne package 
                perplexity=30)

#experimental; requires installation of umap-learn module for python
#nbt=RunUMAP(nbt,dims.use = 1:25, 
               # reduction.use = "pca", 
               # n_neighbors = 30, 
               # min_dist = 0.3)

#Find clusters using shared nearest neighbor method, each cell is assigned an identity corresponding to the cluster it is in 

nbt=FindClusters(nbt, genes.use = var.genes, 
                 reduction.type = "tsne",      #use pca or tsne; PCA takes way longer
                 dims.use = 1:2,               #as many dimensions as computed for pca or tsne respectively
                 resolution = 0.3,             #0.8 if pca, 0.3 if tsne (the higher the more fine clusters are found)
                 plot.SNN = T)                 



#usefuel functions
##SetIdent()       #Sets identity of given cells to a value you choose; you can create own clusters or collapse existing clusters
##RenameIdent()    #rename existing identities
##FindMarkers()    #find genes differentially expressed between two or more clusters
##PCElbowPlot()    #plots standard deviation of PCs


#find biomarkers for certain clusters
markers<-FindMarkers(nbt,13,21,                    #use identities according to clusters you want to compare
                     logfc.threshold=1.5)          #log foldchange cutoff

#plot examples
## "GLUL","AQP4","RLBP1","APOE" are markers for Müller cells, a major retinal cell type
TSNEPlot(nbt, do.label = T)

FeaturePlot(nbt, features.plot ="GLUL",            #which feature (e.g gene/cells) to plot
            cols.use = c("lightblue", "red"))


FeaturePlot(nbt, features.plot =c("AQP4","GLUL"),
                 overlay = T, 
                 cols.use = c("lightblue","red","yellow"))

VlnPlot(nbt,  features.plot = c("AQP4","GLUL"), 
              ident.include = c(13,21),             #use identities according to clusters you want to compare
              do.sort = T) 

VlnPlot(nbt,  features.plot = c("GLUL","AQP4","RLBP1","APOE"), 
              ident.include = c(13,21),             #use identities according to clusters you want to compare
              do.sort = T) 





