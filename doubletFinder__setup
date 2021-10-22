.libPaths("/home/rstudio/R/x86_64-pc-linux-gnu-library/4.0")
#-------------------------------------------------------------------------------------------------------------
library("Seurat")
library("glmGamPoi")
library("harmony")
#-------------------------------------------------------------------------------------------------------------
sobj = readRDS("/home/rstudio/rstudio-docker-sessions/01__mecfs/data/post__cellCycleScoring.RDS")

sobj.list <- SplitObject(sobj, split.by="orig.ident")

sobj.list <- lapply(X = sobj.list, FUN = function(x){
  
  SCTransform(object = x, method = "glmGamPoi", vars.to.regress = "percent.mt", return.only.var.genes = FALSE)
  
  RunPCA(x, verbose = FALSE, dims = 1:50)
  RunUMAP(x, dims = 1:50)
  
  FindNeighbors(x, dims = 1:50)
  FindClusters(x, resolution = 0.6)
  
  })

#-------------------------------------------------------------------------------------------------------------
m= 7.6652E-6
b= -3.57897E-07
sobj.list = lapply(sobj.list, function(x){
  AddMetaData(object = x, metadata = round(m * nrow(x@meta.data) + b, 2), col.name = "expectedRate"  )
})

#-------------------------------------------------------------------------------------------------------------
save.image("INPUT_DOUBLET_FINDER_split__by__orig.ident.Rdata")

