process_seurat<-function(obj, method, ref_datasets = NULL, k.anchor=5,k.weight=100,
                         res=NULL, features=NULL, dims=NULL, 
                         batch=NULL, return_model = F, cluster=T, 
                         type = "seur", nfeats = 2000, neighbor=F) {
  
  if(type == "sce") {
    obj <- CreateSeuratObject(counts = counts(obj))
  } else {
    obj <- obj
  }
  
  
  if(!is.null(features)) {
    features <- rownames(obj)[!grepl(features, rownames(obj))]
    obj <- subset(obj, features = features)
  } 
  
  if(method == "integrate") {
    obj <- .integrate_seurat(obj, split = batch, nfeats = nfeats, ref_datasets=ref_datasets, k.anchor = k.anchor, k.weight = k.weight)
    DefaultAssay(obj) <- "integrated"
    obj <- ScaleData(obj) %>% RunPCA(.)
  } else if (method == "log") {
    DefaultAssay(obj) <- "RNA"
    obj <-
      NormalizeData(obj) %>%
      FindVariableFeatures(., selection.method = "vst", nfeatures = nfeats) %>%
      ScaleData(., vars.to.regress=batch) %>%
      RunPCA(.)
  } else if (method == "glm") {
    obj <-
      SCTransform(obj, method = "glmGamPoi", batch_var=batch, variable.features.n = nfeats) %>%
      RunPCA(.)
  } else if (method == "qpoisson") {
    obj <-
      SCTransform(obj, method = "qpoisson", variable.features.n = nfeats, vars.to.regress = batch) %>%
      RunPCA(.)
  }
  
  if(cluster==T & method == "integrate"){
    obj <-
      obj %>%
      RunUMAP(., dims = seq(dims), return.model=return_model) %>%
      FindNeighbors(., dims = seq(dims), return.neighbor=neighbor) %>%
      FindClusters(., resolution = res)
  } else if(cluster==T & method != "integrate") {
    obj <-
      obj %>%
      RunUMAP(., dims = seq(dims), return.model=return_model) %>%
      FindNeighbors(., dims = seq(dims), return.neighbor=neighbor) %>%
      FindClusters(., resolution = res)
  }
  
  return(obj)
}