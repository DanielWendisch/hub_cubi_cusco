#input variables

# resolution_clustering_1
# resolution_clustering_2
# dim_number_pca_and_neighbors
# 
# basic_cluster_name
# basic_cluster_name_2
#leiden_or_louvain

# TODO change assay names to variables 

print("run script basic_seurat_processing_for_each_assay.R")


for (assay_iter in Assays(seurat_obj)) {
  print("new assay #######################################################################")
  print("          #######################################################################")
  print(assay_iter)
  
  if (assay_iter=="hashtag_oligos") {dimensions <- 1:length(convert_vector)}else{dimensions <- c(1:dim_number_pca_and_neighbors)}
  
  DefaultAssay(seurat_obj) <- assay_iter
  seurat_obj <- seurat_obj |>
    NormalizeData()
  
  print("scale data")
  seurat_obj <- seurat_obj |>
    FindVariableFeatures(assay = assay_iter) |>
    ScaleData()
  
  print("runPCA")
  seurat_obj[[assay_iter]]$scale.data  <- seurat_obj[[assay_iter]]$scale.data %>% write_matrix_memory(compress=FALSE)
  seurat_obj <- seurat_obj |>
    RunPCA(assay = assay_iter, reduction.name = paste0("pca.", assay_iter),npcs =max(dimensions))
  
  print("neighbors")
  seurat_obj <- seurat_obj |> FindNeighbors(
    dims = dimensions,
    reduction = paste0("pca.", assay_iter),
    assay = assay_iter)
  gc()
  
  print("clustering no 1")
  seurat_obj <- seurat_obj |> FindClusters(
    resolution = resolution_clustering_1,
    verbose = TRUE,
    cluster.name = paste(assay_iter, "clusters", leiden_or_louvain, "res." , resolution_clustering_1, sep = "_"), 
    algorithm = ifelse(leiden_or_louvain=="louvain",1,4)
    # graph.name = paste0("kNN_","pca.",assay_iter)
  )
  
  print("clustering no 2")
  seurat_obj <- seurat_obj |> FindClusters(
    resolution = resolution_clustering_2,
    verbose = TRUE,
    cluster.name = paste(assay_iter, "clusters", leiden_or_louvain, "res." , resolution_clustering_2, sep = "_"), 
    algorithm = ifelse(leiden_or_louvain=="louvain",1,4)
    # graph.name = paste0("kNN_","pca.",assay_iter)
  )
  
  
  seurat_obj <- seurat_obj |>
    RunUMAP(
      dims = dimensions,
      reduction = paste0("pca.", assay_iter),
      assay = assay_iter,
      reduction.name = paste0("umap.", assay_iter)
    )
  
  gc()
}

DefaultAssay(seurat_obj) <- "RNA"