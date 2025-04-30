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
break

for (assay_it in Assays(seurat_obj)) {
  print("new assay #######################################################################")
  print("          #######################################################################")
  print(assay_it)
  
  if (assay_it=="hashtag_oligos") {dimensions <- 1:length(convert_vector)}else{dimensions <- c(1:dim_number_pca_and_neighbors)}
  
  DefaultAssay(seurat_obj) <- assay_it
  seurat_obj <- seurat_obj |>
    NormalizeData()
  
  print("scale data")
  seurat_obj <- seurat_obj |>
    FindVariableFeatures(assay = assay_it) |>
    ScaleData()
  
  print("runPCA")
  seurat_obj[[assay_it]]$scale.data  <- seurat_obj[[assay_it]]$scale.data %>% write_matrix_memory(compress=FALSE)
  seurat_obj <- seurat_obj |>
    RunPCA(assay = assay_it, reduction.name = paste0("pca.", assay_it),npcs =max(dimensions))
  
  print("neighbors")
  seurat_obj <- seurat_obj |> FindNeighbors(
    dims = dimensions,
    reduction = paste0("pca.", assay_it),
    assay = assay_it)
  gc()
  
  print("clustering no 1")
  seurat_obj <- seurat_obj |> FindClusters(
    resolution = resolution_clustering_1,
    verbose = TRUE,
    cluster.name = paste(assay_it, basic_cluster_name,sep = "_"), 
    algorithm = ifelse(leiden_or_louvain=="louvain",1,4)
    # graph.name = paste0("kNN_","pca.",assay_it)
  )
  
  print("clustering no 1")
  seurat_obj <- seurat_obj |> FindClusters(
    resolution = resolution_clustering_2,
    verbose = TRUE,
    cluster.name = paste(assay_it, basic_cluster_name_2,sep = "_"), 
    algorithm = ifelse(leiden_or_louvain=="louvain",1,4)
    # graph.name = paste0("kNN_","pca.",assay_it)
  )
  
  
  seurat_obj <- seurat_obj |>
    RunUMAP(
      dims = dimensions,
      reduction = paste0("pca.", assay_it),
      assay = assay_it,
      reduction.name = paste0("umap.", assay_it)
    )
  
  gc()
}

DefaultAssay(seurat_obj) <- "RNA"