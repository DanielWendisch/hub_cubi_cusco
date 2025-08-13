#library(sceasy)
library(reticulate)
library(here)
library(tidyverse)
library(Seurat)
library(SeuratDisk)
library(Matrix)

output_path <- here("output", "prepare_for_cellxgene")
dir.create(output_path)


convert_all_layers_to_dgCMatrix <- function(seurat_obj){
  for (Assay in Assays(seurat_obj)) {
    cat(Assay)
    cat("counts")
    seurat_obj[[Assay]]$counts <- as(object = seurat_obj[[Assay]]$counts, Class = "dgCMatrix",)
    cat("scale.data")
    seurat_obj[[Assay]]$scale.data <- as(object = seurat_obj[[Assay]]$scale.data, Class = "dgCMatrix")
    cat("data")
    seurat_obj[[Assay]]$data <- as(object = seurat_obj[[Assay]]$data, Class = "dgCMatrix")
    
  }
  return(seurat_obj) 
}



# use_condaenv('EnvironmentName')
# loompy <- reticulate::import('loompy')


hub_01_10_11[["cellbender_RNA"]] <- NULL
hub_01_10_11[["hashtag_oligos"]] <- NULL
hub_01_10_11[["mnn.reconstructed"]] <- NULL


#hub_01_10_11 <- convert_all_layers_to_dgCMatrix(hub_01_10_11)

hub_01_10_11 <- convert_all_layers_to_dgCMatrix(hub_01_10_11)
hub_01_10_11[["RNA"]]$counts %>% str()
hub_01_10_11 %>% str()
#hub_01_10_11[["RNA"]] <- as(object = hub_01_10_11[["RNA"]], Class = "Assay")

hub_01_10_11 <- UpdateSeuratObject(hub_01_10_11)


library(tidyseurat)
hub_01_10_11$organoid <-  hub_01_10_11 |> separate(cell_line_cellranger, sep = "_", into = "organoid") |> pull(organoid)

hub_01_10_11$cell_line <-  hub_01_10_11 |> separate(cell_line_cellranger, sep = "_", into = c("organoid", "cell_line")) |> pull(cell_line)

hub_01_10_11 |> write_rds(here(output_path, "seurat_obj_hub_01_10_11_harmony_mnn_noBPCell.rds"))


# This is done in a separate environment on wsl
# SaveH5Seurat(hub_01_10_11, filename =here(output_path, "hub_01_10_11.h5Seurat"))
# Convert(here(output_path, "hub_01_10_11.h5Seurat"), dest = "h5ad")
# 
# 
# 

from_script <- "lung_and_intestinal_organoid_integration"
hub_01_10_11 <-  read_rds(here("output", from_script, "seurat_obj_hub_01_10_11_harmony_mnn.rds"))

hub_01_10_11 <- hub_01_10_11 |> JoinLayers()

hub_01_10_11[["RNA"]]$counts <- as(object = hub_01_10_11[["RNA"]]$counts, Class = "matrix")
hub_01_10_11[["RNA"]]$data <- as(object = hub_01_10_11[["RNA"]]$data, Class = "matrix")
hub_01_10_11[["RNA"]]$scale.data <- as(object = hub_01_10_11[["RNA"]]$scale.data, Class = "matrix")



hub_01_10_11 <- hub_01_10_11 |> JoinLayers()
hub_01_10_11[["cellbender_RNA"]] <- NULL
hub_01_10_11[["hashtag_oligos"]] <- NULL
hub_01_10_11[["mnn.reconstructed"]] <- NULL

library(scCustomize)
DefaultAssay(hub_01_10_11)
data_v3 <- Convert_Assay(hub_01_10_11, convert_to = "V3")






