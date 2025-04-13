library(here)
library(Seurat)
library(BPCells)
library(scCustomize)
library(reticulate)
# py_install(AnnData)
py_require("addata")


#############################
### FUCTIONS##################
convert_.rds_name_to_.h5ad <- function(filename) {
  if (grepl("\\.rds$", filename)) {
    sub("\\.rds$", ".h5ad", filename)
  } else {
    stop("Filename does not end with '.rds'")
  }
}


convert_from_bpcells_format <- function(seurat_obj) {
  #as(seurat_obj[[assay]]@layers[[layer]], "Assay")
  for (assay in Assays(seurat_obj)) {
    
    if (length(Layers(seurat_obj[[assay]]))>1) {
      for (layer in Layers(seurat_obj[[assay]])) {
        #browser()
        seurat_obj[[assay]]@layers[[layer]] <- as(GetAssayData(seurat_obj, assay = assay, layer = layer), Class = "dgCMatrix")
        print(paste(assay, layer))
        print("_________________________")
      }
    }
      else { # some not so good handling when assay has no layers
        print(paste("nolayers in", assay, layer))
        if (Layers(seurat_obj[[assay]])[[1]]=="data") {
          print("data")
          browser()
          seurat_obj[[assay]]@data <- as(GetAssayData(seurat_obj, assay = assay, slot = "data"), Class = "dgCMatrix")
          print("_________________________")
        }
    }
    
  }
  return(seurat_obj)
}

###############################
###############################

# read in .rds from intermediate data directory and make BPCells related data access to normal matrix
filenames <- c("hub_01_filtered_seurat_obj_processed.rds",
              "hub_01_integrated_.rds",
              "hub_02_filtered_seurat_obj_processed.rds")

for (filename in filenames) {
  path_in <- here("intermediate_data", filename)
  directory_out <- here("intermediate_data", "anndata")
  path_out<- here(directory_out, convert_.rds_name_to_.h5ad(filename))
  if (file.exists(path_out)) {
    print(paste(path_out, "exists","________________________","skipped to next"))
    next}
  obj <- readRDS(path_in) 
  obj <- obj |> convert_from_bpcells_format()
  
  # convert to anndata and save in intermediate/anndata with .h5ad file name extension
  as.anndata(obj, file_path = directory_out, file_name = convert_.rds_name_to_.h5ad(filename))
}




hub_01_integrated_$


Layers(obj["RNA"])

Layers(obj[["RNA"]])
Layers(obj[["uncorrected_RNA"]])
Layers(obj[["mnn.reconstructed"]])


Assays(obj)



slot()obj# saveRDS(hub_01,file=here("intermediate_data",converted_filename))
