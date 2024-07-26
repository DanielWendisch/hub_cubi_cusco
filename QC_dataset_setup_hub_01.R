# This script makes a Seurat object from data set with BPCell processing of data matrices(memory efficient)
# the Seurat dataset contains both raw (cellranger output) derived counts as well as cellbender corrected counts
# as separate assays
# 



# set to calc_bp_cells to TRUE if BPCells matrices should be recalculated and saved
calc_bp_cells <- FALSE
###----------------------------------------------------------------------------------------------------
library(readr)
library(dplyr)
library(tidyr)
library(purrr)
library(magrittr)
library(ggplot2)
library(Seurat)
library(tidyseurat)
library(BPCells)
#library(Azimuth)
library(stringr)
library(ggridges)
library(forcats)
library(rhdf5)
###----------------------------------------------------------------------------------------------------








#set paths for dataset-associated files
###----------------------------------------------------------------------------------------------------
dataset_name <- "hub_01"
path_raw_data_file <- "..\\..\\raw_data\\cubi_cusco"

path_cellranger_output <- paste(sep = "\\", path_raw_data_file, "hub_01_outs\\multi\\count\\raw_feature_bc_matrix.h5")
path_cellbender_output <-  paste(sep = "\\", path_raw_data_file,"cellbender\\hub_01\\hub_01_cellbender_corrected_filtered_seurat.h5")

## set hashtag conversion
convert_vector <- c("CMO302"="hIO_BIHi005-A",
                    "CMO301"="hIO_BIHi001-B",
                    "CMO303"="hIO_BIHi250-A",
                    "CMO304"="hIO_UCSFi001-A")

tag_calls_summary <- read_csv(paste(sep = "\\", path_raw_data_file,  "hub_01_outs/multi/multiplexing_analysis/tag_calls_summary.csv"))
tag_calls_per_cell <- read_csv(paste(sep = "\\", path_raw_data_file, "hub_01_outs/multi/multiplexing_analysis/tag_calls_per_cell.csv"))
assignment_confidence_table <- read_csv(paste(sep = "\\",path_raw_data_file,  "hub_01_outs/multi/multiplexing_analysis/assignment_confidence_table.csv"))



# write BPCells object, if not yet done
######################################################################################
# set to calc_bp_cells to TRUE if BPCells matrices should be recalculated and saved
# BPcells is needed for faster data handling on computers with little RAM
# execute only when no BP cell files have been saved, eg. first run on machine
# if matrices are alreadz stored, they will have to be moved before recalculation

if (calc_bp_cells) {
  
  if (!file.exists("BPcell_matrices")) {
    dir.create("BPcell_matrices")
  }
  
  # load uncorrected data
  ###----------------------------------------------------------------------------------------------------
  dat <- Read10X_h5(
    path_cellranger_output,
    use.names = TRUE)
  #separate the matrix in hashtag and gene reads
  
  cmo_dat <- dat$`Multiplexing Capture`
  dat <- dat$`Gene Expression`
  
  # Write the matrix to a directory
  write_matrix_dir(
    mat = dat,
    dir = paste("BPcell_matrices/", dataset_name))
  
  write_matrix_dir(
    mat = cmo_dat,
    dir = paste0("BPcell_matrices/", "cmo_", dataset_name))
  
  ###--------------------------------------------------------------------------------------------------
  
  
  # load cellbender corrected data and save BPcells matrix
  ###----------------------------------------------------------------------------------------------
  
  dat_cellbender <- Read10X_h5(filename = path_cellbender_output,unique.features = TRUE)
  
  setdiff(rownames(dat_cellbender)[!(grepl("ENSG", rownames(dat_cellbender)))],
          rownames(dat)) |> length()
  
  write_matrix_dir(dat_cellbender, dir = paste0("BPcell_matrices/" , dataset_name, "/cellbender"))
  }
#############################################################################################



# Create Seurat Matrix from stored BPcell matrices
###----------------------------------------------------------------------------------------------------

dat_cellbender <- open_matrix_dir(dir = paste0("BPcell_matrices/" , dataset_name, "/cellbender"))

dat <- open_matrix_dir(dir = paste("BPcell_matrices/", dataset_name))

intersecting_cells <- intersect(colnames(dat_cellbender), colnames(dat))
seurat_obj <- CreateSeuratObject(dat_cellbender[,intersecting_cells])

dat_raw <- CreateAssay5Object(dat[,intersecting_cells])

seurat_obj[["uncorrected_RNA"]] <- dat_raw

dat_cellbender <- NULL
dat<- NULL
gc()
###-----------------------------------------------------------------------------------





# basic seurat processing for each assay
###--------------------------------------------------------------------------------------------
for (assay_it in Assays(seurat_obj)) {
  DefaultAssay(seurat_obj) <- assay_it
  seurat_obj <- seurat_obj |> 
    NormalizeData() 
  
  seurat_obj <- seurat_obj |>
    FindVariableFeatures(assay = assay_it) |> 
    ScaleData() 
  #print("runPCA")
  seurat_obj <- seurat_obj |>
    RunPCA(assay = assay_it,reduction.name = paste0("pca.",assay_it)) 
  #print("neighbors")
  seurat_obj <- seurat_obj |> FindNeighbors(dims = 1:30,
                                                  reduction = paste0("pca.",assay_it),
                                                  assay = assay_it) 
  
  #print("clusters")
  seurat_obj <- seurat_obj |>  FindClusters(
    resolution = 0.8,
    verbose = FALSE,
    cluster.name =paste0(assay_it, "_clusters_leiden_res0.8"),algorithm=4
    #graph.name = paste0("kNN_","pca.",assay_it)
  ) 
  
  seurat_obj <- seurat_obj |>
    RunUMAP(dims = 1:30, reduction = paste0("pca.",assay_it),
            assay = assay_it,
            reduction.name =paste0("umap.",assay_it) )
}
###-----------------------------------------------------------------------------------




#add metadata
###-----------------------------------------------------------------------------------

# make metadata column: cell_line_cellranger 
#-- match CMOxx to cell line name, adjusted feature_call column from cellranger multi
cell_barcodes <- tag_calls_per_cell |>
  right_join(tibble(cell_barcode=colnames(seurat_obj)), keep=T)  |>
  mutate(cell_line_cellranger=str_replace_all(feature_call,convert_vector) |>
           replace_na("no cell"))


# make metadata column: singlet_doublet_cellranger 
#-- cellranger calls with levels doublet, singlet no cell(not called a cell from cellranger and therefore not in file per_sample_outs but in file multi>count)
#adjusted feature_call column from cellranger multi
cell_barcodes <-  cell_barcodes |> mutate(singlet_doublet_cellranger = ifelse(grepl("\\|", cell_line_cellranger), "doublet", 
                                                                              ifelse(cell_line_cellranger=="no cell","no cell","singlet"
                                                                              )))
# make metadata column: cell_line_doublet_cellranger 
#--- like singlet_doublet_cellranger but with cell line names
cell_barcodes <-  cell_barcodes |> mutate(cell_line_doublet_cellranger = ifelse(grepl("\\|", cell_line_cellranger), "doublet",cell_line_cellranger))

# add relevant dataframe columns to Seurat Metadata
cell_barcodes <- cell_barcodes |>
  as.data.frame()
rownames(cell_barcodes) <- cell_barcodes$cell_barcode.y
cell_barcodes <- cell_barcodes[,c("cell_line_cellranger","singlet_doublet_cellranger","cell_line_doublet_cellranger")]
seurat_obj <- AddMetaData(seurat_obj,cell_barcodes)

# saveRDS(
#   object = seurat_obj,
#   file = "intermediate_data/seurat_obj_hub_01_raw_and_corrected.Rds"
# )

# add cellbender cell probabilities
corrected = H5Fopen("C:\\Users\\Danne\\raw_data\\cubi_cusco\\cellbender\\hub_01\\hub_01_cellbender_corrected.h5")
cell_prob_tbl <- tibble(
  cell=corrected$metadata$barcodes_analyzed,
  cell_probability= corrected$droplet_latents$cell_probability)

H5Fclose(corrected)

cell_prob_tbl <- cell_prob_tbl |> filter(cell %in% colnames(seurat_obj))
cell_prob_vec <- cell_prob_tbl |> pull(cell_probability)
names(cell_prob_vec) <- cell_prob_tbl |> pull(cell)
seurat_obj <- AddMetaData(seurat_obj, col.name = "cellbender_prob_to_be_cell", cell_prob_vec)
###-----------------------------------------------------------------------------------




















