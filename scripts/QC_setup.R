# This script makes a Seurat object from data set with BPCell processing of data matrices(memory efficient)
# the Seurat dataset contains both raw (cellranger output) derived counts as well as cellbender corrected counts
# as separate assays
#

## inputs:
#     path_cellranger_output <- paste(sep = "\\", path_raw_data_file, "hub_10_outs\\multi\\count\\raw_feature_bc_matrix.h5")
#     path_cellbender_output <-  paste(sep = "\\", path_raw_data_file,"cellbender\\hub_10\\hub_10_cellbender_corrected_filtered_seurat.h5")
#     convert_vector <- c("CMO308"="EC_UCSFi001-A",
#                     "CMO307"="EC_BIHi250-A",
#                     "CMO306"="EC_BIHi005-A",
#                     "CMO305"="EC_BIHi001-B",
#                     "CMO309"="EC_BIHi001-A")
#
#     tag_calls_summary <- read_csv(paste(sep = "\\", path_raw_data_file,  "hub_10_outs/multi/multiplexing_analysis/tag_calls_summary.csv"))
#     tag_calls_per_cell <- read_csv(paste(sep = "\\", path_raw_data_file, "hub_10_outs/multi/multiplexing_analysis/tag_calls_per_cell.csv"))
#     assignment_confidence_table <- read_csv(paste(sep = "\\",path_raw_data_file,  "hub_10_outs/multi/multiplexing_analysis/assignment_confidence_table.csv"))

## outputs:
# BPCells matrix
#
# write_rds(seurat_obj.markers,
#           file = paste0("output/",
#                         dataset_name,
#                         "QC_cluster_markers_min.pct_",min_pct,
#                         "_logfc.threshold_",logfc_threshold,
#                         "_max.cells.per.ident_",max_cells_,
#                         ".rds"))
#
# write_csv(seurat_obj.markers,
#           file = paste0("output/",
#                         dataset_name,
#                         "QCmarkers_min.pct_",min_pct,
#                         "_logfc.threshold_",logfc_threshold,
#                         "_max.cells.per.ident_",max_cells_,
#                         ".csv"))
# write_rds(seurat_obj, "intermediate_data/QC_dataset_setup_hub_10.rds")


# set to calc_bp_cells to TRUE if BPCells matrices should be recalculated and saved


### libraries----------------------------------------------------------------------------------------------------
#renv::use_python()
library(reticulate)
library(readr)
library(dplyr)
library(tidyr)
library(purrr)
library(ggplot2)
library(Seurat)
library(tidyseurat)
library(BPCells)
# library(Azimuth)
library(stringr)
library(ggridges)
library(forcats)
library(rhdf5)
library(here)
library(patchwork)
#library(future)

# plan("multisession", workers = 4)
# plan()

#renv::use_python(type="auto")

### set paths for dataset-associated files 
#----------------------------------------------------------------------------------------------------

## set hashtag conversion
#TODO move these fully to master script

# convert_vector <- c("CMO301"= "hLOA_BIHi001-B",
#                     "CMO302"= "hLOA_BIHi005-A",
#                     "CMO303"= "hLOA_BIHi250-A",
#                     "CMO304"= "hLOA_UCSFi001-A")


# input data fiels, dir paths specified in master
path_cellranger_output <- here( path_cellranger_output_dir,"count","raw_feature_bc_matrix.h5")
path_cellbender_output <- here(path_cellbender_output_dir,paste0(dataset_name, "_cellbender_corrected_filtered_seurat.h5"))
path_vireo_output <-  here(path_vireo_output_dir, paste0(dataset_name,"_vireo_donor_ids.tsv"))

#read-in csv
tag_calls_summary <- read_csv(here(path_cellranger_output_dir,"multiplexing_analysis","tag_calls_summary.csv"))
tag_calls_per_cell <- read_csv(here(path_cellranger_output_dir,"multiplexing_analysis","tag_calls_per_cell.csv"))
assignment_confidence_table <- read_csv(here(path_cellranger_output_dir,"multiplexing_analysis","assignment_confidence_table.csv"))

path_bp_cells_cellranger <- here("BPcell_matrices", paste0("cellranger_", dataset_name))
path_bp_cells_cellranger_cmo <- here("BPcell_matrices", paste0("cellranger_cmo_", dataset_name))
path_bp_cells_cellbender <- here("BPcell_matrices", paste0("cellbender_", dataset_name))


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
  ### ----------------------------------------------------------------------------------------------------
  dat <- Read10X_h5(
    path_cellranger_output,
    use.names = TRUE
  )
  # separate the matrix in hashtag and gene reads

  cmo_dat <- dat$`Multiplexing Capture`
  dat <- dat$`Gene Expression`

  # Write the matrix to a directory
  write_matrix_dir(
    mat = dat,
    dir = path_bp_cells_cellranger,
    overwrite = T
  )

  write_matrix_dir(
    mat = cmo_dat,
    dir = path_bp_cells_cellranger_cmo,
    overwrite = T
  )

  ### --------------------------------------------------------------------------------------------------


  # load cellbender corrected data and save BPcells matrix
  ### ----------------------------------------------------------------------------------------------

  dat_cellbender <- Read10X_h5(filename = path_cellbender_output, unique.features = TRUE)

  write_matrix_dir(dat_cellbender,#$`Gene Expression`,
                   path_bp_cells_cellbender, overwrite = T) # TODO harmomize between dataset $Gene expression
}
#############################################################################################








### # Create Seurat Matrix from stored BPcell matrices ---------------------------------------------------------------
# Load data
dat_cellbender <- open_matrix_dir(dir = path_bp_cells_cellbender)
dat_ranger <- open_matrix_dir(dir = path_bp_cells_cellranger)
dat_cmo <- open_matrix_dir(dir = path_bp_cells_cellranger_cmo)




## ###########plots for library evaluation

###
rna_reads_per_cell <- dat_ranger |> colSums()

knee_plot_rna <- plot_read_count_knee(rna_reads_per_cell, cutoff = 2e3)+
  plot_read_count_knee(rna_reads_per_cell, cutoff = 1e3)+
  plot_read_count_knee(rna_reads_per_cell, cutoff = 5e2)+
  plot_read_count_knee(rna_reads_per_cell, cutoff = 1e2)+
  plot_annotation(title="RNA reads per cell")


cmo_reads_per_cell <- dat_cmo |> colSums()

knee_plot_cmo <- plot_read_count_knee(cmo_reads_per_cell, cutoff = 2e3)+
  plot_read_count_knee(cmo_reads_per_cell, cutoff = 1e3)+
  plot_read_count_knee(cmo_reads_per_cell, cutoff = 5e2)+
  plot_read_count_knee(cmo_reads_per_cell, cutoff = 1e2)+
  plot_annotation(title="hashtag reads per cell")


###
ggsave(knee_plot_rna, filename = here("output",dataset_name, "qc_process", "plots_and_plot_data", paste(dataset_name ,"rna_read_count_knee.png", sep = "_")))
ggsave(knee_plot_cmo, filename = here("output",dataset_name,"qc_process", "plots_and_plot_data", paste(dataset_name ,"hashtag_read_count_knee.png", sep = "_")))


############ #inital trimming of data
dat_ranger <- dat_ranger[,(rna_reads_per_cell>counts_per_cell_pre_filter)]

#cells that are present in celranger counts &cellranger hashtag & cellbender output
intersecting_cells <- intersect(colnames(dat_cellbender), colnames(dat_ranger))
intersecting_cells <- intersect(intersecting_cells, colnames(dat_cmo))


seurat_obj <- CreateSeuratObject(dat_ranger[, intersecting_cells])
assay_cellbender <- CreateAssay5Object(dat_cellbender[, intersecting_cells])

seurat_obj[["cellbender_RNA"]] <- assay_cellbender

dat_cmo <- dat_cmo[, intersecting_cells]
dat_cmo_tbl <- as.data.frame(t(as.matrix(dat_cmo)))
dat_cmo <- CreateAssay5Object(dat_cmo)

seurat_obj[["hashtag_oligos"]] <- dat_cmo

# add hashtag counts to metadata
dat_cmo_tbl <- dat_cmo_tbl[,names(convert_vector)]
colnames(dat_cmo_tbl) <- paste("hashtag_id",convert_vector, sep="_")
seurat_obj <- AddMetaData(object = seurat_obj, dat_cmo_tbl)

#clean up memory
dat_cellbender <- NULL
dat_ranger <- NULL
dat_cmo <- NULL
gc()


### vireo ---------------------------------------------------------------------------------------------------------


if (add_vireo) {
  vireo_ids <- read_tsv(path_vireo_output)
  
  vireo_ids <- vireo_ids |> select(cell, donor_id) |>
    right_join(tibble(cell=colnames(seurat_obj))) |>
    mutate(donor_id=if_else(is.na(donor_id), "no vireo data",donor_id))
  
  vireo_ids_vec <- pull(vireo_ids, donor_id)
  names(vireo_ids_vec) <- pull(vireo_ids, cell)
  seurat_obj <- AddMetaData(seurat_obj,
                            col.name ="vireo_snp_demux",
                            vireo_ids_vec) 
  
  #seurat_obj <- seurat_obj |> filter(vireo_snp_demux  %in% c("donor0","donor1","donor2","donor3"))
  
  seurat_obj <- seurat_obj |> mutate(vireo=case_when(vireo_snp_demux == "unassigned" ~ "vireo.unassigned" ,
                         is.na(vireo_snp_demux) ~ "vireo.no data",
                         TRUE ~ "vireo.cell"))
} else { # some weird error when using else
  seurat_obj$vireo <- "not_computed"
  seurat_obj$vireo_snp_demux <- "not_computed"
  
}






# basic seurat processing for each assay
### -------------------------------------------------------------------------------------------------------

resolution_clustering_1 <- resolution_clustering_1
resolution_clustering_2 <- resolution_clustering_2
dim_number_pca_and_neighbors <- dim_number_pca_and_neighbors

source(here("scripts", "basic_seurat_processing_for_each_assay.R"))

### -------------------------------------------------------------------------------------------------




# add metadata
### ---------------------------------------------------------------------------------------------------

# make metadata column: cell_line_cellranger
#-- match CMOxx to cell line name, adjusted feature_call column from cellranger multi
cell_barcodes <- tag_calls_per_cell |>
  right_join(tibble(cell_barcode = colnames(seurat_obj)), keep = T) |>
  mutate(cell_line_cellranger = str_replace_all(feature_call, convert_vector) |>
           replace_na("no cell"))


# METADATA: singlet_doublet_cellranger
#-- cellranger calls with levels doublet, singlet no cell(not called a cell from cellranger and therefore not in file per_sample_outs but in file multi>count)
# adjusted feature_call column from cellranger multi

cell_barcodes <- cell_barcodes |> mutate(singlet_doublet_cellranger = ifelse(grepl("\\|", cell_line_cellranger), "doublet",
                                                                             ifelse(cell_line_cellranger == "no cell", "no cell", "singlet")))
# METADATA: cell_line_doublet_cellranger
#--- like singlet_doublet_cellranger but with cell line names
cell_barcodes <- cell_barcodes |> mutate(cell_line_doublet_cellranger = ifelse(grepl("\\|", cell_line_cellranger), "doublet", cell_line_cellranger))

# add relevant dataframe columns to Seurat Metadata
cell_barcodes <- cell_barcodes |>
  as.data.frame()
rownames(cell_barcodes) <- cell_barcodes$cell_barcode.y
cell_barcodes <- cell_barcodes[, c("cell_line_cellranger", "singlet_doublet_cellranger", "cell_line_doublet_cellranger")]
seurat_obj <- AddMetaData(seurat_obj, cell_barcodes)



### add cellbender cell probabilities ----
#--------------------------------------------------------

corrected <- H5Fopen(here(path_cellbender_output_dir, paste0(dataset_name,"_cellbender_corrected.h5")))
cell_prob_tbl <- tibble(
  cell = corrected$metadata$barcodes_analyzed,
  cell_probability = corrected$droplet_latents$cell_probability)

H5Fclose(corrected)

cell_prob_tbl <- cell_prob_tbl |> filter(cell %in% colnames(seurat_obj))
cell_prob_vec <- cell_prob_tbl |> pull(cell_probability)
names(cell_prob_vec) <- cell_prob_tbl |> pull(cell)

# METADATA cellbender_prob_to_be_cell
# from cellbender output  H5Fopen cell = corrected$metadata$barcodes_analyzed, cell_probability = corrected$droplet_latents$cell_probability
seurat_obj <- AddMetaData(seurat_obj, col.name = "cellbender_prob_to_be_cell", cell_prob_vec)

# METADATA cellbender_p95
# METADATA vireo
# METADATA cellranger

seurat_obj <- seurat_obj |> 
  mutate(cellbender_p95=if_else(cellbender_prob_to_be_cell > 0.95, "cellbender.higher0.95" , "cellbender.lower0.95")) |> 

  mutate(cellranger=paste("cellranger", singlet_doublet_cellranger, sep="."))

# METADATA cellranger percent_mito
seurat_obj[["percent_mito"]] <- PercentageFeatureSet(object = seurat_obj, pattern = "^MT-", assay = "RNA")



seurat_obj@misc$filtering_steps <- list()
seurat_obj@misc$filtering_steps$qc_setup <-paste("counts_per_cell_pre_filter >" ,as.character(counts_per_cell_pre_filter))


seurat_obj$organoid <-  seurat_obj |> separate(cell_line_cellranger, sep = "_", into = "organoid")|> pull(organoid)

seurat_obj$cell_line <-  seurat_obj |> separate(cell_line_cellranger, sep = "_", into = c("organoid", "cell_line")) |> pull(cell_line)





