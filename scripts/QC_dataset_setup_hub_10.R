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
calc_bp_cells <- TRUE
add_vireo <- TRUE
leiden_or_louvain <- "louvain"
basic_cluster_name <- paste("clusters", leiden_or_louvain, "res0.8", sep="_")
basic_cluster_name_2 <- paste("clusters", leiden_or_louvain, "res0.4", sep="_")
dataset_name <- "hub_10"

### ----------------------------------------------------------------------------------------------------
#renv::use_python()
library(reticulate)
library(readr)
library(dplyr)
library(tidyr)
library(purrr)
library(magrittr)
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
library(future)

# plan("multisession", workers = 4)
# plan()

#renv::use_python(type="auto")
### ----------------------------------------------------------------------------------------------------








# set paths for dataset-associated files
### ----------------------------------------------------------------------------------------------------

path_raw_data_file <- "..\\..\\raw_data\\cubi_cusco"

path_cellranger_output <- paste(sep = "\\", path_raw_data_file, "cellranger\\hub_10_outs\\multi\\count\\raw_feature_bc_matrix.h5")
path_cellbender_output <- paste(sep = "\\", path_raw_data_file, "cellbender\\hub_10\\hub_10_cellbender_corrected_filtered_seurat.h5")

path_vireo_output <- paste(sep = "\\", path_raw_data_file, "vireo", paste(dataset_name, "vireo_donor_ids.tsv", sep="_"))



## set hashtag conversion
convert_vector <- c("CMO301"= "hLOA_BIHi001-B",
                    "CMO302"= "hLOA_BIHi005-A",
                    "CMO303"= "hLOA_BIHi250-A",
                    "CMO304"= "hLOA_UCSFi001-A")

tag_calls_summary <- read_csv(paste(sep = "\\", path_raw_data_file, "cellranger\\hub_10_outs\\multi\\multiplexing_analysis/tag_calls_summary.csv"))
tag_calls_per_cell <- read_csv(paste(sep = "\\", path_raw_data_file, "cellranger\\hub_10_outs\\multi\\multiplexing_analysis/tag_calls_per_cell.csv"))
assignment_confidence_table <- read_csv(paste(sep = "\\", path_raw_data_file, "cellranger\\hub_10_outs\\multi\\multiplexing_analysis\\assignment_confidence_table.csv"))



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
    dir = paste("BPcell_matrices/", dataset_name),overwrite = T
  )

  write_matrix_dir(
    mat = cmo_dat,
    dir = paste0("BPcell_matrices/", "cmo_", dataset_name),overwrite = T
  )

  ### --------------------------------------------------------------------------------------------------


  # load cellbender corrected data and save BPcells matrix
  ### ----------------------------------------------------------------------------------------------

  dat_cellbender <- Read10X_h5(filename = path_cellbender_output, unique.features = TRUE)

  write_matrix_dir(dat_cellbender$`Gene Expression`, dir = paste0("BPcell_matrices/", dataset_name, "/cellbender"), ,overwrite = T)
}
#############################################################################################



# Create Seurat Matrix from stored BPcell matrices
### ----------------------------------------------------------------------------------------------------

dat_cellbender <- open_matrix_dir(dir = paste0("BPcell_matrices/", dataset_name, "/cellbender"))
dat <- open_matrix_dir(dir = paste("BPcell_matrices/", dataset_name))
dat_cmo <- open_matrix_dir(dir = paste0("BPcell_matrices/", "cmo_", dataset_name))


reads_per_cell <- dat |> colSums()

plot_read_count_knee(reads_per_cell, cutoff = 2e3)+
  plot_read_count_knee(reads_per_cell, cutoff = 1e3)+
  plot_read_count_knee(reads_per_cell, cutoff = 5e2)+
  plot_read_count_knee(reads_per_cell, cutoff = 1e2)+
  plot_annotation(title="RNA reads per cell")

ggsave(filename = here("output/",dataset_name, c(paste0(dataset_name,"_read_count_knee.png"))), width = 22)


reads_per_cell <- dat_cmo |> colSums()

plot_read_count_knee(reads_per_cell, cutoff = 2e3)+
  plot_read_count_knee(reads_per_cell, cutoff = 1e3)+
  plot_read_count_knee(reads_per_cell, cutoff = 5e2)+
  plot_read_count_knee(reads_per_cell, cutoff = 1e2)+
  plot_annotation(title="hashtag reads per cell")

ggsave(filename = here("output/",dataset_name, c(paste0(dataset_name,"_hashtag_count_knee.png"))), width = 22)



#dat <- dat[,reads_per_cell>1e3]

intersecting_cells <- intersect(colnames(dat_cellbender), colnames(dat))
intersecting_cells <- intersect(intersecting_cells, colnames(dat_cmo))

seurat_obj <- CreateSeuratObject(dat_cellbender[, intersecting_cells])

dat_raw <- CreateAssay5Object(dat[, intersecting_cells])

seurat_obj[["uncorrected_RNA"]] <- dat_raw

dat_cmo <- dat_cmo[, intersecting_cells]
dat_cmo_tbl <- as.data.frame(t(as.matrix(dat_cmo)))
dat_cmo <- CreateAssay5Object(dat_cmo)

seurat_obj[["hashtag_oligos"]] <- dat_cmo

dat_cmo_tbl <- dat_cmo_tbl[,names(convert_vector)]
colnames(dat_cmo_tbl) <- paste("hashtag_id",convert_vector, sep="_")
seurat_obj <- AddMetaData(object = seurat_obj, dat_cmo_tbl)

dat_cellbender <- NULL
dat <- NULL
gc()
### ---------------------------------------------------------------------------------------------------------



# basic seurat processing for each assay
### -------------------------------------------------------------------------------------------------------
for (assay_it in Assays(seurat_obj)) {
  DefaultAssay(seurat_obj) <- assay_it
  seurat_obj <- seurat_obj |>
    NormalizeData()
  
  seurat_obj <- seurat_obj |>
    FindVariableFeatures(assay = assay_it) |>
    ScaleData()
  print("runPCA")
  seurat_obj <- seurat_obj |>
    RunPCA(assay = assay_it, reduction.name = paste0("pca.", assay_it),npcs =30)
  print("neighbors")
  seurat_obj <- seurat_obj |> FindNeighbors(
    dims = 1:15,
    reduction = paste0("pca.", assay_it),
    assay = assay_it
  )
  gc()
  
  print("clusters")
  
  seurat_obj <- seurat_obj |> FindClusters(
    resolution = 0.8,
    verbose = TRUE,
    cluster.name = paste(assay_it, basic_cluster_name,sep = "_"), 
    algorithm = ifelse(leiden_or_louvain=="louvain",1,4)
    # graph.name = paste0("kNN_","pca.",assay_it)
  )
  
  seurat_obj <- seurat_obj |>
    RunUMAP(
      dims = 1:15, reduction = paste0("pca.", assay_it),
      assay = assay_it,
      reduction.name = paste0("umap.", assay_it)
    )
  
  gc()
}
### -------------------------------------------------------------------------------------------------




# add metadata
### ---------------------------------------------------------------------------------------------------

# make metadata column: cell_line_cellranger
#-- match CMOxx to cell line name, adjusted feature_call column from cellranger multi
cell_barcodes <- tag_calls_per_cell |>
  right_join(tibble(cell_barcode = colnames(seurat_obj)), keep = T) |>
  mutate(cell_line_cellranger = str_replace_all(feature_call, convert_vector) |>
           replace_na("no cell"))


# make metadata column: singlet_doublet_cellranger
#-- cellranger calls with levels doublet, singlet no cell(not called a cell from cellranger and therefore not in file per_sample_outs but in file multi>count)
# adjusted feature_call column from cellranger multi
cell_barcodes <- cell_barcodes |> mutate(singlet_doublet_cellranger = ifelse(grepl("\\|", cell_line_cellranger), "doublet",
                                                                             ifelse(cell_line_cellranger == "no cell", "no cell", "singlet")
))
# make metadata column: cell_line_doublet_cellranger
#--- like singlet_doublet_cellranger but with cell line names
cell_barcodes <- cell_barcodes |> mutate(cell_line_doublet_cellranger = ifelse(grepl("\\|", cell_line_cellranger), "doublet", cell_line_cellranger))

# add relevant dataframe columns to Seurat Metadata
cell_barcodes <- cell_barcodes |>
  as.data.frame()
rownames(cell_barcodes) <- cell_barcodes$cell_barcode.y
cell_barcodes <- cell_barcodes[, c("cell_line_cellranger", "singlet_doublet_cellranger", "cell_line_doublet_cellranger")]
seurat_obj <- AddMetaData(seurat_obj, cell_barcodes)



# add cellbender cell probabilities
#### TODO: change path to here()

corrected <- H5Fopen(here(path_raw_data_file,"cellbender","hub_10","hub_10_cellbender_corrected.h5"))

cell_prob_tbl <- tibble(
  cell = corrected$metadata$barcodes_analyzed,
  cell_probability = corrected$droplet_latents$cell_probability
)
#plot cell bender probabilities
cell_prob_tbl |> ggplot(aes( -rank(cell_probability),cell_probability))+geom_point()+theme_bw()+ggtitle("CellBender probability to be a cell")
ggsave(filename = here("output/",dataset_name, c(paste0(dataset_name,"CellBender probability.png"))))

H5Fclose(corrected)

cell_prob_tbl <- cell_prob_tbl |> filter(cell %in% colnames(seurat_obj))
cell_prob_vec <- cell_prob_tbl |> pull(cell_probability)
names(cell_prob_vec) <- cell_prob_tbl |> pull(cell)
seurat_obj <- AddMetaData(seurat_obj, col.name = "cellbender_prob_to_be_cell", cell_prob_vec)




write_rds(seurat_obj, "intermediate_data/QC_dataset_setup_hub_10.rds")
### -----------------------------------------------------------------------------------------------------




# find markers for every cluster compared to all remaining cells, report only the positive
# ones
### ------------------------------------------------------------------------------------------------------
min_pct <- 0.4
logfc_threshold <- 0.25
max_cells_ <- 300

Idents(seurat_obj) <- paste0("RNA", basic_cluster_name,"_")
seurat_obj.markers <- FindAllMarkers(seurat_obj,
                                     only.pos = TRUE,
                                     min.pct = min_pct,
                                     logfc.threshold = logfc_threshold,
                                     max.cells.per.ident = max_cells_,
) |>
  as_tibble()

write_rds(seurat_obj.markers,
          file = paste0(
            "output/",
            dataset_name,
            "QC_cluster_markers_min.pct_", min_pct,
            "_logfc.threshold_", logfc_threshold,
            "_max.cells.per.ident_", max_cells_,
            ".rds"
          )
)

write_csv(seurat_obj.markers,
          file = paste0(
            "output/",
            dataset_name,
            "QCmarkers_min.pct_", min_pct,
            "_logfc.threshold_", logfc_threshold,
            "_max.cells.per.ident_", max_cells_,
            ".csv"
          )
)



###VIREO -----------------------------------------------------------------------------------------------------

seurat_obj <- read_rds( "intermediate_data/QC_dataset_setup_hub_10_vireo_filtered.rds")


if (add_vireo) {
  vireo_ids <- read_tsv(path_vireo_output)
  
  vireo_ids |> ggplot(aes(donor_id)) + geom_histogram(stat="count")+ theme_bw()
  
  vireo_ids |> ggplot(aes(prob_max))+geom_histogram()+facet_wrap(vars(donor_id))
  

  
  

  
  vireo_ids <- vireo_ids |> select(cell, donor_id) |> right_join(tibble(cell=colnames(seurat_obj)))
  vireo_ids$donor_id |> unique()
  vireo_ids_vec <- pull(vireo_ids, donor_id)
  names(vireo_ids_vec) <- pull(vireo_ids, cell)
  seurat_obj <- AddMetaData(seurat_obj,col.name ="vireo_snp_demux",vireo_ids_vec) 
  
  
  seurat_obj |>  select(vireo_snp_demux, contains("hashtag_id") ) |>
    pivot_longer(cols = contains("hashtag_id"), names_to = "hashtag") |> 
    ggplot(aes(hashtag, log(value+1)))+geom_violin()+
    facet_wrap(vars(vireo_snp_demux))
  
  seurat_obj |> select(vireo_snp_demux)
  #seurat_obj@meta.data[["vireo_snp_demux_celltype"]] 
  
  seurat_obj |> group_by(vireo_snp_demux) |> count(cell_line_doublet_cellranger)
  seurat_obj |> group_by(cell_line_doublet_cellranger) |> count(vireo_snp_demux)
  seurat_obj |> ggplot(aes(vireo_snp_demux, log(nCount_RNA)))+ geom_violin()
  
  
  seurat_obj_vireo_filtered <- seurat_obj |> filter(vireo_snp_demux  %in% c("donor0","donor1","donor2","donor3"))

  
  
  for (assay_it in Assays(seurat_obj_vireo_filtered)) {
    
    if (assay_it=="hashtag_oligos") {dimensions <- 1:length(convert_vector)}else{dimensions <- c(1:15)}
    
    DefaultAssay(seurat_obj_vireo_filtered) <- assay_it
    seurat_obj_vireo_filtered <- seurat_obj_vireo_filtered |>
      NormalizeData()
    
    seurat_obj_vireo_filtered <- seurat_obj_vireo_filtered |>
      FindVariableFeatures(assay = assay_it) |>
      ScaleData()
    print("runPCA")
    seurat_obj_vireo_filtered <- seurat_obj_vireo_filtered |>
      RunPCA(assay = assay_it, reduction.name = paste0("pca.", assay_it),npcs =max(dimensions))
    
    print("neighbors")
    seurat_obj_vireo_filtered <- seurat_obj_vireo_filtered |> FindNeighbors(
      dims = dimensions,
      reduction = paste0("pca.", assay_it),
      assay = assay_it
    )
    gc()
    
    print("clusters")
    
    seurat_obj_vireo_filtered <- seurat_obj_vireo_filtered |> FindClusters(
      resolution = 0.8,
      verbose = TRUE,
      cluster.name = paste(assay_it, basic_cluster_name,sep = "_"), 
      algorithm = ifelse(leiden_or_louvain=="louvain",1,4)
      # graph.name = paste0("kNN_","pca.",assay_it)
    )
    
    
    seurat_obj_vireo_filtered <- seurat_obj_vireo_filtered |> FindClusters(
      resolution = 0.2,
      verbose = TRUE,
      cluster.name = paste(assay_it, basic_cluster_name_2,sep = "_"), 
      algorithm = ifelse(leiden_or_louvain=="louvain",1,4)
      # graph.name = paste0("kNN_","pca.",assay_it)
    )
    
    
    seurat_obj_vireo_filtered <- seurat_obj_vireo_filtered |>
      RunUMAP(
        dims = dimensions,
        reduction = paste0("pca.", assay_it),
        assay = assay_it,
        reduction.name = paste0("umap.", assay_it)
      )
    
    gc()
  }
  
  seurat_obj_vireo_filtered |> DimPlot(group.by ="vireo_snp_demux", reduction = "pca.hashtag_oligos")
  seurat_obj_vireo_filtered |> FeaturePlot("hashtag_id_hLOA_BIHi001-B")
  seurat_obj_vireo_filtered |> VlnPlot("hashtag_id_hLOA_BIHi001-B")
  
  write_rds(seurat_obj_vireo_filtered, "intermediate_data/QC_dataset_setup_hub_10_vireo_filtered.rds")
  
  
}


seurat_obj_vireo_filtered |> VlnPlot("hashtag_id_hLOA_BIHi001-B")

DefaultAssay(seurat_obj_vireo_filtered) <- "hashtag_oligos"
seurat_obj_vireo_filtered |> DimPlot(group.by ="vireo_snp_demux" ) +
  seurat_obj_vireo_filtered |> DimPlot(group.by ="cell_line_doublet_cellranger" )
                                     
DefaultAssay(seurat_obj_vireo_filtered) <- "uncorrected_RNA"

seurat_obj_vireo_filtered |> DimPlot(group.by ="vireo_snp_demux" )

seurat_obj_vireo_filtered |> FeaturePlot("hashtag_id_hLOA_BIHi001-B" )


seurat_obj_vireo_filtered |> DimPlot(group.by ="RNA_clusters_louvain_res0.4", label = T)


seurat_obj_vireo_filtered |> FeaturePlot("EPCAM")

seurat_obj_vireo_filtered |> FeaturePlot("SOX11")

seurat_obj_vireo_filtered |> FeaturePlot("CDH7")

seurat_obj_vireo_filtered |> FeaturePlot("RPL15")


seurat_obj_vireo_filtered |> FeaturePlot("KCND2")


seurat_obj_vireo_filtered |> FeaturePlot("MUC2")




markers_04 <- seurat_obj_vireo_filtered |> FindAllMarkers(group.by ="RNA_clusters_louvain_res0.4" )


markers_04 |> as_tibble() |> slice_min(p_val,by = cluster, n = 10) |> View()











