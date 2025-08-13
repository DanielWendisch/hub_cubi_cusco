
dataset_name <- "hub_10"

seurat_obj <- read_rds(here("intermediate_data",paste0("QC_dataset_setup_", dataset_name  ,".rds")))

# to fix the double clusterinfo issue
# seurat_obj@meta.data <- seurat_obj@meta.data[,-13]
# seurat_obj |> write_rds(here("intermediate_data",paste0("QC_dataset_setup_", dataset_name  ,".rds")))

seurat_obj@meta.data <- seurat_obj@meta.data[,(names(seurat_obj@meta.data) |> unique())]

#plots <- readRDS "output/plots/hub_01_QC_plots.rds")
#seurat_obj <- read_rds(here("intermediate_data", "QC_dataset_setup_hub_10_vireo_filtered.rds"))


leiden_or_louvain <- "louvain"
basic_cluster_name <- paste("RNA_clusters", leiden_or_louvain, "res0.8", sep="_")
basic_cluster_name_2 <- paste("clusters", leiden_or_louvain, "res0.4", sep="_")

uncorrected_cluster_name <- paste("RNA_clusters", leiden_or_louvain, "res0.8", sep="_")


#plots
#######################################



plot_specs <- tribble(
  ~default_assay,   ~list_int,    ~meta_data_column,    ~plote_title,                                         ~order,                           ~subtitle,
  
  "cellbender_RNA",    100,   basic_cluster_name,  " clusters (0.8)",                                             NULL,          "background corrected counts (cellbender)",
  
  "cellbender_RNA",    1,    "cell_line_cellranger",          "cell line (cellranger)",                                    NULL,      "background corrected counts (cellbender)",
  
  "cellbender_RNA",    2,    "singlet_doublet_cellranger",    "singlets, doublets, empty droplets (hashtag call by cellranger)", "doublet",  "background corrected counts (cellbender)",
  
  "cellbender_RNA",    3,    "cell_line_doublet_cellranger",  "cell line, singlets, doublets (hashtag call by cellranger)", NULL,    "background corrected counts (cellbender)",
  
  "cellbender_RNA",    200,   uncorrected_cluster_name,  " clusters (0.8)",                                             NULL,          "background uncorrected counts",
  
  "RNA",    21,    "cell_line_cellranger",          "cell line (cellranger)",                                    NULL,     "background uncorrected counts",
  
  "RNA",    22,    "singlet_doublet_cellranger",    "singlets, doublets, empty droplets  (hashtag call by cellranger)", "doublet",  "background uncorrected counts",
  
  "RNA",    23,    "cell_line_doublet_cellranger",  "cell line, singlets, doublets (hashtag call by cellranger)", NULL,    "background uncorrected counts",
  
  "RNA",    300,   uncorrected_cluster_name,  " clusters on uncorrected RNA(0.8)",                                             NULL,          "background corrected counts (cellbender)",
  
)



dataset_content <- "Alveolar Lung Organoids"
lower_nCount_RNA_thresh_quantile <- NA
lower_nCount_RNA_thresh <- 3000
upper_nCount_RNA_thresh_quantile<- NA
upper_nCount_RNA_thresh<- 42000
upper_mito_thresh_quantile<- NA
upper_mito_thresh <- 7.5
deleted_clusters<- c(1) #c(11,12)


##
log_n_count_rna_plot_limits <- c(0.001,log(1.05*(max(seurat_obj$nCount_RNA))))


n_count_rna_plot_limits <- c(1,1.05*(max(seurat_obj$nCount_RNA)))
