

seurat_obj <- read_rds("intermediate_data/QC_dataset_setup_hub_02.rds")
#plots <- readRDS "output/plots/hub_01_QC_plots.rds")



#plots
#######################################



plot_specs <- tribble(
  ~default_assay,   ~list_int,    ~meta_data_column,    ~plote_title,                                         ~order,                           ~subtitle,
  
  "RNA",    100,   "RNA_clusters_leiden_res0.8",  "Leiden clusters (0.8)",                                             NULL,          "background corrected counts (cellbender)",
  
  "RNA",    1,    "cell_line_cellranger",          "cell line (cellranger)",                                    NULL,      "background corrected counts (cellbender)",
  
  "RNA",    2,    "singlet_doublet_cellranger",    "singlets, doublets, empty droplets (hashtag call by cellranger)", "doublet",  "background corrected counts (cellbender)",
  
  "RNA",    3,    "cell_line_doublet_cellranger",  "cell line, singlets, doublets (hashtag call by cellranger)", NULL,    "background corrected counts (cellbender)",
  
  "uncorrected_RNA",    200,   "uncorrected_RNA_clusters_leiden_res0.8",  "Leiden clusters (0.8)",                                             NULL,          "background uncorrected counts",
  
  "uncorrected_RNA",    21,    "cell_line_cellranger",          "cell line (cellranger)",                                    NULL,     "background uncorrected counts",
  
  "uncorrected_RNA",    22,    "singlet_doublet_cellranger",    "singlets, doublets, empty droplets  (hashtag call by cellranger)", "doublet",  "background uncorrected counts",
  
  "uncorrected_RNA",    23,    "cell_line_doublet_cellranger",  "cell line, singlets, doublets (hashtag call by cellranger)", NULL,    "background uncorrected counts",
  "RNA",    300,   "uncorrected_RNA_clusters_leiden_res0.8",  "Leiden clusters on uncorrected RNA(0.8)",                                             NULL,          "background corrected counts (cellbender)",
  
)


lower_nCount_RNA_thresh_quantile <- NA
lower_nCount_RNA_thresh <- 3000
upper_nCount_RNA_thresh_quantile<- NA
upper_nCount_RNA_thresh<- 36000
upper_mito_thresh_quantile<- NA
upper_mito_thresh <- 7.5
deleted_clusters<- c(11,12)