plot_specs <- tribble(
  ~default_assay,   ~list_int,    ~meta_data_column,   ~title_short,                         ~plote_title,                                         ~order,                           ~subtitle,
  
  "cellbender_RNA",    100,   basic_cluster_name,     "clus0.8_cellbender_emedding"  ,      " clusters (0.8)",                                             NULL,          "background corrected counts (cellbender)",
  
  "cellbender_RNA",    1,    "cell_line_cellranger",      "cell_line_cellbender_embedding"  ,     "cell line (cellranger)",                                    NULL,      "background corrected counts (cellbender)",
  
  "cellbender_RNA",    2,    "singlet_doublet_cellranger", "singlets_doublets_cellbender_embedding"  ,    "singlets, doublets, empty droplets (hashtag call by cellranger)", "doublet",  "background corrected counts (cellbender)",
  
  "cellbender_RNA",    3,    "cell_line_doublet_cellranger","celline_doublets_cellbender_embedding",   "cell_line_singlets_doublets_cellbender_embedding", NULL,    "background corrected counts (cellbender)",
  
  "cellbender_RNA",    200,   uncorrected_cluster_name,   "clus0.8_cellbender_embedding",     " clusters (0.8)",                                             NULL,          "background uncorrected counts",
  
  "RNA",    21,    "cell_line_cellranger",        "cell_line",          "cell line (cellranger)",                                    NULL,     "background uncorrected counts",
  
  "RNA",    22,    "singlet_doublet_cellranger",  "singelts_doublets",  "singlets, doublets, empty droplets  (hashtag call by cellranger)", "doublet",  "background uncorrected counts",
  
  "RNA",    23,    "cell_line_doublet_cellranger", "cell_line", "cell line, singlets, doublets (hashtag call by cellranger)", NULL,    "background uncorrected counts",
  
  "RNA",    300,   uncorrected_cluster_name, "clus0.8",  " clusters on uncorrected RNA(0.8)",                                             NULL,          "background corrected counts (cellbender)",
  
)
