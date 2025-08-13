clustering_fine_RNA <- paste("RNA", "clusters", leiden_or_louvain, "res." , resolution_clustering_1, sep = "_")
clustering_rough_RNA <- paste("RNA", "clusters", leiden_or_louvain, "res." , resolution_clustering_2, sep = "_")
clustering_fine_cellbender <- paste("cellbender_RNA", "clusters", leiden_or_louvain, "res." , resolution_clustering_1, sep = "_")
clustering_rough_cellbender <- paste("cellbender_RNA", "clusters", leiden_or_louvain, "res." , resolution_clustering_2, sep = "_")






plot_specs <- tribble(
  ~default_assay,   ~list_int,    ~meta_data_column,   ~title_short,                         ~plote_title,                                         ~order,                           ~subtitle,
  
  "cellbender_RNA",    100,   clustering_rough_RNA,     "clus0.8_cellbender_emedding"  ,      " clusters (0.8)",                                             NULL,          "background corrected counts (cellbender)",
  
  "cellbender_RNA",    1,    "cell_line_cellranger",      "cell_line_cellbender_embedding"  ,     "cell line (cellranger)",                                    NULL,      "background corrected counts (cellbender)",
  
  "cellbender_RNA",    2,    "singlet_doublet_cellranger", "singlets_doublets_cellbender_embedding"  ,    "singlets, doublets, empty droplets (hashtag call by cellranger)", "doublet",  "background corrected counts (cellbender)",
  
  "cellbender_RNA",    3,    "cell_line_doublet_cellranger","celline_doublets_cellbender_embedding",   "cell_line_singlets_doublets_cellbender_embedding", NULL,    "background corrected counts (cellbender)",
  
  "cellbender_RNA",    200,   clustering_rough_cellbender,   "clus0.8_cellbender_embedding",     " clusters (0.8)",                                             NULL,          "background uncorrected counts",
  
  "RNA",    21,    "cell_line_cellranger",        "cell_line",          "cell line (cellranger)",                                    NULL,     "background uncorrected counts",
  
  "RNA",    22,    "singlet_doublet_cellranger",  "singelts_doublets",  "singlets, doublets, empty droplets  (hashtag call by cellranger)", "doublet",  "background uncorrected counts",
  
  "RNA",    23,    "cell_line_doublet_cellranger", "cell_line", "cell line, singlets, doublets (hashtag call by cellranger)", NULL,    "background uncorrected counts",
  
  "RNA",    300,   clustering_rough_cellbender, "clus0.8",  " clusters on uncorrected RNA(0.8)",                                             NULL,          "background corrected counts (cellbender)",
  
)
