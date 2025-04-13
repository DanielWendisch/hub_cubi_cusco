import os

# Define variables
datasets = "hub_01"
RAW_DATA_PATH = os.path.join("..", "..", "raw_data", "cubi_cusco")
# output_files = [f"QC_dashboard_{dataset}.html" for dataset in datasets]

# Additional parameters
min_pct = 0.1
logfc_threshold = 0.5
max_cells_ = 100

# Rule to render QC_dashboard.qmd for each dataset
quit
    output:
        markers_rds=os.path.join(
            "output",
            f"{datasets}QC_cluster_markers_min.pct_{min_pct}_logfc.threshold_{logfc_threshold}_max.cells.per.ident_{max_cells_}.rds"
        ),
        seurat_obj=os.path.join("intermediate_data", f"QC_dataset_setup_{datasets}.rds")
    shell:
        """
        1_QC_dataset_setup.R
        """
