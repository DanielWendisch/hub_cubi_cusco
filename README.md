![workflow](workflow.PNG)
# hub_cubi_cusco

## meta_data
QC_setup.R
### cell_line_cellranger -- match CMOxx to cell line name, adjusted feature_call column from cellranger multi
### singlet_doublet_cellranger -- cellranger calls with levels doublet, singlet no cell(not called a cell from cellranger and therefore not in file per_sample_outs but in file multi>count)
### adjusted feature_call column from cellranger multi
### cell_line_doublet_cellranger --- like singlet_doublet_cellranger but with cell line names



# METADATA: singlet_doublet_cellranger
#-- cellranger calls with levels doublet, singlet no cell(not called a cell from cellranger and therefore not in file per_sample_outs but in file multi>count)
# adjusted feature_call column from cellranger multi

cell_barcodes <- cell_barcodes |> mutate(singlet_doublet_cellranger = ifelse(grepl("\\|", cell_line_cellranger), "doublet",
                                                                             ifelse(cell_line_cellranger == "no cell", "no cell", "singlet")))
# METADATA: cell_line_doublet_cellranger
#--- like singlet_doublet_cellranger but with cell line names
cell_barcodes <- cell_barcodes |> mutate(cell_line_doublet_cellranger = ifelse(grepl("\\|", cell_line_cellranger), "doublet", cell_line_cellranger))


 METADATA cellbender_prob_to_be_cell
# from cellbender output  H5Fopen cell = corrected$metadata$barcodes_analyzed, cell_probability = corrected$droplet_latents$cell_probability
seurat_obj <- AddMetaData(seurat_obj, col.name = "cellbender_prob_to_be_cell", cell_prob_vec)

# METADATA cellbender_p95
# METADATA vireo
# METADATA cellranger

seurat_obj <- seurat_obj |> 
  mutate(cellbender_p95=if_else(cellbender_prob_to_be_cell>0.95,"cellbender.higher0.95","cellbender.lower0.95")) |> 
  mutate(vireo=case_when(vireo_snp_demux == "unassigned" ~ "vireo.unassigned" ,
                         is.na(vireo_snp_demux) ~ "vireo.no data",
                         TRUE ~ "vireo.cell")) |> 
  mutate(cellranger=paste("cellranger",singlet_doublet_cellranger,sep="."))

# METADATA cellranger percent_mito
seurat_obj[["percent_mito"]] <- PercentageFeatureSet(object = seurat_obj, pattern = "^MT-")

# METADATA cellranger dead
seurat_obj$"dead" <- ifelse(pull(seurat_obj,
                                 percent_mito) > upper_mito_thresh,"dead","alive")






----
# METADATA cellranger dead
seurat_obj$"dead" <- ifelse(pull(seurat_obj,
                                 percent_mito) > upper_mito_thresh,"dead","alive")
