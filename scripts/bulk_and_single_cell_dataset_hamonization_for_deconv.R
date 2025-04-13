# reads in integrated single cell data set and matched bulkRNASeq data set
# hormonizes names of samples
# changes gene names of single cell data set to ensemble ids
# keep only genes (ensemble ids) that are shared in bulk and sc sample
# output: harmonized bulk matrix and harmonized *single cell experiment*, genes are nemad with ensemble ids
# 

library(Seurat)
library(tidyseurat)
library(tidyverse)
library(SingleCellExperiment)
library(tidySingleCellExperiment)
#library(biomaRt) needed but called directly because of conflict

dataset_name <- "hub_01"
path_raw_data_file <-  "..\\..\\raw_data\\cubi_cusco"

obj <- read_rds(paste0("intermediate_data/", dataset_name, "_integrated_", ".rds"))#sc data
combined_bulk<- read.csv(paste(path_raw_data_file,"hub_01_02_bulk\\counts_table_hub_01_02.csv", sep = "\\")) # bulk counts
hub_01_meta_data <-  read_csv( "intermediate_data/hub_01_meta_data.csv")


## select  bulk samples of experiment according to sodar meta data

hub_01_bulk_library_names <- hub_01_meta_data |> pull(sample) |> str_replace_all("-",".")

combined_bulk <- combined_bulk |>  select(-gene,-X) #leave only counts and ensemble id
bulk_mtx <- combined_bulk[,hub_01_bulk_library_names] |> as.matrix() |>  t() 
colnames(bulk_mtx) <- combined_bulk |> pull(ensembl_gene_id)
bulk_mtx <- t(bulk_mtx)


## chose cell state labels, remove clusters

### what is this for?
cell.state.labels <- obj |> pull(harmony_clusters_res0.4) |> as.character()
table(cell.state.labels)

## make single cell mtx with ensemble ids from seurat object

sc_mtx <- obj@assays$uncorrected_RNA@layers$counts |> as.matrix() 
colnames(sc_mtx) <- colnames(obj) #cells
rownames(sc_mtx) <- rownames(obj) #features

# make named vactor for conversion of esnemble ids and hgnc gene names
if (!(file.exists("intermediate_data/gene_ensg_mapping_vec.rds")&
    file.exists("intermediate_data/ensg_gene_mapping_vec.rds"))) {
  # Use the Ensembl database for humans (change the dataset for different species)
  ensembl <- biomaRt::useMart("ensembl", dataset = "hsapiens_gene_ensembl")
  # Query Ensembl to get the Ensembl Gene IDs for your list of gene names
  gene_mapping <- biomaRt::getBM(
    filters = "hgnc_symbol",    # Filter based on HGNC gene symbols
    attributes = c("hgnc_symbol", "ensembl_gene_id"),  # Retrieve HGNC symbols and Ensembl Gene IDs
    values = rownames(sc_mtx),        # The list of gene names to map
    mart = ensembl
  )
  # make named vactor for conversion
  gene_ensg_mapping_vec <- gene_mapping$ensembl_gene_id
  names(gene_ensg_mapping_vec)<- gene_mapping$hgnc_symbol
  
  ensg_gene_mapping_vec <- gene_mapping$hgnc_symbol
  names(ensg_gene_mapping_vec)<- gene_mapping$ensembl_gene_id
  
  
  #write to disk for later 
  write_rds(gene_ensg_mapping_vec, "intermediate_data/gene_ensg_mapping_vec.rds")
  write_rds(ensg_gene_mapping_vec, "intermediate_data/ensg_gene_mapping_vec.rds")
}else  {
  gene_ensg_mapping_vec <-  read_rds( "intermediate_data/gene_ensg_mapping_vec.rds")
  ensg_gene_mapping_vec <- read_rds("intermediate_data/ensg_gene_mapping_vec.rds")
  }


rownames(sc_mtx) <- ifelse(rownames(sc_mtx) %in% names(gene_ensg_mapping_vec),
                           gene_ensg_mapping_vec[rownames(sc_mtx)],
                           rownames(sc_mtx))

# remove gene names that are not in conversion vector
sc_mtx <- sc_mtx[rownames(sc_mtx) %in% gene_ensg_mapping_vec,]
#sc_mtx <- t(sc_mtx)



bulk_mtx <- bulk_mtx[rowSums(bulk_mtx)>20,]
# keep only genes (ensemble ids) that are shared in bulk and sc sample
ensg_intersect_to_keep <- intersect(rownames(sc_mtx),rownames(bulk_mtx))
bulk_mtx <- bulk_mtx[ensg_intersect_to_keep,]
sc_mtx <- sc_mtx[ensg_intersect_to_keep,]

meta_data_arranged <- obj@meta.data[colnames(sc_mtx),]
meta_data_arranged <- meta_data_arranged |> bind_cols(obj@reductions$umap.mnn@cell.embeddings[colnames(sc_mtx),]) 

sce <- SingleCellExperiment(assays = list(counts = sc_mtx), colData=meta_data_arranged)


sce |> write_rds(here("intermediate_data/" , paste0("single_cell_experiment_", "RNA_", "counts_" ,dataset_name,".rds")  ))
bulk_mtx |> write_rds(here("intermediate_data" ,paste0( "bulk_mtx_", "RNA_", "counts_" ,dataset_name,".rds")  ))




