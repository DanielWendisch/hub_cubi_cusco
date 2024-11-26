#this script calculates all additional meta data of dataset
dataset_name <- "hub_01"
#after QC

library(Seurat)
library(tidyseurat)
library(tidyverse)


path_raw_data_file <-  "..\\..\\raw_data\\cubi_cusco"

obj <- read_rds(paste0("intermediate_data/", dataset_name, "_integrated_", ".rds"))#sc data


obj <- obj |> mutate(m_e_class =case_when(
  harmony_clusters_res0.1 == 1 ~ "M",
  harmony_clusters_res0.1 == 2 ~ "M",
  harmony_clusters_res0.1 == 3 ~ "M",
  harmony_clusters_res0.1 == 4 ~ "M",
  harmony_clusters_res0.1 == 0 ~ "E"
  
))



# module score function # see also script hub_x_processing_after_qc.rmd

add_named_module_score <- function(seurat_object,named_gene_list){
  
  seurat_object <- AddModuleScore(seurat_object, features=named_gene_list, name = "m.s.")
  
  names_vec <-paste0("m.s.", names(named_gene_list))
  names(names_vec) <- paste0("m.s.", 1:length(names_vec))
  
  names(seurat_object@meta.data) <- str_replace_all(names(seurat_object@meta.data),names_vec)
  
  return(seurat_object)
}

yu_gramp_Cell_2021_tHIO <- read_excel("../../raw_data/cubi_cusco/genesets_from_literature/yu...gramp_Cell_2021_tHIO.xlsx")


vector_list <- yu_gramp_Cell_2021_tHIO %>%
  group_by(group) %>%
  summarise(value_vector = list(feature)) %>%
  deframe()
vector_list <- clean_names(vector_list)
names(vector_list) <-    paste0("t_hio__",names(vector_list)  )

#named_gene_list <- vector_list

obj<- add_named_module_score(obj,vector_list)


obj |>  write_rds(paste0("intermediate_data/", dataset_name, "_integrated_", ".rds"))#sc data