# script to read in tables with aligned reads from seasnap pipeline 
# and saves these in a combined tibble

library(purrr)
library(tidyverse)
library(vroom)
library(here)


path_raw_data_file <- "..\\..\\raw_data\\cubi_cusco"
path_raw_data_file <- here("../../raw_data/cubi_cusco")
path_input<- paste(sep = "\\", path_raw_data_file, "hub_01_02_bulk")

ensembl <- biomaRt::useMart("ensembl", dataset = "hsapiens_gene_ensembl")

#the following path contains folders with subfolders libname---all_lanes----all_mates----out in which the counts (all_mates.feature_counts) and the seqeuncing summary are stored
# store all count files of subfolders in object"all_files"
all_file_paths <- list.files(path="../../raw_data/cubi_cusco/hub_01_02_bulk", recursive = T, pattern = "all_mates.feature_counts", full.names = T)

# Step 2: Filter out files that contain "summary" in their name
all_file_paths <- grep("summary", all_file_paths, invert = TRUE, value = TRUE)

# Step 3: Use purrr::map to read the filtered files into a list
#file_list <- map(all_file_paths, ~ read_table(.x))

bulk_table_list <- map(all_file_paths, ~ vroom(.x, skip=1))

# Now file_list contains the data frames from each file
combined_table <- map(bulk_table_list, ~ .x[[7]]) |>
  bind_cols()

names(combined_table) <- sub("feature_counts\\.([^.]+)\\..*", "\\1", basename(all_file_paths))

# add ensemble id column 
combined_table$ensembl_gene_id <-  bulk_table_list[[1]] |> pull(Geneid )

# Add gene names

# Query to convert Ensembl IDs to gene names
gene_info <- biomaRt::getBM(
  attributes = c("ensembl_gene_id", "external_gene_name"),
  filters = "ensembl_gene_id",
  values = combined_table$ensembl_gene_id ,
  mart = ensembl
)

colnames(gene_info)[2] <- "gene"
combined_table <- as_tibble(combined_table)
combined_table <- left_join(combined_table,gene_info) |> relocate(gene,ensembl_gene_id)
combined_table |> write.csv(paste(path_raw_data_file,"hub_01_02_bulk\\counts_table_hub_01_02.csv", sep = "\\"))
