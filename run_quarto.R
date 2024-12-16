library(quarto)
library(tidyverse)
library(here)

setwd(here("scripts"))
dataset_params <- c("hub_01","hub_02")

for (param in dataset_params) {
  print(param)
  
  outpute_file_name <- paste0(param, "_processing_after_qc",".html")
  param_title <- paste0(param, "processing after QC")
  quarto_render(
    input = "hub_0x_processing_after_qc.qmd",
    output_file = outpute_file_name,
    execute_params = list(dataset_name=param, title=param_title))
  file.rename(outpute_file_name, paste0("../docs/",outpute_file_name))# move file to docs
  
}

setwd(here())



