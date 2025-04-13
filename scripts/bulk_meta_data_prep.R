library(here)
library(tidyverse)
library(readxl)
library(janitor)
## hub_01

## meta data from 3p Sequencing, hub_01

# meta data from 3p Sequencing

raw_data_path <- here("../../raw_data/cubi_cusco")#raw_data sits in folder outside of project 
intermediate_data_path <- here("intermediate_data")

meta_data <- read_xlsx(paste0(raw_data_path,"/hub_01_02_bulk/a_HUB_model_characterisation_GEX_3pseq_transcription_profiling.xlsx"))

names(meta_data) <- meta_data[1,]
meta_data <- meta_data[-1,]
meta_data <- meta_data |> clean_names()

meta_data <- meta_data |> select(name_4,name,name_2,experiment_identifier, tissue_type,rna_amount_ng,tapestation_qc_rin,perform_date_4,differentiation_duration, condition, treatment, treatment_duration)

colnames(meta_data) <- c("sample","cell_line","name_2","experiment_identifier", "tissue_type","rna_amount_ng","tapestation_qc_rin","perform_date_4","differentiation_duration", "condition", "treatment", "treatment_duration"            )

hub_01_meta_data_3p <- meta_data |> filter(experiment_identifier=="HUB_01")
hub_01_02_meta_data_3p <- meta_data |> filter(experiment_identifier  %in% c("HUB_01", "HUB_02"))

## meta data from polyA Sequencing, hub_01


meta_data <- read_xlsx(paste0(raw_data_path,"/hub_01_02_bulk/a_HUB_model_characterisation_GEX_polyA_transcription_profiling.xlsx"))

names(meta_data) <- meta_data[1,]
meta_data <- meta_data[-1,]
meta_data <- meta_data |> clean_names()

meta_data <- meta_data |> select(name_4,name,name_2,experiment_identifier, tissue_type,rna_amount_ng,tapestation_qc_rin,perform_date_4,differentiation_duration, condition, treatment, treatment_duration)

colnames(meta_data) <- c("sample","cell_line","name_2","experiment_identifier", "tissue_type","rna_amount_ng","tapestation_qc_rin","perform_date_4","differentiation_duration", "condition", "treatment", "treatment_duration"            )

hub_01_meta_data_polyA <- meta_data |> filter(experiment_identifier=="HUB_01")
hub_01_02_meta_data_polyA <- meta_data |> filter(experiment_identifier  %in% c("HUB_01", "HUB_02"))

hub_01_meta_data <- bind_rows(hub_01_meta_data_polyA,hub_01_meta_data_3p)
hub_01_02_meta_data <- bind_rows(hub_01_02_meta_data_polyA,hub_01_02_meta_data_3p)
write_csv(hub_01_meta_data, here("intermediate_data","hub_01_meta_data.csv"))

write_csv(hub_01_02_meta_data, here("intermediate_data","hub_01_02_meta_data.csv"))



