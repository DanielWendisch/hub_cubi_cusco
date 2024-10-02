# load conventionallz without BPcells
```{r eval=FALSE, include=FALSE}
dat <- Read10X("C:\\Users\\Danne\\raw_data\\cubi_cusco\\hub_01_outs\\multi\\count\\raw_feature_bc_matrix\\")
hub_01 <- CreateSeuratObject(counts = dat$`Gene Expression`, project = "hub_01")

hashtags <- CreateAssay5Object(counts=dat$`Multiplexing Capture`)
#check if column names are equal
if (all.equal(colnames(hashtags), colnames(hub_01))) {
  hub_01[["hashtags"]] <- hashtags
  hub_01 |> saveRDS(file = "intermediate_data/hub_01_raw_seurate.rds")
}else
  warning("cell names are not alike")

dat <- NULL
hashtags <- NULL
```


```{r}

hub_01 <- BPCells::open_matrix_10x_hdf5(
  path = "C:/Users/Danne/raw_data/cubi_cusco/hub_01_outs/multi/count/raw_feature_bc_matrix.h5",
  
)


hub_01 <- Read10X_h5(
  "C:/Users/Danne/raw_data/cubi_cusco/hub_01_outs/multi/count/raw_feature_bc_matrix.h5",use.names = TRUE
  
)
#separate the matrix in hashtag and gene reads
rows_with_CMO <- grepl("CMO", rownames(hub_01))
cmo_hub_01 <- hub_01[rows_with_CMO,]
hub_01 <- hub_01[!rows_with_CMO,]
hub_01 <- Azimuth:::ConvertEnsembleToSymbol(mat = hub_01,
                                            species = "human")

# Write the matrix to a directory
write_matrix_dir(
  mat = hub_01,
  dir = 'BPcell_matrices/hub_01')

write_matrix_dir(
  mat = cmo_hub_01,
  dir = 'BPcell_matrices/cmo_hub_01')
```

# setdiff(rownames(hub_01),
#         rownames(hub_01_cellbender)[!(grepl("ENSG", rownames(hub_01_cellbender)))]) |> length()
# 
# setdiff(rownames(hub_01_cellbender),
# rownames(hub_01)) |> length()
# 
# hub_01_cellbender <- scCustomize::Read_CellBender_h5_Mat("C:\\Users\\Danne\\raw_data\\cubi_cusco\\cellbender\\hub_01\\hub_01_cellbender_corrected_filtered.h5")
# 
# hub_01_cellbender <- Azimuth:::ConvertEnsembleToSymbol(mat = hub_01_cellbender, species = "human")



#| title: dings
#| echo: false
#| eval: false

filtered_seurat_obj <- seurat_obj |>
  filter(
    !(seurat_clusters %in% deleted_clusters) &
      nCount_RNA < upper_nCount_RNA_thresh &
      nCount_RNA > lower_nCount_RNA_thresh &
      percent_mito < upper_mito_thresh &
      sampletag_multiplets == "single_hashtag"
  )

filtered_seurat_obj |> DimPlot(label = T)

```




|                     | **lower**                              | **upper**                                    |
  |-----------------|------------------------|--------------------------------|
  | **nCount quantile** | `{r} lower_nCount_RNA_thresh_quantile` | `{r} upper_nCount_RNA_thresh_quantile`       |
  | **nCount RNA**      | `{r} lower_nCount_RNA_thresh`          | `{r} upper_nCount_RNA_thresh_non_scientific` |
  
  |                         | **lower** | **upper**                        |
  |-------------------------|-----------|----------------------------------|
  | **quantile**            | none      | `{r} upper_mito_thresh_quantile` |
  | **mitochondrial genes** | none      | `{r} upper_mito_thresh`          |
  
  |                      |                        |
  |----------------------|------------------------|
  | **deleted clusters** | `{r} deleted_clusters` |

