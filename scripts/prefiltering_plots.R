library(Seurat)
library(tidyverse)
library(ggplotify)
library(tidyseurat)
library(ggExtra)
library(ggridges)
library(gridExtra)
library(gt)

## remove enviroment links from ggplot
# clean_plot <- function(plot) {
#   plot$env <- new.env(parent = emptyenv())
#   plot
# }

DefaultAssay(seurat_obj) <- "RNA"

plot_list_cellranger_cells <- list()
ggextra_plot_list_cellranger_cells<- list()
clean_plot_list_cellranger_cells<- list()


###############

p1 <- seurat_obj |> 
  ggplot(aes(x = log(rank(-nCount_RNA)),
             y = log(nCount_RNA),
             color = vireo_snp_demux)) +
  geom_point() +
  theme_classic()


p1_marginal <- ggExtra::ggMarginal(p1 +
                            theme(legend.position = "bottom") ,
                          type = "histogram",
                          groupColour = TRUE,
                          groupFill = TRUE,
                          size = 3,
                          margins = "y")


plot_list_cellranger_cells[["nCount_RNA_rank__vs__RNA_rank"]] <- p1
ggextra_plot_list_cellranger_cells[["nCount_RNA_rank__vs__RNA_rank"]] <- as.ggplot(p1_marginal)




############

px <- seurat_obj |>
  filter(!is.na(vireo_snp_demux)) |> 
  ggplot(aes(log(rank(-nCount_RNA)),
             nCount_RNA,color=vireo)) +
  geom_point() +
  # scale_x_log10(labels = scales::comma_format()) +
  scale_y_log10() +
  expand_limits(y = 0) +
  ggplot2::annotation_logticks() +
  ggplot2::theme_classic() +
  ggtitle("without vireo.NA")

px_marginal <- ggExtra::ggMarginal(px +
                      theme(legend.position = "bottom") ,
                    type = "histogram",
                    groupColour = TRUE,
                    groupFill = TRUE,
                    size = 3,
                    margins = "y")

plot_list_cellranger_cells[["nCount_RNA_rank__vs__RNA_rank_no_NA"]] <- px
ggextra_plot_list_cellranger_cells[["nCount_RNA_rank__vs__RNA_rank_no_NA"]] <- as.ggplot(px_marginal) #|> ggplotGrob()



#################

p2 <- ggplot(seurat_obj, aes(x = log(rank(-nCount_RNA)), 
                             y = log(nCount_RNA),
                             color = cellbender_p95)) +
  geom_point()+
  theme_classic()


p2_marginal <- ggExtra::ggMarginal(p2+
                            theme(legend.position = "bottom") ,
                          type = "histogram",
                          groupColour = TRUE,
                          groupFill = TRUE,
                          size = 3,
                          margins = "y")



plot_list_cellranger_cells[["nCount_RNA_rank__vs__RNA_rank_color_cellbender"]] <- p2
ggextra_plot_list_cellranger_cells[["nCount_RNA_rank__vs__RNA_rank_color_cellbender"]] <- as.ggplot(p2_marginal)




################

plot_list_cellranger_cells[["vireo_calls_barchart_grid_cellbender_vireo"]] <- seurat_obj |>   
  select(cellbender_p95, vireo_snp_demux, cellranger) |> 
  group_by(cellbender_p95, vireo_snp_demux, cellranger) |> 
  count() |> 
  ggplot(aes(vireo_snp_demux,n, fill=vireo_snp_demux))+
  geom_col() +
  facet_grid(cellranger ~ cellbender_p95)+
  theme_bw()





###############

p1b <-   seurat_obj |> 
  ggplot( aes(x = log(rank(-nCount_RNA)), 
              y= cellbender_prob_to_be_cell,
              color = cellbender_p95)) +
  geom_point() +
  theme_classic()

p1b_marginal<- ggExtra::ggMarginal(p1b +
                      theme(legend.position = "bottom") ,
                    type = "histogram",
                    groupColour = TRUE,
                    groupFill = TRUE,
                    size = 3,
                    margins = "y")

plot_list_cellranger_cells[["cellbender_prob__vs__RNA_rank_color_cellbender"]] <- p1b
ggextra_plot_list_cellranger_cells[["cellbender_prob__vs__RNA_rank_color_cellbender"]] <- as.ggplot(p1b_marginal)


#########

plot_list_cellranger_cells[["cellbender_prob_histogram"]] <- seurat_obj|> 
  ggplot( aes(x = log(nCount_RNA) , y=cellbender_p95, fill = cellbender_p95)) +
  geom_density_ridges(alpha=0.6, stat="binline")



#########

DefaultAssay(seurat_obj) <- "hashtag_oligos"

plot_list_cellranger_cells[["hashtag_dimplot_vireo"]] <- seurat_obj |> DimPlot(group.by ="vireo_snp_demux" )
plot_list_cellranger_cells[["hashtag_dimplot_cellranger"]] <- seurat_obj |> DimPlot(group.by ="cell_line_doublet_cellranger" )
plot_list_cellranger_cells[["hashtag_dimplot_cellbender"]] <- seurat_obj |> DimPlot(group.by ="cell_line_doublet_cellranger" )


DefaultAssay(seurat_obj) <- "RNA"

plot_list_cellranger_cells[["RNA_dimplot_vireo"]] <- seurat_obj |> DimPlot(group.by ="vireo_snp_demux" )
plot_list_cellranger_cells[["RNA_dimplot_cellranger"]] <- seurat_obj |> DimPlot(group.by ="cell_line_doublet_cellranger" )
plot_list_cellranger_cells[["RNA_dimplot_cellbender"]] <- seurat_obj |> DimPlot(group.by ="cell_line_doublet_cellranger" )



DefaultAssay(seurat_obj) <- "cellbender_RNA"

plot_list_cellranger_cells[["cellbender_RNA_dimplot_vireo"]] <- seurat_obj |> DimPlot(group.by ="vireo_snp_demux" )
plot_list_cellranger_cells[["cellbender_RNA_dimplot_cellranger"]] <- seurat_obj |> DimPlot(group.by ="cell_line_doublet_cellranger" )
plot_list_cellranger_cells[["cellbender_RNA_dimplot_cellbender"]] <- seurat_obj |> DimPlot(group.by ="cell_line_doublet_cellranger" )




########


text_list <- seurat_obj@misc$filtering_steps
# Combine into one long-format table
combined_table <- imap_dfr(text_list, function(vec, name) {
  tibble(value = vec, group = name)
}, .id = "row_id") |> 
  select(-row_id)

# Create a gt table with row groups
gt_combined <- combined_table |> 
  gt(groupname_col = "group") |> 
  cols_label(value = "") |> 
  tab_options(
    row_group.font.weight = "bold",
    table.font.size = px(14)
  )



# Optionally, display them together
#grid.arrange(grobs = plot_list,ncol=1 )
#filter_annotations <- arrangeGrob(grobs = gt_combined, ncol = 1)






#############
#normal ggplots

plot_list_grob <- map(plot_list_cellranger_cells, ~ ggplotGrob((.x)))



###
#remove enviornemt liks from ggplot
#clean_plot_list_grob <- map(plot_list_cellranger_cells, ~ clean_plot((.x)))
#clean_plot_list_grob |> write_rds(here("output", "hub_10","clean_plot_list.rds"))

####
# ggextra as.ggplots
plotified_ggextra_plot_list_cellranger_cells <- map(ggextra_plot_list_cellranger_cells, ~ as.grob((.x)))



plot_list_grob$text_filter_annotations <- gt_combined















