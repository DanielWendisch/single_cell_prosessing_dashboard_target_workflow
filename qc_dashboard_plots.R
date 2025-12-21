#makes the following lists of plots objects an saves them as grobs and ggplots with linekd data (.rds)\

plots_list <- list()
ggextra_plot_list <- list()

assay_iter <- "RNA"
cluster.name = paste(assay_iter, "clusters", leiden_or_louvain, "res." , resolution_clustering_1, sep = "_") # calculated in basic_seurat_processing_for...R
basic_cluster_name <- cluster.name
#######################set-up
library(Seurat)
library(tidyverse)
library(ggplotify)
library(ggridges)

my_mad <- function(x, constant = 1) {
  median(abs(x - median(x, na.rm = TRUE)), na.rm = TRUE) * constant
}


theme_1 <-   theme( plot.title = element_text(hjust = 0.5),
                    plot.subtitle = element_text(hjust = 0.5))


######################
# TODO why are these not in prespecified tibble in plot_specs...R





# ############
# 
# for (i in 1:nrow(plot_specs)) {
#   DefaultAssay(seurat_obj) <- plot_specs$default_assay[i]
#   
#   
#   plots_list[[plot_specs$list_int[i]]] <- seurat_obj |>
#     DimPlot(group.by = plot_specs$meta_data_column[i], order = plot_specs$order[[i]][1])  + 
#     theme_1+
#     ggtitle(label=plot_specs$plote_title[i], subtitle = plot_specs$subtitle[i])
#   
# }


# make dim plot according to plot_specs
plots_list <- list()
assays <- unique(plot_specs$default_assay)

for (j in seq_along(assays)) {
  message(assays[j])
  
  
  plot_specs_iter <- plot_specs |> filter(default_assay==assays[j])
  DefaultAssay(seurat_obj) <- assays[j]
  list_iter <- list()
  for (i in seq_along(unique(plot_specs_iter$title_short))){
    
    message(plot_specs_iter$title_short[[i]])
    
    
    
    list_iter[[plot_specs_iter$title_short[i]]] <- seurat_obj |>
      
    DimPlot(group.by = plot_specs_iter$meta_data_column[i], order = plot_specs_iter$order[[i]][1])  + 
    theme_1+
    ggtitle(label=plot_specs_iter$plote_title[i], subtitle = plot_specs_iter$subtitle[i])
    }
  plots_list[[assays[j]]] <- list_iter
  
  }


DefaultAssay(seurat_obj) <- "RNA"

plots_list$RNA[["cellbender_cell_prob_embedding_RNA"]] <- seurat_obj |> 
  FeaturePlot("cellbender_prob_to_be_cell", min.cutoff = 0.98, cols = c("red","black"))+
  theme_1+
  ggtitle(label="probability to be cell (cellbender)",
          subtitle = "background corrected counts)")


plots_list$RNA[["cellbender_cell_prob_embedding_RNA"]]

DefaultAssay(seurat_obj) <- "cellbender_RNA"
plots_list$RNA[["cellbender_cell_prob_embedding_cellbender_RNA"]] <- seurat_obj |> 
  FeaturePlot("cellbender_prob_to_be_cell", min.cutoff = 0.98, cols = c("red","black"))+
  theme_1+
  ggtitle(label="probability to be cell (cellbender)",
          subtitle = "background corrected counts (cellbender)")


DefaultAssay(seurat_obj) <- "RNA"




plots_list$RNA[[clustering_rough_RNA]] <- seurat_obj |> DimPlot(group.by = clustering_rough_RNA)


plots_list$RNA[["bar_cell_line_doublet_cellranger"]] <- seurat_obj |> ggplot(aes(.data[[basic_cluster_name]],fill=cell_line_doublet_cellranger))+
geom_bar()+
theme_minimal()


DefaultAssay(seurat_obj) <- "cellbender_RNA"

plots_list$cellbender_RNA[[clustering_rough_cellbender]] <- seurat_obj |> DimPlot(group.by = clustering_rough_cellbender)


plots_list$cellbender[["bar_cell_line_doublet_cellranger"]] <- seurat_obj |> ggplot(aes(.data[[clustering_rough_RNA]],fill=cell_line_doublet_cellranger))+
  geom_bar()+
  theme_minimal()

#############################################################################################
#############################################################################################
# QC metrics, after background correction




plots_list$RNA[["percent_mito_RNA_embedding"]] <- seurat_obj |> ggplot(aes(umapRNA_1, umapRNA_2, color=percent_mito)) + geom_point() +theme_1
plots_list$RNA[["mito_group_RNA_embedding"]] <- seurat_obj |> ggplot(aes(umapRNA_1, umapRNA_2, color=mito_group)) + geom_point() +theme_1



DefaultAssay(seurat_obj) <- "RNA"
plots_list$RNA[["perc_mito_RNA_embedding"]] <-seurat_obj |>
  FeaturePlot("percent_mito", reduction = "umap.RNA", max.cutoff= 10) +
  ggtitle("percent mitochondrial reads, raw RNA counts.")


DefaultAssay(seurat_obj) <- "cellbender_RNA"
plots_list$cellbender_RNA[["perc_mito_RNA_embedding"]] <- seurat_obj |> FeaturePlot("percent_mito", reduction = "umap.cellbender_RNA", max.cutoff= 10) + ggtitle("percent mitochondrial reads, with backgr. correc.")





plots_list$RNA[["binned_nCount_RNA_singlet_doublet_cellranger"]] <- seurat_obj |> 
  mutate(quantile_nCount_RNA = ntile(nCount_RNA,100)) |>  # Create groups in steps of 1000
  group_by(quantile_nCount_RNA,singlet_doublet_cellranger) |> 
  summarise(n = n(), .groups = 'drop', min_of_bin_nCount_RNA=min(nCount_RNA)) |> 
  group_by(quantile_nCount_RNA) |> 
  mutate(freq=n/sum(n)) |> 
  ggplot(aes(min_of_bin_nCount_RNA,freq,color=singlet_doublet_cellranger)) +
  geom_point()+
  geom_smooth()+
  geom_hline(yintercept = 0.1, color = "red", linetype = "dashed", size = 1)+
  geom_hline(yintercept = 0.9, color = "red", linetype = "dashed", size = 1)+
  geom_vline(xintercept = lower_nCount_RNA_thresh, color = "black", linetype = "dashed", size = 1)+
  geom_vline(xintercept = upper_nCount_RNA_thresh, color = "black", linetype = "dashed", size = 1)+
  scale_x_continuous(breaks = seq(0, max(seurat_obj$nCount_RNA, na.rm = TRUE), by = 2000))+
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


DefaultAssay(seurat_obj) <- "RNA"
# TODO FIX
# plots_list$RNA[["Dimplot_percent_mito"]] <- seurat_obj |> VlnPlot(features = "percent_mito")


DefaultAssay(seurat_obj) <- "RNA"
plots_list$RNA[["Dimplot_nCount_RNA"]] <- seurat_obj |> VlnPlot("nCount_RNA")

################ Mean absolute deviation of perc mito
#Do I nbeed this??
mt_mad_tbl <- tibble(percent_mito = seurat_obj$percent_mito,
                     seurat_clusters= seurat_obj$seurat_clusters,
                     cell_line_doublet_cellranger= seurat_obj$cell_line_doublet_cellranger)


mt_mad_tbl <- mt_mad_tbl |>
  group_by(cell_line_doublet_cellranger) |>
  mutate(MAD_percent_mito =my_mad(x=percent_mito)) |>
  mutate(dead_by_mad=ifelse(MAD_percent_mito>upper_mito_thresh, "dead","alive"))



p <- mt_mad_tbl |> ggplot(aes(percent_mito, cell_line_doublet_cellranger, fill = dead_by_mad)) +
  geom_density_ridges()+ theme(axis.text.x = element_text(angle = 45, hjust=1)) +
  xlim(0,20)+
  theme_minimal()


##################

plots_list$RNA[["nCount_RNA_histo"]] <- seurat_obj |> ggplot(aes(nCount_RNA,  fill = cell_line_doublet_cellranger))  +
  geom_histogram(bins = 200)+
  theme(axis.text.x = element_text(angle = 45, hjust=1))+
  facet_wrap(~cell_line_doublet_cellranger,ncol = 1, scales="free" )+
  scale_x_continuous(breaks = seq(0, max(seurat_obj$nCount_RNA, na.rm = TRUE), by = 5000))+
  xlim(0,70000) +
  theme_minimal()




#make grob of all plots in sublists
dim_plot_grob_list <- list()
for (sub_list_name in names(plots_list)) {
  print(sub_list_name)
  dim_plot_grob_list[[sub_list_name]] <-map(plots_list[[sub_list_name]], ~ ggplotGrob((.x)))
}





####################
## ggExtra has to be converted to grob with library(ggplotify)


pmultiplet_n_count_RNA_n_Feature<- seurat_obj |>
  ggplot(aes(nCount_RNA, nFeature_RNA, color=singlet_doublet_cellranger)) +  
  geom_point(size=0.5) +
  theme_bw() +
  #xlim(c(-100,70000)) +
  #ylim(c(-1000,8000)) +
  #scale_color_manual(values = c( "grey","red","green"))  +  
  geom_vline(aes(xintercept = lower_nCount_RNA_thresh),linetype = "dashed") +
  geom_vline(aes(xintercept = upper_nCount_RNA_thresh),linetype = "dashed")+
  guides(colour = guide_legend(override.aes = list(size=7)))+
  scale_x_continuous(breaks = seq(0, 150000, by = 5000))+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))



ggextra_plot_list$nFeature_vs_nCount_RNA <- ggExtra::ggMarginal(pmultiplet_n_count_RNA_n_Feature +
                                                                      theme(legend.position = "bottom") ,
                                                                    type = "density",
                                                                    groupColour = TRUE,
                                                                    groupFill = TRUE,
                                                                    size = 3) 






ggextra_plot_list_grob <- as.ggplot(ggextra_plot_list$nFeature_vs_nCount_RNA) |> ggplotGrob()


dim_plot_grob_list$RNA$ggextra_plots <- ggextra_plot_list_grob 





