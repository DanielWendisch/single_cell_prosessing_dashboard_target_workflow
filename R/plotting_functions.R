
generate_plot_grobs <- function(step_no, dataset = dataset_name, obj) {
    
  #library(Seurat) library(tidyverse) library(ggplotify) library(tidyseurat) library(ggExtra) library(ggridges) library(gridExtra) library(gt)
    


  ## remove enviroment links from ggplot
  # clean_plot <- function(plot) {
  #   plot$env <- new.env(parent = emptyenv())
  #   plot
  # }
  browser()
  DefaultAssay(obj) <- "RNA"
  
  plot_list_cellranger_cells <- list()
  ggextra_plot_list_cellranger_cells<- list()
  clean_plot_list_cellranger_cells<- list()
  
  
  ###############
  
  p1 <- obj |> 
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
  
  px <- obj |>
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
  
  p2 <- ggplot(obj, aes(x = log(rank(-nCount_RNA)), 
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
  
  plot_list_cellranger_cells[["vireo_calls_barchart_grid_cellbender_vireo"]] <- obj |>   
    select(cellbender_p95, vireo_snp_demux, cellranger) |> 
    group_by(cellbender_p95, vireo_snp_demux, cellranger) |> 
    count() |> 
    ggplot(aes(vireo_snp_demux,n, fill=vireo_snp_demux))+
    geom_col() +
    facet_grid(cellranger ~ cellbender_p95)+
    theme_bw()
  
  
  
  
  
  ###############
  
  p1b <-   obj |> 
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
  
  plot_list_cellranger_cells[["cellbender_prob_histogram"]] <- obj|> 
    ggplot( aes(x = log(nCount_RNA) , y=cellbender_p95, fill = cellbender_p95)) +
    geom_density_ridges(alpha=0.6, stat="binline")
  
  
  
  #########
  
  DefaultAssay(obj) <- "hashtag_oligos"
  
  plot_list_cellranger_cells[["hashtag_dimplot_vireo"]] <- obj |> DimPlot(group.by ="vireo_snp_demux" )
  plot_list_cellranger_cells[["hashtag_dimplot_cellranger"]] <- obj |> DimPlot(group.by ="cell_line_doublet_cellranger" )
  plot_list_cellranger_cells[["hashtag_dimplot_cellbender"]] <- obj |> DimPlot(group.by ="cell_line_doublet_cellranger" )
  
  
  DefaultAssay(obj) <- "RNA"
  
  plot_list_cellranger_cells[["RNA_dimplot_vireo"]] <- obj |> DimPlot(group.by ="vireo_snp_demux" )
  plot_list_cellranger_cells[["RNA_dimplot_cellranger"]] <- obj |> DimPlot(group.by ="cell_line_doublet_cellranger" )
  plot_list_cellranger_cells[["RNA_dimplot_cellbender"]] <- obj |> DimPlot(group.by ="cell_line_doublet_cellranger" )
  
  
  
  DefaultAssay(obj) <- "cellbender_RNA"
  
  plot_list_cellranger_cells[["cellbender_RNA_dimplot_vireo"]] <- obj |> DimPlot(group.by ="vireo_snp_demux" )
  plot_list_cellranger_cells[["cellbender_RNA_dimplot_cellranger"]] <- obj |> DimPlot(group.by ="cell_line_doublet_cellranger" )
  plot_list_cellranger_cells[["cellbender_RNA_dimplot_cellbender"]] <- obj |> DimPlot(group.by ="cell_line_doublet_cellranger" )
  
  
  
  
  ########
  
  
  text_list <- obj@misc$filtering_steps
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
  
  
  path_grob_filtered <- here(
    path_output_plots_root,
    paste(step_no, "plot_list_grob", dataset, "rds", sep = ".")
  )
  
  path_ggextra_filtered <- here(
    path_output_plots_root,
    paste(step_no, "plotified_plot_list_grob", dataset, "rds", sep = ".")
  )
  
  plot_list_grob |> write_rds(path_grob_filtered)
  #plotified_ggextra_plot_list_cellbender_cells |> write_rds(path_ggextra_filtered)
  plotified_ggextra_plot_list_cellranger_cells |> write_rds(path_ggextra_filtered)
  ## TODO remove side action
  return(path_grob_filtered)
  

}




make_dashboard_from_grobs <- function(step_no, dataset , path_grob_filtered) {
  path_grob_filtered <- here(
    path_output_plots_root,
    paste(step_no, "plot_list_grob", dataset, "rds", sep = ".")
  )
  

  path_ggextra_filtered <- here(
    path_output_plots_root,
    paste(step_no, "plotified_plot_list_grob", dataset, "rds", sep = ".")
  )
  
  quarto_title <- paste0(" \"title=\" ", paste(dataset_name, step_no),"  \" ")
  quarto_title <- paste0("title=", paste0(dataset_name,": ",  step_no))
  
  quarto_render_move(here("older workflow", "vireo.qmd"),
                     execute_params = list(
                       param_plot_list_grob = path_grob_filtered ,
                       param_ggextra_plot_list_grob = path_ggextra_filtered
                     ),
                     output_file = paste(dataset_name, step_no, "html", sep = "."),
                     output_dir = path_output_dashboards,
                     quarto_args = c("--metadata", quarto_title),
                     metadata= list(title=quarto_title)
  )
}


quarto_render_move <- function( #https://github.com/jhelvy/jph/blob/master/R/quarto_render_move.R
  input,
  output_file = NULL,
  output_dir = NULL,
  ...
) {
  
  # Get all the input / output file names and paths
  x <- quarto::quarto_inspect(input)
  output_format <- names(x$formats)
  output <- x$formats[[output_format]]$pandoc$`output-file`
  if (is.null(output_file)) { output_file <- output }
  input_dir <- dirname(input)
  if (is.null(output_dir)) { output_dir <- input_dir }
  output_path_from <- file.path(input_dir, output)
  output_path_to <- file.path(output_dir, output_file)
  
  # Render qmd file to input_dir
  quarto::quarto_render(input = input, ... = ...)
  
  # If output_dir is different from input_dir, copy the rendered output
  # there and delete the original file
  if (input_dir != output_dir) {
    
    # Try to make the folder if it doesn't yet exist
    if (!dir.exists(output_dir)) { dir.create(output_dir) }
    
    # Now move the output to the output_dir and remove the original output
    file.copy(
      from = output_path_from,
      to = output_path_to,
      overwrite = TRUE
    )
    file.remove(output_path_from)
    
    # If the output_dir is the same as input_dir, but the output_file
    # has a different name from the input file, then just rename it
  } else if (output_file != output) {
    file.rename(from = output_path_from, to = output_path_to)
  }
}



make_path_seurat_object <- function(step, dataset=dataset_name){
  path_seurat_obj<- here(path_intermediate_data_root,
                         paste("seurat_obj",
                               step,
                               dataset,
                               "rds",
                               sep = "."))
  return(path_seurat_obj)
  
}


