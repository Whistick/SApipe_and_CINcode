
adata$Sub_CellType <- factor(
  adata$Sub_CellType,
  levels = c("Basal_Cells", "Sq_01", "Sq_02", "Mal_01", "Mal_02","Mal_03" )
)
my_colors <- c(
  "Basal_Cells" = "#FDBF6F",
  "Sq_01" = "#A6CEE3",
  "Sq_02" = "#1F78B4",
  "Mal_01" = "#B2DF8A",
  "Mal_02" = "#33A02C",
  "Mal_03" = "#FF7F00"
)

plots <- DynamicPlot(
  srt = adata, 
  lineages = paste0("Lineage", 1:3), 
  group.by = "Sub_CellType",
  features = c("CDKN2A","MKI67"),
  compare_lineages = FALSE, 
  compare_features = TRUE,
  ncol = 3,
  combine = FALSE
)

y_min <- min(sapply(plots, function(p) min(ggplot_build(p)$data[[1]]$y, na.rm=TRUE)))
y_max <- max(sapply(plots, function(p) max(ggplot_build(p)$data[[1]]$y, na.rm=TRUE)))

plots_fixed <- lapply(plots, function(p) {
  p + coord_cartesian(ylim=c(0, 1)) + scale_color_manual(values = my_colors) + scale_fill_manual(values = my_colors)
})

library(patchwork)
combined <- wrap_plots(plots_fixed, ncol=3)
combined