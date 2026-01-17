# 1. read data
metadata <- read.csv("metadata.csv")
metadata$CIN01 <- factor(metadata$CIN01, levels = c("NOR", "LSIL", "HSIL", "SCC"))

cell_map <- c(
  "Basal_Cell" = "Epithelial", "Columnar" = "Epithelial", "Glandular" = "Epithelial", "Squamous" = "Epithelial",
  "Malignants" = "Malignant",
  "B_Cell" = "Immune", "CD4T_Cell" = "Immune", "CD8T_Cell" = "Immune", "Dendritic_Cell" = "Immune", 
  "Marcophage" = "Immune", "Mast_Cell" = "Immune", "Neutrophil" = "Immune", "NK_Cell" = "Immune", "Plasma_Cell" = "Immune",
  "Endothelial" = "Stromal", "Fibroblast" = "Stromal", "Schwann_Cell" = "Stromal", "SMCs" = "Stromal",
  "Doublet" = "Unassigned", "low-quality" = "Unassigned"
)
metadata$type <- cell_map[metadata$CellType]
table(metadata$type)
df <- data.frame(
  CIN = metadata$CIN01,
  Celltype = metadata$type
)
df <- df[df$Celltype %in% c('Epithelial','Malignant','Stromal','Immune'),]
df$Celltype <- factor(df$Celltype ,levels = c('Immune','Stromal','Epithelial','Malignant'))
colors <- c(
  "Immune" = "#FCCDE5",
  "Stromal" = "#A6CEE2",
  "Epithelial" = "#dadaeb",
  "Malignant" = "#FB8072"
)

plot_data <- df %>%
  group_by(CIN, Celltype) %>%
  summarise(count = n(), .groups = "drop") %>%
  group_by(CIN) %>%
  mutate(freq = count / sum(count)) %>%
  ungroup()
ggplot(plot_data,
       aes(x = CIN, stratum = Celltype, alluvium = Celltype, y = freq)) +
  geom_stratum(aes(fill = Celltype), width = 0.6, color = "grey60") +
  geom_flow(aes(fill = Celltype), width = 0.3, alpha = 0.4) +
  scale_fill_manual(values = colors) +
  scale_x_discrete(expand = c(0, 0.15)) + 
  scale_y_continuous(labels = scales::percent_format(),expand = c(0, 0)) +
  theme_few() +
  labs(x = "", y = "", fill = "Cell Type") +
  theme(
    panel.background = element_blank(),
    plot.background = element_blank(),
    axis.line = element_line(color = "black", linewidth = 0.2),
    axis.ticks = element_line(color = "black", linewidth = 0.2),
    axis.ticks.length = unit(0.25, "cm"), 
    axis.title.x = element_text(size = 10, color = "black"),
    axis.title.y = element_text(size = 10, color = "black"),
    axis.text.x = element_text(color = "black", size = 8, angle = 45, hjust = 1),
    axis.text.y = element_text(color = "black", size = 8),
    plot.title = element_text(size = 10, color = "black", hjust = 0.5),
    panel.grid = element_blank(),
    # panel.border = element_rect(color = "black", fill = NA, linewidth = 0.2),
    panel.border = element_blank(),
    # legend.position = "none",
    # aspect.ratio = 1.5,
  )