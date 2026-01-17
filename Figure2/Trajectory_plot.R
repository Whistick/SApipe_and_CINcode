############# Branch SlingShot ##########
df <- FetchData(adata,c("CIN01","BranchID"))
df_long <- df %>%
  mutate(BranchID_list = str_split(as.character(BranchID), ",")) %>%
  rowwise() %>%
  mutate(weight = 1 / length(BranchID_list)) %>%
  unnest(BranchID_list) %>%
  select(-BranchID) %>% 
  rename(BranchID = BranchID_list) %>%
  mutate(BranchID = as.character(BranchID)) %>%
  ungroup()
plot_data <- df_long %>%
  group_by(CIN01, BranchID) %>%
  summarise(freq = sum(weight), .groups = "drop")
plot_data <- plot_data %>%
  group_by(CIN01) %>%
  mutate(prop = freq / sum(freq)) %>%
  ungroup() %>%
  rename(CIN = CIN01, Celltype = BranchID)
colors <- c("1" = "#1b9e77", "2" = "#d95f02", "3" = "#7570b3")
ggplot(plot_data,
       aes(x = CIN, stratum = Celltype, alluvium = Celltype, y = prop)) +
  geom_stratum(aes(fill = Celltype), width = 0.6, color = "grey60") +
  geom_flow(aes(fill = Celltype), width = 0.3, alpha = 0.4) +
  scale_fill_manual(values = colors) +
  scale_x_discrete(expand = c(0, 0.15)) + 
  scale_y_continuous(labels = scales::percent_format(), expand = c(0, 0)) +
  theme_few() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(x = "", y = "", fill = "Trajectory")
