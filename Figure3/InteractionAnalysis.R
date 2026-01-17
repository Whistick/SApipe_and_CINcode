##################### LIANA ##########################
testdata  <- qread(qs_path)
liana_path <- system.file(package = "liana")
testdata <- NormalizeData(testdata, normalization.method = "LogNormalize")  # generate data slot 
# testdata$`cell_labels` <- testdata$`clusters` 
testdata$`cell_labels` <- testdata$`Sub_CellType`  

# filter
testdata <- testdata[,testdata$`cell_labels` %in% lbl_all]
testdata <- SetIdent(testdata, value = testdata$cell_labels) # set ident
testdata$cell_labels <- droplevels(testdata$cell_labels)
table(testdata$cell_labels)

# Run LIANA
liana_test <- liana_wrap(testdata)
liana_test <- liana_test %>% liana_aggregate()
liana_test %>% dplyr::glimpse()


#### Visualize ####
# Heatmap Overview
liana_trunc <- liana_test %>%
  filter(aggregate_rank <= 0.01) # note that these pvals are already corrected
heat_freq(liana_trunc)

mph <- c("Mph_01_MMP12",'Mph_02_PDGFC','Mph_03_KCNQ3','Mph_04_VCAN')
epi <- c('Malignant_Atypical_hyperplasia','Malignant_Differationated','Squamous_Differationated','Malignant_ImmuneResponse','Basal_Cells','Squamous_Proliferation')
liana_dotplot(liana_test,source_groups = epi,target_groups = mph,ntop = 70) + theme(
  axis.text.x = element_text(size = 12,angle=45, hjust = 1,vjust=1), 
  axis.text.y = element_text(size = 10),                   
  plot.title = element_text(size = 18),       
  legend.text = element_text(size = 12),             
  strip.text = element_text(size = 12),                
  strip.text.x = element_text(size = 14)
)

