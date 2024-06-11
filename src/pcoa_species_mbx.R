##################################################
#R program for creating PCoA plots
#Extended Data Figure 1A
#Extended Data Figure 2A
##################################################

library(dplyr)
library(tidyverse)
library(viridis)
library(ggplot2)
library(vegan)
library(gridExtra)
library(RColorBrewer)

setwd("~/b2b")
unfiltered_species <- read.delim('input/species_v4.tsv',row.names=1) %>% t() %>% as.data.frame() %>% rownames_to_column() %>% rename(barcode_metagenomics = rowname) 
final_metadata <- read.delim('input/meta_df.tsv',row.names=1)
final_metadata <- final_metadata[!duplicated(final_metadata$barcode_metagenomics),]
df_w_meta <- left_join(unfiltered_species,final_metadata,by="barcode_metagenomics") 

nafld_data<-df_w_meta %>% filter(cohort=="NAFLD" | case==0)

nafld_data_all <- nafld_data 
nafld_data_meta <- nafld_data %>% select(!starts_with('s__')) %>% rownames_to_column("sample_id") 
nafld_data_species <- nafld_data %>% select(starts_with('s__')) 
nafld_data_species_id <- nafld_data %>% select("barcode_metagenomics",starts_with('s__')) 

df.min <- (min(nafld_data_species[nafld_data_species > 0])/2)
nafld_data_species[nafld_data_species == 0] <- df.min 

#Read in taxonomic data 
tax<-read_tsv("metaphlan_taxonomic_profiles.tsv")
names(tax) = gsub("_taxonomic_profile", "", names(tax))
bactfirm.data<-tax %>% filter(`# taxonomy`=="k__Bacteria|p__Bacteroidetes" | `# taxonomy`=="k__Bacteria|p__Firmicutes") %>% select(-`# taxonomy`) %>% t() %>% as.data.frame() %>% rownames_to_column() %>% rename(firm.abund=V1,bact.abund=V2,barcode_metagenomics=rowname)

nafld_data_species_w_bactfirm <- left_join(nafld_data_species_id,bactfirm.data,by="barcode_metagenomics") 
nafld_data_species_w_bactfirm_input <- nafld_data_species_w_bactfirm %>% select(starts_with('s__')) 
nafld_data_species_w_bactfirm_meta <- nafld_data_species_w_bactfirm %>% select(-c(starts_with('s__'))) %>% rownames_to_column("sample_id")

bray <- nafld_data_species_w_bactfirm_input %>% vegdist(., "bray") 
pc = capscale(bray~1, comm = nafld_data_species_w_bactfirm_input)
pc.summary<-summary(pc)

pcl.bray <- as.data.frame(pc$CA$u) %>%
  select(MDS1, MDS2)

pcl.bray <- pcl.bray %>%
  rownames_to_column("sample_id") %>%
  inner_join(nafld_data_species_w_bactfirm_meta, by = "sample_id") %>%
  select(colnames(nafld_data_species_w_bactfirm_meta), everything()) 

# save axes r2
pco1.r2 <- paste("PCo1 (", round(pc.summary$cont$importance[2,1]*100, digits = 1), "%)", sep = '')
pco2.r2 <- paste("PCo2 (", round(pc.summary$cont$importance[2,2]*100, digits = 1), "%)", sep = '')
rm(pc, pc.summary)
#Bacteroidetes
ggplot(pcl.bray, 
       aes(MDS1, MDS2)) + 
  geom_point(aes(color=bact.abund/100), alpha = 0.7, size = 2, stroke = 1) + 
  scale_color_viridis_c(option = "E",limits=c(0,1),breaks=seq(0,1,0.2)) +
  coord_fixed() + 
  theme_bw(base_size=24) + 
  ggtitle("Bacteroidetes") + 
  labs(x = pco1.r2, 
       y = pco2.r2,
       fill = "") +
  labs(color = "Relative abundance") 
ggsave(
  file.path("output.mp4", "pco_bacteroidetes.pdf"),
  dpi = 300, width=10, height=6
)
#Firmicutes
ggplot(pcl.bray, 
       aes(MDS1, MDS2)) + 
  geom_point(aes(color=firm.abund/100), alpha = 0.7, size = 2, stroke = 1) + 
  scale_color_viridis_c(option = "E",limits=c(0,1),breaks=seq(0,1,0.2)) +
  coord_fixed() + 
  theme_bw(base_size=24) + 
  ggtitle("Firmicutes") + 
  labs(x = pco1.r2, 
       y = pco2.r2,
       fill = "") +
  labs(color = "Relative abundance") 
ggsave(
  file.path("output.mp4", "pco_firmicutes.pdf"),
  dpi = 300, width=10, height=6
)

###Extended Data Figure 1A
###Regular PCoA -- NAFLD                                         
bray <- nafld_data_species %>% vegdist(., "bray") 
pc = capscale(bray~1, comm = nafld_data_species)
pc.summary<-summary(pc)

pcl.bray <- as.data.frame(pc$CA$u) %>%
  select(MDS1, MDS2)

pcl.bray <- pcl.bray %>%
  rownames_to_column("sample_id") %>%
  inner_join(nafld_data_meta, by = "sample_id") %>%
  select(colnames(nafld_data_meta), everything()) 


# save axes r2
pco1.r2 <- paste("PCo1 (", round(pc.summary$cont$importance[2,1]*100, digits = 1), "%)", sep = '')
pco2.r2 <- paste("PCo2 (", round(pc.summary$cont$importance[2,2]*100, digits = 1), "%)", sep = '')
rm(pc, pc.summary)

ggplot(pcl.bray, 
       aes(MDS1, MDS2)) + 
  geom_point(aes(color = as.factor(case)), alpha = 0.7, size = 2, stroke = 1) + 
  scale_color_manual(values = c("0" = "#999999", "1" = "#E69F00")) +
  coord_fixed() + 
  theme_bw(base_size=24) + 
  ggtitle("NAFLD") + 
  labs(x = pco1.r2, 
       y = pco2.r2,
       fill = "") +
  labs(color = "NAFLD case") 

ggsave(
  file.path("output.mp4", "pco.pdf"),
  dpi = 300, width=10, height=6
  )

###Metabolites
unfiltered_mbx_old <- read.delim("annotated_metabolites.tsv",row.names=1) %>% t() %>% as.data.frame() %>% rownames_to_column() %>% rename(barcode_metabolomics = rowname) 
unfiltered_mbx <- read.delim("annotated_metabolites_w_methods.tsv",row.names=1) %>% select(starts_with("X")) %>% t() %>% as.data.frame() %>% rownames_to_column() %>% rename(barcode_metabolomics = rowname) 

###Normalize the mbx by methods. Divide TF by median of the ratio and pick the ones that is the most abundant among columns (hilic, c18)
#first, check the median of the ratio: (TF+1)/(QI+1)
unfiltered_mbx_test <- read.delim("annotated_metabolites_w_methods.tsv",row.names=1) %>% select(starts_with("X")) %>% rownames_to_column("metabolite")
df <- separate(unfiltered_mbx_test, metabolite, into = c("column1", "column2", "column3"), sep = "_")
df <- df[duplicated(df$column1) | duplicated(df$column1, fromLast = TRUE), ]
filtered_df <- subset(df, column1 %in% unique(column1[column3 == "TF" & column1 %in% column1[column3 == "QI"]]))
#26 metabolites measured by both TF and QI. Glycocholic acid has 3
filtered_df <- filtered_df %>%
  filter(!(column1 == "Glycocholic acid" & column2 == "HILIC-neg"))
grouped_data <- split(filtered_df, filtered_df$column1)
# Function to calculate the ratios of TF to QI for each group
calculate_ratios <- function(group) {
  tf_rows <- group[group$column3 == "TF", ]
  qi_rows <- group[group$column3 == "QI", ]
  
  ratios <- (tf_rows[, 4:ncol(group)]+1) / (qi_rows[, 4:ncol(group)]+1)
  
  return(ratios)
}
# Create an empty dataframe to store the ratios
ratios_df <- data.frame()

for (metabolite in names(grouped_data)) {
  group <- grouped_data[[metabolite]]
  
  # Calculate the ratios for the current metabolite
  ratios <- calculate_ratios(group)
  
  # Add the ratios to the dataframe
  ratios_df <- rbind(ratios_df, ratios)
}
values <- unlist(ratios_df)
summary(values) #median is 41

# Divide columns that end with _TF by 41
unfiltered_mbx[, grep("_TF$", names(unfiltered_mbx))] <- unfiltered_mbx[, grep("_TF$", names(unfiltered_mbx))] / 41
# Create a list of unique prefixes (metabolites) in the column names
prefixes <- unique(sapply(names(unfiltered_mbx[-1]), function(x) strsplit(x, "_")[[1]][1]))
df_final <- data.frame(matrix(ncol = 0, nrow = nrow(unfiltered_mbx)))
# Loop through each prefix
for (p in prefixes) {
  # Find all columns that start with the prefix and store their names in a vector
  cols <- grep(p, names(unfiltered_mbx), value = TRUE, fixed = TRUE)
  
  # If there is only one column with this prefix, keep it as is
  if (length(cols) == 1) {
    df_final <- cbind(df_final, unfiltered_mbx[, cols, drop = FALSE])
  } else if (length(cols) > 1) {
    # If there are multiple columns with this prefix, find the column with the highest value
    max_col <- apply(unfiltered_mbx[, cols], 1, max) %>% as.data.frame() 
    names(max_col)[names(max_col) == "."] <- p
    df_final <- cbind(df_final, max_col)
  }
}
# Rename the columns to remove the prefix and suffix
names(df_final) <- sapply(names(df_final), function(x) strsplit(x, "_")[[1]][1])
rownames(df_final) <- unfiltered_mbx$barcode_metabolomics
df_final_join <- df_final %>% rownames_to_column("barcode_metabolomics")

final_metadata <- read.delim('input/meta_df.tsv',row.names=1)
final_metadata <- final_metadata[!duplicated(final_metadata$barcode_metabolomics),]
df_w_meta <- left_join(df_final_join,final_metadata,by="barcode_metabolomics") 
mbx_list<-names(df_final)
mbx.data<-df_w_meta %>% select(all_of(mbx_list))

nafld.data<-df_w_meta %>% filter(cohort=="NAFLD") 
mbx.nafld.data<-nafld.data %>% select(all_of(mbx_list)) 

nafld_data<-df_w_meta %>% filter(cohort=="NAFLD" | case==0)
nafld_data_all <- nafld_data 
nafld_data_meta <- nafld_data %>% select(!all_of(mbx_list)) %>% rownames_to_column("sample_id")  
nafld_data_mbx <- nafld_data %>% select(all_of(mbx_list)) 
nafld_data_mbx_id <- nafld_data %>% column_to_rownames("barcode_metabolomics") %>% select(all_of(mbx_list)) 

###Extended Data Figure 2A
###Regular PCoA -- NAFLD                                         
bray <- nafld_data_mbx %>% vegdist(., "bray") 
pc = capscale(bray~1, comm = nafld_data_mbx)
pc.summary<-summary(pc)

pcl.bray <- as.data.frame(pc$CA$u) %>%
  select(MDS1, MDS2)

pcl.bray <- pcl.bray %>%
  rownames_to_column("sample_id") %>%
  inner_join(nafld_data_meta, by = "sample_id") %>%
  select(colnames(nafld_data_meta), everything()) 


# save axes r2
pco1.r2 <- paste("PCo1 (", round(pc.summary$cont$importance[2,1]*100, digits = 1), "%)", sep = '')
pco2.r2 <- paste("PCo2 (", round(pc.summary$cont$importance[2,2]*100, digits = 1), "%)", sep = '')
rm(pc, pc.summary)

ggplot(pcl.bray, 
       aes(MDS1, MDS2)) + 
  geom_point(aes(color = as.factor(case)), alpha = 0.7, size = 2, stroke = 1) + 
  scale_color_manual(values = c("0" = "#999999", "1" = "#E69F00")) +
  coord_fixed() + 
  theme_bw(base_size=24) + 
  ggtitle("MASLD") + 
  labs(x = pco1.r2, 
       y = pco2.r2,
       fill = "") +
  labs(color = "MASLD case") 
#10*6

ggsave(
  file.path("output_mbx", "pco.pdf"),
  dpi = 300
)
