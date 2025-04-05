##################################################
#R program for creating acylcarnitine heatmaps
#Figure 4C & Extended Data Figure 2
##################################################

library(ggplot2)
library(ggpubr)
library(dplyr)
library(tidyverse)
library(stringr)
library(readr)
library(reshape2)
library(gridExtra)
library(grid)
library(pheatmap)

#####################################################################
###                           species                            ####
#####################################################################
setwd("~/b2b")

tax<-read_tsv("metaphlan_taxonomic_profiles.tsv")

#metadata
final_metadata <- read.delim('input/meta_df.tsv',row.names=1)
final_metadata <- final_metadata[!duplicated(final_metadata$barcode_metagenomics),]

#process metaphlan output
tax <- tax %>% rename("taxonomy"="# taxonomy")
names(tax) = gsub("_taxonomic_profile", "", names(tax))
tax_t<-tax[grepl("t__", tax$taxonomy), ] #only use the ones with t__ level
levels <- c("Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species", "Species2")
tree_order<-tax_t %>% separate(taxonomy, sep = "\\|", into = levels, remove = FALSE) %>% select("Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species", "Species2") %>% rename("SGB"="Species2")
tax_species<-tax_t %>% separate(taxonomy, sep = "\\|", into = levels, remove = FALSE) %>% 
  select(-c("taxonomy","Domain", "Phylum", "Class", "Order", "Family", "Genus")) %>%
  mutate(Species=gsub("SGB\\d+\\s?", "", Species)) %>%
  unite("speciesname",Species:Species2,remove=TRUE) %>% 
  mutate(speciesname = gsub("s__","",speciesname)) %>% mutate(speciesname = gsub("t__","",speciesname)) %>%
  mutate(speciesname = gsub("__","_",speciesname)) %>%
  column_to_rownames("speciesname") %>% 
  t() 
tax_species<-tax_species/100
tax_species_for_merge<-tax_species%>% 
  as.data.frame() %>%
  rownames_to_column() %>% 
  rename("barcode_metagenomics"="rowname")

df_w_meta <- left_join(tax_species_for_merge,final_metadata,by="barcode_metagenomics") %>% 
  distinct(barcode_metagenomics, .keep_all = TRUE) 

nafld_data<-df_w_meta %>% filter(cohort=="NAFLD" | case==0) %>% 
  mutate(lean_vs_nonlean_case = case_when(bmi17v >= 25 & case==1 ~ 'nonlean_nafld', bmi17v < 25 & case==1 ~ 'lean_nafld')) #nonlean case is 1, lean case is 0

nafld_data_species <- nafld_data %>% column_to_rownames("alias_id") %>% select(matches('SGB')) 
taxa_t <- nafld_data_species / 100

# change to numeric
for(i in 1:ncol(taxa_t))
{
  taxa_t[ , i] <- as.numeric(as.character(taxa_t[, i]))
}

AST <- function(x) {
  y <- sign(x) * asin(sqrt(abs(x)))
  if(any(is.na(y))) {
    logging::logerror(
      paste0("AST transform only valid for values between -1 and 1"))
    stop()
  }
  return(y)
}
taxa_t_ast <- apply(taxa_t,2,AST)

#####################################################################
###                          metabolites                         ####
#####################################################################

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

nafld_data<-df_w_meta %>% filter(cohort=="NAFLD" | case==0) %>% 
  mutate(lean_vs_nonlean_case = case_when(bmi17v >= 25 & case==1 ~ 'nonlean_nafld', bmi17v < 25 & case==1 ~ 'lean_nafld')) #nonlean case is 1, lean case is 0

nafld_data_mbx <- nafld_data %>% select(all_of(mbx_list)) 
#took average of 10 technically duplicated samples
nafld_data_mbx_id <- nafld_data %>%
  group_by(alias_id) %>%
  summarize(across(all_of(mbx_list), mean)) %>%
  column_to_rownames("alias_id")

mbx_t<-nafld_data_mbx_id
# change to numeric
for(i in 1:ncol(mbx_t))
{
  mbx_t[ , i] <- as.numeric(as.character(mbx_t[, i]))
}
#note, mbx data doesnt add up to 1 and it is already normalized based on median
# Log transformation
LOG <- function(x) {
  y <- replace(x, x == 0, min(x[x>0]) / 2)
  return(log2(y))
}

mbx_t_log <- apply(mbx_t, 2, LOG)

################################################################################
###                         merge species with metabolites                  ####
################################################################################
mbx_taxa <- merge(taxa_t_ast,mbx_t_log,by='row.names')
rownames(mbx_taxa) <- mbx_taxa[,1]
mbx_taxa <- mbx_taxa[,-1]

#################################################################
###                           metadata                       ####
#################################################################
meta <- read.delim('input/meta_df.tsv',row.names=1) %>%
  mutate(lean_vs_nonlean_case = case_when(bmi17v >= 25 & case==1 ~ 'nonlean_nafld', bmi17v < 25 & case==1 ~ 'lean_nafld')) #nonlean case is 1, lean case is 0
meta <- meta[!duplicated(meta$alias_id),]

rpmeat_and_fiber <- read_csv('input/redmeat.fiber.aliasid.csv') 
meta <- left_join(meta, rpmeat_and_fiber, by="alias_id")

meta <- meta[,colnames(meta) %in% c("alias_id","lean_vs_nonlean_case",'case','aheiv2010_15','rpmeats15v', 'aofib15v')]
rownames(meta) <- meta$alias_id

all <- merge(mbx_taxa, meta, by = 'row.names')
rownames(all) <- all[, 1]
all <- all[, -1]

##### CAR and MASLD bugs
significant_results <- read.delim("~/b2b/output.mp4/significant_results.tsv", header = TRUE, sep = "\t")
significant_bugs <- significant_results %>%
  filter(value == "case") %>%
  pull(feature) %>%
  unique()

all_control <- all %>% filter(case==0) %>% select(-c(alias_id,case,lean_vs_nonlean_case))
all_nonlean <- all %>% filter(lean_vs_nonlean_case=="nonlean_nafld") %>% select(-c(alias_id,case,lean_vs_nonlean_case))
all_lean <- all %>% filter(lean_vs_nonlean_case=="lean_nafld") %>% select(-c(alias_id,case,lean_vs_nonlean_case))
all_case <- all %>% filter(lean_vs_nonlean_case=="nonlean_nafld" | lean_vs_nonlean_case=="lean_nafld") %>% select(-c(alias_id,case,lean_vs_nonlean_case))

######################################################
###        main heatmap (MASLD Bugs vs. CAR)       ###
######################################################

species_of_interest <- intersect(colnames(all_case), significant_bugs)

carnitines <- colnames(mbx_t_log)[grepl("^CAR \\d+:\\d+", colnames(mbx_t_log))]
carnitines_chain_order <- carnitines[order(as.numeric(gsub("CAR (\\d+):.*", "\\1", carnitines)))]

enrichment_status <- significant_results %>%
  filter(value == "case") %>%
  filter(feature %in% significant_bugs) %>%
  mutate(MASLD = ifelse(coef > 0, "Enriched", "Depleted")) %>%
  select(feature, MASLD) %>%
  distinct() %>%
  column_to_rownames("feature")

calculate_spearman <- function(data, species, metabolites) {
  results <- data.frame(Species = character(), Metabolite = character(), Correlation = numeric(), P_value = numeric())
  
  for (met in metabolites) {
    for (sp in species) {
      cor_test <- cor.test(data[[sp]], data[[met]], method = "spearman", use = "pairwise.complete.obs")
      results <- rbind(results, data.frame(
        Species = sp,
        Metabolite = met,
        Correlation = cor_test$estimate,
        P_value = cor_test$p.value
      ))
    }
  }
  
  #multiple testing correction
  results <- results %>%
    group_by(Metabolite) %>%
    mutate(FDR = p.adjust(P_value, method = "BH")) %>%
    ungroup() %>%
    mutate(Significance = ifelse(FDR < 0.20, "*", ""))
  
  return(results)
}

results <- calculate_spearman(all_case, species_of_interest, carnitines_chain_order)

heatmap_data <- results %>%
  select(Species, Metabolite, Correlation) %>%
  pivot_wider(names_from = Metabolite, values_from = Correlation) %>%
  column_to_rownames("Species")

significance_matrix <- results %>%
  select(Species, Metabolite, Significance) %>%
  pivot_wider(names_from = Metabolite, values_from = Significance) %>%
  column_to_rownames("Species")

rownames(heatmap_data) <- gsub("_", " ", rownames(heatmap_data))
rownames(significance_matrix) <- gsub("_", " ", rownames(significance_matrix))

############################################################
###        left strip (AHEI & fiber for species)         ###
############################################################

ahei_results <- calculate_spearman(all_case, species_of_interest, c("aheiv2010_15")) %>%
  select(Species, Correlation) %>%
  rename(AHEI = Correlation)

fiber_results <- calculate_spearman(all_case, species_of_interest, c("aofib15v")) %>%
  select(Species, Correlation) %>%
  rename(Fiber = Correlation)

enrichment_status <- enrichment_status %>%
  rownames_to_column(var = "Species") 

species_annotations <- left_join(ahei_results, fiber_results, by = "Species") %>%
  left_join(enrichment_status, by = "Species") %>%
  column_to_rownames("Species")

rownames(species_annotations) <- gsub("_", " ", rownames(species_annotations))

#################################################
###      top strip (AHEI & fiber for CAR)     ###
#################################################

ahei_results_car <- calculate_spearman(all_case, carnitines_chain_order, c("aheiv2010_15")) %>%
  select(Species, Correlation) %>%
  rename(AHEI = Correlation)

fiber_results_car <- calculate_spearman(all_case, carnitines_chain_order, c("aofib15v")) %>%
  select(Species, Correlation) %>%
  rename(Fiber = Correlation)

carnitine_annotations <- left_join(ahei_results_car, fiber_results_car, by = "Species") %>%
  column_to_rownames("Species")

rownames(carnitine_annotations) <- gsub("_", " ", rownames(carnitine_annotations))

annotation_colors <- list(
  MASLD = c("Enriched" = "#E69F00", "Depleted" = "#999999"), 
  AHEI = colorRampPalette(c("blue", "white", "red"))(50),
  Fiber = colorRampPalette(c("blue", "white", "red"))(50)
)

pheatmap(
  heatmap_data,
  display_numbers = significance_matrix,  #asterisks for significant correlations
  main = "Among MASLD cases (N=209)",
  fontsize_number = 15,
  color = colorRampPalette(c("blue", "white", "red"))(100),  
  border_color = "grey",
  cluster_rows = TRUE,
  cluster_cols = FALSE,
  annotation_row = species_annotations,   
  annotation_col = carnitine_annotations, 
  annotation_colors = annotation_colors,
  annotation_names_row = TRUE,
  angle_col=90
)

###Extended Data Figure 2B
##among controls
species_of_interest <- intersect(colnames(all_control), significant_bugs)

carnitines <- colnames(mbx_t_log)[grepl("^CAR \\d+:\\d+", colnames(mbx_t_log))]
carnitines_chain_order <- carnitines[order(as.numeric(gsub("CAR (\\d+):.*", "\\1", carnitines)))]

enrichment_status <- significant_results %>%
  filter(value == "case") %>%
  filter(feature %in% significant_bugs) %>%
  mutate(MASLD = ifelse(coef > 0, "Enriched", "Depleted")) %>%
  select(feature, MASLD) %>%
  distinct() %>%
  column_to_rownames("feature")

calculate_spearman <- function(data, species, metabolites) {
  results <- data.frame(Species = character(), Metabolite = character(), Correlation = numeric(), P_value = numeric())
  
  for (met in metabolites) {
    for (sp in species) {
      cor_test <- cor.test(data[[sp]], data[[met]], method = "spearman", use = "pairwise.complete.obs")
      results <- rbind(results, data.frame(
        Species = sp,
        Metabolite = met,
        Correlation = cor_test$estimate,
        P_value = cor_test$p.value
      ))
    }
  }
  
  #multiple testing correction 
  results <- results %>%
    group_by(Metabolite) %>%
    mutate(FDR = p.adjust(P_value, method = "BH")) %>%
    ungroup() %>%
    mutate(Significance = ifelse(FDR < 0.20, "*", ""))
  
  return(results)
}

results <- calculate_spearman(all_control, species_of_interest, carnitines_chain_order)

heatmap_data <- results %>%
  select(Species, Metabolite, Correlation) %>%
  pivot_wider(names_from = Metabolite, values_from = Correlation) %>%
  column_to_rownames("Species")

significance_matrix <- results %>%
  select(Species, Metabolite, Significance) %>%
  pivot_wider(names_from = Metabolite, values_from = Significance) %>%
  column_to_rownames("Species")

rownames(heatmap_data) <- gsub("_", " ", rownames(heatmap_data))
rownames(significance_matrix) <- gsub("_", " ", rownames(significance_matrix))

ahei_results <- calculate_spearman(all_control, species_of_interest, c("aheiv2010_15")) %>%
  select(Species, Correlation) %>%
  rename(AHEI = Correlation)

fiber_results <- calculate_spearman(all_control, species_of_interest, c("aofib15v")) %>%
  select(Species, Correlation) %>%
  rename(Fiber = Correlation)

enrichment_status <- enrichment_status %>%
  rownames_to_column(var = "Species") 

species_annotations <- left_join(ahei_results, fiber_results, by = "Species") %>%
  left_join(enrichment_status, by = "Species") %>%
  column_to_rownames("Species")

rownames(species_annotations) <- gsub("_", " ", rownames(species_annotations))

ahei_results_car <- calculate_spearman(all_control, carnitines_chain_order, c("aheiv2010_15")) %>%
  select(Species, Correlation) %>%
  rename(AHEI = Correlation)

fiber_results_car <- calculate_spearman(all_control, carnitines_chain_order, c("aofib15v")) %>%
  select(Species, Correlation) %>%
  rename(Fiber = Correlation)

carnitine_annotations <- left_join(ahei_results_car, fiber_results_car, by = "Species") %>%
  column_to_rownames("Species")

rownames(carnitine_annotations) <- gsub("_", " ", rownames(carnitine_annotations))

annotation_colors <- list(
  MASLD = c("Enriched" = "#E69F00", "Depleted" = "#999999"), 
  AHEI = colorRampPalette(c("blue", "white", "red"))(50),
  Fiber = colorRampPalette(c("blue", "white", "red"))(50)
)

pheatmap(
  heatmap_data,
  display_numbers = significance_matrix,  #asterisks for significant correlations
  main = "Among controls (N=478)",
  fontsize_number = 15,
  color = colorRampPalette(c("blue", "white", "red"))(100), 
  border_color = "grey",
  cluster_rows = TRUE,
  cluster_cols = FALSE,
  annotation_row = species_annotations,  
  annotation_col = carnitine_annotations, 
  annotation_colors = annotation_colors,
  annotation_names_row = TRUE,
  angle_col=90
)
