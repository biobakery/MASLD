##################################################
#R program for creating Figure 1B
##################################################

# Load necessary libraries
library(vegan)
library(dplyr)
library(tidyverse)
library(stringr)
library(readr)
library(Maaslin2)
library(ggplot2)
library(ggrepel)
library(viridis)
library(ggsignif)
library(lme4)
library(reshape2)

setwd("~/b2b")
unfiltered_species <- read.delim('input/species_v4.tsv',row.names=1) %>% t() %>% as.data.frame() %>% rownames_to_column() %>% rename(barcode_metagenomics = rowname) 
final_metadata <- read.delim('input/meta_df.tsv',row.names=1)
final_metadata <- final_metadata[!duplicated(final_metadata$barcode_metagenomics),]
df_w_meta <- left_join(unfiltered_species,final_metadata,by="barcode_metagenomics") %>%
  mutate(lean_vs_nonlean_case = case_when(bmi17v >= 25 & case==1 ~ 1, bmi17v < 25 & case==1 ~ 0)) #nonlean case is 1, lean case is 0

species.data<-df_w_meta %>% select(starts_with('s__'))
species.data<-species.data/100
species.data <- species.data %>% t() %>% as.data.frame() %>% rownames_to_column() %>%
  mutate(rowname=gsub("SGB(\\d+)_SGB\\1", "SGB\\1", rowname)) %>%
  mutate(rowname = gsub("s__","",rowname)) %>% column_to_rownames("rowname") %>% t()

nafld.data<-df_w_meta %>% filter(cohort=="NAFLD") 
species.nafld.data<-nafld.data %>% select(starts_with('s__')) 
species.nafld.data <- species.nafld.data %>% t() %>% as.data.frame() %>% rownames_to_column() %>%
  mutate(rowname=gsub("SGB(\\d+)_SGB\\1", "SGB\\1", rowname)) %>%
  mutate(rowname = gsub("s__","",rowname)) %>% column_to_rownames("rowname") %>% t()

#filtering: 
species_maaslin_allresults <- read_tsv("all_results_species_v4.tsv")
features<-species_maaslin_allresults$feature
species.data.filtered<-species.data %>% t() %>% as.data.frame() 
species.data.filtered<-species.data.filtered%>% filter(row.names(species.data.filtered) %in% features) %>% t() %>% as.matrix() 

nafld_data<-df_w_meta %>% filter(cohort=="NAFLD" | case==0) %>% 
  mutate(lean_vs_nonlean_case = case_when(bmi17v >= 25 & case==1 ~ 1, bmi17v < 25 & case==1 ~ 0)) #nonlean case is 1, lean case is 0

sp_w_meta <- df_w_meta %>% filter(cohort=="NAFLD" | case==0)
names(sp_w_meta) <- gsub(names(sp_w_meta), pattern = "SGB(\\d+)_SGB\\1", replacement = "SGB\\1")  
names(sp_w_meta) <- gsub(names(sp_w_meta), pattern = "s__", replacement = "")  

sp_df <- species.data.filtered
sp_NAFLD_df<-sp_w_meta %>% mutate(obesity = case_when(bmi17v >= 30 ~ 1, bmi17v <30 ~ 0)) %>% mutate(by_bmi = case_when(bmi17v >= 30 & case==1 ~ 'Obese NAFLD', bmi17v <30 & case==1 ~ 'Lean NAFLD', bmi17v >= 30 & case==0 ~ 'Obese control', bmi17v <30 & case==0 ~ 'Lean control')) 

#####mbx
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
df_w_meta <- left_join(df_final_join,final_metadata,by="barcode_metabolomics") %>%
  mutate(lean_vs_nonlean_case = case_when(bmi17v >= 25 & case==1 ~ 1, bmi17v < 25 & case==1 ~ 0)) #nonlean case is 1, lean case is 0

mbx_list<-names(df_final)
mbx.data<-df_w_meta %>% select(all_of(mbx_list))

nafld.data<-df_w_meta %>% filter(cohort=="NAFLD") 
mbx.nafld.data<-nafld.data %>% select(all_of(mbx_list)) 

mbx.data.filtered<-mbx.data %>% t() %>% as.data.frame

#no filtering
mbx.data.filtered<-mbx.data.filtered%>% t() 
mbx_filtered_list <- colnames(mbx.data.filtered)

nafld_data<-df_w_meta %>% filter(cohort=="NAFLD" | case==0) %>% 
  mutate(lean_vs_nonlean_case = case_when(bmi17v >= 25 & case==1 ~ 1, bmi17v < 25 & case==1 ~ 0)) #nonlean case is 1, lean case is 0

nafld_data_mbx <- nafld_data %>% select(all_of(mbx_list)) 

mbx_w_meta <- df_w_meta
mbx_df <- mbx.data.filtered
mbx_NAFLD_df<-mbx_w_meta %>% filter(cohort=="NAFLD" | case==0) %>% mutate(obesity = case_when(bmi17v >= 30 ~ 1, bmi17v <30 ~ 0)) %>% mutate(by_bmi = case_when(bmi17v >= 30 & case==1 ~ 'Obese NAFLD', bmi17v <30 & case==1 ~ 'Lean NAFLD', bmi17v >= 30 & case==0 ~ 'Obese control', bmi17v <30 & case==0 ~ 'Lean control'))

#####Virome
virome_profile <- read.delim("MGXBAQLaVa_VGB_table.tsv",row.names=1) %>% t() %>% as.data.frame() %>% rownames_to_column() %>% rename(barcode_metagenomics = rowname) 
virome_profile$barcode_metagenomics <- gsub("_Abundance.RPKs", "", virome_profile$barcode_metagenomics)

virome_profile <- virome_profile %>% column_to_rownames("barcode_metagenomics")
virome_profile$total_rpks <- rowSums(virome_profile)

# Convert RPK values to relative abundance
rpks_df_relative <- virome_profile / virome_profile$total_rpks
virome_profile_new <- rpks_df_relative %>% select(-total_rpks) %>% rownames_to_column("barcode_metagenomics") 

virome_list<-virome_profile_new %>% select(-barcode_metagenomics)

final_metadata <- read.delim('input/meta_df.tsv',row.names=1)
final_metadata <- final_metadata[!duplicated(final_metadata$barcode_metagenomics),]
df_w_meta <- left_join(virome_profile_new,final_metadata,by="barcode_metagenomics") %>%
  mutate(lean_vs_nonlean_case = case_when(bmi17v >= 25 & case==1 ~ 1, bmi17v < 25 & case==1 ~ 0)) #nonlean case is 1, lean case is 0

virome.data<-df_w_meta %>% select(c(names(virome_list)))

nafld.data<-df_w_meta %>% filter(cohort=="NAFLD") 
virome.nafld.data<-nafld.data %>% select(c(names(virome_list)))

virome_maaslin_allresults <- read_tsv("output_virome_v0.3/all_results.tsv")
virome<-virome_maaslin_allresults$feature
virome.data.filtered<-virome.data %>% t() %>% as.data.frame 
virome.data.filtered<-virome.data.filtered%>% filter(row.names(virome.data.filtered) %in% virome) %>% t() 
virome.data.filtered<-virome.data.filtered%>% as.data.frame()

virome.nafld.data<-df_w_meta %>% filter(cohort=="NAFLD" | case==0) %>% select(c(names(virome_list),"alias_id")) %>% column_to_rownames("alias_id")
virome.nafld.data<-virome.nafld.data[order(row.names(virome.nafld.data)), ]

nafld_data<-df_w_meta %>% filter(cohort=="NAFLD" | case==0) %>% 
  mutate(lean_vs_nonlean_case = case_when(bmi17v >= 25 & case==1 ~ 1, bmi17v < 25 & case==1 ~ 0)) #nonlean case is 1, lean case is 0

nafld_data_virome <- nafld_data %>% select(c(names(virome_list)))
unclassified_virome <- nafld_data_virome %>% select(contains("Unclassified")) 

virome_w_meta <- df_w_meta
virome_df <- virome.data.filtered
virome_NAFLD_df<-virome_w_meta %>% filter(cohort=="NAFLD" | case==0) %>% mutate(obesity = case_when(bmi17v >= 30 ~ 1, bmi17v <30 ~ 0)) %>% mutate(by_bmi = case_when(bmi17v >= 30 & case==1 ~ 'Obese NAFLD', bmi17v <30 & case==1 ~ 'Lean NAFLD', bmi17v >= 30 & case==0 ~ 'Obese control', bmi17v <30 & case==0 ~ 'Lean control'))

#####MGX pathways
mgx_pathway_file <- read.delim('pathabundance_unstrat.tsv',row.names=1) %>% t() %>% as.data.frame() %>% rownames_to_column() %>% rename(barcode_metagenomics = rowname)
mgx_pathway_file$barcode_metagenomics <- gsub("_Abundance", "", mgx_pathway_file$barcode_metagenomics)

df_w_meta <- left_join(mgx_pathway_file, final_metadata, by = "barcode_metagenomics") %>%
  mutate(lean_vs_nonlean_case = case_when(bmi17v >= 25 & case==1 ~ 1, bmi17v < 25 & case==1 ~ 0)) #nonlean case is 1, lean case is 0

pathway_list<-read.delim('pathabundance_unstrat.tsv',row.names=1) %>% rownames_to_column()

mgx_pathway_prevalence <- df_w_meta %>%
  select(pathway_list$rowname) 

mgx_pathway_prevalence_summary <- mgx_pathway_prevalence %>%
  summarise_all(~ mean(. > 0))

mgx_pathway_prevalence_summary <- as.data.frame(t(mgx_pathway_prevalence_summary)) 
colnames(mgx_pathway_prevalence_summary) <- c("prevalence") 

mgx_pathway_prevalence_summary <- mgx_pathway_prevalence_summary %>%
  rownames_to_column(var = "feature")

# Filter features with at least 10% prevalence
pathway_features_to_keep <- mgx_pathway_prevalence_summary %>%
  filter(prevalence >= 0.1) %>%
  pull(feature) 

mgx_pathway_filtered <- df_w_meta %>%
  select(all_of(pathway_features_to_keep), barcode_metagenomics) 

nafld_data <- df_w_meta %>% 
  filter(cohort == "NAFLD" | case == 0) %>% 
  mutate(lean_vs_nonlean_case = ifelse(bmi17v >= 25 & case == 1, 1, 0)

mgx_pathway_w_meta <- df_w_meta %>% filter(cohort=="NAFLD" | case==0)

mgx_pathway_df <- mgx_pathway_filtered %>% select(-barcode_metagenomics)
mgx_pathway_NAFLD_df<-mgx_pathway_w_meta %>% mutate(obesity = case_when(bmi17v >= 30 ~ 1, bmi17v <30 ~ 0)) %>% mutate(by_bmi = case_when(bmi17v >= 30 & case==1 ~ 'Obese NAFLD', bmi17v <30 & case==1 ~ 'Lean NAFLD', bmi17v >= 30 & case==0 ~ 'Obese control', bmi17v <30 & case==0 ~ 'Lean control')) 

#####MTX pathways
mtxpathway<-read_tsv("pathabundance_relab_nospecial.tsv")
names(mtxpathway) = gsub("_Abundance", "", names(mtxpathway))
mtx_pathway<-mtxpathway %>% column_to_rownames("# Pathway") %>% t() %>% as.data.frame() 
samples <- mtx_pathway %>% rownames_to_column()
samplenames = samples$rowname 
mtx_pathway_unstratified<-mtx_pathway[,!grepl("\\|", colnames(mtx_pathway)) ] 
mtx_pathway_for_join<-mtx_pathway_unstratified %>% rownames_to_column("barcode_metagenomics")
pathways_list <- colnames(mtx_pathway_unstratified)

df_w_meta <- left_join(mtx_pathway_for_join, final_metadata, by = "barcode_metagenomics") %>%
  mutate(lean_vs_nonlean_case = case_when(bmi17v >= 25 & case==1 ~ 1, bmi17v < 25 & case==1 ~ 0)) #nonlean case is 1, lean case is 0

mtx_pathway_prevalence <- df_w_meta %>%
  select(pathways_list) 

mtx_pathway_prevalence_summary <- mtx_pathway_prevalence %>%
  summarise_all(~ mean(. > 0))

mtx_pathway_prevalence_summary <- as.data.frame(t(mtx_pathway_prevalence_summary)) 
colnames(mtx_pathway_prevalence_summary) <- c("prevalence") 

mtx_pathway_prevalence_summary <- mtx_pathway_prevalence_summary %>%
  rownames_to_column(var = "feature")

# Filter features with at least 10% prevalence
pathway_features_to_keep <- mtx_pathway_prevalence_summary %>%
  filter(prevalence >= 0.1) %>%
  pull(feature) 

mtx_pathway_filtered <- df_w_meta %>%
  select(all_of(pathway_features_to_keep), barcode_metagenomics) # Retain metadata column if needed

nafld_data <- df_w_meta %>% 
  filter(cohort == "NAFLD" | case == 0) %>% 
  mutate(lean_vs_nonlean_case = ifelse(bmi17v >= 25 & case == 1, 1, 0)

mtx_pathway_w_meta <- df_w_meta %>% filter(cohort=="NAFLD" | case==0)

mtx_pathway_df <- mtx_pathway_filtered %>% select(-barcode_metagenomics)
mtx_pathway_NAFLD_df<-mtx_pathway_w_meta %>% mutate(obesity = case_when(bmi17v >= 30 ~ 1, bmi17v <30 ~ 0)) %>% mutate(by_bmi = case_when(bmi17v >= 30 & case==1 ~ 'Obese NAFLD', bmi17v <30 & case==1 ~ 'Lean NAFLD', bmi17v >= 30 & case==0 ~ 'Obese control', bmi17v <30 & case==0 ~ 'Lean control')) 

# Define datasets
datasets <- list(
  MTX_functions = mtx_pathway_NAFLD_df, 
  MGX_functions = mgx_pathway_NAFLD_df,
  Metabolites = mbx_NAFLD_df,
  Viruses = virome_NAFLD_df,
  Bacteria = sp_NAFLD_df
)

# Define covariates
covariates <- c("case", "lean_vs_nonlean_case", "age", "bmi17v", "db17", "act17v", "aheiv2010_15")

features <- list(
  MTX_functions = mtx_pathway_df,
  MGX_functions = mgx_pathway_df,
  Metabolites = mbx_df,
  Viruses = virome_df,
  Bacteria = sp_df
)

R2_matrix <- matrix(NA, nrow = length(covariates), ncol = length(datasets))
rownames(R2_matrix) <- covariates
colnames(R2_matrix) <- names(datasets)

P_matrix <- R2_matrix

# Loop over datasets and covariates
for (dataset_name in names(datasets)) {
  dataset <- datasets[[dataset_name]]
  feature_columns <- colnames(features[[dataset_name]])
  
  numeric_data <- dataset[, colnames(dataset) %in% feature_columns]
  
  for (covariate in covariates) {
    if (covariate %in% colnames(dataset)) {
      # Remove rows with NA in the covariate or numeric data
      data_clean <- dataset[!is.na(dataset[[covariate]]), ]
      numeric_data_clean <- numeric_data[rownames(data_clean), , drop = FALSE]
      numeric_data_clean <- numeric_data_clean[rowSums(numeric_data_clean, na.rm = TRUE) > 0, , drop = FALSE]
      
      # Ensure numeric data and metadata align
      if (nrow(numeric_data_clean) > 1 && ncol(numeric_data_clean) > 0) {
        # Create a distance matrix
        distance_matrix <- vegdist(numeric_data_clean, method = "bray")
        
        common_rows <- rownames(as.matrix(distance_matrix))
        data_clean <- data_clean[common_rows, , drop = FALSE]
        
        # Run PERMANOVA
        result <- tryCatch({
          adonis2(
            distance_matrix ~ data_clean[[covariate]],
            data = data_clean,
            permutations = 999,
            by = "margin"
          )
        }, error = function(e) {
          message(paste("Error with covariate:", covariate, "in dataset:", dataset_name, "-", e$message))
          return(NULL)
        })
        
        if (!is.null(result)) {
          R2_matrix[covariate, dataset_name] <- result$R2[1]
          P_matrix[covariate, dataset_name] <- result$`Pr(>F)`[1]
        } else {
          R2_matrix[covariate, dataset_name] <- NA
          P_matrix[covariate, dataset_name] <- NA
        }
      } else {
        message(paste("Insufficient numeric data in dataset:", dataset_name, "for covariate:", covariate))
      }
    } else {
      message(paste("Covariate not found:", covariate, "in dataset:", dataset_name))
    }
  }
}


P_matrix_adj <- apply(P_matrix, 2, p.adjust, method = "fdr")

heatmap_data <- melt(R2_matrix)
colnames(heatmap_data) <- c("Covariate", "Dataset", "R2")
heatmap_data$P <- melt(P_matrix_adj)$value

heatmap_data$R2_scaled <- ifelse(heatmap_data$R2 <= 0.03, heatmap_data$R2, 0.03 + (heatmap_data$R2 - 0.03) / 3)

heatmap_data$Significance <- cut(
  heatmap_data$P,
  breaks = c(-Inf, 0.001, 0.01, 0.05, Inf),
  labels = c("***", "**", "*", "")
)

ggplot(heatmap_data, aes(x = Covariate, y = Dataset, fill = R2_scaled)) +
  geom_tile(color = "white") +
  geom_text(
    aes(
      label = sprintf("%.1f%% %s", R2 * 100, Significance),  
      color = ifelse(R2 > 0.015, "white", "black")  
    ),
    size = 3
  ) +
  scale_fill_gradientn(
    colors = c("white", "#D7E4F2", "#A6BDD6", "#6E97C7", "#2962A1"),
    name = "Variance Explained (%)",
    limits = c(0, max(heatmap_data$R2, na.rm = TRUE)),
    trans = "sqrt",  
    breaks = c(0, 0.005, 0.01, 0.02, 0.03, max(heatmap_data$R2, na.rm = TRUE)), 
    labels = scales::percent_format(accuracy = 0.1),
    na.value = "white"
  ) +
  scale_color_manual(values = c("black" = "black", "white" = "white"), guide = "none") +
  labs(title = "PERMANOVA Heatmap", x = "Dataset", y = "Covariate") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.text.y = element_text(size = 10),
    plot.title = element_text(hjust = 0.5, face = "bold", size = 14)
  )
