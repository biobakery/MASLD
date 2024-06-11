#!/usr/bin/env Rscript
##################################################
#R program for performing PERMANOVA
#Figure 5A
#Extended Data Figure 1A
#Extended Data Figure 2A
##################################################

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

setwd("~/b2b")
unfiltered_species <- read.delim('input/species_v4.tsv',row.names=1) %>% t() %>% as.data.frame() %>% rownames_to_column() %>% rename(barcode_metagenomics = rowname)
final_metadata <- read.delim('input/meta_df.tsv',row.names=1)
final_metadata <- final_metadata[!duplicated(final_metadata$barcode_metagenomics),]
df_w_meta <- left_join(unfiltered_species,final_metadata,by="barcode_metagenomics")

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

###########
###NAFLD###
###########
#211 nafld cases and 502 controls (193 matched and 309 unmatched controls)
nafld_data<-df_w_meta %>% filter(cohort=="NAFLD" | case==0) %>%
  mutate(obesity = case_when(bmi17v >= 30 ~ 1, bmi17v <30 ~ 0)) %>%
  mutate(lean = case_when(bmi17v < 25 ~ 1, bmi17v >= 25 ~ 0)) %>%
  mutate(lean_nafld_binary = ifelse(bmi17v < 25 & case==1, 1, 0)) %>%
  mutate(nonlean_nafld_binary = ifelse(bmi17v >= 25 & case==1, 1, 0)) %>%
  mutate(lean_nafld_lean_control = case_when(bmi17v < 25 & case==1 ~ 1, bmi17v < 25 & case==0 ~ 0)) %>%
  mutate(nonlean_nafld_nonlean_control = case_when(bmi17v >= 25 & case==1 ~ 1, bmi17v >= 25 & case==0 ~ 0)) %>%
  mutate(lean_vs_nonlean_case = case_when(bmi17v >= 25 & case==1 ~ 1, bmi17v < 25 & case==1 ~ 0)) #nonlean case is 1, lean case is 0

sp_w_meta <- df_w_meta %>% filter(cohort=="NAFLD" | case==0)
names(sp_w_meta) <- gsub(names(sp_w_meta), pattern = "SGB(\\d+)_SGB\\1", replacement = "SGB\\1")
names(sp_w_meta) <- gsub(names(sp_w_meta), pattern = "s__", replacement = "")

sp_df <- species.data.filtered
sp_NAFLD_df<-sp_w_meta

#Extended Data Figure 1A
sp_NAFLD_perm <- adonis2(sp_NAFLD_df[, colnames(sp_NAFLD_df) %in% colnames(sp_df)] ~case + age + aheiv2010_15 + act17v + bmi17v+ db17 + plate, data = sp_NAFLD_df, permutations = 999, method = "bray", sqrt.dist = FALSE, add = FALSE, by = 'margin')
sp_NAFLD_perm

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
df_w_meta <- left_join(df_final_join,final_metadata,by="barcode_metabolomics")

mbx_list<-names(df_final)
mbx.data<-df_w_meta %>% select(all_of(mbx_list))

nafld.data<-df_w_meta %>% filter(cohort=="NAFLD")
mbx.nafld.data<-nafld.data %>% select(all_of(mbx_list))

mbx_maaslin_allresults <- read_tsv("output_mbx/nonlog/all_results.tsv")
mbx<-mbx_maaslin_allresults$feature
mbx.data.filtered<-mbx.data %>% t() %>% as.data.frame
row.names(mbx.data.filtered) <- gsub("^(\\(|\\d)", "X\\1", row.names(mbx.data.filtered))
row.names(mbx.data.filtered) = gsub("-", ".", row.names(mbx.data.filtered))
row.names(mbx.data.filtered) = gsub(":", ".", row.names(mbx.data.filtered))
row.names(mbx.data.filtered) = gsub("\\(", ".", row.names(mbx.data.filtered))
row.names(mbx.data.filtered) = gsub("\\)", ".", row.names(mbx.data.filtered))
row.names(mbx.data.filtered) = gsub("\\[", ".", row.names(mbx.data.filtered))
row.names(mbx.data.filtered) = gsub("\\]", ".", row.names(mbx.data.filtered))
row.names(mbx.data.filtered) = gsub(" ", ".", row.names(mbx.data.filtered))
row.names(mbx.data.filtered) = gsub("'", ".", row.names(mbx.data.filtered))
row.names(mbx.data.filtered) = gsub("/", ".", row.names(mbx.data.filtered))
row.names(mbx.data.filtered) = gsub(";", ".", row.names(mbx.data.filtered))
row.names(mbx.data.filtered) = gsub(",", ".", row.names(mbx.data.filtered))
mbx.data.filtered<-mbx.data.filtered%>% filter(row.names(mbx.data.filtered) %in% mbx) %>% t()
mbx_filtered_list <- colnames(mbx.data.filtered)
###########
###NAFLD###
###########
nafld_data<-df_w_meta %>% filter(cohort=="NAFLD" | case==0) %>%
  mutate(obesity = case_when(bmi17v >= 30 ~ 1, bmi17v <30 ~ 0)) %>%
  mutate(lean = case_when(bmi17v < 25 ~ 1, bmi17v >= 25 ~ 0)) %>%
  mutate(lean_nafld_binary = ifelse(bmi17v < 25 & case==1, 1, 0)) %>%
  mutate(nonlean_nafld_binary = ifelse(bmi17v >= 25 & case==1, 1, 0)) %>%
  mutate(lean_nafld_lean_control = case_when(bmi17v < 25 & case==1 ~ 1, bmi17v < 25 & case==0 ~ 0)) %>%
  mutate(nonlean_nafld_nonlean_control = case_when(bmi17v >= 25 & case==1 ~ 1, bmi17v >= 25 & case==0 ~ 0))

nafld_data_mbx <- nafld_data %>% select(all_of(mbx_list))

mbx_w_meta <- df_w_meta
mbx_df <- mbx.data.filtered
mbx_NAFLD_df<-mbx_w_meta %>% filter(cohort=="NAFLD" | case==0)

#Extended Data Figure 1A
mbx_NAFLD_perm <- adonis2(mbx_NAFLD_df[, colnames(mbx_NAFLD_df) %in% colnames(mbx_df)]~case + age + aheiv2010_15 + act17v + bmi17v + db17 + plate, data = mbx_NAFLD_df, permutations = 999, method = "bray", sqrt.dist = FALSE, add = FALSE, by = 'margin')
mbx_NAFLD_perm

#####Virome
virome_profile <- read.delim("MGXBAQLaVa_VGB_table.tsv",row.names=1) %>% t() %>% as.data.frame() %>% rownames_to_column() %>% rename(barcode_metagenomics = rowname)
virome_profile$barcode_metagenomics <- gsub("_Abundance.RPKs", "", virome_profile$barcode_metagenomics)

virome_profile <- virome_profile %>% column_to_rownames("barcode_metagenomics")
# Calculate total RPKs for each sample
virome_profile$total_rpks <- rowSums(virome_profile)

# Convert RPK values to relative abundance
rpks_df_relative <- virome_profile / virome_profile$total_rpks
virome_profile_new <- rpks_df_relative %>% select(-total_rpks) %>% rownames_to_column("barcode_metagenomics")

virome_list<-virome_profile_new %>% select(-barcode_metagenomics)

final_metadata <- read.delim('input/meta_df.tsv',row.names=1)
final_metadata <- final_metadata[!duplicated(final_metadata$barcode_metagenomics),]
df_w_meta <- left_join(virome_profile_new,final_metadata,by="barcode_metagenomics")

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
  mutate(obesity = case_when(bmi17v >= 30 ~ 1, bmi17v <30 ~ 0)) %>%
  mutate(lean = case_when(bmi17v < 25 ~ 1, bmi17v >= 25 ~ 0)) %>%
  mutate(lean_nafld_binary = ifelse(bmi17v < 25 & case==1, 1, 0)) %>%
  mutate(nonlean_nafld_binary = ifelse(bmi17v >= 25 & case==1, 1, 0)) %>%
  mutate(lean_nafld_lean_control = case_when(bmi17v < 25 & case==1 ~ 1, bmi17v < 25 & case==0 ~ 0)) %>%
  mutate(nonlean_nafld_nonlean_control = case_when(bmi17v >= 25 & case==1 ~ 1, bmi17v >= 25 & case==0 ~ 0)) %>%
  mutate(lean_vs_nonlean_case = case_when(bmi17v >= 25 & case==1 ~ 1, bmi17v < 25 & case==1 ~ 0)) #nonlean case is 1, lean case is 0

nafld_data_virome <- nafld_data %>% select(c(names(virome_list)))
unclassified_virome <- nafld_data_virome %>% select(contains("Unclassified"))

virome_w_meta <- df_w_meta
virome_df <- virome.data.filtered
virome_NAFLD_df<-virome_w_meta %>% filter(cohort=="NAFLD" | case==0)

#Figure 5A
virome_NAFLD_perm <- adonis2(virome_NAFLD_df[, colnames(virome_NAFLD_df) %in% colnames(virome_df)]~case + age + aheiv2010_15 + act17v + bmi17v + db17 + plate, data = virome_NAFLD_df, permutations = 999, method = "bray", sqrt.dist = FALSE, add = FALSE, by = 'margin')
virome_NAFLD_perm
