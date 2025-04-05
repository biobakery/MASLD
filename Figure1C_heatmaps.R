##################################################
#R program for creating Figure 1C
##################################################

library(dplyr)
library(tidyverse)
library(stringr)
library(readr)
library(Maaslin2)
library(ggplot2)
library(ggrepel)
library(pheatmap)
library(RColorBrewer)
library(viridis)
library(gridExtra)
library(grid)
library(lattice)

setwd("~/b2b")
###species###
unfiltered_species <- read.delim('input/species_v4.tsv',row.names=1) %>% t() %>% as.data.frame() %>% rownames_to_column() %>% rename(barcode_metagenomics = rowname)
final_metadata <- read.delim('input/meta_df.tsv',row.names=1)
final_metadata <- final_metadata[!duplicated(final_metadata$barcode_metagenomics),]
df_w_meta <- left_join(unfiltered_species,final_metadata,by="barcode_metagenomics")

species.data<-df_w_meta %>% select(starts_with('s__'))

nafld.data<-df_w_meta %>% filter(cohort=="NAFLD")
species.nafld.data<-nafld.data %>% select(starts_with('s__'))

avg.abundance.data<-as.data.frame(colMeans(species.nafld.data)) %>% rownames_to_column()
top5<-avg.abundance.data %>% top_n(5)
top5<-top5[order(top5$`colMeans(species.nafld.data)`,decreasing=TRUE), ] #sort
top5species<-top5$rowname

nafld_data<-df_w_meta %>% filter(cohort=="NAFLD" | case==0) %>%
  mutate(lean_vs_nonlean_case = case_when(bmi17v >= 25 & case==1 ~ 2, bmi17v < 25 & case==1 ~ 1, case==0 ~ 0)) #nonlean case is 2, lean case is 1

nafld_case<-nafld_data %>% filter(case==1)
nafld_control<-nafld_data %>% filter(case==0)
species.nafld.data <- nafld_data %>% column_to_rownames("alias_id") %>% select(starts_with('s__'))
species.nafld.data <- species.nafld.data / 100

species.nafld.data<-species.nafld.data[order(row.names(species.nafld.data)), ]

avg.abundance.data<-as.data.frame(colMeans(species.nafld.data)) %>% rownames_to_column()
top5<-avg.abundance.data %>% top_n(5)
top5<-top5[order(top5$`colMeans(species.nafld.data)`,decreasing=TRUE), ] #sort
top5species<-top5$rowname

selected_top5 <- species.nafld.data %>% t() %>% as.data.frame()
selected_top5_species <- selected_top5 %>% filter(row.names(selected_top5) %in% top5species) %>% as.matrix()
selected_top5_species <- selected_top5_species %>% t() %>% as.data.frame() %>% rownames_to_column()
casecontrol<-nafld_data %>% select(alias_id, case) %>% mutate(rowname=alias_id) %>% select(-alias_id)
groupcase <- left_join(selected_top5_species,casecontrol,by="rowname") %>% group_by(case) %>% arrange(-case) %>% ungroup()
casebar <- groupcase %>% select(rowname,case) %>% column_to_rownames("rowname")
selected_top5_species<-groupcase %>% column_to_rownames("rowname") %>% select(starts_with('s__')) %>% t()
selected_top5_species<-selected_top5_species[top5species, ] #sort
df.min <- (min(selected_top5_species[selected_top5_species > 0])/2) # assign 1/2 of minimum to 0 so we can log transform
log_nafld_data_species<-selected_top5_species
log_nafld_data_species[log_nafld_data_species == 0] <- df.min
log_nafld_data_species <- log10(log_nafld_data_species)
ann_colors = list(
  case = c("0" = "#999999", "1" = "#E69F00"))
# Customize the magma color palette to start from light gray
magma_colors <- viridis::viridis(100, option = "A")
magma_colors[1] <- "#D3D3D3"  # Change the first color to light gray
#without rownames and colnames
plot_species <- pheatmap(log_nafld_data_species, 
                         cluster_rows = FALSE, 
                         cluster_cols = FALSE, 
                         show_rownames = FALSE, 
                         show_colnames = FALSE, 
                         annotation_col = casebar, 
                         annotation_colors = ann_colors, 
                         legend = FALSE, 
                         color = magma_colors)
# Reorder the columns based on the clustered_species
column_order <- clustered_species$tree_col$order
log_nafld_data_species_reordered <- log_nafld_data_species[, column_order]
# Create the reordered heatmap
clustered_species_reordered <- pheatmap(log_nafld_data_species_reordered, 
                                        cluster_rows = FALSE, 
                                        cluster_cols = FALSE, 
                                        show_rownames = FALSE, 
                                        show_colnames = FALSE, 
                                        legend = FALSE, 
                                        color = magma_colors)

## Get annotation bar for nonlean, lean, control
#selected_top5 <- species.nafld.data %>% t() %>% as.data.frame()
#selected_top5_species <- selected_top5 %>% filter(row.names(selected_top5) %in% top5species) %>% as.matrix()
#selected_top5_species <- selected_top5_species %>% t() %>% as.data.frame() %>% rownames_to_column()
#casecontrol<-nafld_data %>% select(alias_id, lean_vs_nonlean_case) %>% mutate(rowname=alias_id) %>% select(-alias_id)
#groupcase <- left_join(selected_top5_species,casecontrol,by="rowname") %>% group_by(lean_vs_nonlean_case) %>% arrange(-lean_vs_nonlean_case) %>% ungroup()
#casebar <- groupcase %>% select(rowname,lean_vs_nonlean_case) %>% column_to_rownames("rowname")
#nonleancasebar <- groupcase %>% select(rowname,lean_vs_nonlean_case)
#selected_top5_species<-groupcase %>% column_to_rownames("rowname") %>% select(starts_with('s__')) %>% t()
#selected_top5_species<-selected_top5_species[top5species, ] #sort
#df.min <- (min(selected_top5_species[selected_top5_species > 0])/2) # assign 1/2 of minimum to 0 so we can log transform
#log_nafld_data_species<-selected_top5_species
#log_nafld_data_species[log_nafld_data_species == 0] <- df.min
#log_nafld_data_species <- log10(log_nafld_data_species)
#ann_colors = list(
#  lean_vs_nonlean_case = c("0" = "#999999", "1" = "blue", "2" = "red"))
#clustered_species<-pheatmap(log_nafld_data_species,cluster_rows=FALSE, cluster_cols=T, show_rownames = F, show_colnames = F, annotation_col=casebar, annotation_colors = ann_colors,legend=T,color = magma_colors)

###virome###
virome_file <- read.delim("~/b2b/viromeprofile_v3.tsv",row.names=1) %>% t() %>% as.data.frame() %>% rownames_to_column() %>% rename(barcode_metagenomics = rowname)
final_metadata <- read.delim('input/meta_df.tsv',row.names=1)
final_metadata <- final_metadata[!duplicated(final_metadata$barcode_metagenomics),]
df_w_meta <- left_join(virome_file,final_metadata,by="barcode_metagenomics")

virome.nafld.data<-df_w_meta %>% filter(cohort=="NAFLD" | case==0) %>% select(c(contains("VGB"),"alias_id")) %>% column_to_rownames("alias_id")
virome.nafld.data<-virome.nafld.data[order(row.names(virome.nafld.data)), ]

samples_not_in_virome <- species.nafld.data %>% filter(!(row.names(species.nafld.data) %in% row.names(virome.nafld.data))) %>% row.names()
samples_not_in_species <- virome.nafld.data %>% filter(!(row.names(virome.nafld.data) %in% row.names(species.nafld.data))) %>% row.names()

avg.abundance.data<-as.data.frame(colMeans(virome.nafld.data)) %>% rownames_to_column()
top5<-avg.abundance.data %>% top_n(5)
top5<-top5[order(top5$`colMeans(virome.nafld.data)`,decreasing=TRUE), ] #sort
top5virome<-top5$rowname

selected_top5 <- virome.nafld.data %>% t() %>% as.data.frame() %>% mutate(!!!setNames(rep(list(NA), length(samples_not_in_virome)), samples_not_in_virome))
selected_top5<-selected_top5[ , order(names(selected_top5))] #sort by column
selected_top5_virome <- selected_top5 %>% filter(row.names(selected_top5) %in% top5virome) %>% as.matrix()
selected_top5_virome <- selected_top5_virome %>% t() %>% as.data.frame() %>% rownames_to_column()
groupcase <- left_join(selected_top5_virome,casecontrol,by="rowname") %>% group_by(case) %>% arrange(-case) %>% ungroup()
casebar <- groupcase %>% select(rowname,case) %>% column_to_rownames("rowname")
selected_top5_virome<-groupcase %>% column_to_rownames("rowname") %>% select(contains("VGB")) %>% t()
selected_top5_virome<-selected_top5_virome[top5virome, ] #sort
df.min <- (min(selected_top5_virome[selected_top5_virome > 0],na.rm=T)/2) # assign 1/2 of minimum to 0 so we can log transform. Note: here, na.rm=T
log_nafld_data_virome<-selected_top5_virome
log_nafld_data_virome[log_nafld_data_virome == 0] <- df.min
log_nafld_data_virome <- log10(log_nafld_data_virome)

#without rownames and colnames
plot_virome<-pheatmap(log_nafld_data_virome,cluster_rows=FALSE, cluster_cols=FALSE, show_rownames = F, show_colnames = F, na_col="white",annotation_col=casebar, annotation_colors = ann_colors,legend=F)

column_order <- clustered_species$tree_col$order
log_nafld_data_virome_reordered <- log_nafld_data_virome[, column_order]
clustered_virome_reordered<-pheatmap(log_nafld_data_virome_reordered, cluster_rows = FALSE, cluster_cols = F, show_rownames = F, show_colnames = F, na_col = "white", legend=F,color = magma_colors)

###pathway###
pathway_file <- read.delim("pathabundance_unstrat.tsv",row.names=1) %>% t() %>% as.data.frame() %>% rownames_to_column() %>% rename(barcode_metagenomics = rowname)
pathway_file$barcode_metagenomics = gsub("_Abundance", "", pathway_file$barcode_metagenomics)
final_metadata <- read.delim('input/meta_df.tsv',row.names=1)
final_metadata <- final_metadata[!duplicated(final_metadata$barcode_metagenomics),]
df_w_meta <- left_join(pathway_file,final_metadata,by="barcode_metagenomics")

pathway.nafld.data<-df_w_meta %>% filter(cohort=="NAFLD" | case==0) %>% select(c(contains(": "),"alias_id")) %>% column_to_rownames("alias_id")
pathway.nafld.data<-pathway.nafld.data[order(row.names(pathway.nafld.data)), ]

samples_not_in_pathway <- species.nafld.data %>% filter(!(row.names(species.nafld.data) %in% row.names(pathway.nafld.data))) %>% row.names()
samples_not_in_species <- pathway.nafld.data %>% filter(!(row.names(pathway.nafld.data) %in% row.names(species.nafld.data))) %>% row.names()

avg.abundance.data<-as.data.frame(colMeans(pathway.nafld.data)) %>% rownames_to_column()
top5<-avg.abundance.data %>% top_n(5)
top5<-top5[order(top5$`colMeans(pathway.nafld.data)`,decreasing=TRUE), ] #sort
top5pathway<-top5$rowname

selected_top5 <- pathway.nafld.data %>% t() %>% as.data.frame() %>% mutate(!!!setNames(rep(list(NA), length(samples_not_in_pathway)), samples_not_in_pathway))
selected_top5<-selected_top5[ , order(names(selected_top5))] #sort by column
selected_top5_pathway <- selected_top5 %>% filter(row.names(selected_top5) %in% top5pathway) %>% as.matrix()
selected_top5_pathway <- selected_top5_pathway %>% t() %>% as.data.frame() %>% rownames_to_column()
groupcase <- left_join(selected_top5_pathway,casecontrol,by="rowname") %>% group_by(case) %>% arrange(-case) %>% ungroup()
casebar <- groupcase %>% select(rowname,case) %>% column_to_rownames("rowname")
selected_top5_pathway<-groupcase %>% column_to_rownames("rowname") %>% select(contains(": ")) %>% t()
selected_top5_pathway<-selected_top5_pathway[top5pathway, ] #sort
df.min <- (min(selected_top5_pathway[selected_top5_pathway > 0],na.rm=T)/2) # assign 1/2 of minimum to 0 so we can log transform. Note: here, na.rm=T
log_nafld_data_pathway<-selected_top5_pathway
log_nafld_data_pathway[log_nafld_data_pathway == 0] <- df.min
log_nafld_data_pathway <- log10(log_nafld_data_pathway)

#without rownames and colnames
log_nafld_data_pathway_only_five<-log_nafld_data_pathway[1:(nrow(log_nafld_data_pathway) - 2), ]
plot_pathway<-pheatmap(log_nafld_data_pathway_only_five,cluster_rows=FALSE, cluster_cols=FALSE, show_rownames = F, show_colnames = F, na_col="white",annotation_col=casebar, annotation_colors = ann_colors,legend=F)
log_nafld_data_pathway_only_five_reordered <- log_nafld_data_pathway_only_five[, column_order]
clustered_pathway_reordered<-pheatmap(log_nafld_data_pathway_only_five_reordered, cluster_rows = FALSE, cluster_cols = F, show_rownames = F, show_colnames = F, na_col = "white", legend=F,color = magma_colors)

###MTX pathway###
mtxpathway<-read_tsv("pathabundance_relab_nospecial.tsv")
names(mtxpathway) = gsub("_Abundance", "", names(mtxpathway))
mtx_pathway<-mtxpathway %>% column_to_rownames("# Pathway") %>% t() %>% as.data.frame()
samples <- mtx_pathway %>% rownames_to_column()
samplenames = samples$rowname
mtx_pathway_unstratified<-mtx_pathway[,!grepl("\\|", colnames(mtx_pathway)) ] %>% rownames_to_column("barcode_metagenomics")

final_metadata <- read.delim('input/meta_df.tsv',row.names=1)
final_metadata <- final_metadata[!duplicated(final_metadata$barcode_metagenomics),]
df_w_meta <- left_join(mtx_pathway_unstratified,final_metadata,by="barcode_metagenomics")

pathway.nafld.data<-df_w_meta %>% filter(cohort=="NAFLD" | case==0) %>% select(c(contains(": "),"alias_id")) %>% column_to_rownames("alias_id")
pathway.nafld.data<-pathway.nafld.data[order(row.names(pathway.nafld.data)), ]

samples_not_in_mtx <- species.nafld.data %>% filter(!(row.names(species.nafld.data) %in% row.names(pathway.nafld.data))) %>% row.names()
samples_not_in_species <- pathway.nafld.data %>% filter(!(row.names(pathway.nafld.data) %in% row.names(species.nafld.data))) %>% row.names()

avg.abundance.data<-as.data.frame(colMeans(pathway.nafld.data)) %>% rownames_to_column()
top5<-avg.abundance.data %>% top_n(5)
top5<-top5[order(top5$`colMeans(pathway.nafld.data)`,decreasing=TRUE), ] #sort
top5pathway<-top5$rowname

selected_top5 <- pathway.nafld.data %>% t() %>% as.data.frame() %>% mutate(!!!setNames(rep(list(NA), length(samples_not_in_mtx)), samples_not_in_mtx))
selected_top5<-selected_top5[ , order(names(selected_top5))] #sort by column
selected_top5_pathway <- selected_top5 %>% filter(row.names(selected_top5) %in% top5pathway) %>% as.matrix()
selected_top5_pathway <- selected_top5_pathway %>% t() %>% as.data.frame() %>% rownames_to_column()
groupcase <- left_join(selected_top5_pathway,casecontrol,by="rowname") %>% group_by(case) %>% arrange(-case) %>% ungroup()
casebar <- groupcase %>% select(rowname,case) %>% column_to_rownames("rowname")
selected_top5_pathway<-groupcase %>% column_to_rownames("rowname") %>% select(contains(": ")) %>% t()
selected_top5_pathway<-selected_top5_pathway[top5pathway, ] #sort
df.min <- (min(selected_top5_pathway[selected_top5_pathway > 0],na.rm=T)/2) # assign 1/2 of minimum to 0 so we can log transform. Note: here, na.rm=T
log_nafld_data_pathway_mtx<-selected_top5_pathway
log_nafld_data_pathway_mtx[log_nafld_data_pathway_mtx == 0] <- df.min
log_nafld_data_pathway_mtx <- log10(log_nafld_data_pathway_mtx)

#without rownames and colnames
plot_mtx<-pheatmap(log_nafld_data_pathway_mtx,cluster_rows=FALSE, cluster_cols=FALSE, show_rownames = F, show_colnames = F,na_col = "white",annotation_col=casebar, annotation_colors = ann_colors,legend=F)
log_nafld_data_pathway_mtx_reordered <- log_nafld_data_pathway_mtx[, column_order]
clustered_mtx_reordered<-pheatmap(log_nafld_data_pathway_mtx_reordered, cluster_rows = FALSE, cluster_cols = F, show_rownames = F, show_colnames = F, na_col = "white", legend=F,color = magma_colors)

###mbx###
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
nafld_data<-df_w_meta %>% filter(cohort=="NAFLD" | case==0) %>%
  mutate(obesity = case_when(bmi17v >= 30 ~ 1, bmi17v <30 ~ 0)) %>%
  mutate(lean = case_when(bmi17v < 25 ~ 1, bmi17v >= 25 ~ 0)) %>%
  mutate(lean_nafld_binary = ifelse(bmi17v < 25 & case==1, 1, 0)) %>%
  mutate(nonlean_nafld_binary = ifelse(bmi17v >= 25 & case==1, 1, 0)) %>%
  mutate(lean_nafld_lean_control = case_when(bmi17v < 25 & case==1 ~ 1, bmi17v < 25 & case==0 ~ 0)) %>%
  mutate(nonlean_nafld_nonlean_control = case_when(bmi17v >= 25 & case==1 ~ 1, bmi17v >= 25 & case==0 ~ 0)) %>%
  mutate(lean_vs_nonlean_case = case_when(bmi17v >= 25 & case==1 ~ 1, bmi17v < 25 & case==1 ~ 0)) #nonlean case is 1, lean case is 0

nafld_data_mbx <- nafld_data %>% select(all_of(mbx_list))
#took average of 10 technically duplicated samples
nafld_data_mbx_id <- nafld_data %>%
  group_by(alias_id) %>%
  summarize(across(all_of(mbx_list), mean)) %>%
  column_to_rownames("alias_id")

mbx.nafld.data<-nafld_data_mbx_id[order(row.names(nafld_data_mbx_id)), ]

samples_not_in_mbx <- species.nafld.data %>% filter(!(row.names(species.nafld.data) %in% row.names(mbx.nafld.data))) %>% row.names()
samples_not_in_species <- mbx.nafld.data %>% filter(!(row.names(mbx.nafld.data) %in% row.names(species.nafld.data))) %>% row.names()

avg.abundance.data<-as.data.frame(colMeans(mbx.nafld.data)) %>% rownames_to_column()
top5<-avg.abundance.data %>% top_n(5)
top5<-top5[order(top5$`colMeans(mbx.nafld.data)`,decreasing=TRUE), ] #sort
top5mbx<-top5$rowname

mbx.nafld.data.same.sample<-mbx.nafld.data %>% filter(!(row.names(mbx.nafld.data) %in% samples_not_in_species))
selected_top5 <- mbx.nafld.data.same.sample %>% t() %>% as.data.frame() %>% mutate(!!!setNames(rep(list(NA), length(samples_not_in_mbx)), samples_not_in_mbx))
selected_top5<-selected_top5[ , order(names(selected_top5))] #sort by column
selected_top5_mbx <- selected_top5 %>% filter(row.names(selected_top5) %in% top5mbx) %>% as.matrix()
selected_top5_mbx <- selected_top5_mbx %>% t() %>% as.data.frame() %>% rownames_to_column()
groupcase <- left_join(selected_top5_mbx,casecontrol,by="rowname") %>% group_by(case) %>% arrange(-case) %>% ungroup()
casebar <- groupcase %>% select(rowname,case) %>% column_to_rownames("rowname")
selected_top5_mbx<-groupcase %>% column_to_rownames("rowname") %>% select(-case) %>% t()
selected_top5_mbx<-selected_top5_mbx[top5mbx, ] #sort
df.min <- (min(selected_top5_mbx[selected_top5_mbx > 0],na.rm=T)/2) # assign 1/2 of minimum to 0 so we can log transform. Note: here, na.rm=T
log_nafld_data_mbx<-selected_top5_mbx
log_nafld_data_mbx[log_nafld_data_mbx == 0] <- df.min
log_nafld_data_mbx <- log2(log_nafld_data_mbx)

#without rownames and colnames
plot_mbx<-pheatmap(log_nafld_data_mbx,cluster_rows=FALSE, cluster_cols=FALSE, show_rownames = F, show_colnames = F,na_col = "white", annotation_col=casebar, annotation_colors = ann_colors,legend=F)
log_nafld_data_mbx_reordered <- log_nafld_data_mbx[, column_order]
clustered_mbx_reordered<-pheatmap(log_nafld_data_mbx_reordered, cluster_rows = FALSE, cluster_cols = F, show_rownames = F, show_colnames = F, na_col = "white", legend=F,color = magma_colors)

clustered_list=list()
clustered_list[['clustered_species_reordered']]=clustered_species_reordered[[4]]
clustered_list[['clustered_mbx_reordered']]=clustered_mbx_reordered[[4]]
clustered_list[['clustered_virome_reordered']]=clustered_virome_reordered[[4]]
clustered_list[['clustered_pathway_reordered']]=clustered_pathway_reordered[[4]]
clustered_list[['clustered_mtx_reordered']]=clustered_mtx_reordered[[4]]
grid.arrange(grobs=clustered_list,
             ncol=1,nrow=5) 
dev.off()
