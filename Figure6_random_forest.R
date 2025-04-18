##################################################
#R program for creating Figure 6
##################################################

library(dplyr)
library(tidyverse)
library(stringr)
library(readr)
library(ggplot2)
library(ggrepel)
library(ggpubr)
library(SIAMCAT)

################################################################################
## species
################################################################################
tax<-read_tsv("metaphlan_taxonomic_profiles.tsv")

#metadata
final_metadata <- read.delim('meta_df.tsv',row.names=1)
final_metadata <- final_metadata[!duplicated(final_metadata$barcode_metagenomics),]

#process metaphlan output
tax <- tax %>% rename("taxonomy"="# taxonomy")
names(tax) = gsub("_taxonomic_profile", "", names(tax))
tax_t<-tax[grepl("t__", tax$taxonomy), ] #only use the ones with t__ level
levels <- c("Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species", "Species2")
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

species.data<-df_w_meta %>% select(matches("SGB|EUK"))

oralspecies <- read_csv('oral_vs_gut.csv')

nafld_data<-df_w_meta %>% filter(cohort=="NAFLD" | case==0) %>%
  mutate(obesity = case_when(bmi17v >= 30 ~ 1, bmi17v <30 ~ 0)) %>%
  mutate(lean = case_when(bmi17v < 25 ~ 1, bmi17v >= 25 ~ 0)) %>%
  mutate(lean_nafld_binary = ifelse(bmi17v < 25 & case==1, 1, 0)) %>%
  mutate(nonlean_nafld_binary = ifelse(bmi17v >= 25 & case==1, 1, 0)) %>%
  mutate(lean_nafld_lean_control = case_when(bmi17v < 25 & case==1 ~ 1, bmi17v < 25 & case==0 ~ 0)) %>%
  mutate(nonlean_nafld_nonlean_control = case_when(bmi17v >= 25 & case==1 ~ 1, bmi17v >= 25 & case==0 ~ 0)) %>%
  mutate(lean_vs_nonlean_case = case_when(bmi17v >= 25 & case==1 ~ 1, bmi17v < 25 & case==1 ~ 0)) %>% #nonlean case is 1, lean case is 0
  mutate(nonlean_vs_control = case_when(bmi17v >= 25 & case==1 ~ 1, case==0 ~ 0)) %>%
  mutate(lean_vs_control = case_when(bmi17v < 25 & case==1 ~ 1, case==0 ~ 0))
nafld_data_species <- nafld_data %>% column_to_rownames("alias_id") %>% select(matches("SGB|EUK"))

################################################################################
## metabolites
################################################################################
unfiltered_mbx <- read.delim("annotated_metabolites_w_methods.tsv",row.names=1) %>% select(starts_with("X")) %>% t() %>% as.data.frame() %>% rownames_to_column() %>% rename(barcode_metabolomics = rowname)

###Normalize the mbx by methods.
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

final_metadata_mbx <- read.delim('meta_df.tsv',row.names=1)
final_metadata_mbx <- final_metadata_mbx[!duplicated(final_metadata_mbx$barcode_metabolomics),]
df_w_meta_mbx <- left_join(df_final_join,final_metadata_mbx,by="barcode_metabolomics")

mbx_list<-names(df_final)
mbx.data<-df_w_meta_mbx %>% select(all_of(mbx_list))

###########
###NAFLD###
###########
nafld_data_mbx<-df_w_meta_mbx %>% filter(cohort=="NAFLD" | case==0) %>%
  mutate(obesity = case_when(bmi17v >= 30 ~ 1, bmi17v <30 ~ 0)) %>%
  mutate(lean = case_when(bmi17v < 25 ~ 1, bmi17v >= 25 ~ 0)) %>%
  mutate(lean_nafld_binary = ifelse(bmi17v < 25 & case==1, 1, 0)) %>%
  mutate(nonlean_nafld_binary = ifelse(bmi17v >= 25 & case==1, 1, 0)) %>%
  mutate(lean_nafld_lean_control = case_when(bmi17v < 25 & case==1 ~ 1, bmi17v < 25 & case==0 ~ 0)) %>%
  mutate(nonlean_nafld_nonlean_control = case_when(bmi17v >= 25 & case==1 ~ 1, bmi17v >= 25 & case==0 ~ 0)) %>%
  mutate(lean_vs_nonlean_case = case_when(bmi17v >= 25 & case==1 ~ 1, bmi17v < 25 & case==1 ~ 0)) %>% #nonlean case is 1, lean case is 0
  mutate(nonlean_vs_control = case_when(bmi17v >= 25 & case==1 ~ 1, case==0 ~ 0)) %>%
  mutate(lean_vs_control = case_when(bmi17v < 25 & case==1 ~ 1, case==0 ~ 0))
#took average of 10 technically duplicated samples
nafld_data_mbx_id <- nafld_data_mbx %>%
  group_by(alias_id) %>%
  summarize(across(all_of(mbx_list), mean)) %>%
  column_to_rownames("alias_id")

################################################################################
## MTX pathways
################################################################################
mtxpathway<-read_tsv("pathabundance_relab_nospecial.tsv")
names(mtxpathway) = gsub("_Abundance", "", names(mtxpathway))
mtx_pathway<-mtxpathway %>% column_to_rownames("# Pathway") %>% t() %>% as.data.frame() 
samples <- mtx_pathway %>% rownames_to_column()
samplenames = samples$rowname 
mtx_pathway_unstratified<-mtx_pathway[,!grepl("\\|", colnames(mtx_pathway)) ] %>% rownames_to_column("barcode_metagenomics")

final_metadata_mtx <- read.delim('meta_df.tsv',row.names=1)
final_metadata_mtx <- final_metadata_mtx[!duplicated(final_metadata_mtx$barcode_metagenomics),]
df_w_meta_mtx <- left_join(mtx_pathway_unstratified,final_metadata_mtx,by="barcode_metagenomics") 

nafld_data_mtx_pathway<-df_w_meta_mtx %>% filter(cohort=="NAFLD" | case==0) %>% select(c(contains(": "),"alias_id")) %>% column_to_rownames("alias_id")
nafld_data_mtx_pathway<-nafld_data_mtx_pathway[order(row.names(nafld_data_mtx_pathway)), ]

################################################################################
## Virome
################################################################################
virome_profile <- read.delim("viromeprofile_v3.tsv",row.names=1) %>% t() %>% as.data.frame() %>% rownames_to_column() %>% rename(barcode_metagenomics = rowname)

virome_list<-virome_profile %>% select(-barcode_metagenomics)

final_metadata_vir <- read.delim('meta_df.tsv',row.names=1)
final_metadata_vir <- final_metadata_vir[!duplicated(final_metadata_vir$barcode_metagenomics),]
df_w_meta_vir <- left_join(virome_profile,final_metadata_vir,by="barcode_metagenomics")

nafld_data_virome<-df_w_meta_vir %>% filter(cohort=="NAFLD" | case==0) %>% select(c(names(virome_list),"alias_id")) %>% column_to_rownames("alias_id")
nafld_data_virome<-nafld_data_virome[order(row.names(nafld_data_virome)), ]

################################################################################
## ML MODELING
################################################################################

bugs.s.siam <- t(nafld_data_species)
mbx.s.siam <- t(nafld_data_mbx_id)

col_sums <- colSums(mbx.s.siam, na.rm = TRUE)
norm.mbx.s.siam <- t(apply(mbx.s.siam, 1, function(col) col / col_sums))

meta.siam <- nafld_data %>%
  select('alias_id', 'case', 'age', 'db17', 'bmi17v', 'act17v', 'aheiv2010_15') %>%
  column_to_rownames('alias_id')
# label cases vs. controls
meta.siam_label <- create.label(meta = meta.siam, label = 'case', case = '1')

# create a siamcat object
siam.all <- siamcat (feat = bugs.s.siam, label=meta.siam_label, meta = meta.siam)
siam.all.mbx <- siamcat (feat = norm.mbx.s.siam, label=meta.siam_label, meta = meta.siam)

show(siam.all)

# feature selection
siam.all <- filter.features(siam.all, filter.method='abundance', cutoff=0.001)
siam.all.mbx <- filter.features(siam.all.mbx, filter.method='abundance', cutoff=0.001)

siam.all.confounders <- check.confounders(
  siam.all,
  fn.plot = './confounder_plots.pdf',
  meta.in = c('age',
              'db17',
              'bmi17v',
              'act17v',
              'aheiv2010_15'),
  feature.type = 'filtered'
)

# data normalizing
siam.all <- normalize.features(siam.all, norm.method = 'log.std',
                               norm.param = list(log.n0=1e-06, sd.min.q=0))
siam.all.mbx <- normalize.features(siam.all.mbx, norm.method = 'log.std',
                                   norm.param = list(log.n0=1e-06, sd.min.q=0))
set.seed(1123)
# prepare cross-validation
siam.all <-  create.data.split(
  siam.all,
  num.folds = 5,
  num.resample = 2
)
set.seed(1123)
siam.all.mbx <-  create.data.split(
  siam.all.mbx,
  num.folds = 5,
  num.resample = 2
)

set.seed(1123)
siam.all <- train.model(
  siam.all,
  method = "randomForest"
)
set.seed(1123)
siam.all.mbx <- train.model(
  siam.all.mbx,
  method = "randomForest"
)

# get information about the model type
model_type(siam.all)

# access the models
models <- models(siam.all)
models[[1]]

# make predictions using the data-split and the models trained in previous step
siam.all <- make.predictions(siam.all)
siam.all.mbx <- make.predictions(siam.all.mbx)
siam.predict <- pred_matrix(siam.all)

# model evaluation and interpretation
siam.all <-  evaluate.predictions(siam.all)
siam.all.mbx <-  evaluate.predictions(siam.all.mbx)

# evaluation plot 
model.evaluation.plot(siam.all,
                      fn.plot = './eval_plot_RF.pdf')
model.evaluation.plot(siam.all.mbx,
                      fn.plot = './eval_plot_mbx_RF.pdf')

### mtx pathway
mtx.s.siam<-t(nafld_data_mtx_pathway)

# create a siamcat object
siam.all.mtx <- siamcat (feat = mtx.s.siam, label = meta.siam_label, meta = meta.siam)

# feature selection
siam.all.mtx <- filter.features(siam.all.mtx, filter.method='abundance', cutoff=0.001)

# data normalizing
siam.all.mtx <- normalize.features(siam.all.mtx, norm.method = 'log.std',
                                   norm.param = list(log.n0=1e-06, sd.min.q=0))

set.seed(1123)
siam.all.mtx <-  create.data.split(
  siam.all.mtx,
  num.folds = 5,
  num.resample = 2
)

set.seed(1123)
siam.all.mtx <- train.model(
  siam.all.mtx,
  method = "randomForest"
)

siam.all.mtx <- make.predictions(siam.all.mtx)

siam.all.mtx <-  evaluate.predictions(siam.all.mtx)

model.evaluation.plot(siam.all.mtx,
                      fn.plot = './eval_plot.mtx_RF.pdf')

### Virome
vir.s.siam<-t(nafld_data_virome)

# create a siamcat object
siam.all.vir <- siamcat (feat = vir.s.siam, label = meta.siam_label, meta = meta.siam)

# feature selection
siam.all.vir <- filter.features(siam.all.vir, filter.method='abundance', cutoff=0.001)

# data normalizing
siam.all.vir <- normalize.features(siam.all.vir, norm.method = 'log.std',
                                   norm.param = list(log.n0=1e-06, sd.min.q=0))
set.seed(1123)
siam.all.vir <-  create.data.split(
  siam.all.vir,
  num.folds = 5,
  num.resample = 2
)

set.seed(1123)
siam.all.vir <- train.model(
  siam.all.vir,
  method = "randomForest"
)

siam.all.vir <- make.predictions(siam.all.vir)

siam.all.vir <-  evaluate.predictions(siam.all.vir)

model.evaluation.plot(siam.all.vir,
                      fn.plot = './eval_plot.vir_RF.pdf')


########################################################
######## add metadata #########
########################################################
meta.siam <- meta.siam %>% mutate(db17 = ifelse(db17 ==0, 0, 1))

common_samples <- Reduce(intersect, lapply(list(bugs.s.siam, mtx.s.siam, norm.mbx.s.siam, vir.s.siam), colnames))
bugs.s.siam_common <- bugs.s.siam[, common_samples]
mtx.s.siam_common <- mtx.s.siam[, common_samples]
norm.mbx.s.siam_common <- norm.mbx.s.siam[, common_samples]
vir.s.siam_common <- vir.s.siam[, common_samples]
all_features <- rbind(bugs.s.siam_common, mtx.s.siam_common, norm.mbx.s.siam_common, vir.s.siam_common)

# create a siamcat object
siam.all.vl <- siamcat (feat = all_features, label = meta.siam_label, meta = meta.siam)

# feature selection
siam.all.vl <- filter.features(siam.all.vl, filter.method='abundance', cutoff=0.001)

# data normalizing
siam.all.vl <- normalize.features(siam.all.vl, norm.method = 'log.std',
                                  norm.param = list(log.n0=1e-06, sd.min.q=0))

set.seed(1123)
# prepare cross-validation
siam.all.vl <-  create.data.split(
  siam.all.vl,
  num.folds = 5,
  num.resample = 2
)

siam.all.vl <- add.meta.pred(siam.all.vl,
                             pred.names=c('age',
                                          'db17',
                                          'bmi17v',
                                          'act17v',
                                          'aheiv2010_15'))

set.seed(1123)
siam.all.vl <- train.model(
  siam.all.vl,
  method = "randomForest"
)

siam.all.vl <- make.predictions(siam.all.vl)

siam.all.vl <-  evaluate.predictions(siam.all.vl)

# model interpretation plot for all data + metadata
model.interpretation.plot(
  siam.all.vl,
  color.scheme = 'RdBu',
  fn.plot = './interpretation_plot_all_metadata_RF.pdf',
consens.thres = 0.001
)

# evaluation plot
model.evaluation.plot('Taxa' = siam.all,
                      'MBX' = siam.all.mbx,
                      'MTX pathways' = siam.all.mtx,
                      'Virome' = siam.all.vir,
                      'All + metadata' = siam.all.vl,
                      fn.plot = './eval_plot.full_RF.pdf')


############################
############################
######nonlean vs lean case
############################
############################
# label cases vs. controls
new.meta.siam <- nafld_data %>%
  select('alias_id', 'lean_vs_nonlean_case', 'age', 'db17', 'bmi17v', 'act17v', 'aheiv2010_15') %>%
  column_to_rownames('alias_id')
meta.siam_label <- create.label(meta = new.meta.siam, label = 'lean_vs_nonlean_case', case = '1')

# create a siamcat object
siam.all_nonleanlean <- siamcat (feat = bugs.s.siam, label=meta.siam_label, meta = new.meta.siam)
siam.all.mbx_nonleanlean <- siamcat (feat = norm.mbx.s.siam, label=meta.siam_label, meta = new.meta.siam)

show(siam.all_nonleanlean)

# feature selection
siam.all_nonleanlean <- filter.features(siam.all_nonleanlean, filter.method='abundance', cutoff=0.001)
siam.all.mbx_nonleanlean <- filter.features(siam.all.mbx_nonleanlean, filter.method='abundance', cutoff=0.001)

# data normalizing
siam.all_nonleanlean <- normalize.features(siam.all_nonleanlean, norm.method = 'log.std',
                               norm.param = list(log.n0=1e-06, sd.min.q=0))
siam.all.mbx_nonleanlean <- normalize.features(siam.all.mbx_nonleanlean, norm.method = 'log.std',
                                   norm.param = list(log.n0=1e-06, sd.min.q=0))
# prepare cross-validation
set.seed(1123)
siam.all_nonleanlean <-  create.data.split(
  siam.all_nonleanlean,
  num.folds = 5,
  num.resample = 2
)
set.seed(1123)
siam.all.mbx_nonleanlean <-  create.data.split(
  siam.all.mbx_nonleanlean,
  num.folds = 5,
  num.resample = 2
)
set.seed(1123)
siam.all_nonleanlean <- train.model(
  siam.all_nonleanlean,
  method = "randomForest"
)
set.seed(1123)
siam.all.mbx_nonleanlean <- train.model(
  siam.all.mbx_nonleanlean,
  method = "randomForest"
)

model_type(siam.all_nonleanlean)

# access the models
models <- models(siam.all_nonleanlean)
models[[1]]

# make predictions using the data-split and the models trained in previous step
siam.all_nonleanlean <- make.predictions(siam.all_nonleanlean)
siam.all.mbx_nonleanlean <- make.predictions(siam.all.mbx_nonleanlean)
siam.predict_nonleanlean <- pred_matrix(siam.all_nonleanlean)

# model evaluation and interpretation
siam.all_nonleanlean <-  evaluate.predictions(siam.all_nonleanlean)
siam.all.mbx_nonleanlean <-  evaluate.predictions(siam.all.mbx_nonleanlean)

### mtx pathway
mtx.s.siam<-t(nafld_data_mtx_pathway)

# create a siamcat object
siam.all.mtx_nonleanlean <- siamcat (feat = mtx.s.siam, label = meta.siam_label, meta = new.meta.siam)

# feature selection
siam.all.mtx_nonleanlean <- filter.features(siam.all.mtx_nonleanlean, filter.method='abundance', cutoff=0.001)

# data normalizing
siam.all.mtx_nonleanlean <- normalize.features(siam.all.mtx_nonleanlean, norm.method = 'log.std',
                                   norm.param = list(log.n0=1e-06, sd.min.q=0))

# prepare cross-validation
set.seed(1123)
siam.all.mtx_nonleanlean <-  create.data.split(
  siam.all.mtx_nonleanlean,
  num.folds = 5,
  num.resample = 2
)

# model training, default setting
set.seed(1123)
siam.all.mtx_nonleanlean <- train.model(
  siam.all.mtx_nonleanlean,
  method = "randomForest"
)

# make predictions using the data-split and the models trained in previous step
siam.all.mtx_nonleanlean <- make.predictions(siam.all.mtx_nonleanlean)

# model evaluation and interpretation
siam.all.mtx_nonleanlean <-  evaluate.predictions(siam.all.mtx_nonleanlean)

### Virome
vir.s.siam<-t(nafld_data_virome)

# create a siamcat object
siam.all.vir_nonleanlean <- siamcat (feat = vir.s.siam, label = meta.siam_label, meta = new.meta.siam)

# feature selection
siam.all.vir_nonleanlean <- filter.features(siam.all.vir_nonleanlean, filter.method='abundance', cutoff=0.001)

# data normalizing
siam.all.vir_nonleanlean <- normalize.features(siam.all.vir_nonleanlean, norm.method = 'log.std',
                                   norm.param = list(log.n0=1e-06, sd.min.q=0))

# prepare cross-validation
set.seed(1123)
siam.all.vir_nonleanlean <-  create.data.split(
  siam.all.vir_nonleanlean,
  num.folds = 5,
  num.resample = 2
)

# model training, default setting
set.seed(1123)
siam.all.vir_nonleanlean <- train.model(
  siam.all.vir_nonleanlean,
  method = "randomForest"
)

# make predictions using the data-split and the models trained in previous step
siam.all.vir_nonleanlean <- make.predictions(siam.all.vir_nonleanlean)

# model evaluation and interpretation
siam.all.vir_nonleanlean <-  evaluate.predictions(siam.all.vir_nonleanlean)


########################################################
######## add metadata #########
########################################################
new.meta.siam <- new.meta.siam %>% mutate(db17 = ifelse(db17 ==0, 0, 1))

#combining all datasets
common_samples <- Reduce(intersect, lapply(list(bugs.s.siam, mtx.s.siam, norm.mbx.s.siam, vir.s.siam), colnames))
# Subset each matrix to include only common samples
bugs.s.siam_common <- bugs.s.siam[, common_samples]
mtx.s.siam_common <- mtx.s.siam[, common_samples]
norm.mbx.s.siam_common <- norm.mbx.s.siam[, common_samples]
vir.s.siam_common <- vir.s.siam[, common_samples]
# Combine matrices using rbind
all_features <- rbind(bugs.s.siam_common, mtx.s.siam_common, norm.mbx.s.siam_common, vir.s.siam_common)

# create a siamcat object
siam.all.vl_nonleanlean <- siamcat (feat = all_features, label = meta.siam_label, meta = new.meta.siam)

# feature selection
siam.all.vl_nonleanlean <- filter.features(siam.all.vl_nonleanlean, filter.method='abundance', cutoff=0.001)

# data normalizing
siam.all.vl_nonleanlean <- normalize.features(siam.all.vl_nonleanlean, norm.method = 'log.std',
                                  norm.param = list(log.n0=1e-06, sd.min.q=0))

set.seed(1123)
# prepare cross-validation
siam.all.vl_nonleanlean <-  create.data.split(
  siam.all.vl_nonleanlean,
  num.folds = 5,
  num.resample = 2
)

siam.all.vl_nonleanlean <- add.meta.pred(siam.all.vl_nonleanlean,
                             pred.names=c('age',
                                          'db17',
                                          #'bmi17v',
                                          'act17v',
                                          'aheiv2010_15'))

# model training, default setting
set.seed(1123)
siam.all.vl_nonleanlean <- train.model(
  siam.all.vl_nonleanlean,
  method = "randomForest"
)

# make predictions using the data-split and the models trained in previous step
siam.all.vl_nonleanlean <- make.predictions(siam.all.vl_nonleanlean)

# model evaluation and interpretation
siam.all.vl_nonleanlean <-  evaluate.predictions(siam.all.vl_nonleanlean)

# evaluation plot for non-lean case vs. lean case
model.evaluation.plot('Taxa' = siam.all_nonleanlean,
                      'MBX' = siam.all.mbx_nonleanlean,
                      'MTX pathways' = siam.all.mtx_nonleanlean,
                      'Virome' = siam.all.vir_nonleanlean,
                      fn.plot = './eval_plot.full.nonleancase_leancase_RF.pdf')


############################
############################
######nonlean vs control
############################
############################
# label cases vs. controls
new.meta.siam <- nafld_data %>%
  select('alias_id', 'nonlean_vs_control', 'age', 'db17', 'bmi17v', 'act17v', 'aheiv2010_15') %>%
  column_to_rownames('alias_id')
meta.siam_label <- create.label(meta = new.meta.siam, label = 'nonlean_vs_control', case = '1')

siam.all_nonleancontrol <- siamcat (feat = bugs.s.siam, label=meta.siam_label, meta = new.meta.siam)
siam.all.mbx_nonleancontrol <- siamcat (feat = norm.mbx.s.siam, label=meta.siam_label, meta = new.meta.siam)

show(siam.all_nonleancontrol)

siam.all_nonleancontrol <- filter.features(siam.all_nonleancontrol, filter.method='abundance', cutoff=0.001)
siam.all.mbx_nonleancontrol <- filter.features(siam.all.mbx_nonleancontrol, filter.method='abundance', cutoff=0.001)

siam.all_nonleancontrol <- normalize.features(siam.all_nonleancontrol, norm.method = 'log.std',
                               norm.param = list(log.n0=1e-06, sd.min.q=0))
siam.all.mbx_nonleancontrol <- normalize.features(siam.all.mbx_nonleancontrol, norm.method = 'log.std',
                                   norm.param = list(log.n0=1e-06, sd.min.q=0))
set.seed(1123)
siam.all_nonleancontrol <-  create.data.split(
  siam.all_nonleancontrol,
  num.folds = 5,
  num.resample = 2
)
set.seed(1123)
siam.all.mbx_nonleancontrol <-  create.data.split(
  siam.all.mbx_nonleancontrol,
  num.folds = 5,
  num.resample = 2
)
set.seed(1123)
siam.all_nonleancontrol <- train.model(
  siam.all_nonleancontrol,
  method = "randomForest"
)
set.seed(1123)
siam.all.mbx_nonleancontrol <- train.model(
  siam.all.mbx_nonleancontrol,
  method = "randomForest"
)

model_type(siam.all_nonleancontrol)

models <- models(siam.all_nonleancontrol)
models[[1]]

siam.all_nonleancontrol <- make.predictions(siam.all_nonleancontrol)
siam.all.mbx_nonleancontrol <- make.predictions(siam.all.mbx_nonleancontrol)
siam.predict_nonleancontrol <- pred_matrix(siam.all_nonleancontrol)

siam.all_nonleancontrol <-  evaluate.predictions(siam.all_nonleancontrol)
siam.all.mbx_nonleancontrol <-  evaluate.predictions(siam.all.mbx_nonleancontrol)

### mtx pathway
mtx.s.siam<-t(nafld_data_mtx_pathway)

siam.all.mtx_nonleancontrol <- siamcat (feat = mtx.s.siam, label = meta.siam_label, meta = new.meta.siam)
siam.all.mtx_nonleancontrol <- filter.features(siam.all.mtx_nonleancontrol, filter.method='abundance', cutoff=0.001)
siam.all.mtx_nonleancontrol <- normalize.features(siam.all.mtx_nonleancontrol, norm.method = 'log.std',
                                   norm.param = list(log.n0=1e-06, sd.min.q=0))

set.seed(1123)
siam.all.mtx_nonleancontrol <-  create.data.split(
  siam.all.mtx_nonleancontrol,
  num.folds = 5,
  num.resample = 2
)

set.seed(1123)
siam.all.mtx_nonleancontrol <- train.model(
  siam.all.mtx_nonleancontrol,
  method = "randomForest"
)

siam.all.mtx_nonleancontrol <- make.predictions(siam.all.mtx_nonleancontrol)
siam.all.mtx_nonleancontrol <-  evaluate.predictions(siam.all.mtx_nonleancontrol)

### Virome
vir.s.siam<-t(nafld_data_virome)

siam.all.vir_nonleancontrol <- siamcat (feat = vir.s.siam, label = meta.siam_label, meta = new.meta.siam)
siam.all.vir_nonleancontrol <- filter.features(siam.all.vir_nonleancontrol, filter.method='abundance', cutoff=0.001)
siam.all.vir_nonleancontrol <- normalize.features(siam.all.vir_nonleancontrol, norm.method = 'log.std',
                                   norm.param = list(log.n0=1e-06, sd.min.q=0))
set.seed(1123)
siam.all.vir_nonleancontrol <-  create.data.split(
  siam.all.vir_nonleancontrol,
  num.folds = 5,
  num.resample = 2
)

set.seed(1123)
siam.all.vir_nonleancontrol <- train.model(
  siam.all.vir_nonleancontrol,
  method = "randomForest"
)

siam.all.vir_nonleancontrol <- make.predictions(siam.all.vir_nonleancontrol)
siam.all.vir_nonleancontrol <-  evaluate.predictions(siam.all.vir_nonleancontrol)

########################################################
######## add metadata #########
########################################################
new.meta.siam <- new.meta.siam %>% mutate(db17 = ifelse(db17 ==0, 0, 1))

common_samples <- Reduce(intersect, lapply(list(bugs.s.siam, mtx.s.siam, norm.mbx.s.siam, vir.s.siam), colnames))
bugs.s.siam_common <- bugs.s.siam[, common_samples]
mtx.s.siam_common <- mtx.s.siam[, common_samples]
norm.mbx.s.siam_common <- norm.mbx.s.siam[, common_samples]
vir.s.siam_common <- vir.s.siam[, common_samples]
all_features <- rbind(bugs.s.siam_common, mtx.s.siam_common, norm.mbx.s.siam_common, vir.s.siam_common)

siam.all.vl_nonleancontrol <- siamcat (feat = all_features, label = meta.siam_label, meta = new.meta.siam)
siam.all.vl_nonleancontrol <- filter.features(siam.all.vl_nonleancontrol, filter.method='abundance', cutoff=0.001)
siam.all.vl_nonleancontrol <- normalize.features(siam.all.vl_nonleancontrol, norm.method = 'log.std',
                                  norm.param = list(log.n0=1e-06, sd.min.q=0))

set.seed(1123)
siam.all.vl_nonleancontrol <-  create.data.split(
  siam.all.vl_nonleancontrol,
  num.folds = 5,
  num.resample = 2
)

siam.all.vl_nonleancontrol <- add.meta.pred(siam.all.vl_nonleancontrol,
                             pred.names=c('age',
                                          'db17',
                                          #'bmi17v',
                                          'act17v',
                                          'aheiv2010_15'))

set.seed(1123)
siam.all.vl_nonleancontrol <- train.model(
  siam.all.vl_nonleancontrol,
  method = "randomForest"
)

siam.all.vl_nonleancontrol <- make.predictions(siam.all.vl_nonleancontrol)
siam.all.vl_nonleancontrol <-  evaluate.predictions(siam.all.vl_nonleancontrol)

# evaluation plot for non-lean case vs. control
model.evaluation.plot('Taxa' = siam.all_nonleancontrol,
                      'MBX' = siam.all.mbx_nonleancontrol,
                      'MTX pathways' = siam.all.mtx_nonleancontrol,
                      'Virome' = siam.all.vir_nonleancontrol,
                      fn.plot = './eval_plot.full.nonleancase_control_RF.pdf')

############################
############################
######lean vs control
############################
############################
# label cases vs. controls
new.meta.siam <- nafld_data %>%
  select('alias_id', 'lean_vs_control', 'age', 'db17', 'bmi17v', 'act17v', 'aheiv2010_15') %>%
  column_to_rownames('alias_id')
meta.siam_label <- create.label(meta = new.meta.siam, label = 'lean_vs_control', case = '1')

siam.all_leancontrol <- siamcat (feat = bugs.s.siam, label=meta.siam_label, meta = new.meta.siam)
siam.all.mbx_leancontrol <- siamcat (feat = norm.mbx.s.siam, label=meta.siam_label, meta = new.meta.siam)

show(siam.all_leancontrol)

siam.all_leancontrol <- filter.features(siam.all_leancontrol, filter.method='abundance', cutoff=0.001)
siam.all.mbx_leancontrol <- filter.features(siam.all.mbx_leancontrol, filter.method='abundance', cutoff=0.001)

siam.all_leancontrol <- normalize.features(siam.all_leancontrol, norm.method = 'log.std',
                               norm.param = list(log.n0=1e-06, sd.min.q=0))
siam.all.mbx_leancontrol <- normalize.features(siam.all.mbx_leancontrol, norm.method = 'log.std',
                                   norm.param = list(log.n0=1e-06, sd.min.q=0))
set.seed(1123)
siam.all_leancontrol <-  create.data.split(
  siam.all_leancontrol,
  num.folds = 5,
  num.resample = 2
)
set.seed(1123)
siam.all.mbx_leancontrol <-  create.data.split(
  siam.all.mbx_leancontrol,
  num.folds = 5,
  num.resample = 2
)
set.seed(1123)
siam.all_leancontrol <- train.model(
  siam.all_leancontrol,
  method = "randomForest"
)
set.seed(1123)
siam.all.mbx_leancontrol <- train.model(
  siam.all.mbx_leancontrol,
  method = "randomForest"
)

model_type(siam.all_leancontrol)

models <- models(siam.all_leancontrol)
models[[1]]

siam.all_leancontrol <- make.predictions(siam.all_leancontrol)
siam.all.mbx_leancontrol <- make.predictions(siam.all.mbx_leancontrol)
siam.predict_leancontrol <- pred_matrix(siam.all_leancontrol)

siam.all_leancontrol <-  evaluate.predictions(siam.all_leancontrol)
siam.all.mbx_leancontrol <-  evaluate.predictions(siam.all.mbx_leancontrol)

### mtx pathway
mtx.s.siam<-t(nafld_data_mtx_pathway)

siam.all.mtx_leancontrol <- siamcat (feat = mtx.s.siam, label = meta.siam_label, meta = new.meta.siam)
siam.all.mtx_leancontrol <- filter.features(siam.all.mtx_leancontrol, filter.method='abundance', cutoff=0.001)
siam.all.mtx_leancontrol <- normalize.features(siam.all.mtx_leancontrol, norm.method = 'log.std',
                                   norm.param = list(log.n0=1e-06, sd.min.q=0))
set.seed(1123)
siam.all.mtx_leancontrol <-  create.data.split(
  siam.all.mtx_leancontrol,
  num.folds = 5,
  num.resample = 2
)

set.seed(1123)
siam.all.mtx_leancontrol <- train.model(
  siam.all.mtx_leancontrol,
  method = "randomForest"
)

siam.all.mtx_leancontrol <- make.predictions(siam.all.mtx_leancontrol)
siam.all.mtx_leancontrol <-  evaluate.predictions(siam.all.mtx_leancontrol)

### Virome
vir.s.siam<-t(nafld_data_virome)

siam.all.vir_leancontrol <- siamcat (feat = vir.s.siam, label = meta.siam_label, meta = new.meta.siam)
siam.all.vir_leancontrol <- filter.features(siam.all.vir_leancontrol, filter.method='abundance', cutoff=0.001)
siam.all.vir_leancontrol <- normalize.features(siam.all.vir_leancontrol, norm.method = 'log.std',
                                   norm.param = list(log.n0=1e-06, sd.min.q=0))
set.seed(1123)
siam.all.vir_leancontrol <-  create.data.split(
  siam.all.vir_leancontrol,
  num.folds = 5,
  num.resample = 2
)

set.seed(1123)
siam.all.vir_leancontrol <- train.model(
  siam.all.vir_leancontrol,
  method = "randomForest"
)

siam.all.vir_leancontrol <- make.predictions(siam.all.vir_leancontrol)
siam.all.vir_leancontrol <-  evaluate.predictions(siam.all.vir_leancontrol)

########################################################
######## add metadata #########
########################################################
new.meta.siam <- new.meta.siam %>% mutate(db17 = ifelse(db17 ==0, 0, 1))

common_samples <- Reduce(intersect, lapply(list(bugs.s.siam, mtx.s.siam, norm.mbx.s.siam, vir.s.siam), colnames))
bugs.s.siam_common <- bugs.s.siam[, common_samples]
mtx.s.siam_common <- mtx.s.siam[, common_samples]
norm.mbx.s.siam_common <- norm.mbx.s.siam[, common_samples]
vir.s.siam_common <- vir.s.siam[, common_samples]
all_features <- rbind(bugs.s.siam_common, mtx.s.siam_common, norm.mbx.s.siam_common, vir.s.siam_common)

siam.all.vl_leancontrol <- siamcat (feat = all_features, label = meta.siam_label, meta = new.meta.siam)
siam.all.vl_leancontrol <- filter.features(siam.all.vl_leancontrol, filter.method='abundance', cutoff=0.001)
siam.all.vl_leancontrol <- normalize.features(siam.all.vl_leancontrol, norm.method = 'log.std',
                                  norm.param = list(log.n0=1e-06, sd.min.q=0))

set.seed(1123)
siam.all.vl_leancontrol <-  create.data.split(
  siam.all.vl_leancontrol,
  num.folds = 5,
  num.resample = 2
)

siam.all.vl_leancontrol <- add.meta.pred(siam.all.vl_leancontrol,
                             pred.names=c('age',
                                          'db17',
                                          #'bmi17v',
                                          'act17v',
                                          'aheiv2010_15'))

set.seed(1123)
siam.all.vl_leancontrol <- train.model(
  siam.all.vl_leancontrol,
  method = "randomForest"
)

siam.all.vl_leancontrol <- make.predictions(siam.all.vl_leancontrol)
siam.all.vl_leancontrol <-  evaluate.predictions(siam.all.vl_leancontrol)

# evaluation plot for lean case vs. control
model.evaluation.plot('Taxa' = siam.all_leancontrol,
                      'MBX' = siam.all.mbx_leancontrol,
                      'MTX pathways' = siam.all.mtx_leancontrol,
                      'Virome' = siam.all.vir_leancontrol,
                      fn.plot = './eval_plot.full.leancase_control_RF.pdf')
