##################################################
#R program for creating alpha diversity plots
#Extended Data Figure 1BC
##################################################

library(vegan)
library(dplyr)
library(tidyverse)
library(stringr)
library(readr)
library(Maaslin2)
library(ggplot2)
library(ggrepel)

setwd("~/b2b")

unfiltered_species <- read.delim('input/species_v4.tsv',row.names=1) %>% t() %>% as.data.frame() %>% rownames_to_column() %>% rename(barcode_metagenomics = rowname) 
final_metadata <- read.delim('input/meta_df.tsv',row.names=1)
final_metadata <- final_metadata[!duplicated(final_metadata$barcode_metagenomics),]
df_w_meta <- left_join(unfiltered_species,final_metadata,by="barcode_metagenomics") 

nafld_data<-df_w_meta %>% filter(cohort=="NAFLD" | case==0) %>% 
  mutate(lean_nonlean_control = case_when(bmi17v >= 25 & case==1 ~ 'Non-leanMASLD', bmi17v <25 & case==1 ~ 'LeanMASLD', case==0 ~ 'Controls')) 

################################################################################
## shannon index
################################################################################

alpha <- nafld_data %>% column_to_rownames("barcode_metagenomics") %>% select(starts_with('s__')) 
alpha <- alpha/100
alpha2 <- as.data.frame(vegan::diversity(alpha, index = "shannon")) %>%
  select(shannon = `vegan::diversity(alpha, index = \"shannon\")`) %>% rownames_to_column("barcode_metagenomics")

alpha3 <- nafld_data %>%
  left_join(alpha2, by = "barcode_metagenomics")

###Extended Data Figure 1B
library(ggsignif)
ggplot(data = alpha3,
       aes(x = as.factor(case), y = shannon, fill = as.factor(case))) +
  scale_fill_manual(values = c("0" = "#999999", "1" = "#E69F00")) +
  geom_boxplot(notch = FALSE) +
  geom_jitter(shape = 16, position = position_jitter(0.1)) +
  theme_classic(base_size = 18) +
  theme(legend.position = "none") +
  ylab("Alpha diversity\n(Shannon Index)") +
  geom_signif(comparisons = list(c("0", "1")),
              test=wilcox.test,
              tip_length = 0,
              map_signif_level = TRUE) 

###Extended Data Figure 1C
# nonlean case vs lean case vs control
ggplot(data = alpha3, aes(x = as.factor(lean_nonlean_control), y = shannon, fill = as.factor(lean_nonlean_control))) +
  scale_fill_manual(values = c("Non-leanMASLD" = "#FF0000", "LeanMASLD" = "#0000FF", "Controls" = "#999999")) +
  geom_boxplot(notch = FALSE) +
  geom_jitter(shape = 16, position = position_jitter(0.1)) +
  theme_classic(base_size = 18) +
  theme(legend.position = "none") +
  xlab("Group") +
  ylab("Alpha diversity\n(Shannon Index)") +
  geom_signif(comparisons = list(c("Non-leanMASLD", "LeanMASLD"), c("LeanMASLD", "Controls"), c("Non-leanMASLD", "Controls")),
              tip_length = 0.02,
              y_position = c(4.8, 4.8, 5), 
              map_signif_level = TRUE)
