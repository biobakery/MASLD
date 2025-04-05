##################################################
#R program for creating Figure 5
##################################################

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
library(gridExtra)
library(vegan)

setwd("~/b2b")
virome_profile <- read.delim("MGXBAQLaVa_VGB_table.tsv",row.names=1) %>% t() %>% as.data.frame() %>% rownames_to_column() %>% rename(barcode_metagenomics = rowname)
virome_profile$barcode_metagenomics <- gsub("_Abundance.RPKs", "", virome_profile$barcode_metagenomics)

virome_profile <- virome_profile %>% column_to_rownames("barcode_metagenomics")
virome_profile$total_rpks <- rowSums(virome_profile)

# Convert RPK values to relative abundance
rpks_df_relative <- virome_profile / virome_profile$total_rpks
virome_profile_new <- rpks_df_relative %>% select(-total_rpks) %>% rownames_to_column("barcode_metagenomics")

virome_list<-virome_profile_new %>% select(-barcode_metagenomics)

final_metadata_vir <- read.delim('~/b2b/input/meta_df.tsv',row.names=1)
final_metadata_vir <- final_metadata_vir[!duplicated(final_metadata_vir$barcode_metagenomics),]
df_w_meta_vir <- left_join(virome_profile_new,final_metadata_vir,by="barcode_metagenomics") %>%
  distinct(barcode_metagenomics, .keep_all = TRUE)

alias_key<-final_metadata_vir %>% select(alias_id, barcode_metagenomics)

VGBnames <- read.delim("VGB_taxonomy.tsv",row.names=1) %>% as.data.frame() %>% rownames_to_column() %>% rename(VGB = rowname)
indices <- which(VGBnames$Reference.Species == "")
VGBnames$Reference.Species[indices] <- VGBnames$VGB[indices]
VGBnames$fullname <- ifelse(VGBnames$VGB == VGBnames$Reference.Species,
                            paste("Unclassified species", VGBnames$VGB),
                            paste(VGBnames$Reference.Species, VGBnames$VGB))
VGBnames_list <- VGBnames %>% select(VGB,fullname)

###########
###NAFLD###
###########
nafld_data<-df_w_meta_vir %>% filter(cohort=="NAFLD" | case==0) %>%
  mutate(lean_nonlean_control = case_when(bmi17v >= 25 & case==1 ~ 'Non-leanMASLD', bmi17v <25 & case==1 ~ 'LeanMASLD', case==0 ~ 'Controls')) %>%
  mutate(lean_nafld_lean_control = case_when(bmi17v < 25 & case==1 ~ 1, bmi17v < 25 & case==0 ~ 0)) %>%
  mutate(nonlean_nafld_nonlean_control = case_when(bmi17v >= 25 & case==1 ~ 1, bmi17v >= 25 & case==0 ~ 0)) %>%
  mutate(lean_vs_nonlean_case = case_when(bmi17v >= 25 & case==1 ~ 1, bmi17v < 25 & case==1 ~ 0)) #nonlean case is 1, lean case is 0

nafld_data_species <- nafld_data %>% select(c(names(virome_list)))

# assign 1/2 of minimum to 0 so we can log transform
df.min <- (min(nafld_data_species[nafld_data_species > 0])/2)
log_nafld_data_species<-nafld_data_species
log_nafld_data_species[log_nafld_data_species == 0] <- df.min
log_nafld_data_species <- log10(log_nafld_data_species)

###Figure 5B
###alpha diversity for virome
alpha <- nafld_data %>% column_to_rownames("barcode_metagenomics") %>% select(contains("VGB"))
alpha2 <- as.data.frame(vegan::diversity(alpha, index = "shannon")) %>%
  select(shannon = `vegan::diversity(alpha, index = \"shannon\")`) %>% rownames_to_column("barcode_metagenomics")

alpha3 <- nafld_data %>%
  left_join(alpha2, by = "barcode_metagenomics")

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

###Extended Data Figure 5A
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
              y_position = c(5.8, 5.8, 6),
              map_signif_level = TRUE)

#adjusted for age, most recent diabetes, cumulative avg of physical activity, cum avg of bmi, cum avg of AHEI
fit_data <- Maaslin2(
  input_data = nafld_data %>% select(c(names(virome_list))), 
  input_metadata = nafld_data %>% select(!c(names(virome_list))), #metadata
  output="output_virome_v0.3",
  normalization = "none", #our data already normalized
  transform = "LOG",
  analysis_method = "LM",
  max_significance = 0.20, # q-value threshold for significance. default is 0.25
  random_effects = NULL,
  fixed_effects = c('case', 'age', 'db17', 'act17v', 'bmi17v', 'aheiv2010_15'),
  correction = "BH",
  standardize = TRUE,
  cores = 30,
  plot_heatmap = TRUE,
  plot_scatter = TRUE,
  heatmap_first_n = 50)

###Figure 5C
###Volcano plot
results <- fit_data[["results"]]

results.nafld <- results %>%
  filter(metadata == "case") %>%
  mutate(adjlog2 = exp(coef)) %>%
  mutate(color = case_when(
    # changed this to 0 from 2
    coef > 0.25 & qval <= 0.2 ~ "Up",
    coef < -0.25 & qval <=0.2 ~ "Down",
    TRUE ~ "Stable")) %>%
  mutate(color2 = case_when(
    adjlog2 > 1 & qval <= 0.25 ~ "Up",
    adjlog2 < -1 & qval <=0.25 ~ "Down",
    TRUE ~ "Stable")
  )

# label top bugs
top <- 5

results_with_names <- merge(results.nafld, VGBnames_list, by.x = "feature", by.y = "VGB", all.x = TRUE)
top_bugs <- bind_rows(
  results_with_names %>%
    filter(color == 'Up') %>%
    arrange(qval, desc(abs(coef))) %>%
    head(top),
  results_with_names %>%
    filter(color == 'Down') %>%
    arrange(qval, desc(abs(coef))) %>%
    head(top)
)

options(ggrepel.max.overlaps = Inf)

ggplot(results_with_names, aes(x = coef,
                               y = -log(qval, 10),
                               colour= color)) + 
  geom_point(alpha=0.4, size=3.5) +
  xlab(expression("B-coefficient")) +
  ylab(expression("-log"[10] * "(FDR p-value)")) +
  xlim(c(-1, 1)) +
  scale_color_manual(values=c("blue", "grey","red")) +
  geom_vline(xintercept=c(-0.25, 0.25),lty=4, col="black",lwd=0.8) +
  geom_hline(yintercept = -log10(0.2),lty=4, col="black",lwd=0.8) +
  theme_classic(base_size = 18)+
  theme(plot.title = element_text(hjust = 0.5),
        legend.position="right",
        legend.title = element_blank()) +
  geom_label_repel(data = top_bugs,
                   mapping = aes(coef, -log(qval, 10),
                                 label = fullname),
                   show.legend = FALSE,
                   size = 6,
                   min.segment.length = 0.2,
                   color = "black",
                   fontface='bold')

###Figure 5A
nafld_data_meta <- nafld_data %>% select(!c(names(virome_list))) %>% rownames_to_column("sample_id")

bray <- nafld_data_species %>% vegdist(., "bray")
pc = capscale(bray~1, comm = nafld_data_species)
pc.summary<-summary(pc)

pcl.bray <- as.data.frame(pc$CA$u) %>%
  select(MDS1, MDS2)

pcl.bray <- pcl.bray %>%
  rownames_to_column("sample_id") %>%
  inner_join(nafld_data_meta, by = "sample_id") %>%
  select(colnames(nafld_data_meta), everything())

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

###Figure 5D
#lean nafld vs lean control
fit_data_lean <- Maaslin2(
  input_data = nafld_data %>% select(c(names(virome_list))), 
  input_metadata = nafld_data %>% select(!c(names(virome_list))), #metadata
  output="output_virome_v0.3/lean",
  normalization = "none",
  transform = "LOG",
  analysis_method = "LM",
  max_significance = 0.2, # q-value threshold for significance. default is 0.25
  random_effects = NULL,
  fixed_effects = c('lean_nafld_lean_control', 'age', 'db17', 'act17v', 'aheiv2010_15'),
  correction = "BH",
  standardize = TRUE,
  cores = 30,
  plot_heatmap = TRUE,
  plot_scatter = TRUE,
  heatmap_first_n = 50)

results <- fit_data_lean[["results"]]

leanresults <- read_tsv('output_virome_v0.3/lean/all_results.tsv') %>% filter(metadata=="lean_nafld_lean_control")
leanresults_sig <- read_tsv('output_virome_v0.3/lean/significant_results.tsv') %>% filter(metadata=="lean_nafld_lean_control")
leanresults_sig_with_names <- merge(leanresults_sig, VGBnames_list, by.x = "feature", by.y = "VGB", all.x = TRUE)

#nonlean nafld vs nonlean control
fit_data_nonlean <- Maaslin2(
  input_data = nafld_data %>% select(c(names(virome_list))),
  input_metadata = nafld_data %>% select(!c(names(virome_list))), #metadata
  output="output_virome_v0.3/nonlean",
  normalization = "none",
  transform = "LOG",
  analysis_method = "LM",
  max_significance = 0.20, # q-value threshold for significance. default is 0.25
  random_effects = NULL,
  fixed_effects = c('nonlean_nafld_nonlean_control', 'age', 'db17', 'act17v', 'aheiv2010_15'),
  correction = "BH",
  standardize = TRUE,
  cores = 30,
  plot_heatmap = TRUE,
  plot_scatter = TRUE,
  heatmap_first_n = 50)

results <- fit_data_nonlean[["results"]]

nonleanresults <- read_tsv('output_virome_v0.3/nonlean/all_results.tsv') %>% filter(metadata=="nonlean_nafld_nonlean_control")
nonleanresults_sig <- read_tsv('output_virome_v0.3/nonlean/significant_results.tsv') %>% filter(metadata=="nonlean_nafld_nonlean_control")
nonleanresults_sig_with_names <- merge(nonleanresults_sig, VGBnames_list, by.x = "feature", by.y = "VGB", all.x = TRUE)

#scatter plot for lean and nonlean
lean_and_nonlean<-full_join(leanresults,nonleanresults,by="feature")
lean_and_nonlean <- merge(lean_and_nonlean, VGBnames_list, by.x = "feature", by.y = "VGB", all.x = TRUE)

lean_and_nonlean$color <- ifelse(grepl("^Unclassified", lean_and_nonlean$fullname), "gray", "black")

ggplot(lean_and_nonlean, aes(x=coef.x, y=coef.y, color=color)) +
  xlab("Beta coefficient for lean MASLD vs controls") +
  ylab("Beta coefficient for nonlean MASLD vs controls") +
  xlim(-0.95, 0.95) +
  ylim(-0.6, 0.6) +
  geom_point() +
  scale_color_identity() +  
  geom_text_repel(data = subset(lean_and_nonlean, !grepl("^Unclassified", fullname)),
                  aes(label = gsub(" ", "_", gsub("VGB_", "VGB", fullname))),
                  fontface='bold',
                  size=5,
                  box.padding = 0.35,
                  point.padding = 0.5,
                  max.overlaps = 30) +
  geom_vline(xintercept=0, linetype = "dashed") +
  geom_hline(yintercept=0, linetype="dashed") +
  geom_label(aes(x = -0.8, y = 0.55, label = "Nonlean MASLD increased\nLean MASLD decreased"),
             label.padding = unit(6, "mm"), fill = "lightgrey", color="black", fontface='bold') +
  geom_label(aes(x = -0.8, y = -0.55, label = "Nonlean MASLD decreased\nLean MASLD decreased"),
             label.padding = unit(6, "mm"), fill = "lightgrey", color="black", fontface='bold') +
  geom_label(aes(x = 0.8, y = 0.55, label = "Nonlean MASLD increased\nLean MASLD increased"),
             label.padding = unit(6, "mm"), fill = "lightgrey", color="black", fontface='bold') +
  geom_label(aes(x = 0.8, y = -0.55, label = "Nonlean MASLD decreased\nLean MASLD increased"),
             label.padding = unit(6, "mm"), fill = "lightgrey", color="black", fontface='bold') +
  theme_bw() +
  theme(panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"))
