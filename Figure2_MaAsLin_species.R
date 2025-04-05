##################################################
#R program for creating Figure 2
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

setwd("~/b2b")
tax<-read_tsv("metaphlan_taxonomic_profiles.tsv")

#metadata
final_metadata <- read.delim('input/meta_df.tsv',row.names=1)
final_metadata <- final_metadata[!duplicated(final_metadata$barcode_metagenomics),]

#other metadata (medication)
other_metadata <- read_csv('input/b2b_medication.csv') %>% select(c(alias_id, contains("medication")))

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
  distinct(barcode_metagenomics, .keep_all = TRUE) %>% column_to_rownames("barcode_metagenomics")

new_data<-df_w_meta %>% rownames_to_column("barcode_metagenomics")
check_medication <- left_join(new_data,other_metadata,by="alias_id") %>% column_to_rownames("barcode_metagenomics")

species.data<-df_w_meta %>% select(matches("SGB|EUK"))

nafld.data<-df_w_meta %>% filter(cohort=="NAFLD")
species.nafld.data<-nafld.data %>% select(matches("SGB|EUK"))

oralspecies <- read_csv('oral_vs_gut.csv')

#alias id key
alias_key<-final_metadata %>% select(alias_id, barcode_metagenomics)

###########
###NAFLD###
###########
#211 nafld cases and 502 controls (193 matched and 309 unmatched controls)
nafld_data<-df_w_meta %>% filter(cohort=="NAFLD" | case==0) %>%
  mutate(lean_nafld_lean_control = case_when(bmi17v < 25 & case==1 ~ 1, bmi17v < 25 & case==0 ~ 0)) %>%
  mutate(nonlean_nafld_nonlean_control = case_when(bmi17v >= 25 & case==1 ~ 1, bmi17v >= 25 & case==0 ~ 0)) %>%
  mutate(lean_vs_nonlean_case = case_when(bmi17v >= 25 & case==1 ~ 1, bmi17v < 25 & case==1 ~ 0)) #nonlean case is 1, lean case is 0

#check medication
nafld_data_check<-check_medication %>% filter(cohort=="NAFLD" | case==0) %>%
  mutate(lean_vs_nonlean_case = case_when(bmi17v >= 25 & case==1 ~ 1, bmi17v < 25 & case==1 ~ 0)) #nonlean case is 1, lean case is 0
###let us get the lean nafld cases that do not fit into the cardiometabolic diagnosis criteria
lean_nafld <- nafld_data_check %>%
  select(!matches("SGB|EUK")) %>%
  filter(bmi17v < 25 & case==1 ) %>%
  select(c(hbp17, chol17, db17, contains("medication")))

lean_nafld_without_cardio <- lean_nafld %>%
  mutate(across(c(hbp17, chol17, db17), ~as.integer(. == "yes"))) %>%
  filter(rowSums(select(., hbp17, chol17, db17)) == 0)
lean_nafld_without_cardio_no_medication <- lean_nafld_without_cardio %>% filter(rowSums(select(., contains("medication")) == 1, na.rm = TRUE) == 0)
lean_nafld_without_cardio_samples<-rownames(lean_nafld_without_cardio_no_medication)

nafld_data_species <- nafld_data %>% select(matches("SGB|EUK")) 

# assign 1/2 of minimum to 0 so we can log transform
df.min <- (min(nafld_data_species[nafld_data_species > 0])/2)
log_nafld_data_species<-nafld_data_species
log_nafld_data_species[log_nafld_data_species == 0] <- df.min
log_nafld_data_species <- log10(log_nafld_data_species)

#MaAsLin
#adjusted for age, most recent diabetes, cumulative avg of physical activity, cum avg of bmi, cum avg of AHEI
fit_data <- Maaslin2(
  input_data = nafld_data %>% select(matches("SGB|EUK")), #species
  input_metadata = nafld_data %>% select(!matches("SGB|EUK")), #metadata
  output="output.mp4",
  normalization = "TSS", #our data already normalized
  transform = "LOG",
  analysis_method = "LM",
  max_significance = 0.20, # q-value threshold for significance. default is 0.25
  random_effects = NULL,
  fixed_effects = c('case', 'age', 'db17', 'act17v', 'bmi17v', 'aheiv2010_15'),
  correction = "BH",
  standardize = TRUE,
  cores = 1,
  plot_heatmap = TRUE,
  plot_scatter = TRUE,
  heatmap_first_n = 50)

###Figure 2A
#Volcano plot
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
top_bugs <- bind_rows(
  results.nafld %>%
    filter(color == 'Up') %>%
    arrange(qval, desc(abs(coef))) %>%
    head(top),
  results.nafld %>%
    filter(color == 'Down') %>%
    arrange(qval, desc(abs(coef))) %>%
    head(top)
)

options(ggrepel.max.overlaps = Inf)

ggplot(results.nafld, aes(x = coef,
                             y = -log(qval, 10),
                             colour= color)) + # -log10 conversion
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
                                 label = feature),
                   show.legend = FALSE,
                   size = 5.5,
                   min.segment.length = 0.2,
                   color = "black",
                   fontface='bold')

#oral species
sigresults <- read_tsv('output.mp4/significant_results.tsv') %>% filter(metadata=="case")
allresults <- read_tsv('output.mp4/all_results.tsv')
number_of_features_after_prev_filter<-unique(allresults$feature)

nafld_data_species_and_case <- nafld_data %>% select(c(matches("SGB|EUK"),case))
casecontrol<-nafld_data_species_and_case %>% select(case)
oralspecies_selected<-oralspecies %>% filter(class=="oral")
oral_dataframe<-nafld_data_species_and_case[, colnames(nafld_data_species_and_case) %in% oralspecies_selected$speciesname]
sigresults_oral<-sigresults %>% filter(sigresults$feature %in% oralspecies_selected$speciesname)

sumabund<-rowSums(oral_dataframe) %>% as.data.frame() %>% rename("sumabund" = ".")
sumabund_and_case <- merge(sumabund,casecontrol,by = 'row.names')

###Figure 2B
oral<-ggplot(data = sumabund_and_case,
       aes(x = as.factor(case), y = log10(sumabund), fill = as.factor(case))) +
  scale_fill_manual(values = c("0" = "#999999", "1" = "#E69F00")) +
  geom_boxplot(notch = FALSE) +
  geom_jitter(shape = 16, position = position_jitter(0.1)) +
  theme_classic(base_size = 18) +
  theme(legend.position = "none") +
  ylab("log10(abundance)") +
  geom_signif(comparisons = list(c("0", "1")),
              test=wilcox.test,
              tip_length = 0,
              map_signif_level = TRUE) +
  ggtitle("Oral SGBs")

###Figure 2C
#genus boxplots using metaphlan taxonomic profiles
#Read in taxonomic data
tax<-read_tsv("metaphlan_taxonomic_profiles.tsv")
names(tax) = gsub("_taxonomic_profile", "", names(tax))
strepveilloactino.data<-tax %>% filter(`# taxonomy`=="k__Bacteria|p__Firmicutes|c__Bacilli|o__Lactobacillales|f__Streptococcaceae|g__Streptococcus" | `# taxonomy`=="k__Bacteria|p__Firmicutes|c__Negativicutes|o__Veillonellales|f__Veillonellaceae|g__Veillonella" | `# taxonomy` == "k__Bacteria|p__Actinobacteria|c__Actinobacteria|o__Actinomycetales|f__Actinomycetaceae|g__Actinomyces") %>% select(-`# taxonomy`) %>% t() %>% as.data.frame() %>% rownames_to_column() %>% rename(strep.abund=V1,veillo.abund=V2,actino.abund=V3,barcode_metagenomics=rowname)
nafld_data_id_leannonlean <- nafld_data %>% select(lean_vs_nonlean_case) %>% rownames_to_column("barcode_metagenomics")
nafld_data_species_w_strepveilloactino <- left_join(nafld_data_id_leannonlean,strepveilloactino.data,by="barcode_metagenomics") %>% na.omit()
strep<-ggplot(data = nafld_data_species_w_strepveilloactino,
       aes(x = as.factor(lean_vs_nonlean_case), y = log10(strep.abund/100), fill = as.factor(lean_vs_nonlean_case))) +
  scale_fill_manual(values = c("0" = "blue", "1" = "red")) +
  geom_boxplot(notch = FALSE) +
  geom_jitter(shape = 16, position = position_jitter(0.1)) +
  theme_classic(base_size = 18) +
  theme(legend.position = "none") +
  ylab("log10(abundance)") +
  geom_signif(comparisons = list(c("0", "1")),
              tip_length = 0,
              map_signif_level = TRUE) +
  ggtitle("Streptococcus spp.")

###Figure 2D
#lean nafld vs lean control MaAsLin
fit_data_lean <- Maaslin2(
  input_data = nafld_data %>% select(matches("SGB|EUK")), #species
  input_metadata = nafld_data %>% select(!matches("SGB|EUK")), #metadata
  output="output.mp4/lean",
  normalization = "TSS",
  transform = "LOG",
  analysis_method = "LM",
  max_significance = 0.2, # q-value threshold for significance. default is 0.25
  random_effects = NULL,
  fixed_effects = c('lean_nafld_lean_control', 'age', 'db17', 'act17v', 'aheiv2010_15'),
  correction = "BH",
  standardize = TRUE,
  cores = 1,
  plot_heatmap = TRUE,
  plot_scatter = TRUE,
  heatmap_first_n = 50)

results <- fit_data_lean[["results"]]

leanresults <- read_tsv('output.mp4/lean/all_results.tsv') %>% filter(metadata=="lean_nafld_lean_control")

#nonlean nafld vs nonlean control MaAsLin
fit_data_nonlean <- Maaslin2(
  input_data = nafld_data %>% select(matches("SGB|EUK")), #species
  input_metadata = nafld_data %>% select(!matches("SGB|EUK")), #metadata
  output="output.mp4/nonlean",
  normalization = "TSS",
  transform = "LOG",
  analysis_method = "LM",
  max_significance = 0.20, # q-value threshold for significance. default is 0.25
  random_effects = NULL,
  fixed_effects = c('nonlean_nafld_nonlean_control', 'age', 'db17', 'act17v', 'aheiv2010_15'),
  correction = "BH",
  standardize = TRUE,
  cores = 1,
  plot_heatmap = TRUE,
  plot_scatter = TRUE,
  heatmap_first_n = 50)

results <- fit_data_nonlean[["results"]]

nonleanresults <- read_tsv('output.mp4/nonlean/all_results.tsv') %>% filter(metadata=="nonlean_nafld_nonlean_control")

#scatter plot for lean and nonlean
lean_and_nonlean<-full_join(leanresults,nonleanresults,by="feature")
oral_lean_and_nonlean<-lean_and_nonlean %>% filter(feature %in% oralspecies_selected$speciesname)
#scatter plot for oral SGBs for lean and nonlean
model<-lm(coef.y~coef.x,data=oral_lean_and_nonlean)
#with v=0 and h=0 lines
ggplot(oral_lean_and_nonlean, aes(x=coef.x, y=coef.y)) +
  xlab("Beta coefficient for lean MASLD vs controls")+
  ylab("Beta coefficient for nonlean MASLD vs controls") +
  xlim(-0.5,0.75) +
  ylim(-0.2,0.65) +
  geom_point()+
  geom_text_repel(aes(label=feature),fontface='bold',size=4)+
  geom_vline(xintercept=0,linetype = "dashed")+
  geom_hline(yintercept=0,linetype="dashed")+
  geom_label(aes(x = -0.35, y = 0.65, label = "Nonlean MASLD increased\nLean MASLD decreased"),
               label.padding = unit(4, "mm"),  fill = "lightgrey", color="black",fontface='bold')+
  geom_label(aes(x = -0.35, y = -0.2, label = "Nonlean MASLD decreased\nLean MASLD decreased"),
             label.padding = unit(4, "mm"),  fill = "lightgrey", color="black",fontface='bold')+
  geom_label(aes(x = 0.6, y = 0.65, label = "Nonlean MASLD increased\nLean MASLD increased"),
             label.padding = unit(4, "mm"),  fill = "lightgrey", color="black",fontface='bold')+
  geom_label(aes(x = 0.6, y = -0.2, label = "Nonlean MASLD decreased\nLean MASLD increased"),
             label.padding = unit(4, "mm"),  fill = "lightgrey", color="black",fontface='bold')+
  theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))

###Supp table 13
#MaAsLin without the cases that did not meet the new MASLD diagnosis criteria
nafld_data_sens <- nafld_data[!rownames(nafld_data) %in% lean_nafld_without_cardio_samples, ]
fit_data <- Maaslin2(
  input_data = nafld_data_sens %>% select(matches("SGB|EUK")), #species
  input_metadata = nafld_data_sens %>% select(!matches("SGB|EUK")), #metadata
  output="output.mp4.without.diagnosis",
  normalization = "TSS", #our data already normalized
  transform = "LOG",
  analysis_method = "LM",
  max_significance = 0.20, # q-value threshold for significance. default is 0.25
  random_effects = NULL,
  fixed_effects = c('case', 'age', 'db17', 'act17v', 'bmi17v', 'aheiv2010_15'),
  correction = "BH",
  standardize = TRUE,
  cores = 1,
  plot_heatmap = TRUE,
  plot_scatter = TRUE,
  heatmap_first_n = 50)
