##################################################
#R program for creating Figure 3
##################################################

library(dplyr)
library(tidyverse)
library(stringr)
library(readr)
library(Maaslin2)
library(ggplot2)
library(ggrepel)
library(stats)
library(gridExtra)

setwd("~/b2b")
unfiltered_species <- read.delim('input/species_unfiltered.tsv',row.names=1) %>% t() %>% as.data.frame() %>% rownames_to_column() %>% rename(barcode_metagenomics = rowname) 
final_metadata <- read.delim('input/meta_df.tsv',row.names=1)
final_metadata <- final_metadata[!duplicated(final_metadata$barcode_metagenomics),]
df_w_meta <- left_join(unfiltered_species,final_metadata,by="barcode_metagenomics")

species.data<-df_w_meta %>% select(starts_with('s__'))

nafld.data<-df_w_meta %>% filter(cohort=="NAFLD") 
species.nafld.data<-nafld.data %>% select(starts_with('s__')) 
species.nafld.data<-species.nafld.data/100

#get oral species
oralspecies <- read_csv('input/oralvsgut.csv') 
oralspecies_selected<-oralspecies %>% filter(major_site=="oral")

#get the filtered (10%) list of species
species_maaslin_allresults <- read_tsv("output/all_results.tsv")
filtered_species_result<-species_maaslin_allresults$feature
filtered_species_list<-gsub("s__","",filtered_species_result)

#pathway code starts here
pathway_file <- read.delim('pathabundance_unstrat.tsv',row.names=1) %>% t() %>% as.data.frame() %>% rownames_to_column() %>% rename(barcode_metagenomics = rowname) 
pathway_file$barcode_metagenomics = gsub("_Abundance", "", pathway_file$barcode_metagenomics)
final_metadata <- read.delim('input/meta_df.tsv',row.names=1)
final_metadata <- final_metadata[!duplicated(final_metadata$barcode_metagenomics),]
df_w_meta <- left_join(pathway_file,final_metadata,by="barcode_metagenomics") 

pathway_list<-read.delim('pathabundance_unstrat.tsv',row.names=1) %>% rownames_to_column()
path.data<-df_w_meta %>% select(pathway_list$rowname)

nafld.data<-df_w_meta %>% filter(cohort=="NAFLD") 
path.nafld.data<-nafld.data %>% select(pathway_list$rowname) 

#stratified
stratified_pathway_file <- read.delim('input/pathabundance_NAFLD.tsv',row.names=1) %>% t() %>% as.data.frame() %>% rownames_to_column() %>% rename(barcode_metagenomics = rowname) 
names_pathway<-stratified_pathway_file %>% select(-c(barcode_metagenomics,NAFLD))
colnames_pathway<-colnames(names_pathway)
df_w_meta_stratified <- left_join(stratified_pathway_file,final_metadata,by="barcode_metagenomics") %>%
  mutate(lean_nafld = case_when(bmi17v >= 25 & case==1 ~ 'Nonlean NAFLD', bmi17v <25 & case==1 ~ 'Lean NAFLD', bmi17v >= 25 & case==0 ~ 'Control', bmi17v <25 & case==0 ~ 'Control')) %>%
  select(c(colnames_pathway),case,lean_nafld,barcode_metagenomics) %>% column_to_rownames("barcode_metagenomics") %>% t()

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

###Figure 3A
#boxplot
require(reshape2)
require(graphics)
selected_pathways <- c("GLUTORN-PWY: L-ornithine biosynthesis I",
                       "PWY-2941: L-lysine biosynthesis II",
                       "P4-PWY: superpathway of L-lysine, L-threonine and L-methionine biosynthesis I",
                       "ARGSYN-PWY: L-arginine biosynthesis I (via L-ornithine)",
                       "ARGSYNBSUB-PWY: L-arginine biosynthesis II (acetyl cycle)",
                       "PWY-7977: L-methionine biosynthesis IV",
                       "PWY-702: L-methionine biosynthesis II")

nafld_data_path_w_case <- nafld_data %>% select(all_of(selected_pathways) | case | alias_id) 
sig.path.m <- melt(nafld_data_path_w_case, id = c("alias_id","case"))

ggplot(data = sig.path.m, aes(x = variable, y = log10(value))) +
  geom_boxplot(aes(fill = as.factor(case)), outlier.size = 0.3) +
  coord_flip() +
  scale_fill_manual(name = "MASLD", values = c("#999999", "#E69F00"),
                    labels = c("control", "case")) +
  theme_bw() +
  theme(axis.text.y = element_text(size = 12),
        legend.position="bottom") +
  ylab("log10(relative abundance)") +
  xlab("DNA")  #10*5

#plot the oral ones
bug_colors <- c(
  # Oral bugs (green)
  "Strep"      = "blue3",
  "Non-Strep"  = "lightblue",
  "Non-oral"   = "orange")
plot_stratified_distribution <- function(name) {

  pcl=t(df_w_meta_stratified)
  
  pwy <- grepl(name, colnames(pcl), fixed=T)
  sum_stratified <- pcl[,pwy,drop=F] 
  t_sum_stratified <- as.data.frame(sum_stratified)
  
  colnames(t_sum_stratified) <- sub(".*\\|", "", colnames(t_sum_stratified))
  colnames(t_sum_stratified) <- gsub(".*s__","", colnames(t_sum_stratified))

  #transform to relative abundance by sum normalization and filtering
  t_sum_stratified[ is.na(t_sum_stratified) ] <- NA
  # make numeric
  for(i in 1:ncol(t_sum_stratified)) 
  { 
    t_sum_stratified[ ,i] <- as.numeric(as.character(t_sum_stratified[,i])) 
  }
  dim(t_sum_stratified) 

  # sum normalize - transform to relative abundance and filtering
  t_sum_stratified_sweep <- sweep(t_sum_stratified, 1, rowSums(t_sum_stratified), `/`) #mean of that species (divide by sum of all the abundance for each sample)
  t_sum_stratified_sweep[is.na(t_sum_stratified_sweep)] <- 0
  stratified <- as.data.frame(t(t_sum_stratified_sweep)) 

  #sort by mean relative abundance
  mns <- rowMeans(stratified, na.rm=TRUE)
  order(-mns)
  stratified <- stratified[order(-mns),]
  rowMeans(stratified)

  #oral species vs non-oral -- strep vs non-strep
  oral_strep<-oralspecies_selected[grepl("Strep", oralspecies_selected$feature), ]
  oral_nonstrep<-oralspecies_selected[!grepl("Strep", oralspecies_selected$feature), ]
  strep <- as.data.frame(t(colSums(stratified[rownames(stratified) %in% c(oral_strep$feature),])))
  rownames(strep) <- "Strep"
  nonstrep <- as.data.frame(t(colSums(stratified[rownames(stratified) %in% c(oral_nonstrep$feature),])))
  rownames(nonstrep) <- "Non-Strep"
  stratified_new <- rbind(strep,nonstrep)
  
  #take difference for cases and controls
  #include case status
  df_case <- pcl %>% as.data.frame %>% select("case") %>% t() %>% as.data.frame
  stratified_case <- rbind(stratified_new, df_case)  # Combine df1 and df_case using rbind()
  
  # Identify the case/control status
  case_control <- stratified_case[nrow(stratified_case), ]
  # Calculate the mean abundance for case and control
  mean_case <- rowMeans(stratified_new[, case_control == 1])
  mean_control <- rowMeans(stratified_new[, case_control == 0])
  # Calculate the difference in mean abundance between case and control
  diff_mean <- (mean_case - mean_control)*100
  diff.abundance<-diff_mean %>% as.data.frame() %>% rename('diff'='.') %>% rownames_to_column() %>% mutate(ec = name)
  
  dna_pwy_bar<-ggplot(diff.abundance, aes(x = ec, y = diff, fill = rowname)) +
    geom_bar(stat = "identity") +
    theme_void() +  # Use theme_void() to remove the outline box and tick marks
    scale_y_continuous(breaks = c(-0.5, 0, 0.5), limits = c(-2, 2)) +
    geom_hline(yintercept = 0, linetype = "dashed", size = 0.5) +
    scale_fill_manual(values = bug_colors) +
    coord_flip() +
    theme(legend.position = "none")
  print(dna_pwy_bar)
  
}
plot7<-plot_stratified_distribution("GLUTORN-PWY: L-ornithine biosynthesis I")
plot6<-plot_stratified_distribution("PWY-2941: L-lysine biosynthesis II")
plot5<-plot_stratified_distribution("P4-PWY: superpathway of L-lysine, L-threonine and L-methionine biosynthesis I")
plot4<-plot_stratified_distribution("ARGSYN-PWY: L-arginine biosynthesis I (via L-ornithine)")
plot3<-plot_stratified_distribution("ARGSYNBSUB-PWY: L-arginine biosynthesis II (acetyl cycle)")
plot2<-plot_stratified_distribution("PWY-7977: L-methionine biosynthesis IV")
plot1<-plot_stratified_distribution("PWY-702: L-methionine biosynthesis II")

grid.arrange(plot1,plot2,plot3,plot4,plot5,plot6,plot7,
             ncol=1,nrow=7) #10*6.5

###Figure 3B
###### Nonlean vs lean case for MGX (DNA)
#boxplot
require(reshape2)
require(graphics)

nafld_data_new <- nafld_data
nafld_data_new$lean_vs_nonlean_case<-as.factor(nafld_data_new$lean_vs_nonlean_case)
nafld_data_new$lean_vs_nonlean_case<-as.factor(nafld_data_new$lean_vs_nonlean_case)
nafld_data_new <- nafld_data_new %>% filter(!is.na(lean_vs_nonlean_case))

fit_data <- Maaslin2(
  input_data = nafld_data_new %>% select(pathway_list$rowname), #nafld_data_pathway, 
  input_metadata =nafld_data_new %>% select(!pathway_list$rowname), #metadata
  output="output_pathway_notfiltered_nonleanlean",
  normalization = "TSS", 
  transform = "LOG",
  analysis_method = "LM",
  max_significance = 0.20, # q-value threshold for significance. default is 0.25
  random_effects = NULL,
  fixed_effects = c('lean_vs_nonlean_case', 'age', 'db17', 'act17v', 'aheiv2010_15'),
  correction = "BH",
  standardize = TRUE,
  cores = 30,
  plot_heatmap = TRUE,
  plot_scatter = TRUE,
  heatmap_first_n = 50)

#get significant pathways for nonlean vs lean case
dataframe_sig <-read.table(file = '~/b2b/output_pathway_notfiltered_nonleanlean/significant_results.tsv', sep = '\t', header = TRUE) %>% filter(metadata=="lean_vs_nonlean_case") %>% select(c(feature, coef, qval))
sig_pathways <- c("COLANSYN-PWY: colanic acid building blocks biosynthesis",
                  "PWY-6588: pyruvate fermentation to acetone",
                  "PWY-7323: superpathway of GDP-mannose-derived O-antigen building blocks biosynthesis",
                  "P4-PWY: superpathway of L-lysine, L-threonine and L-methionine biosynthesis I",
                  "PWY-4984: urea cycle",
                  "PWY0-781: aspartate superpathway",
                  "PWY-7560: methylerythritol phosphate pathway II",
                  "PWY-7197: pyrimidine deoxyribonucleotide phosphorylation",
                  "PWY-1861: formaldehyde assimilation II (assimilatory RuMP Cycle)",
                  "P185-PWY: formaldehyde assimilation III (dihydroxyacetone cycle)",
                  "PWY-6270: isoprene biosynthesis I",
                  "PWY0-1479: tRNA processing",
                  "PWY-6895: superpathway of thiamine diphosphate biosynthesis II",
                  "PWY-6749: CMP-legionaminate biosynthesis I")
#for nafld_data_new, 1 is nonlean case and 0 is lean case
nafld_data_path_w_case <- nafld_data_new %>% select(all_of(sig_pathways) | lean_vs_nonlean_case | alias_id) 
sig.path.m <- melt(nafld_data_path_w_case, id = c("alias_id","lean_vs_nonlean_case"))

df_w_meta_stratified_nonleanlean <- left_join(stratified_pathway_file,final_metadata,by="barcode_metagenomics") %>%
  mutate(lean_nafld = case_when(bmi17v >= 25 & case==1 ~ 'Nonlean NAFLD', bmi17v <25 & case==1 ~ 'Lean NAFLD', bmi17v >= 25 & case==0 ~ 'Control', bmi17v <25 & case==0 ~ 'Control')) %>%
  mutate(lean_vs_nonlean_case = case_when(bmi17v >= 25 & case==1 ~ 1, bmi17v < 25 & case==1 ~ 0)) %>% #nonlean case is 1, lean case is 0
  filter(!is.na(lean_vs_nonlean_case)) %>%
  select(c(colnames_pathway),case,lean_vs_nonlean_case,barcode_metagenomics) %>% column_to_rownames("barcode_metagenomics") %>% t()

ggplot(data = sig.path.m, aes(x = variable, y = log10(value))) +
  geom_boxplot(aes(fill = as.factor(lean_vs_nonlean_case)), outlier.size = 0.3) +
  coord_flip() +
  scale_fill_manual(name = "MASLD", values = c("blue", "red"),
                    labels = c("Lean case", "Nonlean case")) +
  theme_bw() +
  theme(axis.text.y = element_text(size = 12),
        legend.position="bottom") +
  ylab("log10(relative abundance)") +
  xlab("DNA")  #10*5

#plot the oral ones
bug_colors <- c(
  # Oral bugs (green)
  "Strep"      = "blue3",
  "Non-Strep"  = "lightblue",
  "Non-oral"   = "orange")
plot_stratified_distribution <- function(name) {
  pcl=t(df_w_meta_stratified_nonleanlean)
  
  pwy <- grepl(name, colnames(pcl), fixed=T)
  sum_stratified <- pcl[,pwy,drop=F] 
  t_sum_stratified <- as.data.frame(sum_stratified)
  
  colnames(t_sum_stratified) <- sub(".*\\|", "", colnames(t_sum_stratified))
  colnames(t_sum_stratified) <- gsub(".*s__","", colnames(t_sum_stratified))
  
  #transform to relative abundance by sum normalization and filtering
  t_sum_stratified[ is.na(t_sum_stratified) ] <- NA
  # make numeric
  for(i in 1:ncol(t_sum_stratified)) 
  { 
    t_sum_stratified[ ,i] <- as.numeric(as.character(t_sum_stratified[,i])) 
  }
  dim(t_sum_stratified) 

  # sum normalize - transform to relative abundance and filtering
  t_sum_stratified_sweep <- sweep(t_sum_stratified, 1, rowSums(t_sum_stratified), `/`) #mean of that species (divide by sum of all the abundance for each sample)
  t_sum_stratified_sweep[is.na(t_sum_stratified_sweep)] <- 0
  stratified <- as.data.frame(t(t_sum_stratified_sweep)) 

  #sort by mean relative abundance
  mns <- rowMeans(stratified, na.rm=TRUE)
  order(-mns)
  stratified <- stratified[order(-mns),]
  rowMeans(stratified)
  
  #oral species vs non-oral -- strep vs non-strep
  oral_strep<-oralspecies_selected[grepl("Strep", oralspecies_selected$feature), ]
  oral_nonstrep<-oralspecies_selected[!grepl("Strep", oralspecies_selected$feature), ]
  strep <- as.data.frame(t(colSums(stratified[rownames(stratified) %in% c(oral_strep$feature),])))
  rownames(strep) <- "Strep"
  nonstrep <- as.data.frame(t(colSums(stratified[rownames(stratified) %in% c(oral_nonstrep$feature),])))
  rownames(nonstrep) <- "Non-Strep"
  stratified_new <- rbind(strep,nonstrep)
  
  #take difference for cases and controls
  #include case status
  df_case <- pcl %>% as.data.frame %>% select("lean_vs_nonlean_case") %>% t() %>% as.data.frame
  stratified_case <- rbind(stratified_new, df_case)  # Combine df1 and df_case using rbind()
  
  # Identify the case/control status
  case_control <- stratified_case[nrow(stratified_case), ]
  # Calculate the mean abundance for case and control
  mean_case <- rowMeans(stratified_new[, case_control == 1])
  mean_control <- rowMeans(stratified_new[, case_control == 0])
  # Calculate the difference in mean abundance between case and control
  diff_mean <- (mean_case - mean_control)*100
  diff.abundance<-diff_mean %>% as.data.frame() %>% rename('diff'='.') %>% rownames_to_column() %>% mutate(ec = name)
  
  dna_pwy_bar<-ggplot(diff.abundance, aes(x = ec, y = diff, fill = rowname)) +
    geom_bar(stat = "identity") +
    theme_void() +  # Use theme_void() to remove the outline box and tick marks
    scale_y_continuous(breaks = c(-0.5, 0, 0.5), limits = c(-2, 2)) +
    geom_hline(yintercept = 0, linetype = "dashed", size = 0.5) +
    scale_fill_manual(values = bug_colors) +
    coord_flip() +
    theme(legend.position = "none")
  print(dna_pwy_bar)
  
}
plot14<-plot_stratified_distribution("COLANSYN-PWY: colanic acid building blocks biosynthesis")
plot13<-plot_stratified_distribution("PWY-6588: pyruvate fermentation to acetone")
plot12<-plot_stratified_distribution("PWY-7323: superpathway of GDP-mannose-derived O-antigen building blocks biosynthesis")
plot11<-plot_stratified_distribution("P4-PWY: superpathway of L-lysine, L-threonine and L-methionine biosynthesis I")
plot10<-plot_stratified_distribution("PWY-4984: urea cycle")
plot9<-plot_stratified_distribution("PWY0-781: aspartate superpathway")
plot8<-plot_stratified_distribution("PWY-7560: methylerythritol phosphate pathway II")
plot7<-plot_stratified_distribution("PWY-7197: pyrimidine deoxyribonucleotide phosphorylation")
plot6<-plot_stratified_distribution("PWY-1861: formaldehyde assimilation II (assimilatory RuMP Cycle)")
plot5<-plot_stratified_distribution("P185-PWY: formaldehyde assimilation III (dihydroxyacetone cycle)")
plot4<-plot_stratified_distribution("PWY-6270: isoprene biosynthesis I")
plot3<-plot_stratified_distribution("PWY0-1479: tRNA processing")
plot2<-plot_stratified_distribution("PWY-6895: superpathway of thiamine diphosphate biosynthesis II")
plot1<-plot_stratified_distribution("PWY-6749: CMP-legionaminate biosynthesis I")

grid.arrange(plot1,plot2,plot3,plot4,plot5,plot6,plot7,
             plot8,plot9,plot10,plot11,plot12,plot13,plot14,
             ncol=1,nrow=14) #10*6.5

###Figure 3A MTX
############
##### Do the same for MTX (RNA)
############
#####pathway
mtxpathway<-read_tsv("pathabundance_relab_nospecial.tsv")
names(mtxpathway) = gsub("_Abundance", "", names(mtxpathway))
mtx_pathway<-mtxpathway %>% column_to_rownames("# Pathway") %>% t() %>% as.data.frame() 
samples <- mtx_pathway %>% rownames_to_column()
samplenames = samples$rowname 
mtx_pathway_unstratified<-mtx_pathway[,!grepl("\\|", colnames(mtx_pathway)) ] 

final_metadata <- read.delim('input/meta_df.tsv',row.names=1)
final_metadata <- final_metadata[!duplicated(final_metadata$barcode_metagenomics),]
mtx_pathway_for_join<-mtx_pathway_unstratified %>% rownames_to_column("barcode_metagenomics")
pathways_list <- colnames(mtx_pathway_unstratified)
df_w_meta <- left_join(mtx_pathway_for_join,final_metadata,by="barcode_metagenomics") 
pathway.nafld.data<-df_w_meta %>% filter(cohort=="NAFLD" | case==0) %>% column_to_rownames("barcode_metagenomics")
pathway.nafld.data<-pathway.nafld.data[order(row.names(pathway.nafld.data)), ]

mtx_input<-pathway.nafld.data%>%select(all_of(pathways_list))
metadata_input<-pathway.nafld.data%>%select(!all_of(pathways_list))
num_rows <- sum(rowSums(mtx_input == 0) == ncol(mtx_input))
final_metadata <- metadata_input %>% select(c(case,age,db17,act17v,bmi17v,aheiv2010_15,cohort))
new_final_metadata <- final_metadata %>% rownames_to_column("barcode_metagenomics")
#stratified
stratified_pathway_file<-read_tsv("pathabundance_relab_nospecial.tsv")
names(stratified_pathway_file) = gsub("_Abundance", "", names(stratified_pathway_file))
stratified_pathway_file<-stratified_pathway_file %>% column_to_rownames("# Pathway") %>% t() %>% as.data.frame() 
colnames_pathway<-colnames(stratified_pathway_file)
new_stratified_pathway_file<-stratified_pathway_file %>% rownames_to_column("barcode_metagenomics")

df_w_meta_stratified <- left_join(new_stratified_pathway_file,new_final_metadata,by="barcode_metagenomics") %>%
  filter(cohort=="NAFLD" | case==0) %>%
  mutate(lean_nafld = case_when(bmi17v >= 25 & case==1 ~ 'Nonlean NAFLD', bmi17v <25 & case==1 ~ 'Lean NAFLD', bmi17v >= 25 & case==0 ~ 'Control', bmi17v <25 & case==0 ~ 'Control')) %>%
  select(c(colnames_pathway),case,lean_nafld,barcode_metagenomics) %>% column_to_rownames("barcode_metagenomics") %>% t()

#get oral species
oralspecies <- read_csv('input/oralvsgut.csv') 
oralspecies_selected<-oralspecies %>% filter(major_site=="oral")
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

#boxplot
require(reshape2)
require(graphics)
selected_pathways <- c("GLUTORN-PWY: L-ornithine biosynthesis I",
                       "PWY-2941: L-lysine biosynthesis II",
                       "PWY-724: superpathway of L-lysine, L-threonine and L-methionine biosynthesis II",
                       "ARGSYN-PWY: L-arginine biosynthesis I (via L-ornithine)",
                       "ARGSYNBSUB-PWY: L-arginine biosynthesis II (acetyl cycle)",
                       "PWY-7977: L-methionine biosynthesis IV",
                       "PWY-702: L-methionine biosynthesis II")

nafld_data_path_w_case <- nafld_data %>% select(all_of(selected_pathways) | case | alias_id) 
sig.path.m <- melt(nafld_data_path_w_case, id = c("alias_id","case"))

ggplot(data = sig.path.m, aes(x = variable, y = log10(value))) +
  geom_boxplot(aes(fill = as.factor(case)), outlier.size = 0.3) +
  coord_flip() +
  scale_fill_manual(name = "MASLD", values = c("#999999", "#E69F00"),
                    labels = c("control", "case")) +
  theme_bw() +
  theme(axis.text.y = element_text(size = 12),
        legend.position="bottom") +
  ylab("log10(relative abundance)") +
  xlab("RNA")  #10*5

###Figure 3B MTX
###### Nonlean vs lean case for MTX (RNA)
#boxplot
require(reshape2)
require(graphics)

nafld_data_new <- nafld_data
nafld_data_new$lean_vs_nonlean_case<-as.factor(nafld_data_new$lean_vs_nonlean_case)
nafld_data_new$lean_vs_nonlean_case<-as.factor(nafld_data_new$lean_vs_nonlean_case)
nafld_data_new <- nafld_data_new %>% 
  filter(!is.na(lean_vs_nonlean_case))

sig_pathways <- c("COLANSYN-PWY: colanic acid building blocks biosynthesis",
                  "PWY-6588: pyruvate fermentation to acetone",
                  "PWY-7323: superpathway of GDP-mannose-derived O-antigen building blocks biosynthesis",
                  "P4-PWY: superpathway of L-lysine, L-threonine and L-methionine biosynthesis I",
                  "PWY-4984: urea cycle",
                  "PWY0-781: aspartate superpathway",
                  "PWY-7560: methylerythritol phosphate pathway II",
                  "PWY-7197: pyrimidine deoxyribonucleotide phosphorylation",
                  "PWY-1861: formaldehyde assimilation II (assimilatory RuMP Cycle)",
                  "P185-PWY: formaldehyde assimilation III (dihydroxyacetone cycle)",
                  "PWY-6270: isoprene biosynthesis I",
                  "PWY0-1479: tRNA processing",
                  "PWY-6895: superpathway of thiamine diphosphate biosynthesis II",
                  "PWY-6749: CMP-legionaminate biosynthesis I")

# Add missing pathways as NA
for (pathway in sig_pathways) {
  if (!(pathway %in% colnames(nafld_data_new))) {
    nafld_data_new[[pathway]] <- NA
  }
}

#for nafld_data_new, 1 is nonlean case and 0 is lean case
nafld_data_path_w_case <- nafld_data_new %>% 
  select(all_of(sig_pathways), lean_vs_nonlean_case, alias_id) 

sig.path.m <- melt(nafld_data_path_w_case, id = c("alias_id","lean_vs_nonlean_case"))

df_w_meta_stratified_nonleanlean <- left_join(new_stratified_pathway_file,new_final_metadata,by="barcode_metagenomics") %>%
  mutate(lean_nafld = case_when(bmi17v >= 25 & case==1 ~ 'Nonlean NAFLD', bmi17v <25 & case==1 ~ 'Lean NAFLD', bmi17v >= 25 & case==0 ~ 'Control', bmi17v <25 & case==0 ~ 'Control')) %>%
  mutate(lean_vs_nonlean_case = case_when(bmi17v >= 25 & case==1 ~ 1, bmi17v < 25 & case==1 ~ 0)) %>% #nonlean case is 1, lean case is 0
  filter(!is.na(lean_vs_nonlean_case)) %>%
  select(c(colnames_pathway),case,lean_vs_nonlean_case,barcode_metagenomics) %>% column_to_rownames("barcode_metagenomics") %>% t()

ggplot(data = sig.path.m, aes(x = variable, y = log10(value))) +
  geom_boxplot(aes(fill = as.factor(lean_vs_nonlean_case)), outlier.size = 0.3) +
  coord_flip() +
  scale_fill_manual(name = "MASLD", values = c("blue", "red"),
                    labels = c("Lean case", "Nonlean case")) +
  theme_bw() +
  theme(axis.text.y = element_text(size = 12),
        legend.position="bottom") +
  ylab("log10(relative abundance)") +
  xlab("RNA")  #10*5
