##################################################
#R program for creating Figure 4
##################################################

library(dplyr)
library(tidyverse)
library(stringr)
library(readr)
library(Maaslin2)
library(ggplot2)
library(ggrepel)

setwd("~/b2b")
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

###########
###NAFLD###
###########
#209 NAFLD cases
nafld_data<-df_w_meta %>% filter(cohort=="NAFLD" | case==0) %>% 
  mutate(lean_vs_nonlean_case = case_when(bmi17v >= 25 & case==1 ~ 1, bmi17v < 25 & case==1 ~ 0)) #nonlean case is 1, lean case is 0

#alias id key
alias_key<-final_metadata %>% select(alias_id, barcode_metabolomics)

nafld_data_mbx <- nafld_data %>% select(all_of(mbx_list)) 

#assign 1/2 of minimum to 0 so we can log transform 
df.min <- (min(nafld_data_mbx[nafld_data_mbx > 0])/2)
log_nafld_data_mbx<-nafld_data_mbx
log_nafld_data_mbx[log_nafld_data_mbx == 0] <- df.min 
log_nafld_data_mbx <- log2(log_nafld_data_mbx)

fit_data <- Maaslin2(
  input_data = nafld_data_mbx, #mbx
  input_metadata =nafld_data %>% select(!all_of(mbx_list)), #metadata
  output="output_mbx/nonlog",
  normalization = "NONE", #the data is already median normalized data 
  transform = "LOG",
  analysis_method = "LM",
  max_significance = 0.20, # q-value threshold for significance. default is 0.25
  random_effects = NULL,
  fixed_effects = c('case', 'age', 'db17', 'bmi17v', 'act17v', 'aheiv2010_15'),
  min_prevalence = 0,
  min_abundance = 0,
  correction = "BH",
  standardize = TRUE,
  cores = 1,
  plot_heatmap = TRUE,
  plot_scatter = TRUE,
  heatmap_first_n = 50)

results <- fit_data[["results"]] 

results.nafld <- results %>%
  filter(metadata == "case") %>%
  mutate(adjlog2 = exp(coef)) %>%
  mutate(color = case_when(
    # changed this to 0 from 2
    coef > 0.25 & qval <= 0.20 ~ "Up", 
    coef < -0.25 & qval <=0.20 ~ "Down",
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

###Figure 4A
ggplot(results.nafld, aes(x = coef, 
                          y = -log(qval, 10),
                          colour= color)) + # -log10 conversion
  geom_point(alpha=0.4, size=3.5) +
  xlab(expression("B-coefficient")) +
  ylab(expression("-log"[10] * "(FDR p-value)")) + 
  xlim(c(-1.5, 1.5)) +
  scale_color_manual(values=c("blue", "grey","red")) +
  geom_vline(xintercept=c(-0.25, 0.25),lty=4, col="black",lwd=0.8) +
  geom_hline(yintercept = -log10(0.20),lty=4, col="black",lwd=0.8) +
  theme_classic(base_size = 18)+
  theme(plot.title = element_text(hjust = 0.5), 
        legend.position="right", 
        legend.title = element_blank()) +
  geom_label_repel(data = top_bugs,
                   mapping = aes(coef, -log(qval, 10), 
                                 label = feature),
                   show.legend = FALSE, 
                   size = 4,
                   min.segment.length = 0.2,
                   color = "black") 

###Figure 4B
#CAR
car_metabolites_data <- nafld_data %>%
  select(matches("^CAR \\d+:\\d+"), case, lean_vs_nonlean_case) %>%
  mutate(across(matches("^CAR \\d+:\\d+"), ~ {
    col_min_half <- min(.[. > 0], na.rm = TRUE) / 2
    .[. == 0] <- col_min_half
    log2(.)
  }))

boxplot_data_long_cars <- car_metabolites_data %>%
  pivot_longer(
    cols = matches("^CAR \\d+:\\d+"),  
    names_to = "Metabolite",
    values_to = "Abundance"
  )

boxplot_data_long_cars <- boxplot_data_long_cars %>%
  mutate(
    Chain_Length = as.numeric(gsub("CAR (\\d+):.*", "\\1", Metabolite))  
  ) %>%
  arrange(Chain_Length) %>%  # Arrange by chain length
  mutate(
    Metabolite = factor(Metabolite, levels = unique(Metabolite))  
  )

wilcoxon_results <- boxplot_data_long_cars %>%
  group_by(Metabolite) %>%
  summarize(
    p_value = wilcox.test(
      Abundance[case == 0], 
      Abundance[case == 1]
    )$p.value
  )

wilcoxon_results <- wilcoxon_results %>%
  mutate(
    q_value = p.adjust(p_value, method = "BH"),  # Benjamini-Hochberg FDR correction
    significance = case_when(
      q_value < 0.001 ~ "**",
      q_value < 0.01 ~ "*",
      TRUE ~ " "  # Not significant
    )
  )

boxplot_data_long_cars <- boxplot_data_long_cars %>%
  left_join(wilcoxon_results, by = "Metabolite")

ggplot(boxplot_data_long_cars, aes(x = Metabolite, y = Abundance)) +
  geom_boxplot(aes(fill = factor(case)), alpha = 0.6) +
  scale_fill_manual(name = "MASLD", values = c("#999999", "#E69F00"),
                    labels = c("Control", "MASLD")) +
  labs(x = "Metabolite", y = "log2(abundance)") +
  geom_text(
    data = wilcoxon_results,
    aes(x = Metabolite, y = max(boxplot_data_long_cars$Abundance, na.rm = TRUE), label = significance),
    inherit.aes = FALSE,
    size = 5,
    vjust = -0.5,
    fontface = "bold"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
    axis.text.y = element_text(size = 10)
  )

