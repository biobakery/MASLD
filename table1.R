#!/usr/bin/env Rscript
##################################################
#R program for creating Table 1
##################################################

library(tableone)
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
new_meta <- read_csv('input/table1variables.csv') %>%
  mutate(pkyr17 = ifelse(pkyr17==998,0,pkyr17)) %>%
  mutate(pkyr17 = ifelse(pkyr17==999,NA,pkyr17)) %>% #998 is never smokers; 999 is missing
  select(alias_id, pkyr17)
df_w_meta <- left_join(df_w_meta,new_meta,by="alias_id")

nafld_data<-df_w_meta %>% filter(cohort=="NAFLD" | case==0) %>%
  mutate(obesity = case_when(bmi17v >= 30 ~ 1, bmi17v <30 ~ 0)) %>%
  mutate(lean = case_when(bmi17v < 25 ~ 1, bmi17v >= 25 ~ 0)) %>%
  mutate(hbp17 = ifelse(hbp17 =="yes", 1, 0))

nafld_data_sens<-df_w_meta %>% filter(cohort=="NAFLD") %>%
  mutate(obesity = case_when(bmi17v >= 30 ~ 1, bmi17v <30 ~ 0)) %>%
  mutate(lean = case_when(bmi17v < 25 ~ 1, bmi17v >= 25 ~ 0)) %>%
  mutate(hbp17 = ifelse(hbp17 =="yes", 1, 0))

meta_data <- nafld_data %>% select(!starts_with('s__'))
meta_data_sens <- nafld_data_sens %>% select(!starts_with('s__'))

# by subject
table1 <- meta_data %>%
  select(
    alias_id,
    Age = age,
    NAFLD = case,
    Hypertension = hbp17, #high_blood_pressure
    Hypercholesterolemia = chol17, #high_cholesterol
    `Diabetes Mellitus` =  db17,
    `Pack-Years` = pkyr17,
    `Physical activity` = act17v,
    `Body Mass Index (BMI)` = bmi17v,
    `Postmenopausal Hormone Use` = nhor15,
    `Aspirin Use` = aspu15,
    AHEI = aheiv2010_15,
    `Total Calories` = calor15v,
    `Western Diet` = western_diet)

t1 <- CreateTableOne(
  data = table1,
  vars = c(
    "Age",
    "Hypertension",
    "Hypercholesterolemia",
    "Diabetes Mellitus",
    "Pack-Years",
    "Physical activity",
    "Body Mass Index (BMI)",
    "Postmenopausal Hormone Use",
    "Aspirin Use",
    "AHEI",
    "Total Calories",
    "NAFLD",
    "Western Diet"),
  strata = c("NAFLD"),
  factorVars = c(
    "Hypertension",
    "Hypercholesterolemia",
    "Diabetes Mellitus",
    "Postmenopausal Hormone Use",
    "Aspirin Use",
    "NAFLD")
)

table1.subject <- print(
  t1,
  showAllLevels = TRUE,
  quote = FALSE,
  noSpaces = TRUE,
  printToggle = FALSE,
  contDigits = 1
)

table1.subject

write.csv(table1.subject,
          file = file.path("output", "table1.csv"))

rm(table1.subject)

#sensitivity (only including NAFLD matched controls)
table1_sens <- meta_data_sens %>%
  select(
    alias_id,
    Age = age,
    NAFLD = case,
    Hypertension = hbp17,
    Hypercholesterolemia = chol17,
    `Diabetes Mellitus` =  db17,
    `Pack-Years` = pkyr17,
    `Physical activity` = act17v,
    `Body Mass Index (BMI)` = bmi17v,
    `Postmenopausal Hormone Use` = nhor15,
    `Aspirin Use` = aspu15,
    AHEI = aheiv2010_15,
    `Total Calories` = calor15v,
    `Western Diet` = western_diet)

t1_sens <- CreateTableOne(
  data = table1_sens,
  vars = c(
    "Age",
    "Hypertension",
    "Hypercholesterolemia",
    "Diabetes Mellitus",
    "Pack-Years",
    "Physical activity",
    "Body Mass Index (BMI)",
    "Postmenopausal Hormone Use",
    "Aspirin Use",
    "AHEI",
    "Total Calories",
    "NAFLD",
    "Western Diet"),
  strata = c("NAFLD"),
  factorVars = c(
    "Hypertension",
    "Hypercholesterolemia",
    "Diabetes Mellitus",
    "Postmenopausal Hormone Use",
    "Aspirin Use",
    "NAFLD")
)

table1.subject_sens <- print(
  t1_sens,
  showAllLevels = TRUE,
  quote = FALSE,
  noSpaces = TRUE,
  printToggle = FALSE,
  contDigits = 1
)

table1.subject_sens

write.csv(table1.subject_sens,
          file = file.path("output", "table1_sensitivity.csv"))


rm(table1.subject_sens)
