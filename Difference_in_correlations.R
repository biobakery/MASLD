##################################################
#R program for creating difference in correlation figures
#Extended Data Figure 3A & Extended Data Figure 4
##################################################

library(tibble)
library(dplyr)
library(tidyr)
library(ggplot2)
library(ggcorrplot)
library(pheatmap)
library(gplots)
library(plotrix)

#Read in the two files from halla
lean_association_oral<-read.table('all_associations_oral_lean.txt', header=T, sep = "")
nonlean_association_oral<-read.table('all_associations_oral_nonlean.txt', header=T, sep = "")
#join the datasets and get absolute difference and signs
lean_nonlean_merged_oral <- left_join(nonlean_association_oral,lean_association_oral,by=c("X_features","Y_features")) %>%
  mutate(absdiff=abs(association.x-association.y)) %>% mutate(sign=association.x*association.y) %>%
  filter(!is.na(absdiff)) %>% 
  mutate(group=case_when(association.x>=0 & association.y>=0 ~"+,+", association.x<=0 & association.y<=0 ~"-,-",
                         association.x>=0 & association.y<=0 ~"+,-", association.x<=0 & association.y>=0 ~"-,+")) %>%
  mutate(groupinnum=case_when(association.x>=0 & association.y>=0 ~1, association.x<=0 & association.y<=0 ~2,
                              association.x>=0 & association.y<=0 ~3, association.x<=0 & association.y>=0 ~4)) 

#significance calculation
library(psych)
for(i in 1:nrow(lean_nonlean_merged_oral)) 
{ 
  summary_test<-r.test(n = 174, r12 = lean_nonlean_merged_oral$association.x[i], 
                       n2 = 37, r34 = lean_nonlean_merged_oral$association.y[i])
  lean_nonlean_merged_oral$difference_p[i]<-as.numeric(summary_test$p)
}

lean_nonlean_merged_oral$fdr_qval<-p.adjust(lean_nonlean_merged_oral$difference_p,method="fdr")
to_spread<-lean_nonlean_merged_oral %>% select(X_features,Y_features,groupinnum)
to_cluster <- spread(to_spread, Y_features, groupinnum) %>% column_to_rownames("X_features") %>% as.matrix()
rowclust1 = hclust(dist(to_cluster))
colclust1 = hclust(dist(t(to_cluster)))
clustered_df = to_cluster[rowclust1$order,colclust1$order]

lean_nonlean_merged_oral_sig<-lean_nonlean_merged_oral %>% filter(difference_p<0.05) 

###Extended Data Figure 3A
#Oral-typical bacteria and MASLD-associated metabolites
ggplot(lean_nonlean_merged_oral, aes(x = X_features, y = Y_features, color = group)) + 
  geom_count(aes(size = absdiff)) +
  scale_size_area(max_size = 3) + 
  coord_flip() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size = 12),
        axis.text.y = element_text(size = 12)) + 
  scale_y_discrete(limits = colnames(clustered_df)) + 
  scale_x_discrete(limits = rownames(clustered_df))

###Extended Data Figure 4
#Bacteria and acylcarnitines
#Read in the two files from halla
lean_association_oral<-read.table('all_associations_CAR_lean.txt', header=T, sep = "")
nonlean_association_oral<-read.table('all_associations_CAR_nonlean.txt', header=T, sep = "")
#join the datasets and get absolute difference and signs
lean_nonlean_merged_oral <- left_join(nonlean_association_oral,lean_association_oral,by=c("X_features","Y_features")) %>%
  mutate(absdiff=abs(association.x-association.y)) %>% mutate(sign=association.x*association.y) %>%
  filter(!is.na(absdiff)) %>% 
  mutate(group=case_when(association.x>=0 & association.y>=0 ~"+,+", association.x<=0 & association.y<=0 ~"-,-",
                         association.x>=0 & association.y<=0 ~"+,-", association.x<=0 & association.y>=0 ~"-,+")) %>%
  mutate(groupinnum=case_when(association.x>=0 & association.y>=0 ~1, association.x<=0 & association.y<=0 ~2,
                              association.x>=0 & association.y<=0 ~3, association.x<=0 & association.y>=0 ~4)) 

##filtering so that there are at least 4 big absolute difference in correlation (absdiff > 0.3) for species
threshold <- 0.3
filtered_features <- lean_nonlean_merged_oral %>%
  filter(absdiff > threshold) %>%
  count(X_features) %>%
  filter(n >= 4) %>%
  pull(X_features)
lean_nonlean_filtered <- lean_nonlean_merged_oral %>%
  filter(X_features %in% filtered_features)

ggplot(lean_nonlean_filtered, aes(x = X_features, y = Y_features, color = group)) + 
  geom_count(aes(size = absdiff)) +
  scale_size_area(max_size = 3) + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size = 12),
        axis.text.y = element_text(size = 12)) + 
  scale_y_discrete(limits = colnames(clustered_df)) + 
  scale_x_discrete(limits = intersect(rownames(clustered_df), filtered_features))
