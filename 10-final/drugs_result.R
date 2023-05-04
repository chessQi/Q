library(tidyverse)
setwd("D:/Desktop/毕设/liver  cancer/10-final")

rm(list = ls())
proximity_1 <- read.csv("proximity2_result_1.csv", sep = ",", header = T, stringsAsFactors = F)
proximity_2 <- read.csv("proximity2_result_2.csv", sep = ",", header = T, stringsAsFactors = F)
proximity_3 <- read.csv("proximity2_result_3.csv", sep = ",", header = T, stringsAsFactors = F)

drug_1 <- proximity_1[which(proximity_1$zscore < -2 & proximity_1$pvalue <= 0.05), 1]
drug_2 <- proximity_2[which(proximity_2$zscore < -2 & proximity_2$pvalue <= 0.05), 1]
drug_3 <- proximity_3[which(proximity_3$zscore < -2 & proximity_3$pvalue <= 0.05), 1]

drugs <- intersect(drug_1, drug_2) %>% intersect(drug_3)

drugs_df <- data.frame()

for (drug in drugs) {
  drugs_df <- rbind(drugs_df, proximity_1[which(proximity_1$Drug == drug), ]) %>% 
    rbind(proximity_2[which(proximity_2$Drug == drug), ]) %>% 
    rbind(proximity_3[which(proximity_3$Drug == drug), ])
}

drugs_result <- aggregate(drugs_df[, -1], by = list(Drug=drugs_df$Drug), mean) %>% dplyr::arrange(zscore, pvalue)

write.csv(drugs_result, "Drugs_result.csv", row.names = F)