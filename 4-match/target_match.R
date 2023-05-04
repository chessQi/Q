library(tidyverse)
setwd("D://Desktop//毕设//liver  cancer//4-match")

rm(list = ls())
protein_net <- read.table("D:/Desktop/毕设/liver  cancer/1_WCN/top10%_CoreNet.txt", stringsAsFactors = F)
node_module <- read.table("D:/Desktop/毕设/liver  cancer/1_WCN/node_Module.txt", header = T, stringsAsFactors = F)
DB_info <- read.delim("D:/Desktop/毕设/liver  cancer/2_drugBank/DrugBank_TDI.csv", sep = ",", header = T, stringsAsFactors = F)
TTD_info <- read.delim("D:/Desktop/毕设/liver  cancer/3_TTD/TTD_TDI.csv", sep = ",", header = T, stringsAsFactors = F)[, -1:-2]

protein <- distinct(as.data.frame(c(protein_net[, 1], protein_net[, 2])))
#protein是所有的基因名称总和
colnames(protein)[1] <- "GeneName"

DB_merge <- merge(protein, DB_info, by = "GeneName", sort = T) %>% dplyr::select(c(2, 1, 3:7))
TTD_merge <- merge(protein, TTD_info, by = "GeneName", sort = T) %>% dplyr::select(c(2, 1, 3:7))

nDB <- 0
nTTD <- 0
target_list <- c()
for (i in protein$GeneName) {
  if (i %in% DB_merge$GeneName == T) {
    nDB <- nDB + 1
    target_list <- append(target_list, i)
  }
  else if (i %in% TTD_merge$GeneName == T) {
    nTTD <- nTTD + 1
    target_list <- append(target_list, i)
  }
}
#target_list是Genename的子集
target <- data.frame(GeneName = target_list)
non_target <- dplyr::setdiff(protein, target)

node_module <- dplyr::select(node_module, 1)
names(node_module) <- "GeneName"
target_module <- dplyr::intersect(target, node_module)
non_target_module <- dplyr::intersect(non_target, node_module)

target_drug <- rbind(merge(target, DB_merge[, 2:3]), merge(target, TTD_merge[, 2:3]))
target_drug_module <- rbind(merge(target_module, DB_merge[, 2:3]), merge(target_module, TTD_merge[, 2:3]))

protein_type <- data.frame(GeneName = c(target[, 1], non_target[, 1]), 
                           Type = c(rep("T", nrow(target)), rep("N", nrow(non_target))))
protein_type_module <- data.frame(GeneName = c(target_module[, 1], non_target_module[, 1]), 
                                  Type = c(rep("T", nrow(target_module)), rep("N", nrow(non_target_module))))

write.csv(DB_merge, "DrugBank_Target.csv", row.names = F)
write.csv(TTD_merge, "TTD_Target.csv", row.names = F)
write.csv(target, "coretarget.csv", row.names = F)
write.csv(non_target, "non-coretarget.csv", row.names = F)
write.csv(target_module, "coretarget_module.csv", row.names = F)
write.csv(non_target_module, "non-coretarget_module.csv", row.names = F)
write.csv(protein_type, "protein_type.csv", row.names = F)
write.csv(protein_type_module, "protein_type_module.csv", row.names = F)
write.csv(target_drug, "Target_Drug.csv", row.names = F)
write.csv(target_drug_module, "Target_Drug_module.csv", row.names = F)