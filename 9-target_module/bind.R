library(tidyverse)
setwd("D:/Desktop/毕设/liver  cancer/9-target_module")

rm(list = ls())
DB_info <- read.delim("D:/Desktop/毕设/liver  cancer/4-match/DrugBank_Target.csv", sep = ",", header = T, stringsAsFactors = F)
TTD_info <- read.delim("D:/Desktop/毕设/liver  cancer/4-match/TTD_Target.csv", sep = ",", header = T, stringsAsFactors = F)
module_label <- read.table("D:/Desktop/毕设/liver  cancer/1_WCN/node_Module.txt", header = T)
corenet <- read.table("D:/Desktop/毕设/liver  cancer/1_WCN/top10%_CoreNet.txt", header = F, stringsAsFactors = F)

colnames(module_label)[1] <- "GeneName"
# corenet <- corenet[, -3]

# Define modules to analysis
analysis_module <- c(2)

# Get duplicate drugs
DB_info$InChI <- gsub("InChI=", "", DB_info$InChI)
same <- dplyr::union(dplyr::intersect(DB_info$DrugName, TTD_info$DrugName), 
                     dplyr::intersect(DB_info$InChI, TTD_info$InChI))
for (i in same) {
  ttd_srow <- which(TTD_info$DrugName == i | TTD_info$InChI == i)
  db_srow <- which(DB_info$DrugName == i | DB_info$InChI == i)
  TTD_info$DrugID[ttd_srow] <- unique(DB_info$DrugID[db_srow])
}
TD_all <- dplyr::distinct(rbind(dplyr::select(DB_info, c(GeneName, DrugID)), 
                                dplyr::select(TTD_info, c(GeneName, DrugID))))


# ttd_row <- c()
# for (i in same) {
#   ttd_row <- append(ttd_row, which(TTD_info$DrugName == i))
# }
# 
# # Keep the duplicate drugs in DrugBank, and delete the ones in TTD
# TD_all <- dplyr::distinct(rbind(dplyr::select(DB_info, c(1, 3)), dplyr::select(TTD_info[-ttd_row, ], c(1, 3))))

module_protein <- data.frame()
for (gr in 1:nrow(corenet)) {
  if (corenet[gr, 1] %in% module_label$GeneName == T && 
      corenet[gr, 2] %in% module_label$GeneName == T) {
    module_protein <- rbind(module_protein, corenet[gr, ])
  }
}

# # Retain targets related to the analysis module
# rowid <- c()
# for (i in TD_all$GeneName) {
#   if (i %in% module_label[which(module_label$module %in% analysis_module), 1] == T) {
#     rowid <- c(rowid, which(TD_all$GeneName == i))
#   }
# }
# rowid <- unique(rowid)
# TD_net <- TD_all[rowid, ]

# Retain targets related to the analysis module
TD_net <- merge(module_label[which(module_label$module %in% analysis_module), ] %>% dplyr::select(1), 
                TD_all, by = "GeneName", all.x = T) %>% na.omit()

TD_net$label <- "D"
TD_net <- merge(TD_net, module_label, by = "GeneName")
write.csv(dplyr::select(TD_net, c(2:4)), "drug_label.csv", row.names = F, quote = F)
write.csv(TD_net, "module_DTN.csv", row.names = F, quote = F)

# Combine drugs on the original module network to build a new net
# In order to do rbind, name the colum 2 "DrugID"
# colnames(corenet_module) <- c("GeneName", "DrugID")
# newnet <- rbind(corenet_module, TD_net)
# 
# TD_all$label <- "D"
# write.csv(dplyr::select(TD_all, c(2, 3)), "drug_label.csv", row.names = F)
# write.table(newnet, "module_DTN.txt", row.names = F, quote = F)