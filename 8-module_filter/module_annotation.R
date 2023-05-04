library(tidyverse)
library(do)
library(org.Hs.eg.db)
library(clusterProfiler)
library(ggplot2)
library(ggpubr)
setwd("D:/Desktop/毕设/liver  cancer/8-module_filter")

rm(list = ls())
protein_type <- read.csv("D:/Desktop/毕设/liver  cancer/4-match/protein_type_module.csv", header = T)
# target <- read.csv("../1-3_target_match/coretarget.csv", header = T)
# non_target <- read.csv("../1-3_target_match/non-coretarget.csv", header = T)
module_label <- read.table("D:/Desktop/毕设/liver  cancer/1_WCN/node_Module.txt", header = T)
TFCs <- read.table("D:/Desktop/毕设/liver  cancer/1_WCN/WCN.txt", header = T, stringsAsFactors = F)
edge_module <- read.table("D:/Desktop/毕设/liver  cancer/1_WCN/edge_Module.txt", header = T, stringsAsFactors = F)

# target$Type <- "T"
# non_target$Type <- "N"
names(module_label)[1] <- "GeneName"
# protein_class <- rbind(target, non_target)

layout <- merge(module_label, protein_type, by = "GeneName")
module_table <- as.data.frame(table(layout$module))


# Set up analysis and plot parameters
# Module to analysis
analysis_module <- unique(module_label$module)

# Calculate target proportion of each module
proportion <- module_table
for (i in proportion$Var1) {
  Tcount <- length(which(layout$module == i & layout$Type == "T"))
  Psum <- module_table[(which(module_table$Var1 == i)), 2]
  pp <- Tcount/Psum
  proportion[which(proportion$Var1 == i), 3] <- Tcount
  proportion[which(proportion$Var1 == i), 4] <- pp
}
names(proportion) <- c("Module", "Protein_Count", "Target_Count", "Proportion")

# # If the network is too large, the algorithm can be optimizated.
# TFCs_module <- data.frame()
# for (gr in 1:nrow(TFCs)) {
#   if (TFCs[gr, 1] %in% module_label$GeneName == T && 
#       TFCs[gr, 2] %in% module_label$GeneName == T && 
#       module_label[which(module_label$GeneName == TFCs[gr, 1]), 2] == 
#       module_label[which(module_label$GeneName == TFCs[gr, 2]), 2]) {
#     TFCs[gr, 4] <- module_label[which(module_label$GeneName == TFCs[gr, 1]), 2]
#     TFCs_module <- rbind(TFCs_module, TFCs[gr, ])
#   }
# }
# names(TFCs_module)[4] <- "module"

TFCs_module <- merge(edge_module, TFCs, by = c("name1", "name2"))
TFCs_mean <- c()
for (m in analysis_module) {
  TFCs_temp <- TFCs_module[which(TFCs_module$module == m), ]
  TFCs_mean_temp <- mean(TFCs_temp$TFCs)
  TFCs_mean <- c(TFCs_mean, TFCs_mean_temp)
}
proportion$TFCs_mean <- TFCs_mean

write.csv(proportion, "proportion_TFCs.csv", row.names = F)

# Plot the TFC score Difference between Modules

# pattern <- c("^1$:M1", "^2$:M2", "^3$:M3", "^4$:M4", "^5$:M5", "^6$:M6", "^7$:M7", 
#              "^8$:M8", "^9$:M9", "^10$:M10", "^11$:M11", "^12$:M12", "^13$:M13", 
#              "^14$:M14", "^15$:M15", "^16$:M16", "^17$:M17", "^18$:M18", 
#              "^19$:M19", "^20$:M20", "^21$:M21", "^22$:M22")

# Equivalent to:
rule_prefix <- "M"
pattern <- paste0("^", analysis_module, "$:", rule_prefix, analysis_module)

TFCs_module$module <- Replace(data = TFCs_module$module, pattern = pattern)
png("TFCs_difference.png")
ggplot(TFCs_module, aes(x = module, y = TFCs, fill = module)) + geom_boxplot() + 
  stat_compare_means(label = "p.signif", method = "t.test", ref.group = ".all.", hide.ns = TRUE) + 
  ggtitle("TFC score Difference between Modules")
dev.off()

# Caculate degree of each module.
# The degree between modules is defined as the sum of the degrees of all nodes in module A and all nodes in module B.
module_degree <- data.frame(module = analysis_module, freq = 0)
diff_module_link <- dplyr::setdiff(TFCs[, 1:2], TFCs_module[, 1:2])
# for (dgr in 1:nrow(diff_module_link)) {
#   if (diff_module_link[dgr, 1] %in% module_label$GeneName == T && 
#       diff_module_link[dgr, 2] %in% module_label$GeneName == T) {
#     mid1 <- module_label[which(module_label$GeneName == diff_module_link[dgr, 1]), 2]
#     mid2 <- module_label[which(module_label$GeneName == diff_module_link[dgr, 2]), 2]
#     rid <- which(module_degree$module == mid1 | module_degree$module == mid2)
#     module_degree[rid, 2] <- module_degree[rid, 2] + 1
#   }
# }
for (dgr in 1:nrow(diff_module_link)) {
  if (diff_module_link[dgr, 1] %in% module_label$GeneName == T && 
      diff_module_link[dgr, 2] %in% module_label$GeneName == T) {
    mid1 <- module_label[which(module_label$GeneName == diff_module_link[dgr, 1]), 2]
    mid2 <- module_label[which(module_label$GeneName == diff_module_link[dgr, 2]), 2]
    if (mid1 != mid2) {
      rid <- which(module_degree$module == mid1 | module_degree$module == mid2)
      module_degree[rid, 2] <- module_degree[rid, 2] + 1
    }
  }
}
module_degree$module <- Replace(data = module_degree$module, pattern = pattern)
write.csv(module_degree, "module_degree.csv", row.names = F)