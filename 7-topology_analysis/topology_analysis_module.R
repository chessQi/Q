library(tidyverse)
library(ggplot2)
library(ggpubr)
library(do)
setwd("D:/Desktop/毕设/liver  cancer/7-topology_analysis")

rm(list = ls())
topology <- read.delim("module_topology_parameters.csv", sep = ",", header = T, stringsAsFactors = F)
target0 <- read.csv("D:/Desktop/毕设/liver  cancer/4-match/coretarget_module.csv", header = T, stringsAsFactors = F)
non_target0 <- read.csv("D:/Desktop/毕设/liver  cancer/4-match/non-coretarget.csv", header = T, stringsAsFactors = F)
module_label <- read.table("D:/Desktop/毕设/liver  cancer/1_WCN/node_Module.txt", header = T, stringsAsFactors = F)

topo_para <- c("name", "AverageShortestPathLength", "BetweennessCentrality", "ClosenessCentrality", 
               "ClusteringCoefficient", "Degree", "NeighborhoodConnectivity", "Radiality", "Stress", 
               "TopologicalCoefficient")
plot_para <- dplyr::select(topology, all_of(topo_para))
names(plot_para)[1] <- "GeneName"
names(module_label)[1] <- "GeneName"

target <- merge(target0, module_label, by = "GeneName")
non_target <- merge(non_target0, module_label, by = "GeneName")

target_topo <- merge(target, plot_para, by = "GeneName")
non_target_topo <- merge(non_target, plot_para, by = "GeneName")

# # Calculated the average topological parameters of target and non-target proteins
# target_topopara_mean <- c()
# ntarget_topopara_mean <- c()
# for (i in 3:ncol(target_topo)) {
#   target_topopara_mean <- c(target_topopara_mean, mean(target_topo[, i]))
#   ntarget_topopara_mean <- c(ntarget_topopara_mean, mean(non_target_topo[, i]))
# }
# target_non_mean <- as.data.frame(rbind(target_topopara_mean, ntarget_topopara_mean))
# colnames(target_non_mean) <- colnames(target_topo)[c(-1, -2)]
# write.csv(target_non_mean, "target_non_topopara_mean.csv")

all_topo <- rbind(target_topo, non_target_topo) %>% .[order(.$module), ]
# all_topo$module <- do::Replace(data = all_topo$module, 
#                                pattern = c("^1$:M1", "^2$:M2", "^3$:M3", "^4$:M4", "^5$:M5", "^6$:M6", "^7$:M7", 
#                                            "^8$:M8", "^9$:M9", "^10$:M10", "^11$:M11", "^12$:M12", "^13$:M13", 
#                                            "^14$:M14", "^15$:M15", "^16$:M16", "^17$:M17", "^18$:M18", 
#                                            "^19$:M19", "^20$:M20", "^21$:M21", "^22$:M22"))
# Equivalent to:
analysis_module <- unique(module_label$module)
rule_prefix <- "M"
pattern <- paste0("^", analysis_module, "$:", rule_prefix, analysis_module)
all_topo$module <- do::Replace(data = all_topo$module, pattern = pattern)

# Calculated the average topological parameters of all proteins in each module.
module <- unique(all_topo$module)
plotdata <- data.frame()
parameter_mean <- data.frame(module)
for (cn in colnames(all_topo)[c(-1, -2)]) {
  plotdata_part <- dplyr::select(all_topo, c(2, which(colnames(all_topo) == cn)))
  plotdata_part$Parameter <- cn
  colnames(plotdata_part)[2] <- "Value"
  plotdata <- rbind(plotdata, plotdata_part)
  
  parameter_mean_str <- c()
  for (m in module) {
    parameter_data_temp <- all_topo[which(all_topo$module == m), ]
    parameter_mean_temp <- mean(parameter_data_temp[, which(colnames(parameter_data_temp) == cn)])
    parameter_mean_str <- c(parameter_mean_str, parameter_mean_temp)
  }
  parameter_mean[, ncol(parameter_mean)+1] <- parameter_mean_str
}
colnames(parameter_mean)[-1] <- colnames(all_topo)[c(-1, -2)]
write.csv(parameter_mean, "module_topoparameter_mean.csv", row.names = F)

# # Plot the difference between the topological parameters of all modules
# png("all_topo_difference.png", height = 1080, width = 1920)
# ggplot(plotdata, aes(x = module, y = Value, fill = module)) + 
#   geom_boxplot() + facet_wrap(~Parameter, scale = "free") + 
#   stat_compare_means(label = "p.signif", method = "t.test", ref.group = ".all.", hide.ns = TRUE)
# dev.off()

###############################################################################
# But what I need is only "NeighborhoodConnectivity".
###############################################################################

png("NC_difference.png")
ggplot(plotdata[which(plotdata$Parameter == "NeighborhoodConnectivity"), ], 
       aes(x = module, y = Value, fill = module)) + geom_boxplot() + 
  stat_compare_means(label = "p.signif", method = "t.test", ref.group = ".all.", hide.ns = TRUE) + 
  ggtitle("NeighborhoodConnectivity Difference between Modules")
dev.off()


# # Plot the difference between the topological parameters of the target module
# special_module <- c("M2", "M7", "M13", "M14", "M15")
# special_module_topo <- data.frame()
# for (sm in special_module) {
#   special_module_topo <- rbind(special_module_topo, all_topo[which(all_topo$module == sm), ])
# }
# sm_plotdata <- data.frame()
# my_comparisons = list(c("M2","M7"), c("M2","M13"), c("M2","M14"), c("M2","M15"), 
#                       c("M7","M13"), c("M7","M14"), c("M7","M15"), 
#                       c("M13","M14"), c("M13","M15"), 
#                       c("M14","M15"))
# for (smcn in colnames(special_module_topo)[c(-1, -2)]) {
#   sm_plotdata_part <- dplyr::select(special_module_topo, c(2, which(colnames(special_module_topo) == smcn)))
#   sm_plotdata_part$Parameter <- smcn
#   colnames(sm_plotdata_part)[2] <- "Value"
#   sm_plotdata <- rbind(sm_plotdata, sm_plotdata_part)
# }
# png("special_module_topo_difference.png", height = 1080, width = 1920)
# ggplot(sm_plotdata, aes(x = module, y = Value, fill = module)) + 
#   geom_boxplot() + facet_wrap(~Parameter, scale = "free") +
#   stat_compare_means(comparisons = my_comparisons) + 
#   stat_compare_means(method = "anova")
# dev.off()

# # Plot the difference between the topological parameters of special modules and all other modules.
# special_non_topo <- all_topo
# special_non_topo[which(all_topo$module != special_module), 2] <- "Other_Mod"
# sm_non_plotdata <- data.frame()
# my_comparisons_avg = list(c("M2","Other_Mod"), c("M7","Other_Mod"), c("M13","Other_Mod"), 
#                           c("M14","Other_Mod"), c("M15","Other_Mod"))
# for (smncn in colnames(special_non_topo)[c(-1, -2)]) {
#   sm_non_plotdata_part <- dplyr::select(special_non_topo, c(2, which(colnames(special_non_topo) == smncn)))
#   sm_non_plotdata_part$Parameter <- smncn
#   colnames(sm_non_plotdata_part)[2] <- "Value"
#   sm_non_plotdata <- rbind(sm_non_plotdata, sm_non_plotdata_part)
# }
# png("special_non_topo_difference.png", height = 1080, width = 1920)
# ggplot(sm_non_plotdata, aes(x = module, y = Value, fill = module)) + 
#   geom_boxplot() + facet_wrap(~Parameter, scale = "free") +
#   stat_compare_means(comparisons = my_comparisons_avg) + 
#   stat_compare_means(method = "anova")
# dev.off()