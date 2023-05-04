library(tidyverse)
library(ggplot2)
library(ggpubr)
setwd("D:/Desktop/毕设/liver  cancer/5-originnetwork")

rm(list = ls())
topo_para <- c("name", "AverageShortestPathLength", "BetweennessCentrality", "ClosenessCentrality", 
               "ClusteringCoefficient", "Degree", "NeighborhoodConnectivity", "Radiality", "Stress", 
               "TopologicalCoefficient")
topology <- read.delim("original_topology_parameters.csv", sep = ",", header = T, stringsAsFactors = F) %>% 
  dplyr::select(all_of(topo_para))
target <- read.csv("D:/Desktop/毕设/liver  cancer/4-match/coretarget.csv", header = T, stringsAsFactors = F)
non_target <- read.csv("D:/Desktop/毕设/liver  cancer/4-match/non-coretarget.csv", header = T, stringsAsFactors = F)

names(topology)[1] <- "GeneName"
target_topo <- merge(target, topology, by = "GeneName", all.x = T)
non_target_topo <- merge(non_target, topology, by = "GeneName", all.x = T)

target_topo$GeneName <- "Target"
non_target_topo$GeneName <- "Non-Target"
type_bind <- rbind(target_topo, non_target_topo)
colnames(type_bind)[1] <- "Type"

plotdata <- data.frame()
for (cn in colnames(type_bind)[-1]) {
  plotdata_part <- dplyr::select(type_bind, c(1, which(colnames(type_bind) == cn)))
  plotdata_part$Parameter <- cn
  colnames(plotdata_part)[2] <- "Value"
  plotdata <- rbind(plotdata, plotdata_part)
}

png("original_topology_analysis_boxplot.png", width = 1507, height = 822)
ggplot(plotdata, aes(x = Parameter, y = Value, fill = Type))+ geom_boxplot() + facet_wrap(~Parameter, scale = "free") + 
  stat_compare_means(method = "t.test")
dev.off()

# ggplot(plotdata, aes(x = Type, y = Value, fill = Parameter))+ geom_boxplot() + facet_wrap(~Parameter, scale = "free") + 
#   stat_compare_means(method = "t.test")
