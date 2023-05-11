library(ne.PCA)
setwd("D://Desktop//±œ…Ë/liver  cancer/0_prep/")

rm(list = ls())
nodes <- read.csv("nodes.csv", sep = ",")
edges <- read.csv("links.csv", sep = ",")
ne.PCA(nodes, edges)
