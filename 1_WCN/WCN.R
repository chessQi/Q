library(ne.PCA)
setwd("D://Desktop//����/liver  cancer/0_prep/")

rm(list = ls())
nodes <- read.csv("nodes.csv", sep = ",")
edges <- read.csv("links.csv", sep = ",")
ne.PCA(nodes, edges)
