library(org.Hs.eg.db)
library(clusterProfiler)
library(tidyverse)
library(ggplot2)
setwd("D:/Desktop/毕设/liver  cancer/5-originnetwork")

rm(list = ls())
CN <- read.table("D:/Desktop/毕设/liver  cancer/1_WCN/top10%_CoreNet.txt", stringsAsFactors = F)
protein_all <- dplyr::distinct(as.data.frame(c(CN[, 1], CN[, 2])))
colnames(protein_all) <- "GeneName"

# data(geneList, package = "DOSE") # Background gene set for enrichment analysis
# geneid <- bitr(target_symbol$GeneName, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = "org.Hs.eg.db")

# GO
ego <- enrichGO(gene = protein_all$GeneName, 
                # universe = names(geneList), # Background gene set, need geneid to use this parameter
                OrgDb = org.Hs.eg.db, 
                keyType = "SYMBOL",
                ont = "all", # One of CC, BP, MF, all
                pAdjustMethod = "BH", # Correction methods: "holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none"
                pvalueCutoff = 1, 
                qvalueCutoff = 1) # If using geneid, add a parameter: readable = TRUE

go <- as.data.frame(ego)
go_CC <- go[go$ONTOLOGY == "CC",][1:10, ]
go_MF <- go[go$ONTOLOGY == "MF",][1:10, ]
go_BP <- go[go$ONTOLOGY == "BP",][1:10, ]
go_enrich_df <- data.frame(ID = c(go_CC$ID, go_MF$ID, go_BP$ID), 
                           Description = c(go_CC$Description, go_MF$Description, go_BP$Description), 
                           GeneNumber = c(go_CC$Count, go_MF$Count, go_BP$Count), 
                           p.adjust = c(go_CC$p.adjust, go_MF$p.adjust, go_BP$p.adjust),
                           group = factor(c(rep("Cellular Component", 10), 
                                           rep("Molecular Function", 10), 
                                           rep("Biological Process", 10)), 
                                         levels = c("Cellular Component", "Molecular Function", "Biological Process")))
go_enrich_df$number <- factor(rev(1:nrow(go_enrich_df)))

color <- c("#8DA1CB", "#FD8D62", "#66C3A5")
ggplot(data = go_enrich_df, aes(x = number, y = GeneNumber, fill = group)) + 
  geom_bar(stat = "identity", width = 0.8) + coord_flip() + scale_fill_manual(values = color) + 
  scale_x_discrete(labels = rev(go_enrich_df$Description)) + xlab("GO term") + theme_bw() + 
  theme(axis.text = element_text(face = "bold", color = "gray50", size = 12), 
        legend.title = element_text(size = 14), legend.text = element_text(size = 13)) + 
  labs(title = "The Most Enriched GO Terms of All Proteins") + geom_text(aes(label = p.adjust), hjust = 1)

# KEGG
# In KEGG enrichment analysis, symbol does not support human species (hsa) and needs to be converted to geneid
geneid <- bitr(protein_all$GeneName, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = "org.Hs.eg.db")
ekegg <- enrichKEGG(gene = geneid$ENTREZID, 
                    organism = "hsa", # Shorthand for human species:"hsa"
                    pvalueCutoff = 1)
barplot(ekegg, showCategory = 20, title = "The Most Enriched KEGG pathway of All Proteins") + 
  geom_text(aes(label = p.adjust), hjust = 1, color = "white")