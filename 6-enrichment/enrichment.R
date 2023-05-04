library(clusterProfiler)
library(org.Hs.eg.db)
library(ggplot2)
setwd("D:/Desktop/毕设/liver  cancer/6-enrichment")

rm(list = ls())
target_symbol <- read.table("D:/Desktop/毕设/liver  cancer/4-match/coretarget_module.csv", header = T, stringsAsFactors = F)
non_target_symbol <- read.table("D:/Desktop/毕设/liver  cancer/4-match/non-coretarget_module.csv", header = T, stringsAsFactors = F)

# data(geneList, package = "DOSE") # Background gene set for enrichment analysis
# geneid <- bitr(target_symbol$GeneName, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = "org.Hs.eg.db")

gene_list <- list(target_symbol, non_target_symbol)

go_title_list <- list("The Most Enriched GO Terms of Target", "The Most Enriched GO Terms of Non-Target")
go_fname_list <- list("GO_target.png", "GO_non_target.png")

# GO
for (i in 1:2) {
  ego <- enrichGO(gene = gene_list[[i]]$GeneName, 
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
                             group = factor(c(rep("Cellular Component", 10), 
                                              rep("Molecular Function", 10), 
                                              rep("Biological Process", 10)), 
                                            levels = c("Cellular Component", "Molecular Function", "Biological Process")))
  go_enrich_df$number <- factor(rev(1:nrow(go_enrich_df)))
  
  color <- c("#8DA1CB", "#FD8D62", "#66C3A5")
  png(go_fname_list[[i]], width = 1540, height = 850)
  print(ggplot(data = go_enrich_df, aes(x = number, y = GeneNumber, fill = group)) + 
          geom_bar(stat = "identity", width = 0.8) + coord_flip() + scale_fill_manual(values = color) + 
          scale_x_discrete(labels = rev(go_enrich_df$Description)) + xlab("GO term") + theme_bw() + 
          theme(axis.text = element_text(face = "bold", color = "gray50", size = 12), 
                legend.title = element_text(size = 14), legend.text = element_text(size = 13)) + 
          labs(title = go_title_list[[i]]))
  dev.off()
}

# KEGG
# In KEGG enrichment analysis, symbol does not support human species (hsa) and needs to be converted to geneid.
kegg_title_list <- list("The Most Enriched KEGG pathway of Target", "The Most Enriched KEGG pathway of Non-Target")
kegg_fname_list <- list("KEGG_target.png", "KEGG_non_target.png")
for (j in 1:2) {
  geneid <- bitr(gene_list[[j]]$GeneName, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = "org.Hs.eg.db")
  ekegg <- enrichKEGG(gene = geneid$ENTREZID, 
                      organism = "hsa", # Shorthand for human species:"hsa"
                      pvalueCutoff = 1)
  png(kegg_fname_list[[j]], width = 1250, height = 710)
  print(barplot(ekegg, showCategory = 20, title = kegg_title_list[[j]], font.size = 18))
  dev.off()
}
