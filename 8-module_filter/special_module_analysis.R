library(tidyverse)
library(do)
library(org.Hs.eg.db)
library(clusterProfiler)
library(ggplot2)
library(ggpubr)
setwd("D:/Desktop/毕设/liver  cancer/8-module_filter")

rm(list = ls())
module_label <- read.table("D:/Desktop/毕设/liver  cancer/1_WCN/node_Module.txt", header = T)


# Set up analysis and plot parameters
# Module to analysis
# analysis_module <- unique(module_label$module)
# or
analysis_module <- c(1,2,4,6,7,17)
# Number of Gene Ontology Enrichment Analysis
go_num <- 10
# Number of KEGG Pathway Enrichment Analysis
kegg_num <- 10


# Drawing area division
if (length(analysis_module) < 16) {
  pic_nrow <- floor(sqrt(length(analysis_module)))
  pic_ncol <- ceiling(length(analysis_module) / pic_nrow)
} else {
  pic_ncol <- 4
  pic_nrow <- ceiling(length(analysis_module) / pic_ncol)
}

# Enrichment analysis of each module
go_CC <- data.frame()
go_MF <- data.frame()
go_BP <- data.frame()
# About 7 mins for 20 modules
for (md in analysis_module) {
  ego_temp <- enrichGO(gene = module_label[which(module_label$module == md), 1], 
                       OrgDb = org.Hs.eg.db, 
                       keyType = "SYMBOL",
                       ont = "all", # One of CC, BP, MF, all
                       pAdjustMethod = "BH", # Correction methods: "holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none"
                       pvalueCutoff = 1, 
                       qvalueCutoff = 1) # If using geneid, add a parameter: readable = TRUE
  go_temp <- as.data.frame(ego_temp)
  go_temp$Module <- md
  go_CC_temp <- go_temp[go_temp$ONTOLOGY == "CC",][1:go_num, ]
  go_MF_temp <- go_temp[go_temp$ONTOLOGY == "MF",][1:go_num, ]
  go_BP_temp <- go_temp[go_temp$ONTOLOGY == "BP",][1:go_num, ]
  
  go_CC <- rbind(go_CC, go_CC_temp)
  go_MF <- rbind(go_MF, go_MF_temp)
  go_BP <- rbind(go_BP, go_BP_temp)
}

# write.csv(go_CC, "go_CC.csv")
# write.csv(go_MF, "go_MF.csv")
# write.csv(go_BP, "go_BP.csv")
# go_CC <- read.csv("go_CC.csv", header = T, row.names = 1, stringsAsFactors = F)
# go_MF <- read.csv("go_MF.csv", header = T, row.names = 1, stringsAsFactors = F)
# go_BP <- read.csv("go_BP.csv", header = T, row.names = 1, stringsAsFactors = F)

go_enrich_df <- data.frame(ID = c(go_CC$ID, go_MF$ID, go_BP$ID), 
                           Description = c(go_CC$Description, go_MF$Description, go_BP$Description), 
                           GeneNumber = c(go_CC$Count, go_MF$Count, go_BP$Count), 
                           group = factor(c(rep("Cellular Component", go_num*length(analysis_module)), 
                                            rep("Molecular Function", go_num*length(analysis_module)), 
                                            rep("Biological Process", go_num*length(analysis_module))), 
                                          levels = c("Cellular Component", "Molecular Function", "Biological Process")), 
                           p.adjust = c(go_CC$p.adjust, go_MF$p.adjust, go_BP$p.adjust), 
                           Module = c(go_CC$Module, go_MF$Module, go_BP$Module))
go_enrich_df <- na.omit(go_enrich_df)
go_enrich_rank <- dplyr::arrange(go_enrich_df, Module, group, p.adjust)
write.csv(go_enrich_rank, "GO_Enrichment_data.csv", row.names = F)

png("Special_Module_Enrichment_GO.png", width = (90 - 3*(pic_ncol - 1))*go_num*pic_ncol, height = 30*go_num*3*pic_nrow + 60)
go_title <- paste0("The Most Enriched GO Terms of Module ", paste0(analysis_module, collapse = ","))
color <- c("#8DA1CB", "#FD8D62", "#66C3A5")
pgo_list <- list()
pgo_index <- 0
for (mdid in analysis_module) {
  plotdata <- go_enrich_rank[which(go_enrich_rank$Module == mdid), ]
  plotdata$Description <- factor(plotdata$Description, levels = rev(plotdata$Description))
  pgo_index <- pgo_index + 1
  pgo_list[[pgo_index]] <- ggplot(data = plotdata, aes(x = Description, y = GeneNumber, fill = group)) + 
    geom_bar(stat = "identity", width = 0.8) + coord_flip() + scale_fill_manual(values = color) + 
    scale_x_discrete(labels = function(x) str_wrap(x, width = 6*go_num)) + 
    xlab(NULL) + ylab(NULL) + theme_bw() + 
    theme(axis.text = element_text(face = "bold", size = 15 + go_num - 10),
          legend.title = element_text(size = 14), legend.text = element_text(size = 13))
}

# pgo <- ggpubr::ggarrange(pgo_list[[1]], pgo_list[[2]], pgo_list[[3]], pgo_list[[4]], pgo_list[[5]], 
#                          pgo_list[[6]], pgo_list[[7]], pgo_list[[8]], pgo_list[[9]], pgo_list[[10]], 
#                          pgo_list[[11]], pgo_list[[12]], pgo_list[[13]], pgo_list[[14]], pgo_list[[15]], 
#                          pgo_list[[16]], pgo_list[[17]], pgo_list[[18]], pgo_list[[19]], pgo_list[[20]], 
#                          pgo_list[[21]], pgo_list[[22]], nrow = 6, ncol = 4, labels = c(1:22), 
#                          font.label = list(color = 'red'), common.legend = T)

# Equivalent to:
pgo <- do.call(ggpubr::ggarrange, c(pgo_list, list(nrow = pic_nrow, ncol = pic_ncol, labels = analysis_module, 
                                                   font.label = list(color = 'red'), common.legend = T)))
pgo <- annotate_figure(pgo, top = text_grob(go_title, size = 25), 
                       left = text_grob("Go Terms", size = 20, rot = 90), 
                       bottom = text_grob("GeneNumber", size = 20))
# ggsave报错，暂无解决方法
# ggsave("test.png", pgo, width = (90 - 3*(pic_ncol - 1))*go_num*pic_ncol, height = 30*go_num*3*pic_nrow + 60, limitsize = F)
print(pgo)
dev.off()


# KEGG
# In KEGG enrichment analysis, symbol does not support human species (hsa) and needs to be converted to geneid.
png("Special_Module_Enrichment_KEGG.png", width = 50*kegg_num*pic_ncol, height = 35*kegg_num*pic_nrow + 50)
kegg_title <- paste0("The Most Enriched KEGG Pathway of Module ", paste0(analysis_module, collapse = ","))
pkegg_list <- list()
pkegg_index <- 0
# About 6 mins for 20 modules
kegg_enrich_df <- data.frame()
for (mdid in analysis_module) {
  geneid <- bitr(module_label[which(module_label$module == mdid), 1], 
                 fromType = "SYMBOL", toType = "ENTREZID", 
                 OrgDb = "org.Hs.eg.db")
  ekegg <- enrichKEGG(gene = geneid$ENTREZID, 
                      organism = "hsa", # Shorthand for human species:"hsa"
                      pvalueCutoff = 1)
  kegg_temp <- ekegg@result[1:kegg_num, ]
  kegg_temp$Module <- mdid
  kegg_enrich_df <- rbind(kegg_enrich_df, kegg_temp)
  pkegg_index <- pkegg_index + 1
  pkegg_list[[pkegg_index]] <- barplot(ekegg, showCategory = kegg_num)
}
write.csv(kegg_enrich_df, "KEGG_Enrichment_data.csv", row.names = F)

# pk <- ggpubr::ggarrange(pkegg_list[[1]], pkegg_list[[2]], pkegg_list[[3]], pkegg_list[[4]], pkegg_list[[5]], 
#                         pkegg_list[[6]], pkegg_list[[7]], pkegg_list[[8]], pkegg_list[[9]], pkegg_list[[10]], 
#                         pkegg_list[[11]], pkegg_list[[12]], pkegg_list[[13]], pkegg_list[[14]], pkegg_list[[15]], 
#                         pkegg_list[[16]], pkegg_list[[17]], pkegg_list[[18]], pkegg_list[[19]], pkegg_list[[20]], 
#                         pkegg_list[[21]], pkegg_list[[22]], nrow = 6, ncol = 4, labels = c(1:22), 
#                         font.label = list(color = 'red'))

# Equivalent to:
pk <- do.call(ggpubr::ggarrange, c(pkegg_list, list(nrow = pic_nrow, ncol = pic_ncol, labels = analysis_module, 
                                                    font.label = list(color = 'red'))))
pk <- annotate_figure(pk, top = text_grob(kegg_title, size = 25))
print(pk)
dev.off()
