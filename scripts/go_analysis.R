library(org.Dm.eg.db)
library(clusterProfiler)


genes <- read.table("D:\\microRNA\\scripts\\Mluteus_up_microRNA_targets\\Mluteus_up_microRNA_targets.txt")

g_pull <- dplyr::pull(genes, 'V1')


go_an <- enrichGO(g_pull, org.Dm.eg.db, 
                  pvalueCutoff= 0.05,
                  pAdjustMethod = "none",
                  qvalueCutoff = 1)


kegg <- enrichKEGG(g_pull,
                          organism = 'dme',
                          keyType = 'ncbi-geneid',
                          pvalueCutoff= 0.05,
                          pAdjustMethod = "none",
                          qvalueCutoff = 1)

# Show results in a separate window
barplot(go_an)
barplot(kegg)

# Save GO analysis results to a *.csv file
write.csv(go_an[,2:9], "D:\\microRNA\\scripts\\Mluteus_up_microRNA_targets\\Mluteus_up_microRNA_targets_GO.csv")

# Save KEGG analysis results to a *.csv file
write.csv(kegg[,1:11], "D:\\microRNA\\scripts\\Mluteus_up_microRNA_targets\\Mluteus_up_microRNA_targets_KEGG.csv")
