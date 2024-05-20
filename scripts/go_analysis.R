library(org.Dm.eg.db)
library(clusterProfiler)


genes <- read.table("D:\\microRNA\\scripts\\Mluteus_down_microRNA_targets\\Mluteus_down_microRNA_targets.txt")

print(genes)

g_pull <- dplyr::pull(genes, 'V1')


go_an <- enrichGO(g_pull, org.Dm.eg.db, pvalueCutoff= 1,
                  pAdjustMethod = "none",
                  qvalueCutoff = 1)
dotplot(go_an)
