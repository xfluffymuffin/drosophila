# Load necessary libraries
library(biomaRt)
library(clusterProfiler)
library(org.Dm.eg.db)
library(enrichplot)


load_data <- function(path_to_all_pathogen_peaks, path_to_unique_pathogen_peaks) 
{
  # Choose the database and the dataset
  mart_obj <- useEnsembl(biomart = "genes",
                         dataset = "dmelanogaster_gene_ensembl"
  )
  # Read the file with unique peaks
  cols <- read.csv(file = path_to_unique_pathogen_peaks, sep = ';')
  
  # Read the file with ALL peaks
  s <- read.csv(path_to_all_pathogen_peaks, header = FALSE, sep = "\t")
  
  # Convert first 4 columns to vectors
  chr <- as.vector(s[,1])[-1]
  peak_start <- as.vector(s[,2])[-1]
  peak_end <- as.vector(s[,3])[-1]
  log_c <- as.vector(s[,13])[-1]

  # Create dataframes for all single pathogen peaks
  df <- data.frame(chr, peak_start, peak_end, log_c)
  
  # Extract only up- and down-regulated peaks
  up_reg <-  df[df$log_c > 0,]
  down_reg <- df[df$log_c < 0,]
  
  # Filter up- and down-regulated peaks by uniqueness
  up_reg_unique <- up_reg[up_reg$peak_start %in% cols$start,]
  down_reg_unique <- down_reg[down_reg$peak_start %in% cols$start,]
  
  # Preparing unique up-regulated peaks
  up_genes <- getBM(attributes = 'entrezgene_id',
                    filters = c('chromosome_name','start','end'),
                    values = list(up_reg_unique$chr, 
                                  up_reg_unique$peak_start, 
                                  up_reg_unique$peak_end),
                    mart = mart_obj
  )
  pull_up <<- dplyr::pull(up_genes, 'entrezgene_id')
  
  # Preparing unique down-regulated peaks
  down_genes <- getBM(attributes = 'entrezgene_id',
                      filters = c('chromosome_name','start','end'),
                      values = list(down_reg_unique$chr, 
                                    down_reg_unique$peak_start, 
                                    down_reg_unique$peak_end),
                      mart = mart_obj
  )
  pull_down <<- dplyr::pull(down_genes, 'entrezgene_id')
  
  # Preparing all unique peaks of the pathogen
  genes <- getBM(attributes = 'entrezgene_id',
                 filters = c('chromosome_name','start','end'),
                 values = list(cols$seqnames, cols$start, cols$end),
                 mart = mart_obj
  ) 
  
  pull <<- dplyr::pull(genes, 'entrezgene_id')
}

#  Analyse data for one of the two pathogens
#load_data('Ecoli_3_vs_Control_1_peak_compare_table.tsv', 'e_coli_unique.csv')
load_data('Mlut_5_vs_Control_1_peak_compare_table.tsv', 'm_luteus_unique.csv')



# ENRICHMENT ANALYSIS



# 1) Analyse ALL unique peaks of the chosen pathogen


  # a) GO enrichment
enr_go <- enrichGO(pull, org.Dm.eg.db)
dotplot(enr_go)
dev.print(pdf, 'enr_GO_all_unique_peaks_dotplot.pdf')

barplot(enr_go)
dev.print(pdf, 'enr_GO_all_unique_peaks_barplot.pdf')

enr_go_pw <- pairwise_termsim(enr_go)
emapplot(enr_go_pw)
dev.print(pdf, 'enr_GO_all_unique_peaks_emapplot.pdf')

cnetplot(enr_go)
dev.print(pdf, 'enr_GO_all_unique_peaks_cnetplot.pdf')


  # b) KEGG enrichment
enr_kegg <- enrichKEGG(pull,
                       organism = 'dme',
                       keyType = 'ncbi-geneid',
                       pvalueCutoff= 1,
                       pAdjustMethod = "none",
                       qvalueCutoff = 1
                       )
dotplot(enr_kegg)
dev.print(pdf, 'enr_KEGG_all_unique_peaks_dotplot.pdf')

barplot(enr_kegg)
dev.print(pdf, 'enr_KEGG_all_unique_peaks_barplot.pdf')

enr_kegg_pw <- pairwise_termsim(enr_kegg)
emapplot(enr_kegg_pw)
dev.print(pdf, 'enr_KEGG_all_unique_peaks_emapplot.pdf')

cnetplot(enr_kegg)
dev.print(pdf, 'enr_KEGG_all_unique_peaks_cnetplot.pdf')


  # c) WikiPathways enrichment
enr_wiki <- enrichWP(pull,
                     organism = 'Drosophila melanogaster',
                     pvalueCutoff= 1,
                     pAdjustMethod = "none",
                     qvalueCutoff = 1
                     )
dotplot(enr_wiki)
dev.print(pdf, 'enr_WP_all_unique_peaks_dotplot.pdf')

barplot(enr_wiki)
dev.print(pdf, 'enr_WP_all_unique_peaks_barplot.pdf')

enr_wiki_pw <- pairwise_termsim(enr_wiki)
emapplot(enr_wiki_pw)
dev.print(pdf, 'enr_WP_all_unique_peaks_emapplot.pdf')

cnetplot(enr_wiki)
dev.print(pdf, 'enr_WP_all_unique_peaks_cnetplot.pdf')


  # d) compareCluster
comp <- compareCluster(as.data.frame(pull), fun = "enrichGO", OrgDb = org.Dm.eg.db)
dotplot(comp)
dev.print(pdf, 'compareCluster_all_unique_peaks_dotplot.pdf')

comp_sm <- pairwise_termsim(comp)
emapplot(comp_sm)
dev.print(pdf, 'compareCluster_all_unique_peaks_emapplot.pdf')

cnetplot(comp)
dev.print(pdf, 'compareCluster_all_unique_peaks_cnetplot.pdf')



# 2) Analyse all up-regulated peaks of the chosen pathogen


# a) GO enrichment
enr_go_up <- enrichGO(pull_up, org.Dm.eg.db)
dotplot(enr_go_up)
dev.print(pdf, 'enr_GO_UP_unique_peaks_dotplot.pdf')

barplot(enr_go_up)
enr_go_pw_up <- pairwise_termsim(enr_go_up)
dev.print(pdf, 'enr_GO_UP_unique_peaks_barplot.pdf')

emapplot(enr_go_pw_up)
dev.print(pdf, 'enr_GO_UP_unique_peaks_emapplot.pdf')

cnetplot(enr_go_up)
dev.print(pdf, 'enr_GO_UP_unique_peaks_cnetplot.pdf')


# b) KEGG enrichment
enr_kegg_up <- enrichKEGG(pull_up,
                       organism = 'dme',
                       keyType = 'ncbi-geneid',
                       pvalueCutoff= 1,
                       pAdjustMethod = "none",
                       qvalueCutoff = 1
                       )
dotplot(enr_kegg_up)
dev.print(pdf, 'enr_KEGG_UP_unique_peaks_dotplot.pdf')

barplot(enr_kegg_up)
dev.print(pdf, 'enr_KEGG_UP_unique_peaks_barplot.pdf')

enr_kegg_pw_up <- pairwise_termsim(enr_kegg_up)
emapplot(enr_kegg_pw_up)
dev.print(pdf, 'enr_KEGG_UP_unique_peaks_emapplot.pdf')

cnetplot(enr_kegg_up)
dev.print(pdf, 'enr_KEGG_UP_unique_peaks_cnetplot.pdf')


# c) WikiPathways enrichment
enr_wiki_up <- enrichWP(pull_up,
                     organism = 'Drosophila melanogaster',
                     pvalueCutoff= 1,
                     pAdjustMethod = "none",
                     qvalueCutoff = 1
                     )
dotplot(enr_wiki_up)
dev.print(pdf, 'enr_WP_UP_unique_peaks_dotplot.pdf')

barplot(enr_wiki_up)
dev.print(pdf, 'enr_WP_UP_unique_peaks_barplot.pdf')

enr_wiki_pw_up <- pairwise_termsim(enr_wiki_up)
emapplot(enr_wiki_pw_up)
dev.print(pdf, 'enr_WP_UP_unique_peaks_emapplot.pdf')

cnetplot(enr_wiki_up)
dev.print(pdf, 'enr_WP_UP_unique_peaks_cnetplot.pdf')


# d) compareCluster
comp_up <- compareCluster(as.data.frame(pull_up), fun = "enrichGO", OrgDb = org.Dm.eg.db)
dotplot(comp_up)
dev.print(pdf, 'compareCluster_UP_unique_peaks_dotplot.pdf')

comp_sm_up <- pairwise_termsim(comp_up)
emapplot(comp_sm_up)
dev.print(pdf, 'compareCluster_UP_unique_peaks_emapplot.pdf')

cnetplot(comp_up)
dev.print(pdf, 'compareCluster_UP_unique_peaks_cnetplot.pdf')



# 3) Analyse all down-regulated peaks of the chosen pathogen


# a) GO enrichment
enr_go_down <- enrichGO(pull_down, org.Dm.eg.db)
dotplot(enr_go_down)
dev.print(pdf, 'enr_GO_DOWN_unique_peaks_dotplot.pdf')

barplot(enr_go_down)
dev.print(pdf, 'enr_GO_DOWN_unique_peaks_barplot.pdf')

enr_go_pw_down <- pairwise_termsim(enr_go_down)
emapplot(enr_go_pw_down)
dev.print(pdf, 'enr_GO_DOWN_unique_peaks_emapplot.pdf')

cnetplot(enr_go_down)
dev.print(pdf, 'enr_GO_DOWN_unique_peaks_cnetplot.pdf')


# b) KEGG enrichment
enr_kegg_down <- enrichKEGG(pull_down,
                          organism = 'dme',
                          keyType = 'ncbi-geneid',
                          pvalueCutoff= 1,
                          pAdjustMethod = "none",
                          qvalueCutoff = 1
                          )
dotplot(enr_kegg_down)
dev.print(pdf, 'enr_KEGG_DOWN_unique_peaks_dotplot.pdf')

barplot(enr_kegg_down)
dev.print(pdf, 'enr_KEGG_DOWN_unique_peaks_barplot.pdf')

enr_kegg_pw_down <- pairwise_termsim(enr_kegg_down)
emapplot(enr_kegg_pw_down)
dev.print(pdf, 'enr_KEGG_DOWN_unique_peaks_emapplot.pdf')

cnetplot(enr_kegg_down)
dev.print(pdf, 'enr_KEGG_DOWN_unique_peaks_ccnetplot.pdf')


# c) WikiPathways enrichment
enr_wiki_down <- enrichWP(pull_down,
                        organism = 'Drosophila melanogaster',
                        pvalueCutoff= 1,
                        pAdjustMethod = "none",
                        qvalueCutoff = 1
                        )
dotplot(enr_wiki_down)
dev.print(pdf, 'enr_WP_DOWN_unique_peaks_dotplot.pdf')

barplot(enr_wiki_down)
dev.print(pdf, 'enr_WP_DOWN_unique_peaks_barplot.pdf')

enr_wiki_pw_down <- pairwise_termsim(enr_wiki_down)
emapplot(enr_wiki_pw_down)
dev.print(pdf, 'enr_WP_DOWN_unique_peaks_emapplot.pdf')

cnetplot(enr_wiki_down)
dev.print(pdf, 'enr_WP_DOWN_unique_peaks_cnetplot.pdf')


# d) compareCluster
comp_down <- compareCluster(as.data.frame(pull_down), fun = "enrichGO", OrgDb = org.Dm.eg.db)
dotplot(comp_down)
dev.print(pdf, 'compareCluster_DOWN_unique_peaks_dotplot.pdf')

comp_sm_down <- pairwise_termsim(comp_down)
emapplot(comp_sm_down)
dev.print(pdf, 'compareCluster_DOWN_unique_peaks_emapplot.pdf')

cnetplot(comp_down)
dev.print(pdf, 'compareCluster_DOWN_unique_peaks_cnetplot.pdf')
