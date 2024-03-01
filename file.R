# Загрузка библиотек
library(biomaRt)
library(clusterProfiler)
library(org.Dm.eg.db)

# Выбор базы данных и датасета
mart_obj <- useEnsembl(biomart = "genes", 
                      dataset = "dmelanogaster_gene_ensembl"
                      )
  
# Чтение csv-файла с пиками, необходимо указать актуальный путь
cols <- read.csv(file = 'e_coli_unique.csv', sep = ';')
  
# Получение Entrez ID для генов, ассоциированных с пиками
genes <- getBM(attributes = 'entrezgene_id',  
                filters = c('chromosome_name','start','end'),
                values = list(cols$seqnames, cols$start, cols$end),
                mart = mart_obj
                )

# Преобразование единственного столбца датафрейма genes в формат,
# воспринимаемый нижеиспользуемыми фунцкциями
pull <- dplyr::pull(genes, 'entrezgene_id')
  
# GO enrichment + график
enr_go <- enrichGO(pull, org.Dm.eg.db)
dotplot(enr_go)

# Пока график на экране, сохраняем его в pdf (1)
dev.print(pdf, 'enr_GO.pdf')
  
# KEGG enrichment + график
enr_kegg <- enrichKEGG(pull, 
                  organism = 'dme', 
                  keyType = 'ncbi-geneid', 
                  pvalueCutoff= 1, 
                  pAdjustMethod = "none", 
                  qvalueCutoff = 1
                  )
dotplot(enr_kegg)

# Пока график на экране, сохраняем его в pdf (2)
dev.print(pdf, 'enr_KEGG.pdf')
