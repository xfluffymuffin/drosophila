library(ChIPpeakAnno)

# csv-файлы на вход
s <- read.csv("Ecoli_3_vs_Control_1_peak_compare_table.tsv", 
              header = FALSE, sep = "\t")
snew <- read.csv("Mlut_5_vs_Control_1_peak_compare_table.tsv", 
                 header = FALSE, sep = "\t")

# перевод трех столбцов в векторы
chr <- as.vector(s[,1])[-1]
chromosome <- as.vector(snew[,1])[-1]

peak_start <- as.vector(s[,2])[-1]
start <- as.vector(snew[,2])[-1]

peak_end <- as.vector(s[,3])[-1]
end <- as.vector(snew[,3])[-1]

log_c <- as.vector(s[,13])[-1]
log_c2 <- as.vector(snew[,13])[-1]


# создание датафреймов для двух наборов пиков
df <- data.frame(chr, peak_start, peak_end, log_c)
deef <- data.frame(chromosome, start, end, log_c2)

# создание объектов класса GRanges
gr_e_coli <- makeGRangesFromDataFrame(unique(df))
gr_m_luteus <- makeGRangesFromDataFrame(unique(deef))

#modif1 <- rbind(as.data.frame(gr_e_coli), gr_e_coli$log_c)


# нахождение общих и уникальных пиков, построение диаграммы
ol <- findOverlapsOfPeaks(gr_e_coli, gr_m_luteus)
makeVennDiagram(ol, fill=c("#009E73", "#F0E442"), # circle fill color
                col=c("#D55E00", "#0072B2"), #circle border color
                cat.col=c("#D55E00", "#0072B2"))

# помещение пиков в отдельные файлы
Unique_e_coli <- as(ol$peaklist[["gr_e_coli"]], "data.frame")
Unique_e_coli <- apply(Unique_e_coli,2,as.character)
write.csv2(Unique_e_coli, "e_coli_unique.csv")

Unique_m_luteus <- as(ol$peaklist[["gr_m_luteus"]], "data.frame")
Unique_m_luteus <- apply(Unique_m_luteus,2,as.character)
write.csv2(Unique_m_luteus, "m_luteus_unique.csv")

Overlap <- as(ol$peaklist[["gr_e_coli///gr_m_luteus"]], "data.frame")
Overlap <- apply(Overlap,2,as.character)
write.csv2(Overlap, "overlap.csv")
