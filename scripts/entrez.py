import glob
from Bio import Entrez
import pandas as pd
from pandas import read_excel


def get_records(query):
    Entrez.email = "yakrit2013@yandex.ru"

    # Поиск в базе NCBI Gene, сохранение результатов в rec
    handle = Entrez.esearch(db="gene",
                            term=query,
                            retmax=10000000)
    rec = Entrez.read(handle)

    with open(f"entrez_id_{outp}s.csv", "a") as outfile:
        try:
            outfile.write(f"{query[1:].rsplit('[')[0]},{rec["IdList"][0]}\n")
        except IndexError:
            add_id = pd.read_csv(r"D:\microRNA\scripts\all_problematic_genes_FIXED.csv")
            add_ids = []

            for index, row in add_id.iterrows():
                add_ids.append(f"{row['gene_name']}={row['entrez_id']}")

            for id in add_ids:
                if query[1:].rsplit('[')[0] in id:
                    outfile.write(f"{query[1:].rsplit('[')[0]},{id.split("=")[1]}\n")


def get_gene_names(inp):
    global outp
    genes_list = []
    df = pd.read_excel(inp)
    outp = inp.split("__")[1].rstrip(".xlsx")

    for index, row in df.iterrows():
        genes_list.append(row['Target gene'])

    with open(f"entrez_id_{outp}s.csv", "w") as file:
        file.write(f"gene_name,entrez_id\n")

    for name in genes_list:
        get_records(f"({name}[Gene Name]) AND Drosophila melanogaster[Organism]")

# Программа анализирует все *.xlsx файлы в указанной папке
target_scan_dir = r"D:\microRNA\analysis_filtered\TargetScan_data\Mluteus\down\*.xlsx"

for filepath in glob.iglob(target_scan_dir):
    get_gene_names(filepath)
