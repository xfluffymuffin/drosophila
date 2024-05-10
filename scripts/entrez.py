from Bio import Entrez


def get_records(query):
    # Загружаем записи из GenBank по запросу из переменной 'query'
    Entrez.email = "yakrit2013@yandex.ru"

    # Поиск в базе NCBI Gene, сохранение результатов в record.
    handle = Entrez.esearch(db="gene",
                            term=query,
                            retmax=10000000)
    rec = Entrez.read(handle)

    with open("entrez_id.csv", "a") as outfile:
        try:
            outfile.write(f"{query[1:].rsplit('[')[0]},{rec["IdList"][0]}\n")
        except IndexError:
            outfile.write(f"{query[1:].rsplit('[')[0]},NA\n")


def get_gene_names(file):
    genes_list = []
    with open(file, "r") as genes:
        for gene in genes:
            genes_list.append(gene.strip("\n"))

    with open("entrez_id.csv", "w") as file:
        file.write(f"gene_name,entrez_id\n")

    for name in genes_list:
        get_records(f"({name}[Gene Name]) AND Drosophila melanogaster[Organism]")



get_gene_names("gene_names.txt")
