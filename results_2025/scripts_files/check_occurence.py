import pandas as pd

file_path = r"D:\studying\mirna_proj\scripts_files\matches_mirbase.csv"

df = pd.read_csv(file_path, sep=";")

gene_ids = [
    "FBgn0040323", "FBgn0030310", "FBgn0035806", "FBgn0051217", "FBgn0039494", "FBgn0030051", "FBgn0030774",
    "FBgn0052383", "FBgn0052382", "FBgn0039102", "FBgn0003495", "FBgn0262473", "FBgn0026760", "FBgn0033402",
    "FBgn0003882", "FBgn0000250", "FBgn0260632", "FBgn0011274", "FBgn0038134", "FBgn0025574"]

# Фильтруем строки по gene ID и группируем данные
result = df[df["ensembl_gene_id"].isin(gene_ids)].groupby("ensembl_gene_id").agg(
    count=("ensembl_gene_id", "count"),
    mirna=("mirna", lambda x: ", ".join(set(x)))  # Собираем уникальные miRNA
).reset_index()

result.to_csv('result.csv')
print(result.to_string())

