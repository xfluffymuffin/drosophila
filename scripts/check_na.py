# This program checks for NA values in the ID column and prints corresponding gene names.
# Those gene names can be used for manual ID search

import os
import pandas as pd
from pandas import read_csv


def na_check(file):
    na = []
    df = pd.read_csv(file)
    sp = ["E.coli" if "Ecoli" in file else "M.luteus"][0]
    reg = ["down" if "down" in file else "up"][0]
    mi = file.split("id_")[1].split(".pre")[0]

    for index, row in df[df.isna().any(axis=1)].iterrows():
        na.append(row["gene_name"])

    with open("na_ids.txt", "a") as na_ids:
        if na:
            na_ids.write(f"Species: {sp}\nRegulation type: {reg}\nmiRNA: {mi}\nPath to file: {file}\n"
                         f"IDs missing for following genes:\n\t{'\n\t'.join(na)}\n\n")
        else:
            na_ids.write(f"Species: {sp}\nRegulation type: {reg}\nmiRNA: {mi}\nPath to file: {file}\n"
                         f"No IDs missing\n\n")


def list_no_duplicates():
    with open("na_ids.txt", "r") as na_ids:
        genes_list = []
        lines = na_ids.readlines()

        with open("all_problematic_genes.txt", "w") as all_prob:
            for line in lines:
                if "\t" in line:
                    genes_list.append(line.strip("\t\n"))

            for gene in set(genes_list):
                all_prob.write(gene + "\n")

    # Checking if our list of "problematic" genes contains duplicates (it does); conversion to set removes them
    # print(len(id_list) != len(set(id_list)))


for root, dirs, files in os.walk(r"D:\microRNA\scripts\Ecoli_down_microRNA_targets"):
    for file in files:
        na_check(rf"{root}\{file}")

list_no_duplicates()
