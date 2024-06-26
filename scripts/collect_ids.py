import glob
import pandas as pd


def collect_id(path):
    with open(f"{scan_dir.split('\\')[3]}_dupl.txt", "a") as ids:
        for line in pd.read_csv(path)["entrez_id"]:
            ids.write(str(line) + "\n")


def del_duplicates():
    clean_list = []

    with open(f"{scan_dir.split('\\')[3]}_dupl.txt", "r") as dupl:
        for line in dupl.readlines():
            clean_list.append(line.strip("\n"))

    with open(f"{scan_dir.split('\\')[3]}.txt", "w") as fin:
        for id in set(clean_list):
            fin.write(id+"\n")


scan_dir = r"D:\microRNA\scripts\Ecoli_down_microRNA_targets\*.csv"

for filepath in glob.iglob(scan_dir):
    collect_id(filepath)

del_duplicates()
