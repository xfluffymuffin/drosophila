import pandas as pd
import urllib.request
import gzip
import shutil
import argparse
import os
from tqdm import tqdm

# get a list with FlyBase ID (genes\transcripts) we'd like to check + a separate miRNAs list -->
# --> create files with them --> feed them to this script --> find matches -->
# --> maybe filter/prioritize by interaction score


def download_progress(block_num, block_size, total_size):

    downloaded = block_num * block_size
    percent = min(downloaded / total_size * 100, 100)
    print(f"\rDownloading requested database: {percent:.2f}% ({downloaded}/{total_size} bytes)", end="")


def download_db(db):

    script_dir = os.path.dirname(os.path.abspath(__file__))

    download_dir = os.path.join(script_dir, "downloads")
    os.makedirs(download_dir, exist_ok=True)

    db_file = f"mres_drosophila.microT.{db}.txt"
    save_path = os.path.join(download_dir, db_file)
    if os.path.exists(f"{save_path}.gz"):
        pass
    else:
        urllib.request.urlretrieve(
            f"https://dianalab.e-ce.uth.gr/microt_webserver_downloads/flat_files/{db_file}.gz",
            f"{save_path}.gz", download_progress)
        print("\nDownload complete.")

    with gzip.open(f"{save_path}.gz", "rb") as f_in:
        with open(save_path, "wb") as f_out:
            shutil.copyfileobj(f_in, f_out)


def find_matches(db, fbid_list, mirna_list, output):

    script_dir = os.path.dirname(os.path.abspath(__file__))

    table = pd.read_csv(f"{script_dir}\\downloads\\mres_drosophila.microT.{db}.txt.gz", sep="\t")
    df = pd.DataFrame(table)
    with open(fbid_list, 'r') as fbid:
        with open(mirna_list, 'r') as mirna:
            fbids, mirnas = [x.strip() for x in fbid.readlines()], [x.strip() for x in mirna.readlines()]
            df_list = df.values.tolist()
            results = []

            for row in tqdm(df_list, desc="Looking for matches", unit=" lines", ncols=100):

                modif_row = [str(item).lower() for item in row]

                for m in mirnas:
                    # print(f'Finding targets for {m}')
                    modif_m = m[:m.rfind('-')] + '_' + m[m.rfind('-') + 1:] if '-' in m else m

                    for f in fbids:
                        # print(f'Checking {f}')
                        if (m.lower().strip("*") in str(modif_row)) and (f.lower() in str(modif_row)):
                            results.append(row)
                        if (modif_m.lower().strip("*") in str(modif_row)) and (f.lower() in str(modif_row)):
                            results.append(row)

    with open(f'{script_dir}\\downloads\\mres_drosophila.microT.{db}.txt', 'r') as copy:
        header = copy.readline()

    if os.path.exists(f'{script_dir}\\matches_{db}.csv'):
        os.remove(f'{script_dir}\\matches_{db}.csv')

    if output:
        script_dir = output

    with open(f'{script_dir}\\matches_{db}.csv', 'a') as matches:
        matches.write(header.replace("\t", ","))

        unique_results = set()

        for result in results:
            formatted_result = str(result)[1:-1].replace(", ", ',').replace("'", "")
            unique_results.add(formatted_result)  # Добавляем строку в множество и избавляемся от дубликатов

        for unique_result in unique_results:
            matches.write(unique_result + "\n")

    final_df = pd.read_csv(f'{script_dir}\\matches_{db}.csv')

    if final_df.empty:
        print('No matches found.')
    else:
        print(f"Results saved in {script_dir}")
        print(final_df.to_string())


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-fb", "--FBID_list", type=str,
                        help="Path to the potential targets list, i.e. the list of FlyBase gene IDs (FBgnxxxxxxx) "
                             "corresponding to genes which encode proteins user would like to check "
                             "AND/OR the list of FlyBase transcript IDs (FBtrxxxxxxx).",
                        required=True),
    parser.add_argument("-mi", "--miRNA_list", type=str,
                        help="Path to the list of miRNAs user is looking to find targets for",
                        required=True),
    parser.add_argument("-db", "--database", type=str,
                        help="Name of the database to find miRNA+FBID matches in. "
                             "Type 'mirbase' for MiRBase, 'mirgenedb' for MirGeneDB",
                        required=True)
    parser.add_argument("-o", "--output", type=str,
                        help="Saving directory. Folder created if doesn't exist. Script directory by default")

    args = parser.parse_args()

download_db(args.database)
find_matches(args.database, args.FBID_list, args.miRNA_list, args.output)
