import argparse
import sys
import os
import pandas as pd
from Bio import Entrez
Entrez.email = "zhcong.pku@gmail.com"
import time


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input", type=str, required=True,
                        help="A file with accessions to download")
    parser.add_argument("-o", "--output_dir", type=str, required=True,
                        help="The directory to write downloaded files to")
    args = parser.parse_args()
    return vars(args)

def write_seq(in_file, dir):
    dat = pd.read_csv(in_file, sep = "\t", header = 0)
    for index, row in dat.iterrows():
        strain = row['Strain']
        Replicons = row["Replicons"]
        reps = Replicons.split(";")
        seq={}
        print(strain+"\n")
        print(dir)
        with open(os.path.join(dir,  strain+".fasta"), "w") as output:
            for x in reps:
                if (x.split(":")[0].find("plasmid") > -1  ):
                    handle = Entrez.efetch(db="nucleotide", id=x.split(":")[1].split("/")[0], rettype="fasta", retmode="text")
                    record = handle.read()
                    output.write(record.rstrip('\n'))
        time.sleep(2)


def main():
    c_args = parse_args()
    in_file = os.path.abspath(c_args["input"])
    op_dir = os.path.abspath(c_args["output_dir"])
    if not os.path.exists(op_dir):
        os.makedirs(op_dir)
    dat = pd.read_csv(in_file, sep = "\t", header = 0)
    print("#######")
    print(in_file)
    write_seq(in_file, op_dir)
    


if __name__ == "__main__":
    main()

