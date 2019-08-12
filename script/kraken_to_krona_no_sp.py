#!/usr/bin/env python2.7

import argparse

from natsort import humansorted

from collections import defaultdict

import numpy as np

import pandas as pd



def parse_args():

    parser = argparse.ArgumentParser(description='Parse BAM component filter arguments.')

    parser.add_argument('--kraken_file', '-k', type=str, required=True,

                        help='translated kraken file (get from metawrap)')

    parser.add_argument('--cov_file', '-c', type=str, required=True,

                        help='coverage file generate by metabat')

    parser.add_argument('--out_file', '-o', required=True, type=str, default='output',

                        help='krona file')

    args = parser.parse_args()

    return vars(args)
# This script takes in a translated kraken file of either contigs (from SPAdes) or reads, and parses it into a format for ktImportText to produce a kronachart.

def merge_kv(kraken_file, cov_file, outfile):

	krak = pd.read_csv(kraken_file, sep = "\t", names=["contig", "tax"])
	#krak['tax'] = krak['tax'].str.replace(';', '\t')
	cov = pd.read_csv(cov_file, sep = "\t", names=["contig", "lenth", "AvgDepth", "bam", "bam-ver"])
	merged = pd.merge(krak, cov, on = 'contig', how = 'right')
	merged.to_csv("merged_kraken_cov.csv", sep = "\t", index = False)

def k2k(kraken_file, cov_file, outfile):
	data={}
	lend={}
	for line in open("merged_kraken_cov.csv"):
		cut=line.strip().split("\t")
		print(cut[0])
		if cut[0] == "contig" or cut[0] == "contigName":
			continue
		name=cut[0]
		tax="\t".join(cut[1].split(";"))
		lens=cut[2]
		dep=cut[3]

		if tax in data: 
			data[tax]+=float(lens) * float(dep)
			lend[tax]+=float(lens)
		else:
			data[tax]=float(lens) * float(dep)
			lend[tax]=float(lens)
	with open(outfile, 'w') as out:
		for tax in data:
			out.write(str(data[tax]/lend[tax]) +"\t" +tax + "\n")

def main():

    c_args = parse_args()

    kraken_file = c_args["kraken_file"]
    cov_file = c_args["cov_file"]
    outfile=c_args["out_file"]+".krona"

    merge_kv(kraken_file, cov_file, outfile)
    k2k(kraken_file, cov_file, outfile)




if __name__ == "__main__":

    main()

