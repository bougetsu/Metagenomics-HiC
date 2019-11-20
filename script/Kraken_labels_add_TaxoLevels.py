#!/usr/bin/env python

import sys,os,re
import numpy as np

import pandas as pd

kraken_report = "/ms/11/cong/project/HiC/ASSEMBLY/MEGAHIT/CoAsm_150PE/Taxify/Kraken/final.contigs.krareport_NoIncident"
krak = pd.read_csv(kraken_report, sep = "\t", header = 0, skiprows= 3)
print(krak.loc[[1]])


#level = pd.Series(krak.rank,index=krak.taxName).to_dict()
level = dict(zip(krak["taxName"],krak["rank"]))

print("###############")

print(level["Firmicutes"])
print("###############")

krak_label = "/ms/11/cong/project/HiC/ASSEMBLY/MEGAHIT/CoAsm_150PE/Taxify/Kraken/final.contigs.kraLabel"

f = open(krak_label+"_AddRank", "w")


for line in open(krak_label, "r"):
		#print(line)
		cut=line.strip().split("\t")
		if(len(cut)>1):
			tax=cut[1].split(";")
                        rank_tax = ""
			for i in range(len(tax)):
				#print(t)
				if tax[i] in level and (level[tax[i]] != "no rank"):
					#print(tax[i])
					tax[i] = level[tax[i]].capitalize() + "__" + tax[i]
                                        rank_tax = rank_tax+tax[i]+";"
                        if(rank_tax != ""):
                            rank_tax = rank_tax[:-1]
			    f.write(cut[0]+"\t"+rank_tax+"\n")
f.close()

