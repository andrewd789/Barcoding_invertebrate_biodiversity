# -*- coding: utf-8 -*-
"""
Created on Tue Aug 21 16:01:34 2018

@author: dopheidea
"""
import os, glob, re
from Bio import SeqIO

# Exclude any Hauturu barcode sequences from downloaded BOLD database sequences
# Hauturu barcodes have Genbank IDS in the range KP420745 - KP422464
# http://gamon.webfactional.com/regexnumericrangegenerator/
reg = "KP4(2074[5-9]|207[5-9][0-9]|20[89][0-9]{2}|21[0-9]{3}|22[0-3][0-9]{2}|224[0-5][0-9]|2246[0-4])$"

os.chdir('G:/Documents/PhD/Sanger_OTUs_analysis/Sanger_BOLD_454_matching/BOLD_NZ_seqs_2018/')
files = glob.glob("*.fas")
for f in files:
    keep = list()
    exclude = list()
    with open(f, "r") as infile:
        label = f.split(".fas")[0]
        with open("{0}_keep.fasta".format(label), "a") as outfile:
            for seq in SeqIO.parse(infile, "fasta"):
                r = re.search(reg, seq.id)
                if r is not None:
                    exclude.append(seq.id)
                    print("excluded:{0}".format(seq.id))
                elif r is None:
                    keep.append(seq.id)
                    SeqIO.write(seq, outfile, "fasta")