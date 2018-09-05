# -*- coding: utf-8 -*-

import os
import re
from Bio import SeqIO

filepath = ("H:/My Documents/PhD Research PFR folder/")


#os.chdir("H:/My Documents/PhD Research PFR folder/COI-Sanger-454_matching/")

OTUs = "CO1_366_OTU_centroids.fasta"
taxatable = "CO1_366_OTU_centroid_taxatable.txt"
taxdict = {}

OTUs_dict = SeqIO.to_dict(SeqIO.parse(OTUs, "fasta"))

taxa = ["Arthropoda", "Arachnida", "Malacostraca", "Myriapoda", "Insecta", \
"Coleoptera", "Diptera", "Lepidoptera", "Hymenoptera", "Hemiptera", "Orthoptera", \
"Mollusca", "Annelida", "Onychophora"]

for t in taxa:
    print(t)
    outfile = open("COI_%s_OTU_centroids.fasta" % t, "w")
    for row in open(taxatable, "r"):
        if t in row:
            print("%s found" % t)
            OTU = row.split("\t")[0]
            print(OTU)
            rec = OTUs_dict.get(OTU)
            #print(rec)
            if rec is not None:
                SeqIO.write(rec, outfile, "fasta")
    outfile.close()

#################################
outfile = open("COI_Malacostraca_OTU_centroids.fasta", "w")
for record in SeqIO.parse(infile, "fasta"):
    if "Amphipoda" in record.id or "Isopoda" in record.id or "Malacostraca" in record.id:
        print(record)
        SeqIO.write(record, outfile, "fasta")
outfile.close()

outfile = open("COI_Gastropoda_OTU_centroids.fasta", "w")
for record in SeqIO.parse(infile, "fasta"):
    if "Eupulmonata" in record.id:
        print(record)
        SeqIO.write(record, outfile, "fasta")
outfile.close()
                
outfile = open("COI_Annelida_OTU_centroids.fasta", "w")
for record in SeqIO.parse(infile, "fasta"):
    if "Haplotaxida" in record.id:
        print(record)
        SeqIO.write(record, outfile, "fasta")
outfile.close()

outfile = open("COI_Myriapoda_OTU_centroids.fasta", "w")
for record in SeqIO.parse(infile, "fasta"):
    if "Myriapoda" in record.id or "Pauropoda" in record.id or"Symphyla" in record.id \
    or"Chilopoda" in record.id or"Geophilomorpha" in record.id or"Lithobiomorpha" in record.id \
    or "Diplopoda" in record.id or "Polydesmida" in record.id or "Spirostreptida" in record.id \
    or "Glomerida" in record.id:
        print(record)
        SeqIO.write(record, outfile, "fasta")
outfile.close()
                                     
outfile = open("COI_Arachnida_OTU_centroids.fasta", "w")
for record in SeqIO.parse(infile, "fasta"):
    if "Araneae" in record.id or "Opiliones" in record.id or "Pseudoscorpiones" in record.id \
    or "Arachnida" in record.id:
        print(record)
        SeqIO.write(record, outfile, "fasta")
outfile.close()

#####################################
