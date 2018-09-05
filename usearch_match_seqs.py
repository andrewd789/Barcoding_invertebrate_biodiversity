# -*- coding: utf-8 -*-
"""
Created on Fri May 16 16:17:04 2014

@author: adop001

Matches COI OTUs between various data sets
"""
import glob
import math
import os
import subprocess
from Bio import SeqIO

#usearchpath = "G:/Documents/Bioinformatics/Usearch/usearch7.0.1090_i86linux64"
#usearchpath = "G:/Documents/Bioinformatics/Usearch/usearch9.0.2132_win32.exe"
usearchpath = "G:/Documents/Bioinformatics/Vsearch/vsearch-2.4.4-win-x86_64/vsearch.exe"
startupinfo = subprocess.STARTUPINFO() # Prevents cmd windows from opening
startupinfo.dwFlags |= subprocess.STARTF_USESHOWWINDOW # Prevents cmd windows from opening

def match_seqs (f1, f2, label, tx, log):
    searchReads = "{0} -usearch_global {1} \
        -db {2} -strand plus -id 0.97 -uc {3}_match-97_readmap.uc".format(usearchpath, f1, f2, label)
    subprocess.call(searchReads.split(), stdout = log, stderr = subprocess.STDOUT)
    log.write(("Match COI seqs against {0} at 97 percent\n".format(label)) + str(searchReads) + '\n')
#    searchReads2 = "%s -usearch_global CO1_366_OTU_centroids.fasta \
#        -db %s -strand plus -id 0.90 -uc %s_match-90_readmap.uc" % (usearchpath, input_file, label)
#    subprocess.call(searchReads2.split(), stdout = log, stderr = subprocess.STDOUT)
#    log.write(('\nMatch CO1 seqs against %s at 90 percent\n' % label) + str(searchReads) + '\n')    

os.chdir("G:/Documents/PhD/Sanger_OTUs_analysis/Sanger_BOLD_454_matching/")
files1 = glob.glob("Sanger_OTUs_by_taxa/*.fasta")
#files2 = glob.glob("COI-454_OTUs_by_taxa/COI-454_*.fasta") # for matching invertebrate barcode OTUs with soil metabarcoding COI OTUs
files2 = glob.glob("BOLD_NZ_seqs_2018/*centroids_2.fasta") # for matching invertebrate barcode OTUs with BOLD database COI OTUs

nOTUs = dict()

#with open("COI-Sanger-454-OTUs_matching_log.txt", "a") as log:
with open("COI-Sanger-BOLD-OTUs_matching_log.txt", "a") as log:
    for f1 in files1:
        tx = f1.strip("Sanger_OTUs_by_taxa").split("_")[1]
        #f2 = glob.glob("COI-454_OTUs_by_taxa/COI-454_{0}*.fasta".format(tx))
        f2 = glob.glob("BOLD_NZ_seqs_2018/*{0}*centroids_2.fasta".format(tx))
        if len(f2) == 1: # Matching dataset file found
            f2 = f2[0]
            print(tx)
            # Get total number of OTUs per file
            s1 = [seq.id for seq in SeqIO.parse(f1, "fasta")]
            s2 = [seq.id for seq in SeqIO.parse(f2, "fasta")]
            nOTUs[tx] = len(s1), len(s2)
            #label = "{0}_Sanger_OTUs_vs_454".format(tx)
            label = "{0}_Sanger_OTUs_vs_BOLD_NZ".format(tx)
            match_seqs(f1, f2, label, tx, log)


#with open("COI-Sanger-454-OTUs_matching_summary.txt", "a") as summary:
with open("COI-Sanger-BOLD-OTUs_matching_summary.txt", "a") as summary:
    summary.write("group\tbarcode OTUs\tx OTUs\tmatching\testimate\n")
    readmaps = glob.glob("*readmap.uc")
    for f in readmaps:
        tx = f.split("_")[0]
        #print(tx)
        K = nOTUs.get(tx)[0] # Number of OTUs in first dataset
        n = nOTUs.get(tx)[1] # Number of OTUs in second dataset
        matched = set()
        #matched = dict()
        unmatched = set()
        for row in open(f):
            #print(row)
            fields = row.strip().split("\t")
            if fields[0] is "H":
                matched.add(fields[8])
                #matched[fields[8]] = fields[9]
            elif fields[0] is "N":
                unmatched.add(fields[8])
        k = len(matched) # Number of matching OTUs
        u = len(unmatched) # Number of non-matching OTUs
        N = (((K+1)*(n+1))/(k+1))-1 # Chapman estimator
        #N = ((K-1)*(n-1)/(k-2)) # Bayesian estimator
        print("{0}: {1} barcode OTUs; {2} other OTUs; {3} matched; {4:.0f} estimated species".format(tx, K, n, k, N))
        summary.write("{0}\t{1}\t{2}\t{3}\t{4:.0f}\n".format(tx, K, n, k, N))
