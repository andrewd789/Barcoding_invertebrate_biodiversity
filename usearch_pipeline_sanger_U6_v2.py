# -*- coding: utf-8 -*-
"""
Created on Tue Sep 09 13:02:56 2014

@author: adop001

Usearch pipeline with re-centering strategy
(R. Edgar's recommended approach prior to UPARSE)

"""
import glob
import os
import subprocess

#usearchpath = "../Usearch7/usearch7.0.1090_i86linux64"
#usearchpath = 'G:/Documents/Bioinformatics/Usearch/usearch9.0.2132_win32.exe'
usearchpath = 'G:/Documents/Bioinformatics/Vsearch/vsearch-2.4.4-win-x86_64/vsearch.exe'
startupinfo = subprocess.STARTUPINFO() # Prevents cmd windows from opening
startupinfo.dwFlags |= subprocess.STARTF_USESHOWWINDOW # Prevents cmd windows from opening

def get_usearch_Sanger_OTUs(label, input_file, log):

    derep = "%s -derep_fulllength %s \
        -output %s_uniques.fasta -sizeout" % (usearchpath, input_file, label)
    subprocess.call(derep.split(), stdout=log, stderr=subprocess.STDOUT, startupinfo=startupinfo)
    log.write('\nDereplicate filtered/trimmed reads:\n' + str(derep) + '\n')

    clusterOTUs = ("%s -cluster_fast %s_uniques.fasta -id 0.97 -centroids %s_OTUs_centroids.fasta \
                   -consout %s_OTUs_consensus.fasta -sizeout -uc %s_OTUs_readmap.uc"
                   % (usearchpath, label, label, label, label))
    subprocess.call(clusterOTUs.split(), stdout=log, stderr=subprocess.STDOUT, startupinfo=startupinfo)
    log.write('\ncluster_fast into OTUs at 97%:\n' + str(clusterOTUs) + '\n')

    sortbysize = ("%s -sortbysize %s_OTUs_consensus.fasta \
                  -output %s_OTUs_consensus_sorted.fasta"
                  % (usearchpath, label, label))
    subprocess.call(sortbysize.split(), stdout=log, stderr=subprocess.STDOUT, startupinfo=startupinfo)
    log.write('\nSort OTU consensus seqs by size:\n' + str(sortbysize) + '\n')

    clusterOTUs2 = ("%s -cluster_smallmem %s_OTUs_consensus_sorted.fasta \
                    -id 0.97 -usersort -centroids %s_OTUs_centroids_2.fasta \
                    -uc %s_OTUs_readmap_2.uc"
                    % (usearchpath, label, label, label))
    subprocess.call(clusterOTUs2.split(), stdout=log, stderr=subprocess.STDOUT, startupinfo=startupinfo)
    log.write('\nCluster consensus seqs into OTUs at 97%:\n' +
              str(clusterOTUs2) + '\n')

#    OTUtable = open(("%s_OTUtable.txt" % label), "w")
#    makeOTUtable = ("python ../Usearch7/uc2otutab_andrew.py \
#     %s_OTU_readmap.uc > %s_OTUtable.txt" % (label, label))
#    subprocess.call(makeOTUtable.split(), stdout=OTUtable,
#                    stderr=subprocess.STDOUT)
#    OTUtable.close()

#listing = os.listdir("./")
#for infile in listing:
#    if (str("BOLD_NZ") in str(infile) and str(infile).endswith(".fas")):
#        input_file = infile
#        label = str.split(infile, ".fas")[0]
#        log = open(("%s_usearch_log.txt" % label), "a+")
#        get_usearch_Sanger_OTUs_U6(label, input_file, log)

#os.chdir('G:/Documents/BOLD_arthropod_seqs_data/')
#label = 'BOLD_arthropod_seqs'
#input_file = 'BOLD_24180_arthropod_seqs.fas'
os.chdir('G:/Documents/PhD/Sanger_OTUs_analysis/Sanger_BOLD_454_matching/BOLD_NZ_seqs_2018/')
infiles = glob.glob("BOLD_NZ*keep.fasta")
for f in infiles:
    label = f.split(".fas")[0]
    with open(('%s_usearch_log.txt' % label), 'a+') as log:
        get_usearch_Sanger_OTUs(label, f, log)
