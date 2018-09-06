# -*- coding: utf-8 -*-

import glob, os, re, subprocess

#usearch = "usearch7.0.1090_i86linux64"
#usearch = "G:/Documents/Bioinformatics/Usearch/usearch8.0.1623_win32.exe"
usearch = "G:/Documents/Bioinformatics/Vsearch/vsearch-2.4.4-win-x86_64/vsearch.exe" # Vsearch gives much the same results
startupinfo = subprocess.STARTUPINFO() # Prevents cmd windows from opening
startupinfo.dwFlags |= subprocess.STARTF_USESHOWWINDOW # Prevents cmd windows from opening

# OTU clustering using UCLUST and recentering approach
# from https://www.drive5.com/usearch/manual7/recenter.html
def get_usearch_Sanger_OTUs(label, f, log):
    print("Clustering otus from {0} sequences...".format(label))
    # Dereplicate the sequences
    derep = "{0} -derep_fulllength {1} -output {2}_uniques.fasta -sizeout".format(usearch, f, label)
    subprocess.call(derep.split(), stdout=log, stderr=subprocess.STDOUT, startupinfo=startupinfo)
    log.write('\nDereplicate filtered/trimmed reads:\n' + str(derep) + '\n')
    # Cluster the dereplicated sequences into OTUs (automatically sorted by length)
    clusterOTUs = "{0} -cluster_fast {1}_uniques.fasta -id 0.97 -centroids {1}_OTUs_centroids.fasta \
                  -consout {1}_OTUs_consensus.fasta".format(usearch, label)
    subprocess.call(clusterOTUs.split(), stdout=log, stderr=subprocess.STDOUT, startupinfo=startupinfo)
    log.write('\ncluster_fast into OTUs at 97%:\n' + str(clusterOTUs) + '\n')
    # Sort the consensus sequences by abundance
    sortbysize = "{0} -sortbysize {1}_OTUs_consensus.fasta -output {1}_OTUs_consensus_sorted.fasta".format(usearch, label)
    subprocess.call(sortbysize.split(), stdout=log, stderr=subprocess.STDOUT, startupinfo=startupinfo)
    log.write('\nSort OTU consensus seqs by size:\n' + str(sortbysize) + '\n')
    # Cluster the abundance-sorted consensus sequences again, to remove any redundant OTUs
    clusterOTUs2 = "{0} -cluster_smallmem {1}_OTUs_consensus_sorted.fasta -usersort -id 0.97 \
                   -centroids {1}_OTUs_centroids_2.fasta".format(usearch, label)
    subprocess.call(clusterOTUs2.split(), stdout=log, stderr=subprocess.STDOUT, startupinfo=startupinfo)
    log.write('\nCluster consensus seqs into OTUs at 97%:\n' + str(clusterOTUs2) + '\n')
    # Assign all the barcode sequences to the OTU centroids
    searchReads = "{0} -usearch_global {1} -db {2}_OTUs_centroids_2.fasta -strand plus \
                  -id 0.97 -uc {2}_OTUs_readmap.uc -output_no_hits".format(usearch, f, label)
    subprocess.call(searchReads.split(), stdout=log, stderr=subprocess.STDOUT)
    log.write("\nMap the filtered/trimmed reads to the OTUs\n" + searchReads + "\n")

# Function to make an OTU table from usearch readmap output
def make_otutable_usearch(readmap):
    samples = []
    otus = []
    hits = {}
    # Get samples and otu centroids from readmap
    print("Getting samples and otus from readmap...")
    for row in readmap:
        if row.startswith("H"):
            query = row.strip().split("\t")[8]
            sample = "{0}_{1}".format(query.split("|")[2], query.split("|")[4])
            #sample = query.split("|")[2]
            if sample not in samples:
                samples.append(sample)
            hit = row.strip().split("\t")[9]
            #hit = re.split("^centroid=", hit)[1].strip()
            if hit not in otus:
                otus.append(hit)
            hits[query] = hit
        # Include any sequences without hits to OTU centroids
        elif row.startswith("S"):
            seq = row.split("\t")[8].strip()
            #print(seq)
            if seq not in otus:
                otus.append(seq)
            hits[seq] = "S"                
    
    # Set up bins and tables for otu counts
    bins = [[0 for row in range(len(samples))] for col in range(len(otus))]
    # Get OTU counts from readmap
    print("Getting OTU counts...")
    for query, hit in hits.items():
        #sample = "{0}_{1}".format(query.split("|")[2], query.split("|")[4])
        sample = query.split("|")[2]
        if hit is not "S":
            bins[otus.index(hit)][samples.index(sample)] += 1
        else:
            bins[otus.index(query)][samples.index(sample)] += 1

    # Write counts to table        
    with open("{0}_otutable_U6.txt".format(label), "w") as otutable:
        for sample in samples:
            # Make the two highest elevation sample names consistent with the others
            sample = re.sub("CM30c30", "9", sample)
            sample = re.sub("LB1", "10", sample)
            otutable.write("\t{0}".format(sample))
        for i in range(len(bins)):
            otutable.write("\n{0}".format(otus[i]))
            for item in bins[i]:
                otutable.write("\t{0}".format(item))
    print("Finished {0}".format(label))

# Vsearch readmap output is a bit different, includes linebreaks, requiring different treatment:
def make_otutable_vsearch(readmap):
    samples = []
    otus = []
    hits = {}
    # Get samples and otu centroids from readmap
    print("Getting samples and otus from readmap...")
    for row in readmap:
        if row.startswith("H"):
            query = row.strip().split("\t")[8]
            sample = "{0}_{1}".format(query.split("|")[2], query.split("|")[4])
            #sample = query.split("|")[2]
            if sample not in samples:
                samples.append(sample)
            hit = readmap.readline().strip() # Load the next line
            #otu = readmap.readline().split("\t")[1]
            hit = re.split("^centroid=", hit)[1].strip()
            if hit not in otus:
                otus.append(hit)
            hits[query] = hit
        # Include any sequences without hits to OTU centroids
        elif row.startswith("N"):
            seq = row.split("\t")[8].strip()
            #print(seq)
            if seq not in otus:
                otus.append(seq)
            hits[seq] = "None"                
    
    # Set up bins and tables for otu counts
    bins = [[0 for row in range(len(samples))] for col in range(len(otus))]
    # Get OTU counts from readmap
    print("Getting OTU counts...")
    for query, hit in hits.items():
        sample = "{0}_{1}".format(query.split("|")[2], query.split("|")[4])
        #sample = query.split("|")[2]
        if hit is not "None":
            bins[otus.index(hit)][samples.index(sample)] += 1
        else:
            bins[otus.index(query)][samples.index(sample)] += 1

    # Write counts to table        
    with open("{0}_otutable.txt".format(label), "w") as otutable:
        for sample in samples:
            # Make the two highest elevation sample names consistent with the others
            sample = re.sub("CM30c30", "9", sample)
            sample = re.sub("LB1", "10", sample)
            otutable.write("\t{0}".format(sample))
        for i in range(len(bins)):
            otutable.write("\n{0}".format(otus[i]))
            for item in bins[i]:
                otutable.write("\t{0}".format(item))



# To cluster COI and 28S barcode sequences into OTUs:
os.chdir("G:/Documents/GitHub/Barcoding_invertebrate_biodiversity/Invert_DNA_barcode_data/")
os.chdir("./Invert_DNA_barcode_data/")
infiles = glob.glob("*barcode_sequences.fasta")
for f in infiles:
    label = f.split("_barcode_")[0]
    print(label)
    with open("{0}_usearch_log.txt".format(label), "a") as log:
        get_usearch_Sanger_OTUs(label, f, log)
    with open("{0}_OTUs_readmap.uc".format(label), "r") as readmap:
        #make_otutable_usearch(readmap)
        make_otutable_vsearch(readmap)
    print("Finished {0}".format(label))

# To cluster BOLD database sequences into OTUs:
os.chdir("./BOLD_NZ_seqs_2018/")
infiles = glob.glob("BOLD_NZ*keep.fasta")
for f in infiles:
    label = f.split(".fa")[0]
    with open("{0}_usearch_log.txt".format(label), "a") as log:
        get_usearch_Sanger_OTUs(label, f, log)
    print("Finished {0}".format(label))