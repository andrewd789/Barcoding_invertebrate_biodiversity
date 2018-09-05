# Analysis of Hauturu DNA barcoded invertebrates OTU tables

library(ape)
library(ggplot2)
library(gridExtra)
library(phyloseq)
library(reshape2)
#library(RColorBrewer)
library(scales)
library(vegan)
#library(vegetarian)

theme_set(theme_bw(base_size = 10))

###############################################################################
# Load all the data 
#setwd("G:/Documents/PhD/Sanger_vs_454_COI_stuff/Sanger_OTUs_analysis")
setwd("G:/Documents/GitHub/Barcoding_invertebrate_biodiversity/")

# Invertebrate DNA barcoding summary data:
# Successful barcoding results:
barcodes <- read.table("Invert_DNA_barcode_data/28S+COI_sequence_taxonomy.txt", 
                       sep = "\t", header = TRUE, check.names = FALSE)
# Unsuccesful barcoding attempts
fails <- read.table("Invert_DNA_barcode_data/Barcode_failures.txt",
                    sep = "\t", header = TRUE, check.names = FALSE)

# OTU tables:
OTUtable.COI <- read.table("Invert_DNA_barcode_data/COI_366_OTUtable.txt", 
                           header = TRUE, row.names = 1, sep = "\t", check.names=FALSE)
OTUtable.COI <- t(OTUtable.COI)
OTUtable.COI <- OTUtable.COI[order(rownames(OTUtable.COI)),]

OTUtable.28S <- read.table("Invert_DNA_barcode_data/28S_247_OTUtable.txt", 
                           header = TRUE, row.names = 1, sep = "\t", check.names=FALSE)
OTUtable.28S <- t(OTUtable.28S)
OTUtable.28S <- OTUtable.28S[order(rownames(OTUtable.28S)),]

# OTU taxonomy data:
taxatable.COI <- read.table("Invert_DNA_barcode_data/COI_366_OTU_centroid_taxatable.txt", 
                            sep="\t", header=T, row.names=1)
taxatable.28S <- read.table("Invert_DNA_barcode_data/28S_247_OTU_centroid_taxatable.txt", 
                            sep="\t", header=T, row.names=1)

Orders <- c("Protostomia","Arthropoda","Arachnida","Araneae","Opiliones","Pseudoscorpiones",
            "Myriapoda","Pauropoda","Symphyla","Geophilomorpha","Lithobiomorpha","Diplopoda",
            "Polydesmida","Spirostreptida",
            "Amphipoda","Isopoda","Diplura","Collembola","Archaeognatha",
            "Endopterygota","Neoptera","Lepidoptera","Coleoptera","Diptera","Hymenoptera","Neuroptera",
            "Siphonaptera","Blattodea","Orthoptera","Phasmatodea","Hemiptera","Psocoptera","Thysanoptera",
            "Haplotaxida","Gastropoda","Onychophora") 

taxatable.COI$Order <- factor(taxatable.COI$Order, levels = Orders, ordered = TRUE)
taxatable.28S$Order <- factor(taxatable.28S$Order, levels = Orders, ordered = TRUE)

# Phylo trees:
tree.COI <- read.tree("Invert_DNA_barcode_data/RAxML_bestTree.COI_good_LCA_366_OTUs_tree_GTRGAMMA")
tree.28S <- read.tree("Invert_DNA_barcode_data/RAxML_bestTree.28S_good_LCA_247_OTUs_tree_GTRGAMMA")
#str(tree.COI)
#str(tree.28S)

# Sample metadata:
samples <- read.table("Invert_DNA_barcode_data/LBI_Pilot_chem_elevation_data.txt", 
                      sep="\t", header=T, row.names=1)
samples <- samples[order(rownames(samples)),]
samples$PlotN <- gsub("Plot-[A-Z]", "", samples$Plot)
samples$PlotN <- factor(samples$PlotN, ordered = TRUE,
                        levels = c("01","02","03","04","05","06","07","08","09","10"))
samples$Plot_Method <- paste(samples$Method, samples$Plot)
samples$Plot_Method <- factor(samples$Plot_Method, ordered = TRUE)

# OTU top BLAST matches:
bl <- read.table("Invert_DNA_barcode_data/28S+COI_OTUs_top_blast_matches_organised.txt", 
                 header = TRUE, row.names = 1, sep = "\t", check.names = FALSE)

# 454 soil DNA OTU table and taxonomy for biodiversity estimates
OTUtable.454 <- read.table("Soil_DNA_454_data/COI_454_inverts_OTUtable_byPlot.txt", 
                           sep = "\t", header = TRUE, row.names = 1)
#OTUtable.454 <- t(OTUtable.454)
taxonomy.454 <- read.table("Soil_DNA_454_data/COI_454_inverts_taxonomy.txt", 
                     sep = "\t", header = TRUE, row.names = 1)

###############################################################################
# Make a summary table of DNA barcoding results (a bit complicated!)
d <- dcast(barcodes, Phylum + Class + Order ~ Gene)
d$x28S <- d$`28S and COI`+ d$`28S only`
d$COI <- d$`28S and COI`+ d$`COI only`
d$`28S only` <- NULL
d$`COI only` <- NULL
d$total <- d$x28S + d$COI - d$`28S and COI`
colnames(d) <- gsub("28S and COI", "Both_loci", colnames(d))

# Make an extra row with totals for each column:
totals <- data.frame(matrix(nrow = 1, ncol = ncol(d)))
colnames(totals) <- colnames(d)
totals[1, ] <- c("total","total","total", colSums(d[, 4:7]))

# Add failed barcoding attempt details:
fails$`Combined order` <- as.character(fails$`Combined order`)
fails$`Combined order` <- gsub("Empty vial|0", "unknown", fails$`Combined order`)
f <- dcast(fails, `Combined order` ~ Gene)
f$x28S <- f$`28S` + f$Both
f$`28S` <- NULL
f$COI <- f$COI + f$Both
colnames(f) <- gsub("Both", "Both_loci", colnames(f))
x <- colSums(f[,2:4])
totals[, 4:6] <- paste0(totals[, 4:6], " (", x[match(colnames(totals[, 4:6]), names(x))], ")")

# Several groups only had failures, so need to be added separately to barcoding results summary:
x <- f$`Combined order`[!(f$`Combined order` %in% d$Order)] # "(Arachnida)" "Glomerida"   "unknown"
y <- data.frame(matrix(nrow = length(x), ncol = ncol(d)))
colnames(y) <- colnames(d)
y$Order <- x
y$Phylum <- c("Arthropoda","Arthropoda","Unknown")
y$Class <- c("Arachnida","Diplopoda","Unknown")
d <- rbind(d, y) # Add to the summary table

# Add details of unsuccessful barcoding attempts to the summary:
d$x28S <- paste0(d$x28S, " (", f$x28S[match(d$Order, f$`Combined order`)], ")")
d$COI <- paste0(d$COI, " (", f$COI[match(d$Order, f$`Combined order`)], ")")
d$Both_loci <- paste0(d$Both_loci, " (", f$Both[match(d$Order, f$`Combined order`)], ")")
d$x28S <- gsub("^0 |^NA |\\([0NA]*\\)", "", d$x28S)
d$COI <- gsub("^0 |^NA |\\([0NA]*\\)", "", d$COI)
d$Both_loci <- gsub("^0 |^NA |\\([0NA]*\\)", "", d$Both_loci)

# Add collection method details to the summary:
m <- dcast(barcodes, Phylum + Class + Order ~ Method) 
d$Leaf_litter <- m$`Leaf-Litter`[match(d$Order, m$Order)]
d$Pitfall_trap <- m$`Pitfall-Trap`[match(d$Order, m$Order)]

# Get OTU counts for each Order:
OTUs_28S <- aggregate(taxatable.28S$Order, by = list(taxatable.28S$Order), length)
OTUs_COI <- aggregate(taxatable.COI$Order, by = list(taxatable.COI$Order), length)

# Make the OTU Orders match the summary table. (Put non-Orders in parentheses):
not.Orders <- paste0("(Protostomia|Arthropoda|Arachnida|Myriapoda|Chilopoda|Diplopoda|",
                     "Malacostraca|Hexapoda|Diplura|Insecta|Endopterygota|Neoptera|Gastropoda|Onychophora)")
OTUs_28S$Group.1 <- gsub(not.Orders, "\\(\\1\\)", OTUs_28S$Group.1)
OTUs_COI$Group.1 <- gsub(not.Orders, "\\(\\1\\)", OTUs_COI$Group.1)

# Some barcodes were identified as Endopterygota (Infraclass), whereas corresponding OTUs were 
# identified as Neoptera (Subclass; one semi-rank lower, hardly more informative). 
# So for consistency, change the Neoptera ids to Endopterygota: 
OTUs_COI$Group.1[!OTUs_COI$Group.1 %in% d$Order]
OTUs_COI$Group.1 <- gsub("Neoptera", "Endopterygota", OTUs_COI$Group.1)

# Now add the OTU counts to the summary table:
d$x28S_OTUs <- OTUs_28S$x[match(d$Order, OTUs_28S$Group.1)]
d$COI_OTUs <- OTUs_COI$x[match(d$Order, OTUs_COI$Group.1)]

# Get totals for methods and OTUs, then add to summary
d[is.na(d)] <- 0 # Convert NAs to zeros
x <- as.data.frame(colSums(d[, 8:11]))
totals <- cbind(totals, t(x)) # Join all the totals together
d <- rbind(d, totals) # Join totals row to the summary table
d <- d[, c(1:3, 5:6, 4, 7:11)] # Re-order the columns
d <- d[order(d$Phylum, d$Class, d$Order), ] # Re-order the rows
d[d == 0] <- "" # Replace zeros with blank

# Output to file (The final table will need further re-organization by taxonomy):
write.table(d, file = "Invert_barcoding_summary_table.txt", sep = "\t", row.names = FALSE, quote = FALSE)

###############################################################################
# Barplot of top BLAST matches to OTU sequences

cols <- c("Gastropoda" = "#FFAA55",          
          "Myriapoda" = "#7F55FF",
          "Diplopoda"	= "#7F55FF",
          "Chilopoda"  = "#7F55FF",
          "Arachnida" = "#65BA65", 
          "Malacostraca" = "#AA7FFF", 
          "Insecta"	= "#7FAAFF",
          "Hexapoda" = "#0055FF",
          "Diplura" = "#0055FF",
          "Ellipura" = "#0055FF",
          "Arthropoda" = "#C3C3C3",
          "Oligochaeta" = "#FF557F",
          "Onychophora" = "#2D82AC",
          "Protostomia" = "#5A5A5A")

p1 <- ggplot(bl, aes(x = `Identical matches %`, fill = Taxon)) + 
  geom_histogram(binwidth = 1) + facet_wrap( ~ Gene) + 
  #ggtitle("Best BLAST matches for 28S and COI OTUs") +
  scale_x_continuous(breaks = pretty_breaks(5)) +
  xlab("Identical matches (%)") + ylab("Number of OTUs") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        strip.background = element_blank()) +
  scale_fill_manual(values = cols)            
p1
ggsave(p1, filename = "28S_COI_OTUs_BLAST_ID_barplot.pdf", width = 210, height = 120, units = "mm")

# Table of top BLAST matches;
bl <- bl[order(bl$Gene, bl$`Identical matches %`, bl$`Alignment length`), ]
bl.summary <- bl[bl$`Identical matches %` >= 97.0 & bl$`Alignment length` >= 450, ]
table(bl.summary$Gene) # 24 28S and 10 COI OTUs 
bl.summary$OTU <- rownames(bl.summary)
bl.summary <- bl.summary[, c(1,18,8,9,7,3:4,6)] # Bah, species ids are missing

###############################################################################
# Multivariate analyses of DNA barcoding OTU tables

# Generate Jaccard MDS ordination plots:
get_labels <- function(pts_mds){
  pts_mds$plots <- sapply(strsplit(rownames(pts_mds), "_"), "[", 1)
  pts_mds$plots <- gsub("-", " ", pts_mds$plots)
  pts_mds$plots <- ordered(pts_mds$plots, levels = c("Plot 1", "Plot 2", "Plot 3", "Plot 4", "Plot 5",
                                                     "Plot 6", "Plot 7", "Plot 8", "Plot 9", "Plot 10"))
  pts_mds$methods <- sapply(strsplit(rownames(pts_mds), "_"), "[[", 2)
  pts_mds$methods <- gsub("-", " ", pts_mds$methods)
  return(pts_mds)
}

# COI:
df.dist.COI <- vegdist(OTUtable.COI, method = "jaccard", binary = TRUE)
mds.COI <- metaMDS(df.dist.COI, k = 2)
pts_mds.COI <- as.data.frame(mds.COI$points)
stress_mds.COI <- mds.COI$stress
pts_mds.COI <- get_labels(pts_mds.COI)

# Plot COI MDS ordination:
p1 <- ggplot(pts_mds.COI, aes(x = MDS1, y = MDS2, color = plots, shape = methods)) + 
  geom_point(size = 4) + #scale_colour_discrete() +
  geom_text(aes(label = pts_mds.COI$plots), size = 3, vjust = 2) +
  theme(panel.grid = element_blank(), aspect.ratio = 1, legend.text=element_text(size=8)) +
  scale_colour_discrete(name = "Sample plot") + scale_shape(name = "Collection method") +
  scale_x_reverse() +
  ggtitle(paste("a. COI OTUs, Jaccard distance (stress:", round(stress_mds.COI, 2),")"))
p1

# 28S:
df.dist.28S <- vegdist(OTUtable.28S, method = "jaccard", binary = TRUE)
mds.28S <- metaMDS(df.dist.28S, k = 2)
pts_mds.28S <- as.data.frame(mds.28S$points)
stress_mds.28S <- mds.28S$stress
pts_mds.28S <- get_labels(pts_mds.28S)

# Plot 28S MDS ordination:
p2 <- ggplot(pts_mds.28S, aes(x = MDS1, y = MDS2, color = plots, shape = methods)) + 
  geom_point(size = 4) + #scale_colour_discrete() +
  geom_text(aes(label = pts_mds.28S$plots), size = 3, vjust = 2) +
  scale_colour_discrete(name = "Sample plot") + scale_shape(name = "Collection method") +
  theme(panel.grid = element_blank(), aspect.ratio = 1, legend.text=element_text(size=8)) +
  scale_x_reverse() +
  ggtitle(paste("b. 28S OTUs, Jaccard distance (stress:", round(stress_mds.28S, 2),")"))
p2

# Use Phyloseq to get UniFrac distances:
# Convert data to phyloseq format:
samples1 <- sample_data(samples) 

OTUtable.COI <- otu_table(OTUtable.COI, taxa_are_rows=FALSE)
taxatable.COI <- tax_table(as.matrix(taxatable.COI))
phylo_COI <- phyloseq(OTUtable.COI, samples1, tree.COI, taxatable.COI)
phylo_COI

OTUtable.28S <- otu_table(OTUtable.28S, taxa_are_rows=FALSE)
taxatable.28S <- tax_table(as.matrix(taxatable.28S))
phylo_28S <- phyloseq(OTUtable.28S, samples1, tree.28S, taxatable.28S)
phylo_28S

# COI UniFrac:
set.seed(345)
unifrac.COI <- UniFrac(phylo_COI, weighted=FALSE, normalized=TRUE)
mds.COI.unifrac <- metaMDS(unifrac.COI)
pts.COI.unifrac <- as.data.frame(mds.COI.unifrac$points)
stress_mds.COI.unifrac <- mds.COI.unifrac$stress
pts.COI.unifrac <- get_labels(pts.COI.unifrac)

# COI UniFrac MDS ordination:
p3 <- ggplot(pts.COI.unifrac, aes(x = MDS1, y = MDS2, color = plots, shape = methods)) + 
  geom_point(size = 4) + #scale_colour_discrete() +
  geom_text(aes(label = pts.COI.unifrac$plots), size = 3, vjust = 2) +
  scale_colour_discrete(name = "Sample plot") + scale_shape(name = "Collection method") +
  theme(panel.grid = element_blank(), aspect.ratio = 1, legend.text=element_text(size=8)) +
  scale_x_reverse() +
  ggtitle(paste("c. COI OTUs, unweighted UniFrac distance (stress:", round(stress_mds.COI.unifrac, 2),")"))
p3

# 28S UniFrac:
set.seed(456)
unifrac.28S <- UniFrac(phylo_28S, weighted=FALSE, normalized=TRUE)
mds.28S.unifrac <- metaMDS(unifrac.28S)
pts.28S.unifrac <- as.data.frame(mds.28S.unifrac$points)
stress_mds.28S.unifrac <- mds.28S.unifrac$stress
pts.28S.unifrac <- get_labels(pts.28S.unifrac)

# 28S UniFrac MDS ordination:
p4 <- ggplot(pts.28S.unifrac, aes(x = MDS1, y = MDS2, color = plots, shape = methods)) + 
  geom_point(size = 4) + #scale_colour_discrete() +
  geom_text(aes(label = pts.28S.unifrac$plots), size = 3, vjust = 2) +
  scale_colour_discrete(name = "Sample plot") + scale_shape(name = "Collection method") +
  theme(panel.grid = element_blank(), aspect.ratio = 1, legend.text=element_text(size=8)) +
  #scale_x_reverse() +
  ggtitle(paste("d. 28S OTUs, unweighted UniFrac distance (stress:", round(stress_mds.28S.unifrac, 2),")"))
p4

# Output all four MDS plots as one figure:
pdf("COI_28S_OTU_MDS_plots_x4.pdf", width = 29.7/2.54, height = 21/2.54, useDingbats = FALSE)
grid.arrange(p1, p2, p3, p4, ncol = 2)
dev.off()

###############################################################################
# Permanova tests:

adonis(df.dist.COI ~ Method * Elevation, samples)
adonis(df.dist.28S ~ Method * Elevation, samples)
adonis(unifrac.COI ~ Method * Elevation, samples)
adonis(unifrac.28S ~ Method * Elevation, samples)

###############################################################################
# Heatmap figures of OTU distribution using phyloseq
# Unable to order the y-axis correctly, so need to reorganise afterwards in Illustrator/Inkscape. (Grrr).
# Also adjust size and alignment of axis and facet labels

hm.COI <- plot_heatmap(phylo_COI, method = "MDS", distance = "bray", sample.label = "PlotN", 
          low="#66CCFF", high="#000033", na.value="white", sample.order = "PlotN") +
          #xlab("Sample plot") + # Doesn't work
          facet_grid(Order ~ Method, scales = "free", space = "free") +
          theme(strip.text.y = element_text(angle = 0), strip.background = element_blank(),
                axis.text.x = element_text(angle = 0, hjust = 0), axis.text.y = element_blank())

ggsave(hm.COI, file = "COI_OTUs_heatmap.pdf", width = 210, height = 366*0.8, units = "mm")

hm.28S <- plot_heatmap(phylo_28S, method = "MDS", distance = "bray", sample.label = "PlotN", 
          low="#66CCFF", high="#000033", na.value="white", sample.order = "PlotN") +
          #xlab("Sample plot") + # Doesn't work
          facet_grid(Order ~ Method, scales = "free", space = "free") +
          theme(strip.text.y = element_text(angle = 0), strip.background = element_blank(),
                axis.text.x = element_text(angle = 0, hjust = 0), axis.text.y = element_blank())

ggsave(hm.28S, file = "28S_OTUs_heatmap.pdf", width = 210, height = 247*0.8, units = "mm")

###############################################################################
# Procrustes comparison of ordination patterns

plot(mds.COI, type = "text")
plot(mds.28S, type = "text")
plot(mds.COI.unifrac, type = "text")
plot(mds.28S.unifrac, type = "text")

# procrustes(x = target matrix, y = rotated matrix)
# arrows point from rotated (y) to target matrix (x) (to.target = TRUE)
pro1 <- procrustes(mds.28S, mds.COI, symmetric = T, scale = T) # target = 28S, rotated = COI
pt1 <- protest(mds.28S, mds.COI, permutations = 9999)
plot(pro1)
text(pro1, display = c("target"))
text(pro1, display = c("rotated"))

pro2 <- procrustes(mds.28S.unifrac, mds.COI.unifrac, symmetric = T, scale = T)
pt2 <- protest(mds.28S.unifrac, mds.COI.unifrac, permutations = 9999)

pro3 <- procrustes(mds.COI.unifrac, mds.COI, symmetric = T, scale = T)
pt3 <- protest(mds.COI.unifrac, mds.COI, permutations = 9999)

pro4 <- procrustes(mds.28S.unifrac, mds.28S, symmetric = T, scale = T)
pt4 <- protest(mds.28S.unifrac, mds.28S, permutations = 9999)

# Output Procrustes results
pt.results <- list(pt1, pt2, pt3, pt4)
pt.data <- list()
n <- 1
for(pt in pt.results){
  call <- pt$call
  ss <- pt$ss
  t0 <- pt$t0
  sig <- pt$signif
  pt.bits <- c(call, ss, t0, sig)
  pt.data[[n]] <- pt.bits
  n <- n + 1
}
pt <- do.call("rbind", pt.data)
colnames(pt) <- c("call", "ss", "t0", "signif")
write.table(pt, file = "COI_28S_jaccard_unifrac_procrustes.txt", sep = "\t", quote = FALSE, row.names = FALSE)

# Plot procrustes errors
pro_plot <- function(pro = pro1, title = "y vs. x"){
  pts_mds <- data.frame(MDS1=pro$Yrot[,1], MDS2=pro$Yrot[,2], # rotated
                        xMDS1=pro$X[,1], xMDS2=pro$X[,2]) # target
  pts_mds <- get_labels(pts_mds)
  r1 = acos(pro$rotation[1,1])
  r2 = r1 + (pi/2)
  pro.plot <- ggplot(pts_mds) +
    geom_point(aes(x=MDS1, y=MDS2, colour = plots, shape = methods), size = 4) + # rotated (y)
    #geom_point(aes(x=xMDS1, y=xMDS2, colour = plots, shape = methods), size = 4) + # target (x)
    geom_segment(aes(x=MDS1, y=MDS2, xend=xMDS1, yend=xMDS2, colour = plots), 
                 arrow = arrow(length = unit(3, "mm")), show.legend = FALSE) + # Arrows from y (rotated) to x (target) 
    geom_hline(yintercept = 0, linetype = "dashed", color = "#aaaaaa") + geom_vline(xintercept = 0, linetype = "dashed", color = "#aaaaaa") +
    geom_abline(intercept = 0, slope = tan(r1), color = "#aaaaaa") + geom_abline(intercept = 0, slope = tan(r2), color = "#aaaaaa") +
    #geom_text(aes(x=MDS1, y=MDS2, label=rownames(pts_mds), size = 2)) + 
    theme(panel.grid = element_blank(), aspect.ratio = 1, legend.text=element_text(size=8)) +
    scale_colour_discrete(name = "Sample plot") + scale_shape(name = "Collection method") +
    ggtitle(title) + coord_fixed()
}

pro1.plot <- pro_plot(pro1, "a. COI Jaccard vs. 28S Jaccard distances")
pro1.plot <- pro1.plot + scale_x_reverse()
pro2.plot <- pro_plot(pro2, "b. COI UniFrac vs. 28S UniFrac distances")
pro2.plot <- pro2.plot + scale_y_reverse()
pro3.plot <- pro_plot(pro3, "c. COI Jaccard vs. UniFrac distances")
pro3.plot <- pro3.plot + scale_x_reverse()
pro4.plot <- pro_plot(pro4, "d. 28S Jaccard vs. UniFrac distances")
pro4.plot <- pro4.plot + scale_y_reverse()

pdf("COI_28S_OTU_procrustes_plots_x4.pdf", width = 29.7/2.54, height = 21/2.54, useDingbats = FALSE)
grid.arrange(pro1.plot, pro2.plot, pro3.plot, pro4.plot, ncol = 2)
dev.off()

###############################################################################
# Calculation of biodiversity estimates for DNA barcoded invertebrate OTUs 
# and 454 soil DNA invertebrate OTUs using SPECRICH2:
# Requires the number of species observed in exactly n sites, 
# and the number of species observed at each site. 
# Then the values can be entered into the online estimator at 
# https://www.mbr-pwrc.usgs.gov/software/specrich2.shtml

get_occurrences <- function(OTUtable = OTUtable.COI.byPlot, taxatable = taxatable.COI, 
                            tx = "Arthropoda", txrank = "Phylum"){
  # Set up an empty data frames to hold occurrences in 1 to 10 sites
  occs_res <- data.frame(matrix(ncol = 2, nrow = 10)) 
  rownames(occs_res) <- seq(1:10)
  colnames(occs_res) <- c("occs","nOTUs")
  
  # Limit to taxonomic group
  x <- taxatable[taxatable[txrank] ==  tx, ]
  z <- OTUtable[match(rownames(x), rownames(OTUtable)), ]
  
  # Get number of OTUs at each site
  nOTUs <- apply(z, 2, function(x) sum(x > 0)) 
  names(nOTUs) <- gsub("^Plot[-|_]|Plot", "", names(nOTUs)) 
  occs_res$nOTUs <- nOTUs[match(rownames(occs_res), names(nOTUs))]
  # Get number of OTUs in n sites
  z$occ <- apply(z, 1, function(x) sum(x > 0)) 
  occs <- as.data.frame(as.matrix(table(z$occ)))
  occs_res$occs <- occs$V1[match(rownames(occs_res), rownames(occs))]
  occs_res <- as.data.frame(t(occs_res))
  occs_res$group <- tx
  occs_res$rank <- txrank
  return(occs_res)
}

tx.groups <- list(c("Arthropoda","Phylum"),c("Mollusca","Phylum"),c("Annelida","Phylum"),c("Onychophora","Phylum"),
                  c("Arachnida","Class"),c("Malacostraca","Class"),c("Myriapoda", "Class"),c("Insecta","Class"),
                  c("Coleoptera","Order"),c("Diptera","Order"),c("Lepidoptera","Order"),c("Hymenoptera","Order"),
                  c("Hemiptera","Order"),c("Orthoptera","Order"))

# Combine the results from leaf litter and pitfall traps at each site:
combine_plots <- function(OTUtable){
  x <- melt(OTUtable, varnames = c("Sample","OTU"))
  x$Plot <- sapply(strsplit(as.character(x$Sample), "_"), "[[", 1)
  OTUtable <- dcast(x, OTU ~ Plot, value.var = "value", fun.aggregate = sum)
  rownames(OTUtable) <- OTUtable$OTU
  OTUtable$OTU <- NULL
  return(OTUtable)
}

# COI OTU occurrences:
occs_res <- data.frame(matrix(nrow = 0, ncol = 12))
colnames(occs_res) <- c(seq(1:10), "group", "rank")
OTUtable.COI.byPlot <- combine_plots(OTUtable.COI)
for(x in tx.groups){
  occs <- get_occurrences(OTUtable.COI.byPlot, taxatable.COI, tx = x[[1]], txrank = x[[2]])
  occs_res <- rbind(occs_res, occs)
}
occs_res[is.na(occs_res)] <- 0
occs_res$var <- gsub("[0-9]*", "", rownames(occs_res))
write.table(occs_res, file = "COI_OTU_occurrences.txt", sep = "\t", row.names = FALSE, quote = FALSE)

# 28S OTU occurrences:
occs_res <- data.frame(matrix(nrow = 0, ncol = 12))
colnames(occs_res) <- c(seq(1:10), "group", "rank")
OTUtable.28S.byPlot <- combine_plots(OTUtable.28S)
for(x in tx.groups){
  occs <- get_occurrences(OTUtable.28S.byPlot, taxatable.28S, tx = x[[1]], txrank = x[[2]])
  occs_res <- rbind(occs_res, occs)
}
occs_res[is.na(occs_res)] <- 0
occs_res$var <- gsub("[0-9]*", "", rownames(occs_res))
write.table(occs_res, file = "28S_OTU_occurrences.txt", sep = "\t", row.names = FALSE, quote = FALSE)

# 454 soil DNA metabarcoding COI OTU occurrences
occs_res <- data.frame(matrix(nrow = 0, ncol = 12))
colnames(occs_res) <- c(seq(1:10), "group", "rank")
for(x in tx.groups){
  occs <- get_occurrences(OTUtable.454, taxonomy.454, tx = x[[1]], txrank = x[[2]])
  occs_res <- rbind(occs_res, occs)
}
occs_res[is.na(occs_res)] <- 0
occs_res$var <- gsub("[0-9]*", "", rownames(occs_res))
write.table(occs_res, file = "COI_454_soil_DNA_OTU_occurrences.txt", sep = "\t", row.names = FALSE, quote = FALSE)

# Then enter the values and get estimates, one group at a time, at https://www.mbr-pwrc.usgs.gov/software/specrich2.shtml

###############################################################################
# Alpha, beta, gamma estimates ------------------------------------------------
#
#outfile <- "Sanger_CO1_28S_alpha_stats.txt"
#
# df <- t(OTUtable_CO1)
# 
# a_0 <- apply(df, 1, d, lev = "alpha", q = 0)
# a_1 <- apply(df, 1, d, lev = "alpha", q = 1)
# a_2 <- apply(df, 1, d, lev = "alpha", q = 2)
# g_0 <- d(df, lev = "gamma", q = 0)
# g_1 <- d(df, lev = "gamma", q = 1)
# g_2 <- d(df, lev = "gamma", q = 2)
# beta <- betadiver(df, "w")
# 
# alpha <- as.data.frame(cbind(a_0, a_1, a_2))
# write.table(alpha, file = outfile, append = TRUE, sep = "\t", quote = FALSE)
# gamma <-  as.data.frame(cbind(g_0, g_1, g_2)) 
# write.table(gamma, file = outfile, append = TRUE, sep = "\t", quote = FALSE)
# beta <- as.matrix(beta)
# write.table(beta, file = outfile, append = TRUE, sep = "\t", quote = FALSE)
# 
# df <- read.table("CO1_28S_alpha_vs_elevation.txt", 
#                 header = TRUE, sep = "\t", check.names=FALSE)
# df1 <- df[,1:8]
# df2 <- df[,c(1:5, 9:11)]
# 
# df.l <- melt(df2, id.vars = c("Sample", "Plot", "Method", "Gene", "Elevation"))
# ggplot(df.l, aes(x = Elevation, y = value, colour = variable)) + 
#       geom_point(shape = 1, aes(group = variable)) + 
#       geom_smooth(method=lm, fill = NA) +
#       facet_wrap(Gene~Method) +
#       theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
