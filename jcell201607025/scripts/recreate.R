# Recreate published figures as a sanity check

load('../data/annot.RData')
load('../data/fpkm.RData')

# get 226 cells from 12 and 13 week postconception human neocortex specimens
head(cell_info)
subcell_info <- cell_info[cell_info$source_name_s == "fetal human cortex 1" | cell_info$source_name_s == "fetal human cortex 2",]
head(subcell_info)
rownames(subcell_info) %in% rownames(cell_name)
subcell_names <- as.character(cell_name[rownames(subcell_info),])

# convert to same names as in fpkm matrix
head(fpkm)
subcell_names_nice <- subcell_names
subcell_names_nice <- gsub("fetal1_", "", subcell_names_nice)
subcell_names_nice <- gsub("fetal2_", "", subcell_names_nice)
subcell_names_nice <- gsub("fetal3_", "", subcell_names_nice)
#subcell_names_nice <- gsub("_c1", "", subcell_names_nice)
#subcell_names_nice <- gsub("_c2", "", subcell_names_nice)
subcell_names_nice

#foo <- colnames(fpkm)
#foo <- gsub("_c1", "", foo)
#foo <- gsub("_c2", "", foo)
#colnames(fpkm) <- foo

subcell_names_nice[which(!(subcell_names_nice %in% colnames(fpkm)))]
colnames(fpkm)[grepl("F5", colnames(fpkm))]
colnames(fpkm)
# weird, there are no F5 fetal 13wpc
subcell_names_nice <- intersect(subcell_names_nice, colnames(fpkm))

sub_fpkm <- fpkm[, subcell_names_nice]
head(sub_fpkm)

# Look for highly variable genes (var > 0.5) expressed in more than 2 cells
vi <- rowSums(sub_fpkm > 0) > 2
table(vi)
sub_fpkm <- sub_fpkm[vi,]
v <- apply(sub_fpkm, 1, var)
vi <- v > 0.5
table(vi)
sub_fpkm <- sub_fpkm[vi,]

# Recreating Figure S1
mat <- log10(t(sub_fpkm)+1)
base.pca <- prcomp(mat)
batch <- lapply(rownames(mat), function(x) {
    y <- strsplit(x, '_')[[1]]
    y <- y[-1]
    y <- y[-1]
    paste(y, collapse="_")
})
batch <- as.factor(unlist(batch))
batchcol <- c("yellowgreen", "seagreen", "steelblue")[batch]
names(batchcol) <- names(batch)
plot(base.pca$x[,1], base.pca$x[,2], pch=16, col=batchcol, main="Recreation of Fig. S1A", xlab="PC1", ylab="PC2")
abline(v=2, lty = 2, lwd = 2, col="red")

vi <- unlist(lapply(rownames(mat), function(n) { which(grepl(n, as.character(cell_name[,1]))) }))
sra_info <- cell_info[rownames(cell_name)[vi],]$Run_s

names(batch) <- names(batchcol) <- as.character(sra_info)
save(batch, batchcol, file="batch.RData")

# Which are neurons which are NPCs?
neurons <- names(which(base.pca$x[,1] < 2))
npcs <- names(which(base.pca$x[,1] > 2))
groups <- c(rep('neurons', length(neurons)), rep('npcs', length(npcs)))
names(groups) <- c(neurons, npcs)
groupcol <- c("coral4", "olivedrab3", "skyblue2", "slateblue3")[as.integer(as.factor(groups))]
names(groupcol) <- names(groups)
plot(base.pca$x[,1], base.pca$x[,2], pch=16, col=groupcol[rownames(base.pca$x)])

# Get SRA names
head(cell_info)
vi <- unlist(lapply(neurons, function(n) { which(grepl(n, as.character(cell_name[,1]))) }))
neuron_info <- cell_info[rownames(cell_name)[vi],]

vi <- unlist(lapply(npcs, function(n) { which(grepl(n, as.character(cell_name[,1]))) }))
length(npcs)
length(vi)
npc_info <- cell_info[rownames(cell_name)[vi],]

neuron_sra <- as.character(neuron_info$Run_s)
npc_sra <- as.character(npc_info$Run_s)

save(neuron_sra, npc_sra, neurons, npcs, file="../data/cell_subgroup.RData")

load("../data/cell_subgroup.RData")
dput(paste("samtools merge neurons.bam", paste(paste0(neuron_sra, '.bam'), collapse=" "), collapse=" "))
dput(paste("samtools merge npcs.bam", paste(paste0(npc_sra, '.bam'), collapse=" "), collapse=" "))


load("../data/cd.RData")
library(Rsamtools)
seqnames <- paste("../tophat/", colnames(cd), '.bam', sep="")
#seqnames <- paste("../tophat/", npc_sra, '.bam', sep="")
what <- "qwidth"
param <- ScanBamParam(what=what)
qlen <- mclapply(seqnames, function(x) {
    print(x)
    b <- BamFile(x)
    aln <- scanBam(b, param=param)
    aln <- aln[[1]]
    q <- aln$qwidth[1]
    print(q)
    return(q)
}, mc.cores=10)
# get only the ones of length 100
qlen <- unlist(qlen)
names(qlen) <- colnames(cd)
save(qlen, file="qlen.RData")
load("qlen.RData")

table(cutree(hc,2))
npc_sra_100 = names(which(qlen[npc_sra]==100))
neuron_sra_100 = names(which(qlen[neuron_sra] == 100))
length(npc_sra_100)
length(neuron_sra_100)
# keep all NPCs, want same number of neurons
load('tam.RData')
plot(hc)
o <- 1:length(hc$labels)
names(o) <- hc$labels[hc$order]
npc_sra_100_sub <- names(sort(o[npc_sra_100], decreasing=FALSE))[1:45]
neuron_sra_100_sub <- names(sort(o[neuron_sra_100], decreasing=TRUE)[1:45])
neuron_imm_sra_100_sub <- names(sort(o[neuron_sra_100], decreasing=FALSE)[1:45])

dput(paste("samtools merge neurons100sub45.bam", paste(paste0(neuron_sra_100_sub, '.bam'), collapse=" "), collapse=" "))
dput(paste("samtools merge neuronsimm100sub45.bam", paste(paste0(neuron_imm_sra_100_sub, '.bam'), collapse=" "), collapse=" "))
dput(paste("samtools merge npcs100sub45.bam", paste(paste0(npc_sra_100_sub, '.bam'), collapse=" "), collapse=" "))





# Get just organoid data

load('../data/annot.RData')
load('../data/fpkm.RData')

head(cell_info)
table(cell_info$tissue_s)
dput(as.character(cell_info[cell_info$tissue_s == "Microdissected cortical-like ventricle from cerebral organoid", "Run_s"]))
