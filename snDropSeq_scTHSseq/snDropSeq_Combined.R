plot(1,1)
## ###################### read in data
## require(Matrix)
## x <- read.table("../data-raw/Occ_06062017_all_counts_PC_no-MT.tsv.gz",header=T,sep="\t",stringsAsFactors=F)
## countMatrix <- Matrix(as.matrix(x),sparse=T)
## rm(x)
## head(countMatrix)
## save(countMatrix, file="../data/Occ_06062017_all_counts_PC_no-MT.RData")
## x <- read.table("../data-raw/Fcx_06062017_all_counts_PC_no-MT.tsv.gz",header=T,sep="\t",stringsAsFactors=F)
## countMatrix <- Matrix(as.matrix(x),sparse=T)
## rm(x)
## head(countMatrix)
## save(countMatrix, file="../data/Fcx_06062017_all_counts_PC_no-MT.RData")

###################### combine all data
load("../data/Occ_06062017_all_counts_PC_no-MT.RData")
cd.Occ <- countMatrix
load("../data/Fcx_06062017_all_counts_PC_no-MT.RData")
cd.Fcx <- countMatrix
load("../data/Cer_05042017_all_counts_PC_no-MT.RData")
cd.Cer <- countMatrix

require(Matrix)
genes.int <- intersect(intersect(rownames(cd.Occ), rownames(cd.Fcx)), rownames(cd.Cer))
length(genes.int)
countMatrix <- cbind(cd.Occ[genes.int,], cd.Fcx[genes.int,], cd.Cer[genes.int,])
dim(countMatrix)
sample <- factor(c(
    rep('Occ', ncol(cd.Occ)),
    rep('Fcx', ncol(cd.Fcx)),
    rep('Cer', ncol(cd.Cer))
))
names(sample) <- c(colnames(cd.Occ), colnames(cd.Fcx), colnames(cd.Cer))
head(sample)

## annotations from individual analysis
load("Occ_results_06082017.RData")
head(Occ.annot)
load("Fcx_results_06082017.RData")
head(Fcx.annot)
load("Cer_results_06092017-2.RData")
head(Cer.annot)
annot <- c(as.character(Occ.annot), as.character(Fcx.annot), as.character(Cer.annot))
names(annot) <- c(names(Occ.annot), names(Fcx.annot), names(Cer.annot))
annot <- factor(annot)
new.annot <- annot

## Batch annotations
countMatrix.batch <- as.factor(gsub("_.*","",colnames(countMatrix)))
names(countMatrix.batch) <- colnames(countMatrix)

## QC 
## t0 <- 5000
## t1 <- 300
## vc <- colSums(countMatrix>0)>t1 & colSums(countMatrix)<t0
## table(vc)
## countMatrix <- countMatrix[, vc]
## countMatrix.batch <- countMatrix.batch[vc]

vc <- names(na.omit(annot))
countMatrix <- countMatrix[, vc]
countMatrix.batch <- countMatrix.batch[vc]
dim(countMatrix)

########### PAGODA
library("largeVis", lib.loc="/usr/local/lib/R/site-library")
library(pagoda2)

# Generate a new pagoda2 object
dim(countMatrix)
p2 <- Pagoda2$new(x = countMatrix, n.cores = 4, trim=10, batch=countMatrix.batch)
#p2 <- Pagoda2$new(x = countMatrix, n.cores = 4, trim=20, batch=sample)

# Adjust the variance
p2$adjustVariance(plot = T, gam.k = 10)

# Calculate a PCA reduction with the number of PCs specified by nPCs
# and using only the n.odgenes overdispersed genes -- in this case 2000
p2$calculatePcaReduction(nPcs = 150, n.odgenes = 2.e3, maxit=10000)

# Generate K-nearest neighbour graph
p2$makeKnnGraph(k = 15, type = 'PCA', center = T, weight.type = 'none', n.cores = 4, distance = 'cosine')

# Generate an embedding with tSNE on the basis of the PCA reduction
p2$getEmbedding(type = 'PCA', embeddingType = 'tSNE', verbose=TRUE, perplexity=30)
#p2$getEmbedding(type='PCA',M=1/500,perplexity=500,gamma=1)

# Identify clusters using the infomap.community method
# on the basis of the reduction called 'PCA' (generated above)
# Save the resulting clustering as 'infomap'
p2$getKnnClusters(method = infomap.community, type = 'PCA', name = 'infomap')

# Do an independent identification of clusters using the
# multilevel community algorithm again using the PCA reduction
# and save it as 'multilevel'. This does not overwrite the
# previous clustering
p2$getKnnClusters(method = multilevel.community, type = 'PCA', name='multilevel')

# Do yet another clustering
p2$getKnnClusters(method = walktrap.community, type = 'PCA', name='walktrap')

new.annot <- annot
load("../R/Occ_results_04052017.RData")
Occ.annot <- DS.annot
head(Occ.annot)
load("../R/Fcx_results_04052017.RData")
head(Fcx.annot)
load("Cer_results_06092017.RData")
head(Cer.annot)
annot <- c(as.character(Occ.annot), as.character(Fcx.annot), as.character(Cer.annot))
names(annot) <- c(names(Occ.annot), names(Fcx.annot), names(Cer.annot))
old.annot <- factor(annot)

final <- final.annot <- p2$clusters$PCA$infomap
final <- as.character(final)
final[final.annot %in% which(table(final.annot)<50)] <- NA
final[final.annot == 2] <- 'Oli'
final[final.annot == 1] <- 'Ex1a'
final[final.annot == 8] <- 'Ex1b'
final[final.annot == 24] <- 'Ex1c'
final[final.annot == 13] <- 'Ex1d'
final[final.annot == 16] <- 'Ex2a'
final[final.annot == 36] <- 'Ex2b'
final[final.annot == 21] <- 'Ex3a'
final[final.annot == 38] <- 'Ex3b'
final[final.annot == 5] <- 'Ex4'
final[final.annot == 12] <- 'Ex5a'
final[final.annot == 27] <- 'Ex5e'
final[final.annot == 22] <- 'Ex5d'
final[final.annot == 19] <- 'Ex5b'
final[final.annot == 42] <- 'Ex5f'
final[final.annot == 26] <- 'Ex5c'
final[final.annot == 30] <- 'Ex6a'
final[final.annot == 25] <- 'Ex8'
final[final.annot == 10] <- 'Ex6b'
final[final.annot == 33] <- 'Ex3b'

final[final.annot == 14] <- 'In1c'
final[final.annot == 15] <- 'In3'
final[final.annot == 28] <- 'In1b'
final[final.annot == 37] <- 'In1a'
final[final.annot == 45] <- 'In2'
final[final.annot == 29] <- 'In4'
final[final.annot == 17] <- 'In5'
final[final.annot == 31] <- 'In6a'
final[final.annot == 6] <- 'In6b'
final[final.annot == 7] <- 'In8'
final[final.annot == 20] <- 'In7'

final[final.annot == 4] <- 'Ast'
final[final.annot == 35] <- 'Ast'
final[final.annot == 3] <- 'Gran'
final[final.annot == 11] <- 'Purk'
final[final.annot == 23] <- 'Ast_cer'
final[final.annot == 9] <- 'Opc'
final[final.annot == 18] <- 'Mic'
final[final.annot == 40] <- 'Mic'
final[final.annot == 2] <- 'Oli'
final[final.annot == 39] <- 'Oli'
final[final.annot == 34] <- 'Opc_cer'
final[final.annot == 38] <- 'End'
final[final.annot == 32] <- 'Per'

final[final.annot == 41] <- 'Ex6c'
final[final.annot == 44] <- NA
final[final.annot == 46] <- NA
final[final.annot == 47] <- NA
final[final.annot == 48] <- NA
final[final.annot == 49] <- NA
final[final.annot == 50] <- NA
final[final.annot == 51] <- NA
final[final.annot == 52] <- NA
final[final.annot == 43] <- NA

names(final) <- names(final.annot)
table(final)

final <- factor(final)
table(final)

                

pdf('test.pdf', width=30, height=5)
## # Plot
## par(mfrow=c(1,3))
## p2$plotEmbedding(type = 'PCA',
##                  clusterType = 'infomap',
##                  embeddingType = 'tSNE',
##                  mark.clusters = T,
##                  main="infomap",
##                  min.group.size=100,
##                  alpha=0.1)

## p2$plotEmbedding(type = 'PCA',
##                  clusterType = 'multilevel',
##                  embeddingType = 'tSNE',
##                  mark.clusters = T,
##                  main="multilevel",
##                  min.group.size=100,
##                  alpha=0.1)

## p2$plotEmbedding(type = 'PCA',
##                  clusterType = 'walktrap',
##                  embeddingType = 'tSNE',
##                  mark.clusters = T,
##                  main="walktrap",
##                  min.group.size=100,
##                  alpha=0.1)

## Final plot
par(mfrow=c(1,6))
## Color by time point
p2$plotEmbedding(type = 'PCA',
                 groups=sample,
                 embeddingType = 'tSNE',
                 mark.clusters = T,
                 main="Sample",
                 mark.cluster.cex=1,
                 alpha=0.1)
p2$plotEmbedding(type = 'PCA',
                 groups=countMatrix.batch,
                 embeddingType = 'tSNE',
                 mark.clusters = T,
                 main="Batch",
                 mark.cluster.cex=1,
                 alpha=0.1)
p2$plotEmbedding(type = 'PCA',
                 groups=old.annot,
                 embeddingType = 'tSNE',
                 mark.clusters = T,
                 main="Old Annot",
                 mark.cluster.cex=1,
                 alpha=0.1)
p2$plotEmbedding(type = 'PCA',
                 groups=new.annot,
                 embeddingType = 'tSNE',
                 mark.clusters = T,
                 main="Individual Annot",
                 mark.cluster.cex=1,
                 alpha=0.1)
p2$plotEmbedding(type = 'PCA',
                 groups=final,
                 embeddingType = 'tSNE',
                 mark.clusters = T,
                 main="Final Annot",
                 mark.cluster.cex=1,
                 alpha=0.1)
p2$plotEmbedding(type = 'PCA',
                 colors=p2$depth,
                 embeddingType = 'tSNE',
                 main='depth',
                 alpha=0.1)

dev.off()

save(final, p2, file='final.annot.RData')

save.image('comb.06102017.RData')
