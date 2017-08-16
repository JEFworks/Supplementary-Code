## load data
load("../data/Occ_06062017_all_counts_PC_no-MT.RData")
dim(countMatrix)

## annotations from individual analysis
load("../R/Occ_results_04052017.RData")
annot <- DS.annot

## Batch annotations
countMatrix.batch <- as.factor(gsub("_.*","",colnames(countMatrix)))
names(countMatrix.batch) <- colnames(countMatrix)

## QC 
t0 <- 5000
t1 <- 300
vc <- colSums(countMatrix>0)>t1 & colSums(countMatrix)<t0
table(vc)
countMatrix <- countMatrix[, vc]
countMatrix.batch <- countMatrix.batch[vc]

########### PAGODA
library("largeVis", lib.loc="/usr/local/lib/R/site-library")
library(pagoda2)

# Generate a new pagoda2 object
dim(countMatrix)
p2 <- Pagoda2$new(x = countMatrix, n.cores = 4, trim=10, batch=countMatrix.batch)

# Adjust the variance
p2$adjustVariance(plot = T, gam.k = 10)

# Calculate a PCA reduction with the number of PCs specified by nPCs
# and using only the n.odgenes overdispersed genes -- in this case 2000
p2$calculatePcaReduction(nPcs = 150, n.odgenes = 2.e3, maxit=1000)

# Generate K-nearest neighbour graph
p2$makeKnnGraph(k = 15, type = 'PCA', center = T, weight.type = 'none', n.cores = 4, distance = 'cosine')

# Generate an embedding with tSNE on the basis of the PCA reduction
p2$getEmbedding(type = 'PCA', embeddingType = 'tSNE', verbose=TRUE, perplexity=30)
#p2$getEmbedding(type='PCA',M=5,perplexity=500,gamma=1)

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

# Plot
par(mfrow=c(1,3))
p2$plotEmbedding(type = 'PCA',
                  clusterType = 'infomap',
                 embeddingType = 'tSNE',
                 mark.clusters = T,
                 main="infomap",
                 min.group.size=50)

p2$plotEmbedding(type = 'PCA',
                 clusterType = 'multilevel',
                 embeddingType = 'tSNE',
                 mark.clusters = T,
                 main="multilevel",
                 min.group.size=50)

p2$plotEmbedding(type = 'PCA',
                 clusterType = 'walktrap',
                 embeddingType = 'tSNE',
                 mark.clusters = T,
                 main="walktrap",
                 min.group.size=50)

cl0 <- p2$clusters$PCA$infomap
cl0[cl0 %in% which(table(cl0)<50)] <- NA
retest <- cl0

## Final plot
par(mfrow=c(1,4))
# Color by time point
p2$plotEmbedding(type = 'PCA',
                 groups=countMatrix.batch,
                 embeddingType = 'tSNE',
                 mark.clusters = T,
                 main="Batch")
p2$plotEmbedding(type = 'PCA',
                 groups=annot,
                 embeddingType = 'tSNE',
                 mark.clusters = T,
                 main="Annot")
p2$plotEmbedding(type = 'PCA',
                 groups=retest,
                     embeddingType = 'tSNE',
                     mark.clusters = T,
                     main="multilevel")
p2$plotEmbedding(type = 'PCA',
                 colors=p2$depth,
                 embeddingType = 'tSNE',
                 mark.clusters = T,
                 main="Depth")

############## Now retest
good.cells <- na.omit(retest)
length(good.cells)
dim(countMatrix)

Occ <- Pagoda2$new(x = countMatrix[, names(good.cells)], n.cores = 4, trim=10, batch=countMatrix.batch[names(good.cells)])

# Adjust the variance
Occ$adjustVariance(plot = T, gam.k = 10)

# Calculate a PCA reduction with the number of PCs specified by nPCs
# and using only the n.odgenes overdispersed genes -- in this case 2000
Occ$calculatePcaReduction(nPcs = 150, n.odgenes = 2.e3, maxit=1000)

# Generate K-nearest neighbour graph
Occ$makeKnnGraph(k = 15, type = 'PCA', center = T, weight.type = 'none', n.cores = 4, distance = 'cosine')

# Generate an embedding with tSNE on the basis of the PCA reduction
Occ$getEmbedding(type = 'PCA', embeddingType = 'tSNE', verbose=TRUE, perplexity=30)
#Occ$getEmbedding(type='PCA',M=5,perplexity=100,gamma=1)

# Identify clusters using the infomap.community method
# on the basis of the reduction called 'PCA' (generated above)
# Save the resulting clustering as 'infomap'
Occ$getKnnClusters(method = infomap.community, type = 'PCA', name = 'infomap')

# Do an independent identification of clusters using the
# multilevel community algorithm again using the PCA reduction
# and save it as 'multilevel'. This does not overwrite the
# previous clustering
Occ$getKnnClusters(method = multilevel.community, type = 'PCA', name='multilevel')

# Do yet another clustering
Occ$getKnnClusters(method = walktrap.community, type = 'PCA', name='walktrap')

# Plot
par(mfrow=c(1,3))
Occ$plotEmbedding(type = 'PCA',
                 clusterType = 'infomap',
                 embeddingType = 'tSNE',
                 mark.clusters = T,
                 main="infomap",
                 min.group.size=50)

Occ$plotEmbedding(type = 'PCA',
                 clusterType = 'multilevel',
                 embeddingType = 'tSNE',
                 mark.clusters = T,
                 main="multilevel",
                 min.group.size=50)

Occ$plotEmbedding(type = 'PCA',
                 clusterType = 'walktrap',
                 embeddingType = 'tSNE',
                 mark.clusters = T,
                 main="walktrap",
                 min.group.size=50)

newannot <- Occ$clusters$PCA$infomap
newannot[newannot %in% which(table(newannot)<30)] <- NA

## Final plot
par(mfrow=c(1,4))
# Color by time point
Occ$plotEmbedding(type = 'PCA',
                 groups=countMatrix.batch,
                 embeddingType = 'tSNE',
                 mark.clusters = T,
                 main="Batch",
                 mark.cluster.cex = 1,
                 alpha=0.1) 
Occ$plotEmbedding(type = 'PCA',
                 groups=annot,
                 embeddingType = 'tSNE',
                 mark.clusters = T,
                 main="Annot",
                 mark.cluster.cex = 1,
                 alpha=0.1) 
Occ$plotEmbedding(type = 'PCA',
                 groups=newannot,
                     embeddingType = 'tSNE',
                     mark.clusters = T,
                     main="multilevel",
                 mark.cluster.cex = 1,
                 alpha=0.1) 
Occ$plotEmbedding(type = 'PCA',
                 colors=Occ$depth,
                 embeddingType = 'tSNE',
                 mark.clusters = T,
                 main="Depth",
                 mark.cluster.cex = 1,
                 alpha=0.1) 

save.image('Occ.06082017.RData')  

Occ.annot <- as.character(newannot)
names(Occ.annot) <- names(newannot)
Occ.annot[newannot == 1] <- 'Ex1a'
Occ.annot[newannot == 2] <- 'Oli'
Occ.annot[newannot == 4] <- 'Asta'
Occ.annot[newannot == 29] <- 'Astb'
Occ.annot[newannot == 9] <- 'Ex1b'
Occ.annot[newannot == 21] <- 'Ex3c'
Occ.annot[newannot == 22] <- 'Ex3b'
Occ.annot[newannot == 5] <- 'Ex3a'
Occ.annot[newannot == 3] <- 'Ex4'
Occ.annot[newannot == 10] <- 'Ex5a2'
Occ.annot[newannot == 19] <- 'Ex5a1'
Occ.annot[newannot == 11] <- 'Ex5b'
Occ.annot[newannot == 12] <- 'Ex6b'
Occ.annot[newannot == 20] <- 'Ex8'
Occ.annot[newannot == 25] <- 'Ex6a'
Occ.annot[newannot == 28] <- 'Ex6c'
Occ.annot[newannot == 23] <- 'Ina'
Occ.annot[newannot == 6] <- 'In6b'
Occ.annot[newannot == 16] <- 'In7'
Occ.annot[newannot == 7] <- 'In8'
Occ.annot[newannot == 14] <- 'In5'
Occ.annot[newannot == 24] <- 'In4'
Occ.annot[newannot == 18] <- 'In2'
Occ.annot[newannot == 15] <- 'In1'
Occ.annot[newannot == 13] <- 'In3'
Occ.annot[newannot == 27] <- 'End'
Occ.annot[newannot == 26] <- 'Per'
Occ.annot[newannot == 8] <- 'Opc'
Occ.annot[newannot == 17] <- 'Mic'
Occ.annot <- factor(Occ.annot)
head(Occ.annot)

Occ$plotEmbedding(type = 'PCA',
                  groups=Occ.annot,
                  embeddingType = 'tSNE',
                  mark.clusters = T,
                  mark.cluster.cex = 1,
                  alpha=0.1) 

save(Occ, Occ.annot, file='Occ_results_06082017.RData')



#load('Occ.06082017.RData')
load('Occ_results_06082017.RData')
library("largeVis", lib.loc="/usr/local/lib/R/site-library")
library(pagoda2)
load('final.annot.RData')

final.occ <- final[rownames(Occ$counts)]
final.occ[final.occ %in% names(which(table(final.occ) < 30))] <- NA
table(final.occ)


Occ$plotEmbedding(type = 'PCA',
                  groups=final.occ,
                  embeddingType = 'tSNE',
                  mark.clusters = T,
                  main="Final Annot",
                  mark.cluster.cex = 1,
                  alpha=0.1) 
set.seed(3)
Occ$plotEmbedding(type = 'PCA',
                  groups=final.occ,
                  embeddingType = 'tSNE',
                  mark.clusters = T,
                  main="Final Annot",
                  mark.cluster.cex = 1,
                  alpha=0.1, shuffle.colors=TRUE) 
