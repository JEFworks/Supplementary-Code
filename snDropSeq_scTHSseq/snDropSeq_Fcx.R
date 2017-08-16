## load data
load("../data/Fcx_06062017_all_counts_PC_no-MT.RData")
dim(countMatrix)

## annotations from individual analysis
load("../R/Fcx_results_04052017.RData")
annot <- Fcx.annot

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
p2$makeKnnGraph(k = 30, type = 'PCA', center = T, weight.type = 'none', n.cores = 4, distance = 'cosine')

# Generate an embedding with tSNE on the basis of the PCA reduction
p2$getEmbedding(type = 'PCA', embeddingType = 'tSNE', verbose=TRUE, perplexity=100)
#p2$getEmbedding(type='PCA',M=5,perplexity=30,gamma=1)

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
## Color by time point
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
                 main="New Annot")
p2$plotEmbedding(type = 'PCA',
                 colors=p2$depth,
                 embeddingType = 'tSNE',
                 mark.clusters = T,
                 main="Depth")

############## Now retest
good.cells <- na.omit(retest)
length(good.cells)
dim(countMatrix)

Fcx <- Pagoda2$new(x = countMatrix[, names(good.cells)], n.cores = 4, trim=10, batch=countMatrix.batch[names(good.cells)])

# Adjust the variance
Fcx$adjustVariance(plot = T, gam.k = 10)

# Calculate a PCA reduction with the number of PCs specified by nPCs
# and using only the n.odgenes overdispersed genes -- in this case 2000
Fcx$calculatePcaReduction(nPcs = 150, n.odgenes = 2.e3, maxit=1000)

# Generate K-nearest neighbour graph
Fcx$makeKnnGraph(k = 15, type = 'PCA', center = T, weight.type = 'none', n.cores = 4, distance = 'cosine')

# Generate an embedding with tSNE on the basis of the PCA reduction
Fcx$getEmbedding(type = 'PCA', embeddingType = 'tSNE', verbose=TRUE, perplexity=30)
#Fcx$getEmbedding(type='PCA',M=5,perplexity=100,gamma=1)

# Identify clusters using the infomap.community method
# on the basis of the reduction called 'PCA' (generated above)
# Save the resulting clustering as 'infomap'
Fcx$getKnnClusters(method = infomap.community, type = 'PCA', name = 'infomap')

# Do an independent identification of clusters using the
# multilevel community algorithm again using the PCA reduction
# and save it as 'multilevel'. This does not overwrite the
# previous clustering
Fcx$getKnnClusters(method = multilevel.community, type = 'PCA', name='multilevel')

# Do yet another clustering
Fcx$getKnnClusters(method = walktrap.community, type = 'PCA', name='walktrap')

# Plot
par(mfrow=c(1,3))
Fcx$plotEmbedding(type = 'PCA',
                 clusterType = 'infomap',
                 embeddingType = 'tSNE',
                 mark.clusters = T,
                 main="infomap",
                 min.group.size=30)

Fcx$plotEmbedding(type = 'PCA',
                 clusterType = 'multilevel',
                 embeddingType = 'tSNE',
                 mark.clusters = T,
                 main="multilevel",
                 min.group.size=30)

Fcx$plotEmbedding(type = 'PCA',
                 clusterType = 'walktrap',
                 embeddingType = 'tSNE',
                 mark.clusters = T,
                 main="walktrap",
                 min.group.size=30)

newannot <- Fcx$clusters$PCA$infomap
newannot[newannot %in% which(table(newannot)<30)] <- NA

## Final plot
par(mfrow=c(1,4))
# Color by time point
Fcx$plotEmbedding(type = 'PCA',
                 groups=countMatrix.batch,
                 embeddingType = 'tSNE',
                 mark.clusters = T,
                 main="Batch",
                 mark.cluster.cex = 1,
                 alpha=0.1)
Fcx$plotEmbedding(type = 'PCA',
                 groups=annot,
                 embeddingType = 'tSNE',
                 mark.clusters = T,
                 main="Annot",
                 mark.cluster.cex = 1,
                 alpha=0.1)
Fcx$plotEmbedding(type = 'PCA',
                 groups=newannot,
                 embeddingType = 'tSNE',
                 mark.clusters = T,
                 main="multilevel",
                 mark.cluster.cex = 1,
                 alpha=0.1)
Fcx$plotEmbedding(type = 'PCA',
                 colors=Fcx$depth,
                 embeddingType = 'tSNE',
                 mark.clusters = T,
                 main="Depth",
                 mark.cluster.cex = 1,
                 alpha=0.1)



save.image('Fcx.06082017.RData')  

Fcx.annot <- as.character(newannot)
names(Fcx.annot) <- names(newannot)
Fcx.annot[newannot == 1] <- 'Ex1'
Fcx.annot[newannot == 12] <- 'Ex2a'
Fcx.annot[newannot == 11] <- 'Ex2b'
Fcx.annot[newannot == 3] <- 'Ex4'
Fcx.annot[newannot == 5] <- 'Ex5'
Fcx.annot[newannot == 16] <- 'Ex8'
Fcx.annot[newannot == 18] <- 'Ex6a'
Fcx.annot[newannot == 6] <- 'Ex6b'
Fcx.annot[newannot == 9] <- 'In3'
Fcx.annot[newannot == 14] <- 'In1'
Fcx.annot[newannot == 17] <- 'In4'
Fcx.annot[newannot == 13] <- 'In5'
Fcx.annot[newannot == 19] <- 'In6a'
Fcx.annot[newannot == 8] <- 'In6b'
Fcx.annot[newannot == 22] <- 'In7'
Fcx.annot[newannot == 7] <- 'In8'
Fcx.annot[newannot == 2] <- 'Oli'
Fcx.annot[newannot == 15] <- 'Mic'
Fcx.annot[newannot == 20] <- 'End'
Fcx.annot[newannot == 21] <- 'Astb'
Fcx.annot[newannot == 4] <- 'Asta'
Fcx.annot[newannot == 10] <- 'Opc'
Fcx.annot <- factor(Fcx.annot)
head(Fcx.annot)

Fcx$plotEmbedding(type = 'PCA',
                  groups=Fcx.annot,
                 embeddingType = 'tSNE',
                 mark.clusters = T,
                 mark.cluster.cex = 1,
                 alpha=0.1)

save(Fcx, Fcx.annot, file='Fcx_results_06082017.RData')



load('Fcx.06082017.RData')  


load('Fcx_results_06082017.RData')
library("largeVis", lib.loc="/usr/local/lib/R/site-library")
library(pagoda2)
load('final.annot.RData')

final.fcx <- final[rownames(Fcx$counts)]
final.fcx[final.fcx %in% names(which(table(final.fcx) < 30))] <- NA
table(final.fcx)


Fcx$plotEmbedding(type = 'PCA',
                                    groups=final.fcx,
                                    embeddingType = 'tSNE',
                                    mark.clusters = T,
                                    main="Final Annot",
                                    mark.cluster.cex = 1,
                                    alpha=0.1)
set.seed(3)
Fcx$plotEmbedding(type = 'PCA',
                                    groups=final.fcx,
                                    embeddingType = 'tSNE',
                                    mark.clusters = T,
                                    main="Final Annot",
                                    mark.cluster.cex = 1,
                                    alpha=0.1, shuffle.colors=TRUE)
