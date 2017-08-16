## load data
load("../data/Cer_05042017_all_counts_PC_no-MT.RData")
dim(countMatrix)

## annotations from individual analysis
load("../R3/Cer_results_05052017.RData")
annot <- Cer.annot

## Batch annotations
countMatrix.batch <- as.factor(gsub("_.*","",colnames(countMatrix)))
names(countMatrix.batch) <- colnames(countMatrix)
par(mfrow=c(3, length(levels(countMatrix.batch))/3), mar=rep(5,4))
lapply(levels(countMatrix.batch), function(x) {
    vi <- countMatrix.batch==x
    hist(log10(colSums(countMatrix[,vi]>0)+1), main=x)
    #abline(v=median(log10(colSums(countMatrix[,vi]>0)+1)), col='red')
    abline(v=log10(350+1), col='red')
})

## QC 
t0 <- 5000
t1 <- 100
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
#p2 <- Pagoda2$new(x = countMatrix, n.cores = 4, trim=25, batch=countMatrix.batch)

# Adjust the variance
p2$adjustVariance(plot = T, gam.k = 10)

# Calculate a PCA reduction with the number of PCs specified by nPCs
# and using only the n.odgenes overdispersed genes -- in this case 2000
p2$calculatePcaReduction(nPcs = 100, n.odgenes = 2.e3, maxit=1000)

# Generate K-nearest neighbour graph
p2$makeKnnGraph(k = 50, type = 'PCA', center = T, weight.type = 'none', n.cores = 4, distance = 'cosine')

# Generate an embedding with tSNE on the basis of the PCA reduction
#p2$getEmbedding(type = 'PCA', embeddingType = 'tSNE', verbose=TRUE, perplexity=50)
p2$getEmbedding(type='PCA',M=5,perplexity=50,gamma=1)

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
#                 embeddingType = 'tSNE',
                 mark.clusters = T,
                 main="infomap",
                  min.group.size=50)

p2$plotEmbedding(type = 'PCA',
                 clusterType = 'multilevel',
#                 embeddingType = 'tSNE',
                 mark.clusters = T,
                 main="multilevel",
                 min.group.size=50)

p2$plotEmbedding(type = 'PCA',
                 clusterType = 'walktrap',
#                 embeddingType = 'tSNE',
                 mark.clusters = T,
                 main="walktrap",
                 min.group.size=50)

cl0 <- p2$clusters$PCA$walktrap
cl0[cl0 %in% which(table(cl0)<50)] <- NA
cl0[cl0 %in% c(7,4)] <- NA
retest <- cl0
table(retest)

## Final plot
par(mfrow=c(1,4))
## Color by time point
p2$plotEmbedding(type = 'PCA',
                 groups=countMatrix.batch,
#                 embeddingType = 'tSNE',
                 mark.clusters = T,
                 main="Batch")
p2$plotEmbedding(type = 'PCA',
                 groups=annot,
#                 embeddingType = 'tSNE',
                 mark.clusters = T,
                 main="Annot")
p2$plotEmbedding(type = 'PCA',
                 groups=retest,
#                 embeddingType = 'tSNE',
                 mark.clusters = T,
                 main="New Annot")
p2$plotEmbedding(type = 'PCA',
                 colors=p2$depth,
#                 embeddingType = 'tSNE',
                 mark.clusters = T,
                 main="Depth")

############## Now retest
good.cells <- na.omit(retest)
length(good.cells)
dim(countMatrix)

Cer <- Pagoda2$new(x = countMatrix[, names(good.cells)], n.cores = 4, trim=10, batch=countMatrix.batch[names(good.cells)])

# Adjust the variance
Cer$adjustVariance(plot = T, gam.k = 10)

# Calculate a PCA reduction with the number of PCs specified by nPCs
# and using only the n.odgenes overdispersed genes -- in this case 2000
Cer$calculatePcaReduction(nPcs = 150, n.odgenes = 2.e3, maxit=1000)

# Generate K-nearest neighbour graph
Cer$makeKnnGraph(k = 15, type = 'PCA', center = T, weight.type = 'none', n.cores = 4, distance = 'cosine')

# Generate an embedding with tSNE on the basis of the PCA reduction
#Cer$getEmbedding(type = 'PCA', embeddingType = 'tSNE', verbose=TRUE, perplexity=30)
Cer$getEmbedding(type='PCA',M=5,perplexity=50,gamma=1)

# Identify clusters using the infomap.community method
# on the basis of the reduction called 'PCA' (generated above)
# Save the resulting clustering as 'infomap'
Cer$getKnnClusters(method = infomap.community, type = 'PCA', name = 'infomap')

# Do an independent identification of clusters using the
# multilevel community algorithm again using the PCA reduction
# and save it as 'multilevel'. This does not overwrite the
# previous clustering
Cer$getKnnClusters(method = multilevel.community, type = 'PCA', name='multilevel')

# Do yet another clustering
Cer$getKnnClusters(method = walktrap.community, type = 'PCA', name='walktrap')

# Plot
par(mfrow=c(1,3))
Cer$plotEmbedding(type = 'PCA',
                 clusterType = 'infomap',
#                 embeddingType = 'tSNE',
                 mark.clusters = T,
                 main="infomap",
                 min.group.size=30)

Cer$plotEmbedding(type = 'PCA',
                 clusterType = 'multilevel',
#                 embeddingType = 'tSNE',
                 mark.clusters = T,
                 main="multilevel",
                 min.group.size=30)

Cer$plotEmbedding(type = 'PCA',
                 clusterType = 'walktrap',
#                 embeddingType = 'tSNE',
                 mark.clusters = T,
                 main="walktrap",
                 min.group.size=30)

newannot <- Cer$clusters$PCA$infomap
newannot[newannot %in% which(table(newannot)<30)] <- NA

## Final plot
par(mfrow=c(1,4))
# Color by time point
Cer$plotEmbedding(type = 'PCA',
                 groups=countMatrix.batch,
#                 embeddingType = 'tSNE',
                 mark.clusters = T,
                 main="Batch",
                 mark.cluster.cex = 1,
                 alpha=0.1)
Cer$plotEmbedding(type = 'PCA',
                 groups=annot,
#                 embeddingType = 'tSNE',
                 mark.clusters = T,
                 main="Annot",
                 mark.cluster.cex = 1,
                 alpha=0.1)
Cer$plotEmbedding(type = 'PCA',
                 groups=newannot,
#                 embeddingType = 'tSNE',
                 mark.clusters = T,
                 main="multilevel",
                 mark.cluster.cex = 1,
                 alpha=0.1)
Cer$plotEmbedding(type = 'PCA',
                 colors=Cer$depth,
#                 embeddingType = 'tSNE',
                 mark.clusters = T,
                 main="Depth",
                 mark.cluster.cex = 1,
                 alpha=0.1)


#save.image('Cer.06092017.RData')  

Cer.annot <- as.character(newannot)
names(Cer.annot) <- names(newannot)
Cer.annot[newannot == 1] <- 'Gran'
Cer.annot[newannot == 2] <- 'PurkA'
Cer.annot[newannot == 4] <- 'PurkB'
Cer.annot[newannot == 6] <- 'Oli'
Cer.annot[newannot == 5] <- 'Opc_cerB'
Cer.annot[newannot == 8] <- 'Opc_cerA'
Cer.annot[newannot == 9] <- 'Ast_cerA'
Cer.annot[newannot == 3] <- 'Ast_cerB'
Cer.annot[newannot == 10] <- 'Per'
Cer.annot[newannot == 11] <- 'End'
Cer.annot[newannot == 12] <- 'End'
Cer.annot[newannot == 7] <- 'Mic'
Cer.annot <- factor(Cer.annot)
head(Cer.annot)

Cer$plotEmbedding(type = 'PCA',
                  groups=Cer.annot,
#                 embeddingType = 'tSNE',
                 mark.clusters = T,
                 mark.cluster.cex = 1,
                 alpha=0.1)

save(Cer, Cer.annot, file='Cer_results_06092017-2.RData')



load('Cer_results_06092017-2.RData')
library("largeVis", lib.loc="/usr/local/lib/R/site-library")
library(pagoda2)
load('final.annot.RData')

final.cer <- final[rownames(Cer$counts)]
final.cer[grepl('In', final.cer)] <- NA
final.cer[grepl('Ex', final.cer)] <- NA
table(final.cer)


Cer$plotEmbedding(type = 'PCA',
                  groups=final.cer,
#                  embeddingType = 'tSNE',
                  mark.clusters = T,
                  main="Final Annot",
                  mark.cluster.cex = 1,
                  alpha=0.1)
set.seed(3)
Cer$plotEmbedding(type = 'PCA',
                  groups=Cer.annot,
                  ##embeddingType = 'tSNE',
                  mark.clusters = T,
                  main="Final Annot",
                  mark.cluster.cex = 1,
                  alpha=0.1, shuffle.colors=TRUE)

Cer.annot.final <- as.character(Cer.annot)
Cer.annot.final[Cer.annot=='Ast_cerA'] <- 'Ast'
Cer.annot.final[Cer.annot=='Ast_cerB'] <- 'Ast_cer'
Cer.annot.final[Cer.annot=='Opc_cerA'] <- 'Opc'
Cer.annot.final[Cer.annot=='Opc_cerB'] <- 'Opc_cer'
names(Cer.annot.final) <- names(Cer.annot)
Cer.annot.final <- factor(Cer.annot.final)

Cer$getEmbedding(type = 'PCA', embeddingType = 'tSNE', verbose=TRUE, perplexity=30)
Cer$plotEmbedding(type = 'PCA',
                  groups=Cer.annot.final,
                  embeddingType = 'tSNE',
                  mark.clusters = T,
                  main="Final Annot",
                  mark.cluster.cex = 1,
                  alpha=0.1, shuffle.colors=TRUE)
save(Cer.annot.final, Cer, file="Cer.annot.final.RData")

load("Cer.annot.final.RData")

