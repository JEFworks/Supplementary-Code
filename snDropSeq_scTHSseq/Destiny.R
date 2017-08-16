###################
## Diffusion maps to study Opc and Oli trajectory
## @author: Jean Fan
## @email: jeanfan@fas.harvard.edu
###################

#source("https://bioconductor.org/biocLite.R")
#biocLite("destiny")
library(destiny)
library(Biobase)
library("Rcpp", lib.loc="/usr/local/lib/R/site-library")
library("dbscan", lib.loc="/usr/local/lib/R/site-library")
library("largeVis", lib.loc="/usr/local/lib/R/site-library")
library(Cairo)
library(parallel)
library(WGCNA)
library(Matrix)
library("pagoda2", lib.loc="/home/barkasn/R/x86_64-pc-linux-gnu-library/3.3/")
source("/home/pkharchenko/m/p2/schromoclust.r") 

## Load snDrop-seq Results
load("~/Projects/Kun_snDropSeq/Rcomb/fin.blue.06152017.RData")
Occ$plotEmbedding(type='PCA',show.legend=F,mark.clusters=T,groups=fin.blue, mark.cluster.cex=1,alpha=0.1, main="prediction",embeddingType='tSNE')

## Select Occ OPCs and Olis
cells <- names(fin.blue)[which(fin.blue %in% c('OPC', 'OPC_Cer', 'Oli'))]
vi <- grepl('occ', p2$batch[cells])
cells <- cells[vi]
cells.annot <- factor(fin.blue[cells])
levels(cells.annot)
Occ$plotEmbedding(type='PCA',show.legend=F,mark.clusters=T,groups=cells.annot, mark.cluster.cex=1,alpha=0.1, main="Occ used cells",embeddingType='tSNE')

## Get normalized counts
raw_ct <- p2$counts[cells,]
annot_ct <- factor(fin.blue[cells])
dim(raw_ct)
length(annot_ct)
raw_ct <- as.data.frame(as.matrix(raw_ct))
head(raw_ct)
ct <- as.ExpressionSet(raw_ct)
ct$annot <- annot_ct

## Run Destiny
dm <- DiffusionMap(ct, k=100)
#save(dm, file="dm_allgenes_OccOpcOliOnly.RData")

## Plot trajectory
library(RColorBrewer)
par(mfrow=c(2,2))
plot(dm, pch = 20, col = rainbow(2)[ct$annot])
plot(dm, 1:2, pch = 20, col = rainbow(2)[ct$annot])
colcol <- rainbow(length(levels(p2$batch)))
names(colcol) <- levels(p2$batch)
colcol <- colcol[p2$batch[cells]]
plot(dm, pch = 20, col = colcol)
plot(dm, 1:2, pch = 20, col = colcol)

## Order cells by first eigenvector
slotNames(dm)
head(slot(dm, 'eigenvectors'))
ev <- slot(dm, 'eigenvectors')
cell.order <- order(ev[,1])
names(cell.order) <- cells

## Identify differentially upregulated genes
selected <- rep(NA, length(cell.order))
names(selected) <- names(cell.order)
selected <- selected[cell.order]
selected[1:400] <- 1 ## OPCs
selected[700:1100] <- 3 ## Immature Oli
selected[2664:3064] <- 4 ## Mature Oli
selected <- factor(selected)

selected.annot <- rep(NA, nrow(Occ$counts))
names(selected.annot) <- rownames(Occ$counts)
selected.annot[names(selected)] <- selected
table(selected.annot)
selected.annot <- factor(selected.annot)
direct.dg <- Occ$getDifferentialGenes(upregulated.only = TRUE, groups=selected.annot, z.threshold=1.96)

direct.dg.genes <- lapply(direct.dg, function(x) {
    rownames(x[x$highest,])
})
head(direct.dg.genes)

mat <- t(raw_ct[cell.order,unlist(direct.dg.genes)])
mat[is.na(mat)] <- 0
mat <- t(scale(t(mat)))
range(mat)
mat[mat < -1.5] <- -1.5
mat[mat > 1.5] <- 1.5
heatmap(mat, col=colorRampPalette(c('blue', 'white', 'red'))(100), Rowv=NA, Colv=NA, scale="none", ColSideColors=rainbow(5)[selected])
heatmap(mat, col=colorRampPalette(c('blue', 'white', 'red'))(100), Rowv=NA, Colv=NA, scale="none", ColSideColors=rainbow(2)[annot_ct[colnames(mat)]])

## Visualize select known markers
blue.genes <- c('MOBP', 'MBP', 'PLP1', 'VCAN', 'SOX6', 'LUZP2', 'MEG3', 'SNTG1', 'GRIK2', 'FGF14', 'RORA')
blue.genes <- intersect(rownames(mat), blue.genes )
heatmap(mat[blue.genes,], col=colorRampPalette(c('blue', 'white', 'red'))(100), Rowv=NA, Colv=NA, scale="none", ColSideColors=rainbow(5)[selected])


################ Use genes to order cells by accessibility
load('~/Projects/Kun_Epigenetics/RJuly2017/Occ/Occ_refine_results_07062016-2.RData')
acc.r <- Occ.r
acc.annot <- Occ.refine
acc.cells <- acc.annot %in% c('Opc', 'Oli'); table(acc.cells)
acc.cells <- names(acc.annot)[acc.cells]
head(acc.cells)
table(acc.annot[acc.cells])

acc.cells.annot <- factor(acc.annot[acc.cells])
acc.r$plotEmbedding(type='PCA',show.legend=F,mark.clusters=F,groups=acc.cells.annot, mark.cluster.cex=1,alpha=0.1, main="Occ used cells",embeddingType='tSNE')

sites <- colnames(acc.r$counts)
acc.mat <- t(acc.r$misc$rawCounts)[, acc.cells]
acc.annot <- factor(acc.annot[acc.cells])
levels(acc.annot)

## Annotate peaks with genes
peaks <- sites
peak.df <- data.frame(do.call(rbind,strsplit(peaks,":|-")),stringsAsFactors=F); colnames(peak.df) <- c('chr','start','end');
peak.df$mid <- (as.numeric(peak.df$end)+as.numeric(peak.df$start))/2
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
seqlevels(txdb)
library(ChIPseeker)
peak.gr <- with(peak.df, GRanges(as.character(chr), IRanges(as.numeric(as.character(start)), as.numeric(as.character(end)))))
p2g <- annotatePeak(peak.gr, TxDb=txdb)

p2g.detail <- cbind(geneId=p2g@anno$geneId, distanceToTSS=p2g@anno$distanceToTSS, p2g@detailGenomicAnnotation)
rownames(p2g.detail) <- peaks

require(annotate)
library(org.Hs.eg.db)
p2g.detail$symbol <- getSYMBOL(p2g.detail$geneId, data='org.Hs.eg')
head(p2g.detail)

genes <- c(
    direct.dg.genes[[1]],
    direct.dg.genes[[2]],
    direct.dg.genes[[3]]
)
genes.annot <- c(
    rep('Opc genes', length(direct.dg.genes[[1]])),
    rep('Immature Oli genes', length(direct.dg.genes[[2]])),
    rep('Mature Oli genes', length(direct.dg.genes[[3]]))
)
names(genes.annot) <- genes
genes.annot <- factor(genes.annot)

mat <- t(raw_ct[cell.order,genes])
mat[is.na(mat)] <- 0
mat <- t(scale(t(mat)))
range(mat)
mat[mat < -1.5] <- -1.5
mat[mat > 1.5] <- 1.5
heatmap(mat, col=colorRampPalette(c('blue', 'white', 'red'))(100), Rowv=NA, Colv=NA, scale="none", ColSideColors=rainbow(5)[selected])

## Calculate joint accessibility for each group of genes
csdf <- do.call(rbind, lapply(levels(genes.annot), function(ga) {
    gt <- names(genes.annot)[genes.annot==ga]
    vi <- p2g.detail$symbol %in% gt; table(vi)

    p2g.detail.subset <- p2g.detail[vi,]
    vi <- !p2g.detail.subset$distal_intergenic; table(vi)
    p2g.detail.subset <- p2g.detail.subset[vi,]
    
    sites <- na.omit(rownames(p2g.detail.subset))
    print(length(sites))
    x = as.character(sites)
    if(length(sites)>1) {
        ## tpm
        data.use = t(colSums(acc.mat[sites,]))/colSums(acc.mat)*1e6/length(sites)
    } else {
        data.use = t(acc.mat[sites,])
    }
    data.use
}))
rownames(csdf) <- levels(genes.annot)

## Cluster cells
csdf.scaled <- csdf
csdf.scaled <- t(scale(t(csdf.scaled))) ## row scale
csdf.scaled <- scale(csdf.scaled)
csdf.scaled <- t(scale(t(csdf.scaled))) ## row scale
csdf.scaled[is.nan(csdf.scaled)] <- 0
csdf.scaled[is.na(csdf.scaled)] <- 0
range(csdf.scaled)
csdf.scaled[csdf.scaled < -2] <- -2
csdf.scaled[csdf.scaled > 2] <- 2
bar <- t(csdf.scaled)
bar <- cbind(bar, as.numeric(acc.annot[rownames(bar)]))
hc <- hclust(dist(bar), method='ward.D2')

## Visualize
library(gplots)
dend <- as.dendrogram(hc)
temp <- dend[[2]][[1]]
dend[[2]][[1]] <- dend[[2]][[2]]
dend[[2]][[2]] <- temp
heatmap.2(csdf.scaled[c(3,1,2),], Rowv=NA, Colv=dend, scale="none", col=colorRampPalette(c('blue', 'white', 'red'))(100), ColSideColors=rainbow(2)[acc.annot[hc$labels]], trace="none", xlab=NA,margins=c(5,20),labCol=FALSE)
heatmap.2(csdf[c(3,1,2),], Rowv=NA, Colv=dend, scale="row", col=colorRampPalette(c('blue', 'white', 'red'))(100), ColSideColors=rainbow(2)[acc.annot[hc$labels]], trace="none", xlab=NA,margins=c(5,20),labCol=FALSE)

## Double check
new.groups <- as.factor(cutree(as.hclust(dend),3))
heatmap.2(csdf[c(3,1,2),], Rowv=NA, Colv=dend, scale="row", col=colorRampPalette(c('blue', 'white', 'red'))(100), ColSideColors=rainbow(3)[new.groups], trace="none", xlab=NA,margins=c(5,20),labCol=FALSE)
csdf <- do.call(rbind, lapply(levels(genes.annot), function(ga) {
    gt <- names(genes.annot)[genes.annot==ga]
    vi <- p2g.detail$symbol %in% gt; table(vi)

    p2g.detail.subset <- p2g.detail[vi,]
    #vi <- p2g.detail.subset$Promoter; table(vi)
    #p2g.detail.subset <- p2g.detail.subset[vi,]
    vi <- !p2g.detail.subset$distal_intergenic; table(vi)
    p2g.detail.subset <- p2g.detail.subset[vi,]
    
    sites <- na.omit(rownames(p2g.detail.subset))
    print(length(sites))
    x = as.character(sites)
    if(length(sites)>1) {
        ## tpm
        data.use = t(colSums(acc.mat[sites,]))/colSums(acc.mat)*1e6/length(sites)
    } else {
        data.use = t(acc.mat[sites,])
    }
    #data.use

    ident.use = na.omit(new.groups[colnames(data.use)])

    cs <- unlist(lapply(levels(ident.use), function(g) {
        ch <- na.omit(names(ident.use)[ident.use==g])
        sum(data.use[, ch])/length(ch) 
    }))
    names(cs) <- levels(ident.use)
    
    cs
    
}))
csdf

## for each gene
csdf <- do.call(rbind, lapply(1:length(genes), function(i) {
    gene.name <- genes[i]
    type.name <- names(genes)[i]
    print(gene.name)
    vi <- which(p2g.detail$symbol==gene.name)
    p2g.detail.subset <- p2g.detail[vi,]
    
    #vi <- which(p2g.detail.subset$Promoter)
    #p2g.detail.subset <- p2g.detail.subset[vi,]
    vi <- !p2g.detail.subset$distal_intergenic; table(vi)
    p2g.detail.subset <- p2g.detail.subset[vi,]

    sites <- na.omit(rownames(p2g.detail.subset))
    print(length(sites))
    x = as.character(sites)
    
    if(length(sites)==0) {
        return(NA)
    }
    if(length(sites)>1) {
        ## tpm
        data.use = t(colSums(acc.mat[sites,]))/colSums(acc.mat)*1e6/length(sites)
    } else {
        data.use = t(acc.mat[sites,])
    }
    data.use
}))
rownames(csdf) <- genes

mat <- csdf[, hc$labels]
mat <- na.omit(mat)
mat <- log10(mat+1)
#mat <- t(scale(t(mat)))
mat[mat > 0] <- 1
heatmap.2(mat, Rowv=NA, Colv=dend, scale="none", col=colorRampPalette(c('white', 'black'))(100), ColSideColors=rainbow(2)[acc.annot[hc$labels]], trace="none", xlab=NA,margins=c(5,20),labCol=FALSE, labRow=FALSE, RowSideColors=rainbow(3)[genes.annot[rownames(mat)]])

save.image('destiny.RData')

save(direct.dg, dm, cell.order, selected, genes, genes.annot, new.groups, file="destiny_forBlue.RData")
