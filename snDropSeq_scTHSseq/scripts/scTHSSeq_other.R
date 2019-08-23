####################
## Chromoclust analysis on all epigenetic data
## for non-Occ datasets projected onto axis established by Occ data
##
## @date 4/7/2017
## @author: Jean Fan
## @email: jeanfan@fas.harvard.edu
###################

## load Occ results
load("hR10_vc_comb_fixed_withopcend.RData")
hR10_vc_r <- r
hR10_vc_annot <- annot

## load non-Occ data for projection
load("hR7_fc9_cov_peaks_spp_RMMMcombined_RMrepeatmask100_bandwidth500_step100_thr5_span10_fdr1e-07.RData")
#load("hR11_cer_cov_peaks_spp_RMMMcombined_RMrepeatmask100_bandwidth500_step100_thr5_span10_fdr1e-07.RData")
#load("hR11_hip_cov_peaks_spp_RMMMcombined_RMrepeatmask100_bandwidth500_step100_thr5_span10_fdr1e-07.RData")
library(Matrix)
## convert by slice
cov <- do.call(cbind,apply(embed(pmin(ncol(cov),seq(0,ceiling(ncol(cov)/1e3))*1e3),2),1,function(ii) {
    Matrix(cov[,seq(ii[2]+1,ii[1])],sparse=T)
}))
cov.init <- cov
#sites <- intersect(rownames(cov), colnames(hR10_vc_r$counts))
#cov <- cov[sites,]
dim(cov)

quantile(Matrix::rowSums(cov>0))
vi <- Matrix::rowSums(cov>0)>10; table(vi)
cov <- cov[vi, ]
dim(cov)

quantile(Matrix::colSums(cov>0))
vi <- Matrix::colSums(cov>0)>10; table(vi)
cov <- cov[, vi]
dim(cov)

## binarize
cov[cov>0] <- 1


############################# chromoclust
library("Rcpp", lib.loc="/usr/local/lib/R/site-library")
library("dbscan", lib.loc="/usr/local/lib/R/site-library")
library("largeVis", lib.loc="/usr/local/lib/R/site-library")
library(Cairo)
library(parallel)
library(WGCNA)
library(Matrix)
source("/home/pkharchenko/m/p2/schromoclust.r")

#batch <- as.factor(gsub(".*\\|(.*)","\\1",colnames(cov))); names(batch) <- colnames(cov)
batch <- sapply(colnames(cov), function(x) strsplit(strsplit(x, '.unique.')[[1]][2], '_')[[1]][1])
r <- Pagoda2$new(cov,trim=0,n.cores=20,batch=batch,min.cells.per.gene=0)
#r <- Pagoda2$new(cov,trim=0,n.cores=20)
x <- r$getRefinedLibSizes(verbose=T)

#pdf('test.pdf')
#par(mfrow=c(2,2), mar = c(3.5,3.5,2.0,0.5), mgp = c(2,0.65,0), cex = 1.0);
#plot(rowSums(r$misc$rawCounts),x$depth,col=adjustcolor(1,alpha=0.1),cex=0.3); abline(a=0,b=1,lty=2,col=2)
#p <- colMeans(r$misc$rawCounts/rowSums(r$misc$rawCounts)); p <- p/sum(p);
#plot(p,x$p[1,],col=adjustcolor(1,alpha=0.1),cex=0.3); abline(a=0,b=1,lty=2,col=2)
#x <- r$adjustVariance(plot=T,do.par=F,gam.k=20,rescale.mat=T)
x <- r$carefullyNormalizedReduction(sparse=T)
#dev.off()

## Project
opc <- hR10_vc_r$misc$PCA
opc$nv <- opc$v[match(colnames(r$reductions$normalized),rownames(opc$v)),]
opc$nv[is.na(opc$nv)] <- 0
r$reductions$PCA <- as.matrix(r$reductions$normalized %*% opc$nv)
#r$calculatePcaReduction(nPcs=30,type='normalized')

## predict
lvec <- colSumByFac(hR10_vc_r$misc[['rawCounts']], hR10_vc_annot)[-1,] + 1
lvec <- t(lvec/pmax(1,rowSums(lvec)))
colnames(lvec) <- levels(hR10_vc_annot)
rownames(lvec) <- colnames(hR10_vc_r$misc[['rawCounts']])
str(lvec)
ld <- jsDist(lvec); colnames(ld) <- rownames(ld) <- colnames(lvec)
colnames(lvec)
length(levels(hR10_vc_annot))
levels(hR10_vc_annot)

int <- intersect(rownames(lvec), colnames(r$misc[['rawCounts']]))
lvec <- lvec[int,]
ll <- as.matrix(t(r$misc[['rawCounts']][, int] %*% log(lvec))) - colSums(lvec)
ll <- t(ll)-apply(ll,2,max)+10; ll <- exp(ll); ll <- ll/rowSums(ll)
cl1 <- as.factor(colnames(ll)[apply(ll,1,which.max)]);
names(cl1) <- rownames(r$counts)
pred.final <- cl1

## check correlation with sample factor
#apply(r$reductions$PCA,2,cor,as.numeric(r$batch))[1:10]
## check correlation with depth .. we expect to have some
apply(r$reductions$PCA,2,cor,r$depth)[1:10]
r$makeKnnGraph(k=30,type='PCA',center=T,distance='cosine');
r$getKnnClusters(method=multilevel.community,type='PCA',name='multilevel')
r$getKnnClusters(method=walktrap.community,type='PCA',name='walktrap')
r$getKnnClusters(method=infomap.community,type='PCA',name='infomap')
r$getEmbedding(type='PCA',M=5,perplexity=10,gamma=0.9)

## plot
par(mfrow=c(2,2));
r$plotEmbedding(type='PCA',show.legend=F,mark.clusters=T,group=pred.final, mark.cluster.cex=1,alpha=0.1, main="prediction")
r$plotEmbedding(type='PCA',show.legend=F,clusterType='multilevel',mark.clusters=T,mark.cluster.cex=1,alpha=0.1, main="multilevel")
r$plotEmbedding(type='PCA',show.legend=F,clusterType='walktrap',mark.clusters=T,mark.cluster.cex=1,alpha=0.1, main="walktrap")
r$plotEmbedding(type='PCA',show.legend=F,clusterType='infomap',mark.clusters=T,mark.cluster.cex=1,alpha=0.1, main="infomap")
#r$plotEmbedding(type='PCA',show.legend=F,mark.clusters=T,mark.cluster.cex=1,colors=r$depth,alpha=0.1,main="depth")
#r$plotEmbedding(type='PCA',show.legend=F,mark.clusters=T,mark.cluster.cex=1,group=r$batch,alpha=0.1,main="batch")

#pdf('tsne.pdf')
## tSNE embedding
#r$getEmbedding(type='PCA',perplexity=30,verbose=T,embeddingType='tSNE')
r$getEmbedding(type='PCA',perplexity=50,verbose=T,embeddingType='tSNE')
#r$getEmbedding(type='PCA',perplexity=60,verbose=T,embeddingType='tSNE')

par(mfrow=c(2,2));
r$plotEmbedding(type='PCA',show.legend=F,mark.clusters=T,group=pred.final, mark.cluster.cex=1,alpha=0.1, main="prediction",embeddingType='tSNE')
r$plotEmbedding(type='PCA',show.legend=F,clusterType='multilevel',mark.clusters=T,mark.cluster.cex=1,alpha=0.1,embeddingType='tSNE', main="multilevel")
r$plotEmbedding(type='PCA',show.legend=F,clusterType='walktrap',mark.clusters=T,mark.cluster.cex=1,alpha=0.1,embeddingType='tSNE', main="walktrap")
r$plotEmbedding(type='PCA',show.legend=F,clusterType='infomap',mark.clusters=T,mark.cluster.cex=1,alpha=0.1,embeddingType='tSNE', main="infomap", min.group.size=50)
#r$plotEmbedding(type='PCA',show.legend=F,mark.clusters=T,mark.cluster.cex=1,colors=r$depth,alpha=0.1,main="depth",embeddingType='tSNE')
#r$plotEmbedding(type='PCA',show.legend=F,mark.clusters=T,mark.cluster.cex=1,group=r$batch,alpha=0.1,main="batch",embeddingType='tSNE')
#r$plotEmbedding(type='PCA',show.legend=F,mark.clusters=T,mark.cluster.cex=1,group=sample,alpha=0.1,main="sample",embeddingType='tSNE')
#dev.off()

## consolidate groupings
#cl0 <- r$clusters$PCA$multilevel[rownames(r$counts)]
cl0 <- r$clusters$PCA$infomap[rownames(r$counts)]
min.cluster.size = 50
cl0[cl0 %in% which(table(cl0)<min.cluster.size)] <- NA

cl0[cl0==6] <- NA
cl0[cl0==7] <- NA
cl0[cl0==4] <- NA

group <- as.character(cl0)
group[cl0==1] <- 'Oli'
group[cl0==5] <- 'Ast'
group[cl0==2] <- 'Ex'
group[cl0==3] <- 'In'
#group[cl0==8] <- 'Mic'
group[cl0==8] <- NA
group <- factor(group)
names(group) <- names(cl0)

#pdf('plots.pdf')
par(mfrow=c(2,2))
#r$plotEmbedding(type='PCA',groups=cl0, mark.clusters=TRUE, alpha=0.1)
#r$plotEmbedding(type='PCA',groups=cl0, mark.clusters=TRUE, alpha=0.1,embeddingType='tSNE', main="tSNE")
r$plotEmbedding(type='PCA',show.legend=F,mark.clusters=T,group=pred.final, mark.cluster.cex=1,alpha=0.1, main="prediction",embeddingType='tSNE')
r$plotEmbedding(type='PCA',groups=group, mark.clusters=TRUE, alpha=0.1,embeddingType='tSNE', main="groups")
r$plotEmbedding(type='PCA',show.legend=F,mark.clusters=T,mark.cluster.cex=1,colors=r$depth,alpha=0.1,main="depth",embeddingType='tSNE')
r$plotEmbedding(type='PCA',show.legend=F,mark.clusters=T,mark.cluster.cex=1,group=r$batch,alpha=0.1,main="batch",embeddingType='tSNE')
#dev.off()

lvec <- colSumByFac(r$misc[['rawCounts']],group)[-1,] + 1
#lvec <- t(lvec)
lvec <- t(lvec/pmax(1,rowSums(lvec))) 
#lvec <- t(lvec/pmax(1,rowSums(lvec>0))) 
#lvec <- t(lvec/colSums(lvec>0))
#lvec <- t(lvec/colSums(lvec))
#colnames(lvec) <- paste('cl',which(table(cl0)>0),sep='')
colnames(lvec) <- levels(group)
rownames(lvec) <- colnames(r$misc[['rawCounts']])
str(lvec)
ld <- jsDist(lvec); colnames(ld) <- rownames(ld) <- colnames(lvec)
#hctree <- stats::hclust(as.dist(ld),method='ward')
hctree <- stats::hclust(as.dist(ld))
plot(hctree,axes=F,sub="",ylab="",cex=0.8,main='initial cluster dendrogram')

## calculate cluster assignment probabilities for each cell
#ll <- as.matrix(t(r$misc[['rawCounts']] %*% log(lvec))) - colSums(lvec)
#ll <- t(ll)-apply(ll,2,max)+10; ll <- exp(ll); ll <- ll/rowSums(ll)
#cl1 <- as.factor(apply(ll,1,which.max)); 
#names(cl1) <- rownames(r$counts)
#table(cl1==cl0)
#cl1[which(cl1!=cl0)]
#cl0[which(cl0!=cl1)]


############ differential peaks

dg <- r$getDifferentialGenes(upregulated.only = TRUE,groups=group)
par(mfrow=c(3,4))
lapply(1:length(levels(group)), function(i) {
    r$plotEmbedding(type='PCA',embedding='tSNE',show.legend=F,colors=rowMeans(r$counts[,names(dg[[i]])]),alpha=0.1, main=levels(group)[i])
})

save(r, hctree, group, dg, file="results/hR7_fc9/hR7_fc9_results_04072017.RData")

load("results/hR7_fc9/hR7_fc9_results_04072017.RData")
embed <- get('embeddings', slot(r,'.xData'))[['PCA']][['tSNE']]
data <- cbind(embed, as.character(group[rownames(embed)]))
head(data)
colnames(data) <- c('tSNE1', 'tSNE2', 'group')
plot(data[,1], data[,2], col=data[,3])
write.table(data, file="results/hR7_fc9/hR7_fc9_forbrandon_data.txt", quote=FALSE, sep="\t")
### write bed
direct.dg <- dg
lapply(names(direct.dg),function(i) {
        peaks <- names(direct.dg[[i]]);
        peak.df <- data.frame(do.call(rbind,strsplit(peaks,":|-")),stringsAsFactors=F); colnames(peak.df) <- c('chr','start','end');
        peak.df$name <- paste('peak',seq(1,nrow(peak.df)),sep='');
        peak.df$score <- round(pmin(10,as.numeric(direct.dg[[i]]))*1e2)
        fname <- paste0('results/hR7_fc9/hR7_fc9_forbrandon_direct.',i,'.diffPeaks.bed',sep='')
        write(paste('track name="Direct',i,'" description="Direct ',i,' differentially accessible peaks (upregulated only)" visibility=2 useScore=1',sep=''),file=fname) # header
        write.table(peak.df,append=T,file=fname,quote=F,row.names=F,col.names=F,sep='\t')
})
