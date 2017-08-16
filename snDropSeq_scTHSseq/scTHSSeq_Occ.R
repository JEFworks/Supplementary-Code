####################
# Chromoclust analysis on all epigenetic data
# See spp_comb.R for previous work to get coverage matrix
#
# visual cortex only
#
# @date 2/28/2017
####################

## load data
#load("hR10_vc_cov_peaks_spp_RMMM_RMrepeatmask100_bandwidth500_step100_thr5_span10_fdr1e-07.RData")
load("hR10_vc_cov_peaks_spp_RMMMcombined_RMrepeatmask100_bandwidth500_step100_thr5_span10_fdr1e-07.RData")
bams <- colnames(cov)
dim(cov)

library(Matrix)
#cov <- as(cov, "sparseMatrix")
## convert by slice
cov <- do.call(cbind,apply(embed(pmin(ncol(cov),seq(0,ceiling(ncol(cov)/1e3))*1e3),2),1,function(ii) {
    Matrix(cov[,seq(ii[2]+1,ii[1])],sparse=T)
}))
dim(cov)

## filter for common only
quantile(rowSums(cov>0))
vi <- rowSums(cov>0)>0; table(vi)
cov <- cov[vi,]
dim(cov)

## seems like there are some bad cells causing issues
quantile(Matrix::colSums(cov>0))
#vi <- Matrix::colSums(cov>0)>100; table(vi)
#cov <- cov[, vi]

## binarize
cov[cov>0] <- 1

## clean up name
cn <- colnames(cov)
cn <- gsub('../data-raw5/hR10_vc_noAlt/ubam_rm_mm/', '', cn)
cn <- gsub('.unique.', '|', cn)
cn <- gsub('.rm_mm.bam', '', cn)
cn <- gsub('.hR10_vc_noAlt', '', cn)
cn <- paste0('hR10_vc_', cn)
head(cn)
colnames(cov) <- cn
names(bams) <- cn

## Peter's annotations
rold <- readRDS("/home/pkharchenko/kun/ths/combined/vc.p2.rds")
cell.groups.peter <- as.factor(rold$clusters$PCA$walktrap);

gc()

############################# Chromoclust
library("Rcpp", lib.loc="/usr/local/lib/R/site-library")
library("dbscan", lib.loc="/usr/local/lib/R/site-library")
library("largeVis", lib.loc="/usr/local/lib/R/site-library")
library(Cairo)
library(parallel)
library(WGCNA)
library(Matrix)
source("/home/pkharchenko/m/p2/schromoclust.r")

batch <- as.factor(gsub(".*\\|(.*)","\\1",colnames(cov))); names(batch) <- colnames(cov)
r <- Pagoda2$new(cov,trim=0,n.cores=20,batch=batch)
#x <- r$getRefinedLibSizes(verbose=T)

#pdf('test.pdf')
#par(mfrow=c(2,2), mar = c(3.5,3.5,2.0,0.5), mgp = c(2,0.65,0), cex = 1.0);
#plot(rowSums(r$misc$rawCounts),x$depth,col=adjustcolor(1,alpha=0.1),cex=0.3); abline(a=0,b=1,lty=2,col=2)
#p <- colMeans(r$misc$rawCounts/rowSums(r$misc$rawCounts)); p <- p/sum(p);
#plot(p,x$p[1,],col=adjustcolor(1,alpha=0.1),cex=0.3); abline(a=0,b=1,lty=2,col=2)
#x <- r$adjustVariance(plot=T,do.par=F,gam.k=20,rescale.mat=T)
x <- r$carefullyNormalizedReduction(sparse=T)
#dev.off()

#pdf('test.pdf')
r$calculatePcaReduction(nPcs=30,type='normalized')
## show first 10 components
#par(mfrow=c(3,3));
#invisible(apply(matrix(1:10,ncol=2,byrow=T),1,function(ii) {
#   r$embeddings$PCA$two <- r$reductions$PCA[,c(ii[1],ii[2])]
#   r$plotEmbedding(type='PCA',show.legend=F,groups=cell.groups.peter,embeddingType='two',mark.clusters=T,main=pas#te(ii,collapse=" vs. "),cex=0.4)
#}))
#dev.off()

## check correlation with sample factor
#apply(r$reductions$PCA,2,cor,as.numeric(r$batch))[1:10]
## check correlation with depth .. we expect to have some
#apply(r$reductions$PCA,2,cor,r$depth)[1:10]
r$makeKnnGraph(k=10,type='PCA',center=T,distance='cosine');
#r$makeKnnGraph(k=50,type='PCA',center=T,distance='L2');
r$getKnnClusters(method=multilevel.community,type='PCA',name='multilevel')
#r$getKnnClusters(method=infomap.community,type='PCA',name='infomap')
r$getKnnClusters(method=walktrap.community,type='PCA',name='walktrap')
#r$getEmbedding(type='PCA',M=5,perplexity=30,gamma=0.8)
r$getEmbedding(type='PCA',M=5,perplexity=50,gamma=0.8)
#r$getEmbedding(type='PCA',M=5,perplexity=100,gamma=0.8)
#r$embeddings$PCA$two <- NULL


#pdf('pca.pdf')
## plot
par(mfrow=c(2,2));
r$plotEmbedding(type='PCA',show.legend=F,clusterType='multilevel',mark.clusters=T,mark.cluster.cex=1,alpha=0.1, main="multilevel")
r$plotEmbedding(type='PCA',show.legend=F,clusterType='walktrap',mark.clusters=T,mark.cluster.cex=1,alpha=0.1, main="walktrap")
#r$plotEmbedding(type='PCA',show.legend=F,clusterType='infomap',mark.clusters=T,mark.cluster.cex=1,alpha=0.1, main="infomap")
r$plotEmbedding(type='PCA',show.legend=F,mark.clusters=T,mark.cluster.cex=1,groups=cell.groups.peter,alpha=0.1,main="Peter's clusters")
r$plotEmbedding(type='PCA',show.legend=F,mark.clusters=T,mark.cluster.cex=1,colors=r$depth,alpha=0.1,main="depth")
r$plotEmbedding(type='PCA',show.legend=F,mark.clusters=T,mark.cluster.cex=1,group=r$batch,alpha=0.1,main="batch")
#r$plotEmbedding(type='PCA',show.legend=F,mark.clusters=T,mark.cluster.cex=1,group=sample,alpha=0.1,main="sample")
#dev.off()

#pdf('tsne.pdf')
## tSNE embedding
r$getEmbedding(type='PCA',perplexity=50,verbose=T,embeddingType='tSNE')
#r$getEmbedding(type='PCA',perplexity=30,verbose=T,embeddingType='tSNE')
#r$getEmbedding(type='PCA',perplexity=100,verbose=T,embeddingType='tSNE')
#r$plotEmbedding(type='PCA',embeddingType='tSNE',groups=cl0,alpha=0.1,mark.clusters=T,mark.cluster.cex=1)
#r$plotEmbedding(type='PCA',embeddingType='tSNE',groups=cell.groups,alpha=0.1,mark.clusters=T,mark.cluster.cex=1)
#r$plotEmbedding(type='PCA',embeddingType='tSNE',colors=r$depth,alpha=0.1,mark.clusters=T,mark.cluster.cex=1)
par(mfrow=c(2,2));
r$plotEmbedding(type='PCA',show.legend=F,clusterType='multilevel',mark.clusters=T,mark.cluster.cex=1,alpha=0.1,embeddingType='tSNE', main="multilevel")
r$plotEmbedding(type='PCA',show.legend=F,clusterType='walktrap',mark.clusters=T,mark.cluster.cex=1,alpha=0.1,embeddingType='tSNE', main="walktrap")
r$plotEmbedding(type='PCA',show.legend=F,mark.clusters=T,mark.cluster.cex=1,groups=cell.groups.peter,alpha=0.1,embeddingType='tSNE',main="Peter's clusters")
r$plotEmbedding(type='PCA',show.legend=F,mark.clusters=T,mark.cluster.cex=1,colors=r$depth,alpha=0.1,main="depth",embeddingType='tSNE')
r$plotEmbedding(type='PCA',show.legend=F,mark.clusters=T,mark.cluster.cex=1,group=r$batch,alpha=0.1,main="batch",embeddingType='tSNE')
#r$plotEmbedding(type='PCA',show.legend=F,mark.clusters=T,mark.cluster.cex=1,group=sample,alpha=0.1,main="sample",embeddingType='tSNE')
#dev.off()



## consolidate groupings
cl0 <- r$clusters$PCA$multilevel[rownames(r$counts)]
cl1 <- r$clusters$PCA$walktrap[rownames(r$counts)]

cl0[cl0==9] <- NA
cl0[cl0==10] <- NA
cl0[cl0==8] <- NA
cl0[cl0==4] <- NA
cl0[cl1==2] <- NA
cl0[cl0==1 & cl1==3] <- NA

cl0 <- factor(cl0)
              
#pdf('plots.pdf')
r$plotEmbedding(type='PCA',groups=cl0, mark.clusters=TRUE, alpha=0.1)
r$plotEmbedding(type='PCA',groups=cl0, mark.clusters=TRUE, alpha=0.1,embeddingType='tSNE')
#dev.off()


lvec <- colSumByFac(r$misc[['rawCounts']],as.integer(cl0))[-1,] + 1
#lvec <- t(lvec)
lvec <- t(lvec/pmax(1,rowSums(lvec))) 
#lvec <- lvec/rowSums(lvec>0)
colnames(lvec) <- paste('cl',which(table(cl0)>0),sep='')
rownames(lvec) <- colnames(r$misc[['rawCounts']])
str(lvec)
ld <- jsDist(lvec); colnames(ld) <- rownames(ld) <- colnames(lvec)
hctree <- stats::hclust(as.dist(ld),method='ward')
plot(hctree,axes=F,sub="",ylab="",cex=0.8,main='initial cluster dendrogram')

## calculate cluster assignment probabilities for each cell
#ll <- as.matrix(t(r$misc[['rawCounts']] %*% log(lvec))) - colSums(lvec)
#ll <- t(ll)-apply(ll,2,max)+10; ll <- exp(ll); ll <- ll/rowSums(ll)
#cl1 <- as.factor(apply(ll,1,which.max)); 
#names(cl1) <- rownames(r$counts)
#table(cl1==cl0)
#cl1[which(cl1!=cl0)]
#cl0[which(cl0!=cl1)]


############ recursive differential peaks

## dend <- as.dendrogram(hctree)
## groups <- as.character(cl0)
## groups <- paste0('cl', groups)
## names(groups) <- names(cl0)

## dgall <- list()
## recursiveDf <- function(dend) {
##     g1 <- labels(dend[[1]])
##     if(length(g1)>1) { recursiveDf(dend[[1]]) }
##     g2 <- labels(dend[[2]])
##     if(length(g2)>1) { recursiveDf(dend[[2]]) }
##     g1cells <- names(groups[groups %in% g1])
##     g2cells <- names(groups[groups %in% g2])
##     difftest <- c(rep('g1', length(g1cells)), rep('g2', length(g2cells)))
##     names(difftest) <- c(g1cells, g2cells)
##     difftest <- factor(difftest)
##     dg <- r$getDifferentialGenes(upregulated.only = TRUE,groups=difftest)
##     dgall[[paste0(paste0(g1, collapse=":"), ' vs ', paste0(g2, collapse=":"))]] <<- dg
## }
## recursiveDf(dend)
## length(dgall)

## dgpeaks <- unlist(lapply(rev(names(dgall)), function(x) {
##     ## top 50 only
##     c(
##         names(sort(dgall[[x]][[1]], decreasing=TRUE))[1:10],
##         names(sort(dgall[[x]][[2]], decreasing=TRUE))[1:10]
##     )
## #    c(
## #        names(sort(dgall[[x]][[1]], decreasing=TRUE)),
## #        names(sort(dgall[[x]][[2]], decreasing=TRUE))
## #    )
## }))
## #dgpeaks <- na.omit(dgpeaks)
## dgpeaks <- na.omit(unique(dgpeaks))
## length(dgpeaks)

## ## consolidate cov matrix into groups to look at dgpeaks
## matsum2 <- log10(lvec+1)
## hist(matsum2)
## #matsum2 <- t(scale(t(matsum2)))
## matsum2 <- matsum2[dgpeaks,]
## matsum2 <- scale(matsum2)
## range(matsum2)
## matsum2[matsum2>2] <- 2
## matsum2[matsum2< -2] <- -2

## library(gplots)
## #heatmap.2(matsum2[1:5,], Rowv=NA, Colv=dend, col=colorRampPalette(c('blue', 'white', 'red'))(100), trace="none", scale="none", labRow=FALSE)
## #vi <- intersect(rownames(matsum2), names(dgall[[1]][[1]]))
## #length(vi)
## #heatmap.2(matsum2[vi,], Rowv=NA, Colv=dend, col=colorRampPalette(c('blue', 'white', 'red'))(100), trace="none", scale="none")
## heatmap.2(matsum2, Rowv=NA, Colv=dend, col=colorRampPalette(c('blue', 'white', 'red'))(100), trace="none", scale="none", labRow=FALSE)
## heatmap.2(matsum2, col=colorRampPalette(c('blue', 'white', 'red'))(100), trace="none", scale="none", labRow=FALSE)


## ## visualize peaks
## par(mfrow=c(2,5))
## r$plotEmbedding(type='PCA',groups=cl0, mark.clusters=TRUE, alpha=0.1,embeddingType='tSNE')
## lapply(1:length(dgall), function(i) {
##     dgup <- unlist(names(dgall[[i]][[1]]))
##     dgdown <- unlist(names(dgall[[i]][[2]]))
##     #v <- log2(rowMeans(r$counts[,dgup])/rowMeans(r$counts[,dgdown]))
##     v <- log2(rowMeans(r$counts[,dgdown])/rowMeans(r$counts[,dgup]))
##     #v <- rowMeans(r$counts[,dgdown]) - rowMeans(r$counts[,dgup])
##     v[is.infinite(v)] <- 0
##     v[is.na(v)] <- 0
##     range(v)
##     v[v>1] <- 1
##     v[v< -1] <- -1
##     ## highlight just groups
##     gs <- unlist(lapply(strsplit(names(dgall[i]), ' vs ')[[1]], function(x) strsplit(x, ":")))
##     v[names(groups)[!(groups %in% gs)]] <- NA
##     r$plotEmbedding(type='PCA',embedding='tSNE', show.legend=F,colors=v,alpha=0.1, main=names(dgall)[i])
## })

## old way
dg <- r$getDifferentialGenes(upregulated.only = TRUE,groups=cl0)
lapply(1:length(levels(cl0)), function(i) {
    r$plotEmbedding(type='PCA',embedding='tSNE',show.legend=F,colors=rowMeans(r$counts[,names(dg[[i]])]),alpha=0.1, main=levels(cl0)[i])
})



########################################################## GSEA
# map sites to genes
p2g <- hg38.peaks2Symbols(colnames(r$counts))
# run GO enrichment of the differential sites
gene2go <- hg38.getSymbols2Go(p2g)

library(GO.db)
# custom gene sets
d <- read.table("/home/pkharchenko/m/kun/ths/vc/cortical.quake.txt",sep="\t",header=T,stringsAsFactors=F)
cus2gene <- tapply(d$gene,as.factor(d$cluster),I); names(cus2gene) <- paste('quake',names(cus2gene),sep='.')
d <- read.table("/home/pkharchenko/m/kun/ths/vc/cortical.all50.txt",sep="\t",header=T,stringsAsFactors = F)
d$cluster <- paste('top50',d$cluster,sep='.')
cus2gene <- c(cus2gene,tapply(d$gene,as.factor(d$cluster),I))
d <- read.table("/home/pkharchenko/m/kun/ths/vc/cortical.deg.txt",sep="\t",header=F,stringsAsFactors = F); colnames(d)<-c('cluster','accession',"gene")
d$cluster <- paste('deg',d$cluster,sep='.')
cus2gene <- c(cus2gene,tapply(d$gene,as.factor(d$cluster),I))
gene2cus <- list2env(invert.string.list(cus2gene))

calculate.go.enrichment <- function(genelist,universe,pvalue.cutoff=1e-3,mingenes=3,env=entrez2GO,subset=NULL,list.genes=F) {
      all.genes <- unique(ls(env));
        # determine sizes
        universe <- unique(c(universe,genelist));
        ns <- length(intersect(genelist,all.genes));
        us <- length(intersect(universe,all.genes));
        #pv <- lapply(go.map,function(gl) { nwb <- length(intersect(universe,gl[[1]])); if(nwb<mingenes) { return(0.5)} else { p <- phyper(length(intersect(genelist,gl[[1]])),nwb,us-nwb,ns); return(ifelse(p>0.5,1.0-p,p)) }});

        # compile count vectors
        stab <- table(unlist(mget(as.character(genelist),env,ifnotfound=NA),recursive=T))
        utab <- table(unlist(mget(as.character(universe),env,ifnotfound=NA),recursive=T))
        if(!is.null(subset)) {
                stab <- stab[names(stab) %in% subset];
                    utab <- utab[names(utab) %in% subset];
                  }

        tabmap <- match(rownames(stab),rownames(utab))

        cv <- data.frame(cbind(utab,rep(0,length(utab)))); names(cv) <- c("u","s");
        cv$s[match(rownames(stab),rownames(utab))] <- as.vector(stab);
        cv <- na.omit(cv);
        cv <- cv[cv$u>mingenes,];
        pv <- phyper(cv$s,cv$u,us-cv$u,ns,lower.tail=F);
        pr <- dhyper(cv$s,cv$u,us-cv$u,ns)
        # correct for multiple hypothesis
        mg <- length(which(cv$u>mingenes));

        if(pvalue.cutoff<1) {
                ovi <- which(pv<0.5 & p.adjust(pr)<=pvalue.cutoff);
                    uvi <- which(pv>0.5 & p.adjust(pr)<=pvalue.cutoff);
                  } else {
                          ovi <- which(pv<0.5 & pr*mg<=pvalue.cutoff);
                              uvi <- which(pv>0.5 & pr*mg<=pvalue.cutoff);
                            }
        ovi <- ovi[order(pr[ovi])];
        uvi <- uvi[order(pr[uvi])];

        #return(list(over=data.frame(t=rownames(cv)[ovi],o=cv$s[ovi],u=cv$u[ovi],p=pr[ovi]*mg),under=data.frame(t=rownames(cv)[uvi],o=cv$s[uvi],u=cv$u[uvi],p=pr[uvi]*mg)))
        if(list.genes) {
                x <- mget(as.character(genelist),env,ifnotfound=NA);
                    df <- data.frame(id=rep(names(x),unlist(lapply(x,function(d) length(na.omit(d))))),go=na.omit(unlist(x)),stringsAsFactors=F)
                    ggl <- tapply(df$id,as.factor(df$go),I)
                    ovg <- as.character(unlist(lapply(ggl[rownames(cv)[ovi]],paste,collapse=" ")))
                    uvg <- as.character(unlist(lapply(ggl[rownames(cv)[uvi]],paste,collapse=" ")))
                    return(list(over=data.frame(t=rownames(cv)[ovi],o=cv$s[ovi],u=cv$u[ovi],p=pr[ovi]*mg,fe=cv$s[ovi]/(ns*cv$u[ovi]/us),genes=ovg),under=data.frame(t=rownames(cv)[uvi],o=cv$s[uvi],u=cv$u[uvi],p=pr[uvi]*mg,fe=cv$s[uvi]/(ns*cv$u[uvi]/us),genes=uvg)))
                  } else {
                          return(list(over=data.frame(t=rownames(cv)[ovi],o=cv$s[ovi],u=cv$u[ovi],p.raw=pr[ovi],fdr=p.adjust(pr)[ovi],p=pr[ovi]*mg,fe=cv$s[ovi]/(ns*cv$u[ovi]/us),fer=cv$s[ovi]/(length(genelist)*cv$u[ovi]/length(universe))),under=data.frame(t=rownames(cv)[uvi],o=cv$s[uvi],u=cv$u[uvi],p.raw=pr[uvi],fdr=p.adjust(pr)[uvi],p=pr[uvi]*mg,fe=cv$s[uvi]/(ns*cv$u[uvi]/us))))
                            }
      }

## test differential sites for enrichment
lapply(dg,function(s) {
    s <- names(s)
    calculate.go.enrichment(na.omit(unique(p2g[s])),na.omit(unique(p2g)),pvalue.cutoff=0.2,list.genes=F,env=gene2cus)$over;
})

## lapply(dgall,function(s) {
##     s <- names(s[[1]])
##     r1 <- calculate.go.enrichment(na.omit(unique(p2g[s])),na.omit(unique(p2g)),pvalue.cutoff=0.2,list.genes=F,env=gene2cus)$over;
##     s <- names(s[[2]])
##     r2 <- calculate.go.enrichment(na.omit(unique(p2g[s])),na.omit(unique(p2g)),pvalue.cutoff=0.2,list.genes=F,env=gene2cus)$over;
##     list(r1,r2)
## })

## calculate GO enrichment
s <- names(dg[[10]])
x <- calculate.go.enrichment(na.omit(unique(p2g[s])),na.omit(unique(p2g)),pvalue.cutoff=0.05,list.genes=F,env=gene2go)$over;
x$desc <- select(GO.db, keys = as.character(x$t), columns = c("TERM"), multiVals = "CharacterList")
x[x$u<2e3,]



##### check immune genes
pdf('quake.pdf')
par(mfrow=c(2,2), mar=rep(5,4))
immune <- names(p2g[which(p2g %in% cus2gene[["quake.mic"]])])
r$plotEmbedding(type='PCA',embedding='tSNE', show.legend=F,colors=rowMeans(r$counts[,immune]),alpha=0.1, main="quake.mic")
immune <- names(p2g[which(p2g %in% cus2gene[["quake.oli"]])])
r$plotEmbedding(type='PCA',embedding='tSNE', show.legend=F,colors=rowMeans(r$counts[,immune]),alpha=0.1, main="quake.oli")
immune <- names(p2g[which(p2g %in% cus2gene[["quake.ast"]])])
r$plotEmbedding(type='PCA',embedding='tSNE', show.legend=F,colors=rowMeans(r$counts[,immune]),alpha=0.1, main="quake.ast")

#immune <- names(p2g[which(p2g %in% c('CD45', 'PTPRC', 'AIF1', 'HLA-DR'))])
immune <- names(p2g[which(p2g %in% 'PTPRC')])
r$plotEmbedding(type='PCA',embedding='tSNE', show.legend=F,colors=rowMeans(r$counts[,immune]),alpha=0.1, main="CD45")
#r$plotEmbedding(type='PCA',embedding='tSNE', show.legend=F,colors=r$counts[,immune],alpha=1, main="CD45")
dev.off()


group <- as.character(cl0)
group[cl0==1] <- 'Mic'
group[cl0==2] <- 'Ex2'
group[cl0==3] <- 'Ast'
group[cl0==5] <- 'Oli'
group[cl0==6] <- 'Ex1'
group[cl0==7] <- 'In'
group <- factor(group)
names(group) <- names(cl0)
head(group)

r$plotEmbedding(type='PCA',groups=group, mark.clusters=TRUE, alpha=0.1)
r$plotEmbedding(type='PCA',groups=group, mark.clusters=TRUE, alpha=0.1,embeddingType='tSNE')

save(r, group, file='hR10_vc_comb_fixed.RData')

## compare with old
group.new <- group
load('hR10_vc_comb.RData')
par(mfrow=c(2,2), mar=rep(5,4))
r$plotEmbedding(type='PCA',groups=group, mark.clusters=TRUE, alpha=0.1)
r$plotEmbedding(type='PCA',groups=group, mark.clusters=TRUE, alpha=0.1,embeddingType='tSNE')
r$plotEmbedding(type='PCA',groups=group.new, mark.clusters=TRUE, alpha=0.1)
r$plotEmbedding(type='PCA',groups=group.new, mark.clusters=TRUE, alpha=0.1,embeddingType='tSNE')



## ######################## Write tracks
## load('multimap.RData')
## table(cell.groups)

## kun.bams <- read.table('../data/old/hR10_vc_used_files.bamlist')[,1]
## kun.bams <- unlist(lapply(kun.bams, function(x) strsplit(x, '/')[[1]][8]))
## cn <- gsub(".separated.sam.sorted.unique.","|", kun.bams)
## cn <- gsub(".VC.bam","",cn)
## cn <- paste0('hR10_vc_', cn)
## names(kun.bams) <- cn

#table(rownames(r$counts) %in% names(kun.bams))
#table(colnames(hR10_vc) %in% names(kun.bams))
#table(names(cell.groups) %in% colnames(hR10_vc))
#table(names(cell.groups) %in% names(kun.bams))
#names(cell.groups)[which(!(names(cell.groups) %in% names(kun.bams)))]

kun.bams <- bams
cols <- group
levels(cols)

lapply(levels(cols), function(i) {
    bamfiles.g <- kun.bams[names(cols)[which(cols==i)]]
    #bamfiles.g

    chrl <- paste('chr',c(1:22, 'X'),sep=''); names(chrl) <- chrl

    #path <- "../data-raw/"
    bamfiles <- na.omit(bamfiles.g)
    require(parallel)
    require(spp)
    bamdata <- mclapply(bamfiles,function(fname) read.bam.tags(fname)$tags,mc.cores=30)
    pdata <- mclapply(chrl,function(chr) list(tags=unlist(lapply(1:length(bamdata),function(i) { na.omit(bamdata[[i]][[chr]])})),cells=unlist(lapply(1:length(bamdata),function(i){ rep(i,length(na.omit(bamdata[[i]][[chr]])))}))),mc.cores=30)
    pdata <- lapply(pdata,function(d) { co <- order(abs(d$tags),decreasing=F); return(list(tags=d$tags[co],cells=d$cells[co]))})
    treads <- sum(unlist(lapply(pdata, function(d) length(d$tags))))

    bandwidth = 500
    step = 100
    thr = 5
    span = 10
    smoothed.density <- lapply(lapply(pdata,function(d) abs(d$tags)), function(d) {
            tc <- spp::window.tag.count(d, window.size=bandwidth, window.step=step)
                x <- seq(tc$x[1], tc$x[2], by=tc$step)
                y <- tc$y
                ## normalize per million reads
                y <- y / (treads/1e6)
                data.frame('x'=x,'y'=y)
            })
    names(smoothed.density) <- chrl

    writetdf <- function(dat,fname,name,description,tmp.dir=".",igvtools="/home/pkharchenko/igv/tools/igvtools",genome="hg38",save.wig=F,zip.wig=T) {
            if(length(grep("\\.tdf$",fname))>0) {
                        tfname <- paste(tmp.dir,gsub("\\.tdf",".wig",fname),sep="/");
                            } else {
                                        tfname <- paste(tmp.dir,paste(fname,"wig",sep="."),sep="/");
                                            }
                writewig(dat,tfname,feature=description,zip=F)
                paste(igvtools,"-f min,max,mean",tfname,fname)
                paste(igvtools,"toTDF -f max,mean",tfname,fname,genome)
                system(paste(igvtools,"toTDF -f max,mean",tfname,fname,genome))
                if(save.wig) {
                            if(zip.wig) {
                                            zf <- paste(fname,"zip",sep=".");
                                                        system(paste("zip \"",zf,"\" \"",tfname,"\"",sep=""));
                                                        system(paste("rm",tfname))
                                                    }
                                } else {
                                            system(paste("rm",tfname))
                                                }
            }
    writetdf(smoothed.density,paste0('results/pooled_g', i, '.tdf'),'scATAC-seq','cov_noAlt_peaks_spp_kunbams_withchrX_RMrepeatmask100_bandwidth500_step100_thr5_span10_fdr1e-07 normalized window peak count per million')

    })



## embed <- get('embeddings', slot(r,'.xData'))[['PCA']][['tSNE']]
## data <- cbind(embed, as.character(cl0[rownames(embed)]))
## head(data)
## colnames(data) <- c('tSNE1', 'tSNE2', 'group')
## plot(data[,1], data[,2], col=data[,3])
## write.table(data, file="results-multimap/forbrandon_data.txt", quote=FALSE, sep="\t")
## ### write bed
## direct.dg <- dg
## lapply(names(direct.dg),function(i) {
##         peaks <- names(direct.dg[[i]]);
##         peak.df <- data.frame(do.call(rbind,strsplit(peaks,":|-")),stringsAsFactors=F); colnames(peak.df) <- c('chr','start','end');
##         peak.df$name <- paste('peak',seq(1,nrow(peak.df)),sep='');
##         peak.df$score <- round(pmin(10,as.numeric(direct.dg[[i]]))*1e2)
##         fname <- paste0('results-multimap/forbrandon_direct.',i,'.diffPeaks.bed',sep='')
##         write(paste('track name="Direct',i,'" description="Direct ',i,' differentially accessible peaks (upregulated only)" visibility=2 useScore=1',sep=''),file=fname) # header
##         write.table(peak.df,append=T,file=fname,quote=F,row.names=F,col.names=F,sep='\t')
## })
## #lapply(1:length(dg), function(i) {
## #    s <- names(dg[[i]])
## #    dspeaks <- data.frame(do.call(rbind,strsplit(s,":|-")),stringsAsFactors=F)
## #    write.table(as.data.frame(dspeaks), file=paste0("results-multimap/forbrandon_pvsig_group", i, "peaksonly.bed"), quote=FALSE, sep="\t", row.names=FALSE, col.names=FALSE)
## #})





## ############ make plots
## cell.groups <- cl0

## par(mfrow=c(1,7), mar=rep(5,4))

## r$plotEmbedding(type='PCA',groups=cell.groups,mark.clusters=TRUE,main="Top PCs")
## r$plotEmbedding(type='PCA',embeddingType = 'tSNE',groups=cell.groups,mark.clusters=TRUE,main="Top PCs tSNE")
## r$plotEmbedding(type='PCA',embeddingType = 'tSNE',groups=cell.groups.peter,mark.clusters=TRUE,main="Peter's Colors")
## r$plotEmbedding(type='PCA',embeddingType = 'tSNE',colors=(r$depth-min(r$depth))/max(r$depth-min(r$depth)),main='depth',zlim=c(0,0.5))
## r$plotEmbedding(type='PCA',embeddingType = 'tSNE',groups=r$batch,main='batch')

## immune <- names(p2g[which(p2g %in% cus2gene[["quake.mic"]])])
## r$plotEmbedding(type='PCA',embedding='tSNE', show.legend=F,colors=rowMeans(r$counts[,immune]),alpha=0.1, main="quake.mic")
## immune <- names(p2g[which(p2g %in% 'PTPRC')])
## r$plotEmbedding(type='PCA',embedding='tSNE', show.legend=F,colors=rowMeans(r$counts[,immune]),alpha=1, main="CD45")


## blue.gs <- c('APBB1IP', 'PLXDC2', 'C10orf11', 'ST6GAL1', 'DOCK8', 'RASGEF1C', 'C3', 'ITPR2', 'DOCK4', 'ADAM28')
## immune <- names(p2g[which(p2g %in% blue.gs)])
## r$plotEmbedding(type='PCA',embedding='tSNE', show.legend=F,colors=rowMeans(r$counts[,immune]),alpha=0.1, main="quake.mic")

