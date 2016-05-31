#####
## Load data
#####

load('../src/MDA1.RData')
MDA1.mut.exp <- MDA1.mut[, grep("MUT", colnames(MDA1.mut))]
MDA1.nrm.exp <- MDA1.mut[, grep("NRM", colnames(MDA1.mut))]

load('../src/MDA2.RData')
MDA2.mut.exp <- MDA2.mut[, grep("MUT", colnames(MDA2.mut))]
MDA2.nrm.exp <- MDA2.mut[, grep("NRM", colnames(MDA2.mut))]

load('../src/MDA3.RData')
MDA3.mut.exp <- MDA3.mut[, grep("MUT", colnames(MDA3.mut))]
MDA3.nrm.exp <- MDA3.mut[, grep("NRM", colnames(MDA3.mut))]

load('../src/NormalB.RData')
NormalB.mut.exp <- NormalB.mut[, grep("MUT", colnames(NormalB.mut))]
NormalB.nrm.exp <- NormalB.mut[, grep("NRM", colnames(NormalB.mut))]

#####
## Prelim visualization
#####

mut.plot <- function(name, mut, nrm) {
    require(ggplot2)
    for (i in 1:ncol(mut)) {
        gene <- colnames(mut)[i]
        gene <- paste(strsplit(gene, " ")[[1]][1:2], collapse=" ")
        x <- mut[,i]
        y <- nrm[,i]
        df <- data.frame(x,y)

        p <- ggplot(df, aes(x=x, y=y)) + geom_point() +
            xlab('mutant allele expression') +
            ylab('normal allele expression') +
            ggtitle(paste('Mutation Calling for ', gene, ' \n in ', name, '\n', sep=''))
        print(p)
    }
}

#pdf('../plots/MDA1_mut_calling.pdf')
#mut.plot('MDA1', MDA1.mut.exp, MDA1.nrm.exp)
#dev.off()

#pdf('../plots/MDA2_mut_calling.pdf')
#mut.plot('MDA2', MDA2.mut.exp, MDA2.nrm.exp)
#dev.off()

#pdf('../plots/MDA3_mut_calling.pdf')
#mut.plot('MDA3', MDA3.mut.exp, MDA3.nrm.exp)
#dev.off()

#pdf('../plots/NormalB_mut_calling.pdf')
#mut.plot('NormalB', NormalB.mut.exp, NormalB.nrm.exp)
#dev.off()

#####
## Consolidate technical replicates by taking
## better value ie. higher value?
## average? (some times one failed assay though...)
#####

consolidate.tr <- function(mat) {
    m <- ncol(mat)
    ind <- (1:(m/2))*2
    mat.tr <- do.call(cbind, lapply(ind, function(i) {
        apply(mat[,c(i-1, i)], 1, mean)
    }))
    colnames(mat.tr) <- colnames(mat[,ind-1])
    return(mat.tr)
}

MDA1.mut.exp <- consolidate.tr(MDA1.mut.exp)
MDA1.nrm.exp <- consolidate.tr(MDA1.nrm.exp)

MDA2.mut.exp <- consolidate.tr(MDA2.mut.exp)
MDA2.nrm.exp <- consolidate.tr(MDA2.nrm.exp)

MDA3.mut.exp <- consolidate.tr(MDA3.mut.exp)
MDA3.nrm.exp <- consolidate.tr(MDA3.nrm.exp)

NormalB.mut.exp <- consolidate.tr(NormalB.mut.exp)
NormalB.nrm.exp <- consolidate.tr(NormalB.nrm.exp)

######
## Normalize by Normal
######

## Use Normal to fit expected mut allele detection based on normal allele detection
mut.nz <- NormalB.mut.exp!=0
lfits.line <- do.call(rbind, lapply(1:ncol(mut.nz), function(i) {
    print(i)
    ind <- mut.nz[,i]
    if(sum(ind) > 1) {
        nrm <- NormalB.nrm.exp[ind,i]
        mut <- NormalB.mut.exp[ind,i]
        fit <- nls(mut~a+b*nrm, start=c(a=-20, b=0), lower=c(a=-20, b=0), algorithm='port')
        lfit <- summary(fit)$coefficients[,1]
        lfit[is.na(lfit)] <- 0

        return(lfit)
    }
    else {
        return(matrix(c(0, 0), 1, 2))
    }
}))
rownames(lfits.line) <- colnames(mut.nz)
lfits.line


## Normalize
mut.normalize <- function(nrm.obs.all, mut.obs.all) {
    mut.normalized <- do.call(cbind, lapply(1:nrow(lfits.line), function(i) {
        lfit <- lfits.line[i,]
        ## observed normal allele expression
        nrm.obs <- nrm.obs.all[,i]
        ## expected mutant allele expression from regression
        pseudo <- 1e-6
        mut.exp <- lfit[1] + lfit[2]*nrm.obs + pseudo
        mut.exp[mut.exp<0] <- 0 ## can't be negative
        ## observed mutant allele expression
        mut.obs <- mut.obs.all[,i]
        ## normalize
        mut.nrm <- mut.obs - mut.exp
        mut.nrm[mut.nrm<0] <- 0 ## can't be negative
        return(mut.nrm)
    }))
    colnames(mut.normalized) <- colnames(mut.obs.all)
    return(mut.normalized)
}

MDA1.mut.normalized <- mut.normalize(MDA1.nrm.exp, MDA1.mut.exp)
MDA2.mut.normalized <- mut.normalize(MDA2.nrm.exp, MDA2.mut.exp)
MDA3.mut.normalized <- mut.normalize(MDA3.nrm.exp, MDA3.mut.exp)
NormalB.mut.normalized <- mut.normalize(NormalB.nrm.exp, NormalB.mut.exp)

j <- 4
plot(MDA1.nrm.exp[,j], MDA1.mut.normalized[,j])
plot(MDA1.nrm.exp[,j], MDA1.mut.exp[,j])
plot(MDA2.nrm.exp[,j], MDA2.mut.normalized[,j])
plot(MDA2.nrm.exp[,j], MDA2.mut.exp[,j])
plot(MDA3.nrm.exp[,j], MDA3.mut.normalized[,j])
plot(MDA3.nrm.exp[,j], MDA3.mut.exp[,j])
plot(NormalB.nrm.exp[,j], NormalB.mut.normalized[,j])
plot(NormalB.nrm.exp[,j], NormalB.mut.exp[,j])
## After normalization, expression of mutant allele:
## Below 2, definitely just Normal
## Between 2 and 5, unclear
## Above 5, likely mutant

mut.mat <- rbind(NormalB.mut.normalized,MDA1.mut.normalized, MDA2.mut.normalized, MDA3.mut.normalized)
nrm.mat <- rbind(NormalB.nrm.exp,MDA1.nrm.exp,MDA2.nrm.exp,MDA3.nrm.exp)
bad.cells <- mut.mat+nrm.mat==0
#mut.mat[bad.cells] <- NA

sg <-  unlist(lapply(rownames(mut.mat), function(x) strsplit(x, '_')[[1]][1]))
names(sg) <- rownames(mut.mat)
icols <- sg
icols[sg=='NB'] <- 'green'
icols[sg=='MDA1'] <- 'yellow'
icols[sg=='MDA2'] <- 'orange'
icols[sg=='MDA3'] <- 'red'

source('/Users/jf154/Dropbox/Grad/PKlab/jfan/Scripts/my_dendrogram.R')
my.heatmap(t(mut.mat),Colv=NA, Rowv=NA, zlim=c(-2, 2), col=colorRampPalette(c("blue","white","red"),space="Lab")(100),labRow=NULL,labCol=NA, scale="row", margin=c(1,15), ColSideColors=icols)

######
## Mutation calling
######

#mutation.call <- function(nrm, mut) {
#    call <- mut
#    ## play around with threshold; use plotting with calls to gage accuracy
#    call[mut > 5] <- 2 ## mutant
#    call[mut <= 5] <- 0 ## unclear
#    call[mut <= 2] <- 1 ## normal
#    call[(mut+nrm)==0] <- NA
#    return(call)
#}
## Updated method
## Mutation calling
mutation.call <- function(nrm, mut) {
    f <- mut/(mut+nrm)
    call <- f
    ## play around with threshold; use plotting with calls to gage accuracy
    call[f > 0.3] <- 2 ## mutant
    call[f <= 0.3] <- 0 ## unclear
    call[f < 0.15] <- 1 ## normal
    #call[(mut+nrm)==0] <- NA
    call[mut+nrm < 6] <- NA
    call[mut < 6 & nrm < 6] <- NA
    return(call)
}

MDA1.mut.call <- mutation.call(MDA1.nrm.exp, MDA1.mut.normalized)
head(MDA1.mut.call)
MDA2.mut.call <- mutation.call(MDA2.nrm.exp, MDA2.mut.normalized)
MDA3.mut.call <- mutation.call(MDA3.nrm.exp, MDA3.mut.normalized)
NormalB.mut.call <- mutation.call(NormalB.nrm.exp, NormalB.mut.normalized)

mut.call.mat <- rbind(NormalB.mut.call, MDA1.mut.call, MDA2.mut.call, MDA3.mut.call)

my.heatmap(t(mut.call.mat),Colv=NA, Rowv=NA, zlim=c(0, 2), col=colorRampPalette(c("grey","blue","red"),space="Lab")(100),labRow=NULL,labCol=NA, scale="none", margin=c(1,15), ColSideColors=icols)

## Lili's recommendations:
## subclonal mutations on top RPS15, EDEM2, PET, STAMBPL1, PLCG2
## clonal mutations XPO1,SF3B1, CCND3, TP53 at the bottom of the figure
## remove the ones without any mutation KMT2D, DGKA, SBNO1, DCK
genes <- c('XPO1','SF3B1', 'CCND3', 'TP53', 'RPS15', 'EDEM2', 'PET', 'STAMBPL1', 'PLCG2')
## Cathy's suggestions: remove less expressed genes
genes <- c('XPO1','SF3B1', 'CCND3', 'TP53', 'PLCG2')
m <- unlist(lapply(genes, function(g) colnames(mut.call.mat)[grep(g, colnames(mut.call.mat))]))
foo <- t(mut.call.mat)[m,]
foo[foo==2] <- -1
foo[is.na(foo)] <- 0
my.heatmap(foo,Colv=NULL, Rowv=NA, zlim=c(-1, 1), col=colorRampPalette(c("red",'white',"blue"),space="Lab")(100),labRow=NULL,labCol=NA, scale="none", margin=c(1,20), ColSideColors=icols)


######
## Plotting with calls
######

mut.plot.lab <- function(name, mut, nrm, label) {
    require(ggplot2)
    for (i in 1:ncol(mut)) {
        gene <- colnames(mut)[i]
        gene <- paste(strsplit(gene, " ")[[1]][1:2], collapse=" ")
        x <- mut[,i]
        y <- nrm[,i]
        l <- label[,i]
        l2 <- l
        l2[l==1] <- 'normal'
        l2[l==0] <- 'unclear'
        l2[l==2] <- 'mutant'
        l2[is.na(l)] <- 'NA'
        df <- data.frame(x,y,l2)

        p <- ggplot(df, aes(x=x, y=y, color=l2)) + geom_point() +
            xlab('mutant allele expression') +
            ylab('normal allele expression') +
            ggtitle(paste('Mutation Calling for ', gene, ' \n in ', name, '\n', sep=''))
        print(p)
    }
}
## Revised to make colors consistent

## Plot
mut.plot.lab <- function(name, mut, nrm, label, lfits.line) {
    require(ggplot2)
    for (i in 1:ncol(mut)) {
        gene <- colnames(mut)[i]
        gene <- paste(strsplit(gene, " ")[[1]][1:2], collapse=" ")
        x <- mut[,i]
        y <- nrm[,i]
        l <- label[,i]
        l2 <- l
        l2[l==1] <- 'normal'
        l2[l==0] <- 'unclear'
        l2[l==2] <- 'mutant'
        l2[is.na(l)] <- 'NA'
        call <- l2
        df <- data.frame(x,y,call)

        group.colors <- c('normal' = "blue", 'unclear' = "grey", 'mutant' ="red", 'NA' = "white")

        ## lfits.line is mut = a + b*nrm
        ## so -a/b + 1/b*mut=nrm
        a <- lfits.line[i,1]
        b <- lfits.line[i,2]

        p <- ggplot(df, aes(x=x, y=y, colour=call)) + scale_colour_manual(values=group.colors) +  geom_point() +
            xlab('log2(mutant allele expression)') +
            ylab('log2(normal allele expression)') +
            ggtitle(paste('Mutation Calling for ', gene, ' \n in ', name, '\n', sep='')) +
            xlim(0, 15) + ylim(0,15) +
            geom_abline(intercept = -a/b, slope=1/b)
        print(p)
    }
}

pdf('../plots/MDA1_mut_calling_lab_withabline.pdf')
mut.plot.lab('MDA1', MDA1.mut.exp, MDA1.nrm.exp, MDA1.mut.call,lfits.line)
dev.off()

pdf('../plots/MDA2_mut_calling_lab_withabline.pdf')
mut.plot.lab('MDA2', MDA2.mut.exp, MDA2.nrm.exp, MDA2.mut.call,lfits.line)
dev.off()

pdf('../plots/MDA3_mut_calling_lab_withabline.pdf')
mut.plot.lab('MDA3', MDA3.mut.exp, MDA3.nrm.exp, MDA3.mut.call,lfits.line)
dev.off()

pdf('../plots/NormalB_mut_calling_lab_withabline.pdf')
mut.plot.lab('NormalB', NormalB.mut.exp, NormalB.nrm.exp, NormalB.mut.call,lfits.line)
dev.off()

save(MDA1.mut.call, MDA2.mut.call, MDA3.mut.call, NormalB.mut.call, file='../src/mut_call.RData')


#####
## To show:
## 1. XPO1 mutation is founding mutation in all MDA samples, not found in B cells
## 2.1. SF3B1p.K667 (A2045C), TP53, CCND3 in one subsample
## 2.2. SF3B1p.E622D (A2033T) in other subsample
## 3. Order by PLCG2 mutations
#####

load('../src/mut_call.RData')
mut.call <- t(rbind(NormalB.mut.call, MDA1.mut.call, MDA2.mut.call, MDA3.mut.call))
head(mut.call)

#genes.order <- c("XPO1 G2439A CLL264.MUT7", "SF3B1 A2033T CLL253.MUT7", "SF3B1 G1914T CLL254.MUT7", "SF3B1 A2045C CLL252.MUT7", "TP53 C1211T CLL267.MUT7", "CCND3 C791T CLL255.MUT7", "PLCG2 T3636G CLL258.MUT7", "PLCG2 C2334T CLL260.MUT7", "PLCG2 G3191C CLL259.MUT7", "PLCG2 T3636A CLL257.MUT7")
#genes.order <- c("PLCG2 T3636A CLL257.MUT7", "PLCG2 G3191C CLL259.MUT7", "PLCG2 C2334T CLL260.MUT7", "PLCG2 T3636G CLL258.MUT7", "CCND3 C791T CLL255.MUT7", "TP53 C1211T CLL267.MUT7", "SF3B1 A2045C CLL252.MUT7", "SF3B1 G1914T CLL254.MUT7", "SF3B1 A2033T CLL253.MUT7", "XPO1 G2439A CLL264.MUT7")
#genes.order <- c("XPO1 G2439A CLL264.MUT7", "SF3B1 A2045C CLL252.MUT7", "SF3B1 A2033T CLL253.MUT7", "SF3B1 G1914T CLL254.MUT7", "PLCG2 T3636G CLL258.MUT7", "PLCG2 C2334T CLL260.MUT7", "PLCG2 G3191C CLL259.MUT7", "PLCG2 T3636A CLL257.MUT7")
genes.order <- c("SF3B1 A2045C CLL252.MUT7", "SF3B1 A2033T CLL253.MUT7", "SF3B1 G1914T CLL254.MUT7", "PLCG2 T3636G CLL258.MUT7", "PLCG2 C2334T CLL260.MUT7", "PLCG2 G3191C CLL259.MUT7", "PLCG2 T3636A CLL257.MUT7")
length(genes.order)

## focus on just these mutations
m <- mut.call[genes.order,]

## restrict to just celsl where we can confidently call all mutations
m[m==0] <- NA
good.cells <- colSums(is.na(m))==0
m <- m[,good.cells]

sg <-  unlist(lapply(colnames(m), function(x) strsplit(x, '_')[[1]][1]))
names(sg) <- colnames(m)
icols <- sg
icols[sg=='NB'] <- 'green'
icols[sg=='MDA1'] <- 'yellow'
icols[sg=='MDA2'] <- 'orange'
icols[sg=='MDA3'] <- 'red'
sg[sg=='NB']='A' ## reorder with NB first

## order by PLCG2 first
cells.order <- order(sg,
                     m[genes.order[1],],
                     m[genes.order[2],],
                     m[genes.order[3],],
                     m[genes.order[4],],
                     m[genes.order[5],],
                     m[genes.order[6],],
                     m[genes.order[7],])
## order by SF3B1 first
#cells.order <- order(
#                     m[genes.order[7],],
#                     m[genes.order[8],],
#                     m[genes.order[9],],
#                     m[genes.order[1],],
#                     m[genes.order[2],],
#                     m[genes.order[3],],
#                     m[genes.order[4],],
#                     m[genes.order[5],],
#                     m[genes.order[6],])

#pdf('../plots/All_mut_repconsolidated_binarized_heatmap_geneordered_SF3B1.pdf')
my.heatmap(m[,cells.order],Colv=NA, Rowv=NA, zlim=c(0, 2), col=colorRampPalette(c("grey","blue","red"),space="Lab")(100),labRow=NULL,labCol=NA, scale="none", margin=c(1,20), ColSideColors=icols[cells.order])
#dev.off()
#pdf('../plots/All_mut_repconsolidated_binarized_heatmap_subset.pdf')
#my.heatmap(m,Colv=NA, Rowv=NA, zlim=c(0, 2), col=colorRampPalette(c("grey","blue","red"),space="Lab")(100),labRow=NULL,labCol=NA, scale="none", margin=c(1,20), ColSideColors=icols)
#dev.off()

## number of cells
length(grep('NB', colnames(m)))
length(grep('MDA1', colnames(m)))
length(grep('MDA2', colnames(m)))
length(grep('MDA3', colnames(m)))

#####
## hclust order
## Done before consolidating technical replicates
#####

m <- t(mut.mat)
m.nb <- m[,sg=='NB']; m.nb[is.na(m.nb)] <- 0
nb.order <- hclust(d=dist(t(m.nb)), method='ward')
m.mda1 <- m[,sg=='MDA1']; m.mda1[is.na(m.mda1)] <- 0
mda1.order <- hclust(d=dist(t(m.mda1)), method='ward')
m.mda2 <- m[,sg=='MDA2']; m.mda2[is.na(m.mda2)] <- 0
mda2.order <- hclust(d=dist(t(m.mda2)), method='ward')
m.mda3 <- m[,sg=='MDA3']; m.mda3[is.na(m.mda3)] <- 0
mda3.order <- hclust(d=dist(t(m.mda3)), method='ward')

m.ordered <- cbind(m[,sg=='NB'][,nb.order$order], m[,sg=='MDA1'][,mda1.order$order], m[,sg=='MDA2'][,mda2.order$order], m[,sg=='MDA3'][,mda3.order$order])

my.heatmap(m.ordered, Colv=NA, Rowv=NA, zlim=c(-2, 2), col=colorRampPalette(c("blue","white","red"),space="Lab")(100),labRow=NULL,labCol=NA, scale="row", margin=c(1,15), ColSideColors=icols)


m <- t(mut.call.mat)
m.nb <- m[,sg=='NB']; m.nb[is.na(m.nb)] <- 0
nb.order <- hclust(d=dist(t(m.nb)), method='ward')
m.mda1 <- m[,sg=='MDA1']; m.mda1[is.na(m.mda1)] <- 0
mda1.order <- hclust(d=dist(t(m.mda1)), method='ward')
m.mda2 <- m[,sg=='MDA2']; m.mda2[is.na(m.mda2)] <- 0
mda2.order <- hclust(d=dist(t(m.mda2)), method='ward')
m.mda3 <- m[,sg=='MDA3']; m.mda3[is.na(m.mda3)] <- 0
mda3.order <- hclust(d=dist(t(m.mda3)), method='ward')

m.call.ordered <- cbind(m[,sg=='NB'][,nb.order$order], m[,sg=='MDA1'][,mda1.order$order], m[,sg=='MDA2'][,mda2.order$order], m[,sg=='MDA3'][,mda3.order$order])

my.heatmap(m.call.ordered,Colv=NA, Rowv=NA, zlim=c(0, 2), col=colorRampPalette(c("grey","blue","red"),space="Lab")(100),labRow=NULL,labCol=NA, scale="none", margin=c(1,15), ColSideColors=icols)


## focus on PLCG2
m <- t(mut.call.mat)
g <- grep('PLCG2', rownames(m))
m <- m[g,]
m.nb <- m[,sg=='NB']; m.nb[is.na(m.nb)] <- 0
nb.order <- hclust(d=dist(t(m.nb)), method='ward')
m.mda1 <- m[,sg=='MDA1']; m.mda1[is.na(m.mda1)] <- 0
mda1.order <- hclust(d=dist(t(m.mda1)), method='ward')
m.mda2 <- m[,sg=='MDA2']; m.mda2[is.na(m.mda2)] <- 0
mda2.order <- hclust(d=dist(t(m.mda2)), method='ward')
m.mda3 <- m[,sg=='MDA3']; m.mda3[is.na(m.mda3)] <- 0
mda3.order <- hclust(d=dist(t(m.mda3)), method='ward')

m.call.ordered <- cbind(m[,sg=='NB'][,nb.order$order], m[,sg=='MDA1'][,mda1.order$order], m[,sg=='MDA2'][,mda2.order$order], m[,sg=='MDA3'][,mda3.order$order])

my.heatmap(m.call.ordered,Colv=NA, Rowv=NA, zlim=c(0, 2), col=colorRampPalette(c("grey","blue","red"),space="Lab")(100),labRow=NULL,labCol=NA, scale="none", margin=c(1,15), ColSideColors=icols)

## stick with just one assay set, ignore replicates
genes <- c('SF3B1 A2045C CLL252.MUT7', 'PLCG2 T3636G CLL258.MUT7','PLCG2 C2334T CLL260.MUT7','PLCG2 G3191C CLL259.MUT7','PLCG2 T3636A CLL257.MUT7')
mco <- t(mut.call.mat)[genes, sg=='MDA3']
mco[mco==0] <- NA
## restrict to just cell where we can detect all 4 P variants
good.cells <- colSums(is.na(mco))==0
mco <- mco[,good.cells]
## restrict to just cells with the SF3B1K66T mutation
good.cells <- mco[genes[1],]==2
mco <- mco[,good.cells]
mco <- as.data.frame(mco)
mco <- mco[,order(mco[genes[1],], mco[genes[2],], mco[genes[3],], mco[genes[4],], mco[genes[5],])]
head(mco)

my.heatmap(as.matrix(mco),Colv=NA, Rowv=NA, zlim=c(0, 2), col=colorRampPalette(c("grey","blue","red"),space="Lab")(100),labRow=NULL,labCol=NA, scale="none", margin=c(1,25))







## Write tables

load('../src/mut_call.RData')
mut.call <- t(rbind(NormalB.mut.call, MDA1.mut.call, MDA2.mut.call, MDA3.mut.call))
head(mut.call)

mut.call.table <- mut.call
mut.call.table[mut.call==1] <- 'Normal'
mut.call.table[mut.call==0] <- 'Unclear'
mut.call.table[mut.call==2] <- 'Mutant'
mut.call.table[is.na(mut.call)] <- 'Unknown'

require('xlsx')
write.xlsx(mut.call.table, file='../mut_call_table/mut_call_table.xlsx')
