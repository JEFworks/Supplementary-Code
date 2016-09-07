library(scde)


## t.read.count.files <- function(path) {
##     files <- list.files(path)
##     nn <- files; nn <- gsub(".count","",files); names(files) <- nn;
##     dd <- lapply(paste(path,files,sep="/"),read.delim,sep="\t",header=F,stringsAsFactors=F);
##     names(dd) <- nn;
##     gis <- dd[[1]][,1]
##     dm <- do.call(cbind,lapply(dd,function(d) d[match(gis,d[,1]),2]))
##     rownames(dm) <- gis;
##     dm <- dm[-grep("no_feature|ambiguous|too_low_aQual|not_aligned|alignment_not_unique",gis),]
##     dm;
## }

## cd <- t.read.count.files("../counts/");
## head(cd)
## dim(cd)
## save(cd, file='../data/cd.RData')


load('../data/cd.RData')
dim(cd)
load('../data/cell_subgroup.RData')
sg <- c(rep('npc', length(npc_sra)), rep('neuron', length(neuron_sra)))
names(sg) <- c(npc_sra, neuron_sra)
sg <- as.factor(sg)
sg
l2cols <- c("coral4", "olivedrab3", "skyblue2", "slateblue3")[as.integer(sg)]
names(l2cols) <- names(sg)

# filter
cd <- clean.counts(cd)
dim(cd)

## # error model
## knn <- knn.error.models(cd, groups=sg, k = ncol(cd)/4, n.cores = 10, min.count.threshold = 2, min.nonfailed = 5, max.model.plots = 10, verbose=1)
## save(knn, file='../data/knn.RData')

load('../data/knn.RData')
valid.cells <- knn$corr.a > 0
table(valid.cells)

## varinfo <- pagoda.varnorm(knn, counts = cd, trim = 3/ncol(cd), max.adj.var = 5, n.cores = 10, plot = TRUE)
## varinfo <- pagoda.subtract.aspect(varinfo, colSums(cd[, rownames(knn)]>0))
## save(varinfo, file='varinfo.RData')
load('varinfo.RData')

library(liger)
library(GO.db)
go.env <- org.Hs.GO2Symbol.list
desc <- AnnotationDbi::select(GO.db, keys = names(go.env), columns = c("TERM"), multiVals = "CharacterList")
names(go.env) <- paste(names(go.env), desc$TERM)
go.env <- list2env(go.env)

## pwpca <- pagoda.pathway.wPCA(varinfo, go.env, n.components = 1, n.cores = 10)
## save(pwpca, file='pwpca.RData')
## clpca <- pagoda.gene.clusters(varinfo, trim = 7.1/ncol(varinfo$mat), n.clusters = 50, n.cores = 10, plot = TRUE)
## save(clpca, file='clpca.RData')

load('pwpca.RData')
load('clpca.RData')
## #clpca <- NULL
## # get full info on the top aspects
tam <- pagoda.top.aspects(pwpca, clpca, n.cells = NULL, z.score = qnorm(0.01/2, lower.tail = FALSE))
hc <- pagoda.cluster.cells(tam, varinfo)
# determine overall cell clustering
tamr <- pagoda.reduce.loading.redundancy(tam, pwpca, clpca)
tamr2 <- pagoda.reduce.redundancy(tamr, distance.threshold = 0.6, plot = TRUE, cell.clustering = hc, labRow = NA, labCol = NA, box = TRUE, margins = c(0.5, 0.5), trim = 0)
#save(tam, tamr, tamr2, hc, file='tam_refined2.RData')
load('tam_refined.RData')

pagoda.view.aspects(tamr2, cell.clustering = hc, box = TRUE, labCol = NA, margins = c(0.5, 20), col.cols = rbind(l2cols[hc$labels]))



library(scde)
#load('tam.RData')
#load('pwpca.RData')
#load('varinfo.RData')
#clpca = NULL
load('qlen.RData')
qlencol <- c('red', 'blue')[as.factor(qlen)]; names(qlencol) <- names(qlen)
# cells that were selected
neurons <- c("SRR2967760", "SRR2967678", "SRR2967621", "SRR2967761", "SRR2967699", "SRR2967731", "SRR2967639", "SRR2967736", "SRR2967649", "SRR2967689", "SRR2967715", "SRR2967648", "SRR2967739", "SRR2967644", "SRR2967728", "SRR2967619", "SRR2967708", "SRR2967666", "SRR2967675", "SRR2967640", "SRR2967611", "SRR2967746", "SRR2967735", "SRR2967740", "SRR2967627", "SRR2967677", "SRR2967714", "SRR2967662", "SRR2967688", "SRR2967660", "SRR2967657", "SRR2967764", "SRR2967722", "SRR2967724", "SRR2967630", "SRR2967742", "SRR2967663", "SRR2967743", "SRR2967629", "SRR2967683", "SRR2967671", "SRR2967733", "SRR2967661", "SRR2967646", "SRR2967710")
npcs <- c("SRR2967741", "SRR2967712", "SRR2967751", "SRR2967713", "SRR2967738", "SRR2967697", "SRR2967737", "SRR2967707", "SRR2967612", "SRR2967658", "SRR2967711", "SRR2967766", "SRR2967720", "SRR2967659", "SRR2967615", "SRR2967667", "SRR2967610", "SRR2967642", "SRR2967654", "SRR2967616", "SRR2967752", "SRR2967638", "SRR2967664", "SRR2967681", "SRR2967680", "SRR2967608", "SRR2967609", "SRR2967669", "SRR2967684", "SRR2967674", "SRR2967690", "SRR2967716", "SRR2967758", "SRR2967651", "SRR2967636", "SRR2967672", "SRR2967771", "SRR2967622", "SRR2967747", "SRR2967694", "SRR2967668", "SRR2967653", "SRR2967691", "SRR2967703", "SRR2967702")
neuronsimm <- c("SRR2967719", "SRR2967706", "SRR2967665", "SRR2967718", "SRR2967637", "SRR2967647", "SRR2967754", "SRR2967682", "SRR2967750", "SRR2967762", "SRR2967768", "SRR2967700", "SRR2967723", "SRR2967623", "SRR2967765", "SRR2967618", "SRR2967632", "SRR2967634", "SRR2967685", "SRR2967744", "SRR2967676", "SRR2967730", "SRR2967652", "SRR2967695", "SRR2967770", "SRR2967693", "SRR2967698", "SRR2967620", "SRR2967625", "SRR2967655", "SRR2967633", "SRR2967705", "SRR2967709", "SRR2967626", "SRR2967756", "SRR2967759", "SRR2967650", "SRR2967614", "SRR2967721", "SRR2967656", "SRR2967613", "SRR2967643", "SRR2967641", "SRR2967631", "SRR2967624")
chosen <- rep(NA, length(hc$labels))
names(chosen) <- hc$labels
chosen[neurons] <- 'neurons'
chosen[neuronsimm] <- 'neuronsimm'
chosen[npcs] <- 'npcs'
chosencol <- rainbow(10)[as.factor(chosen)]
names(chosencol) <- names(chosen)
app <- make.pagoda.app(tamr2, tam, varinfo, go.env, pwpca, clpca, col.cols = rbind(l2cols[hc$labels], qlencol[hc$labels], chosencol[hc$labels], groupcol[hc$labels]), cell.clustering = hc, title = "npcsvsneurons")
show.app(app, "xiaochang", browse = TRUE, port = 1468)
saveRDS(app, file="xiaochang_app.rds")

save.image('scde_run.RData')
load('scde_run.RData')


library(scde)
# Take 2
#plot(hc)
groups <- as.factor(cutree(hc, 10))
library(RColorBrewer)
groupcol <- brewer.pal(10, "Set3")[groups]
names(groups) <- hc$labels
names(groupcol) <- hc$labels

load('qlen.RData')
qlencol <- c('red', 'blue')[as.factor(qlen)]; names(qlencol) <- names(qlen)
load('batch.RData')

pagoda.view.aspects(tamr2, cell.clustering = hc, box = TRUE, labCol = NA, margins = c(0.5, 20), col.cols = rbind(l2cols[hc$labels], groupcol[hc$labels], batchcol[hc$labels], qlencol[hc$labels]))

table(qlen == 100)
groups <- groups[names(qlen)]
groups <- groups[hc$labels[hc$order]]
qlen <- qlen[hc$labels[hc$order]]
# NPC1
table(groups == 3 & qlen == 100)
# NPC3
table(groups == 5 & qlen == 100)
# immature neuron
table((groups == 4 | groups == 6 | groups == 2 | groups == 7) & qlen == 100)
# mature neuron
table(groups == 1 & qlen == 100)
# smallest group is 15

g1 <- names(groups[groups == 3 & qlen == 100][1:19])
g2 <- names(groups[groups == 5 & qlen == 100][1:19])
g3 <- names(groups[(groups == 4 | groups == 6 | groups == 2 | groups == 7) & qlen == 100][1:19])
g4 <- rev(names(groups[groups == 1 & qlen == 100]))[1:19]
chosen <- rep(NA, length(groups))
names(chosen) <- names(groups)
chosen[g1] <- 'npcs1'
chosen[g2] <- 'npcs2'
chosen[g3] <- 'immneur'
chosen[g4] <- 'neur'
chosencol <- c('red', 'orange', 'green', 'blue')[as.factor(chosen)]
names(chosencol) <- names(chosen)

pdf('../heatmaps/pagoda_final.pdf', width=12, height=8)
pagoda.view.aspects(tamr2, cell.clustering = hc, box = TRUE, labCol = NA, margins = c(0.5, 20), col.cols = rbind(groupcol[hc$labels], l2cols[hc$labels], batchcol[hc$labels], qlencol[hc$labels], chosencol[hc$labels]), top=4)
dev.off()

dput(paste("samtools merge neurons100sub19.bam", paste(paste0(g4, '.bam'), collapse=" "), collapse=" "))
dput(paste("samtools merge neuronsimm100sub19.bam", paste(paste0(g3, '.bam'), collapse=" "), collapse=" "))
dput(paste("samtools merge npcs100g1sub19.bam", paste(paste0(g1, '.bam'), collapse=" "), collapse=" "))
dput(paste("samtools merge npcs100g2sub19.bam", paste(paste0(g2, '.bam'), collapse=" "), collapse=" "))






# Take 3 (diff exp)
groups <- cutree(hc, 10)
groups[hc$order]
groups[groups!=3] <- 1
groups <- as.factor(groups)
groups
# estimate gene expression prior
o.prior <- scde.expression.prior(models = knn, counts = cd, length.out = 400, show.plot = FALSE)
# diff exp
ediff <- scde.expression.difference(knn, cd, o.prior, groups  =  groups, n.randomizations  =  100, n.cores  =  10, verbose  =  1)
# top upregulated genes (tail would show top downregulated ones)
head(ediff[order(ediff$Z, decreasing  =  TRUE), ])
ediff['EOMES',]
#ediffsig <- ediff[abs(ediff$cZ) > 1.96,]
#ediffsig <- ediff[abs(ediff$cZ) > 3,]
ediffsigup <- ediff[ediff$cZ < -1.96,]
ediffsigdown <- ediff[ediff$cZ > 3.02,]

#ediffsig <- ediffsigdown
#ediffsig <- ediffsigup
ediffsig <- rbind(ediffsigdown, ediffsigup)
ediffsig <- ediffsig[order(ediffsig$cZ),]
dim(ediffsig)

# heatmap
mat <- log10(cd+1)
mat <- scale(mat)
mat <- t(scale(t(mat)))
#mat <- mat[rev(rownames(ediffsig))[1:20],]
mat <- mat[rownames(ediffsig),]
range(mat)
mat[mat < -1] <- -1
mat[mat > 1] <- 1
heatmap(mat[,hc$labels], Colv=as.dendrogram(hc), Rowv=NA, scale="row", col=colorRampPalette(c("blue", "white", "red"))(100))

# Test for gene set enrichment
library(liger)
library(GO.db)
go.env <- org.Hs.GO2Symbol.list
desc <- AnnotationDbi::select(GO.db, keys = names(go.env), columns = c("TERM"), multiVals = "CharacterList")
names(go.env) <- paste(names(go.env), desc$TERM)

#values <- abs(ediffsigup$mle)
#names(values) <- rownames(ediffsigup)
#values <- abs(ediffsigdown$mle)
#names(values) <- rownames(ediffsigdown)
values <- ediffsig$mle
names(values) <- rownames(ediffsig)
vi <- sapply(go.env, function(x) length(intersect(x, names(values)))) > 2
table(vi)
go.env.sub <- go.env[vi]

gs <- iterative.bulk.gsea(values, set.list=go.env.sub)
gs <- gs[order(gs$q.val, decreasing=FALSE),]
head(gs)
gsup <- gs[gs$sscore > 0 & gs$edge > 0,]
gsdown <- gs[gs$sscore < 0 & gs$edge < 0,]
head(gsup)
head(gsdown)
g <- "GO:0007049 cell cycle"
g <- "GO:0007154 cell communication"
gsea(values, go.env.sub[[g]])

# Test for specific enrichment of RBPs
xiaochang.gs <- read.csv("../data-raw/2016-3-28_human_RBP_MC350-Cell750.csv", header=FALSE, stringsAsFactors=FALSE)
xiaochang.gs <- xiaochang.gs[,1]
gsea(values, xiaochang.gs)
intersect(names(values), xiaochang.gs)

ediffsig['EOMES',]
write.csv(ediffsig, file="ediffsig_npcg1vseverythingelse.csv", quote=FALSE)





# NPC vs. neurons
groups <- cutree(hc, 2)
groups <- as.factor(groups)
names(groups) <- hc$labels
groups[colnames(cd)]
# estimate gene expression prior
o.prior <- scde.expression.prior(models = knn, counts = cd, length.out = 400, show.plot = FALSE)
# diff exp
ediff <- scde.expression.difference(knn, cd, o.prior, groups  =  groups, n.randomizations  =  100, n.cores  =  10, verbose  =  1)
# top upregulated genes (tail would show top downregulated ones)
head(ediff[order(ediff$Z, decreasing  =  TRUE), ])
#ediffsig <- ediff[abs(ediff$cZ) > 1.96,]
ediffsig <- ediff[abs(ediff$cZ) > mean(ediff$cZ)+1.96*sd(ediff$cZ),]
#ediffsigup <- ediff[ediff$cZ < -1.96,]
#ediffsigdown <- ediff[ediff$cZ > 3.02,]

#ediffsig <- ediffsigdown
#ediffsig <- ediffsigup
#ediffsig <- rbind(ediffsigdown, ediffsigup)
#ediffsig <- ediffsig[order(ediffsig$cZ),]
dim(ediffsig)

# heatmap
mat <- log10(cd+1)
mat <- scale(mat)
mat <- t(scale(t(mat)))
#mat <- mat[rev(rownames(ediffsig))[1:20],]
mat <- mat[rownames(ediffsig),]
range(mat)
mat[mat < -1] <- -1
mat[mat > 1] <- 1
heatmap(mat[,hc$labels], Colv=as.dendrogram(hc), Rowv=NA, scale="row", col=colorRampPalette(c("blue", "white", "red"))(100))

g <- "GO:0043005 neuron projection"
g <- "GO:0000278 mitotic cell cycle"







# plot
pdf('../heatmaps/pagoda.pdf', width=12, height=6)
pagoda.view.aspects(tamr2, cell.clustering = hc, box = TRUE, labCol = NA, margins = c(0.5, 20), col.cols = rbind(groupcol[hc$labels], l2cols[hc$labels], batchcol[hc$labels], qlencol[hc$labels], chosencol[hc$labels]), top=4)

pagoda.show.pathways(c("GO:0007399 nervous system development", "GO:0022008 neurogenesis", "GO:0043005 neuron projection"), varinfo, go.env, cell.clustering = hc, margins = c(1,5), show.cell.dendrogram = TRUE, showRowLabels = TRUE, showPC = TRUE)

pagoda.show.pathways(c("GO:0008283 cell proliferation"), varinfo, go.env, cell.clustering = hc, margins = c(1,5), show.cell.dendrogram = TRUE, showRowLabels = TRUE, showPC = TRUE)

pagoda.show.pathways(c("geneCluster.41"), varinfo, list2env(clpca$clusters), cell.clustering = hc, margins = c(1,5), show.cell.dendrogram = TRUE, showRowLabels = TRUE, showPC = TRUE)

pagoda.show.pathways(c("geneCluster.45"), varinfo, list2env(clpca$clusters), cell.clustering = hc, margins = c(1,5), show.cell.dendrogram = TRUE, showRowLabels = TRUE, showPC = TRUE)

dev.off()



#pdf('../heatmaps/markers.pdf', width=12, height=6)
pdf('../heatmaps/markers3.pdf', width=12, height=6)
# biased markers
markers <- c(
    "SCN2A","GRIK3","CDH6","NRCAM",
    "SOX11",
    "SLC24A2", "SOX4",  "DCX", "TUBB3","MAPT",
    "KHDRBS3",     "KHDRBS2",    "KHDRBS1",
    "RBFOX3",
    "CELF6",    "CELF5",    "CELF4",    "CELF3",    "CELF2",    "CELF1",
    "PTBP2", "PTBP1",    "ZFP36L2",
    "HMGN2", "PAX6", "SFRP1",
    "SOX2", "HES1", "NOTCH2", "CLU","HOPX",
    "MKI67","TPX2",
    "EOMES", "NEUROD4","HES6"
    )
# heatmap
mat <- log10(cd+1)
mat <- scale(mat)
mat <- t(scale(t(mat)))
mat <- mat[markers,]
range(mat)
mat[mat < -1] <- -1
mat[mat > 1] <- 1
heatmap(mat[,hc$labels], Colv=as.dendrogram(hc), Rowv=NA, scale="none", col=colorRampPalette(c("blue", "white", "red"))(100))
#heatmap(mat[,hc$labels], Colv=as.dendrogram(hc), scale="none", col=colorRampPalette(c("blue", "white", "red"))(100))
dev.off()

pdf('../heatmaps/rnabinding.pdf', width=12, height=6)
markers <- c(
    "KHDRBS3",     "KHDRBS2",    "KHDRBS1",
    "RBFOX3",
    "CELF6",    "CELF5",    "CELF4",    "CELF3",    "CELF2",    "CELF1",
    "PTBP2", "PTBP1",    "ZFP36L2"
    )
# heatmap
mat <- log10(cd+1)
mat <- scale(mat)
mat <- t(scale(t(mat)))
mat <- mat[markers,]
range(mat)
mat[mat < -1] <- -1
mat[mat > 1] <- 1
heatmap(mat[,hc$labels], Colv=as.dendrogram(hc), Rowv=NA, scale="none", col=colorRampPalette(c("blue", "white", "red"))(100))
dev.off()


# Biased
# estimate gene expression prior
o.prior <- scde.expression.prior(models = knn, counts = cd, length.out = 400, show.plot = FALSE)
# NPC vs. neurons
groups <- cutree(hc, 2)
groups <- as.factor(groups)
names(groups) <- hc$labels
groups[colnames(cd)]
# diff exp
ediff <- scde.expression.difference(knn, cd, o.prior, groups  =  groups, n.randomizations  =  100, n.cores  =  10, verbose  =  1)
# top upregulated genes (tail would show top downregulated ones)
g1 <- ediff[with(ediff, order(-ediff$cZ, -ediff$mle)), ][1:10,]
g2 <- ediff[with(ediff, order(ediff$cZ, ediff$mle)), ][1:10,]
markers <- c(rownames(g1), rownames(g2))
# heatmap
mat <- log10(cd+1)
mat <- scale(mat)
mat <- t(scale(t(mat)))
mat <- mat[markers,]
range(mat)
mat[mat < -1] <- -1
mat[mat > 1] <- 1
heatmap(mat[,hc$labels], Colv=as.dendrogram(hc), Rowv=NA, scale="none", col=colorRampPalette(c("blue", "white", "red"))(100))

# Eomes vs. everything else
groups <- cutree(hc, 10)
groups[hc$order]
groups[groups!=3] <- 1
groups <- as.factor(groups)
groups
# diff exp
ediff <- scde.expression.difference(knn, cd, o.prior, groups  =  groups, n.randomizations  =  100, n.cores  =  10, verbose  =  1)
# top upregulated genes (tail would show top downregulated ones)
g3 <- ediff[with(ediff, order(ediff$cZ, ediff$mle)), ][1:10,]
markers <- rownames(g3)
# heatmap
mat <- log10(cd+1)
mat <- scale(mat)
mat <- t(scale(t(mat)))
mat <- mat[markers,]
range(mat)
mat[mat < -1] <- -1
mat[mat > 1] <- 1
heatmap(mat[,hc$labels], Colv=as.dendrogram(hc), Rowv=NA, scale="none", col=colorRampPalette(c("blue", "white", "red"))(100))

# RG vs everything else
groups <- cutree(hc, 10)
groups[hc$order]
groups[groups==6] <- 8
groups[groups!=8] <- 1
groups <- as.factor(groups)
groups
# diff exp
ediff <- scde.expression.difference(knn, cd, o.prior, groups  =  groups, n.randomizations  =  100, n.cores  =  10, verbose  =  1)
# top upregulated genes (tail would show top downregulated ones)
g4 <- ediff[with(ediff, order(ediff$cZ, ediff$mle)), ][1:10,]
markers <- rownames(g4)
# heatmap
mat <- log10(cd+1)
mat <- scale(mat)
mat <- t(scale(t(mat)))
mat <- mat[markers,]
range(mat)
mat[mat < -1] <- -1
mat[mat > 1] <- 1
heatmap(mat[,hc$labels], Colv=as.dendrogram(hc), Rowv=NA, scale="none", col=colorRampPalette(c("blue", "white", "red"))(100))


markers <- c(rownames(g1), rownames(g3), rownames(g4))
# heatmap
mat <- log10(cd+1)
mat <- scale(mat)
mat <- t(scale(t(mat)))
mat <- mat[markers,]
range(mat)
mat[mat < -1] <- -1
mat[mat > 1] <- 1
heatmap(mat[,hc$labels], Colv=as.dendrogram(hc), Rowv=NA, scale="none", col=colorRampPalette(c("blue", "white", "red"))(100))
#dev.off()




#rbps from xiaochang
markers <- c(
        "ELAVL1",
        "ELAVL2",
        "ELAVL3",
        "ELAVL4",
        "RBFOX1",
        "RBFOX2",
        "RBFOX3",
        "PTBP1",
        "PTBP2",
        "QKI",
        "MSI1",
        "ZFP36L1",
        "ZFP36L2",
        "RALYL",
        "CELF1",
        "CELF2",
        "CELF3",
        "CELF4",
        "CELF5",
        "CELF6",
        "KHDRBS1",
        "KHDRBS2",
        "KHDRBS3",
        "NOVA1",
        "NOVA2",
        "SRRM4"
    )  
markers <- intersect(markers, rownames(cd))
# heatmap
mat <- log10(cd+1)
mat <- scale(mat)
mat <- t(scale(t(mat)))
mat <- mat[markers,]
range(mat)
mat[mat < -2] <- -2
mat[mat > 2] <- 2
heatmap(mat[,hc$labels], Colv=as.dendrogram(hc), Rowv=NA, scale="none", col=colorRampPalette(c("blue", "white", "red"))(100))
