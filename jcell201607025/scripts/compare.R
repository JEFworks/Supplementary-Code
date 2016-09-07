# Compare with old results

# Read in old results
library(readxl)
mat0 <- read_excel('../data-raw/2016-3-17 HumanBrainSplicing_CP-VZ_SE_huttner_V2.xlsx')
rownames(mat0) <- mat0$event_name
range(abs(mat0$diff_CP_VZ)) # min diff 0.1
range(mat0$bayes_factor_CP_VZ) # bayes factor 5
# Note this is CP vs. VZ
# Filter more?
#vi <- abs(mat0$bayes_factor_CP_VZ) > 100
#table(vi)
#mat0 <- mat0[vi,]
#vi <- abs(mat0$diff_CP_VZ) > 0.2
#table(vi)
#mat0 <- mat0[vi,]
#vi <- mat0$CP_PSI > 0 & mat0$VZ_PSI > 0
#table(vi)
#mat0 <- mat0[vi,]

# Read in new results
#mat1 <- read.table('../miso3/comparisons/neurons_vs_neuronsimm/bayes-factors/neurons_vs_neuronsimm.miso_bf', sep='\t', header=TRUE)
mat1 <- read.table('../miso3/comparisons/neurons_vs_npcsg2/bayes-factors/neurons_vs_npcsg2.miso_bf', sep='\t', header=TRUE) #RG
#mat1 <- read.table('../miso3/comparisons/neuronsimm_vs_npcsg2/bayes-factors/neuronsimm_vs_npcsg2.miso_bf', sep='\t', header=TRUE)
#mat1 <- read.table('../miso3/comparisons/neuronsimm_vs_npcsg1/bayes-factors/neuronsimm_vs_npcsg1.miso_bf', sep='\t', header=TRUE)
#mat1 <- read.table('../miso3/comparisons/neurons_vs_npcsg1/bayes-factors/neurons_vs_npcsg1.miso_bf', sep='\t', header=TRUE)
#mat1 <- read.table('../miso3/comparisons/npcsg1_vs_npcsg2/bayes-factors/npcsg1_vs_npcsg2.miso_bf', sep='\t', header=TRUE)
rownames(mat1) <- mat1$event_name
dim(mat1) # 11531 events; 9389 events
# Filter for significant
vi <- abs(mat1$diff) > 0.1
table(vi)
mat1 <- mat1[vi,]
vi <- abs(mat1$bayes_factor) > 5
table(vi)
mat1 <- mat1[vi,]
head(mat1)
dim(mat1) # 768/11531; 554
table(mat1$diff > 0)

head(mat1)
head(mat1[order(abs(mat1$bayes_factor), decreasing=TRUE),])
head(mat1[order(abs(mat1$bayes_factor), decreasing=FALSE),])
head(mat1[order(abs(mat1$diff), decreasing=TRUE),])

# Match by event names
vi <- intersect(rownames(mat0), rownames(mat1))
length(vi)
dim(mat0)
dim(mat1)

mat0_sub <- mat0[vi,]
head(mat0_sub)
mat1_sub <- mat1[vi,]
head(mat1_sub)
lmfit <- lm(mat1_sub$diff ~ 0+mat0_sub$diff)
lmfit
summary(lmfit)
plot(mat0_sub$diff, mat1_sub$diff, pch=16, xlim=c(-1,1), ylim=c(-1,1))
abline(lmfit, col="red", pch=2)

mat1_sub$diff
mat0_sub$diff
rownames(mat1_sub)
# Match up event names with annotation
gene <- unlist(lapply(rownames(mat1_sub), function(e) {
    vi <- which(name %in% e)
    gsub("gsymbol=", "", annot[[vi]][5])
}))
paste(gene, rownames(mat1_sub), sep=": ")

pdf('../heatmaps/corr_neuronsvsnpcg2.pdf', height=3, width=3, useDingbats = FALSE)
dat <- data.frame('xvar'=mat0_sub$diff, 'yvar'=mat1_sub$diff)
library(ggplot2)
ggplot(dat, aes(x=xvar, y=yvar)) +
    geom_point() + 
        geom_abline(slope=coef(lmfit), col="blue") +
            theme_bw() + xlim(c(-1,1)) + ylim(c(-1,1)) 
dev.off()


# GSEA

# Annotate with gene names
gff3 <- read.table("/groups/pklab/jfan/Projects/Walsh_NPC_AltSplice/data-raw/miso_annotations_hg19_v2/SE.hg19.gff3", stringsAsFactors = FALSE)
annot <- gff3[,9]
annot <- lapply(annot, function(x) strsplit(x, ";")[[1]])
head(annot)

name <- unlist(lapply(annot, function(x) gsub("Name=", "", x[[1]][1])))
head(name)                 

# Match up event names with annotation
gene <- unlist(lapply(rownames(mat1), function(e) {
    vi <- which(name %in% e)
    gsub("gsymbol=", "", annot[[vi]][5])
}))
#gene <- gene[gene != "NA"]

values <- abs(mat1$diff)
#values <- mat1$diff
names(values) <- gene
head(values)
barplot(sort(values))

#run enrichment
library(liger)
library(GO.db)
go.env <- org.Hs.GO2Symbol.list
desc <- AnnotationDbi::select(GO.db, keys = names(go.env), columns = c("TERM"), multiVals = "CharacterList")
names(go.env) <- paste(names(go.env), desc$TERM)

vi <- sapply(go.env, function(x) length(intersect(x, names(values)))) > 2
table(vi)
go.env.sub <- go.env[vi]
head(go.env.sub)

gs <- iterative.bulk.gsea(values, set.list=go.env.sub, rank=TRUE)
#gs <- gs[gs$q.val < 0.05,]
gs <- gs[order(gs$q.val, decreasing=FALSE),]
head(gs)
dim(gs)
gsup <- gs[gs$sscore > 0 & gs$edge > 0,]
head(gsup, 20)

g <- "GO:0008219 cell death"
g <- "GO:0006414 translational elongation"
g <- "GO:0051052 regulation of DNA metabolic process"
g <- "GO:0051128 regulation of cellular component organization"
g <- "GO:0043228 non-membrane-bounded organelle"
g <- "GO:0043232 intracellular non-membrane-bounded organelle"
g <- "GO:0001558 regulation of cell growth"
gsea(values, go.env.sub[[g]])


g <- "GO:0000922 spindle pole"
g <- "GO:0031175 neuron projection development"
g <- "GO:0040011 locomotion"
g <- "GO:0008092 cytoskeletal protein binding"
gsea(values, go.env.sub[[g]])

values1 <- values
gsup1 <- gsup




# Test for specific enrichment of RBPs
xiaochang.gs <- read.csv("../data-raw/2016-3-28_human_RBP_MC350-Cell750.csv", header=FALSE, stringsAsFactors=FALSE)
xiaochang.gs <- xiaochang.gs[,1]
gsea(values, xiaochang.gs)
vi <- which(names(values) %in%  xiaochang.gs)
rownames(mat1)[vi]
gene[vi]






# Compare with old results
values <- abs(mat0$diff_CP_VZ)
names(values) <- mat0$gsymbol
head(values)

vi <- sapply(go.env, function(x) length(intersect(x, names(values)))) > 2
table(vi)
go.env.sub <- go.env[vi]

gs <- iterative.bulk.gsea(values, set.list=go.env.sub, rank=TRUE)
#gs <- gs[gs$q.val < 0.05,]
gs <- gs[order(gs$q.val, decreasing=FALSE),]
head(gs)
dim(gs)
gsup <- gs[gs$sscore > 0 & gs$edge > 0,]
head(gsup, 50)

g <- "GO:0000922 spindle pole"
g <- "GO:0031175 neuron projection development"
g <- "GO:0040011 locomotion"
g <- "GO:0008092 cytoskeletal protein binding"
gsea(values, go.env.sub[[g]])

gsup0 <- gsup
values0 <- values


gnam <- intersect(rownames(gsup1), rownames(gsup0))
plot(gsup0[gnam,]$sscore, gsup1[gnam,]$sscore)
plot(log10(gsup0[gnam,]$p), log10(gsup1[gnam,]$p))

vi <- which(gsup1[gnam,]$p < 0.05 & gsup0[gnam,]$p < 0.05)
gnam[vi]

g <- "GO:0032879 regulation of localization"
g <- "GO:0010941 regulation of cell death"
gsup0[g,]
gsup1[g,]
gsea(values1, go.env[[g]], rank=TRUE)
gsea(values0, go.env[[g]], rank=TRUE)






e <- "chr7:42974554:42974735:+@chr7:42976763:42976840:+@chr7:42976921:42977453:+"
e <- "chr20:43514344:43514527:+@chr20:43516289:43516383:+@chr20:43530172:43530474:+"
vi <- which(name %in% e)
vi
gsub("gsymbol=", "", annot[[vi]][5])
