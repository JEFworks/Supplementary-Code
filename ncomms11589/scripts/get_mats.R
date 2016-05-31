library(xlsx)

######
## RNA
######
file <- '../data/MD Anderson EX4 RNA Results.xlsx'
MDA1 <- read.xlsx(file, sheetName='MDA1', stringsAsFactors=F, header=T)
MDA2 <- read.xlsx(file, sheetName='MDA2', stringsAsFactors=F, header=T)
MDA3 <- read.xlsx(file, sheetName='MDA3', stringsAsFactors=F, header=T)
NormalB <- read.xlsx(file, sheetName='Normal B', stringsAsFactors=F, header=T)

process.mat <- function(mat) {
    rownames(mat) <- mat[,1]
    mat <- mat[,-1]
    mat[mat==999] <- 28
    mat <- 28-mat
    return(mat)
}
MDA1.rna <- process.mat(MDA1)
MDA2.rna <- process.mat(MDA2)
MDA3.rna <- process.mat(MDA3)
NormalB.rna <- process.mat(NormalB)

MDA1.geneexp <- MDA1.rna[, c(1:(96-8))]
MDA1.altsplice <- MDA1.rna[, c((96-7):96)]
MDA2.geneexp <- MDA2.rna[, c(1:(96-8))]
MDA2.altsplice <- MDA2.rna[, c((96-7):96)]
MDA3.geneexp <- MDA3.rna[, c(1:(96-8))]
MDA3.altsplice <- MDA3.rna[, c((96-7):96)]
NormalB.geneexp <- NormalB.rna[, c(1:(96-8))]
NormalB.altsplice <- NormalB.rna[, c((96-7):96)]

######
## Mutations
######

file <- '../data/MD Anderson Mutation Results.xlsx'
MDA1 <- read.xlsx(file, sheetName='MDA1', stringsAsFactors=F, header=F)
MDA2 <- read.xlsx(file, sheetName='MDA2', stringsAsFactors=F, header=F)
MDA3 <- read.xlsx(file, sheetName='MDA3', stringsAsFactors=F, header=F)
NormalB <- read.xlsx(file, sheetName='Normal B', stringsAsFactors=F, header=F)

process.mat <- function(mat) {
    ## fix colnames
    cn <- mat[c(1,2),-1]
    cn.name.ind <- seq(from=1, to=96, by=4)
    cn.names <- cn[1,cn.name.ind]
    cn.names.col <- unlist(lapply(cn.names, function(x) rep(x, 4)))
    cn.final <- paste(cn.names.col, cn[2,])
    cn.unique <- make.unique(cn.final)
    ## get numeric part
    mat <- mat[-c(1,2),]
    cell.names <- mat[,1]
    mat <- mat[,-1]
    ## set column and row names
    colnames(mat) <- cn.unique
    mat <- as.matrix(sapply(mat, as.numeric))
    rownames(mat) <- cell.names
    ## flip
    mat[mat==999] <- 28
    mat <- 28-mat
    return(mat)
}

MDA1.mut <- process.mat(MDA1)
MDA2.mut <- process.mat(MDA2)
MDA3.mut <- process.mat(MDA3)
NormalB.mut <- process.mat(NormalB)

save(MDA1.geneexp, MDA1.altsplice, MDA1.mut, file='../src/MDA1.RData')
save(MDA2.geneexp, MDA2.altsplice, MDA2.mut, file='../src/MDA2.RData')
save(MDA3.geneexp, MDA3.altsplice, MDA3.mut, file='../src/MDA3.RData')
save(NormalB.geneexp, NormalB.altsplice, NormalB.mut, file='../src/NormalB.RData')
