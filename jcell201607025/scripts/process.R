# Process SRA annotation file

annot <- readLines("../data-raw/GSE75140_series_matrix.txt.gz")
annot <- lapply(annot, strsplit, '\t')
mat <- annot[which(lapply(annot, function(x) length(x[[1]])) == 735)]
mat <- do.call(rbind, unlist(mat, recursive = FALSE))
cell_name <- mat[1,]
cell_name <- gsub("\"", "", cell_name)
id <- mat[2,]
id <- gsub("\"", "", id)
names(cell_name) <- id
cell_name <- cell_name[-1]
cell_name <- data.frame(cell_name)
head(cell_name)

annot <- readLines("../data-raw/SraRunTable_SRP066834.txt")
annot <- lapply(annot, strsplit, '\t')
mat <- annot[which(lapply(annot, function(x) length(x[[1]])) == 28)]
mat <- do.call(rbind, unlist(mat, recursive = FALSE))
colnames(mat) <- mat[1,]
mat <- mat[-1,]
rownames(mat) <- mat[, 'Sample_Name_s']
head(mat)
cell_info <- data.frame(mat)
head(cell_info)

save(cell_name, cell_info, file="../data/annot.RData")

# FPKM matrix
fpkm <- read.table("../data-raw/GSE75140_hOrg.fetal.master.data.frame.txt.gz", sep="\t", header=TRUE)
fpkm <- t(fpkm)
dim(fpkm)
colnames(fpkm) <- fpkm[1,]
fpkm <- fpkm[-1,]
fpkm <- fpkm[1:(nrow(fpkm)-1),] # last row is species for some reason; remove
head(fpkm)
mat <- do.call(rbind, lapply(1:nrow(fpkm), function(i) as.numeric(fpkm[i,]))) # convert to numeric
head(mat)
rownames(mat) <- rownames(fpkm)
colnames(mat) <- colnames(fpkm)
head(mat)
head(fpkm)
fpkm <- mat
save(fpkm, file="../data/fpkm.RData")

