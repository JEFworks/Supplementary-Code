#####################
# @desc Combined analysis of
# human visual cortex (data-raw2)
# (data-raw3)
# human cerebellum
# human hippocampus R11 run (rep2)
# human frontal cortex (fc9)
# human hippcampus R7 run (rep1)
#
# @date Feburary 1, 2017
# @author Jean Fan
######################

######################
# First call speaks using SPP
######################

require(spp)
require(parallel)

# valid chromosmes
chrl <- paste('chr',c(1:22, 'X'),sep=''); names(chrl) <- chrl;

########################################## Pooling data

####### hR10_vc
path <- "../data-raw2"
bamfiles <- list.files(path,pattern='.*bam$')
names(bamfiles) <- gsub(".bam","",bamfiles)
length(bamfiles)
kb <- readLines("~/Projects/Kun_Epigenetics/data/hR10_vc_rmmm_used_files.bamlist")
kb <- unlist(lapply(kb, function(x) { y=strsplit(x, '/'); y[[1]][length(y[[1]])] }))
names(kb) <- gsub(".bam","",kb)
vi <- names(bamfiles) %in% names(kb)
bamfiles <- bamfiles[vi]
bamfiles.hR10_vc <- paste(path, bamfiles, sep="/")
head(bamfiles.hR10_vc)

############# hR11_cer
path <- "../data-raw3/hR11_cer/ubam_rm_mm"
bamfiles <- list.files(path,pattern='.*bam$')
names(bamfiles) <- gsub(".bam","",bamfiles)
length(bamfiles)
kb <- readLines("/home/jfan/Projects/Kun_Epigenetics/data-raw3/hR11_cer/hR11_cer_rm_mm_merged_dfilter/hR11_cer_rm_mm_used_bams.txt")
kb <- unlist(lapply(kb, function(x) { y=strsplit(x, '/'); y[[1]][length(y[[1]])] }))
names(kb) <- gsub(".bam","",kb)
vi <- names(bamfiles) %in% names(kb)
bamfiles <- bamfiles[vi]
## guuuuhhh it looked like they didn't filter anything out....
bamfiles.hR11_cer <- paste(path, names(which(log10(reads+1)>2.5)), sep="/")
head(bamfiles.hR11_cer)

############# hR11_hip
path <- "../data-raw3/hR11_hip/ubam_rm_mm"
bamfiles <- list.files(path,pattern='.*bam$')
names(bamfiles) <- gsub(".bam","",bamfiles)
length(bamfiles)
reads <- unlist(mclapply(bamfiles, function(fname) {
    tryCatch({
        x <- read.bam.tags(paste(path, fname, sep='/'))$tags
        sum(unlist(lapply(x[chrl], length)))
    }, error = function(r) { return(NA) }) ## empty files
}, mc.cores=40))
names(reads) <- bamfiles
reads <- na.omit(reads)
hist(log10(reads+1), breaks=50)
abline(v=3, col="red")
bamfiles.hR11_hip <- paste(path, names(which(log10(reads+1)>3)), sep="/")
head(bamfiles.hR11_hip)

############# hR7_hip
path <- "../data-raw3/hR7_hip/hR7_hip/ubam_rm_mm"
bamfiles <- list.files(path,pattern='.*bam$')
names(bamfiles) <- gsub(".bam","",bamfiles)
length(bamfiles)
reads <- unlist(mclapply(bamfiles, function(fname) {
    tryCatch({
        x <- read.bam.tags(paste(path, fname, sep='/'))$tags
        sum(unlist(lapply(x[chrl], length)))
    }, error = function(r) { return(NA) }) ## empty files
}, mc.cores=40))
names(reads) <- bamfiles
reads <- na.omit(reads)
hist(log10(reads+1), breaks=50)
abline(v=2.5, col="red")
bamfiles.hR7_hip <- paste(path, names(which(log10(reads+1)>2.5)), sep="/")
head(bamfiles.hR7_hip)

############# hR7_fc9
path <- "../data-raw3/hR7_fc9/hR7_fc9/ubam_rm_mm"
bamfiles <- list.files(path,pattern='.*bam$')
names(bamfiles) <- gsub(".bam","",bamfiles)
length(bamfiles)
reads <- unlist(mclapply(bamfiles, function(fname) {
    tryCatch({
        x <- read.bam.tags(paste(path, fname, sep='/'))$tags
        sum(unlist(lapply(x[chrl], length)))
    }, error = function(r) { return(NA) }) ## empty files
}, mc.cores=40))
names(reads) <- bamfiles
reads <- na.omit(reads)
hist(log10(reads+1), breaks=50)
abline(v=2.5, col="red")
bamfiles.hR7_fc9 <- paste(path, names(which(log10(reads+1)>2.5)), sep="/")
head(bamfiles.hR7_fc9)

save(
    bamfiles.hR10_vc,
    bamfiles.hR11_cer,
    bamfiles.hR11_hip,
    bamfiles.hR7_hip,
    bamfiles.hR7_fc9,
    file="bamfiles_used.RData"
)

### all together now
all.bams <- c(
    bamfiles.hR10_vc,
    bamfiles.hR11_cer,
    bamfiles.hR11_hip,
    bamfiles.hR7_hip,
    bamfiles.hR7_fc9
)
length(all.bams)
head(all.bams)

### pool
bamdata <- mclapply(all.bams, function(fname) {
    tryCatch({
        x <- read.bam.tags(fname)$tags
    }, error = function(r) { return(NA) })
}, mc.cores=40)
table(is.na(bamdata))

## pooled data. $tags will represent Tn5 cut positions on each chromosome (sign carries the strand information), $cell will give the id (integer) of a cell giving the read; both will be sorted in the chromosome coordinate order
pdata <- mclapply(chrl,function(chr) list(tags=unlist(lapply(1:length(bamdata),function(i) { na.omit(bamdata[[i]][[chr]])})), cells=unlist(lapply(1:length(bamdata),function(i) { rep(i,length(na.omit(bamdata[[i]][[chr]])))}))), mc.cores=40)
# sort
pdata <- lapply(pdata,function(d) { co <- order(abs(d$tags),decreasing=F); return(list(tags=d$tags[co],cells=d$cells[co]))})
#head(pdata)

########################################## Remove reads in masked regions
pos2GRanges <- function(chr, pos) 
{
    require(GenomicRanges)
    require(IRanges)
    gr <- GenomicRanges::GRanges(seqnames = chr, ranges = IRanges(start = pos, width = 1))
    return(gr)
}

library(BSgenome.Hsapiens.UCSC.hg38.masked)
genome <- BSgenome.Hsapiens.UCSC.hg38.masked

rl <- 100
pdata.filtered <- lapply(chrl, function(chr) { 
    pdatachr <- pdata[[chr]] ## unfiltered
    ## expand mask region to account for read length
    mask <- GenomicRanges::GRanges(seqnames = chr, IRanges(start=masks(genome[[chr]])$RM@start-rl, end=masks(genome[[chr]])$RM@start + masks(genome[[chr]])$RM@width +rl))
    p <- pos2GRanges(chr, abs(pdatachr$tags))
    overlap <- GenomicRanges::findOverlaps(p, mask, ignore.strand=T)
    
    ## double check
    ##p[slot(overlap, 'queryHits')]
    ##mask
    
    ## remove these bad reads
    sd <- setdiff(1:length(pdatachr$tags), slot(overlap, 'queryHits'))
    pdatachr$tags <- pdatachr$tags[sd]
    pdatachr$cells <- pdatachr$cells[sd]
    
    return(pdatachr)
})
names(pdata.filtered) <- chrl


############################## Call Peaks

## call peaks
require(spp)
peaks.c <- function(x,thr=min(x),max.span=1) {
    storage.mode(x) <- storage.mode(thr) <- "double";
    storage.mode(max.span) <- "integer";
    results <- .Call("find_peaks",x,thr,max.span,"spp");
    return(results);
}

bandwidth = 500
step = 100
thr = 5
span = 10
fdr <- 10e-8

##smoothed.density <- get.smoothed.tag.density(lapply(pdata,function(d) abs(d$tags)),bandwidth=bandwidth,step=step,tag.shift=0,scale.by.dataset.size = F)
#smoothed.density <- lapply(lapply(pdata,function(d) abs(d$tags)), function(d) {
smoothed.density <- lapply(lapply(pdata.filtered,function(d) abs(d$tags)), function(d) {
    tc <- window.tag.count(d, window.size=bandwidth, window.step=step)
    x <- seq(tc$x[1], tc$x[2], by=tc$step)
    y <- tc$y
    data.frame('x'=x,'y'=y)
})    
names(smoothed.density) <- chrl
## this will calculate index positions of all local maxima on the density vector $y
peaks.all <- lapply(chrl, function(chr) {
    peak.indices <- peaks.c( smoothed.density[[chr]]$y, thr=thr, max.span=span);
    df <- data.frame(peak.position = smoothed.density[[chr]]$x[peak.indices], peak.magnitude = smoothed.density[[chr]]$y[peak.indices] );    
    return(df)
})
names(peaks.all) <- chrl
sum(unlist(lapply(peaks.all, nrow)))

peaks.filtered <- mclapply(chrl, function(chr) {
    test.x <- smoothed.density[[chr]]$x
    test.y <- smoothed.density[[chr]]$y
    
    ## init peaks
    df <- peaks.all[[chr]]
    fn <- ecdf(df$peak.magnitude)
    
    ## randomize
    set.seed(0)
    tags <- pdata[[chr]]$tags
    shuffle <- sort(runif(length(tags), min=range(tags)[1], max=range(tags)[2]), decreasing=TRUE)
    ## just jitter?
    #shuffle <- jitter(tags)
    ##shuffle.smooth <- get.smoothed.tag.density(list(chr=abs(shuffle)),bandwidth=bandwidth,step=step,tag.shift=0,scale.by.dataset.size = F)
    shuffle.smooth <- lapply(list(chr=abs(shuffle)), function(d) {
        tc <- window.tag.count(d, window.size=bandwidth, window.step=step)
        x <- seq(tc$x[1], tc$x[2], by=tc$step)
        y <- tc$y
        data.frame('x'=x,'y'=y)
    })
    peak.indices.shuffle <- peaks.c( shuffle.smooth[[1]]$y, thr=thr, max.span=span );
    df.shuffle <- data.frame(peak.position = shuffle.smooth[[1]]$x[peak.indices.shuffle], peak.magnitude = shuffle.smooth[[1]]$y[peak.indices.shuffle] );
    fn.shuffle <- ecdf(df.shuffle$peak.magnitude)
    
    #plot(df.shuffle[1:1000,], type="l")

    t <- knots(fn.shuffle)
    #V <- (1-fn.shuffle(t))*nrow(df.shuffle) + fn(t)*nrow(df)
    #S <- (fn.shuffle(t))*nrow(df.shuffle) + (1-fn(t))*nrow(df)

    V <- (1-fn.shuffle(t))*nrow(df.shuffle)
    S <- (1-fn(t))*nrow(df) 
    
    pseudo <- 1e-6
    fdrs <- (V + pseudo) / (V + S + pseudo)
    
    ## optimal threshold
    above <- fdrs > fdr
    intersect.points <- which(diff(above)!=0)
    threshold <- t[intersect.points][1]
    
    ##plot(t, log10(fdrs), type="l", xlim=c(1,100))
    ##abline(h=fdr, col="red")
    ##abline(v=threshold, col="red")
    
    table(df$peak.magnitude > threshold)
    table(df.shuffle$peak.magnitude > threshold)
    
    df <- peaks.all[[chr]]
    df.t <- df[df$peak.magnitude > threshold,]
    #plot(smoothed.density[[chr]], type="l")
    ##abline(v=df.t[,1], col="red")
    
    return(df.t)
}, mc.cores=length(chrl))
names(peaks.filtered) <- chrl
sum(unlist(lapply(peaks.filtered, nrow)))
sum(unlist(lapply(peaks.filtered, function(x) nrow(na.omit(x)))))

# http://jef.works/blog/2016/12/06/mapping-snps-and-peaks-to-genes-in-R/
#' Convert from ranges to GRanges
#'
#' @param df Dataframe with columns as sequence name, start, and end
#'
#' @returns GRanges version
#'
range2GRanges <- function(df) {
    require(GenomicRanges)
    require(IRanges)
    gr <- GenomicRanges::GRanges(
        seqnames = df[,1],
        ranges=IRanges(start = df[,2], end = df[,3])
        )
    return(gr)
}

## spp results
peaks.spp <- do.call(rbind, lapply(chrl, function(chr) {
    df <- peaks.filtered[[chr]]
    #df <- peaks.all[[chr]]
    #data.frame(chr, df-bandwidth/2, df+bandwidth/2)
    data.frame(chr, df[,1]-bandwidth/2, df[,1]+bandwidth/2)
}))
peaks.spp <- range2GRanges(peaks.spp)
names(peaks.spp) <- paste(peaks.spp)
head(peaks.spp)

## write out bed file
gr <- peaks.spp
df <- data.frame(seqnames=seqnames(gr),
                   starts=start(gr)-1,
                   ends=end(gr),
                   names=c(rep(".", length(gr))),
                   scores=c(rep(".", length(gr))),
                   strands=strand(gr))
head(df)
write.table(df, file=paste0("results-2-1-2017/peaks_spp_RMrepeatmask100_bandwidth", bandwidth, "_step", step, "_thr", thr, "_span", span, "_fdr", fdr, ".bed"), quote=F, sep="\t", row.names=F, col.names=F)
save(peaks.spp, file=paste0("results-2-1-2017/peaks_spp_RMrepeatmask100_bandwidth", bandwidth, "_step", step, "_thr", thr, "_span", span, "_fdr", fdr, ".RData"))   


######################## Get matrix

#cov <- do.call(cbind, mclapply(all.bams, function(fname) {
# errors occur if all processed together
# long vectors not supported yet: fork.c:376

curbams <- list(
    bamfiles.hR10_vc,
    bamfiles.hR11_cer,
    bamfiles.hR11_hip,
    bamfiles.hR7_hip,
    bamfiles.hR7_fc9
)
names(curbams) <- c(
    "bamfiles.hR10_vc",
    "bamfiles.hR11_cer",
    "bamfiles.hR11_hip",
    "bamfiles.hR7_hip",
    "bamfiles.hR7_fc9"
)

require(BiocParallel)
#rev(1:length(curbams))
lapply(c(3,2), function(i) {
    cb <- curbams[[i]]
    cbn <- names(curbams)[i]
    print(cbn)
    cov <- do.call(cbind, bplapply(cb, function(fname) {
        pdata.cell <- read.bam.tags(fname)$tags
        pdata.cell <- pdata.cell[chrl]
        smoothed.density.cell <- unlist(lapply(chrl, function(chr) {
            pdata.cell.chr <- pdata.cell[[chr]]
            peaks.chr <- peaks.filtered[[chr]]
            tca <- window.tag.count.around(pdata.cell.chr, window.size=bandwidth, peaks.chr[,1])$y
        }))
        names(smoothed.density.cell) <- names(peaks.spp)
        smoothed.density.cell
    }, BPPARAM=MulticoreParam(workers=40)))
    colnames(cov) <- cb

    save(cov, file=paste0("results-2-1-2017/", cbn, "_cov_peaks_spp_RMrepeatmask100_bandwidth", bandwidth, "_step", step, "_thr", thr, "_span", span, "_fdr", fdr, ".RData"))
    rm(cov);
    gc();
})








##################### Repeat without peak filtering for comparison

smoothed.density <- lapply(lapply(pdata,function(d) abs(d$tags)), function(d) {
    tc <- window.tag.count(d, window.size=bandwidth, window.step=step)
    x <- seq(tc$x[1], tc$x[2], by=tc$step)
    y <- tc$y
    data.frame('x'=x,'y'=y)
})    
names(smoothed.density) <- chrl
## this will calculate index positions of all local maxima on the density vector $y
peaks.all <- lapply(chrl, function(chr) {
    peak.indices <- peaks.c( smoothed.density[[chr]]$y, thr=thr, max.span=span);
    df <- data.frame(peak.position = smoothed.density[[chr]]$x[peak.indices], peak.magnitude = smoothed.density[[chr]]$y[peak.indices] );    
    return(df)
})
names(peaks.all) <- chrl
sum(unlist(lapply(peaks.all, nrow)))

peaks.filtered <- mclapply(chrl, function(chr) {
    test.x <- smoothed.density[[chr]]$x
    test.y <- smoothed.density[[chr]]$y
    
    ## init peaks
    df <- peaks.all[[chr]]
    fn <- ecdf(df$peak.magnitude)
    
    ## randomize
    set.seed(0)
    tags <- pdata[[chr]]$tags
    shuffle <- sort(runif(length(tags), min=range(tags)[1], max=range(tags)[2]), decreasing=TRUE)
    ## just jitter?
    #shuffle <- jitter(tags)
    ##shuffle.smooth <- get.smoothed.tag.density(list(chr=abs(shuffle)),bandwidth=bandwidth,step=step,tag.shift=0,scale.by.dataset.size = F)
    shuffle.smooth <- lapply(list(chr=abs(shuffle)), function(d) {
        tc <- window.tag.count(d, window.size=bandwidth, window.step=step)
        x <- seq(tc$x[1], tc$x[2], by=tc$step)
        y <- tc$y
        data.frame('x'=x,'y'=y)
    })
    peak.indices.shuffle <- peaks.c( shuffle.smooth[[1]]$y, thr=thr, max.span=span );
    df.shuffle <- data.frame(peak.position = shuffle.smooth[[1]]$x[peak.indices.shuffle], peak.magnitude = shuffle.smooth[[1]]$y[peak.indices.shuffle] );
    fn.shuffle <- ecdf(df.shuffle$peak.magnitude)
    
    #plot(df.shuffle[1:1000,], type="l")

    t <- knots(fn.shuffle)
    #V <- (1-fn.shuffle(t))*nrow(df.shuffle) + fn(t)*nrow(df)
    #S <- (fn.shuffle(t))*nrow(df.shuffle) + (1-fn(t))*nrow(df)

    V <- (1-fn.shuffle(t))*nrow(df.shuffle)
    S <- (1-fn(t))*nrow(df) 
    
    pseudo <- 1e-6
    fdrs <- (V + pseudo) / (V + S + pseudo)
    
    ## optimal threshold
    above <- fdrs > fdr
    intersect.points <- which(diff(above)!=0)
    threshold <- t[intersect.points][1]
    
    ##plot(t, log10(fdrs), type="l", xlim=c(1,100))
    ##abline(h=fdr, col="red")
    ##abline(v=threshold, col="red")
    
    table(df$peak.magnitude > threshold)
    table(df.shuffle$peak.magnitude > threshold)
    
    df <- peaks.all[[chr]]
    df.t <- df[df$peak.magnitude > threshold,]
    #plot(smoothed.density[[chr]], type="l")
    ##abline(v=df.t[,1], col="red")
    
    return(df.t)
}, mc.cores=length(chrl))
names(peaks.filtered) <- chrl
sum(unlist(lapply(peaks.filtered, nrow)))
sum(unlist(lapply(peaks.filtered, function(x) nrow(na.omit(x)))))

## spp results
peaks.spp <- do.call(rbind, lapply(chrl, function(chr) {
    df <- peaks.filtered[[chr]]
    #df <- peaks.all[[chr]]
    #data.frame(chr, df-bandwidth/2, df+bandwidth/2)
    data.frame(chr, df[,1]-bandwidth/2, df[,1]+bandwidth/2)
}))
peaks.spp <- range2GRanges(peaks.spp)
names(peaks.spp) <- paste(peaks.spp)
head(peaks.spp)

## write out bed file
gr <- peaks.spp
df <- data.frame(seqnames=seqnames(gr),
                   starts=start(gr)-1,
                   ends=end(gr),
                   names=c(rep(".", length(gr))),
                   scores=c(rep(".", length(gr))),
                   strands=strand(gr))
head(df)
write.table(df, file=paste0("results-2-1-2017/peaks_spp_WITHrepeats_bandwidth", bandwidth, "_step", step, "_thr", thr, "_span", span, "_fdr", fdr, ".bed"), quote=F, sep="\t", row.names=F, col.names=F)
save(peaks.spp, file=paste0("results-2-1-2017/peaks_spp_WITHrepeats_bandwidth", bandwidth, "_step", step, "_thr", thr, "_span", span, "_fdr", fdr, ".RData"))   

require(BiocParallel)
lapply(rev(1:length(curbams)), function(i) {
    cb <- curbams[[i]]
    cbn <- names(curbams)[i]
    print(cbn)
    cov <- do.call(cbind, bplapply(cb, function(fname) {
        pdata.cell <- read.bam.tags(fname)$tags
        pdata.cell <- pdata.cell[chrl]
        smoothed.density.cell <- unlist(lapply(chrl, function(chr) {
            pdata.cell.chr <- pdata.cell[[chr]]
            peaks.chr <- peaks.filtered[[chr]]
            tca <- window.tag.count.around(pdata.cell.chr, window.size=bandwidth, peaks.chr[,1])$y
        }))
        names(smoothed.density.cell) <- names(peaks.spp)
        smoothed.density.cell
    }, BPPARAM=MulticoreParam(workers=40)))
    colnames(cov) <- cb

    save(cov, file=paste0("results-2-1-2017/", cbn, "_cov_peaks_spp_WITHrepeats_bandwidth", bandwidth, "_step", step, "_thr", thr, "_span", span, "_fdr", fdr, ".RData"))
    rm(cov);
    gc();
})
