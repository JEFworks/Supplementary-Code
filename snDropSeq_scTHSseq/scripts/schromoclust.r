Pagoda2 <- setRefClass(
  "Pagoda2",
  fields=c('counts','clusters','graphs','reductions','embeddings','diffgenes','pathways','n.cores','misc','batch','modelType','verbose','depth','batchNorm','mat'),
  methods = list(
    initialize=function(x, ..., modelType='plain',batchNorm='glm',n.cores=30,verbose=TRUE,min.cells.per.gene=30,trim=round(min.cells.per.gene/2),lib.sizes=NULL) {
      # # init all the output lists
      embeddings <<- list();
      graphs <<- list();
      diffgenes <<- list();
      reductions <<-list();
      clusters <<- list();
      pathways <<- list()
      misc <<-list();
      batch <<- NULL;
      counts <<- NULL;
      if(!missing(x) && class(x)=='Pagoda2') {
        callSuper(x, ..., modelType=modelType, batchNorm=batchNorm, n.cores=n.cores);
      } else {
        callSuper(..., modelType=modelType, batchNorm=batchNorm, n.cores=n.cores,verbose=verbose);
        if(!missing(x) && is.null(counts)) { # interpret x as a countMatrix
          setCountMatrix(x,min.cells.per.gene=min.cells.per.gene,trim=trim,lib.sizes=lib.sizes)
        }
      }
    },
    # provide the initial count matrix, and estimate deviance residual matrix (correcting for depth and batch)
    setCountMatrix=function(countMatrix,depthScale=1,min.cells.per.gene=30,max.cell.fraction.per.gene=0.08, trim=round(min.cells.per.gene/2),lib.sizes=NULL) {
      # check names
      if(any(duplicated(rownames(countMatrix)))) {
        stop("duplicate gene names are not allowed - please reduce")
      }
      if(!is.null(batch)) {
        if(!all(colnames(countMatrix) %in% names(batch))) { stop("the supplied batch vector doesn't contain all the cells in its names attribute")}
        colBatch <- as.factor(batch[colnames(countMatrix)])
        batch <<- colBatch;
      }

      depth <<- Matrix::colSums(countMatrix);

      countMatrix@x <- as.numeric(countMatrix@x>0); # binarize

      counts <<- t(countMatrix)
      counts <<- counts[,diff(counts@p)>min.cells.per.gene]
      counts <<- counts[,diff(counts@p)<max.cell.fraction.per.gene*nrow(counts)]

      if(!is.null(lib.sizes)) {
        if(!all(rownames(counts) %in% names(lib.sizes))) { stop("the supplied lib.sizes vector doesn't contain all the cells in its names attribute")}
        lib.sizes <- lib.sizes[rownames(counts)]
        depth <<- lib.sizes/mean(lib.sizes)*mean(depth);
      }

      misc[['rawCounts']] <<- counts;

      cat(nrow(counts),"cells,",ncol(counts),"sites; normalizing ... ")
      # get normalized matrix
      if(modelType=='linearObs') { # this shoudln't work well, since the depth dependency is not completely normalized out

        # winsorize in normalized space first in hopes of getting a more stable depth estimate
        if(trim>0) {
          counts <<- counts/as.numeric(depth);
          inplaceWinsorizeSparseCols(counts,trim);
          counts <<- counts*as.numeric(depth);
          if(is.null(lib.sizes)) {
            depth <<- round(Matrix::rowSums(counts))
          }
        }


        ldepth <- log(depth);

        # rank cells, cut into n pieces
        n.depth.slices <- 20;
        #depth.fac <- as.factor(floor(rank(depth)/(length(depth)+1)*n.depth.slices)+1); names(depth.fac) <- rownames(counts);
        depth.fac <- cut(cumsum(sort(depth)),breaks=seq(0,sum(depth),length.out=n.depth.slices)); names(depth.fac) <- rownames(counts);
        depth.fac <- depth.fac[rank(depth)]
        # dataset-wide gene average
        gene.av <- (Matrix::colSums(counts)+n.depth.slices)/(sum(depth)+n.depth.slices)

        # pooled counts, df for all genes
        tc <- colSumByFac(counts,as.integer(depth.fac))[-1,]
        tc <- log(tc+1)- log(as.numeric(tapply(depth,depth.fac,sum))+1)
        md <- log(as.numeric(tapply(depth,depth.fac,mean)))
        # combined lm
        cm <- lm(tc ~ md)
        colnames(cm$coef) <- colnames(counts)
        # adjust counts
        # predict log(p) for each non-0 entry
        count.gene <- rep(1:counts@Dim[2],diff(counts@p))
        exp.x <- exp(log(gene.av)[count.gene] - cm$coef[1,count.gene] - ldepth[counts@i+1]*cm$coef[2,count.gene])
        counts@x <<- counts@x*exp.x/(depth[counts@i+1]/depthScale); # normalize by depth as well
        # performa another round of trim
        if(trim>0) {
          inplaceWinsorizeSparseCols(counts,trim);
        }


        # regress out on non-0 observations of ecah gene
        #non0LogColLmS(counts,mx,ldepth)
      } else if(modelType=='plain') {
        cat("using plain model ")

        if(!is.null(batch)) {
          cat("batch ... ")

          # dataset-wide gene average
          gene.av <- (Matrix::colSums(counts)+length(levels(batch)))/(sum(depth)+length(levels(batch)))

          # pooled counts, df for all genes
          tc <- colSumByFac(counts,as.integer(batch))[-1,]
          tc <- t(log(tc+1)- log(as.numeric(tapply(depth,batch,sum))+1))
          bc <- exp(tc-log(gene.av))

          # adjust every non-0 entry
          count.gene <- rep(1:counts@Dim[2],diff(counts@p))
          counts@x <<- counts@x/bc[cbind(count.gene,as.integer(batch)[counts@i+1])]
        }

        if(trim>0) {
          cat("winsorizing ... ")
          counts <<- counts/as.numeric(depth);
          inplaceWinsorizeSparseCols(counts,trim);
          counts <<- counts*as.numeric(depth);
          if(is.null(lib.sizes)) {
            depth <<- round(Matrix::rowSums(counts))
          }
        }

        counts <<- counts/(depth/depthScale);
        #counts@x <<- log10(counts@x+1)
      } else {
        stop('modelType ',modelType,' is not implemented');
      }
      misc[['rescaled.mat']] <<- NULL;
      cat("done.\n")
    },

    # adjust variance of the residual matrix, determine overdispersed sites
    adjustVariance=function(gam.k=5, alpha=5e-2, plot=FALSE, rescale.mat=TRUE, use.unadjusted.pvals=FALSE,do.par=T,max.adjusted.variance=20,min.adjusted.variance=1e-3,rowSel=NULL,verbose=TRUE,min.gene.cells=0,cells=NULL,plim=1e-13) {
      persist <- is.null(rowSel)
      if(verbose) cat("calculating variance fit ...")

      # calculate poisson residuals
      ## x <- misc[['rawCounts']];
      ## if(!is.null(batch)) {
      ##   generes <- unlist(mclapply(1:ncol(x), function(i) {
      ##     df <- data.frame('c'=x[,i],'d'=depth,'b'=batch);
      ##     m <- glm(formula=cbind(c,d-c)~d*b,data=df,family=binomial(link='logit'));
      ##     sum(abs(resid(m,type='deviance')))
      ##   },mc.cores=n.cores,mc.preschedule=T))
      ##   names(generes) <- colnames(x);
      ## } else {
      ##   generes <- unlist(mclapply(1:ncol(x), function(i) {
      ##     df <- data.frame('c'=x[,i],'d'=depth)
      ##     m <- glm(formula=cbind(c,d-c)~d,data=df,family=binomial(link='logit'));
      ##     sum(abs(resid(m,type='deviance')))
      ##   },mc.cores=n.cores,mc.preschedule=T))
      ##   names(generes) <- colnames(x);
      ## }

      # quick estimation based on raw counts and refined library parameters
      x <- misc[['rawCounts']];
      rp <- misc[['cpois']]
      # TODO: implement batch support, move to C++
      if(!is.null(batch)) {
        # p is a matrix instead of a vector, with rows specifying batches, and columns sites
        generes <- unlist(mclapply(1:ncol(x), function(i) {
          zi <- x[,i]==0;
          #sum(rp$depth[zi]*(rp$p[,i][as.integer(batch[zi])]))-sum(ppois(x[!zi,i]-1,rp$depth[!zi]*(rp$p[,i][batch[!zi]]),lower.tail=FALSE,log.p=TRUE))
          sum(rp$depth[zi]*(rp$p[,i][as.integer(batch[zi])]))-sum(ppois(x@x[seq(x@p[i]+1,x@p[i+1])]-1,rp$depth[!zi]*(rp$p[,i][batch[!zi]]),lower.tail=FALSE,log.p=TRUE))
        },mc.preschedule=T,mc.cores=n.cores))
      } else {
        generes <- unlist(mclapply(1:ncol(x), function(i) {
          zi <- x[,i]==0;
          sum(rp$depth[zi]*rp$p[i])-sum(ppois(x[!zi,i]-1,rp$depth[!zi]*rp$p[i],lower.tail=FALSE,log.p=TRUE))
        },mc.preschedule=T,mc.cores=n.cores))
      }

      misc[['deviance']] <<- generes;
      df <- colMeanVarS(counts,NULL);
      df$r <- misc[['deviance']];
      #df$m <- log10(df$m);
      df$m <- log10(Matrix::colSums(misc[['rawCounts']]>0)+1)
      #df <- data.frame(m=log10(rowSums(counts>0)+1),r=misc[['deviance']])
      rownames(df)<-colnames(counts);
      vi <- which(is.finite(df$v) & df$nobs>=min.gene.cells);
      m <- mgcv::gam(r ~ s(m, k = gam.k), data = df)
      df$d <- df$r-predict(m)

      require(VarianceGamma)
      cat(" p-values ... ")
      p <- pvg(-1*df$d,sigma=sqrt(2),nu=2/2,lower.tail=T,log.p=F)
      df$p <- p;
      p[p<plim] <- plim; p[p> 1-plim] <- 1-plim
      df$lp <- log(p)
      n.col <- ncol(counts)/2
      df$lpa <- misc[['varp']] <<- p.adjust(p);
      ods <- which(misc[['varp']]<=alpha);
      cat(length(ods),' significant sites ... ' )
      misc[['devsites']] <<- rownames(df)[ods];



      n.col <- nrow(counts)/2; # crude adjustment for dropout
      df$qv <- qchisq(p, n.col-1, lower.tail = FALSE)/n.col

      # rescale mat variance
      if(rescale.mat) {
        cat("rescaling variance ... ")
        if(!is.null(misc[['rescaled.mat']]) &&  misc[['rescaled.mat']]) { stop('adjusted variance rescaling has already been performed!')}

        if(!is.null(misc[['rescaled.mat']])) {
          if(!persist) { stop("encountered rescaled matrix in non-persist state!") }
          error('adjusted variance rescaling has already been performed!');
        }
        df$gsf <- geneScaleFactors <- sqrt(pmax(min.adjusted.variance,pmin(max.adjusted.variance,df$qv))/exp(df$v));
        inplaceColMult(counts,geneScaleFactors,rowSel);  # normalize variance of each gene
        #inplaceColMult(counts,rep(1/mean(Matrix::colSums(counts)),ncol(counts))); # normalize the column sums to be around 1
        if(persist) misc[['rescaled.mat']] <<- geneScaleFactors;
      }
      if(plot) {
        if(do.par) {
          par(mfrow=c(1,2), mar = c(3.5,3.5,2.0,0.5), mgp = c(2,0.65,0), cex = 1.0);
        }
        smoothScatter(df$m,df$r,main='unadjusted',xlab='log10[ magnitude ]',ylab='deviance')
        grid <- seq(min(df$m[vi]),max(df$m[vi]),length.out=1000)
        lines(grid,predict(m,newdata=data.frame(m=grid)),col="blue")
        if(length(ods)>0) {
          points(df$m[ods],df$r[ods],pch='.',col=2,cex=1)
        }
        smoothScatter(df$m[vi],df$qv[vi],xlab='log10[ magnitude ]',ylab='adjusted deviance',main='adjusted')
        abline(h=1,lty=2,col=8)
        if(is.finite(max.adjusted.variance)) { abline(h=max.adjusted.variance,lty=2,col=1) }
        points(df$m[ods],df$qv[ods],col=2,pch='.')
      }

      misc[['varinfo']] <<- df;

      if(verbose) cat("done.\n")
      return(invisible(df));
    },
    # make a Knn graph
    # note: for reproducibility, set.seed() and set n.cores=1
    makeKnnGraph=function(k=30,nrand=1e3,type='counts',weight.type='cauchy',odgenes=NULL,n.cores=.self$n.cores,distance='cosine',center=TRUE,x=NULL,verbose=TRUE,p=NULL) {
      require(igraph)
      if(is.null(x)) {
        x.was.given <- FALSE;
        if(type=='counts') {
          x <- counts;
        } else {
          if(type %in% names(reductions)) {
            x <- reductions[[type]];
          }
        }
        if(!is.null(odgenes)) {
          if(!all(odgenes %in% rownames(x))) { wraning("not all of the provided odgenes are present in the selected matrix")}
          if(verbose) cat("using provided odgenes ... ")
          x <- x[,odgenes]
        }

      } else {
        x.was.given <- TRUE;
      }

      # TODO: enable sparse matrix support for hnsKnn2

      if(distance=='cosine') {
        if(center) {
          x<- x - Matrix::rowMeans(x) # centering for consine distance
        }
        xn <- hnswKnn2(x,k,nThreads=n.cores,verbose=verbose)
      } else if(distance=='JS') {
        x <- x/pmax(1,Matrix::rowSums(x));
        xn <- hnswKnnJS(x,k,nThreads=n.cores)
      } else if(distance=='L2') {
        xn <- hnswKnnLp(x,k,nThreads=n.cores,p=2.0,verbose=verbose)
      } else if(distance=='L1') {
        xn <- hnswKnnLp(x,k,nThreads=n.cores,p=1.0,verbose=verbose)
      } else if(distance=='Lp') {
        if(is.null(p)) stop("p argument must be provided when using Lp distance")
        xn <- hnswKnnLp(x,k,nThreads=n.cores,p=p,verbose=verbose)
      } else {
        stop("unknown distance measure specified")
      }
      xn <- xn[!xn$s==xn$e,]

      if(n.cores==1) { # for reproducibility, sort by node names
        if(verbose) cat("ordering neighbors for reproducibility ... ");
        xn <- xn[order(xn$s+xn$e),]
        if(verbose) cat("done\n");
      }
      df <- data.frame(from=rownames(x)[xn$s+1],to=rownames(x)[xn$e+1],weight=xn$d,stringsAsFactors=F)
      if(weight.type %in% c("cauchy","normal") && ncol(x)>sqrt(nrand)) {
        # generate some random pair data for scaling
        if(distance=='cosine') {
          #rd <- na.omit(apply(cbind(sample(colnames(x),nrand,replace=T),sample(colnames(x),nrand,replace=T)),1,function(z) if(z[1]==z[2]) {return(NA); } else {1-cor(x[,z[1]],x[,z[2]])}))
          rd <- na.omit(apply(cbind(sample(colnames(x),nrand,replace=T),sample(colnames(x),nrand,replace=T)),1,function(z) if(z[1]==z[2]) {return(NA); } else {1-sum(x[,z[1]]*x[,z[2]])/sqrt(sum(x[,z[1]]^2)*sum(x[,z[2]]^2))}))
        } else if(distance=='JS') {
          rd <- na.omit(apply(cbind(sample(colnames(x),nrand,replace=T),sample(colnames(x),nrand,replace=T)),1,function(z) if(z[1]==z[2]) {return(NA); } else {jw.disR(x[,z[1]],x[,z[2]])}))
        } else if(distance=='L2') {
          rd <- na.omit(apply(cbind(sample(colnames(x),nrand,replace=T),sample(colnames(x),nrand,replace=T)),1,function(z) if(z[1]==z[2]) {return(NA); } else {sqrt(sum((x[,z[1]]-x[,z[2]])^2))}))
        } else if(distance=='L1') {
          rd <- na.omit(apply(cbind(sample(colnames(x),nrand,replace=T),sample(colnames(x),nrand,replace=T)),1,function(z) if(z[1]==z[2]) {return(NA); } else {sum(abs(x[,z[1]]-x[,z[2]]))}))
        }
        suppressWarnings(rd.model <- fitdistr(rd,weight.type))
        if(weight.type=='cauchy') {
          df$weight <- 1/pcauchy(df$weight,location=rd.model$estimate['location'],scale=rd.model$estimate['scale'])-1
        } else {
          df$weight <- 1/pnorm(df$weight,mean=rd.model$estimate['mean'],sd=rd.model$estimate['sd'])-1
        }
      }
      df$weight <- pmax(0,df$weight);
      # make a weighted edge matrix for the largeVis as well
      if(x.was.given) {
        return(invisible(as.undirected(graph.data.frame(df))))
      } else {
        misc[['edgeMat']][[type]] <<- cbind(xn,rd=df$weight);
        g <- as.undirected(graph.data.frame(df))
        graphs[[type]] <<- g;
      }
    },
    # calculate KNN-based clusters
    getKnnClusters=function(type='counts',method=fastgreedy.community, name='community', test.stability=FALSE, subsampling.rate=0.8, n.subsamplings=10, cluster.stability.threshold=0.95, n.cores=.self$n.cores, g=NULL, metaclustering.method='ward.D', min.cluster.size=2, persist=TRUE, plot=FALSE, return.details=FALSE, ...) {
      if(is.null(g)) {
        if(is.null(graphs[[type]])) { stop("call makeKnnGraph(type='",type,"', ...) first")}
        g <- graphs[[type]];
      }

      if(is.null(method)) {
        if(length(vcount(g))<2000) {
          method <- infomap.community;
        } else {
          method <- multilevel.community;
        }
      }

      # method <- multilevel.community; n.cores <- 20
      # n.subsamplings <- 10; cluster.stability.dilution <- 1.5; cluster.stability.fraction <- 0.9; subsampling.rate <- 0.8; metaclustering.method<- 'ward.D'
      # g <- r$graphs$PCA
      # x <- r$counts
      # cls <- method(g)

      #library(parallel)
      #x <- mclapply(1:5,function(z) method(g),mc.cores=20)

      cls <- method(g,...)
      cls.groups <- as.factor(membership(cls));


      if(test.stability) {
        # cleanup the clusters to remove very small ones
        cn <- names(cls.groups);
        vg <- which(unlist(tapply(cls.groups,cls.groups,length))>=min.cluster.size);
        cls.groups <- as.integer(cls.groups); cls.groups[!cls.groups %in% vg] <- NA;
        cls.groups <- as.factor(cls.groups);
        names(cls.groups) <- cn;
        # is there more than one cluster?
        if(length(levels(cls.groups))>1) {
          # run subsamplings
          cls.cells <- tapply(1:length(cls.groups),cls.groups,I)

          # if(type=='counts') {
          #   x <- counts;
          # } else {
          #   if(!type %in% names(reductions)) { stop("reduction ",type,' not found')}
          #   x <- reductions[[type]]
          # }
          #x <- counts;
          hcd <- multi2dend(cls,misc[['rawCounts']])
          m1 <- cldend2array(hcd);

          ai <- do.call(cbind,mclapply(1:n.subsamplings,function(i) {
            sg <- g;
            vi <- sample(1:vcount(sg),round(vcount(sg)*(1-subsampling.rate)))
            sg <- delete.vertices(sg,vi)
            scls <- method(sg)
            m2 <- cldend2array(multi2dend(scls,misc[['rawCounts']]))
            m1s <- m1[,colnames(m2)]
            ai <- (m1s %*% t(m2));
            ai <- ai/(outer(Matrix::rowSums(m1s),Matrix::rowSums(m2),"+") - ai)
            ns <- apply(ai,1,max)
          },mc.cores=n.cores))
          stevl <- apply(ai,1,mean); # node stability measure

          require(dendextend)
          hcd <- hcd %>% set("nodes_pch",19) %>% set("nodes_cex",3) %>% set("nodes_col",val2col(stevl,zlim=c(0.9,1)))


          # annotate n cells on the dednrogram
          t.find.biggest.stable.split <- function(l,env=environment()) {
            if(is.leaf(l)) { return(FALSE) } # don't report stable leafs ?
            bss <- mget("biggest.stable.split.size",envir=env,ifnotfound=-1)[[1]]
            if(attr(l,'nCells') <= bss) { return(FALSE) }

            # test current split for stability
            if(min(stevl[unlist(lapply(l,attr,'nodeId'))]) >= cluster.stability.threshold) { # stable
              # record size
              assign("biggest.stable.split.size",attr(l,'nCells'),envir=env)
              assign("biggest.stable.split",l,envir=env)
              return(TRUE);
            } else {
              # look within
              #return(na.omit(unlist(c(t.find.biggest.stable.split,env=env),recursive=F)))
              return(lapply(l,t.find.biggest.stable.split,env=env))
            }
          }
          # find biggest stable cell split
          e <- environment()
          bss.found <- any(unlist(t.find.biggest.stable.split(hcd,e)));
          if(bss.found) {
            # a stable split was found
            bss <- get('biggest.stable.split',envir=e)
            bss.par <- attr(bss,'nodesPar'); bss.par$col <- 'blue'; attr(bss,'nodesPar') <- bss.par;

            # find all untinterrupted stable subsplits
            consecutiveStableSubleafs <- function(l) {
              if(is.leaf(l) || min(stevl[unlist(lapply(l,attr,'nodeId'))]) < cluster.stability.threshold) {
                # either leaf or not sitting on top of a stable split - return own factor
                return(paste(unlist(l),collapse='+'));
              } else {
                # if both children are stable, return combination of their returns
                return(lapply(l,consecutiveStableSubleafs))
              }
            }
            stable.clusters <- unlist(consecutiveStableSubleafs(bss))


            final.groups <- rep("other", length(cls.groups));
            cf <- rep("other",length(cls.cells));
            for(k in stable.clusters) {
              ci <- as.integer(unlist(strsplit(k,'\\+')));
              final.groups[unlist(cls.cells[ci])] <- k;
              cf[ci] <- k;
            }
            final.groups <- as.factor(final.groups); names(final.groups) <- names(cls.groups);

            hcd <- hcd %>% branches_attr_by_clusters(clusters=as.integer(as.factor(cf[order.dendrogram(hcd)]))) %>% set("branches_lwd", 3)

          } else {
            # TODO: check for any stable node, report that

            final.groups <- rep(1, length(cls.groups)); names(final.groups) <- names(cls.groups);
          }

          #TODO: clean up cell-cluster assignment based on the pooled cluster data

          if(plot) {
            #.self$plotEmbedding(type='PCA',groups=cls.groups,show.legend=T)
            par(mar = c(3.5,3.5,2.0,0.5), mgp = c(2,0.65,0), cex = 1.0);
            hcd %>% plot();
            z <- get_nodes_xy(hcd); text(z,labels=round(stevl,2),adj=c(0.4,0.4),cex=0.7)
          }

          # return details

          rl <- list(groups=final.groups,"original.groups"=cls.groups,'hcd'=hcd,"stevl"=stevl,"cls"=cls)
          if(persist) {
            clusters[[type]][[name]] <<- final.groups;
            misc[['community']][[type]][[name]] <<- rl;
          }
          if(return.details) { return(invisible(rl))}
          return(invisible(final.groups))
        } # end: more than one group
        return(invisible(cls.groups))
      }
      if(persist) {
        clusters[[type]][[name]] <<- cls.groups;
        misc[['community']][[type]][[name]] <<- cls;
      }
      return(invisible(cls))
    },

    # calculate density-based clusters
    getDensityClusters=function(type='counts', embeddingType=NULL, name='density', v=0.7, s=1, ...) {
      if(is.null(embeddings[[type]])) { stop("first, generate embeddings for type ",type)}
      if(is.null(embeddingType)) {
        # take the first one
        embeddingType <- names(embeddings[[type]])[1]
        cat("using",embeddingType,"embedding\n")
        emb <- embeddings[[type]][[embeddingType]]

      } else {
        emb <- embeddings[[type]][[embeddingType]]
        if(is.null(emb)) { stop("embedding ",embeddingType," for type ", type," doesn't exist")}
      }
      cl <- dbscan::dbscan(emb, ...)$cluster;
      cols <- rainbow(length(unique(cl)),v=v,s=s)[cl+1];    cols[cl==0] <- "gray70"
      names(cols) <- rownames(emb);
      clusters[[type]][[name]] <<- cols;
      misc[['clusters']][[type]][[name]] <<- cols;
      return(invisible(cols))
    },
    # determine subpopulation-specific genes
    getDifferentialGenes=function(type='counts',clusterType=NULL,groups=NULL,testType='Fisher',name='customClustering', z.threshold=3,correct.z.scores=TRUE, upregulated.only=FALSE,long.form=FALSE) {
      if(is.null(groups)) {
        # look up the clustering based on a specified type
        if(is.null(clusterType)) {
          # take the first one
          cols <- clusters[[type]][[1]]
        } else {
          cols <- clusters[[type]][[clusterType]]
          if(is.null(cols)) { stop("clustering ",clusterType," for type ", type," doesn't exist")}
        }
        cols <- as.factor(cols);
      } else {
        # use clusters information
        if(!all(rownames(counts) %in% names(groups))) { warning("provided cluster vector doesn't list groups for all of the cells")}
        cols <- as.factor(groups);
      }
      cat("running differential expression with ",length(levels(cols))," clusters ... ")
      # use offsets based on the base model

      if(testType=='Fisher') {
        cm <- misc[['rawCounts']];
      } else {
        cm <- counts;
      }
      if(!all(rownames(cm) %in% names(cols))) { warning("cluster vector doesn't specify groups for all of the cells, dropping missing cells from comparison")}
      # determine a subset of cells that's in the cols and cols[cell]!=NA
      valid.cells <- rownames(cm) %in% names(cols)[!is.na(cols)];
      if(!all(valid.cells)) {
        # take a subset of the count matrix
        cm <- cm[valid.cells,]
      }
      # reorder cols
      cols <- as.factor(cols[match(rownames(cm),names(cols))]);

      if(testType=='Fisher') {
        # Fisher's exact test
        # total number of gene observations
        pickedN <- Matrix::colSums(cm)
        # number of gene observations within the group
        pickedRed <- colSumByFac(cm,as.integer(cols))[-1,]
        # total depth within each group
        redN <- Matrix::rowSums(pickedRed);
        # total depth within other cells
        whiteN <- sum(redN)-redN;
        lpv <- do.call(rbind,lapply(1:nrow(pickedRed),function(i) {
          phyper(pickedRed[i,],redN[i],whiteN[i],pickedN+1,lower.tail=FALSE,log.p=TRUE)
        }))
        lpv[pickedRed==0] <- 0;
        lower.lpv.limit <- -100;
        lpv[lpv<lower.lpv.limit] <- lower.lpv.limit;
        lpv <- -1*lpv;
        if(!upregulated.only) {
          llpv <- do.call(rbind,lapply(1:nrow(pickedRed),function(i) {
            phyper(pickedRed[i,],redN[i],whiteN[i],pickedN,lower.tail=TRUE,log.p=TRUE)
          }))
          llpv[llpv<lower.lpv.limit] <- lower.lpv.limit;
          di <- llpv < -1*lpv;
          lpv[di] <- llpv[di];
        }
        x <- t(lpv);
      } else {
        # Wilcoxon rank test
        # calculate rank per-column (per-gene) average rank matrix
        xr <- sparse_matrix_column_ranks(cm);
        # calculate rank sums per group
        grs <- colSumByFac(xr,as.integer(cols))[-1,]
        # calculate number of non-zero entries per group
        xr@x <- numeric(length(xr@x))+1
        gnzz <- colSumByFac(xr,as.integer(cols))[-1,]
        #group.size <- as.numeric(tapply(cols,cols,length));
        group.size <- as.numeric(tapply(cols,cols,length))[1:nrow(gnzz)]; group.size[is.na(group.size)]<-0; # trailing empty levels are cut off by colSumByFac
        # add contribution of zero entries to the grs
        gnz <- (group.size-gnzz)
        # rank of a 0 entry for each gene
        zero.ranks <- (nrow(xr)-diff(xr@p)+1)/2 # number of total zero entries per gene
        ustat <- t((t(gnz)*zero.ranks)) + grs - group.size*(group.size+1)/2
        # standardize
        n1n2 <- group.size*(nrow(cm)-group.size);
        # usigma <- sqrt(n1n2*(nrow(cm)+1)/12) # without tie correction
        # correcting for 0 ties, of which there are plenty
        usigma <- sqrt(n1n2*(nrow(cm)+1)/12)
        usigma <- sqrt((nrow(cm) +1 - (gnz^3 - gnz)/(nrow(cm)*(nrow(cm)-1)))*n1n2/12)
        x <- t((ustat - n1n2/2)/usigma); # standardized U value- z score
      }


      ## # run quick fisher test on the count matrix for each group
      ## lower.lpv.limit <- -100;
      ## x <- do.call(cbind,tapply(1:nrow(cm),cols,function(ii) {
      ##   # total number of matching, non-matching rows
      ##   redN <- length(ii); whiteN <- length(cols)-redN
      ##   pickedN <- diff(cm@p)
      ##   pickedRed <- diff(cm[ii,]@p)
      ##   lpv <- -1*pmax(lower.lpv.limit,phyper(pickedRed,redN,whiteN,pickedN+1,lower.tail=FALSE,log.p=TRUE))
      ##   if(!upregulated.only) {
      ##     lpvl <- pmax(lower.lpv.limit,phyper(pickedRed,redN,whiteN,pickedN,lower.tail=TRUE,log.p=TRUE))
      ##     di <- lpvl < -1*lpv;
      ##     lpv[di] <- lpvl[di];
      ##   }
      ##   lpv
      ## }))

      # TODO: batch correction


      # correct for multiple hypothesis
      cat("adjusting p-values ... ")
      if(correct.z.scores) {
        x <- matrix(qnorm(scde:::bh.adjust(pnorm(as.numeric(abs(x)), lower.tail = FALSE, log.p = TRUE), log = TRUE), lower.tail = FALSE, log.p = TRUE),ncol=ncol(x))*sign(x)
      }
      rownames(x) <- colnames(cm); colnames(x) <- levels(cols)[1:ncol(x)];
      cat("done.\n")
      if(upregulated.only) {
        ds <- apply(x,2,function(z) {vi <- which(z>=z.threshold); r <- z[vi]; names(r) <- rownames(x)[vi]; sort(r,decreasing=T)})
      } else {
        ds <- apply(x,2,function(z) {vi <- which(abs(z)>=z.threshold); r <- z[vi]; names(r) <- rownames(x)[vi]; sort(r,decreasing=T)})
      }
      if(is.null(groups)) {
        if(is.null(clusterType)) {
          diffgenes[[type]][[names(clusters[[type]])[1]]] <<- ds;
        } else {
          diffgenes[[type]][[clusterType]] <<- ds;
        }
      } else {
        diffgenes[[type]][[name]] <<- ds;
      }
      return(invisible(ds))
    },
    plotDiffGeneHeatmap=function(type='mat',clusterType=NULL, groups=NULL, n.genes=100, z.score=2, gradient.range.quantile=0.95, inner.clustering=FALSE, gradientPalette=NULL, v=0.8, s=1, box=TRUE, ... ) {
      if(!is.null(clusterType)) {
        x <- diffgenes[[type]][[clusterType]];
        if(is.null(x)) { stop("differential genes for the specified cluster type haven't been calculated") }
      } else {
        x <- diffgenes[[type]][[1]];
        if(is.null(x)) { stop("no differential genes found for data type ",type) }
      }

      if(is.null(groups)) {
        # look up the clustering based on a specified type
        if(is.null(clusterType)) {
          # take the first one
          cols <- clusters[[type]][[1]]
        } else {
          cols <- clusters[[type]][[clusterType]]
          if(is.null(cols)) { stop("clustering ",clusterType," for type ", type," doesn't exist")}
        }
      } else {
        # use clusters information
        if(!all(colnames(counts) %in% names(groups))) { warning("provided cluster vector doesn't list groups for all of the cells")}
        cols <- as.factor(groups[match(colnames(counts),names(groups))]);
      }
      cols <- as.factor(cols);
      # select genes to show
      if(!is.null(z.score)) {
        x <- lapply(x,function(zv) zv[zv>=z.score])
        if(!is.null(n.genes)) {
          x <- lapply(x,function(zv) { if(length(zv)>0) { zv <- zv[1:min(length(zv),n.genes)] }; zv })
        }
      } else {
        if(!is.null(n.genes)) {
          x <- lapply(x,function(zv) { if(length(zv)>0) { zv <- zv[1:min(length(zv),n.genes)] }; zv })
        }
      }
      x <- lapply(x,names);
      # make expression matrix
      em <- mat[unlist(x),];

      # renormalize rows
      if(all(sign(em)>=0)) {
        if(is.null(gradientPalette)) {
          gradientPalette <- colorRampPalette(c('gray90','red'), space = "Lab")(1024)
        }
        em <- t(apply(em,1,function(x) {
          zlim <- as.numeric(quantile(x,p=c(1-gradient.range.quantile,gradient.range.quantile)))
          if(diff(zlim)==0) {
            zlim <- as.numeric(range(x))
          }
          x[x<zlim[1]] <- zlim[1]; x[x>zlim[2]] <- zlim[2];
          x <- (x-zlim[1])/(zlim[2]-zlim[1])
        }))
      } else {
        if(is.null(gradientPalette)) {
          gradientPalette <- colorRampPalette(c("blue", "grey90", "red"), space = "Lab")(1024)
        }
        em <- t(apply(em,1,function(x) {
          zlim <- c(-1,1)*as.numeric(quantile(abs(x),p=gradient.range.quantile))
          if(diff(zlim)==0) {
            zlim <- c(-1,1)*as.numeric(max(abs(x)))
          }
          x[x<zlim[1]] <- zlim[1]; x[x>zlim[2]] <- zlim[2];
          x <- (x-zlim[1])/(zlim[2]-zlim[1])
        }))
      }

      # cluster cell types by averages
      rowfac <- factor(rep(names(x),unlist(lapply(x,length))),levels=names(x))
      if(inner.clustering) {
        clclo <- hclust(as.dist(1-cor(do.call(cbind,tapply(1:nrow(em),rowfac,function(ii) Matrix::colMeans(em[ii,,drop=FALSE]))))),method='complete')$order
      } else {
        clclo <- 1:length(levels(rowfac))
      }

      if(inner.clustering) {
        # cluster genes within each cluster
        clgo <- tapply(1:nrow(em),rowfac,function(ii) {
          ii[hclust(as.dist(1-cor(t(em[ii,]))),method='complete')$order]
        })
      } else {
        clgo <- tapply(1:nrow(em),rowfac,I)
      }
      if(inner.clustering) {
        # cluster cells within each cluster
        clco <- tapply(1:ncol(em),cols,function(ii) {
          ii[hclust(as.dist(1-cor(em[,ii])),method='complete')$order]
        })
      } else {
        clco <- tapply(1:ncol(em),cols,I)
      }

      cellcols <- fac2col(cols,v=v,s=s)[unlist(clco[clclo])]
      genecols <- rev(rep(fac2col(cols,v=v,s=s,return.level.colors=T),unlist(lapply(clgo,length)[clclo])))
      my.heatmap2(em[rev(unlist(clgo[clclo])),unlist(clco[clclo])],col=gradientPalette,Colv=NA,Rowv=NA,labRow=NA,labCol=NA,RowSideColors=genecols,ColSideColors=cellcols,margins=c(0.5,0.5),ColSideColors.unit.vsize=0.05,RowSideColors.hsize=0.05,useRaster=T, box=box, ...)
      abline(v=cumsum(unlist(lapply(clco[clclo],length))),col=1,lty=3)
      abline(h=cumsum(rev(unlist(lapply(clgo[clclo],length)))),col=1,lty=3)
    },

    # recalculate library sizes using robust regression within clusters
    getRefinedLibSizes=function( ..., rescale.counts=TRUE ) {
      # TODO: add support for groups
      # TODO: implement using RcppGSL
      cpois.em <- function(d,maxit=10,tol=1e-3,n.cores=30,verbose=FALSE) {
        # assumes binary d

        # derivative of the poisson likelihood
        # x - gene observations (in a given cell)
        # p - underlying gene probabilities
        dLdn <- function(n,x,p) {
          zi <- x==0;
          sum(p[!zi]/(exp(n*p[!zi])-1))-sum(p[zi])
        }
        # initial guesses
        dep <- Matrix::rowSums(d);
        p <- Matrix::colSums(d)/sum(dep); p <- p/sum(p);
        cs <- Matrix::colSums(d>0);
        maxp <- cs==nrow(d)
        minp <- cs==0
        vi <- which(!maxp & !minp);
        maxp <- which(maxp); minp <- which(minp);
        for(iteration in 1:maxit) {
          # update depth
          cat("iteration",iteration," ")
          new.dep <- unlist(mclapply(1:nrow(d),function(i) uniroot(dLdn,c(1/10,10)*dep[i],x=d[i,],p=p,extendInt="yes",tol=1e-10)$root,mc.cores=n.cores,mc.preschedule=T))
          dep.tol <- abs(dep-new.dep)/new.dep;
          cat("E (tol=",mean(dep.tol),") ");
          new.pp <- unlist(mclapply(vi,function(i) uniroot(dLdn,c(1e-10,1),x=d[,i],p=new.dep,tol=1e-10,extendInt='yes')$root,mc.cores=30,mc.preschedule=T))

          new.p <- p; new.p[vi] <- new.pp;
          new.p[maxp] <- 1-1e-3;
          new.p[minp] <- 1e-10;
          new.p <- new.p/sum(new.p);
          p.tol <- abs(p-new.p)/new.p;
          cat("M (tol=",mean(new.p),") \n");
          #if(mean(dep.tol)<=tol & mean(p.tol)<=tol) { break; }
          if(mean(dep.tol)<=tol) { break; }
          dep <- new.dep; p <- new.p;
        }
        names(p) <- colnames(d); names(dep) <- rownames(d);
        return(list(p=p,depth=round(dep)));
      }
      if(is.null(batch)) {
        x <- cpois.em(misc[['rawCounts']],...)
      } else {
        x <- tapply(1:nrow(misc[['rawCounts']]),batch,function(ii) {
          cat("batch ",levels(batch)[batch[ii[1]]],":\n")
          cpois.em((misc[['rawCounts']])[ii,],...)
        })

        # merge into a common structure: single depth vector, and a matrix of p (rows - batch, columns - sites)
        cdepth <- unlist(lapply(x,function(z) z$depth));
        names(cdepth) <- unlist(lapply(x,function(z) names(z$depth)))
        cdepth <- cdepth[rownames(misc[['rawCounts']])]
        x <- list(p=do.call(rbind,lapply(x,function(z) z$p)),depth=cdepth)
      }

      misc[['cpois']] <<- x;
      if(rescale.counts) {
        lib.sizes <- x$depth[rownames(misc[['rawCounts']])];
        if(modelType=='plain') {
          cat("rescaling matrix")
          counts <<- counts * (depth/lib.sizes);
          cat(" done\n")
        }
        depth <<- lib.sizes;
      }
      return(invisible(x));
    },

    # plot heatmap for a given set of genes
    plotGeneHeatmap=function(genes, type='mat', clusterType=NULL, groups=NULL, z.score=2, gradient.range.quantile=0.95, cluster.genes=FALSE, inner.clustering=FALSE, gradientPalette=NULL, v=0.8, s=1, box=TRUE, ... ) {

      if(is.null(groups)) {
        # look up the clustering based on a specified type
        if(is.null(clusterType)) {
          # take the first one
          cols <- clusters[[type]][[1]]
        } else {
          cols <- clusters[[type]][[clusterType]]
          if(is.null(cols)) { stop("clustering ",clusterType," for type ", type," doesn't exist")}
        }
      } else {
        # use clusters information
        if(!all(colnames(counts) %in% names(groups))) { warning("provided cluster vector doesn't list groups for all of the cells")}
        cols <- as.factor(groups[match(colnames(counts),names(groups))]);
      }
      cols <- as.factor(cols);
      # make expression matrix
      if(!all(genes %in% rownames(mat))) { warning(paste("the following specified genes were not found in the data: [",paste(genes[!genes %in% rownames(mat)],collapse=" "),"], omitting",sep="")) }
      x <- intersect(genes,rownames(mat));
      if(length(x)<1) { stop("too few genes") }
      em <- mat[x,];

      # renormalize rows
      if(all(sign(em)>=0)) {
        if(is.null(gradientPalette)) {
          gradientPalette <- colorRampPalette(c('gray90','red'), space = "Lab")(1024)
        }
        em <- t(apply(em,1,function(x) {
          zlim <- as.numeric(quantile(x,p=c(1-gradient.range.quantile,gradient.range.quantile)))
          if(diff(zlim)==0) {
            zlim <- as.numeric(range(x))
          }
          x[x<zlim[1]] <- zlim[1]; x[x>zlim[2]] <- zlim[2];
          x <- (x-zlim[1])/(zlim[2]-zlim[1])
        }))
      } else {
        if(is.null(gradientPalette)) {
          gradientPalette <- colorRampPalette(c("blue", "grey90", "red"), space = "Lab")(1024)
        }
        em <- t(apply(em,1,function(x) {
          zlim <- c(-1,1)*as.numeric(quantile(abs(x),p=gradient.range.quantile))
          if(diff(zlim)==0) {
            zlim <- c(-1,1)*as.numeric(max(abs(x)))
          }
          x[x<zlim[1]] <- zlim[1]; x[x>zlim[2]] <- zlim[2];
          x <- (x-zlim[1])/(zlim[2]-zlim[1])
        }))
      }

      # cluster cell types by averages
      clclo <- 1:length(levels(cols))

      if(cluster.genes) {
        # cluster genes within each cluster
        clgo <- hclust(as.dist(1-cor(t(em))),method='complete')$order
      } else {
        clgo <- 1:nrow(em)
      }

      if(inner.clustering) {
        # cluster cells within each cluster
        clco <- tapply(1:ncol(em),cols,function(ii) {
          ii[hclust(as.dist(1-cor(em[,ii])),method='complete')$order]
        })
      } else {
        clco <- tapply(1:ncol(em),cols,I)
      }

      cellcols <- fac2col(cols,v=v,s=s)[unlist(clco[clclo])]
      #genecols <- rev(rep(fac2col(cols,v=v,s=s,return.level.colors=T),unlist(lapply(clgo,length)[clclo])))
      my.heatmap2(em[rev(clgo),unlist(clco[clclo])],col=gradientPalette,Colv=NA,Rowv=NA,labCol=NA,ColSideColors=cellcols,margins=c(0.5,5),ColSideColors.unit.vsize=0.05,RowSideColors.hsize=0.05,useRaster=T, box=box, ...)
      abline(v=cumsum(unlist(lapply(clco[clclo],length))),col=1,lty=3)
      #abline(h=cumsum(rev(unlist(lapply(clgo[clclo],length)))),col=1,lty=3)
    },

    # show embedding
    plotEmbedding=function(type='counts', embeddingType=NULL, clusterType=NULL, groups=NULL, colors=NULL, do.par=T, cex=0.6, alpha=0.4, gradientPalette=NULL, zlim=NULL, s=1, v=0.8, min.group.size=1, show.legend=FALSE, mark.clusters=FALSE, mark.cluster.cex=2, shuffle.colors=F, legend.x='topright', gradient.range.quantile=0.95, quiet=F, unclassified.cell.color='gray70', ...) {
      if(is.null(embeddings[[type]])) { stop("first, generate embeddings for type ",type)}
      if(is.null(embeddingType)) {
        # take the first one
        emb <- embeddings[[type]][[1]]
      } else {
        emb <- embeddings[[type]][[embeddingType]]
      }
      factor.mapping=FALSE;
      if(is.null(colors) && is.null(groups)) {
        # look up the clustering based on a specified type
        if(is.null(clusterType)) {
          # take the first one
          groups <- clusters[[type]][[1]]
        } else {
          groups <- clusters[[type]][[clusterType]]
          if(is.null(groups)) { stop("clustering ",clusterType," for type ", type," doesn't exist")}
        }

        groups <- as.factor(groups[rownames(emb)]);
        if(min.group.size>1) { groups[groups %in% levels(groups)[unlist(tapply(groups,groups,length))<min.group.size]] <- NA; groups <- as.factor(groups); }
        cols <- fac2col(groups,s=s,v=v,shuffle=shuffle.colors,min.group.size=min.group.size)[rownames(emb)]
        factor.mapping=TRUE;
      } else {
        if(!is.null(colors)) {
          # use clusters information
          if(!all(rownames(emb) %in% names(colors))) { warning("provided cluster vector doesn't list colors for all of the cells; unmatched cells will be shown in gray. ")}
          if(all(areColors(colors))) {
            if(!quiet) cat("using supplied colors as is\n")
            cols <- colors[match(rownames(emb),names(colors))]; cols[is.na(cols)] <- unclassified.cell.color;
          } else {
            cols <- val2col(colors,gradient.range.quantile=gradient.range.quantile); cols[is.na(cols)] <- unclassified.cell.color;
            ## if(is.numeric(colors)) { # treat as a gradient
            ##   if(!quiet) cat("treating colors as a gradient")
            ##   if(is.null(gradientPalette)) { # set up default gradients
            ##     if(all(sign(colors)>=0)) {
            ##       gradientPalette <- colorRampPalette(c('gray80','red'), space = "Lab")(1024)
            ##     } else {
            ##       gradientPalette <- colorRampPalette(c("blue", "grey70", "red"), space = "Lab")(1024)
            ##     }
            ##   }
            ##   if(is.null(zlim)) { # set up value limits
            ##     if(all(sign(colors)>=0)) {
            ##       zlim <- as.numeric(quantile(colors,p=c(1-gradient.range.quantile,gradient.range.quantile),na.rm=T))
            ##       if(diff(zlim)==0) {
            ##         zlim <- as.numeric(range(colors))
            ##       }
            ##     } else {
            ##       zlim <- c(-1,1)*as.numeric(quantile(abs(colors),p=gradient.range.quantile,na.rm=T))
            ##       if(diff(zlim)==0) {
            ##         zlim <- c(-1,1)*as.numeric(max(abs(colors)))
            ##       }
            ##     }
            ##   }
            ##   # restrict the values
            ##   colors[colors<zlim[1]] <- zlim[1]; colors[colors>zlim[2]] <- zlim[2];

            ##   if(!quiet) cat(' with zlim:',zlim,'\n')
            ##   colors <- (colors-zlim[1])/(zlim[2]-zlim[1])
            ##   cols <- gradientPalette[colors[match(rownames(emb),names(colors))]*(length(gradientPalette)-1)+1]
            ## } else {
            ##   stop("colors argument must be a cell-named vector of either character colors or numeric values to be mapped to a gradient")
            ## }
          }
        } else {
          if(!is.null(groups)) {
            if(min.group.size>1) { groups[groups %in% levels(groups)[unlist(tapply(groups,groups,length))<min.group.size]] <- NA; groups <- as.factor(groups); }
            groups <- as.factor(groups)[rownames(emb)]
            if(!quiet) cat("using provided groups as a factor\n")
            factor.mapping=TRUE;
            # set up a rainbow color on the factor
            cols <- fac2col(groups,s=s,v=v,shuffle=shuffle.colors,min.group.size=min.group.size,unclassified.cell.color=unclassified.cell.color)
          }
        }
        names(cols) <- rownames(emb)
      }

      if(do.par) {
        par(mar = c(0.5,0.5,2.0,0.5), mgp = c(2,0.65,0), cex = 1.0);
      }
      plot(emb,col=adjustcolor(cols,alpha=alpha),cex=cex,pch=19,axes=F, panel.first=grid(lty=2,nx=5), ...); box();
      if(mark.clusters) {
        if(!is.null(groups)) {
          cent.pos <- do.call(rbind,tapply(1:nrow(emb),groups,function(ii) apply(emb[ii,,drop=F],2,median)))
          #rownames(cent.pos) <- levels(groups);
          cent.pos <- na.omit(cent.pos);
          text(cent.pos[,1],cent.pos[,2],labels=rownames(cent.pos),cex=mark.cluster.cex)
        }
      }
      if(show.legend) {
        if(factor.mapping) {
          legend(x=legend.x,pch=rep(19,length(levels(groups))),col=rainbow(length(unique(groups)),v=v,s=s),legend=levels(groups))
        }
      }

    },
    getOdGenes=function(n.odgenes=NULL,alpha=5e-2,use.unadjusted.pvals=FALSE) {
      if(is.null(misc[['varinfo']])) { stop("please run adjustVariance first")}
      if(is.null(n.odgenes)) { #return according to alpha
        if(use.unadjusted.pvals) {
          rownames(misc[['varinfo']])[misc[['varinfo']]$lp <= log(alpha)]
        } else {
          rownames(misc[['varinfo']])[misc[['varinfo']]$lpa <= log(alpha)]
        }
      } else { # return top n.odgenes sites
        rownames(misc[['varinfo']])[(order(misc[['varinfo']]$lp,decreasing=F)[1:min(ncol(counts),n.odgenes)])]
      }
    },

    recursiveClustering=function(alpha=5e-2,nPcs=30,k=50,cluster.method=NULL,n.cores=.self$n.cores,groups=NULL,cells=NULL,verbose=TRUE,use.unadjusted.pvals=FALSE,gam.k=5,max.adjusted.variance=1e3,min.adjusted.variance=1e-3,dist='L1',weight.type='none',center.knn=TRUE,current.recursion.level=1,max.recursion.level=10,min.cells=nPcs,min.cluster.size=4,label='R') {
      if(is.null(cells)) { cells <- seq(1,nrow(counts)) }

      dendr.node <- function() {
        l <- as.integer(cells[1]);
        rowFac <- rep(NA,nrow(counts)); rowFac[cells] <- 1;
        tc <- colSumByFac(counts,as.integer(rowFac))
        rownames(tc) <- c("other","leaf")
        attributes(l) <- list(members=1,rheight=0,leaf=TRUE,cells=cells,label=label,pooledClusterCounts=tc)
        return(l);
      }



      if(current.recursion.level>=max.recursion.level || length(cells)<min.cells) return(dendr.node()); # base case

      # remove scaling
      if(!is.null(misc[['rescaled.mat']])) {
        if(!is.null(cells)) { stop("scaled matrix encountered inside recursion") }
        if(length(misc[['rescaled.mat']])!=nrow(counts)) { stop("invalid geneScalingFactor state! aborting.") }
        inplaceColMult(counts,1/misc[['rescaled.mat']],NULL);
        misc[['rescaled.mat']] <<- NULL;
      }


      if(is.null(groups)) {
        # perform a clustering round

        cf <- rep(FALSE,nrow(counts)); cf[cells] <- TRUE;
        df <- .self$adjustVariance(rescale.mat=FALSE,plot=FALSE,use.unadjusted.pvals=use.unadjusted.pvals,rowSel=cf,verbose=FALSE,min.gene.cells=min.cluster.size,gam.k=gam.k)

        if(use.unadjusted.pvals) {
          ods <- which(df$lp<log(alpha))
        } else {
          ods <- which(df$lpa<log(alpha))
        }

        if(length(ods)<1) {
          if(verbose) { cat("[",length(cells),":",length(ods),":0]") }
          return(dendr.node()) # base case
        }

        # add more genes
        ods <- order(df$lp,decreasing=F)[1:(length(ods)+min(nrow(df),max(length(ods)*2,1e3)))]

        # determine subgroups
        # PCA
        x <- counts[cells,ods,drop=F];
        gsf <- sqrt(pmax(min.adjusted.variance,pmin(max.adjusted.variance,df$qv[ods]))/exp(df$v[ods]))
        inplaceColMult(x,gsf,NULL);


        npc <- pmin(nPcs,ceiling(ncol(x)/5),ceiling(nrow(x)/5))
        if(npc>nrow(x)/5 || npc>ncol(x)/5) { # no longer a large problem
          pcs <- svd(x,LINPACK=TRUE)
        } else {
          pcs <- irlba(x, nv=npc, nu=0, center=Matrix::colMeans(x), right_only=FALSE)
        }
        pcas <- as.matrix(x %*% pcs$v);
        rownames(pcas) <- rownames(x)
        colnames(pcas) <- paste('PC',seq(ncol(pcas)),sep='')

        # KNN
        g <- .self$makeKnnGraph(k=k,x=pcas,dist=dist,weight.type=weight.type,center=center.knn,n.cores=n.cores,verbose=FALSE)

        # clustering
        groups <- .self$getKnnClusters(method=cluster.method,g=g,test.stability=T,plot=T,persist=FALSE,metaclustering.method = 'ward.D')

      } else {
        # make sure the groups is a factor on cells

        if(is.null(names(groups))) { stop("starting groups must be a named factor on all cells") }
        nm <- match(rownames(counts)[cells],names(groups));
        if(any(is.na(nm))) { stop("starting groups do not explicitly cover all of the cells") }
        groups <- as.factor(groups[nm]);
        if(verbose) cat("using starting groups with",length(levels(groups)),"levels");
      }

      # cleanup the clusters to remove very small ones
      cn <- names(groups);
      vg <- which(unlist(tapply(groups,groups,length))>=min.cluster.size);
      groups <- as.integer(groups); groups[!groups %in% vg] <- NA;
      groups <- as.factor(groups);
      names(groups) <- cn;

      if(length(levels(groups))<2) {
        if(verbose) { cat("[",label,":",length(cells),":",length(ods),":",length(levels(groups)),"]") }
        return(dendr.node()) # base case
      }
      if(verbose) { cat("[",label,":",length(cells),":",length(ods),":",length(levels(groups)),"]->") }
      # recursive calls
      # TODO: future evaluation
      #children <- tapply(cells,groups,function(ci) .self$recursiveClustering(cells=ci,current.recursion.level=current.recursion.level+1,alpha=alpha,nPcs=nPcs,k=k,cluster.method=cluster.method,n.cores=n.cores,verbose=verbose,use.unadjusted.pvals=use.unadjusted.pvals,min.adjusted.variance,weight.type=weight.type,center.knn=center.knn,max.recursion.level=max.recursion.level,groups=NULL,label=paste(label,groups[ci[1]],sep='.')),simplify=FALSE)
      children <- lapply(levels(groups),function(cl) {
        .self$recursiveClustering(cells=cells[which(groups==cl)],current.recursion.level=current.recursion.level+1,alpha=alpha,nPcs=nPcs,k=k,cluster.method=cluster.method,n.cores=n.cores,verbose=verbose,gam.k=gam.k,use.unadjusted.pvals=use.unadjusted.pvals,min.adjusted.variance,dist=dist,weight.type=weight.type,center.knn=center.knn,max.recursion.level=max.recursion.level,groups=NULL,label=paste(label,cl,sep='.'))
      })
      #children <- lapply(children,value); # future attempt - didn't work
      # postprocessing
      attributes(children) <- list(members=sum(unlist(lapply(children,function(x) attr(x,'members')))),rheight=0,cells=cells,names=levels(groups),label=label)
      attr(children,'midpoint') <- attr(children,'members')/2;
      attr(children,'groups') <- groups;
      # calculate aggregate cluster profiles for the groups
      # build a global factor based on the groups
      rowFac <- rep(-1,nrow(counts)); rowFac[cells] <- as.integer(groups);
      tc <- colSumByFac(counts,as.integer(rowFac))
      rownames(tc) <- c("nonclust",levels(groups))
      tc <- rbind("total"=Matrix::colSums(tc),tc)
      attr(children,'pooledClusterCounts') <- tc;
      d <- jsDist(t(((tc/pmax(1,Matrix::rowSums(tc)))))); rownames(d) <- colnames(d) <- rownames(tc)
      #d <- 1-cor(t(tc)); rownames(d) <- colnames(d) <- rownames(tc)
      #d <- 1-cor(t(log10(tc/pmax(1,Matrix::rowSums(tc))*1e4+1))); rownames(d) <- colnames(d) <- rownames(tc)
      attr(children,'groupDist') <- d
      hc <- hclust(as.dist(d[-c(1,2),-c(1,2)]),method='ward.D')
      attr(children,'groupClust') <- hc
      # set the hights of the child nodes
      for(i in levels(groups)) {  attr(children[[i]],'rheight') <- d['total',i]  }

      if(current.recursion.level>1) return(children);

      # determine cumulative heights
      cumh <- function(l,h=0) { attr(l,'height') <-  attr(l,'rheight')+h; if(!is.leaf(l)) { x <- lapply(l,cumh,h=attr(l,'height')); attributes(x) <- attributes(l); return(x)} else { return(l)} }
      children <- cumh(children)
      # now reverse so that the lowest leaf is at 0
      maxh <- function(l) { if(is.leaf(l)) { return(attr(l,'height')) }; return(max(c(attr(l,'height'),unlist(lapply(l,maxh))))) }
      subh <- function(l,c) { attr(l,'height') <- c-attr(l,'height'); if(!is.leaf(l)) { x <- lapply(l,subh,c); attributes(x) <- attributes(l); return(x)} else { return(l) }}
      children <- subh(children,maxh(children));

      # updated node values to be compatible with as.hclust
      pe <- environment();
      assign('v',1,envir=pe);
      nvup <- function(l) {
        if(is.leaf(l)) {
          x <- get('v',envir=pe); assign('v',x+1,envir=pe); attributes(x) <- attributes(l); return(x)
        } else {
          z <- lapply(l,nvup); attributes(z) <- attributes(l); return(z);
        }
      }
      children <- nvup(children);

      class(children) <- "dendrogram";
      misc[['recursiveClusters']] <<- children;
      return(invisible(children))
    },

    getBinarizedClusterDendrogram=function(children=NULL) {
      if(is.null(children)) { children <- misc[['recursiveClusters']] }
      if(is.null(children)) stop("run recursiveClustering() first, or provide pre-existing result")
      # determine relative heights from the standard ones
      setrh <- function(l,h) {
        attr(l,'rheight') <- h-attr(l,'height');
        if(!is.leaf(l)) { x <- lapply(l,setrh,attr(l,'height')); attributes(x) <- attributes(l); return(x) } else { return(l) }
      }

      getbin <- function(l) {
        if(is.leaf(l)) { return(l) }
        d <- as.dendrogram(attr(l,'groupClust'));
        d <- setrh(d,attr(d,'height'));
        attr(d,'rheight') <- attr(l,'rheight')
        # substitute the leaves with the recursive calls
        attachChildren <- function(x) {
          if(is.leaf(x)) {
            z <- getbin(l[[attr(x,'label')]]);
            attr(z,'rheight') <- attr(x,'rheight');
            attr(z,'midpoint') <- 0.5;
            return(z)
          } else {
            z <- lapply(x,attachChildren);
            attributes(z) <- attributes(x);
            return(z);
          }
        }
        d <- attachChildren(d);
        attr(l,'names') <- attr(d,'names');
        attributes(d) <- attributes(l);
        d
      }
      children <- getbin(children);
      # calculate heights
      cumh <- function(l,h=0) { attr(l,'height') <-  attr(l,'rheight')+h; if(!is.leaf(l)) { x <- lapply(l,cumh,h=attr(l,'height')); attributes(x) <- attributes(l); return(x)} else { return(l)} }
      children <- cumh(children)
      maxh <- function(l) { if(is.leaf(l)) { return(attr(l,'height')) }; return(max(c(attr(l,'height'),unlist(lapply(l,maxh))))) }
      subh <- function(l,c) { attr(l,'height') <- c-attr(l,'height'); if(!is.leaf(l)) { x <- lapply(l,subh,c); attributes(x) <- attributes(l); return(x)} else { return(l) }}
      children <- subh(children,maxh(children));
      # adjust members
      amem <- function(l) { if(is.leaf(l)) { attr(l,'members') <- 1; return(l) } else { x <- lapply(l,amem); attributes(x) <- attributes(l); attr(x,'members') <- sum(unlist(lapply(x,attr,'members'))); return(x); }}
      attr(children,'midpoint') <- 0.5;
      children <- amem(children);
    },
    calculateRPReduction=function(nRPs=1e3,type='counts',name='RP', use.odgenes=FALSE, n.odgenes=NULL, odgenes=NULL, scale=T, center=T) {
      if(type=='counts') {
        x <- counts;
      } else {
        if(!type %in% names(reductions)) { stop("reduction ",type,' not found')}
        x <- reductions[[type]]
      }
      if((use.odgenes || !is.null(n.odgenes)) && is.null(odgenes)) {
        if(is.null(misc[['devsites']] )) { stop("please run adjustVariance() first")}
        odgenes <- misc[['devsites']];
        if(!is.null(n.odgenes)) {
          if(n.odgenes>length(odgenes)) {
            #warning("number of specified odgenes is higher than the number of the statistically significant sites, will take top ",n.odgenes,' sites')
            odgenes <- rownames(misc[['varinfo']])[(order(misc[['varinfo']]$lp,decreasing=F)[1:min(ncol(counts),n.odgenes)])]
          } else {
            odgenes <- odgenes[1:n.odgenes]
          }
        }
      }
      if(!is.null(odgenes)) { x <- x[,odgenes] }
      # TODO: sparse version

      x <- as.matrix(t(x)); x <- x-Matrix::rowMeans(x);
      wvec <- c(1,1,1)/3;
      smat <- do.call(rbind,lapply(1:nRPs,function(i) sample(c(-1,0,1),nrow(x),replace=T,prob=wvec)))
      pcas <- t(smat %*% x)

      rownames(pcas) <- colnames(x)
      colnames(pcas) <- paste('RP',seq(ncol(pcas)),sep='')

      reductions[[name]] <<- pcas;
      return(invisible(pcas))
    },

    carefullyNormalizedReduction=function(type='counts',use.odgenes=FALSE, n.odgenes=NULL, odgenes=NULL, scale=F,center=T, name='normalized',n.depth.slices=20,depthScale=1e4,sparse=FALSE) {

      if(type=='counts') {
        #x <- counts;
        x <- misc[['rawCounts']]
      } else {
        if(!type %in% names(reductions)) { stop("reduction ",type,' not found')}
        x <- reductions[[type]]
      }
      if((use.odgenes || !is.null(n.odgenes)) && is.null(odgenes)) {
        if(is.null(misc[['devsites']] )) { stop("please run adjustVariance() first")}
        odgenes <- misc[['devsites']];
        if(!is.null(n.odgenes)) {
          if(n.odgenes>length(odgenes)) {
            #warning("number of specified odgenes is higher than the number of the statistically significant sites, will take top ",n.odgenes,' sites')
            odgenes <- rownames(misc[['varinfo']])[(order(misc[['varinfo']]$lp,decreasing=F)[1:min(ncol(counts),n.odgenes)])]
          } else {
            odgenes <- odgenes[1:n.odgenes]
          }
        }
      }
      if(!is.null(odgenes)) { x <- x[,odgenes] }

      if(is.null(misc[['cpois']])) {
        warning("using naive guesses for library size/abundance. getRefinedLibSizes() is recommended")
        dep <- Matrix::rowSums(x); names(dep) <- rownames(x);
        if(is.null(batch)) {
          p <- Matrix::colSums(x)/sum(dep); p <- p/sum(p); names(p) <- colnames(x)
        } else {
          # in this case, p is a matrix with the number of rows equal to the number of batches
          p <- colSumByFac(x,as.integer(batch))[-1,];
          p <- p/Matrix::rowSums(p);
          p <- p/sum(p);
          colnames(p) <- colnames(x)
          rownames(p) <- levels(batch)
        }
        rp <- list(p=p,depth=dep);
      } else {
        rp <- misc[['cpois']]
      }

      if(sparse) {

        rx <- x
        rx.gene <- rep(1:rx@Dim[2],diff(rx@p))
        if(is.null(batch)) {
          # -log prob of seeing more than 0 observations for all the non-zero observations, given site specific p and cell-specific depth
          rx@x <- -1*ppois(rx@x-1,rp$depth[rownames(rx)[rx@i+1]]*rp$p[colnames(rx)][rx.gene],lower.tail=FALSE,log.p=TRUE)
        } else {
          ps <- rp$p[,colnames(rx)]; gbi <- cbind(as.integer(batch[rownames(rx)[rx@i+1]]), rx.gene); # rows/column index of batch/gene elements in p
          rx@x <- -1*ppois(rx@x-1,rp$depth[rownames(rx)[rx@i+1]]*ps[gbi],lower.tail=FALSE,log.p=TRUE)

          # even out the residuals between batches for each gene
          tc <- colSumByFac(rx,as.integer(batch))[-1,]; # total residuals per batch (row) per gene (column)
          batch.sizes <- as.numeric(tapply(batch,batch,length));
          dataset.average <- Matrix::colSums(tc)/sum(batch.sizes); # dataset-wide averages
          batch.scaling <- t(t(tc/batch.sizes)/dataset.average);
          # find discrepant genes
          mins <- apply(batch.scaling,2,min); maxs <- apply(batch.scaling,2,max);
          batch.scaling[,mins<0.5 | maxs>2] <- 1e10; # take out the genes.

          rx@x <- rx@x/batch.scaling[gbi];
        }
      } else {
        if(!is.null(batch)) { stop("batch mode is not yet supported in the dense calculation")}
        rx <- do.call(cbind,mclapply(colnames(x), function(i) {
          nzi <- x[,i]>0;
          lp <- as.numeric(rp$depth*rp$p[i])
          lp[nzi] <- ppois(0,rp$depth[nzi]*rp$p[i],lower.tail=FALSE,log.p=TRUE)
          return(-1*lp)
        },mc.preschedule=T,mc.cores=n.cores))
      }

      # sdepth <- sort(depth,decreasing=T)
      # rx <- do.call(cbind,mclapply(colnames(x), function(i) {
      #   nzi <- x[,i]>0;
      #   plp <- lp <- as.numeric(rp$depth*rp$p[i])
      #   lp[nzi] <- ppois(0,rp$depth[nzi]*rp$p[i],lower.tail=FALSE,log.p=TRUE)
      #   names(plp) <- rownames(x);
      #   # scale relative to a perfect model where the occurrences are all within the top cells
      #   tn <- names(sdepth)[1:sum(nzi)];
      #   plp[tn] <- ppois(0,sdepth[1:sum(nzi)]*rp$p[i],lower.tail=FALSE,log.p=TRUE)
      #   #lp <- -1*lp; plp <- -1*plp;
      #   tl <- (sum(plp[plp<0]))/sum(lp<0)
      #   lp[lp<0] <- pmin(0,lp[lp<0]-tl);
      #   tl <- (sum(plp[plp>0]))/sum(lp>0)
      #   lp[lp>0] <- pmax(0,lp[lp>0]-tl);
      #   return(-1*lp)
      # },mc.preschedule=T,mc.cores=n.cores))

      rownames(rx) <- rownames(x)
      colnames(rx) <- colnames(x)

      reductions[[name]] <<- rx;

      return(invisible(rx))
    },

    # run PCA analysis on the overdispersed genes
    calculatePcaReduction=function(nPcs=20,type='counts', name='PCA', use.odgenes=FALSE, n.odgenes=NULL, odgenes=NULL, scale=F,center=T, cells=NULL) {

      if(type=='counts') {
        x <- counts;
      } else {
        if(!type %in% names(reductions)) { stop("reduction ",type,' not found')}
        x <- reductions[[type]]
      }
      if((use.odgenes || !is.null(n.odgenes)) && is.null(odgenes)) {
        if(is.null(misc[['devsites']] )) { stop("please run adjustVariance() first")}
        odgenes <- misc[['devsites']];
        if(!is.null(n.odgenes)) {
          if(n.odgenes>length(odgenes)) {
            #warning("number of specified odgenes is higher than the number of the statistically significant sites, will take top ",n.odgenes,' sites')
            odgenes <- rownames(misc[['varinfo']])[(order(misc[['varinfo']]$lp,decreasing=F)[1:min(ncol(counts),n.odgenes)])]
          } else {
            odgenes <- odgenes[1:n.odgenes]
          }
        }
      }
      if(!is.null(odgenes)) { x <- x[,odgenes] }


      require(irlba)
      if(!is.null(cells)) {
        # cell subset is just for PC determination
        pcs <- irlba(x[cells,], nv=nPcs, nu=0, center=Matrix::colMeans(x), right_only=FALSE)
      } else {
        pcs <- irlba(x, nv=nPcs, nu=0, center=Matrix::colMeans(x), right_only=FALSE)
      }
      rownames(pcs$v) <- colnames(x);


      misc$PCA <<- pcs;
      pcas <- as.matrix(x %*% pcs$v);
      #pcas <- scde::winsorize.matrix(pcas,0.05)
      # # control for sequencing depth
      # if(is.null(batch)) {
      #   mx <- model.matrix(x ~ d,data=data.frame(x=1,d=depth))
      # } else {
      #   mx <- model.matrix(x ~ d*b,data=data.frame(x=1,d=depth,b=batch))
      # }
      # # TODO: how to get rid of residual depth effects in the PCA-based clustering?
      # #pcas <- t(t(colLm(pcas,mx,returnResid=TRUE))+Matrix::colMeans(pcas))
      # pcas <- colLm(pcas,mx,returnResid=TRUE)
      rownames(pcas) <- rownames(x)
      colnames(pcas) <- paste('PC',seq(ncol(pcas)),sep='')
      #pcas <- pcas[,-1]
      #pcas <- scde::winsorize.matrix(pcas,0.1)

      reductions[[name]] <<- pcas;
      ## nIcs <- nPcs;
      ## a <- ica.R.def(t(pcas),nIcs,tol=1e-3,fun='logcosh',maxit=200,verbose=T,alpha=1,w.init=matrix(rnorm(nIcs*nPcs),nIcs,nPcs))
      ## reductions[['ICA']] <<- as.matrix( x %*% pcs$v %*% a);
      ## colnames(reductions[['ICA']]) <<- paste('IC',seq(ncol(reductions[['ICA']])),sep='');

      return(invisible(pcas))
    },

    getEmbedding=function(type='counts', embeddingType='largeVis', name=NULL, M=5, gamma=1, perplexity=100, sgd_batches=2e6, ... ) {
      if(type=='counts') {
        x <- counts;
      } else {
        if(!type %in% names(reductions)) { stop("reduction ",type,' not found')}
        x <- reductions[[type]]
      }
      if(is.null(name)) { name <- embeddingType }
      if(embeddingType=='largeVis') {
        xn <- misc[['edgeMat']][[type]];
        if(is.null(xn)) { stop(paste('KNN graph for type ',type,' not found. Please run makeKnnGraph with type=',type,sep='')) }
        #edgeMat <- sparseMatrix(i=xn$s+1,j=xn$e+1,x=xn$rd,dims=c(ncol(x),ncol(x)))
        edgeMat <- sparseMatrix(i=c(xn$s,xn$e)+1,j=c(xn$e,xn$s)+1,x=c(xn$rd,xn$rd),dims=c(nrow(x),nrow(x)))
        require(largeVis)
        #if(!is.null(seed)) { set.seed(seed) }
        wij <- buildWijMatrix(edgeMat,perplexity=perplexity)
        coords <- projectKNNs(wij = wij, M = M, verbose = TRUE,sgd_batches = sgd_batches,gamma=gamma, seed=1, ...)
        colnames(coords) <- rownames(x);
        emb <- embeddings[[type]][[name]] <<- t(coords);
      } else if(embeddingType=='tSNE') {
        require(Rtsne);
        cat("calculating distance ... ")
        #d <- dist(x);
        #d <- as.dist(1-WGCNA::cor(t(x), method = 'pearson', nThreads = n.cores))
        d <- as.dist(1-stats::cor(t(x)))
        cat("done\n")
        emb <- Rtsne(d,is_distance=T, perplexity=perplexity, ...)$Y;
        rownames(emb) <- labels(d)
        embeddings[[type]][[name]] <<- emb;
      } else if(embeddingType=='FR') {
        g <- graphs[[type]];
        if(is.null(g)){ stop(paste("generate KNN graph first (type=",type,")",sep=''))}
        emb <- layout.fruchterman.reingold(g, weights=E(g)$weight)
        rownames(emb) <- colnames(mat); colnames(emb) <- c("D1","D2")
        embeddings[[type]][[name]] <<- emb;
      } else {
        stop('unknown embeddingType ',embeddingType,' specified');
      }

      return(invisible(emb));
     }
  )
);

# a utility function to translate factor into colors
fac2col <- function(x,s=1,v=1,shuffle=FALSE,min.group.size=1,return.level.colors=F,unclassified.cell.color='gray50') {
  x <- as.factor(x);
  if(min.group.size>1) {
    x <- factor(x,exclude=levels(x)[unlist(tapply(rep(1,length(x)),x,length))<min.group.size])
  }
  col <- rainbow(length(levels(x)),s=s,v=v);
  if(shuffle) col <- sample(col);
  if(return.level.colors) { names(col) <- levels(x); return(col); }
  y <- col[as.integer(x)]; names(y) <- names(x);
  y[is.na(y)] <- unclassified.cell.color;
  y
}

val2col <- function(x,gradientPalette=NULL,zlim=NULL,gradient.range.quantile=0.95) {
  if(all(sign(na.omit(x))>=0)) {
    if(is.null(gradientPalette)) {
      gradientPalette <- colorRampPalette(c('gray90','red'), space = "Lab")(1024)
    }
    if(is.null(zlim)) {
      zlim <- as.numeric(quantile(x,p=c(1-gradient.range.quantile,gradient.range.quantile),na.rm=T))
      if(diff(zlim)==0) {
        zlim <- as.numeric(range(x))
      }
    }
    x[x<zlim[1]] <- zlim[1]; x[x>zlim[2]] <- zlim[2];
    x <- (x-zlim[1])/(zlim[2]-zlim[1])

  } else {
    if(is.null(gradientPalette)) {
      gradientPalette <- colorRampPalette(c("blue", "grey90", "red"), space = "Lab")(1024)
    }
    if(is.null(zlim)) {
      zlim <- c(-1,1)*as.numeric(quantile(abs(x),p=gradient.range.quantile,na.rm=T))
      if(diff(zlim)==0) {
        zlim <- c(-1,1)*as.numeric(max(abs(x)))
      }
    }
    x[x<zlim[1]] <- zlim[1]; x[x>zlim[2]] <- zlim[2];
    x <- (x-zlim[1])/(zlim[2]-zlim[1])

  }
  gradientPalette[x*(length(gradientPalette)-1)+1]
}


# note transpose is meant to speed up calculations when neither scaling nor centering is required
fast.pca <- function(m,nPcs=2,tol=1e-10,scale=F,center=F,transpose=F) {
  require(irlba)
  if(transpose) {
    if(center) { m <- m-Matrix::rowMeans(m)}; if(scale) { m <- m/sqrt(Matrix::rowSums(m*m)); }
    a <- irlba(tcrossprod(m)/(ncol(m)-1), nu=0, nv=nPcs,tol=tol);
    a$l <- t(t(a$v) %*% m)
  } else {
    if(scale||center) { m <- scale(m,scale=scale,center=center) }
    #a <- irlba((crossprod(m) - nrow(m) * tcrossprod(Matrix::colMeans(m)))/(nrow(m)-1), nu=0, nv=nPcs,tol=tol);
    a <- irlba(crossprod(m)/(nrow(m)-1), nu=0, nv=nPcs,tol=tol);
    a$l <- m %*% a$v
  }
  a
}

# quick utility to check if given character vector is colors
# thanks, Josh O'Brien: http://stackoverflow.com/questions/13289009/check-if-character-string-is-a-valid-color-representation
areColors <- function(x) {
  is.character(x) & sapply(x, function(X) {tryCatch(is.matrix(col2rgb(X)), error = function(e) FALSE)})
}

jw.disR <- function(x,y) { x <- x+1/length(x)/1e3; y <- y+1/length(y)/1e3; a <- x*log(x)  + y*log(y) - (x+y)*log((x+y)/2); sqrt(sum(a)/2)}

# translate multilevel segmentation into a dendrogram, with the lowest level of the dendrogram listing the cells
multi2dend <- function(cl,counts,deep=F) {
  if(deep) {
    clf <- as.integer(cl$memberships[1,]); # take the lowest level
  } else {
    clf <- as.integer(membership(cl));
  }
  names(clf) <- names(membership(cl))
  clf.size <- unlist(tapply(clf,factor(clf,levels=seq(1,max(clf))),length))
  rowFac <- rep(NA,nrow(counts));
  rowFac[match(names(clf),rownames(counts))] <- clf;
  lvec <- colSumByFac(counts,rowFac)[-1,];
  lvec.dist <- jsDist(t(lvec/pmax(1,Matrix::rowSums(lvec))));
  d <- as.dendrogram(hclust(as.dist(lvec.dist),method='ward.D'))
  # add cell info to the laves
  addinfo <- function(l,env) {
    v <- as.integer(mget("index",envir=env,ifnotfound=0)[[1]])+1;
    attr(l,'nodeId') <- v
    assign("index",v,envir=env)
    attr(l,'nCells') <- sum(clf.size[as.integer(unlist(l))]);
    if(is.leaf(l)) {
      attr(l,'cells') <- names(clf)[clf==attr(l,'label')];
    }
    attr(l,'root') <- FALSE;
    return(l);
  }
  d <- dendrapply(d,addinfo,env=environment())
  attr(d,'root') <- TRUE;
  d
}
# translate cell cluster dendrogram to an array, one row per node with 1/0 cluster membership
cldend2array <- function(d,cells=NULL) {
  if(is.null(cells)) { # figure out the total order of cells
    cells <- unlist(dendrapply(d,attr,'cells'))
  }
  getcellbin <- function(l) {
    if(is.leaf(l)) {
      vi <- match(attr(l,'cells'),cells)
      ra <- sparseMatrix(i=vi,p=c(0,length(vi)),x=rep(1,length(vi)),dims=c(length(cells),1),dimnames=list(NULL,attr(l,'nodeId')))
      return(ra);
    } else { # return rbind of the children arrays, plus your own
      ra <- do.call(cbind,lapply(l,getcellbin))
      ur <- unique(ra@i);
      ra <- cbind(sparseMatrix(ur+1,x=rep(1,length(ur)),p=c(0,length(ur)),dims=c(length(cells),1),dimnames=list(NULL,attr(l,'nodeId'))),ra);
      return(ra)
    }
  }
  a <- getcellbin(d)
  rownames(a) <- cells;
  return(t(a));
}


# peak annotation utils (based on Jean's code)

# hg38-specific utility function
# given a set of peaks ("chr:start-end", returns assignment to the ? gene using ChIPseeker)
hg38.peaks2Symbols <- function(peaks,includeDistal=FALSE, returnEntrezIds=FALSE,return.details=FALSE) {
  peak.df <- data.frame(do.call(rbind,strsplit(peaks,":|-")),stringsAsFactors=F); colnames(peak.df) <- c('chr','start','end');
  require(GenomicRanges)
  peaks.df <- with(peak.df, GRanges(chr, IRanges(as.numeric(start), as.numeric(end)), strand=NULL))

  ## Annotating peaks
  library(TxDb.Hsapiens.UCSC.hg38.knownGene)
  txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
  require(ChIPseeker)
  spp.ann <- ChIPseeker::annotatePeak(peaks.df, TxDb=txdb, verbose = T)
  spp.ann.df <- as.data.frame(spp.ann)
  #rownames(spp.ann.df) <- paste(peaks)

  if(!includeDistal) {
    spp.ann.df$geneId[spp.ann.df$annotation=='Distal Intergenic'] <- NA
  }

  if(returnEntrezIds) {
    sym <- spp.ann.df$geneId;
  } else {
    ## get symbols
    library(org.Hs.eg.db)
    spp.ann.df$symbol <- sym <- select(org.Hs.eg.db,keys=spp.ann.df$geneId,columns=c('SYMBOL'))$SYMBOL
  }
  if(return.details) { return(spp.ann.df)}
  names(sym) <- peaks;
  sym
}

# group nearby peaks (returns a factor given peak names)
group.nearby.peaks <- function(peaks,window.size=1e3,verbose=T) {
  peaks.df <- data.frame(do.call(rbind,strsplit(peaks,":|-")),stringsAsFactors=F); colnames(peaks.df) <- c('chr','start','end');
  peaks.df$pos <- (as.numeric(peaks.df$start)+as.numeric(peaks.df$end))/2
  peak.f <- rep(NA,nrow(peaks.df));

  for(chr in unique(peaks.df$chr)) {
    ii <- which(peaks.df$chr==chr); ii <- ii[order(peaks.df$pos[ii],decreasing=F)];
    cf <- nearbyPointsGreedyCluster(peaks.df$pos[ii],window.size);
    peak.f[ii] <- cf+max(c(0,peak.f),na.rm=TRUE);
  }
  peak.f <- as.factor(peak.f);
  # translate into coordinate names
  peak.cfn <- unlist(tapply(1:nrow(peaks.df),peak.f,function(ii) paste(peaks.df$chr[ii[1]],':',min(peaks.df$pos[ii]),'-',max(peaks.df$pos[ii]),sep='')))
  peak.f <- as.factor(peak.cfn[peak.f]);
  names(peak.f) <- peaks;
  if(verbose) { cat(round(100-length(levels(peak.f))/length(peaks)*100,2),"% reduction\n")}
  return(peak.f);
}

hg38.closestGenes2Peaks <- function(peaks,n=10) {
  # parse out peak names
  peak.df <- data.frame(do.call(rbind,strsplit(peaks,":|-")),stringsAsFactors=F); colnames(peak.df) <- c('chr','start','end');
  peak.df$mid <- (as.numeric(peak.df$end)+as.numeric(peak.df$start))/2

  # get gene annotation
  require(TxDb.Hsapiens.UCSC.hg38.knownGene)
  require(Homo.sapiens) # to look names homo sapiens can read
  # ... 20 packages and numberous overloaded functions later
  txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene;
  seqlevels(txdb) <- unique(peak.df$chr)
  genes <- genes(txdb)
  x <- tapply(1:nrow(peak.df),as.factor(peak.df$chr),function(ii) {
    ii <- ii[order(peak.df$mid[ii],decreasing=F)]
    chr <- peak.df$chr[ii[1]];
    gi <- seqnames(genes)==chr;
    gdf <- data.frame(start=start(genes[gi]),end=end(genes[gi]),entrezId=mcols(genes[gi])$gene_id,strand=as.integer(as.character(strand(genes[gi]))=="+"),stringsAsFactors=F)
    gdf <- gdf[order(gdf$start,decreasing=F),]
    x <- closestNSegmentsToPoints(gdf$start,gdf$end,peak.df$mid[ii],gdf$strand,n)
    x$i <- x$i+1; x$entrezId <- gdf$entrezId;
    rownames(x$i) <- rownames(x$d) <- rownames(x$s) <- peaks[ii];
    x
  })

  # combine joint index and distance matrices
  entrezIds <- unique(unlist(lapply(x,function(z) z$entrezId)))
  genes <- entrezIds;
  genes <- select(Homo.sapiens,keys=entrezIds,column=c("SYMBOL"),keytype="ENTREZID")$SYMBOL
  genes[is.na(genes)] <- entrezIds[is.na(genes)]

  # merge gene indices, translating them to the global index
  closestGenes <- do.call(rbind,lapply(x,function(z) {
    gids <- match(z$entrezId,entrezIds);
    gind <- apply(z$i,2,function(y) gids[y]);
    rownames(gind) <- rownames(z$i);
    return(gind)
  }))
  closestGenes <- closestGenes[peaks,]

  geneDistances <- do.call(rbind,lapply(x,function(z) {  z$d }))
  geneDistances <- geneDistances[peaks,]

  tssDistances <- do.call(rbind,lapply(x,function(z) {  z$s }))
  tssDistances <- tssDistances[peaks,]

  return(list(genes=genes,closestGenes=closestGenes,geneDistances=geneDistances,tssDistances=tssDistances))
}

hg38.closestPeaks2Genes <- function(peaks,n=10) {
  # parse out peak names
  peak.df <- data.frame(do.call(rbind,strsplit(peaks,":|-")),stringsAsFactors=F); colnames(peak.df) <- c('chr','start','end');
  peak.df$mid <- (as.numeric(peak.df$end)+as.numeric(peak.df$start))/2

  # get gene annotation
  require(TxDb.Hsapiens.UCSC.hg38.knownGene)
  require(Homo.sapiens) # to look names homo sapiens can read
  # ... 20 packages and numberous overloaded functions later
  txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene;
  seqlevels(txdb) <- unique(peak.df$chr)
  genes <- genes(txdb)
  x <- tapply(1:nrow(peak.df),as.factor(peak.df$chr),function(ii) {
    ii <- ii[order(peak.df$mid[ii],decreasing=F)]
    chr <- peak.df$chr[ii[1]];
    gi <- seqnames(genes)==chr;
    gdf <- data.frame(start=start(genes[gi]),end=end(genes[gi]),entrezId=mcols(genes[gi])$gene_id,strand=as.integer(as.character(strand(genes[gi]))=="+"),stringsAsFactors=F)
    gdf <- gdf[order(gdf$start,decreasing=F),]
    x <- closestNPointsToSegments(gdf$start,gdf$end,peak.df$mid[ii],gdf$strand,n)
    x$i <- x$i+1;
    x$peaks <- peaks[ii];
    rownames(x$s) <- rownames(x$i) <- rownames(x$d) <- gdf$entrezId;
    x
  })

  # combine joint index and distance matrices
  entrezIds <- unique(unlist(lapply(x,function(z) rownames(z$i))))
  genes <- entrezIds;
  genes <- select(Homo.sapiens,keys=entrezIds,column=c("SYMBOL"),keytype="ENTREZID")$SYMBOL
  genes[is.na(genes)] <- entrezIds[is.na(genes)]

  # merge gene indices, translating them to the global index
  closestPeaks <- do.call(rbind,lapply(x,function(z) {
    gids <- match(z$peaks,peaks);
    gind <- apply(z$i,2,function(y) gids[y]);
    rownames(gind) <- rownames(z$i);
    return(gind)
  }))
  closestPeaks <- closestPeaks[entrezIds,]
  peakDistances <- do.call(rbind,lapply(x,function(z) {  z$d }))
  peakDistances <- peakDistances[entrezIds,]
  peakTssDistances <- do.call(rbind,lapply(x,function(z) {  z$s }))
  peakTssDistances <- peakTssDistances[entrezIds,]
  rownames(closestPeaks) <- rownames(peakDistances) <- rownames(peakTssDistances) <- genes;
  return(list(peaks=peaks,closestPeaks=closestPeaks,peakDistances=peakDistances,peakTssDistances=peakTssDistances))
}


# returns GO -> gene symbol map
hg38.getGo2Symbols <- function(symbols,entrezIds=FALSE,min.go.size=5,max.go.size=Inf,return.list=FALSE) {
  library(org.Hs.eg.db)
  # translate gene names to ids
  ids <- unlist(lapply(mget(unique(na.omit(symbols)),org.Hs.egALIAS2EG,ifnotfound=NA),function(x) x[1]))
  # reverse map
  rids <- names(ids); names(rids) <- ids;
  # list all the ids per GO category
  go.env <- eapply(org.Hs.egGO2ALLEGS,function(x) as.character(na.omit(rids[x])))
  sz <- unlist(lapply(go.env,length));
  go.env <- go.env[sz>=min.go.size & sz<=max.go.size];
  if(return.list) { go.env } else { list2env(go.env)}
}

# invert one-to-many string map (list)
invert.string.list <- function(map) {
  df <- cbind(rep(names(map),unlist(lapply(map,length))),unlist(map))
  tapply(df[,1],as.factor(df[,2]),I)
}

hg38.getSymbols2Go <- function(symbols,entrezIds=FALSE, return.list=FALSE, ...) {
  forward <- hg38.getGo2Symbols(na.omit(unique(symbols)),entrezIds=entrezIds,return.list = TRUE, ...)
  # reverse map
  gene2go <- invert.string.list(forward)
  if(return.list) { gene2go } else {return(list2env(gene2go))}
}

hg38.symbols2Peaks <- function( peaks, ... ) {
  symbols <- hg38.peaks2Symbols( peaks, ... );
  # build a reverse map
  tapply(peaks,as.factor(symbols),I)
}

#' Given mapping to peaks to symbols and symbols to gene sets, maps peaks to gene sets
peaks2GO <- function(peaks2Symbols, go.env, min.size=5, max.size=Inf) {
  go.env.peaks <- lapply(go.env, function(g) {
    na.omit(unique(unlist(symbols2Peaks[g])))
  })
  size <- unlist(lapply(go.env.peaks, length))
  go.env.peaks <- go.env.peaks[size > min.size & size < max.size]
  return(go.env.peaks)
}

hg38.peaks2GO <- function(peaks, includeDistal=F, ... ) {
  peaks2Symbols <- hg38.symbols2Peaks(peaks, includeDistal=includeDistal);
  library(liger)
  go.env <- org.Hs.GO2Symbol.list
  library(GO.db)
  desc <- select(GO.db, keys = names(go.env), columns = c("TERM"), multiVals = "CharacterList")
  stopifnot(all(names(go.env) == desc$GOID))
  names(go.env) <- paste(names(go.env), desc$TERM)
  go.env.peaks <- peaks2GO(peaks2Symbols, go.env, ...)
}




ths.collapse.aspect.clusters <- function(d,ct, scale = TRUE, pick.top = FALSE) {
  xvm <- do.call(rbind, tapply(seq_len(nrow(d)), factor(ct, levels = sort(unique(ct))), function(ii) {
    if(length(ii) == 1) return(d[ii, ])
    if(pick.top) {
      return(d[ii[which.max(apply(d[ii, ], 1, var))], ])
    }
    xp <- pcaMethods::pca(t(d[ii, ]), nPcs = 1, center = TRUE, scale = "none")
    xv <- pcaMethods::scores(xp)[, 1]
    if(sum(abs(diff(xv))) > 0 && cor(xv, Matrix::colMeans(d[ii, ]*abs(pcaMethods::loadings(xp)[, 1])))<0) { xv <- -1*xv }
    #set scale at top pathway?
    if(sum(abs(diff(xv))) > 0) {
      if(scale) {
        xv <- xv*sqrt(max(apply(d[ii, ], 1, var)))/sqrt(var(xv))
      }
      if(sum(abs(xv)) == 0) { xv <- abs(rnorm(length(xv), sd = 1e-6)) }
    } else {
      xv <- abs(rnorm(length(xv), sd = 1e-6))
    }
    #xv <- xv/sqrt(length(ii))
    xv
  }))
  rownames(xvm) <- unlist(tapply(seq_len(nrow(d)), factor(ct, levels = sort(unique(ct))), function(ii) {
    if(length(ii) == 1) return(rownames(d)[ii])
    return(rownames(d)[ii[which.max(apply(d[ii, ], 1, var))]])
  }))
  return(xvm);
}


draw.color.key <- function(x,y,w,h,col=colorRampPalette(c("green","black","red"),space="Lab")(100),n=300,labels=c("low","high"),lab="expression",labcol="white",labcex=1.2,vert=F,standalone=F,minor.ticks=NULL,useRaster=F,cex=1) {
  require(gridBase)
  require(grid)
  if(!standalone) {
    #opar <- par(no.readonly=TRUE)
    #grid.newpage();
    vp <- baseViewports()
    pushViewport(vp$inner,vp$figure,vp$plot)
    pushViewport(viewport(x=x,y=y,width=w,height=h,just=c("left","top")))
    par(plt=gridPLT(),new=T)
  }

  if(vert) {
    las=2; axd=4; srt=90;
    image(z=matrix(seq(-1,1,length=n),nrow=1),col=col,yaxt="n",xaxt="n",useRaster=useRaster);
  } else {
    las=0; axd=1; srt=0;
    image(z=matrix(seq(-1,1,length=n),ncol=1),col=col,yaxt="n",xaxt="n",useRaster=useRaster);
  }
  par(cex=cex)
  axis(axd,at=seq(0,(length(labels)-1))/(length(labels)-1),labels=labels,las=las);
  if(!is.null(minor.ticks)) {
    axis(axd,at=(minor.ticks-labels[1])/diff(range(labels)),labels=FALSE,tcl=-0.2)
  }
  text(0.5,0.1,labels=lab,col=labcol,adj=c(0.5,0.5),cex=labcex,srt=srt)
  if(!standalone) {
    popViewport(4)
    #par(opar);

  }
}
