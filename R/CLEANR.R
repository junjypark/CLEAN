CleanR=function(xmat, ymat, NNmatrix, nperm=10000, alpha=0.05,
                alternative=c("two.sided", "less", "greater"),
                seed=NULL, npartition=50, parallel=F, ncores=1){
  if (length(which(is(NNmatrix)=="sparseMatrix"))==0){
    stop("NN is not a sparse matrix. Please refer the Matrix R package to convert it.")
  }
  if ( ncol(NNmatrix)!=nrow(ymat) ){
    stop("The number of columns of NN and the number of rows of ymat needs to be the same (# voxels).")
  }
  if ( alpha<0 | alpha>1){
    stop("alpha should range between 0 and 1.")
  }
  if (length(alternative)>1){
    alternative="two.sided"
    cat("Conducting the two-sided test as alternative has not been specified...\n")
  }
  if (is.null(seed)){ seed=sample(1e6,1) }
  
  # xmat=t(apply(xmat,1,scale))/ncol(xmat)
  # ymat=t(apply(ymat,1,scale))/ncol(ymat)
  NNList = list()
  len=ceiling(nrow(NNmatrix)/npartition)
  for (i in 1:npartition){
    start=len*(i-1)+1
    end=min(len*i,nrow(NNmatrix))
    NNList[[i]]=NNmatrix[start:end,]
  }
  
  CleanerPerm=CleanerPermC(xmat, ymat, nperm, seed)
  if (isTRUE(parallel)){
    cl=makeCluster(ncores)
    registerDoParallel(cl)
    result=foreach(i=1:npartition, .packages=("Clean"))%dopar%{
      fit=CleanerExpandPermC(CleanerPerm$U, CleanerPerm$permU, NNList[[i]])
      fit$alternative=alternative
      fit$seed=seed
      fit
    }
    stopCluster(cl)
  } else{
    result=list()
    for (i in 1:npartition){
      result[[i]]=CleanerExpandPermC(CleanerPerm$U, CleanerPerm$permU, NNList[[i]])
      result[[i]]$alternative=alternative
      result[[i]]$seed=seed
    }  
    out=combine(result, alpha=alpha)
    out$nlocations=ncol(NNmatrix)
  }
  
  
  out = list(combine_out = combine2(result, alpha=alpha),
             perm_null_dist = do.call("rbind",lapply(result, function(x){x$permNNU})))
  return(out)
}
