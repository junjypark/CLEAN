Clean=function(ymat, NNmatrix=NULL, group=NULL, nperm=10000, alpha=0.05, alternative=c("two.sided","less", "greater"), seed=NULL, 
               partition=T, npartition=NULL, parallel=F, ncores=1){
  if (is.null(NNmatrix)){
    cat("As NNmatrix is not specified, Clean conducts massive univariate analysis.\n")
    if (is.null(group)){
      fit=MassiveMean(ymat=ymat,
                      nperm=nperm, 
                      alpha=alpha, 
                      alternative=alternative,
                      seed=seed,
                      partition=partition, 
                      npartition = npartition,
                      parallel=parallel,
                      ncores=ncores)
      set.seed(NULL)
      return(fit)
    } else{
      fit=MassiveDiff(ymat=ymat,
                      group=group,
                      nperm=nperm, 
                      alpha=alpha, 
                      alternative=alternative,
                      seed=seed,
                      partition=partition, 
                      npartition = npartition,
                      parallel=parallel,
                      ncores=ncores)
      set.seed(NULL)
      return(fit) 
    }
  } else{
    if (is.null(group)){
      fit=CleanMean(ymat=ymat, 
                    NNmatrix=NNmatrix, 
                    nperm=nperm, 
                    alpha=alpha, 
                    alternative=alternative,
                    seed=seed,
                    partition=partition, 
                    npartition = npartition,
                    parallel=parallel,
                    ncores=ncores)
      set.seed(NULL)
      return(fit)
    }
    else{
      fit=CleanDiff(ymat=ymat, 
                    NNmatrix=NNmatrix, 
                    group=group, 
                    nperm=nperm, 
                    alpha=alpha, 
                    alternative=alternative,
                    seed=seed,
                    partition=partition, 
                    npartition = npartition,
                    parallel=parallel, 
                    ncores=ncores)
      set.seed(NULL)
      return(fit)
    }
  }
}


CleanMean=function(ymat, NNmatrix, nperm=10000, alpha=0.05, alternative=c("two.sided", "less", "greater"), seed=NULL, 
                   partition=T, npartition=1, parallel=F, ncores=1){
  if (length(which(is(NNmatrix)=="sparseMatrix"))==0){
    stop("NN is not a sparse matrix. Please refer the Matrix R package to convert it.")
  }
  if ( ncol(NNmatrix)!=nrow(ymat) ){
    stop("The number of columns of NN and the number of rows of ymat needs to be the same (# voxels).")
  }
  if ( alpha<0 |alpha>1){
    stop("alpha should range between 0 and 1.")
  }
  if (length(alternative)>1){ 
    alternative="two.sided"
    cat("Conducting the two-sided test as alternative has not been specified...\n")
  }
  if (is.null(seed)){ seed=sample(1e6,1) }
  if (isTRUE(partition)){
    if (is.null(npartition)){
      npartition=nrow(NNmatrix)%/%10000+1
    }
    else if (npartition==1){
      partition=FALSE
    }
  }
  
  if (isTRUE(partition)){
    NNList=list()
    len=ceiling(nrow(NNmatrix)/npartition)
    for (i in 1:npartition){
      start=len*(i-1)+1
      end=min(len*i,nrow(NNmatrix))
      NNList[[i]]=NNmatrix[start:end,]
    }
    
    if (isTRUE(parallel)){
      cl=makeCluster(ncores)
      registerDoParallel(cl)
      result=foreach(i=1:npartition, .packages=("Clean"),.noexport = "CleanC" )%dopar%{
        fit=CleanMeanC(ymat, NNList[[i]], nperm, seed)
        fit$alternative=alternative
        fit$seed=seed
        fit
      }
      stopCluster(cl)
    } else{
      result=list()
      for (i in 1:npartition){
        result[[i]]=CleanMeanC(ymat, NNList[[i]], nperm, seed)
        result[[i]]$alternative=alternative
        result[[i]]$seed=seed
      }
    }
    
    out=combine(result, alpha=alpha)
    out$nlocations=ncol(NNmatrix)
    return(out)
  } else{
    out=CleanMeanC(ymat, NNmatrix, nperm, seed)
    if (alternative=="less"){
      out$threshold=quantile(out$permMin,alpha)
      out$pvalue=(1+sum(c(out$permMin)<min(out$Tstat,na.rm=T)))/(1+nperm[1])
    } else if (alternative=="greater"){
      threshold=quantile(out$permMax,1-alpha)
      out$pvalue=(1+sum(c(out$permMax)>max(out$Tstat,na.rm=T)))/(1+nperm[1])
    } else {
      perm=pmax(abs(out$permMin),abs(out$permMax))
      out$threshold=quantile(pmax(abs(out$permMin),abs(out$permMax)),1-alpha)
      out$pvalue=(1+sum(c(perm)>max(abs(out$Tstat),na.rm=T)))/(1+nperm[1])
    }
    
    out$seed=seed
    out$alternative=alternative
    out$nlocations=ncol(NNmatrix)
    return(out)    
  }
}


CleanDiff=function(ymat, NNmatrix, group, nperm=10000, alpha=0.05, alternative=c("two.sided", "less", "greater"), seed=NULL, 
                   partition=T, npartition=1, parallel=F, ncores=1){
  if (length(group)!=ncol(ymat)){
    stop("The number of elements in group does not match with the number of columns in ymat.")
  }
  if (length(which(is(NNmatrix)=="sparseMatrix"))==0){
    stop("NN is not a sparse matrix. Please refer the Matrix R package to convert it.")
  }
  if ( ncol(NNmatrix)!=nrow(ymat) ){
    stop("The number of columns of NN and the number of rows of ymat needs to be the same (# voxels).")
  }
  if ( alpha<0 |alpha>1){
    stop("alpha should range between 0 and 1.")
  }
  if (length(alternative)>1){ 
    alternative="two.sided"
    cat("Conducting the two-sided test as alternative has not been specified...\n")
  }
  if (is.null(seed)){ seed=sample(1e6,1) }
  if (isTRUE(partition)){
    if (is.null(npartition)){
      npartition=nrow(NNmatrix)%/%10000+1
    }
    else if (npartition==1){
      partition=FALSE
    }
  }
  
  if (isTRUE(partition)){
    NNList=list()
    len=ceiling(nrow(NNmatrix)/npartition)
    for (i in 1:npartition){
      start=len*(i-1)+1
      end=min(len*i,nrow(NNmatrix))
      NNList[[i]]=NNmatrix[start:end,]
    }
    
    if (isTRUE(parallel)){
      cl=makeCluster(ncores)
      registerDoParallel(cl)
      result=foreach(i=1:npartition, .packages=("Clean"),.noexport = "CleanC" )%dopar%{
        fit=CleanDiffC(ymat, NNList[[i]], group, nperm, seed)
        fit$alternative=alternative
        fit$seed=seed
        fit
      }
      stopCluster(cl)
    } else{
      result=list()
      for (i in 1:npartition){
        result[[i]]=CleanDiffC(ymat, NNList[[i]], group, nperm, seed)
        result[[i]]$alternative=alternative
        result[[i]]$seed=seed
      }
    }
    
    out=combine(result, alpha=alpha)
    out$nlocations=ncol(NNmatrix)
    return(out)
  } else{
    out=CleanDiffC(ymat, NNmatrix, group, nperm, seed)
    if (alternative=="less"){
      out$threshold=quantile(out$permMin,alpha)
      out$pvalue=(1+sum(c(out$permMin)<min(out$Tstat,na.rm=T)))/(1+nperm[1])
    } else if (alternative=="greater"){
      threshold=quantile(out$permMax,1-alpha)
      out$pvalue=(1+sum(c(out$permMax)>max(out$Tstat,na.rm=T)))/(1+nperm[1])
    } else {
      perm=pmax(abs(out$permMin),abs(out$permMax))
      out$threshold=quantile(perm,1-alpha)
      out$pvalue=(1+sum(c(perm)>max(abs(out$Tstat),na.rm=T)))/(1+nperm[1])
    }
    
    out$seed=seed
    out$alternative=alternative
    out$nlocations=ncol(NNmatrix)
    return(out)    
  }
}


MassiveMean=function(ymat, nperm=10000, alpha=0.05, alternative=c("two.sided", "less", "greater"), seed=NULL, 
                     partition=T, npartition=NULL, parallel=F, ncores=1){
  if ( alpha<0 |alpha>1){
    stop("alpha should range between 0 and 1.")
  }
  if (length(alternative)>1){ 
    alternative="two.sided"
    cat("Conducting the two-sided test as alternative has not been specified...\n")
  }
  if (is.null(seed)){ seed=sample(1e6,1) }
  
  if (isTRUE(partition)){
    if (is.null(npartition)){
      npartition=nrow(ymat)%/%10000+1
    }
    
    ymatList=list()
    len=ceiling(nrow(ymat)/npartition)
    for (i in 1:npartition){
      start=len*(i-1)+1
      end=min(len*i,nrow(ymat))
      ymatList[[i]]=ymat[start:end,]
    }
    
    if (isTRUE(parallel)){
      cl=makeCluster(ncores)
      registerDoParallel(cl)
      result=foreach(i=1:npartition, .packages=("Clean"),.noexport = "CleanC" )%dopar%{
        fit=MassiveMeanC(ymatList[[i]], nperm, seed)
        fit$alternative=alternative
        fit$seed=seed
        fit
      }
      stopCluster(cl)
    } else{
      result=list()
      for (i in 1:npartition){
        result[[i]]=MassiveMeanC(ymatList[[i]], nperm, seed)
        result[[i]]$alternative=alternative
        result[[i]]$seed=seed
      }
    }
    
    out=combine(result, alpha=alpha)
    out$nlocations=nrow(ymat)
    return(out)
  } else{
    out=MassiveMeanC(ymat, nperm, seed)
    if (alternative=="less"){
      out$threshold=quantile(out$permMin,alpha)
      out$pvalue=(1+sum(c(out$permMin)<min(out$Tstat,na.rm=T)))/(1+nperm[1])
    } else if (alternative=="greater"){
      threshold=quantile(out$permMax,1-alpha)
      out$pvalue=(1+sum(c(out$permMax)>max(out$Tstat,na.rm=T)))/(1+nperm[1])
    } else {
      perm=pmax(abs(out$permMin),abs(out$permMax))
      out$threshold=quantile(perm,1-alpha)
      out$pvalue=(1+sum(c(perm)>max(abs(out$Tstat),na.rm=T)))/(1+nperm[1])
    }
    
    out$seed=seed
    out$alternative=alternative
    out$nlocations=nrow(ymat)
    return(out)    
  }
}


MassiveDiff=function(ymat, group, nperm=10000, alpha=0.05, alternative=c("two.sided", "less", "greater"), seed=NULL, 
                     partition=T, npartition=NULL, parallel=F, ncores=1){
  if (length(group)!=ncol(ymat)){
    stop("The number of elements in group does not match with the number of columns in ymat.")
  }
  if ( alpha<0 |alpha>1){
    stop("alpha should range between 0 and 1.")
  }
  if (length(alternative)>1){ 
    alternative="two.sided"
    cat("Conducting the two-sided test as alternative has not been specified...\n")
  }
  if (is.null(seed)){ seed=sample(1e6,1) }
  
  if (isTRUE(partition)){
    if (is.null(npartition)){
      npartition=nrow(ymat)%/%10000+1
    }
    
    ymatList=list()
    len=ceiling(nrow(ymat)/npartition)
    for (i in 1:npartition){
      start=len*(i-1)+1
      end=min(len*i,nrow(ymat))
      ymatList[[i]]=ymat[start:end,]
    }
    
    if (isTRUE(parallel)){
      cl=makeCluster(ncores)
      registerDoParallel(cl)
      result=foreach(i=1:npartition, .packages=("Clean"),.noexport = "CleanC" )%dopar%{
        fit=MassiveDiffC(ymatList[[i]], group, nperm, seed)
        fit$alternative=alternative
        fit$seed=seed
        fit
      }
      stopCluster(cl)
    } else{
      result=list()
      for (i in 1:npartition){
        result[[i]]=MassiveDiffC(ymatList[[i]], group, nperm, seed)
        result[[i]]$alternative=alternative
        result[[i]]$seed=seed
      }
    }
    
    out=combine(result, alpha=alpha)
    out$nlocations=nrow(ymat)
    
    return(out)
  } else{
    out=MassiveDiffC(ymat, group, nperm, seed)
    if (alternative=="less"){
      out$threshold=quantile(out$permMin,alpha)
      out$pvalue=(1+sum(c(out$permMin)<min(out$Tstat,na.rm=T)))/(1+nperm[1])
    } else if (alternative=="greater"){
      threshold=quantile(out$permMax,1-alpha)
      out$pvalue=(1+sum(c(out$permMax)>max(out$Tstat,na.rm=T)))/(1+nperm[1])
    } else {
      perm=pmax(abs(out$permMin),abs(out$permMax))
      out$threshold=quantile(pmax(abs(out$permMin),abs(out$permMax)),1-alpha)
      out$pvalue=(1+sum(c(perm)>max(abs(out$Tstat),na.rm=T)))/(1+nperm[1])
    }
    
    out$seed=seed
    out$alternative=alternative
    out$nlocations=nrow(ymat)
    return(out)    
  }
}

