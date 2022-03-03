
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


