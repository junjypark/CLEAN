CleanMean=function(ymat, 
                   distmat, 
                   cortex,
                   sacf,
                   nperm = 5000, 
                   alpha = 0.05, 
                   alternative = c("two.sided", "less", "greater"), 
                   seed = NULL,
                   nngp = T,
                   nngp.J = 50,
                   partition = T, 
                   npartition = 1, 
                   parallel = F, 
                   ncores = 1){
  
  if (!is.null(cortex)){
    ymat = ymat[cortex,]
    distmat = distmat[cortex,cortex]
  }

  ymat.leverage = spLeverage(ymat, distmat, sacf, nngp, nngp.J)
  
  
  if (isTRUE(partition)){
    if (is.null(npartition)){
      npartition = nrow(NNmatrix)%/%10000+1
    }
    else if (npartition == 1){
      partition = FALSE
    }
  }
  
  if (isTRUE(partition)){
    NNList = list()
    len = ceiling(nrow(NNmatrix)/npartition)
    for (i in 1:npartition){
      start = len*(i-1)+1
      end = min(len*i,nrow(NNmatrix))
      NNList[[i]] = NNmatrix[start:end,]
    }
    
    if (isTRUE(parallel)){
      cl = makeCluster(ncores)
      registerDoParallel(cl)
      result = foreach(i = 1:npartition, .packages=("CLEAN"),.noexport = "CleanC" )%dopar%{
        fit = CleanMeanC(ymat.leverage, NNList[[i]], nperm, seed)
        fit$alternative = alternative
        fit$seed = seed
        fit
      }
      stopCluster(cl)
    } else{
      result = list()
      for (i in 1:npartition){
        result[[i]] = CleanMeanC(ymat.leverage, NNList[[i]], nperm, seed)
        result[[i]]$alternative = alternative
        result[[i]]$seed = seed
      }
    }
    
    out = combine(result, alpha=alpha)
    out$nlocations = ncol(NNmatrix)
    return(out)
  } else{
    out=CleanMeanC(ymat.leverage, NNmatrix, nperm, seed)
    if (alternative == "less"){
      out$threshold = quantile(out$permMin,alpha)
      out$pvalue = (1+sum(c(out$permMin) < min(out$Tstat, na.rm = T)))/(1+nperm[1])
    } else if (alternative == "greater"){
      threshold = quantile(out$permMax, 1-alpha)
      out$pvalue = (1+sum(c(out$permMax)>max(out$Tstat,na.rm=T)))/(1+nperm[1])
    } else {
      perm = pmax(abs(out$permMin),abs(out$permMax))
      out$threshold = quantile(pmax(abs(out$permMin),abs(out$permMax)),1-alpha)
      out$pvalue = (1+sum(c(perm)>max(abs(out$Tstat),na.rm=T)))/(1+nperm[1])
    }
    
    out$seed = seed
    out$alternative = alternative
    out$nlocations = ncol(NNmatrix)
    return(out)    
  }
}

