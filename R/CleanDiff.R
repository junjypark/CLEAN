CleanDiff=function(ymat, 
                   distmat, 
                   cortex,
                   mod0,
                   cov.interest, 
                   sacf,
                   max.radius,
                   nperm = 5000, 
                   alpha = 0.05, 
                   alternative = "two.sided", 
                   seed = NULL, 
                   nngp = T,
                   nngp.J = 50,
                   partition = T, 
                   npartition = 1, 
                   parallel = F, 
                   ncores = 1){
  
  V = nrow(ymat)
  if (!is.null(cortex)){
    ymat = ymat[cortex,]
    distmat = distmat[cortex,cortex]
  } else {
    cortex= 1:V
  }
  
  ymat.leverage = spLeverage(data=ymat, distmat=distmat, mod0=mod0, sacf=sacf, nngp=nngp, nngp.J=nngp.J)$out
  NNmatrix = buildNNmatrixDist(distmat, max.radius = max.radius)
  
  if (isTRUE(partition)){
    if (is.null(npartition)){
      npartition = nrow(NNmatrix)%/%10000+1
    }
    else if (npartition == 1){
      partition = FALSE
    }
  }
  
  if (isTRUE(partition)){
    NNList=list()
    len=ceiling(nrow(NNmatrix)/npartition)
    for (i in 1:npartition){
      start = len*(i-1)+1
      end = min(len*i,nrow(NNmatrix))
      NNList[[i]] = NNmatrix[start:end,]
    }
    
    if (isTRUE(parallel)){
      cl = makeCluster(ncores)
      registerDoParallel(cl)
      result = foreach(i = 1:npartition, .packages=("CLEAN"),.noexport = "CleanC" )%dopar%{
        fit = CleanDiffC(ymat.leverage, NNList[[i]], cov.interest, nperm, seed)
        fit$alternative = alternative
        fit$seed = seed
        fit
      }
      stopCluster(cl)
    } else{
      result = list()
      for (i in 1:npartition){
        result[[i]] = CleanDiffC(ymat.leverage, NNList[[i]], cov.interest, nperm, seed)
        result[[i]]$alternative = alternative
        result[[i]]$seed = seed
      }
    }
    
    out = combine(result, alpha = alpha, collapse = T)
    out$nlocations = ncol(NNmatrix)
    
    result_proc = process(out)
    out$Tstat = rep(0, V)
    out$Tstat[cortex]= result_proc$Tstat
    out$Tstat_thresholded = rep(0, V)
    out$Tstat_thresholded[cortex] = result_proc$Tstat_thresholded
    
    return(out)
  } else{
    out = CleanDiffC(ymat.leverage, NNmatrix, cov.interest, nperm, seed)
    if (alternative == "less"){
      out$threshold = quantile(out$permMin,alpha)
    } else if (alternative == "greater"){
      threshold = quantile(out$permMax,1-alpha)
    } else {
      perm = pmax(abs(out$permMin),abs(out$permMax))
      out$threshold = quantile(perm,1-alpha)
    }
    
    out$seed = seed
    out$alternative = alternative
    out$nlocations = ncol(NNmatrix)
    
    out = combine(result, alpha = alpha, collapse = T)
    result_proc = process(out)
    out$Tstat = rep(0, V)
    out$Tstat[cortex]= result_proc$Tstat
    out$Tstat_thresholded = rep(0, V)
    out$Tstat_thresholded[cortex] = result_proc$Tstat_thresholded
    
    return(out)    
  }
}

