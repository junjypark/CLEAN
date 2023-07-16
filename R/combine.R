combine=function(lst, alpha = 0.05, collapse = F){
  n=length(lst)
  seed=do.call("c",lapply(lst, function(x){x$seed}))
  nperm=do.call("c",lapply(lst, function(x){x$nperm}))
  nlocations=do.call("c",lapply(lst, function(x){x$nlocations}))
  alternative=do.call("c",lapply(lst, function(x){x$alternative}))
  
  if (length(unique(seed))>1){
    stop("[CLEAN] Use the same seed for every CLEAN output.")
  } else {
    seed=seed[1]
  }
  
  if (length(unique(alternative))>1){
    stop("[CLEAN] Use the same alternative for every CLEAN output.")
  } else { 
    alternative=alternative[1]
  }
  
  if (length(unique(nperm))>1){
    stop("[CLEAN] Use the same number of permutations for every CLEAN output.")
  } else { 
    nperm=nperm[1]
  }
  
  Tstat=do.call("c",lapply(lst, function(x){x$Tstat}))

  if (alternative=="two.sided"){
    permMin=apply(do.call("cbind",lapply(lst, function(x){x$permMin})),1,min)
    permMax=apply(do.call("cbind",lapply(lst, function(x){x$permMax})),1,max)
    perm=pmax(abs(permMin), abs(permMax))
    threshold=quantile(perm, 1-alpha)
  } else if (alternative=="greater"){
    permMin=NULL
    permMax=apply(do.call("cbind",lapply(lst, function(x){x$permMax})),1,max)
    threshold=quantile(permMax,1-alpha)
  } else if (alternative=="less"){
    permMax=NULL
    permMin=apply(do.call("cbind",lapply(lst, function(x){x$permMin})),1,min)
    threshold=quantile(permMin,alpha)
  } 
  
  if (collapse){
    lst=list(
      threshold = threshold,
      Tstat = Tstat,
      permMin = permMin,
      permMax = permMax,
      seed = seed,
      nperm = nperm,
      nlocations = nlocations,
      alternative = alternative
    )
  } else {
    for (j in 1:n){
      lst[[j]]$threshold = threshold
      lst[[j]]$nlocations = nlocations
    }    
  }
  
  return(lst)
}


