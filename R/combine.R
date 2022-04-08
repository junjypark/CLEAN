combine=function(lst, alpha=0.05){
  n=length(lst)
  seed=do.call("c",lapply(lst, function(x){x$seed}))
  nperm=do.call("c",lapply(lst, function(x){x$nperm}))
  nlocations=do.call("c",lapply(lst, function(x){x$nlocations}))
  alternative=do.call("c",lapply(lst, function(x){x$alternative}))
  
  if (length(unique(seed))>1){stop("[CLEAN] Use the same seed for every CLEAN output.")}
  else{seed=seed[1]}
  
  if (length(unique(alternative))>1){stop("[CLEAN] Use the same alternative for every CLEAN output.")}
  else{alternative=alternative[1]}
  
  if (length(unique(nperm))>1){stop("[CLEAN] Use the same number of permutations for every CLEAN output.")}
  else{nperm=nperm[1]}
  
  Tstat=do.call("c",lapply(lst, function(x){x$Tstat}))
  permMin=apply(do.call("cbind",lapply(lst, function(x){x$permMin})),1,min)
  permMax=apply(do.call("cbind",lapply(lst, function(x){x$permMax})),1,max)
  
  if (alternative=="less"){
    threshold=quantile(permMin,alpha)
    pvalue=(1+sum(c(permMin)<min(Tstat,na.rm=T)))/(1+nperm)
  } else if (alternative=="greater"){
    threshold=quantile(permMax,1-alpha)
    pvalue=(1+sum(c(permMax)>max(Tstat,na.rm=T)))/(1+nperm)
  } else {
    perm=pmax(abs(permMin),abs(permMax))
    threshold=quantile(perm,1-alpha)
    pvalue=(1+sum(c(perm)>max(abs(Tstat),na.rm=T)))/(1+nperm)
  }
  
  return(list(
    threshold=threshold,
    Tstat=Tstat,
    permMin=permMin,
    permMax=permMax,
    pvalue=pvalue,
    seed=seed,
    nperm=nperm,
    nlocations=nlocations,
    alternative=alternative
  ))
}


