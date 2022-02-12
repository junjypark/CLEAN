getSummaryMatrix=function(ymat, X=NULL, mask,
                          parallel=F, ncores=1){
  mask.index=which(mask!=0)
  ymat=ymat[mask.index,]
  p=nrow(ymat); n=ncol(ymat)
  
  if (!is.null(X)){
    n=nrow(X);
    Q=diag(n)-tcrossprod(tcrossprod(X,solve(crossprod(X))),X)
    ymat=tcrossprod(ymat,Q)
  }
  
  if (isTRUE(parallel)){
    cl=makeCluster(ncores)
    registerDoParallel(cl)
    out=foreach(i=1:p, .combine="rbind")%dopar%{
      ymat[i,]/sum(ymat[i,]^2,na.rm=T)
    }
    stopCluster(cl)
  } else{
    out=ymat
    for (i in 1:p){
      out[i,]=out[i,]/sum(out[i,]^2,na.rm=T)
    }
  }
  return(out)
}
