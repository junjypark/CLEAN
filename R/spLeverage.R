spLeverage=function(ymat, distMat=NULL, mod0=NULL, J=50){
  if (!is.null(mod0)){
    data=t(lm(t(data)~mod0)$residuals)
  } else{
    cat("mod0 has not been specified. Assuming the test for the mean.\n")
  }
  
  if (!is.null(distMat)){
    covreg.fit=CovReg(data, distMat)
    NNGPinfo=constructNNGPinfo(distMat, NN=50)
    NNGPprec=buildNNGPmat(distMat, NNGPinfo, covreg.fit)$NNGPprec
    data.new=as.matrix(NNGPprec%*%data)
  } else{
    cat("distMat has not been specified. Obtaining residuals wihtout leveraging...\n")
  }
  
  return(list(out=data.new, params.covariance=covreg.fit,
              NNGPprec=NNGPprec, NNGPinfo=NNGPinfo))
}

# getSummaryMatrix=function(ymat, X=NULL, mask,
#                           parallel=F, ncores=1){
#   mask.index=which(mask!=0)
#   ymat=ymat[mask.index,]
#   p=nrow(ymat); n=ncol(ymat)
#   
#   if (!is.null(X)){
#     n=nrow(X);
#     Q=diag(n)-tcrossprod(tcrossprod(X,solve(crossprod(X))),X)
#     ymat=tcrossprod(ymat,Q)
#   }
#   
#   if (isTRUE(parallel)){
#     cl=makeCluster(ncores)
#     registerDoParallel(cl)
#     out=foreach(i=1:p, .combine="rbind")%dopar%{
#       ymat[i,]/sum(ymat[i,]^2,na.rm=T)
#     }
#     stopCluster(cl)
#   } else{
#     out=ymat
#     for (i in 1:p){
#       out[i,]=out[i,]/sum(out[i,]^2,na.rm=T)
#     }
#   }
#   return(out)
# }
