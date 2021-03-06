spLeverage=function(data, distMat=NULL, mod0=NULL, J=50){
  if (!is.null(mod0)){
    data=t(lm(t(data)~mod0)$residuals)
    q=ncol(mod0)
  } else{
    q=0
    cat("[CLEAN] mod0 has not been specified. Assuming the test for the mean.\n")
  }
  
  if (!is.null(distMat)){
    covreg.fit=CovReg(data, distMat, n.covariates=q)
    NNGPinfo=constructNNGPinfo(distMat, NN=J)
    NNGPprec=buildNNGPmat(distMat, NNGPinfo, covreg.fit)$NNGPprec
    data.new=as.matrix(NNGPprec%*%data)
  } else{
    cat("[CLEAN] distMat has not been specified. Obtaining residuals wihtout leveraging...\n")
  }
  
  return(list(out=data.new, params.covariance=covreg.fit,
              NNGPprec=NNGPprec, NNGPinfo=NNGPinfo))
}
