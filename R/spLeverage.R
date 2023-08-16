spLeverage=function(data, distmat=NULL, mod0=NULL, sacf="exp", nngp=T, nngp.J=50){
  if (!is.null(mod0)){
    data=t(lm(t(data)~mod0)$residuals)
    q=ncol(mod0)
  } else{
    q=0
  }
  
  covreg.fit=CovReg(data, distmat, n.covariates=q, kernel=sacf)
  
  if (nngp & nngp.J<nrow(distmat)){
    NNGPinfo=constructNNGPinfo(distmat, NN=min(nngp.J))
    NNGPprec=buildNNGPmat(distmat, NNGPinfo, covreg.fit, sacf)$NNGPprec
    data.new=as.matrix(NNGPprec%*%data)
  } else{
    
    phi=covreg.fit$phi
    sigma2=covreg.fit$sigma2
    tau2=covreg.fit$tau2
    
    if (sacf=="exp"){
      f.exp=exp(-phi*distmat)
      Sigma=sigma2*f.exp+tau2*diag(nrow(distmat))
    } else if (sacf=="gau"){
      f.gau=exp(-phi*distmat^2/2)
      Sigma=sigma2*f.gau+tau2*diag(nrow(distmat))
    } else if (sacf=="mix"){
      f.exp=exp(-phi*distmat)
      f.gau=exp(-phi*distmat^2/2)
      
      Sigma=sigma2[1]*f.exp+
        sigma2[2]*f.gau+
        tau2*diag(nrow(distmat))
    }
    
    data.new=solve(Sigma, data)
  }
  
  return(list(out = data.new, 
              params.covariance = covreg.fit))
}
