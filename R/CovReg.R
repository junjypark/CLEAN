usethis::use_package("Matrix")

CovReg=function(epsilon,  distMat, kernel="exponential", n.covariates=NULL, sparse=T, qtl=0.5, maxdist=NULL){
  if (is.null(n.covariates)){    
    n.covariates=0
    cat("[CLEAN] Fitting the model assuming no covariate.\n")
    }
  if (!is.null(qtl) & !is.null(maxdist)){ stop("Only one of qtl or maxdist should be specified.")}
  if (sparse & is.null(qtl) & is.null (maxdist)){ stop("One of qtl or maxdist should be specified to enable the sparse option.") }
  
  q=NULL
  if (!is.null(qtl)){  q=quantile(distMat, qtl) }
  if (!is.null(maxdist)){ q=maxdist }
  
  n=ncol(epsilon); p=nrow(epsilon)
  if (kernel%in%c("exponential", "gaussian")){
    if (kernel=="exponential"){corMat.base=exp(-distMat)}
    else if (kernel=="gaussian"){corMat.base=exp(-distMat^2/2)}
    
    if (!is.null(q)){ corMat.base=ifelse(distMat<q, corMat.base, 0) }
    
    if (sparse){
      corMat.base=Matrix(corMat.base,sparse=T)
      phi.hat=optimize(CovRegOptim,interval=c(10^-5, 10),epsilon=epsilon, corMat_base=corMat.base)$`minimum`
      varcomps=ObtainVarComps(phi.hat, epsilon, corMat.base)
    } else{
      phi.hat=optimize(CovRegOptimC,interval=c(10^-5, 10),epsilon=epsilon, corMat_base=corMat.base)$`minimum`
      varcomps=ObtainVarCompsC(phi.hat, epsilon, corMat.base)
    }
  } else if (kernel=="mixture"){
    corMat.base1=exp(-distMat)
    corMat.base2=exp(-distMat^2/2)
    
    if (!is.null(q)){
      corMat.base1=ifelse(distMat<q, corMat.base1, 0)
      corMat.base2=ifelse(distMat<q, corMat.base2, 0)
    }
    
    corMat.base1=Matrix(corMat.base1,sparse=T)
    corMat.base2=Matrix(corMat.base2,sparse=T)
    phi.hat=optim(c(0.01, 0.01), CovRegOptim2,
                  epsilon=epsilon, 
                  corMat_base1=corMat.base1,
                  corMat_base2=corMat.base2)$par
    varcomps=ObtainVarComps2(phi.hat, epsilon, corMat.base1, corMat.base2)
  }
  return(list(sigma2=varcomps$sigma2*(n/(n-n.covariates)), 
              tau2=varcomps$tau2*(n/(n-n.covariates)), 
              phi=phi.hat,
              kernel=kernel))
}
