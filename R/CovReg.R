usethis::use_package("Matrix")

CovReg=function(epsilon,  distmat, kernel="exp", n.covariates=NULL, sparse=T, qtl=0.5, maxdist=NULL){
  if (is.null(n.covariates)){    
    n.covariates=0
  }
  if (!is.null(qtl) & !is.null(maxdist)){ stop("Only one of qtl or maxdist should be specified.")}
  if (sparse & is.null(qtl) & is.null (maxdist)){ stop("One of qtl or maxdist should be specified to enable the sparse option.") }
  
  q=NULL
  if (!is.null(qtl)){  q=quantile(distmat, qtl) }
  if (!is.null(maxdist)){ q=maxdist }
  
  n=ncol(epsilon); p=nrow(epsilon)
  if (kernel%in%c("exp", "gau")){
    if (kernel=="exp"){corMat.base=exp(-distmat)}
    else if (kernel=="gau"){corMat.base=exp(-distmat^2/2)}
    
    if (!is.null(q)){ corMat.base=ifelse(distmat<q, corMat.base, 0) }
    
    if (sparse){
      corMat.base=Matrix(corMat.base,sparse=T)
      phi.hat=optimize(CovRegOptim,interval=c(10^-5, 10),epsilon=epsilon, corMat_base=corMat.base)$`minimum`
      varcomps=ObtainVarComps(phi.hat, epsilon, corMat.base)
    } else{
      phi.hat=optimize(CovRegOptimC,interval=c(10^-5, 10),epsilon=epsilon, corMat_base=corMat.base)$`minimum`
      varcomps=ObtainVarCompsC(phi.hat, epsilon, corMat.base)
    }
  } else if (kernel=="mix"){
    corMat.base1=exp(-distmat)
    corMat.base2=exp(-distmat^2/2)
    
    if (!is.null(q)){
      corMat.base1=ifelse(distmat<q, corMat.base1, 0)
      corMat.base2=ifelse(distmat<q, corMat.base2, 0)
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
