usethis::use_package("Matrix")



buildNNGPmat=function(distMat, NNGPinfo, params, kernel = "exponential"){
  m=nrow(distMat)
  A=matrix(0,m,m)
  D=matrix(0,m,m)
  phi=params$phi
  sigma2=params$sigma2
  tau2=params$tau2
  k=params$K
  
  f.exp=function(phi, d){ exp(-phi*d) }
  f.gau=function(phi, d){ exp(-phi*d^2/2)}
  
  for (i in 1:(nrow(NNGPinfo$NN)-1)){
    nn=na.omit(NNGPinfo$NN[i+1,])
    lnn=length(nn)
    coordip1=nn[lnn]
    
    if (kernel=="exponential"){
      K=sigma2*f.exp(phi,distMat[nn,nn,drop=F])+tau2*diag(lnn)
    } else if (kernel=="gaussian"){
      K=sigma2*f.gau(phi,distMat[nn,nn,drop=F])+tau2*diag(lnn)
    } else if (kernel=="mixture"){
      K=sigma2[1]*f.exp(phi[1],distMat[nn,nn,drop=F])+
        sigma2[2]*f.gau(phi[2],distMat[nn,nn,drop=F])+
        tau2*diag(lnn)
    }
    
    A[coordip1,nn[-lnn]]=solve(K[-lnn,-lnn], K[lnn,-lnn])
    D[coordip1,coordip1]=K[lnn,lnn]-sum(K[lnn,-lnn]*solve(K[-lnn,-lnn], K[lnn,-lnn]))
  }
  D[NNGPinfo$NN[1,1],NNGPinfo$NN[1,1]]=sum(sigma2)+tau2
  
  IA=Matrix(diag(ncol(A))-A,sparse=T)
  sD=Matrix(diag(1/diag(D)),sparse=T)
  NNGPprec=IA%*%sD%*%Matrix::t(IA)
  return(list(A=IA, D=sD, NNGPprec=NNGPprec))
}
