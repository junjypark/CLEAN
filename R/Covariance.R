usethis::use_package("Matrix")


CovReg=function(epsilon,  distMat, kernel="exponential", sparse=T, qtl=0.5, maxdist=NULL){
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
      rho.hat=optimize(CovRegOptim,interval=c(10^-5, 10),epsilon=epsilon, corMat_base=corMat.base)$`minimum`
      varcomps=ObtainVarComps(rho.hat, epsilon, corMat.base)
    } else{
      rho.hat=optimize(CovRegOptimC,interval=c(10^-5, 10),epsilon=epsilon, corMat_base=corMat.base)$`minimum`
      varcomps=ObtainVarCompsC(rho.hat, epsilon, corMat.base)
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
    rho.hat=optim(c(0.01, 0.01), CovRegOptim2,
                  epsilon=epsilon, 
                  corMat_base1=corMat.base1,
                  corMat_base2=corMat.base2)$par
    varcomps=ObtainVarComps2(rho.hat, epsilon, corMat.base1, corMat.base2)
  }
  return(list(sigma2=varcomps$sigma2, 
              tau2=varcomps$tau2, 
              rho=rho.hat,
              kernel=kernel))
}

CovRegOptim=function(rho, epsilon, corMat_base){
  p=nrow(epsilon)
  n=ncol(epsilon)
  corMat=corMat_base^rho
  corMat_norm=sum(corMat^2)
  if (corMat_norm>p+1e-10){
    y1=sum(epsilon*(corMat%*%epsilon))/n
    y2=sum(epsilon^2)/n
    sigma2= (y1-y2)/(corMat_norm-p)
    tau2= (-p*y1+corMat_norm*y2)/(corMat_norm*p-p^2)
    ss= -2*(sigma2*y1*n+tau2*y2*n)+p*tau2*tau2+corMat_norm*sigma2*sigma2+2*sigma2*tau2*p;
  } else{
    sigma2=-1;
    tau2=-1;
    ss=10^10;
  }
  return(ss)
}

CovRegOptim2=function(rho, epsilon, corMat_base1, corMat_base2){
  rho1=rho[1]
  rho2=rho[2]
  p=nrow(epsilon)
  n=ncol(epsilon)
  corMat1=corMat_base1^rho1
  corMat2=corMat_base2^rho2
  corMat1_norm=sum(corMat1^2)
  corMat2_norm=sum(corMat2^2)
  corMat12_norm=sum(corMat1*corMat2)
  Mat=matrix(c(corMat1_norm, corMat12_norm, p, 
               corMat12_norm, corMat2_norm, p,
               p,p,p),3,3)
  if (det(Mat)>0){
    y1=sum(epsilon*(corMat1%*%epsilon))/n
    y2=sum(epsilon*(corMat2%*%epsilon))/n
    y3=sum(epsilon^2)/n
    y=c(y1,y2,y3)
    params=solve(Mat, y)
    sigma1=params[1]
    sigma2=params[2]
    tau2=params[3]
    ss= sigma1^2*corMat1_norm+sigma2^2*corMat2_norm+tau2^2*p+
      2*(sigma1*sigma2*corMat12_norm+sigma1*tau2*p+sigma2*tau2*p)-
      2*(sigma1*y1*n+sigma2*y2*n+tau2*y3*n)
  } else{
    sigma1=-1;
    sigma2=-1;
    tau2=-1;
    ss=10^10;
  }
  return(ss)
}

ObtainVarComps=function(rho, epsilon, corMat_base){
  p=nrow(epsilon)
  n=ncol(epsilon)
  corMat=corMat_base^rho
  corMat_norm=sum(corMat^2)

  y1=sum(epsilon*(corMat%*%epsilon))/n
  y2=sum(epsilon^2)/n
  sigma2= (y1-y2)/(corMat_norm-p)
  tau2= (-p*y1+corMat_norm*y2)/(corMat_norm*p-p^2)

  return(list(sigma2=sigma2, tau2=tau2))
}

ObtainVarComps2=function(rho, epsilon, corMat_base1, corMat_base2){
  rho1=rho[1]
  rho2=rho[2]
  p=nrow(epsilon)
  n=ncol(epsilon)
  corMat1=corMat_base1^rho1
  corMat2=corMat_base2^rho2
  corMat1_norm=sum(corMat1^2)
  corMat2_norm=sum(corMat2^2)
  corMat12_norm=sum(corMat1*corMat2)
  Mat=matrix(c(corMat1_norm, corMat12_norm, p, 
               corMat12_norm, corMat2_norm, p,
               p,p,p),3,3)

  y1=sum(epsilon*(corMat1%*%epsilon))/n
  y2=sum(epsilon*(corMat2%*%epsilon))/n
  y3=sum(epsilon^2)/n
  y=c(y1,y2,y3)
  params=solve(Mat, y)
  sigma2=params[1:2]
  tau2=params[3]

  return(list(sigma2=sigma2, tau2=tau2))
}


constructNNGPinfo=function(distMat, NN=50, start.vertex=NULL){
  orders=NULL
  m=nrow(distMat)
  out=matrix(NA, m,NN)
  if (is.null(start.vertex)){
    sum.dist=apply(distMat,1,sum)
    prev.vertex=which.min(sum.dist)
  }
  else{prev.vertex=start.vertex}
  
  current.indices=1:m
  out[1,1]=prev.vertex
  orders=c(orders, prev.vertex)
  
  for (j in 2:m){
    current.indices=setdiff(current.indices, prev.vertex)
    vertex=current.indices[which.min(distMat[prev.vertex, current.indices])]
    if (j<=NN){
      out[j,1:j]=c(na.omit(out[j-1,]),vertex) 
    } else{
      out[j,1:min(j,NN)]=c(out[j-1, 2:min(j-1,NN)], vertex)
    }
    orders=c(orders, vertex)
    prev.vertex=vertex
  }
  
  return(list(NN=out, orders=orders))
}

buildNNGPmat=function(distMat, NNGPinfo, params, kernel = "exponential"){
  m=nrow(distMat)
  A=matrix(0,m,m)
  D=matrix(0,m,m)
  rho=params$rho
  sigma2=params$sigma2
  tau2=params$tau2
  k=params$K
  
  f.exp=function(rho, d){ exp(-rho*d) }
  f.gau=function(rho, d){ exp(-rho*d^2/2)}

  for (i in 1:(nrow(NNGPinfo$NN)-1)){
    nn=na.omit(NNGPinfo$NN[i+1,])
    lnn=length(nn)
    coordip1=nn[lnn]
    
    if (kernel=="exponential"){
      K=sigma2*f.exp(rho,distMat[nn,nn,drop=F])+tau2*diag(lnn)
    } else if (kernel=="gaussian"){
      K=sigma2*f.gau(rho,distMat[nn,nn,drop=F])+tau2*diag(lnn)
    } else if (kernel=="mixture"){
      K=sigma2[1]*f.exp(rho[1],distMat[nn,nn,drop=F])+
        sigma2[2]*f.gau(rho[2],distMat[nn,nn,drop=F])+
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
