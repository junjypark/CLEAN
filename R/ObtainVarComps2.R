usethis::use_package("Matrix")



ObtainVarComps2=function(phi, epsilon, corMat_base1, corMat_base2){
  phi1=phi[1]
  phi2=phi[2]
  p=nrow(epsilon)
  n=ncol(epsilon)
  corMat1=corMat_base1^phi1
  corMat2=corMat_base2^phi2
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

