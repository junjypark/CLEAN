usethis::use_package("Matrix")


ObtainVarComps=function(phi, epsilon, corMat_base){
  p=nrow(epsilon)
  n=ncol(epsilon)
  corMat=corMat_base^phi
  corMat_norm=sum(corMat^2)
  
  y1=sum(epsilon*(corMat%*%epsilon))/n
  y2=sum(epsilon^2)/n
  sigma2= (y1-y2)/(corMat_norm-p)
  tau2= (-p*y1+corMat_norm*y2)/(corMat_norm*p-p^2)
  
  return(list(sigma2=sigma2, tau2=tau2))
}
