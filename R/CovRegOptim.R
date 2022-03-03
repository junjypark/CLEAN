usethis::use_package("Matrix")

CovRegOptim=function(phi, epsilon, corMat_base){
  p=nrow(epsilon)
  n=ncol(epsilon)
  corMat=corMat_base^phi
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
