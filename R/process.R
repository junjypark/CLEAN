process=function(fit, threshold=NULL){
  if (is.null(threshold)){ threshold=fit$threshold }
  alternative=fit$alternative
  n.locations=fit$nlocations
  cl1=cl2=NULL
  if (alternative=="two.sided"){
    Tmax=apply(matrix(fit$Tstat,n.locations),1,max)
    Tmin=apply(matrix(fit$Tstat,n.locations),1,min)
    cl1=which(Tmax> threshold)
    cl2=which(Tmin< -threshold)
    inter=intersect(cl1,cl2)
    n.inter=length(inter)
    if (n.inter>0){
      for (j in 1:n.inter){
        if (abs(Tmax[inter[j]])>abs(Tmax[inter[j]])){
          cl2=setdiff(cl2, inter[j])
        } else{
          cl1=setdiff(cl1, inter[j])
        }
      }
    }
  } else if (alternative=="greater"){
    cl1=which(apply(matrix(fit$Tstat,n.locations),1,max)> threshold)
  } else if (alternative=="less"){
    cl2=which(apply(matrix(fit$Tstat,n.locations),1,min)< threshold)
  }
  
  return(list(indices.greater=cl1,indices.less=cl2))
}
