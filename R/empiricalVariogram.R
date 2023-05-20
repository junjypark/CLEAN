get.empirical.variogram=function(epsilon, distmat, qtl=0.5, n.bins=500){
  qt.dist=quantile(distmat, qtl)
  up=upper.tri(distmat,diag=F) * (distmat<qt.dist)
  distvec=distmat[up==1]
  distIndex=which(up==1, arr.ind=T) # pairs of locations in upper triangle of distance matrix
  
  bins=seq(min(distvec),qt.dist,length.out=2*n.bins+1)
  points=bins[2*(1:n.bins)]
  start=bins[2*(1:n.bins)-1] # lower distance limit
  end=bins[2*(1:n.bins)+1] # upper distance limit
  result=vector(mode="numeric",length=n.bins)
  for (bin in 1:n.bins){
    index=which(distvec>=start[bin] &  distvec<end[bin]) # indices in distvec with distances between "start" and "end" distances for current bin
    index.sub=distIndex[index,,drop=FALSE] # row/column indices for those positions
    result[bin]=mean(apply((epsilon[,index.sub[,1]]-epsilon[,index.sub[,2]])^2,2,mean))
  }
  
  return(list(empirical=result, points=start))
}
