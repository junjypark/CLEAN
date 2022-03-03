buildNNmatrixDist=function(distMat, nnSet=c(1,5,10*1:10,50*3:10,100*6:10)){
  
  p=nrow(distMat)
  nnSet=unique(pmin(sort(nnSet),p))
  n.nnSet=length(nnSet)
  nnMax=nnSet[n.nnSet]
  
  out=foreach(i=1:p,.combine="rbind")%do%{
    dt=distMat[i,]
    rk=rank(dt)
    ind=which(rk<=nnMax)
    cbind(i,ind, rk[ind], dt[ind])
  }
  
  NNmatrix=foreach(r=1:n.nnSet, .combine="rbind")%do%{
    nn=nnSet[r]
    phi=phiSet[r]
    ind2=which(out[,3]<=nn)
    out2=out[ind2,]
    sp=sparseMatrix(i=out2[,1], j=out2[,2], x=1, dims=c(p,p))
    sp
  }
  
  return(NNmatrix)
}

