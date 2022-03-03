buildNNmatrixDist_radius=function(distMat, max.radius=20){
  p=nrow(distMat)
  dist_range = sort(c(ceiling(min(distMat, na.rm=T)),ceiling(min(c(max(distMat,na.rm=T),max.radius),na.rm=T))))
  dist_range_sequence = seq(dist_range[1], dist_range[2], by=1) # different radii to consider
  
  out=foreach(i=1:p,.combine="rbind")%do%{
    dt=unname(distMat[i,])
    ind = which(dt<= dist_range[2]) 
    cbind(i,ind, dt[ind])
  }
  
  NNmatrix=foreach(r=1:length(dist_range_sequence), .combine="rbind")%do%{
    dist=dist_range_sequence[r]
    ind2=which(out[,3]<=dist)
    out2=out[ind2,]
    sp=sparseMatrix(i=out2[,1], j=out2[,2], x=1, dims=c(p,p))
    sp
  }
  
  return(NNmatrix)
}

