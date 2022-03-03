usethis::use_package("Matrix")



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
      distvec=distMat[orders, vertex]
      elements=sort(distvec)[1:(NN-1)]
      ind.elements=(orders[which(distvec%in%elements)])[1:(NN-1)]
      out[j,]=c(ind.elements, vertex)
      # 
      # out[j,]=c(out[j-1, 2:min(j-1,NN)], vertex)
    }
    orders=c(orders, vertex)
    prev.vertex=vertex
  }
  
  return(list(NN=out, orders=orders))
}
