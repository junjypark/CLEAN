# buildNNmatrixDist=function(distMat, max.radius=20){
#   p=nrow(distMat)
#   dist_range = sort(c(ceiling(min(distMat, na.rm=T)),ceiling(min(c(max(distMat,na.rm=T),max.radius),na.rm=T))))
#   dist_range_sequence = seq(dist_range[1], dist_range[2], by=1) # different radii to consider
#   
#   out=foreach(i=1:p,.combine="rbind")%do%{
#     dt=unname(distMat[i,])
#     ind = which(dt<= dist_range[2]) 
#     cbind(i,ind, dt[ind])
#   }
#   
#   NNmatrix=foreach(r=1:length(dist_range_sequence), .combine="rbind")%do%{
#     dist=dist_range_sequence[r]
#     ind2=which(out[,3]<=dist)
#     out2=out[ind2,]
#     sp=sparseMatrix(i=out2[,1], j=out2[,2], x=1, dims=c(p,p))
#     sp
#   }
#   
#   return(NNmatrix)
# }



buildNNmatrixDist <- function(distMat, max.radius = 20) {
  p <- nrow(distMat)
  
  d_min <- ceiling(min(distMat, na.rm = TRUE))
  d_max <- ceiling(min(max(distMat, na.rm = TRUE), max.radius))
  radius_seq <- seq(d_min, d_max, by = 1)
  
  idx <- which(distMat <= d_max, arr.ind = TRUE)
  dvals <- distMat[idx]
  
  # Track all stacked sparse blocks
  row_offset <- 0
  NNmatrix_rows <- list()
  
  for (r in radius_seq) {
    keep <- which(dvals <= r)
    sub_i <- idx[keep, 1]
    sub_j <- idx[keep, 2]
    mat <- sparseMatrix(i = sub_i, j = sub_j, x = 1, dims = c(p, p))
    
    NNmatrix_rows[[length(NNmatrix_rows) + 1]] <- mat
  }
  
  # rbind all sparse matrices (this increases the number of rows!)
  NNmatrix <- do.call(rbind, NNmatrix_rows)
  return(NNmatrix)
}


buildLNNmatrix <- function(distMat, max.radius = 20) {
  p <- nrow(distMat)
  
  # Apply Gaussian kernel with cutoff
  weights <- matrix(1,p,p)
  weights[distMat > max.radius] <- 0
  
  # Get nonzero entries for sparse matrix construction
  idx <- which(weights > 0, arr.ind = TRUE)
  
  vals <- weights[idx]
  i <- idx[, 1]
  j <- idx[, 2]
  
  LNNmatrix <- sparseMatrix(i = i, j = j, x = 1, dims = c(p, p))
  
  return(LNNmatrix)
}



