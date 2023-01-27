# Temporary version --- use with cautions
CleanV = function(ymat, NNmatrix, Kmatrix, Mmatrix=NULL,nperm=5000,seed=2) {
  
  if (!is.null(Kmatrix)) {
    fit = CleanVarC(ymat, NNmatrix, Matrix(Kmatrix,sparse=T), decomp=F, nperm, seed)
  } else {
    fit = CleanVarC(ymat, NNmatrix, Matrix(Mmatrix,sparse=T), decomp=T, nperm, seed)
  }
  
  return (list(Tstat =fit$Tstat, permMax= fit$permMax))
}
