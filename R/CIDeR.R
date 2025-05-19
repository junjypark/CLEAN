CIDeR=function(m1, m2, cov_df, distmat,
               cortex,
               cov.nuisance = NULL,
               cov.interest = NULL,
               mean_var_mod =list(method="between", mod_mean=NULL, var_formula=NULL, local_variance=T), #for mean and varaince modelling,
               #mean_var_mod =list(method="within", radius=15), # second approach
               #mean_var_mod =list(method="both", radius=15, mod_mean=NULL, var_formula=NULL), #combine two approaches
               spatial = T, #for spatial autocorrelation modelling
               sp.radius=0, #if scale the residual before spatial modelling
               sacf = "mix",
               nngp = T,
               nngp.J = 50,
               parallel=T,
               ncores=3,
               max.radius = 15, #for cluster enhancement
               nperm = 5000, #for permutation
               alpha = 0.05,
               seed=1) {
  V = ncol(m1)
  if (!is.null(cortex)){
    m1 = m1[, cortex]
    m2 = m2[, cortex]
    distmat = distmat[cortex, cortex]
  } else{
    cortex= 1:V
  }
  
  #StageI-step1, adjust the individual mean and variance
  print("stageI-step1")
  if (mean_var_mod$method=="between") {
    res <- MeanVarBetween(m1,m2,cov_df,
                          mean_var_mod$mod_mean,
                          mean_var_mod$var_formula,
                          parallel=parallel,
                          ncores=ncores)
    res1 <- res$res1
    res2 <- res$res2
    print("between")
    #StageI-step2, adjust the spatial autocorrelation
    if (mean_var_mod$local_variance==T) {
      print("pooling local varaince")
      LNNmatrix <- buildLNNmatrix(distmat, max.radius = max.radius)
      res1 <- adjust_data_local_sd(res1, LNNmatrix)
      res2 <- adjust_data_local_sd(res2, LNNmatrix)
    }
  } else if (mean_var_mod$method=="within") {
    LNNmatrix <- buildLNNmatrix(distmat, max.radius = mean_var_mod$radius)
    res1 <- MeanVarWithin(m1, LNNmatrix)
    res2 <- MeanVarWithin(m2, LNNmatrix)
    print("within")
  } else if (mean_var_mod$method=="both") {
    LNNmatrix <- buildLNNmatrix(distmat, max.radius = mean_var_mod$radius)
    res1 <- MeanVarboth(m1, LNNmatrix,
                        cov_df,
                        mean_var_mod$mod_mean,
                        mean_var_mod$var_formula,
                        parallel=parallel,
                        ncores=ncores)
    res2 <- MeanVarboth(m2, LNNmatrix,
                        cov_df,
                        mean_var_mod$mod_mean,
                        mean_var_mod$var_formula,
                        parallel=parallel,
                        ncores=ncores)
    print("both")
  }
    
  
 
  
  if (spatial) {
    print("stageI-step2_spatial")
    # if (scale.sp) {
    #   #scale the residual
    #   # Result: N × V → scaled version still N × V
    #   res1 <- sweep(res1, 2, apply(res1, 2, sd), `/`)
    #   res2 <- sweep(res2, 2, apply(res2, 2, sd), `/`)
    #   print("scale the data")
    # } 
    # we don't need this
   
    res1_sp_out <- safe_spLeverage(res1, distmat, sacf=sacf, nngp, nngp.J)
    res2_sp_out <- safe_spLeverage(res2, distmat, sacf=sacf, nngp, nngp.J)
    print("spatial autocorrelation finished")
    # res1_sp <- res1_sp_out$result$out
    # res2_sp <- res2_sp_out$result$out
    
    res1 <- t(res1_sp_out$result$out)
    res2 <- t(res2_sp_out$result$out)
  } else {
    print("skipping stageI-step2_spatial")
  }
 
  
  if (sp.radius > 0) {
      print("pooling local varaince after spatial")
      LNNmatrix <- buildLNNmatrix(distmat, max.radius = sp.radius)
      res1 <- adjust_data_local_sd(res1, LNNmatrix)
      res2 <- adjust_data_local_sd(res2, LNNmatrix)
  }
  
  print("stageII")
  
  
  rho <- res1*res2
  #StageIII regression under the null
  if (is.null(cov.nuisance)) {
    rho_res <- t(lm(rho~1)$residuals)
    print("stageIII")
  } else {
    rho_res <- t(lm(rho~as.matrix(cov_df[,cov.nuisance,drop = FALSE]))$residuals)
    print("stageIII")
  }
  
  
  
  #StageIV and V obtain test statistics and permutation 
  NNmatrix <- buildNNmatrixDist(distmat, max.radius = max.radius)
  out = Ciderperm(rho_res, as.matrix(cov_df[,cov.interest,drop = FALSE]), NNmatrix, nperm, seed) #V*N, N*p
  out$seed = seed
  out$nlocations = ncol(NNmatrix)
  out$alternative = "two.sided"
  out = combine(list(out), alpha = alpha, collapse = T)
  
  result_proc = process(out)
  
  out$Tstat = rep(0, V)
  out$Tstat[cortex]= result_proc$Tstat
  out$Tstat_thresholded = rep(0, V)
  out$Tstat_thresholded[cortex] = result_proc$Tstat_thresholded
  
  set.seed(NULL)
  print("stageIV")
  if (spatial) {
    return(list(out=out, rho=rho, sp=list(res1_sp=res1_sp_out,res2_sp=res2_sp_out)))  
  } else {
    return(list(out=out, rho=rho, sp=NULL))  
  }
}

