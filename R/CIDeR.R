
CIDeR=function(m1, m2, cov_df, distmat,
               cortex=NULL,
               cov.nuisance = NULL,
               cov.interest = NULL,
               mean_var_mod =list(formula1=~1, sigma.formula1=~1, family1=NO(),
                                  formula2=~1, sigma.formula2=~1, family2=NO()), #for mean and varaince modelling,
               spatial = T, #for spatial autocorrelation modelling
               sp.radius=5, #if local pooling after spatial modelling
               sacf = "mix",
               nngp = F,
               nngp.J = 50,
               parallel=F,
               ncores=1,
               max.radius = 15, #for cluster enhancement
               nperm = 2000, #for permutation
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
  res1 <- MeanVarBetween(m1,cov_df,
                          mean_var_mod$formula1,
                          mean_var_mod$sigma.formula1,
                          mean_var_mod$family1)
  res2 <- MeanVarBetween(m2,cov_df,
                           mean_var_mod$formula2,
                           mean_var_mod$sigma.formula2,
                           mean_var_mod$family2)
  
 
  
  if (spatial) {
    print("stageI-step2_spatial")
   
    res1_sp_out <- safe_spLeverage(res1, distmat, sacf=sacf, nngp, nngp.J)
    res2_sp_out <- safe_spLeverage(res2, distmat, sacf=sacf, nngp, nngp.J)
    print("spatial autocorrelation finished")
    
    
    res1 <- t(res1_sp_out$result$out)
    res2 <- t(res2_sp_out$result$out)
  } else {
    print("skipping stageI-step2_spatial")
  }
 
  
  if (sp.radius > 0 ) {
      print("pooling local varaince")
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

