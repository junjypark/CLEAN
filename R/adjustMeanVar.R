AdjustMeanVar=function(m, cov, formula, sigma.formula, family) {
  df <-  data.frame(y = m, cov)
  fit <- gamlss(formula=formula, sigma.formula = sigma.formula, data = df,
                family = family)
  
  return(list(res_adjust= (fit$y - fitted(fit, what = "mu")) / fitted(fit, what = "sigma"),
              mu_hat = fitted(fit, what = "mu"),
              s=fitted(fit, what = "sigma"),
              coef=coef(fit,what="sigma")))
}


MeanVarBetween <- function(m,
                           cov,                 # data‑frame of covariates
                           formula = NULL,          # formula for the mean model
                           sigma.formula = NULL,
                           family=NO()#formula for the sd model
                           ) {
  stopifnot(is.matrix(m))
  
  V <- ncol(m)  # number of vertices
  N <- nrow(m)
  res_mat <- vapply(seq_len(V), function(v) {
    AdjustMeanVar(m[, v], cov, formula, sigma.formula, family)[["res_adjust"]]
  }, numeric(N))
  
  return(res_mat)
}


