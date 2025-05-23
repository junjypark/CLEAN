AdjustMeanVar=function(m, cov, mod, var_formula) {
  df <-  data.frame(y = m, cov)
  model <- gls(mod,  data = df,
                weights = var_formula)
  
  s <- sigma(model)  # scalar
  
  # Get the weights (inverse relative standard deviations)
  w <- varWeights(model$modelStruct$varStruct)  # vector of length n
  
  # Compute estimated residual variances
  res_var <- (s / w)^2
  
  res <- resid(model)
  return(list(res_adjust=as.vector(res/sqrt(res_var)),
              res = res,
              w = sqrt(res_var),
              model=model))
}


MeanVarBetween <- function(m1, m2,
                           cov,                 # data‑frame of covariates
                           mod = NULL,          # formula for the mean model
                           var_formula = NULL,  # formula for AdjustMeanVar
                           parallel = TRUE,
                           ncores = 1) {
  stopifnot(is.matrix(m1), is.matrix(m2),
            nrow(m1) == nrow(m2),
            ncol(m1) == ncol(m2))
  
  V <- ncol(m1)  # number of vertices
  N <- nrow(m1)
  
  ## ------------------------------------------------------------------ ##
  ## 1. helper that returns a V‑column matrix of residuals/adjustments ##
  ## ------------------------------------------------------------------ ##
  adjust_block <- function(mat) {
    col_job <- function(v) {
      
      if (!is.null(mod)) {
        if (!is.null(var_formula)) {
          tryCatch(
            AdjustMeanVar(mat[, v], cov, mod, var_formula)[["res_adjust"]],
            error = function(e) {
              fit <- lm(mod, data = cbind(cov, y = mat[, v]))
              residuals(fit) / summary(fit)$sigma
            }
          )
        } else {
         
          fit <- lm(mod, data = cbind(cov, y = mat[, v]))
          residuals(fit) / summary(fit)$sigma
        }
      } else {
        mat[, v]
      }
    }
    
    if (parallel) {
      result_list <- parallel::mclapply(seq_len(V), col_job, mc.cores = ncores)
      result_mat <- do.call(cbind, result_list)  # N × V
    } else {
      result_mat <- vapply(seq_len(V), col_job, numeric(N))
    }
    
    return(result_mat)
  }
  
  ## ------------------------------------------------------------------ ##
  ## 2. run for both matrices and return result                        ##
  ## ------------------------------------------------------------------ ##
  list(
    res1 = adjust_block(m1),
    res2 = adjust_block(m2)
  )
}



