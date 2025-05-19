library(nlme)
library(parallel)
#m is N*V matrix
#covariate, N*p matrix, p covariates, 

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



# MeanVarBetween <- function(m1, m2,
#                           cov,                 # data‑frame of covariates
#                           mod  = NULL,   # formula for the mean model
#                           var_formula = NULL)  # formula for AdjustMeanVar
# {
#   stopifnot(is.matrix(m1), is.matrix(m2),
#             nrow(m1) == nrow(m2))
#   
#   V <- ncol(m1)                      # <‑‑ define once
#   
#   ## ------------------------------------------------------------------ ##
#   ## 1. helper that returns a V‑column matrix of residuals/adjustments  ##
#   ## ------------------------------------------------------------------ ##
#   adjust_block <- function(mat) {
#     
#     col_job <- function(v) {
#       if (!is.null(mod)) {                       # mean model supplied
#         if (!is.null(var_formula)) {             # heteroskedastic adj.
#           tryCatch(
#             AdjustMeanVar(mat[, v], cov, mod, var_formula)[["res_adjust"]],
#             error = function(e) {
#               residuals(lm(mod, data = cbind(cov, y = mat[, v])))
#             }
#           )
#         } else {                                 # plain linear residual
#           residuals(lm(mod, data = cbind(cov, y = mat[, v])))
#         }
#       } else {                                   # intercept‑only model
#         # residuals(lm(mat[, v] ~ 1))
#         mat[, v]
#       }
#     }
#     
#   
#       vapply(seq_len(V), col_job, numeric(nrow(mat)))  # serial
#   }
#   
#   ## run for both matrices ------------------------------------------- ##
#   list(res1 = adjust_block(m1),
#        res2 = adjust_block(m2))
# }
# 



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
              residuals(lm(mod, data = cbind(cov, y = mat[, v])))
            }
          )
        } else {
         
          residuals(lm(mod, data = cbind(cov, y = mat[, v])))
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



MeanVarboth <- function(m, NNmatrix, cov, mod = NULL, var_formula = NULL,
                        parallel = TRUE, ncores = 1) {
  stopifnot(is.matrix(m), nrow(NNmatrix) == ncol(m))

  N <- nrow(m)
  V <- ncol(m)

  adjust_block <- function(mat) {
    col_job <- function(v) {
      neigh_idx <- which(NNmatrix[v, ] != 0)
      if (length(neigh_idx) == 0) return(rep(NA_real_, N))  # no neighbors

      y_mat <- mat[, neigh_idx, drop = FALSE]  # N × N_r(v)
      y_vec <- as.vector(y_mat)                # N * N_r(v)

      cov_rep <- cov[rep(seq_len(N), length(neigh_idx)), , drop = FALSE]
      j_v <- which(neigh_idx == v)

      if (!is.null(mod)) {
        if (!is.null(var_formula)) {
          res_vec <- tryCatch(
            AdjustMeanVar(y_vec, cov_rep, mod, var_formula)[["res_adjust"]],
            error = function(e) {
              residuals(lm(mod, data = cbind(cov_rep, y = y_vec)))
            }
          )
        } else {
          res_vec <- residuals(lm(mod, data = cbind(cov_rep, y = y_vec)))
        }
        res_mat <- matrix(res_vec, nrow = N, byrow = F)
        return(res_mat[, j_v])
      } else {
        return(y_mat[, j_v])
      }
    }

    if (parallel) {
      parallel::mclapply(seq_len(V), col_job, mc.cores = ncores) |>
        simplify2array() # each output is length N, so transpose to get N × V
    } else {
      do.call(cbind, lapply(seq_len(V), col_job))
    }
  }

  res <- adjust_block(m)

  res <- adjust_data_local_sd(res, NNmatrix)
  return(res)
}

# MeanVarboth <- function(m, NNmatrix, cov, mod = NULL, var_formula = NULL,
#                         parallel = TRUE, ncores = 1) {
#   stopifnot(is.matrix(m), nrow(NNmatrix) == ncol(m))
#   
#   N <- nrow(m)
#   V <- ncol(m)
#   
#   adjust_block <- function(mat) {
#     col_job <- function(v) {
#       neigh_idx <- which(NNmatrix[v, ] != 0)
#       if (length(neigh_idx) == 0) return(rep(NA_real_, N))  # no neighbors
#       
#       y_mat <- mat[, neigh_idx, drop = FALSE]  # N × N_r(v)
#       y_vec <- as.vector(y_mat)                # N * N_r(v)
#       
#       cov_rep <- cov[rep(seq_len(N), length(neigh_idx)), , drop = FALSE]
#       j_v <- which(neigh_idx == v) 
#       
#       if (!is.null(mod)) {
#         if (!is.null(var_formula)) {
#           res_vec <- tryCatch(
#             AdjustMeanVar(y_vec, cov_rep, mod, var_formula)[["res_adjust"]],
#             error = function(e) {
#               residuals(lm(mod, data = cbind(cov_rep, y = y_vec)))
#             }
#           )
#         } else {
#           res_vec <- residuals(lm(mod, data = cbind(cov_rep, y = y_vec)))
#         }
#         res_mat <- matrix(res_vec, nrow = N, byrow = F)
#         return(res_mat[, j_v])
#       } else {
#         return(y_mat[, j_v])
#       }
#     }
#     
#     if (parallel) {
#       cl <- parallel::makeCluster(ncores)
#       on.exit(parallel::stopCluster(cl))
#       res_list <- parallel::parLapply(cl, seq_len(V), col_job)
#       do.call(cbind, res_list)
#     } else {
#       do.call(cbind, lapply(seq_len(V), col_job))
#     }
#   }
#   
#   res <- adjust_block(m)
#   res <- adjust_data_local_sd(res, NNmatrix)
#   return(res)
# }


MeanVarboth2 <- function(m, NNmatrix, cov, mod = NULL, var_formula = NULL,
                        parallel = TRUE, ncores = 1) {
  stopifnot(is.matrix(m), nrow(NNmatrix) == ncol(m))
  
  N <- nrow(m)
  V <- ncol(m)
  
  adjust_block <- function(mat) {
    process_vertex <- function(v) {
      neigh_idx <- which(NNmatrix[v, ] != 0)
      K <- length(neigh_idx)
      if (K == 0) return(rep(NA_real_, N))
      
      y_mat <- mat[, neigh_idx, drop = FALSE]
      y_vec <- matrix(y_mat, ncol = 1)
      
      cov_rep <- cov[rep(seq_len(N), each = K), , drop = FALSE]
      j_v <- which(neigh_idx == v)
      
      if (!is.null(mod)) {
        if (!is.null(var_formula)) {
          res_vec <- tryCatch(
            AdjustMeanVar(y_vec, cov_rep, mod, var_formula)[["res_adjust"]],
            error = function(e) {
              residuals(lm(mod, data = cbind(cov_rep, y = y_vec)))
            }
          )
        } else {
          print(v)
          res_vec <- residuals(lm(mod, data = cbind(cov_rep, y = y_vec)))
        }
        # reshape back to N × K matrix
        res_mat <- matrix(res_vec, nrow = N, ncol = K, byrow = TRUE)
        return(res_mat[, j_v])
      } else {
        return(y_mat[, j_v])
      }
    }
    
    res_list <- if (parallel) {
      parallel::mclapply(seq_len(V), process_vertex, mc.cores = ncores)
    } else {
      lapply(seq_len(V), process_vertex)
    }
    
    do.call(cbind, res_list)
  }
  
  res <- adjust_block(m)
  res <- adjust_data_local_sd(res, NNmatrix)
  
  return(res)
}



## Sample data
# set.seed(123)
# N <- 10   # subjects
# V <- 5    # vertices
# 
# m <- matrix(rnorm(N * V), nrow = N, ncol = V)
# 
# # Neighborhood matrix: assume each vertex is neighbor to itself and one other
# NNmatrix <- matrix(0, V, V)
# diag(NNmatrix) <- 1
# NNmatrix[cbind(1:(V-1), 2:V)] <- 1
# NNmatrix[cbind(2:V, 1:(V-1))] <- 1
# NNmatrix_sparse <- Matrix(NNmatrix, sparse = TRUE)
# # Covariates: data frame with intercept, age, sex
# cov <- data.frame(intercept = 1, age = rnorm(N), sex = sample(c(0,1), N, TRUE))
# 
# # Now test the function
# res <- MeanVarboth(m, NNmatrix_sparse , cov, mod = as.formula(y ~ age + sex), parallel = FALSE)
# print(res)


#N*V
CiderMeanVar=function(m1,m2,
                      cov, mod, var_formula=NULL,
                      parallel=T){
  if (!is.null(mod)) {
    if (!is.null(var_formula)) {
      V <- ncol(m1)
      if (parallel) {
        cores <- detectCores() - 10
        
        data1 <- do.call(cbind, mclapply(1:V, function(v) (AdjustMeanVar(m1[, v], cov, mod, var_formula))$res_adjust, mc.cores = cores))
        data2 <- do.call(cbind, mclapply(1:V, function(v) (AdjustMeanVar(m2[, v], cov, mod, var_formula))$res_adjust, mc.cores = cores))
      } else {
        data1 <- NULL
        data2 <- NULL
        for(v in 1:V) {
          data1 <- cbind(data1, AdjustMeanVar(m1[, v], cov, mod, var_formula)$res_adjust)
          data2 <- cbind(data2, AdjustMeanVar(m2[, v], cov, mod, var_formula)$res_adjust) 
        }
      }
      
      
    } else {
      data1 <- as.matrix(lm(m1~as.matrix(cov))$residuals, ncol=V)
      data2 <- as.matrix(lm(m2~as.matrix(cov))$residuals, ncol=V)
    }
  } else {
    data1 <-as.matrix(lm(m1~1)$residuals, ncol=V)
    data2 <-as.matrix(lm(m2~1)$residuals, ncol=V)
  }
  return(list(res1=data1, res2=data2))
}






CiderMeanVar2 <- function(m1, m2,
                         cov,                 # data‑frame of covariates
                         mod        = NULL,   # formula for the mean model
                         var_formula = NULL,  # formula for AdjustMeanVar
                         parallel   = TRUE,
                         ncores     = 1)
{
  stopifnot(is.matrix(m1), is.matrix(m2),
            nrow(m1) == nrow(m2))
  
  V <- ncol(m1)                      # <‑‑ define once
  
  ## ------------------------------------------------------------------ ##
  ## 1. helper that returns a V‑column matrix of residuals/adjustments  ##
  ## ------------------------------------------------------------------ ##
  adjust_block <- function(mat) {
    
    col_job <- function(v) {
      if (!is.null(mod)) {                       # mean model supplied
        if (!is.null(var_formula)) {             # heteroskedastic adj.
          AdjustMeanVar(mat[, v], cov, mod, var_formula)[["res_adjust"]]
        } else {                                 # plain linear residual
          residuals(lm(mod, data = cbind(cov, y = mat[, v])))
        }
      } else {                                   # intercept‑only model
        # residuals(lm(mat[, v] ~ 1))
        mat[, v]
      }
    }
    
    ## ----- choose parallel or serial execution ---------------------- ##
    if (parallel && .Platform$OS.type != "windows" && ncores > 1) {
      ncores <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", unset = parallel::detectCores()))
      print(ncores)
      do.call(cbind,
              parallel::mclapply(seq_len(V), col_job,
                                 mc.cores = ncores-8))
    } else {
      vapply(seq_len(V), col_job, numeric(nrow(mat)))  # serial
    }
  }
  
  ## run for both matrices ------------------------------------------- ##
  list(res1 = adjust_block(m1),
       res2 = adjust_block(m2))
}





# set.seed(42)
# 
# # Parameters
# N <- 100    # individuals
# V <- 2000   # vertices
# p1 <- 2    # nuisance covariates
# p2 <- 1    # covariates of interest
# 
# # Covariates
# cov.nuisance <- matrix(rnorm(N * p1), ncol = p1)
# cov.interest <- matrix(rnorm(N * p2), ncol = p2)
# covariate <- cbind(cov.nuisance, cov.interest)
# 
# # Design matrix
# X <- cbind(1, covariate)  # include intercept
# beta <- matrix(rnorm(ncol(X) * V, sd = 0.5), ncol = V)
# 
# # Mean structure
# mu1 <- X %*% beta
# mu2 <- mu1 + matrix(rnorm(N * V, sd = 0.3), nrow = N)  # slightly different signal
# 
# # Variance depends on covariate (e.g., covariate 1)
# sigma <-pmax(0.1, 0.5 + 0.5 * scale(covariate[, 1]))  # heteroskedastic
# 
# # Add noise
# m1 <- mu1 + matrix(rnorm(N * V, sd = rep(sigma, each = V)), nrow = N)
# m2 <- mu2 + matrix(rnorm(N * V, sd = rep(sigma, each = V)), nrow = N)
# 
# # Run adjustment
# microbenchmark(
#   # original = CiderMeanVar(m1, m2, cov.nuisance, cov.interest, mean = TRUE, var = TRUE),
#   fast = 
#   times = 3
# )
# res1 <- CiderMeanVar(m1, m2, cov.nuisance, cov.interest, mean = TRUE, var = TRUE, parallel=T)
# res2 <- CiderMeanVar(m1, m2, cov.nuisance, cov.interest, mean = TRUE, var = TRUE, parallel=F)
# # Output
# sum(res1$res1!=res2$res1)
# sum(res1$res2!=res2$res2)
# res2$res1
# res1$res2 



# -----------------------------------------------------
# Example usage of AdjustMeanVar with two covariates
# -----------------------------------------------------

# # 1) Generate some synthetic data
# set.seed(123)
# n <- 500
# x1 <- runif(n, 0, 1)               # first covariate
# x2 <- runif(n, 0, 1)
# m  <- rnorm(n, mean = 5, sd = sqrt(2*(1+(x1+x2)^2)))    # "response" variable
#               # second covariate
# 
# # Combine covariates into a data.frame
# cov<- cbind(x1,x2)
# 
# df <- data.frame(y = m, cov)
# # print(colnames(df))
# # print(colnames(cov_df))
# # var_formula <- as.formula(
# #   paste("~", paste(colnames(cov_df), collapse = "+"))
# # )
# # W <- varConstProp(form=var_formula)
# # mod <- as.formula(y~.)
# # 
# # AdjustMeanVar(m, cov_df, mean_formula, var_formula = W)
# # 
# # df <-  data.frame(y = m, cov_df)
# # model <- gls(mod,  data = df,
# #              weights = varConstPower(form=~x1+x2))
# # s <- sigma(model)  # scalar
# # 
# # # Get the weights (inverse relative standard deviations)
# # w <- varWeights(model$modelStruct$varStruct)  # vector of length n
# # 
# # # Compute estimated residual variances
# # res_var <- (s / w)^2
# # plot(2*(x1+x2)^2,res_var)
# # AdjustMeanVar(m, cov_df,  as.formula(y~.), varConstPower(form=~x1+x2))
# 
# CiderMeanVar(cbind(m,m), cbind(m,m), cbind(x1,x2), as.formula(y~.), varConstPower(form=~x1+x2))

