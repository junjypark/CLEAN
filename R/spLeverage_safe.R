safe_spLeverage <- function(res, distmat, sacf = c("mix"), nngp, nngp.J) {
  
  try_spLeverage <- function(method) {
    tryCatch({
      out <- spLeverage_op(data = t(res), distmat = distmat, sacf = method, nngp = nngp, nngp.J = nngp.J)
      phi    <- out$params.covariance$phi
      sigma2 <- out$params.covariance$sigma2
      tau2   <- out$params.covariance$tau2
      
      if (any(c(phi, sigma2, tau2) < 0)) {
        message("sacf='", method, "' failed: negative phi/sigma2/tau2")
        NULL
      } else {
        list(result = out, method = method)
      }
    }, error = function(e) {
      message("sacf='", method, "' error: ", e$message)
      NULL
    })
  }
  
  # User explicitly requested "exp": skip "mix"
  if (sacf == "exp") {
    attempt_exp <- try_spLeverage("exp")
    if (!is.null(attempt_exp)) return(attempt_exp)
    
  } else if (sacf == "mix") {
    # Try "mix" first
    attempt_mix <- try_spLeverage("mix")
    if (!is.null(attempt_mix)) return(attempt_mix)
    
    # Then fallback to "exp"
    attempt_exp <- try_spLeverage("exp")
    if (!is.null(attempt_exp)) return(attempt_exp)
  }
  
  # Fallback to raw result
  message("All attempts failed. Returning raw transposed data.")
  list(result = list(out = t(res)), method = "none")
}
