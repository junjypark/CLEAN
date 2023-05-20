process = function(fit, threshold = NULL){
  if (is.null(threshold)){ threshold = abs(fit$threshold[1]) }
  alternative = fit$alternative
  n.locations = fit$nlocations
  Tstat=rep(0, n.locations)

  if (alternative ==" two.sided"){
    Tstat = apply(matrix(fit$Tstat^2, n.locations), 1, max)
    Tstat_thresholded = Tstat
    Tstat_thresholded[Tstat < threshold] = 0
  } else if (alternative == "greater"){
    Tstat = apply(matrix(fit$Tstat, n.locations), 1, max)
    Tstat_thresholded = Tstat
    Tstat_thresholded[Tstat < threshold] = 0
  } else if (alternative == "less"){
    threshold = -threshold
    Tstat = apply(matrix(fit$Tstat,n.locations), 1, min)
    Tstat_thresholded = Tstat
    Tstat_thresholded[Tstat > threshold] = 0
  }
  
  return(list(Tstat = Tstat, 
              Tstat_thresholded = Tstat_thresholded))
}

