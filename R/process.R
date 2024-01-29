process = function(fit, threshold = NULL){
  if (is.null(threshold)){ threshold = fit$threshold }
  alternative = fit$alternative
  n.locations = fit$nlocations

  if (alternative =="two.sided"){
    Tstat_max = apply(matrix(fit$Tstat, n.locations), 1, max)
    Tstat_min = apply(matrix(fit$Tstat, n.locations), 1, min)
    Tstat = ifelse(abs(Tstat_max) > abs(Tstat_min), Tstat_max, Tstat_min)
    Tstat_thresholded = Tstat
    Tstat_thresholded[abs(Tstat) < threshold] = 0    
  } else if (alternative == "greater"){
    Tstat = apply(matrix(fit$Tstat, n.locations), 1, max)
    Tstat_thresholded = Tstat
    Tstat_thresholded[Tstat < threshold] = 0
  } else if (alternative == "less"){
    Tstat = apply(matrix(fit$Tstat,n.locations), 1, min)
    Tstat_thresholded = Tstat
    Tstat_thresholded[Tstat > threshold] = 0
  }
  
  return(list(Tstat = Tstat, 
              Tstat_thresholded = Tstat_thresholded))
}

