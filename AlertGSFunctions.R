## Functions needed in the context of applying the AlertGS methodology to data


# First: likelihood functions for rectified Gauss and Gumbel distributions

# Gauss
# Input:
#   par: vector of the two parameters mu and sigma of a normal distribution
#   x:   vector of data points used for fitting the likelihood
# Output:
#   value of the log likelihood of the rectified Gauss distribution for the observed data points and given parameters
llrectifiedgauss <- function(par, x){
  mu <- par[1]
  sigma <- par[2]
  return(-sum(ifelse(x==0, pnorm(-mu/sigma, log.p=TRUE), dnorm(x, mu, sigma, log=TRUE))))
}

# Gumbel
# Input:
#   par: vector of the two parameters mu and beta of a Gumbel distribution
#   x:   vector of data points used for fitting the likelihood
# Output:
#   value of the log likelihood of the rectified Gumbel distribution for the observed data points and given parameters
library(ordinal)
llrectifiedgumbel <- function(par, x){
  mu <- par[1]
  beta <- par[2]
  return(-sum(log(ifelse(x==0, pgumbel(x, location=mu, scale=beta), dgumbel(x, location=mu, scale=beta)))))
}


# Main function: AlertGS methodology

#ks_alec_test_pdistr: function to calculate the global-p value and the GO-group wise ALEC
# this function is applied to each GO group of interest individually
#input: gogroup:      names of the genes in the GO group of interest
#       namegroup:    name of the GO group
#       universe:     gene universe
#       alecs.sorted: ALECS for all the genes in the universe, named vector
#       timecutoff:   only until which time point should the ALEC be considered (due to the dealing with neg. and inf. alerts) (default 48 for WD mice)
#       alpha:        significance level for both the global test and the alert conc
#       seed:         set a seed for the permutations; allow for "outer" seed when this is done for many groups
#       RR:           number of permutations (default 1000)
#       plot:         should the GO group be plotted

#output: list containing the following elements:
#       pval_perm:       global p-value calculated on the permuted test statistics only
#       pval_gauss:      global p-value calculated based on the fitting of a rectified Gauss distribution to the permutation test statistics
#       pval_gumbel:     global p-value calculated based on the fitting of a rectified Gumbel distribution to the permutation test statistics
#       par_gauss:       parameters of the fitted rectified Gauss distribution
#       par_gumbel:      parameters of the fitted rectified Gumbel distribution
#       alert:           resulting AlertGS for the group
#       argmax:          condition value (here: time) where the maximum of the KS test statistic is attained
#       annotated:       number of genes in the group
#       nralertingroup:  number of genes with an alert in the group

ks_alec_test_pdistr <- function(gogroup,
                                namegroup = NULL,
                                universe,
                                alecs.sorted,
                                timecutoff = 48,
                                alpha = 0.05,
                                seed,
                                RR = 1000,
                                plot=TRUE){
  
  # only set a seed if this function is called individually, should be possible to have an overall seed
  # when this is called several times
  if(!seed=="outer") set.seed(seed)
  
  NN <- length(universe)
  MM <- length(gogroup)
  
  # set up the KS-matrix with entries "-MM" for the non-group members and "NN-MM" for the group members
  ks.alecs <- cbind(alecs.sorted, rep(-sqrt(MM/(NN-MM)), length(alecs.sorted))) # -sqrt(MM/(NN-MM))
  ks.alecs[gogroup,2] <- sqrt((NN-MM)/MM) # sqrt((NN-MM)/MM)
  ks.alecs <- cbind(ks.alecs, cumsum(ks.alecs[,2]))
  
  ## global p-value
  max_ks <- max(ks.alecs[which(ks.alecs[,1]<=timecutoff),3])
  ##the KS statistic always starts at 0, which is neglected in the setup here. So if the maximum is negative, set it to 0
  if(max_ks < 0) max_ks <- 0
  argmax_ks <- ks.alecs[which.max(ks.alecs[which(ks.alecs[,1]<=timecutoff),3]), 1]
  
  cumsum_permuted <- lapply(1:RR, function(k){
    alecs.perm <- matrix(0, ncol=3, nrow=length(alecs.sorted))
    
    alecs.perm[,1] <- ks.alecs[,1]
    alecs.perm[,2] <- sample(ks.alecs[,2])
    alecs.perm[,3] <- cumsum(alecs.perm[,2])
    
    # for the global p-value: which is the maximum value in the time before the cutoff?
    max_value <- max(alecs.perm[which(alecs.perm[,1]<=timecutoff),3])
    ##the KS statistic always starts at 0, which is neglected in the setup here. So if the maximum is negative, set it to 0
    if(max_value < 0) max_value <- 0
    
    # for the alert conc: calculate all the cumulative sums of the permutations (as done above)
    list(max = max_value, alert = alecs.perm[,3])
  })
  
  ## All test statistics from the permutation
  max_perm_teststat <- unlist(lapply(cumsum_permuted, function(set) set$max))
  
  ##Permutation p-value
  pval_global <- sum(max_perm_teststat >= max_ks)/RR
  
  ##Distribution approach: fit a rectified Gauss distribution to the test statistics
  optres_gauss <- try(optim(par=c(10, 20), llrectifiedgauss, x = max_perm_teststat, method="L-BFGS-B", lower=c(-Inf, 0)))
  ##Distribution approach: fit a rectified Gumbel distribution to the test statistics
  optres_gumbel <- try(optim(par=c(10, 10), llrectifiedgumbel, x = max_perm_teststat, method="L-BFGS-B", lower=c(-Inf, 0)))
  
  if(inherits(optres_gauss, "try-error")) { 
    pval_gauss <- NA 
    par_gauss <- NA
  }
  if(! inherits(optres_gauss, "try-error")) {
    pval_gauss <- ifelse(max_ks == 0, 1, pnorm(max_ks, mean=optres_gauss$par[1], sd = optres_gauss$par[2], lower.tail=FALSE))
    par_gauss <- optres_gauss$par
  }
  
  if(inherits(optres_gumbel, "try-error")) { 
    pval_gumbel <- NA 
    par_gumbel <- NA
  }
  if(! inherits(optres_gumbel, "try-error")) {
    pval_gumbel <- ifelse(max_ks == 0, 1, pgumbel(max_ks, location=optres_gumbel$par[1], scale = optres_gumbel$par[2], lower.tail=FALSE))
    par_gumbel <- optres_gumbel$par
  }
  
  cat(c(pval_global, pval_gauss, pval_gumbel))
  cat("\n")
  
  
  ## only calculate an alert conc if any of these p-values is smaller than alpha
  if(all(c(pval_global, pval_gauss, pval_gumbel) > alpha, na.rm=TRUE)){
    cat(paste("No significant global result."))
    cat("\n")
    alertconc_group <- NA
    pvals_ks <- 1 # needed only for the "if" in the plotting
  }
  if(any(c(pval_global, pval_gauss, pval_gumbel) <= alpha, na.rm=TRUE)){
    cat(paste("Significant global result"))
    cat("\n")
    
    # Calculation of the alert concentration
    permuted_ks_matrix <- sapply(cumsum_permuted, function(set){
      set$alert
    })
    
    pvals_ks <- rowSums(permuted_ks_matrix[which(ks.alecs[,1]<=timecutoff),] >= ks.alecs[which(ks.alecs[,1]<=timecutoff),3]) / RR
    
    if(any(pvals_ks < alpha)){
      cat(paste("Significant alert at condition value",ks.alecs[min(which(pvals_ks <= 0.05)), 1]))
      cat("\n")
      alertconc_group <- round(ks.alecs[min(which(pvals_ks <= 0.05)), 1], 2)
    }
    if(!any(pvals_ks < alpha)){
      cat("No alert can be calculated")
      cat("\n")
      alertconc_group <- NA
    }
    
  }
  
  
  if(plot == TRUE){
    plot(ks.alecs[,1], ks.alecs[,3], type="l", 
         main=paste(namegroup, ", M = ", MM, ", Perm. p-value = ", pval_global))
    abline(h=0, lty=2, col="red")
    abline(v=timecutoff+0.1, lty=2, col="blue")
    points(ks.alecs[,1][which.max(ks.alecs[which(ks.alecs[,1]<=timecutoff),3])], 
           max(ks.alecs[which(ks.alecs[,1]<=timecutoff),3]), pch=16, col="red")
    if(any(pvals_ks < alpha)){
      points(alertconc_group, ks.alecs[min(which(pvals_ks <= 0.05)), 3], pch=16, col="blue")
    }
  }
  cat("\n")
  
  
  return(list(pval_perm = pval_global, 
              pval_gauss = pval_gauss,
              pval_gumbel = pval_gumbel,
              par_gauss = par_gauss,
              par_gumbel = par_gumbel,
              alert = alertconc_group,
              argmax = argmax_ks,
              annotated = MM,
              nralertingroup = sum(alecs.sorted[gogroup] < timecutoff)))
}




# ks_results_table_pdistr: function to summarize results from the KS test into a result table

#input:   alertresults: result as obtained from ks_alec_test_pdistr
#         adjustment: adjustment method as passed to "p.adjust", default: fdr

#output:  Results from the KS-test in a data.frame, sorted by alert value and global p-value (permutation-based)

ks_results_table_pdistr <- function(alertresults, adjustment = "fdr"){
  out_df <- data.frame(
    GOId = names(alertresults),
    GOTerm = sapply(names(alertresults), function(go) get(go, GOTERM)@Term),
    AlertGS = sapply(alertresults, function(set) set$alert),
    PermPVal = sapply(alertresults, function(set) set$pval_perm),
    GaussPVal = sapply(alertresults, function(set) set$pval_gauss),
    GumbelPVal = sapply(alertresults, function(set) set$pval_gumbel),
    ArgMax = sapply(alertresults, function(set) set$argmax),
    Annotated = sapply(alertresults, function(set) set$annotated),
    NrAlertGroup = sapply(alertresults, function(set) set$nralertingroup)
  )
  # Adjustment
  out_df$PermPValAdjustment = p.adjust(out_df$PermPVal, method = adjustment)
  out_df$GaussPValAdjustment = p.adjust(out_df$GaussPVal, method = adjustment)
  out_df$GumbelPValAdjustment = p.adjust(out_df$GumbelPVal, method = adjustment)
  
  # sort by alert
  out_df <- out_df[order(out_df$AlertGS, out_df$PermPVal),]
  
  # sort columns of the table
  out_df <- out_df[,c(1,2,8,9,7,3,4,5,6,10,11,12)]
  return(out_df)
}
##Note: Here, an alert is present whenever at least one of the p-values is smaller than 0.05
# so in the analysis, an additional step considering the actual significances (specifically also
# after adjustment!) is required


