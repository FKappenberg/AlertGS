#Function to simulate alerts for a gene expression dataset, where the genes are assigned to 
#   GO groups, in the null situation, i.e. without group-specific alerts

# Input:
#   GOlist:               List of GO groups, with the names of the respectively annotated genes
#   allgenenames:         Vector containing all names of genes for which an alert is to be calculated
#   propunifgenes:        Ratio of genes for which a uniform alert is randomly drawn in the last step (denoted by r in the paper), defaults to 0.1

# Output: List containing the following elements:
#   AlecsSorted:          Sorted vector of simulated Alecs for each gene. Note: this includes also the replacement 
#                         of missing values with values outside the observed range, and the solving of ties
#   Simresults:           Results of the application of the AlertGS methodology
#   OverallNr:            The total number of genes with alert, including those with alert sampled from the uniform distribution


simulation_fun_null <- function(GOlist,
                                allgenenames,
                                propunifgenes = 0.1){
  
  # 0) read out the sizes of the GO groups
  gosizes <- sapply(GOlist, function(set) length(set))
  
  # 0) initialize a vector of the length of all considered genes and assign an alec of NA
  alecs <- rep(NA, length(allgenenames))
  names(alecs) <- allgenenames
  
  # 4) sample uniform genes, across all genes, but only replace if the uniformly sampled alert is smaller
  index_uniform <- sample(names(alecs), ceiling(propunifgenes*length(alecs)))
  alecs[index_uniform] <- pmin(alecs[index_uniform], round(runif(length(index_uniform), min=3, max=48), 4), na.rm=TRUE)
  
  overall_alecnr <- sum(!is.na(alecs))
  
  # 5) replace the NAs with values between 48.1 and 95.9
  alecs[which(is.na(alecs))] <- round(runif(sum(is.na(alecs)), 48.1, 95.9), 4)
  # deal with ties
  alecs[duplicated(alecs)] <- alecs[duplicated(alecs)]  + rnorm(sum(duplicated(alecs)), 0, 0.001)
  
  # 6) sort the alecs
  alecs_sorted <- sort(alecs)
  
  # Apply the AlertGS methodology, implemented in the function ks_alec_test_pdistr
  alertresults_sim <- lapply(1:length(GOlist), function(i){
    print(i)
    ks_alec_test_pdistr(gogroup= GOlist[[i]], 
                        namegroup = names(GOlist)[i],
                        universe = allgenenames,
                        alecs.sorted = alecs_sorted,
                        timecutoff = 48,
                        alpha = 0.05,
                        seed = "outer",
                        RR = 1000,
                        plot=FALSE)
  })
  names(alertresults_sim) <- names(GOlist)
  
  return(list(AlecsSorted = alecs_sorted,
              Simresults = alertresults_sim,
              OverallNr = overall_alecnr))
}
