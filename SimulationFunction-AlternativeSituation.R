#Function to simulate alerts for a gene expression dataset, where the genes are assigned to 
#   GO groups, some of which are assumed to have a specific alert

# Input:
#   GOlist:               List of GO groups, with the names of the respectively annotated genes
#   allgenenames:         Vector containing all names of genes for which an alert is to be calculated
#   nrmeaningfulgroups:   Nr of groups with a simulated alert (denoted by n in the paper), defaults to 10
#   size_min:             Minimal size of GO group that can be selected to be a meaningful group, defaults to 20
#   size_max:             Maximal size of GO group that can be selected to be a meaningful group, defaults to 200
#   propmeaninggenes:     Ratio of genes within the GO group with a meaningful alert (denoted by x in the paper), defaults to 0.3
#   betarange:            Width of the support of the beta distribution used for sampling meaningful alerts, defaults to 4
#   propunifgenes:        Ratio of genes for which a uniform alert is randomly drawn in the last step (denoted by r in the paper), defaults to 0.1
#   variant:              Which of the simulation approaches, "independent" or "iterative" should be chosen
#   anyuniform:           Boolean variable to decide whether uniform alerts should be sampled for the genes in the last step. Defaults to Yes

# Output: List containing the following elements:
#   MeaningfulGroups:     Vector containing the names of the GO groups sampled to be "meaningful" (i.e., have an alert)
#   AlecsSorted:          Sorted vector of simulated Alecs for each gene. Note: this includes also the replacement 
#                         of missing values with values outside the observed range, and the solving of ties
#   Simresults:           Results of the application of the AlertGS methodology
#   BackgInfo:            Matrix containing the following information for each meaningful GO group:
#                         - indivalecmean: sampled true group-wise alert (denoted as g in the paper)
#                         - alecproportion: proportion of genes with alec for each meaningful group
#                         - requiredindex: relevant only in the "independent" approach, gives the number of genes that should be 
#                                          sampled from all genes without alert from this group
#                         - availablenrgenes: relevant only in the "independent" approach, gives the numbers of genes without alerts
#                                             in cases where more alerts are to be sampled than genes are available
#   OverallNrMeaning:     The total number of genes with "meaningful" alert, i.e. not considering genes with an alert purely due
#                         to the uniform sampling
#   OverallNr:            The total number of genes with alert, including those with alert sampled from the uniform distribution
#   MoreIndexThanGenes:   Relevant only in the "Independent" approach, states the number of meaningful groups for which the case occurred
#                         that more genes were supposed to be sampled from the group than were available


simulation_fun <- function(GOlist, 
                           allgenenames, 
                           nrmeaningfulgroups = 10,
                           size_min = 20, 
                           size_max=200,
                           propmeaninggenes = 0.3,
                           betarange = 4,
                           propunifgenes = 0.1,
                           variant,
                           anyuniform = TRUE){
  
  # 0) read out the sizes of the GO groups
  gosizes <- sapply(GOlist, function(set) length(set))
  
  # 0) initialize a vector of the length of all considered genes and assign an alec of NA
  alecs <- rep(NA, length(allgenenames))
  names(alecs) <- allgenenames
  
  # 0) initialize a counter, how often (only relevant in the independent scenario), more genes should obtain an alert than available
  moreindexthangenes <- 0
  
  # Variant: Iterative
  if(variant == "iterative"){
    # 1) randomly choose nrmeaningfulgroups GO groups of sizes between size_min and size_max as `group with meaningful alert'
    meaningfulgroups <- sample(names(GOlist)[which(gosizes >= size_min & gosizes <= size_max)], nrmeaningfulgroups)
    
    # 2) Sample the true alerts and sort them in increasing order
    indivalecmean <- sample(5:46, nrmeaningfulgroups, replace=TRUE)
    indivalecmean <- sort(indivalecmean)
    
    # 3) iteratively, draw propmeaninggenes from the entire GO group. If they do not have an alert yet, draw a Beta(2, 2)-distributed alert
    for(i in 1:nrmeaningfulgroups){
      # consider the genes from the i-th meaningful group
      genes_group <- unlist(GOlist[meaningfulgroups[i]])
      # randomly sample the respective proportion of genes 
      index <- sample(genes_group, size=ceiling(propmeaninggenes*length(genes_group)))
      # find out which do not have an alert yet
      genes_group_noalert <- intersect(index, names(alecs)[is.na(alecs)])
      # sample alecs for these
      alecs[genes_group_noalert] <- 
        round((rbeta(length(genes_group_noalert), 2, 2)*(betarange)+(indivalecmean[i]-betarange/2)), 4)
    }
    
    # find out the proportion of genes with alec for each meaningful group
    alecproportion <- sapply(1:nrmeaningfulgroups, function(i){
      genes_group <- unlist(GOlist[meaningfulgroups[i]])
      round(sum(!is.na(alecs[genes_group]))/length(genes_group), 3)
    })
    
    # find out the overall proportion of genes with meaningful alec
    overall_alecnr_meaning <- sum(!is.na(alecs))
    
    ## these variables are only needed in the "independent"-approach, but in order to be able to save the same objects,
    ##  they are initialized here anyway
    requiredindex <- rep(NA, nrmeaningfulgroups)
    availablenrgenes <- rep(NA, nrmeaningfulgroups)
    
  }
  
  # Variant: Independent
  if(variant == "independent"){
    # 1) randomly choose nrmeaningfulgroups GO groups of sizes between size_min and size_max as `group with meaningful alert'
    meaningfulgroups <- sample(names(GOlist)[which(gosizes >= size_min & gosizes <= size_max)], nrmeaningfulgroups)
    
    # 2) Sample the true alerts, these remain in random order
    indivalecmean <- sample(5:46, nrmeaningfulgroups, replace=TRUE)
    
    # 3) draw propmeaninggenes from the entire GO group, but sample only out of those who do not yet have an alert.
    #    in case the number of genes without alert is not large enough, use all remaining genes
    requiredindex <- rep(NA, nrmeaningfulgroups)
    availablenrgenes <- rep(NA, nrmeaningfulgroups)
    for(i in 1:nrmeaningfulgroups){
      # consider the genes from the i-th meaningful group
      genes_group <- unlist(GOlist[meaningfulgroups[i]])
      # calculate the number of genes that should obtain an alert
      nralerts <- ceiling(propmeaninggenes*length(genes_group))
      # find out which genes do not have an alert yet
      genes_group_noalert <- intersect(genes_group, names(alecs)[is.na(alecs)])
      
      # if the number of genes that should obtain an alert is larger than the number of remaining genes without alert, use all genes
      if(nralerts > length(genes_group_noalert)){
        print("Nr of alerts larger than genes without alert")
        moreindexthangenes <- moreindexthangenes + 1
        requiredindex[i] <- nralerts
        availablenrgenes[i] <- length(genes_group_noalert)
        index <- genes_group_noalert
      }
      # else
      if(nralerts <= length(genes_group_noalert)){
        index <- sample(genes_group_noalert, size=nralerts)
      }
      
      # sample alecs for the selected genes
      alecs[index] <- round((rbeta(length(index), 2, 2)*(betarange)+(indivalecmean[i]-betarange/2)), 4)
    }
    
    # find out the proportion of genes with alec for each meaningful group
    alecproportion <- sapply(1:nrmeaningfulgroups, function(i){
      genes_group <- unlist(GOlist[meaningfulgroups[i]])
      round(sum(!is.na(alecs[genes_group]))/length(genes_group), 3)
    })
    
    # find out the overall proportion of genes with meaningful alec
    overall_alecnr_meaning <- sum(!is.na(alecs))
    
  }
  
  # Additionally have uniform genes
  if(anyuniform){
    # 4) sample uniform genes, across all genes, but only replace if the uniformly sampled alert is smaller
    index_uniform <- sample(names(alecs), ceiling(propunifgenes*length(alecs)))
    alecs[index_uniform] <- pmin(alecs[index_uniform], round(runif(length(index_uniform), min=3, max=48), 4), na.rm=TRUE)
  }
  
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
  
  backgroundinfo_meaningful <- cbind(indivalecmean, alecproportion, requiredindex, availablenrgenes)
  
  
  return(list(MeaningfulGroups = meaningfulgroups, 
              AlecsSorted = alecs_sorted,
              Simresults = alertresults_sim, 
              BackgInfo = backgroundinfo_meaningful,
              OverallNrMeaning = overall_alecnr_meaning,
              OverallNr = overall_alecnr,
              MoreIndexThanGenes = moreindexthangenes))
}