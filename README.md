# AlertGS

This repository contains the necessary functions to run the AlertGS methodology for an appropriately pre-processed dataset, and 
to reproduce the simulations conducted in the corresponding manuscript.

Note: in addition to the calculation of the global p-values via the rectified Gumbel distribution as discussed in the paper, 
the code also contains the version to calculate them via the permutation test statistics only, or via a rectified Gauss distribution.

The functions are contained in the files as follows:

- File `AlertGSFunctions.R' contains the main function for the AlertGS methodology, called ks_alec_test_pdistr. Two helper functions 
(llrectifiedgauss and llrectifiedgumbel) define the log-likelihoods for the rectified Gauss and Gumbel distribution. The function 
ks_results_table_pdistr can be used to store the results of the AlertGS methodology in a dataframe.
- File `SimulationFunction-NullSituation.R' contains a simulation function for a dataset under the null situation, i.e. without
any groups with meaningful group-wise alerts.
- File `SimulationFunction-AlternativeSituation' contains a simulation for a dataset under the alternative situation, i.e. where
a certain number of gene sets have a meaningful alert, which is reflected in the simulated alerts on gene-level. Both the independent
and the iterative approach are implemented in this function.

All implemented functions are commented in detail, and information on the required input parameters and the output values is given.