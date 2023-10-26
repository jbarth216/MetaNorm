# MetaNorm: Incorporating Meta-analytic Priors into Normalization of NanoString nCounter Data

![Logo](/assets/logo.png)

MetaNorm is a normalization procedure for Nanostring nCounter datasets. 

Please refer to our paper for more details: [MetaNorm link here](www.google.com)

Besides the R helper as well as the vignette we provide along with the package, we have also built a detailed [online documentation](https://pmtnet-omni-document.readthedocs.io/en/latest/) where we guide you step-by-step on how to conduct meta analysis as well as normalizing a NanoString nCounter dataset.

## Model Overview 
![Model Overview](/assets/model.png)

## Dependencies 
- R
- Rcpp
- 

## Installation
```shell
intall_github(????)
```

## Quick Start Guide 


MetaNorm is a normalization technique for Nanostring nCounter datasets. The artical associated with this project is currently under review at Bioinformatics.

There are two main projects associated with this analysis. They are described in detail below.

#1) Meta-analysis of Nanostring nCounter datasets

"Meta_Analysis_Clean_20230825.R" contains the code to perform the meta-analysis. Positive probe data for the 13 collected datasets can be found in 
"All_clean2.csv", while estimated coefficients based on this data can be found in "coeffs2.csv". Finally, "armadillo.cpp" contains C++ code used in 
the meta-analysis (via rcpp package) to improve performance speed.

The purpose of the meta-analysis is to provide posterior estimates to plug in to prior distributions of model parameters in the MetaNorm procedure.
The analysis is based on a complex Bayesian hierarchical model, similar to the ones used in MetaNorm and RCRnorm. The model is designed specifically
for these datasets, and while not mean to be reproduced with other data, can certainly be a guide for similar analyses.

#2) MetaNorm

"MetaNorm_Clean_20230825.R" is the R code containing the normalization function. Nested in this is the file "Loop_functions.R", which contains update
functions for various model parameters. Please note that MetaNorm takes 5 inputs, listed in detail below:

1) dat
A list of 4 dataframes containing:
'$pos_dat', a dataframe of positive probe read counts. Must have 6 rows and n columns (where n is the number of samples)
'$neg_dat', a dataframe of negative probe read counts. Can have any number of rows (for varying negative probe amounts) and n columns
'$hk_dat', a dataframe of housekeeping gene read counts (specified by researcher(s)). Can have any number of rows and n columns
'$reg_dat', a dataframe of "regular" (i.e., all other) gene read counts. Can have any number of rows and n columns

Any deviation from this format will result in an error, or faulty results. Please see FFPE_dat from the RCRnorm R package for a properly formatted example.

2) pos_conc = log10(c(128, 32, 8, 2, 0.5, 0.125))
A set of known RNA values for each positive probe. Unless experimenting with different values, these should always be left at the default values.

3) M=15000 
The number of draws for each parameter, inclduing burn-in. The default of 15,000 is conservative, as MetaNorm can produce stable estimates often with
only 2-3k draws total (see article for more details)

4) mm = 3 
The radius (in standard deviations) for parameters with uniform priors (see article for more details)

5) seed=0
Random seed, for reproducibility. If running multiple chains to check convergence/stability, make sure that each chain has a different seed.

OUTPUT: The output of MetaNorm is a list of vectors, dataframes or arrays of parameter draw. At present, all draws are kept (including burn-in). Contained
in the output list are:

'$a': matrix of n columns (one for each intercept parameter) and M rows (for each draw)
'$b': matrix of n columns and M rows
'$cc': vector of size M (represents the "c" parameter used in the negative probe equation. See article)
'$mu_a': vector of size M
'$mu_b': vector of size M

...

For the rest of the output sets, the indices are relatively straightforward and follow the same logic as these 5. For parameters indexed both by sample
and gene (i.e., kappa_reg), the 1st dimension represents the number of rows, the 2nd represents the number of genes, the 3rd represents the number of samples.

Please note that the "normalized" parameter estimates are contained in "kappa_reg" and "kappa_hk" and will need to be aggregated over the first dimension
across all non-burn-in draws.


FUTURE UPDATES: The next iteration of this code will provide posterior estimates automatically, and will give users the options not to keep all draws to 
preserve RAM and minimize the chances of R crashing.

## Citation
The artical associated with this project is currently under review at Bioinformatics.

## Contact 
Jackson Barth, PhD (jackson_barth@baylor.edu)

Assistant Professor

Department of Statistical Science

Baylor University 

