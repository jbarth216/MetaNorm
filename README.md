# MetaNorm: Incorporating Meta-analytic Priors into Normalization of NanoString nCounter Data

![Logo](/assets/logo.png)

MetaNorm is a normalization procedure for Nanostring nCounter datasets. 

Please refer to our paper for more details: [MetaNorm link here](www.google.com)

Besides the R helper as well as the vignette we provide along with the package, we have also built a detailed [online documentation](https://metanorm.readthedocs.io/en/latest/) where we guide you step-by-step on how to conduct meta analysis as well as normalizing a NanoString nCounter dataset.

## Model Overview 
![Model Overview](/assets/model.png)

## Dependencies 
- R>=4.0.2
- Rcpp>=1.0.10
- RcppArmadillo>=0.12.4
- mvtnorm>=1.1
- MASS>=17.3
- truncnorm>=1.0
- progress>=1.2.2

## Installation
```shell
intall_github(????)
```

## Quick Start Guide 
Along with the package, we have provided two example datasets ``meta_analysis_data.RData`` and ``normalization_data.RData``. We will use them to quickly demonstrate how to get started with our package. 

### Meta-analysis of Nanostring nCounter datasets
In the ``meta_analysis_data.RData``, we provide positive probe data for the 13 collected datasets. To curate the data and generate empirical estimated coefficients based on this data simply run the following. 
```shell 
library(MetaNorm)
data("meta_analysis_data")
ds = curate_data(dataset=ds)
results = find_regression_coefs(df=ds)
```
The purpose of the meta-analysis is to provide posterior estimates to plug in to prior distributions of model parameters in the MetaNorm procedure.

The analysis is based on a complex Bayesian hierarchical model, similar to the ones used in MetaNorm and RCRnorm. The model is designed specifically
for these datasets, and while not mean to be reproduced with other data, can certainly be a guide for similar analyses.

To run the Gibbs sampler, use the following code:
```shell
draw = meta_analysis(results$df, results$coeffs2, M=12000, n_keep=-1)
```
The result is a dataframe containing posterior sample draws of major variables. `M` specifies how many samples to draw while `n_keep` tells the program how many samples to keep. A negative number means to save all the samples. 

To run several MCMC chains, simply repeat the previous code several times. 
```shell 
draw1 = meta_analysis(results$df, results$coeffs2, M=12000, n_keep=-1)
draw2 = meta_analysis(results$df, results$coeffs2, M=12000, n_keep=-1)
```

### MetaNorm
To carry out a normalization procedure for a NanoString nCounter dataset, you will need the `MetaNorm` function. We have also included a sample data used by RCRnorm for demonstration `normalization_data.RData`. 

*Please make sure that your datasets conform to our format*. Once you load the dataset
```shell 
data("normalization_data.RData")
```
you will see the data is a list of 4 dataframes containing:
- '$pos_dat', a dataframe of positive probe read counts. Must have 6 rows and n columns (where n is the number of samples)
- '$neg_dat', a dataframe of negative probe read counts. Can have any number of rows (for varying negative probe amounts) and n columns
- '$hk_dat', a dataframe of housekeeping gene read counts (specified by researcher(s)). Can have any number of rows and n columns
- '$reg_dat', a dataframe of "regular" (i.e., all other) gene read counts. Can have any number of rows and n columns

*Any deviation from this format will result in an error, or faulty results.* 

Once you have prepared your dataset, to perform `MetaNorm` simply run:
```shell
draw = MetaNorm(dat=normalization_data, M=5000, n_keep=1000)
```
As with the `meta_analysis` function, `M` is the number of draws for each parameter, inclduing burn-in. The default of 15,000 is conservative, as MetaNorm can produce stable estimates often with
only 2-3k draws total. 

The other two parameters are 
- `mm` (default is 3)
The radius (in standard deviations) for parameters with uniform priors
- `seed` (default is 0)
Random seed, for reproducibility. If running multiple chains to check convergence/stability, make sure that each chain has a different seed.

For more details, please refer to our [online documentation](https://metanorm.readthedocs.io/en/latest/)

## Citation
The artical associated with this project is currently under review at Bioinformatics.

## Contact 
Jackson Barth, PhD (jackson_barth@baylor.edu)

Assistant Professor

Department of Statistical Science

Baylor University 

