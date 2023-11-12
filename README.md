# MetaNorm: Incorporating Meta-analytic Priors into Normalization of NanoString nCounter Data

![Logo](/assets/logo.png)

MetaNorm is a normalization procedure for Nanostring nCounter datasets. 

Please refer to our paper for more details: [MetaNorm link here](www.google.com)

Besides the R helper as well as the vignette we provide along with the package, we have also built a detailed [online documentation](https://metanorm.readthedocs.io/en/latest/) where we guide you step-by-step on how to conduct meta analysis as well as normalizing a NanoString nCounter dataset.

## Model Overview 
![Model Overview](/assets/model.png)

## Dependencies 

- R>=4.0.2

## Installation
```shell
library(devtools)
install_github("Yuqiu-Yang/MetaNorm", subdir="MetaNorm", ref="main")
```
Other than making sure your R version is correct, there is no need for 
you to manually install other packages as they will be automatically 
installed when installing our package.  
<details>
<summary>Other dependencies</summary>

1. Rcpp>=1.0.10
2. RcppArmadillo>=0.12.4
3. mvtnorm>=1.1
4. MASS>=17.3
5. truncnorm>=1.0
6. progress>=1.2.2
</details>


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
Once you plot the posterior draws, you will see something like this:
![MCMC](/docs/source/images/meta_mcmc.png)


### MetaNorm
To carry out a normalization procedure for a NanoString nCounter dataset, you will need the `MetaNorm` function. We have also included a sample data used by RCRnorm for demonstration `normalization_data.RData`. 

*Please make sure that your datasets conform to our format*. Once you load the dataset
```shell 
data("normalization_data")
```
you will see the data is a list of 4 dataframes containing:
- '$pos_dat', a dataframe of positive probe read counts. Must have 6 rows and n columns (where n is the number of samples)
- '$neg_dat', a dataframe of negative probe read counts. Can have any number of rows (for varying negative probe amounts) and n columns
- '$hk_dat', a dataframe of housekeeping gene read counts (specified by researcher(s)). Can have any number of rows and n columns
- '$reg_dat', a dataframe of "regular" (i.e., all other) gene read counts. Can have any number of rows and n columns

*Any deviation from this format will result in an error, or faulty results.* 

Once you have prepared your dataset, to perform `MetaNorm` simply run:
```shell
draw = MetaNorm(dat=normalization_data, M=5000, all_draws=TRUE)
```
As with the `meta_analysis` function, `M` is the number of draws for each parameter, inclduing burn-in. The default of 15,000 is conservative, as MetaNorm can produce stable estimates often with
only 2-3k draws total. `all_draws=TRUE` tells the algorithm to keep all the MCMC draws for all parameters (which is A LOT). 

For most users, the only parameters of interest are the posterior means for `kappa_hk` and `kappa_reg`. Of course, you are welcome to compute the posterior summary statistics by yourself based on the MCMC draws. However, as the data structure is rather convoluted, we provide along with the `MetaNorm` function additional parameters to only return the parameters of interest. To enable this mode of inference, simply set `all_draws=FALSE` and provide some basic information on how many sample are for burn in and how to thin the MCMC draws.  
```shell 
output = MetaNorm(dat=input, M=5000, all_draws=FALSE, burn_in=1000, thin=2)
```
By setting `burn_in=1000`, we are telling the algorithm to discard the first 1000 samples when computing the summary statistics. `thin=2` tells the algorithm to only collect every other sample to reduce auto-correlation. 

<details>
<summary>How do we compute the posterior summary statistics?</summary>
If you are curious about how we summarized the MCMC draws for `kappa_hk` and `kappa_reg` or you simply want to try it yourself, this is how we arrived at the summary statistics you see. 

Starting from the posterior sample list: `draw`, we first discard burnt in samples and thin the rest of the draws. 
```shell 
draws <- seq(burn_in+1, M, by=thin)
```
To compute the summary statistics
```shell
kappa_reg_stat <- apply(draw$kappa_reg[draws,,],c(2,3),mean)
kappa_hk_stat <- apply(draw$kappa_hk[draws,,],c(2,3),mean)
```

</details>

#### Low Memory Mode
Currently, due to the design of the algorithm, even if we are only returning the posterior means, we still store all the MCMC sample draws internally. This would consume lots of memory space which would be a serious issue if we are working with a large dataset. To alleviate the issue, we provide a function: `MetaNorm_LowMemory` that will only keep the most recent MCMC draws for all parameters except for global parameters. Along with the posterior samples for global parameters, we will also compute the summary statistics for `kappa_hk` and `kappa_reg`. To use this function, run the following code 
```shell
output = MetaNorm_LowMemory(dat=dat,  M=5000, burn_in = 1000, thin=2)
```
The output is a list containing two items: ``GlobalParameters_Draws`` (all draws for mu_a, mu_b, sig2_a, sig2_b) and ``Posterior_estimates`` (list of kappa_hk, kappa_reg posterior estimates).

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

