MetaNorm
========================
To carry out a normalization procedure for a 
NanoString nCounter dataset, you will need the
`MetaNorm` function. We have also included a sample
data used by RCRnorm for demonstration `normalization_data.RData`. 

.. note:: 
    Please make sure that your datasets conform to our format

To load the data and carry out `MetaNorm`, simply run the following: 

.. warning:: 
    Running the following code will consume lots of memory. 
    Proceed with caution. If you are working with a large dataset, 
    please proceed directly to `Low Memory Mode`_.

.. code:: bash 
    
    data("normalization_data")
    draw = MetaNorm(dat=normalization_data, M=5000, all_draws=TRUE)

As with the `meta_analysis` function, `M` is the number of draws for each parameter, inclduing burn-in. The default of 15,000 is conservative, as MetaNorm can produce stable estimates often with
only 2-3k draws total. `all_draws=TRUE` tells the algorithm to keep all the MCMC draws for all parameters (which is A LOT). 

The output of MetaNorm is a list of vectors, dataframes or arrays of parameter draw. 
Contained in the output list are:

# a: matrix of n columns (one for each intercept parameter) and M rows (for each draw)
# b: matrix of n columns and M rows
# cc: vector of size M (represents the "c" parameter used in the negative probe equation. See article)
# mu_a: vector of size M
# mu_b: vector of size M

...

For the rest of the output sets, the indices are relatively straightforward and follow the same logic as these 5. 
For parameters indexed both by sample
and gene (i.e., kappa_reg), the 1st dimension represents the number of rows, 
the 2nd represents the number of genes, the 3rd represents the number of samples.

Please note that the "normalized" parameter estimates are contained in "kappa_reg" 
and "kappa_hk" and will need to be aggregated over the first dimension
across all non-burn-in draws, which is detailed as follow:

For most users, the only parameters of interest are the posterior means for `kappa_hk` and `kappa_reg`. 
Of course, you are welcome to compute the posterior summary statistics by yourself based on the MCMC draws. 
However, as the data structure is rather convoluted, we provide along with the `MetaNorm` function additional 
parameters to only return the parameters of interest. To enable this mode of inference, simply set `all_draws=FALSE` 
and provide some basic information on how many sample are for burn in and how to thin the MCMC draws.  

.. code:: bash 

    output = MetaNorm(dat=input, M=5000, all_draws=FALSE, burn_in=1000, thin=2)

By setting `burn_in=1000`, we are telling the algorithm to discard the first 1000 samples when computing the summary statistics. `thin=2` 
tells the algorithm to only collect every other sample to reduce auto-correlation. 

If you are curious about how we summarized the MCMC draws for `kappa_hk` and `kappa_reg` or you simply want to try it yourself, 
this is how we arrived at the summary statistics you see. 

Starting from the posterior sample list: `draw`, we first discard burnt in samples and thin the rest of the draws. 

.. code:: bash 
    
    draws <- seq(burn_in+1, M, by=thin)

To compute the summary statistics

.. code:: bash 
   
    kappa_reg_stat <- apply(draw$kappa_reg[draws,,],c(2,3),mean)
    kappa_hk_stat <- apply(draw$kappa_hk[draws,,],c(2,3),mean)

Low Memory Mode
----------------------
Currently, due to the design of the algorithm, even if we are only 
returning the posterior means, we still store all the MCMC sample draws internally.
This would consume lots of memory space which would be a serious issue if we are working with 
a large dataset. To alleviate the issue, we provide a function: `MetaNorm_LowMemory` that will only 
keep the most recent MCMC draws for all parameters except for global parameters. Along with the posterior 
samples for global parameters, we will also compute the summary statistics for `kappa_hk` and `kappa_reg`. 
To use this function, run the following code 
 
 .. code:: bash 

    data("normalization_data")
    output = MetaNorm_LowMemory(dat=dat,  M=5000, burn_in = 1000, thin=2)

The output is a list containing two items: ``GlobalParameters_Draws`` (all draws 
for mu_a, mu_b, sig2_a, sig2_b) and ``Posterior_estimates`` (list of kappa_hk, 
kappa_reg posterior estimates).


