MetaNorm
========================
To carry out a normalization procedure for a 
NanoString nCounter dataset, you will need the
`MetaNorm` function. We have also included a sample
data used by RCRnorm for demonstration `normalization_data.RData`. 

.. note:: 
    Please make sure that your datasets conform to our format

To load the data and carry out `MetaNorm`, simply run the following: 

.. code:: bash 
    
    data("normalization_data")
    draw = MetaNorm(dat=normalization_data, M=5000, n_keep=1000)

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
across all non-burn-in draws.