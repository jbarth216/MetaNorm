Meta-analysis of Nanostring nCounter datasets
==========================================================
The purpose of the meta-analysis is to provide posterior estimates 
to plug in to prior distributions of model parameters in the MetaNorm 
procedure. The analysis is based on a complex Bayesian hierarchical model, 
similar to the ones used in MetaNorm and RCRnorm. The model is designed 
specifically for these datasets, and while not mean to be reproduced with 
other data, can certainly be a guide for similar analyses.

In the ``meta_analysis_data.RData``, we provide positive probe
 data for the 13 collected datasets. To curate the data and generate
empirical estimated coefficients based on this data simply run the following. 

.. code:: bash 
    library(MetaNorm)
    data("meta_analysis_data")
    ds = curate_data(dataset=ds)
    results = find_regression_coefs(df=ds)

The minimal requirement is that the data contains 5 columns which must be named 

#. DataSet: An unique ID given to each study
#. RNA: The designated mRNA measurement of the positive probes. They must be 128, 32, 8, 2, 0.5, and 0.125. 
#. SampleID: An unique ID given to each patient in each study 
#. Count: The actual measurements 
#. UID: An unique ID given to each **combination** of patient and study

The data curation involves two steps: creating 
indices that are consecutive and performing linear 
regression to get empirical estimates of intercepts, 
slopes, and residuals. 

Once these are done, you are ready to perform meta analysis

.. code:: bash 
    Draws = meta_analysis(ds=results$df,
                      coeffs2=results$coeffs2,
                      M=12000,
                      n_keep=5000)





