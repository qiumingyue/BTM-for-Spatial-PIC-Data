# BTM-for-Spatial-PIC-Data
## Project description

The objective of this project is to build functions and algorithms for the semiparametric transformation model with spatial partly interval-censored data, which is  proposed in the Article ```Bayesian Transformation Model for Spatial Partly Interval-censored Data```. Several Julia functions have been created for using the model. 

## Datasets
```adjacency matrix.csv```: dataset of the adjacency matrix


```simulated data--PH.csv```: simulated dataset with the true value of $\beta$ is $(0.5, 1)'$:

 1.  The first two columns are the covariates
 2.  The third column is the left point of the time interval
 3.  The fourth column is the right point of the time interval
 4.  The fifth column is the censoring index
 5.  The last column is the location


## Files
```Functions.jl``` includes some functions created from the model: spline, sampler, likelihood functions and MCMC algorithms, etc. 


```analysis. jl``` includes some hyper-parameters settings and instructions on how to run the code.
