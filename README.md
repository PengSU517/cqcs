# CQCS: Censored quantile correlation based screening

The `cqcs` package provides the functions to compute the CQC (censored quantile correlation) proposed by Peng Su, Ruoxuan Zhang and Maozai Tian.

We have included a simu (Simulation) in vignettes.


## Key functions
There are two key functions in the package.

function `cqcor` is used to caculate the CQC between two variables, which can be used to screen the variables.

```{r}
cqcor(y, d, x, tau) 
# y: response variable, 
# d: censoring indicator, 
# x: covariate, 
# tau: quantile level
```

function `cqcs` is used to screen the variables based on CQC.

```{r}
cqcs(y, d, x, tau = 0.5, len = 100) 
# y: response variable, 
# d: censoring indicator, 
# x: covariates, 
# tau: quantile level, 
# len: screening length
```



## Installation

You can install the package using:

```r
remotes::install_github("PengSU517/cqcs")
```




