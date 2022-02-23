# CausalHR
 
This is a repository for the CausalHR package, used to estimate the kernel-based and the Cox-based causal hazard ratio estimator. This package accompanies the following paper "Sensitivity analysis for the causal hazard ratio in randomized and observational studies" by Rachel Axelrod and Daniel Nevo (arXiv link pending).

This package is based on a fork of the muhaz package from CRAN. We updated the muhaz R function and Fortran code to allow observation weights to be used when estimation the hazard function. The muhaz package was originally authored by Kenneth Hess, R. Gentleman and David Winsemius, and is available [on CRAN](https://cran.r-project.org/web/packages/muhaz/index.html).

## Installation

For latest development version run:

```{r}
install.packages("devtools")
devtools::install_github("xlrod1/CausalHR")
```

##How to use the CausalHR package?
There is two main function in this package that one should use when analysing real data set: 
1. CausalHR.with.bootstrap- Returns the kernel-based and the Cox-based Causal hazard ratio estimators. It also returns the standard errors estimated by the bootstrap, and 95\% pointwise CIs were calculated by the percentile method.
2. CausalHR.without.bootstrap-Returns the kernel-based and the Cox-based Causal hazard ratio estimators.

##Where can I read more?
"Sensitivity analysis for the causal hazard ratio in randomized and observational studies" by Rachel Axelrod and Daniel Nevo (arXiv link pending).

##Reproducing paper results
A github repository for reproducing paper results paper is found [here](https://github.com/xlrod1/Reproducability_R_code)
