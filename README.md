# Test for sufficient follow-up in survival data with a cure fraction using shape constrained density estimator
Cure models are widely utilized to analyze survival data with a cure fraction. 
A sufficiently long follow-up period is necessary to estimate the proportion of 
'immune' or 'cured' subjects who will never experience failure. 
We introduce a more flexible concept of "practically" sufficient follow-up, 
defined by the distribution's quantiles, and propose a novel nonparametric 
statistical test. This method primarily assumes a non-increasing density 
function in the distribution's tail. The test utilizes a shape-constrained 
density estimator like the Grenander or kernel smoothed Grenander estimator, 
with critical values computed using a bootstrap procedure.
We also extend the test to settings with categorical covariates.


Detailed description of the method without covariate information can be found in 
Yuen and Musta (2024). The extension of the test to settings with categorical covariates
can be found in Yuen, Musta and Van Keilegom (2025).
The implementation is packed as an R package `cureSFUTest` in this repository.

## Installation
Install the package from GitHub via the [**remotes**](https://remotes.r-lib.org) package:
```R
remotes::install_github('tp-yuen/cureSFUTest')
```

## Example {.tabset}
The tabs below contain examples for test without covariate information
and with categorical covariates.

### Without covariate
After installing the package, one can run the following code using a 
breast cancer dataset from 
[**curatedBreastData**](https://bioconductor.org/packages/curatedBreastData) 
as an illustration.
```R
library(cureSFUTest)
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("curatedBreastData")
library(curatedBreastData)

# Read the data
data(curatedBreastDataExprSetList)
d <- curatedBreastDataExprSetList$study_2034_GPL96_all
df <- data.frame(Y = d$RFS_months_or_MIN_months_of_RFS, status = 1 - d$RFS)

# Test setting
tau <- 20 * 12
eps <- 0.01
n.boot <- 1000L

# Test using Smoothed Grenander estimator
sg.test.res <- sfu.test(df$Y, df$status, tau = tau, eps = eps, method = "sg",
                        n.boot = n.boot)
print(sg.test.res)

# Test using Grenander estimator
g.test.res <- sfu.test(df$Y, df$status, tau = tau, eps = eps, method = "g")
print(g.test.res)

# Hypothetical follow-up cutoff at 90 months
cutoff <- 90
df.cutoff <- df
df.cutoff[df.cutoff$Y > cutoff, ]$status <- 0
df.cutoff[df.cutoff$Y > cutoff, ]$Y <- cutoff
# Test using Smoothed Grenander estimator
sg.test.res <- sfu.test(df.cutoff$Y, df.cutoff$status, tau = tau, eps = eps, 
                        method = "sg", n.boot = n.boot)
print(sg.test.res)
```

### With categorical covariates
After installing the package, one can run the following code using a simulated 
dataset from this package as an illustration.
```R
library(cureSFUTest)

# Read the data
data("simdata")

# Test setting
tau <- 20.46742
eps <- 0.01
n.boot <- 1000L

# Test with categorical covariates using Smoothed Grenander estimator
sg.test.res <- sfu.cov.test(simdata$Y, simdata$delta, simdata$X, tau, eps, 
                            n.boot = n.boot)
print(sg.test.res)
```


## References

Planey K (2023). curatedBreastData: Curated breast cancer gene expression data with survival and treatment information. doi:10.18129/B9.bioc.curatedBreastData, R package version 2.30.0, https://bioconductor.org/packages/curatedBreastData.

Yuen, T. P., & Musta, E. (2024). Testing for sufficient follow-up in survival data with a cure fraction. 

Yuen, T. P., Musta, E., & Van Keilegom, I. (2025). Testing for sufficient follow-up in survival data with categorical covariates.
