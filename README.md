
# fire <img src="man/figures/logo.png" align="right" height="120" />

**An R Package for Functional I-prior Regression using RKHS Norm**

[![License: GPL
v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)

## Overview

The **fire** package implements the *Functional I-prior Regression*
(FIRE) model in R.  
Unlike traditional regression models that rely on Euclidean distance,
**fire** uses the reproducing kernel Hilbert space (RKHS) norm to
measure dissimilarity. This approach provides greater flexibility for
modeling nonlinear and structured covariates such as vectors, matrices,
and higher-order tensors.

Key features include:

- Regression with vector, matrix, and tensor covariates  
- RKHS norm kernels (fractional Brownian motion, RBF, polynomial,
  etc.)  
- Efficient hyperparameter estimation via the EM algorithm  
- Built-in cross-validation (`fire_cv`)  
- S3 methods for `print()`, `summary()`, `fitted()`, `predict()`, and
  `plot()`  
- Example datasets: **Manure**, **Housing** and **Sugar**

## Installation

You can install the development version from GitHub:

``` r
install.packages("devtools")
devtools::install_github("ZiqingHo/fire")
library(fire)
```
