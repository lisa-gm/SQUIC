# SQUIC - Sparse Quadratic Inverse Covariance Estimation

This repository hosts the SQUIC algorithm, made available as an R package.

### Description

SQUIC tackles the statistical problem of estimating large sparse inverse covariance matrices. This estimation poses an ubiquitous problem that arises in many applications e.g. coming from the fields mathematical finance, geology and health. SQUIC belongs to the class of second-order L<sub>1</sub>-regularized Gaussian maximum likelihood methods and is especially suitable for high-dimensional datasets with limited number of samples. For further details please see the listed references.

### References

<a id="1">[1]</a> 
Bollhöfer, M., Eftekhari, A., Scheidegger, S. and Schenk, O., 2019. 
Large-scale sparse inverse covariance matrix estimation. 
SIAM Journal on Scientific Computing, 41(1), pp.A380-A401.

<a id="1">[2]</a> 
Eftekhari, A., Bollhöfer, M. and Schenk, O., 2018, November. 
Distributed memory sparse inverse covariance matrix estimation on high-performance computing architectures. 
In SC18: International Conference for High Performance Computing, Networking, Storage and Analysis (pp. 253-264). IEEE.

### Installation

The installation can be done directly from R, e.g. using the devtools package:

library(devtools)
install_github("lisa-gm/SQUIC_R")
