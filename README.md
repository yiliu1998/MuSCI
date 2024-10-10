# MuSCI
`MuSCI` is an R package that implements multi-source conformal inference under distribution shift. Please refer to this main paper [Liu et al. (2024)](https://proceedings.mlr.press/v235/liu24ag.html) published at PMLR, ICML 2024, for the methodology part. Please cite the paper if you use our package. 

## Installation
To install the latest version of the R package from GitHub, please run following commands in R:

```r
if (!require("devtools")) install.packages("devtools")
devtools::install_github("yiliu1998/MuSCI")
```

## Usage
An Rmd file for illustration ([click here](https://github.com/yiliu1998/MuSCI/tree/main/vignettes)) is available to download, which gives an example of the use of our package using a simulated dataset. Please run it on your local after you install the package. 

## Contact 
The R code is maintained by Yi Liu (Please feel free to reach out at yi.liu.biostat@gmail.com, if you have any questions). 

## Reference
Please cite the following paper:

Yi Liu, Alexander Levis, Sharon-Lise Normand, Larry Han. Proceedings of the 41st International Conference on Machine Learning, Proceedings of Machine Learning Research, 235:31344-31382, 2024.
