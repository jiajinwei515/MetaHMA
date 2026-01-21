MetaHMA is a method for mediation analysis with both exposures and mediators being high-dimensional, where the outcome is binary.
Before implementing MetaHMA, it is necessary to load the R functions in helper.

Specifically, "decorr.R" includes the R code for the deconfounded/decorrelated and debiased high-dimensional linear regression proposed by Sun et al. (2024), 
where the source code is available at https://github.com/YinruiSun/MT_confounding.

The R file "ssdblasso.R" includes the R code for the debiased lasso after sample splitting (Vazquez and Nan, 2025),
where the source code is available at  https://github.com/omar-vazquez/samplesplit_glm. We modified the code to be compatible with MetaHMA.

The R file "factor_functions.R" is a function in the R package "cate", where the source code is available at https://github.com/cran/cate/tree/master/R.

References

Sun, Y., Ma, L., & Xia, Y. (2024). A Decorrelating and Debiasing Approach to Simultaneous Inference for High-Dimensional Confounded Models. 
Journal of the American Statistical Association, 119(548), 2857–2868. https://doi.org/10.1080/01621459.2023.2283938

Vazquez, O., & Nan, B. (2025). Debiased lasso after sample splitting for estimation and inference in high‐dimensional generalized linear models. 
Canadian Journal of Statistics, 53(1), e11827.

