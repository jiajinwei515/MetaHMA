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

## You may run MetaHMA by the following steps.

#### Load the required R functions in helper and MetaHMA
#### Set this to the folder where you saved the files: decorr.R, factor_functions.R, MetaHMA.R, ssdblasso.R

helper_dir <- "/YOUR/LOCAL/PATH/TO/helper"  # <-- change this

#### Load the required functions
source(file.path(helper_dir, "decorr.R"))

source(file.path(helper_dir, "factor_functions.R"))

source(file.path(helper_dir, "ssdblasso.R"))

source(file.path(helper_dir, "MetaHMA.R"))

#### Load the required packages
library(paran)

library(doParallel)

library(hdi)

library(glmnet)

library(MASS)

#### Set this to the folder where you saved the CSV files in data

data_dir <- "/YOUR/LOCAL/PATH/TO/data"   # <-- change this

#### Read data (drop first column since it is an index/sample ID)

X <- read.csv(file.path(data_dir, "X_ALR_transformed_viral.csv"),
              header = TRUE, check.names = FALSE)[, -1]

M <- read.csv(file.path(data_dir, "M_ALR_transformed_bacterial.csv"),
              header = TRUE, check.names = FALSE)[, -1]

Y <- read.csv(file.path(data_dir, "Y_disease_outcomes.csv"),
              header = FALSE, check.names = FALSE)[, -1]

#### Check the dimension
dim(X) # 94 205

dim(M) # 94 97

length(Y) # 94

#### Run MetaHMA (Please set your "ncores" for parallel computation)
MetaHMA_results <- MetaHMA(X, Y, M, ncores = 4, B = 500, 
                            sig_cut = 0.1, 
                            ratio_minscreen = 0.1, 
                            p_adjust_method = "BH",
                            target = NULL)

