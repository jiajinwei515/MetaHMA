# Install required R packages (if not already installed):
# - `paran`: Provides parallel analysis for determining the number of factors in PCA/FA.
# - `doParallel`: Enables parallel computing to improve computational efficiency.
# - `hdi`: Implements methods for high-dimensional inference, e.g., de-sparsified lasso.
# - `glmnet`: Fitting penalized regression models, e.g., lasso and elastic net. 
# - `MASS`

# library(paran)
# library(doParallel)
# library(hdi)
# library(glmnet)
# library(MASS)

# Please note that you need to load all the R functions in helper.

# Parameters and data description:
# Input data:
# X: A matrix of size n * p, representing the exposure variables.
# M: A matrix of size n * q, representing the mediator variables.
# Y: A binary outcome vector of length n (values should be 0 or 1).
#
# Additional parameters:
# ncores: number of cores to be used for parallel computing. 
# Set this to the desired number of CPU cores to improve computational efficiency.
# B: number of iterative sample splits used in the debiased lasso with sample splitting.
# sig_cut: the nominal significance level (FDR threshold).
# ratio_minscreen: the ratio of selected candidate indirect effects in the min-screening procedure.
# p_adjust_method: for more details, please refer to the R documentation: help(p.adjust).
#
# target:
# A vector specifying the target mediators to be included in the debiased lasso with sample splitting. 
# The user can specify the indices of mediators of interest. 
# If set to NULL (default), the input of the debiased lasso with sample splitting 
# will be yielded by the de-sparsified lasso.

MetaHMA <- function(X, Y, M, ncores = 11,
                    B = 1000,
                    sig_cut = 0.1,
                    ratio_minscreen = 0.1,
                    p_adjust_method = "BH",
                    target = NULL)
{
  n <- nrow(X)
  p <- ncol(X)
  q <- ncol(M)
  
  colnames(X) <- paste0("X", 1:ncol(X))
  colnames(M) <- paste0("M", 1:ncol(M))
  
  X <- scale(as.matrix(X), center = TRUE, scale = FALSE)
  M <- scale(as.matrix(M), center = TRUE, scale = FALSE)
  
  # Factor analysis
  message("1: Factor analysis \n")
  
  data_factor <- cbind(X, M)
  
  PA <- paran(data_factor, iterations = 500, centile = 0)
  K_hat <- PA$Retained
  fa <- factor.analysis(data_factor, K_hat, method = "ml")  
  Wupdate.t <- fa$Gamma
  Sgm_inv <- solve(diag(fa$Sigma))
  orthg <- t(Wupdate.t) %*% Sgm_inv %*% Wupdate.t/p
  V <- eigen(orthg)$vectors
  Wupdate <- t(Wupdate.t %*% V)
  Uupdate <- t(solve(Wupdate %*% Sgm_inv %*% t(Wupdate)) %*% Wupdate %*% Sgm_inv %*% t(data_factor)) 
  
  Fhat <- Uupdate
  colnames(Fhat) <- paste0("F", 1:ncol(Fhat))
  Uhat <- data_factor - Uupdate %*% Wupdate
  
  message("1: Factor analysis - done.\n")
  
  data_Y <- as.matrix(cbind(Uhat, Fhat))
  
  # Y model
  message("2: Y model \n")
  
  # initial estimation by lasso
  
  penalty_factor = c(rep(1, p+q), rep(0, ncol(data_Y)-(p+q)))
  
  if(is.null(target) == TRUE)
  {
    cv_fit <- cv.glmnet(data_Y, Y, family = "binomial", alpha = 1, 
                        penalty.factor = penalty_factor,
                        nfolds = 10,
                        standardize = FALSE)
    best_lambda <- cv_fit$lambda.min
    beta_lasso <- as.vector(coef(cv_fit, s = best_lambda))[-1]
    names(beta_lasso) <- colnames(data_Y)
    
    dsparlasso_fit <- lasso.proj(data_Y, Y, family = "binomial", standardize = FALSE,
                                 parallel = TRUE, ncores = ncores,
                                 betainit = beta_lasso, sigma = 1)
    
    beta_pval_dsparlasso <- dsparlasso_fit$pval[(p+1):(p+q)]
    names(beta_pval_dsparlasso) <- colnames(M)
    beta_scr <- beta_pval_dsparlasso[beta_pval_dsparlasso < 2*sig_cut]
    beta_scr_indx <- which(colnames(Uhat) %in% names(beta_scr))
    names(beta_scr_indx) <- colnames(Uhat)[beta_scr_indx]
  }
  
  # debiased lasso with sample splitting
  
  # lasso selection for all bootstrap samples
  cl <- makeCluster(ncores)
  registerDoParallel(cl)
  selection_results <- foreach(i = 1:B, 
                               .packages = c("glmnet"),
                               .export = c("multisplit_selection")) %dopar% {
                                 multisplit_selection(X = data_Y, 
                                                      y = Y, 
                                                      penalty.factor = penalty_factor)
                               }
  stopCluster(cl)
  registerDoSEQ()

  # processing the selection results for next steps
  all_selection_results <- select_results_process(selection_results)
  
  # target covariates selected by screening
  if(is.null(target) == TRUE)
  {
    target_covariate_indx <- beta_scr_indx
  }else
  {
    target_covariate_indx <- p + target
    names(target_covariate_indx) <- colnames(M)[target]
  }
  
  cl <- makeCluster(ncores)
  registerDoParallel(cl)
  multisplit_coefficients <- foreach(b = 1:B, .combine = rbind, .packages = c("glmnet")) %dopar% {
    
    subsample_model_indx <- which(all_selection_results$selected_covariates_indx[b, -1])
    subsample_model_indx <- unique(c(target_covariate_indx, subsample_model_indx))
    subsample_data_indx <- all_selection_results$subsample_indx[b, ]
    
    penalty_factor_1 = rep(1, length(subsample_model_indx))
    penalty_factor_1[which(subsample_model_indx %in% (p+q+1):(p+q+K_hat))] <- 0
    subsample_lasso_fit <- try(cv.glmnet(data_Y[!subsample_data_indx, subsample_model_indx],
                                     Y[!subsample_data_indx],
                                     family = "binomial", 
                                     penalty.factor = penalty_factor_1,
                                     standardize = FALSE), silent = TRUE)
    if(inherits(subsample_lasso_fit, "try-error"))
    {
      output <- NA
    }else{
      subsample_lasso_coef <- as.vector(coef(subsample_lasso_fit, s = "lambda.min"))
      # debiased
      glm_fit <- coef(glm(
        Y[!subsample_data_indx] ~ data_Y[!subsample_data_indx, subsample_model_indx],
        family = "binomial",
        control = list(maxit = 1),
        start = subsample_lasso_coef
      ))
      output <- glm_fit[-1][1:length(target_covariate_indx)]
    }
    output
  }
  
  stopCluster(cl)
  registerDoSEQ()
  
  colnames(multisplit_coefficients) <- names(target_covariate_indx)
  
  indx_na <- which(apply(multisplit_coefficients, 1, function(row) all(is.na(row))))
  B_1 <- B
  multisplit_coefficients_1 <- multisplit_coefficients
  all_selection_results_1 <- all_selection_results
  if(length(indx_na)>0)
  { 
    multisplit_coefficients_1 <- as.matrix(multisplit_coefficients[-indx_na, ])
    colnames(multisplit_coefficients_1) <- names(target_covariate_indx)
    all_selection_results_1$subsample_indx <- all_selection_results$subsample_indx[-indx_na, ]
    B_1 <- B - length(indx_na)
  }
  
  beta_est <- beta_estimation(multisplit_coefficients_1, all_selection_results_1, n, B_1)
  # beta_est <- beta_est[beta_est$pvalue < 2*sig_cut, ]
  I <- which(colnames(M) %in% beta_est$id)
  
  beta_pval <- beta_coef <- rep(NA, q)
  beta_pval[I] <- beta_est$pvalue
  beta_coef[I] <- beta_est$coef_est

  message("\n2: Y model - done. \n")
  
  # M model
  
  cl <- makeCluster(ncores)
  registerDoParallel(cl)
  
  A_sel <- foreach(j = 1:q, .combine = cbind, 
                   .export = c("cov_function", "multiple_testing", "nodewise_lasso", 
                                         "test_stat_decorrelation"),
                   .packages = c("glmnet", "RSpectra", "Matrix", "MASS")) %dopar% {
    if(j %in% I)
    {
      decorr_j = test_stat_decorrelation(M[, j], X, 
                                         initial = "decorrelation", centralize = TRUE,
                                         q = NULL, nlambda = 100, rho = 0.5, kmax = 20)
      pval_j <- 2*(1-pnorm(abs(decorr_j$statistics)))
      A_pval_j <- rep(NA, p)
      A_pval_j <- pval_j
      A_coef_j <- rep(NA, p)
      A_coef_j <- decorr_j$beta_debias
      A_sel_j <- c(A_pval_j, A_coef_j)
    }else{
      A_sel_j <- rep(NA, 2*p)
    }
    return(A_sel_j)
  }
  
  stopCluster(cl)
  
  A_pval <- A_sel[1:p, ]
  A_coef <- A_sel[-(1:p), ]
  
  message("3: M model - done. \n")
  
  # Min-screening
  IE_pval_scr  <- matrix(NA, p, q)
  row.names(IE_pval_scr) <-  colnames(X)
  colnames(IE_pval_scr) <- colnames(M)

  for(i in 1:p){
    for(j in I){
      IE_pval_scr[i,j] <- min(A_pval[i,j], beta_pval[j])
    }
  }

  K <- ratio_minscreen*p*length(I)
  non_na_values <- IE_pval_scr[!is.na(IE_pval_scr)]
  top_K_values <- sort(non_na_values, decreasing = FALSE)[1:K]
  IE_pval_scr_filtered <- IE_pval_scr
  IE_pval_scr_filtered[!(IE_pval_scr %in% top_K_values)] <- NA

  # JST + FDR
  IE <- matrix(NA, p, q)
  row.names(IE) <-  colnames(X)
  colnames(IE) <- colnames(M)
  
  for(i in 1:p){
    for(j in I){
      if(is.na(IE_pval_scr_filtered[i,j]) == FALSE)
      {
        IE[i,j] <- A_coef[i,j] * beta_coef[j]
      }
    }
  }
  
  IE_pval  <- matrix(NA, p, q)
  row.names(IE_pval) <-  colnames(X)
  colnames(IE_pval) <- colnames(M)
  
  for(i in 1:p){
    for(j in I){
      if(is.na(IE_pval_scr_filtered[i,j]) == FALSE)
      {
        IE_pval[i,j] <- max(A_pval[i,j], beta_pval[j])
      }
    }
  }
  
  pval_indices <- which(is.na(IE_pval) == F, arr.ind = TRUE)
  pval_vec <- as.vector(IE_pval[is.na(IE_pval) == F])
  pval_adjusted_vec <- p.adjust(pval_vec, method = p_adjust_method)
  IE_pval_adjusted <- matrix(NA, nrow = p, ncol = q)
  IE_pval_adjusted[pval_indices] <- pval_adjusted_vec
  row.names(IE_pval_adjusted) <-   row.names(IE_pval)
  colnames(IE_pval_adjusted) <- colnames(IE_pval)
  
  final_results <- data.frame(which(IE_pval_adjusted<sig_cut, arr.ind = T), 
                              IE[which(IE_pval_adjusted<sig_cut)],
                              A_coef[which(IE_pval_adjusted<sig_cut)],
                              beta_coef[which(IE_pval_adjusted<sig_cut, arr.ind = T)[,2]],
                              IE_pval_adjusted[which(IE_pval_adjusted<sig_cut)])
  colnames(final_results) <- c("X", "M", "IE", "A", "beta", "adjusted_pval")
  rownames(final_results) <- 1 : nrow(final_results)
  
  message("COMPLETE")
  return(final_results)
  
}
