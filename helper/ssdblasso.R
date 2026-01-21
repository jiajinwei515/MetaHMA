multisplit_selection <- function(X, y, penalty.factor){
  n <- nrow(X)
  subsample_indx <- rep(FALSE, n)
  #subsample_indx[sample(1:n, n / 2, replace = FALSE)] <- TRUE
  indx_0 <- which(y == 0)
  indx_1 <- which(y == 1)
  n_sample_0 <- floor(length(indx_0) / 2)
  n_sample_1 <- floor(length(indx_1) / 2)
  sampled_0 <- sample(indx_0, n_sample_0, replace = FALSE)
  sampled_1 <- sample(indx_1, n_sample_1, replace = FALSE)
  subsample_indx[c(sampled_0, sampled_1)] <- TRUE
  selection_lasso_fit <- try(cv.glmnet(X[subsample_indx, ], y[subsample_indx],
                                   family = "binomial", 
                                   penalty.factor = penalty.factor,
                                   standardize = FALSE), silent = TRUE)
  while(inherits(selection_lasso_fit, "try-error"))
  {
    subsample_indx <- rep(FALSE, n)
    indx_0 <- which(y == 0)
    indx_1 <- which(y == 1)
    n_sample_0 <- floor(length(indx_0) / 2)
    n_sample_1 <- floor(length(indx_1) / 2)
    sampled_0 <- sample(indx_0, n_sample_0, replace = FALSE)
    sampled_1 <- sample(indx_1, n_sample_1, replace = FALSE)
    subsample_indx[c(sampled_0, sampled_1)] <- TRUE
    selection_lasso_fit <- try(cv.glmnet(X[subsample_indx, ], y[subsample_indx],
                                         family = "binomial", 
                                         penalty.factor = penalty.factor,
                                         standardize = FALSE), silent = TRUE)
  }
  lasso_coef <- coef(selection_lasso_fit, s = "lambda.min")
  coef_names <- rownames(lasso_coef)
  lasso_coef <- as.data.frame(t(as.vector(lasso_coef)))
  colnames(lasso_coef) <- coef_names
  results <- list(subsample_indx = as.data.frame(t(subsample_indx)),
                  beta_lasso = lasso_coef)
  return(results)
}

# n_cores <- 10
# cl <- makeCluster(n_cores)
# registerDoParallel(cl)
# B <- 1000
# penalty.factor = c(rep(1, p+q), rep(0, ncol(data_Y)-(p+q)))
# selection_results <- foreach(i = seq_len(B), .packages = c("glmnet")) %dopar% {
#   multisplit_selection(X = data_Y, y = as.vector(Y), penalty.factor = penalty.factor)
# }
# stopCluster(cl)

select_results_process <- function(selection_results)
{
  all_selection_results <- selection_results[[1]]
  for (j in seq_along(all_selection_results)){
    if ("data.frame" %in% class(all_selection_results[[j]])){
      tmp <- as.data.frame(matrix(nrow = length(selection_results) - 1, 
                                  ncol = ncol(all_selection_results[[j]])))
      colnames(tmp) <- colnames(all_selection_results[[j]])
      all_selection_results[[j]] <- rbind(all_selection_results[[j]], tmp)
      
    }
  }
  
  for (i in seq_along(selection_results)){
    if (i > 1){
      for (j in seq_along(all_selection_results)){
        all_selection_results[[j]][i, ] <- c(selection_results[[i]][[j]])
      }
    }
  }
  
  all_selection_results$selected_covariates_indx <- all_selection_results$beta_lasso != 0
  return(all_selection_results)
}



# target_covariate_indx <- which(beta_lasso_full[(p+1):(p+q)] != 0) + p

# multisplit_coefficients <- as.data.frame(matrix(0, nrow = B, ncol = length(target_covariate_indx)))
# for (b in seq_len(B)){
#   subsample_model_indx <- which(all_selection_results$selected_covariates_indx[b, -1])
#   subsample_model_indx <- unique(c(target_covariate_indx, subsample_model_indx))
#   subsample_data_indx <- all_selection_results$subsample_indx[b, ]
# 
#   # Fit debiased lasso
#   penalty.factor = rep(1, length(subsample_model_indx))
#   penalty.factor[which(subsample_model_indx %in% (p+q+1):ncol(data_Y))] <- 0
#   subsample_lasso_fit <- cv.glmnet(data_Y[!subsample_data_indx, subsample_model_indx],
#                                    Y[!subsample_data_indx],
#                                    family = "binomial", 
#                                    penalty.factor = penalty.factor,
#                                    standardize = FALSE)
#   subsample_lasso_coef <- as.vector(coef(subsample_lasso_fit, s = "lambda.min"))
#   multisplit_coefficients[b, ] <- coef(glm(Y[!subsample_data_indx] ~
#                                              data_Y[!subsample_data_indx, subsample_model_indx],
#                                            family = "binomial", control = list(maxit = 1),
#                                            start = subsample_lasso_coef))[-1][1:length(target_covariate_indx)]
# }


# # 定义并行计算的核心数
# num_cores <- 11
# cl <- makeCluster(num_cores)
# registerDoParallel(cl)
# 
# # 使用 foreach 进行并行计算
# multisplit_coefficients <- foreach(b = seq_len(B), .combine = rbind, .packages = c("glmnet")) %dopar% {
#   # 提取当前分割的索引
#   subsample_model_indx <- which(all_selection_results$selected_covariates_indx[b, -1])
#   subsample_model_indx <- unique(c(target_covariate_indx, subsample_model_indx))
#   subsample_data_indx <- all_selection_results$subsample_indx[b, ]
# 
#   penalty.factor = rep(1, length(subsample_model_indx))
#   penalty.factor[which(subsample_model_indx %in% (p+q+1):ncol(data_Y))] <- 0
#   subsample_lasso_fit <- cv.glmnet(data_Y[!subsample_data_indx, subsample_model_indx],
#                                    Y[!subsample_data_indx],
#                                    family = "binomial", 
#                                    penalty.factor = penalty.factor,
#                                    standardize = FALSE)
#   subsample_lasso_coef <- as.vector(coef(subsample_lasso_fit, s = "lambda.min"))
# 
#   # 使用 glm 进行校正
#   glm_fit <- coef(glm(
#     Y[!subsample_data_indx] ~ data_Y[!subsample_data_indx, subsample_model_indx],
#     family = "binomial",
#     control = list(maxit = 1),
#     start = subsample_lasso_coef
#   ))
# 
#   # 返回当前分割的系数
#   glm_fit[-1][1:length(target_covariate_indx)]
# }
# 
# # 停止并行计算
# stopCluster(cl)
# 
# colnames(multisplit_coefficients) <- names(target_covariate_indx)

beta_estimation <- function(multisplit_coefficients, all_selection_results, n, B)
{
  coef_ests <- colMeans(multisplit_coefficients)
  coef_stderrs <- vector(length = length(coef_ests))
  centered_sample_indicators <- scale(!all_selection_results$subsample_indx, center = TRUE, scale = FALSE)
  n2 <- ceiling(n/2)
  for (j in seq_along(coef_ests)){
    cov_hat_j <- t(multisplit_coefficients[, j] - coef_ests[j]) %*% centered_sample_indicators / B 
    V_hat_j <- (n * (n - 1)) * sum(cov_hat_j^2) / (n - n2)^2
    coef_stderrs[j] <- sqrt(max(V_hat_j - ((n / B^2) * (n2 / (n - n2)) * 
                                             sum((multisplit_coefficients[, j] - coef_ests[j])^2)),
                                0))
  }
  
  results_table <- data.frame(id = names(coef_ests),
                              coef_est = coef_ests,
                              #stderr = coef_stderrs,
                              #ci_lower = coef_ests - (qnorm(0.975) * coef_stderrs),
                              #ci_upper = coef_ests + (qnorm(0.975) * coef_stderrs),
                              pvalue = 2 * pnorm(abs(coef_ests / coef_stderrs), lower.tail = FALSE))
  
  return(results_table)
}


# target_covariate_indx[which(results_table$pvalue < 0.05)]
# 
# p.adjust(results_table$pvalue, method = "BH")
# 
# unique(ground_truth[,2])
