p <- 100
q <- 100
rho <- 0.5
K <- 3
R_p <- 0.05
n <- 100

# active in both X model and Y model
R_XY <- 0.05

# active only in X model but not in Y model
R_X <- 0.25

# active only in Y model but not in X model
R_Y <- 0.02

R_active_rows_max <- 0.25

# Generate X

E_X <- matrix(0,nrow=p,ncol=p)
for(i in 1:(p-1)){
  for(j in (i+1):p){
    E_X[i,j] <- (rho)**abs(i-j)
  }
}
E_X <- E_X + t(E_X) + diag(p)

U_X <- MASS::mvrnorm(n = n, mu = rep(0, p), Sigma = E_X)

Fa <- matrix(rnorm(n*K, 0, 1), nrow = n, byrow = T)

B_X <- matrix(runif(n = K*p, 0, 1), nrow = K, byrow = T)

X <- Fa %*% B_X + U_X

# Generate M

A <- matrix(0, nrow = p, ncol = q)
active_cols <- sample(1:q, size = ceiling(q*R_XY), replace = FALSE)
active_cols_add <- sample(setdiff(1:q, active_cols), size = ceiling(q*R_X), replace = FALSE)

for (j in union(active_cols, active_cols_add)) {
  num_active_rows <- sample(2: ceiling(p * R_active_rows_max), 1)
  active_rows <- sample(1:p, size = num_active_rows, replace = FALSE)
  A[active_rows, j] <- rnorm(length(active_rows), mean = 1, sd = sqrt(0)) * 
    sample(c(1, -1), length(active_rows), replace = TRUE)
}

B_M <- matrix(runif(K*q, 0, 1), nrow = K, byrow = T)

E <- matrix(rnorm(n*q, 0, 1), nrow = n, byrow = T)

M <- X %*% A + Fa %*% B_M + E

# Generate Y
beta_active_add <- sample(setdiff(1:q, union(active_cols, active_cols_add)),
                          size = ceiling(q*R_Y), replace = FALSE)
beta_active_indx <- union(active_cols, beta_active_add)
beta <- rep(0, q)
beta[beta_active_indx] <- rnorm(length(beta_active_indx), mean = 1, sd = sqrt(0)) * 
  sample(c(1, -1), length(beta_active_indx), replace = TRUE)

gamma_active_indx <- sample(1:p,ceiling(p*R_p))
gamma  <- rep(0,p)
for(i in gamma_active_indx){
  gamma[i] <- rnorm(1, 1, sqrt(0)) * sample(c(1,-1),1)
}

B_Y <- runif(K, 0, 1)
Fy <- X %*% gamma + M %*% beta + Fa %*% as.vector(B_Y)

prob <- 1/(1+exp(-Fy))

Y <- rbinom(n, size = 1, prob = prob)

ground_truth <- which(A %*% diag(beta) != 0, arr.ind = T)

# Try MetaHMA
MetaHMA(X, Y, M, 
        ncores = 11, 
        B = 500, 
        sig_cut = 0.1, 
        ratio_minscreen = 0.1, 
        p_adjust_method = "BH",
        target = NULL)