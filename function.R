
### Functions to generate the simulation data ###

AR_cov <- function(p, acorr){
  cov_mat <- diag(rep(1, p))
  cor_vec <- acorr^c(1:p)
  for (i in 1:(p - 1)){
    cov_mat[i, (i + 1):p] <- cor_vec[1: (p - i)]
    cov_mat[(i + 1):p, i] <- cor_vec[1: (p - i)]
  }
  return(cov_mat)
}

Gen_X_each <- function(n, p, r = 0.4, magn = 0.5, sigma = 0.3, 
                       sign = rep(1, 5)){
  Sigma_X_remain <- AR_cov(p - 50, r)
  Sigma_X_control <- AR_cov(45, r)
  
  X_remain <- mvrnorm(n, rep(0, p - 50), Sigma_X_remain)
  X_control <- mvrnorm(n, rep(0, 45), Sigma_X_control)
  X_pv <- c()
  for (t in 1:5) {
    gamma_coef <- rep(r, 45) * (2 * rbinom(45, 1, 0.5) - 1)
    X_t <- X_control %*% gamma_coef + rnorm(n, 0, 1)
    X_pv <- cbind(X_pv, X_t)
  }
  X_m <- cbind(X_pv, X_control, X_remain)
  
  X_inter <- 0
  for (s in 1:4) {
    X_inter <- X_inter + X_m[,s] * X_m[,s + 1]
  }
  X_main <- (X_m[,1:5] + 0.2 * (X_m[,1:5])^3) %*% (magn + sigma * sign)
  p_vec <- 1 / (1 + exp(- X_main - 0.4 * magn * X_inter))
  Y_m <- rbinom(n, 1, p_vec)
  
  return(list(X = X_m, Y = Y_m))
}

Generate_data <- function(n, p, M, magn_mu, magn_alpha, 
                          r_lst = rep(0.3, M), model = 'correct_sparse'){
  X_lst <- vector('list', M)
  Y_lst <- vector('list', M)
  if (model %in% c('correct_sparse', 'correct_dense')){
    
    # Generate beta
    
    if (model == 'correct_sparse'){
      s <- 6
      sign_lst <- c(1, -1, 1, -1, 1, -1)
      center_support <- c(1:6)
      beta_center <- rep(0, p)
      beta_center[center_support] <- sign_lst * rep(magn_mu, s)
      beta_center <- c(0, beta_center)
      
      center_support <- center_support + 1
      noise_support <- c(4:9)
      
      neg_indx <- sample(c(1:M), as.integer(M / 2))
      beta_delta_lst <- vector('list', M)
      beta_true_lst <- vector('list', M)
      
      sign_lst <-  c(1, 1, 1, -1, -1, -1)
      
      for (m in 1:M) {
        beta_delta_lst[[m]] <- rep(0, (p + 1))
        for (t in 1:length(noise_support)){
          sign_t <- sign_lst[t]
          if (m %in% neg_indx){
            beta_delta_lst[[m]][noise_support[t]] <- - sign_t * magn_alpha
          }else{
            beta_delta_lst[[m]][noise_support[t]] <- sign_t * magn_alpha
          }
        }
        
        beta_true_lst[[m]] <- beta_center + beta_delta_lst[[m]]
      }
      
    }
    if (model == 'correct_dense'){
      s <- 18
      sign_lst <- rep(c(1, -1), s / 2)
      
      #center_support <- sample(1:p, s1)
      center_support <- c(1:s)
      beta_center <- rep(0, p)
      beta_center[center_support] <- sign_lst * rep(magn_mu, s)
      beta_center <- c(0, beta_center)
      noise_support <- c((s / 3 + 1):(4 * s / 3)) + 1
      beta_delta_lst <- vector('list', M)
      beta_true_lst <- vector('list', M)
      
      for (m in 1:M) {
        beta_delta_lst[[m]] <- rep(0, (p + 1))
        if (m == 1){
          alpha_loc <- magn_alpha * c(rep(0, 3), rep(1.8, 3), rep(0.7, 3), 
                                      rep(-1.3, 3), rep(-0.85, 3), rep(1.15, 3))
          beta_delta_lst[[m]][noise_support] <- alpha_loc
        }
        if (m == 2){
          alpha_loc <- magn_alpha * c(rep(0, 3), rep(-1.8, 3), rep(1.3, 3), 
                                      rep(-0.7, 3), rep(-1.15, 3), rep(0.85, 3))
          beta_delta_lst[[m]][noise_support] <- alpha_loc
        }
        if (m == 3){
          alpha_loc <- magn_alpha * c(rep(-1.8, 3), rep(0, 3), rep(-0.85, 3), 
                                      rep(1.15, 3), rep(0.7, 3), rep(-1.3, 3))
          beta_delta_lst[[m]][noise_support] <- alpha_loc
        }
        if (m == 4){
          alpha_loc <- magn_alpha * c(rep(1.8, 3), rep(0, 3), rep(-1.15, 3), 
                                      rep(0.85, 3), rep(1.3, 3), rep(-0.7, 3))
          beta_delta_lst[[m]][noise_support] <- alpha_loc
        }
        if (m == 5){
          alpha_loc <- magn_alpha * c(rep(0, 3), rep(1.5, 3), rep(0.5, 3), 
                                      rep(-1.1, 3), rep(0.8, 3), rep(-1, 3))
          beta_delta_lst[[m]][noise_support] <- alpha_loc
        }
        if (m == 6){
          alpha_loc <- magn_alpha * c(rep(0, 3), rep(-1.5, 3), rep(1.2, 3), rep(-0.6, 3), rep(0.9, 3), rep(-0.7, 3))
          beta_delta_lst[[m]][noise_support] <- alpha_loc
        }
        if (m == 7){
          alpha_loc <- magn_alpha * c(rep(-1.5, 3), rep(0, 3), rep(-0.8, 3), rep(1, 3), rep(-0.5, 3), rep(1.1, 3))
          beta_delta_lst[[m]][noise_support] <- alpha_loc
        }
        if (m == 8){
          alpha_loc <- magn_alpha * c(rep(1.5, 3), rep(0, 3), rep(-0.9, 3), rep(0.7, 3), rep(-1.2, 3), rep(0.6, 3))
          beta_delta_lst[[m]][noise_support] <- alpha_loc
        }
        
        beta_true_lst[[m]] <- beta_center + beta_delta_lst[[m]]
      }
    }
    
    # Generate X and Y.
    
    prop <- 15 / p
    for (m in 1:M){
      r <- r_lst[m]
      gamma_coef <- c()
      for (t in 1:s) {
        gamma_coef <- cbind(gamma_coef, r * rbinom(p - s, 1, prop))
      }
      Sigma_remain <- AR_cov(p - s, r)
      Sigma_X <- rbind(cbind(diag(rep(1, s)) + t(gamma_coef) %*% Sigma_remain %*% gamma_coef,
                             t(Sigma_remain %*% gamma_coef)), 
                       cbind(Sigma_remain %*% gamma_coef, Sigma_remain))
      X <- mvrnorm(n, rep(0, p), Sigma_X)
      odds <- exp(cbind(rep(1, n), X) %*% beta_true_lst[[m]])
      p_vec <- odds / (1 + odds)
      Y <- rbinom(n, 1, p_vec)
      X_lst[[m]] <- X
      Y_lst[[m]] <- Y
    }
  }
  
  if (model == 'wrong'){
    
    beta_delta_lst <- vector('list', M)
    beta_true_lst <- vector('list', M)
    sign_lst <- c(1, -1, 1, -1, 1)
    for (m in 1:M){
      
      # Generate X and Y.
      
      r <- r_lst[m]
      Gen_data_m <- Gen_X_each(n, p, r = r, magn = magn_mu, sigma = magn_alpha, 
                               sign = sign_lst * (-1)^m)
      
      X_lst[[m]] <- Gen_data_m$X
      Y_lst[[m]] <- Gen_data_m$Y
      
      # Estimate the true beta with data of large sample size
      
      Gen_data_large <- Gen_X_each(100000, p, r = r, magn = magn_mu, sigma = magn_alpha, 
                                   sign = sign_lst * (-1)^m)
      X_large <- Gen_data_large$X
      Y_large <- Gen_data_large$Y
      glm.true <- glm.fit(cbind(1, X_large[,1:50]), Y_large, family = binomial(),
                          intercept = F)
    
      beta_true_lst[[m]] <- c(glm.true$coefficients, rep(0, p - 50))
      
    }
    
    beta_center <- 0
    for (m in 1:M){
      beta_center <- beta_center + beta_true_lst[[m]] 
    }
    beta_center <- beta_center / M
    for (m in 1:M){
      beta_delta_lst[[m]] <- beta_true_lst[[m]] - beta_center
    }
    
  }
  return(list(X = X_lst, Y = Y_lst, beta_true = beta_true_lst, 
              beta_center_true = beta_center, beta_delta_true = beta_delta_lst))
}
  

### Function for local fit ###

Local_fit <- function(Y_lst, X_lst, family = 'binomial', lambda_lst = NULL){
  M <- length(X_lst)
  I_lst <- vector('list', M)
  U_lst <- vector('list', M)
  beta_lst <- c()
  
  for (m in 1:M){
    Y <- Y_lst[[m]]
    X <- X_lst[[m]]
    
    if (length(lambda_lst) == 1){
      lambda.cv <- lambda_lst[1]
    }else{
      cv.result <- cv.glmnet(X, Y, family = family, lambda = lambda_lst)
      lambda.cv <- cv.result$lambda.min
    }
    model <- glmnet(X, Y, family = family, lambda = lambda.cv)
    beta_fit <- c(as.vector(model$a0), as.vector(model$beta))
    
    n <- length(Y)
    X_all <- cbind(rep(1, n), X)
    pi_vec <- as.vector(1 / (1 + exp(- X_all %*% beta_fit)))
    grad <- t(X_all) %*% (Y - pi_vec)
    I_mat <- t(X_all) %*% diag(pi_vec * (1- pi_vec)) %*% X_all
    U <- I_mat %*% beta_fit + grad
    I_lst[[m]] <- I_mat
    U_lst[[m]] <- U
    beta_lst <- cbind(beta_lst, beta_fit)
  }
  return(list(I = I_lst, U = U_lst, beta = beta_lst))
}


### Function for SHIR ###


SHIR_fit <- function(H_lst, d_lst, n_lst, lambda_lst, lambda_g_lst, tune = 'BIC'){
  
  options(warn = -1)
  M <- length(d_lst)
  X_lst <- vector('list', M)
  Y_lst <- vector('list', M)
  n <- sum(n_lst)
  p <- length(H_lst[[1]][1, ]) - 1
  initial <- rep(0, (M + 1) * (p + 1))
  
  for (i in 1:M){
    H <- H_lst[[i]]
    d <- d_lst[[i]]
    mat_all <- cbind(H, d)
    mat_all <- rbind(mat_all, t(c(d, max(mat_all) + n)))
    svd_result <- svd(mat_all)
    s_value <- svd_result$d
    s_mat <- diag(sqrt(s_value))[ ,1:(min(p + 1, n_lst[[i]]) + 1)]
    data_all <- svd_result$u %*% s_mat
    X <- t(data_all[-length(data_all[ ,1]),])
    Y <- data_all[length(data_all[ ,1]),]
    X_lst[[i]] <- X
    Y_lst[[i]] <- Y
  }
  
  X_all <- c()
  Y_all <- c()
  for (i in 1:M){
    X_all <- rbind(X_all, X_lst[[i]])
    Y_all <- c(Y_all, Y_lst[[i]])
  }
  
  X_alpha <- X_lst[[1]]
  for (i in 2:M){
    X_alpha <- bdiag(X_alpha, X_lst[[i]])
  }
  X_all <- cbind(X_all, X_alpha)
  
  enlarge_mat <- matrix(0, p + 1, p + 1)
  for (i in 1:M) {
    enlarge_mat <- cbind(enlarge_mat, diag(rep(1, p + 1)))
  }
  
  indx_lst <- c(c(NA, 1:p), rep(c(NA, (p + 1):(2 * p)), M))
  
  # tune with GIC
  min.lambda <- NULL
  min.GIC <- Inf
  min.coef <- NULL
  min.beta <- NULL
  
  total_num <- length(lambda_g_lst) * length(lambda_lst)
  z <- 0
  for (lambda in lambda_lst){
    for (lambda_g in lambda_g_lst){
      
      lambda_g <- lambda_g * lambda
      penscale <- function(x){
        if (x == 1){
          return(1)
        }else{
          return(lambda_g / lambda)
        }
      }
      
      Y_ <- 0
      rho_ <- n
      enlarge_Y <- rep(- Y_ / sqrt(rho_ / 2) / 2, p + 1)
      X_enlarge <- rbind(X_all, sqrt(rho_ / 2) * enlarge_mat)
      Y_enlarge <- c(Y_all, enlarge_Y)
      
      eta <- 2
      for (iter in 1:10) {
        invisible(utils::capture.output(
        fit_result <- grplasso(x = as.matrix(X_enlarge), y = Y_enlarge, 
                               index = indx_lst, standardize = F, 
                               center = FALSE, lambda = lambda, 
                               penscale = penscale, model = LinReg(), 
                               coef.init = initial)))
        fit_coef <- fit_result$coefficients
        initial <- fit_coef
        
        alpha_fit <- c()
        for (i in 2:(M + 1)){
          alpha_fit <- cbind(alpha_fit, fit_coef[((i - 1) * p + i):(i * p + i)])
        }
        
        Y_ <- Y_ + rho_ * rowMeans(alpha_fit)
        rho_ <- rho_ * eta
        X_enlarge <- rbind(X_all, sqrt(rho_ / 2) * enlarge_mat)
        Y_enlarge <- c(Y_all, - Y_ / sqrt(rho_ / 2) / 2 * rep(1, p + 1))
       
      }
      
      result_GIC <- Cal_GIC_pool(d_lst, H_lst, X_enlarge, 
                                 fit_coef, n_lst, lambda_g / 2, type = tune)
      GIC <- result_GIC$GIC
      if (GIC < min.GIC){
        min.lambda <- c(lambda, lambda_g / lambda)
        min.GIC <- GIC
        min.beta <- result_GIC$beta
        min.coef <- fit_coef
      }
      z <- z + 1
      print(paste0('Finishing: ', 100 * round(z / total_num, 4), '%'))
    }
  }
  return(list(min.lambda = min.lambda, min.beta = min.beta,
              min.coef = min.coef, min.GIC = min.GIC))
}


# Function to tune:

Cal_GIC_pool <- function(U_lst, I_lst, X_enlarge, fit_coef, length_lst, lambda_g, 
                         type = 'BIC'){
  M <- length(U_lst)
  p <- length(I_lst[[1]][1, ]) - 1
  mu <- fit_coef[1:(p + 1)]
  beta_fit <- c()
  alpha_fit <- c()
  
  for (i in 2:(M + 1)){
    beta_fit <- cbind(beta_fit, mu + fit_coef[((i - 1) * p + i):(i * p + i)])
    alpha_fit <- cbind(alpha_fit, fit_coef[((i - 1) * p + i):(i * p + i)])
  }
  
  norm_lst <- sqrt(rowSums(alpha_fit^2))
  S_alpha <- which(norm_lst != 0)
  S_mu <- which(mu != 0)
  S_full <- which(fit_coef != 0)
  H_S <- t(X_enlarge[ ,S_full]) %*% X_enlarge[ ,S_full]
  partial <- diag(0, length(S_full), length(S_full))
  if (length(S_alpha) != 0){
    for (t in 1:length(S_alpha)) {
      j <- S_alpha[t]
      partial_j <- t + length(S_mu) + length(S_alpha) * c(0:(M - 1))
      partial[partial_j, partial_j] <- lambda_g / norm_lst[j] * 
        (diag(1, M, M) - alpha_fit[j,] %*% t(alpha_fit[j,]) / norm_lst[j]^2)
    }
  }
  
  df <- sum(diag(solve(H_S + partial) %*% H_S))
  #print(df)
  
  GIC <- 0
  for (i in 1:M){
    GIC <- GIC + t(beta_fit[,i]) %*% I_lst[[i]] %*% beta_fit[,i] - 2 * t(beta_fit[,i]) %*% U_lst[[i]]
  }
  
  if (type == 'BIC'){
    GIC <- GIC / sum(length_lst) + df * log(sum(length_lst)) / sum(length_lst)
  }
  if (type == 'mBIC'){
    GIC <- GIC / sum(length_lst) + df * log(max(exp(1), log(p))) * log(sum(length_lst)) / sum(length_lst)
  }
  if (type == 'AIC'){
    GIC <- GIC / sum(length_lst) + df * 2 / sum(length_lst)
  }
  if (type == 'RIC'){
    GIC <- GIC / sum(length_lst) + df * log(M * p) / sum(length_lst)
  }
  return(list(GIC = GIC, beta = beta_fit))
}



### Function for debias ###

debias_fit <- function(H_lst, g_lst, beta_fit, n_lst, tau_lst = NULL,
                       nu1_lst = NULL, nu2_lst = NULL, tune = 'BIC'){
  
  
  M <- length(H_lst)
  n <- sum(n_lst) / M
  p <- length(H_lst[[1]][1, ]) - 1
  
  if (is.null(tau_lst)){
    tau_lst <- sqrt(log(p) / n)
  }
  if (is.null(nu1_lst)){
    nu1_lst <- sqrt(log(p) / sum(n_lst))
  }
  if (is.null(nu2_lst)){
    nu2_lst <- sqrt((log(p) + M) / n)
  }
  
  X_lst <- vector('list', M)
  Y_lst <- vector('list', M)
  
  for (i in 1:M){
    H <- H_lst[[i]]
    d <- g_lst[[i]]
    mat_all <- cbind(H, d)
    mat_all <- rbind(mat_all, t(c(d, max(mat_all) + n)))
    
    svd_result <- svd(mat_all)
    s_value <- svd_result$d
    s_mat <- diag(sqrt(s_value))[ ,1:(min(p + 1, n_lst[[i]]) + 1)]
    data_all <- svd_result$u %*% s_mat
    
    X <- t(data_all[-length(data_all[ ,1]),])
    Y <- data_all[length(data_all[ ,1]),]
    X_lst[[i]] <- X
    Y_lst[[i]] <- Y
  }
  
  X_all <- c()
  for (i in 1:M){
    X_all <- rbind(X_all, X_lst[[i]])
  }
  
  X_alpha <- X_lst[[1]]
  for (i in 2:M){
    X_alpha <- bdiag(X_alpha, X_lst[[i]])
  }
  X_all <- cbind(X_all, X_alpha)
  
  enlarge_mat <- matrix(0, p + 1, p + 1)
  for (i in 1:M) {
    enlarge_mat <- cbind(enlarge_mat, 10000 * diag(rep(1, p + 1)))
  }
  X_enlarge <- rbind(X_all, enlarge_mat)
  
  min.lambda <- NULL
  min.GIC <- Inf
  min.beta <- NULL
  min.coef <- NULL
  
  beta_deb_mat <- matrix(0, p + 1, M)
  for (j in 1:(p + 1)) {
    ej <- rep(0, p + 1)
    ej[j] <- 1
    
    for (m in 1:M){
      X_mj <- X_lst[[m]][, -j]
      X_j <- X_lst[[m]][,j]
      
      if (length(tau_lst) >= 2){
        model_fit <- cv.glmnet(X_mj, X_j, lambda = tau_lst * sqrt(n / n_lst[m]), intercept = F,
                               standardize = F)
        tau <- model_fit$lambda.min
        model_fit <- glmnet(X_mj, X_j, lambda = tau, intercept = F, standardize = F)
      }else{
        tau <- tau_lst[1]
        model_fit <- glmnet(X_mj, X_j, lambda = tau, intercept = F, standardize = F)
      }
    
      u_j <- rep(1, p + 1)
      u_j[-j] <- - model_fit$beta
      sigma2 <- mean((X_j - X_mj %*% model_fit$beta)^2) + tau * sqrt(n / n_lst[m]) * sum(abs(model_fit$beta))
      u_m <- as.vector(u_j / sigma2)
      
      beta_deb <- beta_fit[j, m] + 
        t(u_m) %*% (g_lst[[m]] - H_lst[[m]] %*% beta_fit[,m]) / n_lst[m]
      beta_deb_mat[j, m] <- beta_deb
    }
  }
  
  for (nu1 in nu1_lst){
    for (nu2 in nu2_lst) {
      
      mu_all <- rowMeans(beta_deb_mat)
      alpha_all <- beta_deb_mat - mu_all
      alpha_norm <- sqrt(rowSums(alpha_all^2))
      
      supp_mu_pos <- setdiff(which(mu_all > nu1), 1)
      mu_all[supp_mu_pos] <- mu_all[supp_mu_pos] - nu1
      supp_mu_neg <- setdiff(which(mu_all < - nu1), 1)
      mu_all[supp_mu_neg] <- mu_all[supp_mu_neg] + nu1
      mu_all[setdiff(2:(p + 1), union(supp_mu_pos, supp_mu_neg))] <- 0
      
      supp_alpha <- setdiff(which(alpha_norm > nu2), 1)
      alpha_all[setdiff(2:(p + 1), supp_alpha), ] <- 0
      alpha_all[supp_alpha, ] <- alpha_all[supp_alpha, ] * (1 - nu2 / alpha_norm[supp_alpha])
      
      fit_coef <- mu_all
      for (m in 1:M){
        fit_coef <- c(fit_coef, alpha_all[,m])
      }
      fit_coef <- as.matrix(fit_coef)
      result_GIC <- Cal_GIC_pool(g_lst, H_lst, X_enlarge, 
                                 fit_coef, n_lst, nu2, type = tune)
      GIC <- result_GIC$GIC
      if (GIC < min.GIC){
        min.lambda <- c(nu1, nu2)
        min.GIC <- GIC
        min.beta <- result_GIC$beta
        min.coef <- fit_coef
      }
    }
  }
  
  return(list(beta = min.beta, coef = min.coef, lambda = min.lambda,
              GIC = min.GIC))
}




### Function for IPD ###

IPD_fit <- function(X_lst_0, Y_lst, n_lst, lambda_lst, lambda_g_lst, tune = 'BIC'){
  
  M <- length(Y_lst)
  options(warn = -1)
  X_lst <- vector('list', M)
  
  for (i in 1:M){
    X_lst[[i]] <- cbind(rep(1, length(X_lst_0[[i]][,1])), X_lst_0[[i]])
  }
  p <- length(X_lst[[1]][1, ]) - 1
  initial <- rep(0, (M + 1) * (p + 1))
  
  X_all <- c()
  Y_all <- c()
  for (i in 1:M){
    X_all <- rbind(X_all, X_lst[[i]])
    Y_all <- c(Y_all, Y_lst[[i]])
  }
  
  X_alpha <- X_lst[[1]]
  for (i in 2:M){
    X_alpha <- bdiag(X_alpha, X_lst[[i]])
  }
  X_all <- cbind(X_all, X_alpha)
  
  enlarge_mat <- matrix(0, p + 1, p + 1)
  for (i in 1:M) {
    enlarge_mat <- cbind(enlarge_mat, diag(rep(1, p + 1)))
  }
  
  indx_lst <- c(c(NA, 1:p), rep(c(NA, (p + 1):(2 * p)), M))
  
  # tune with GIC
  min.lambda <- NULL
  min.GIC <- Inf
  min.coef <- NULL
  min.beta <- NULL
  
  total_num <- length(lambda_g_lst) * length(lambda_lst)
  z <- 0
  
  for (lambda in lambda_lst){
    for (lambda_g in lambda_g_lst){
      lambda_g <- lambda_g * lambda
      penscale <- function(x){
        if (x == 1){
          return(1)
        }else{
          return(lambda_g / lambda)
        }
      }
      
      ind_sample <- 1:sum(n_lst)
      ind_enlarge <- (sum(n_lst) + 1):(sum(n_lst) + p + 1)
      logReg_aL <- function(){
        grpl.model(invlink = function(eta) c(1/(1 + exp(-eta[ind_sample])), eta[ind_enlarge]), 
                   link = function(mu) c(log(mu[ind_sample]/(1 - mu[ind_sample])), mu[ind_enlarge]),
                   nloglik = function(y, eta, weights, ...) 
                     - sum(weights[ind_sample] *
                             (y[ind_sample] * eta[ind_sample] - log(1 + exp(eta[ind_sample])))) +
                      - sum(weights[ind_enlarge] * 
                              (2 * y[ind_enlarge] * eta[ind_enlarge] - eta[ind_enlarge]^2) / 2),
                   ngradient = function(x, y, mu, weights, ...)  -crossprod(x, weights * (y - mu)),
                   nhessian = function(x, mu, weights, ...) 
                     crossprod(x[ind_sample,], weights[ind_sample] * mu[ind_sample] * (1 - mu[ind_sample]) * x[ind_sample,]) + 
                     crossprod(x[ind_enlarge,],  weights[ind_enlarge] * x[ind_enlarge,]), 
                   check = function(y) {T}, name = "augmented Logistic Regression Model", 
                   comment = "")
      } 
      
      Y_ <- 0
      rho_ <- n
      enlarge_Y <- rep(- Y_ / sqrt(rho_), p + 1)
      X_enlarge <- rbind(X_all, sqrt(rho_) * enlarge_mat)
      Y_enlarge <- c(Y_all, enlarge_Y)
      
      eta <- 2
      for (iter in 1:10) {
        #print(iter)
        invisible(utils::capture.output(
        fit_result <- grplasso(x = as.matrix(X_enlarge), y = Y_enlarge, 
                               index = indx_lst, standardize = F, 
                               center = FALSE, lambda = lambda, 
                               penscale = penscale, model = logReg_aL(), coef.init = initial)))
        fit_coef <- fit_result$coefficients
        initial <- fit_coef
        
        alpha_fit <- c()
        for (i in 2:(M + 1)){
          alpha_fit <- cbind(alpha_fit, fit_coef[((i - 1) * p + i):(i * p + i)])
        }
        
        Y_ <- Y_ + rho_ * rowMeans(alpha_fit)
        rho_ <- rho_ * eta
        X_enlarge <- rbind(X_all, sqrt(rho_) * enlarge_mat)
        Y_enlarge <- c(Y_all, - Y_ / sqrt(rho_)* rep(1, p + 1))
        
      }
      
      fit_coef <- fit_result$coefficients
      initial <- fit_coef
    
      result_GIC <- Cal_GIC_pool_llh(Y_lst, X_lst, fit_coef, n_lst, 
                                     lambda_g, type = tune)
      GIC <- result_GIC$GIC
      if (GIC < min.GIC){
        min.lambda <- c(lambda, lambda_g / lambda)
        min.GIC <- GIC
        min.beta <- result_GIC$beta
        min.coef <- fit_coef
      }
      z <- z + 1
      print(paste0('Finishing: ', 100 * round(z / total_num, 4), '%'))
    }
  }
  return(list(min.lambda = min.lambda, min.beta = min.beta,
              min.coef = min.coef, min.GIC = min.GIC))
}


# Function to tune:

Cal_GIC_pool_llh <- function(Y_lst, X_lst, fit_coef, length_lst, lambda_g,
                             type = 'BIC'){
  M <- length(Y_lst)
  p <- length(X_lst[[1]][1, ]) - 1
  mu <- fit_coef[1:(p + 1)]
  beta_fit <- c()
  alpha_fit <- c()
  
  for (i in 2:(M + 1)){
    beta_fit <- cbind(beta_fit, mu + fit_coef[((i - 1) * p + i):(i * p + i)])
    alpha_fit <- cbind(alpha_fit, fit_coef[((i - 1) * p + i):(i * p + i)])
  }

  llh <- 0
  X_weight_lst <- vector('list', M)
  for (i in 1:M){
    Y <- Y_lst[[i]]
    X <- X_lst[[i]]
    n <- length(Y)
    pi_vec <- as.vector(exp(X %*% beta_fit[,i]) / (1 + exp(X %*% beta_fit[,i])))
    X_weight_lst[[i]] <- diag(sqrt(pi_vec * (1 - pi_vec))) %*% X
    llh <- llh + sum(Y * log(pi_vec) + (1 - Y) * log(1 - pi_vec))
  } 
  
  X_all <- c()
  for (i in 1:M){
    X_all <- rbind(X_all, X_weight_lst[[i]])
  }
  
  X_alpha <- X_weight_lst[[1]]
  for (i in 2:M){
    X_alpha <- bdiag(X_alpha, X_weight_lst[[i]])
  }
  X_all <- cbind(X_all, X_alpha)
  
  enlarge_mat <- matrix(0, p + 1, p + 1)
  for (i in 1:M) {
    enlarge_mat <- cbind(enlarge_mat, diag(rep(10000, p + 1)))
  }
  X_enlarge <- rbind(X_all, enlarge_mat)
  
  norm_lst <- sqrt(rowSums(alpha_fit^2))
  S_alpha <- which(norm_lst != 0)
  S_mu <- which(mu != 0)
  S_full <- which(fit_coef != 0)
  H_S <- t(X_all[ ,S_full]) %*% X_all[ ,S_full]
  partial <- diag(0, length(S_full), length(S_full))
  if (length(S_alpha) != 0){
    for (t in 1:length(S_alpha)) {
      j <- S_alpha[t]
      partial_j <- t + length(S_mu) + length(S_alpha) * c(0:(M - 1))
      partial[partial_j, partial_j] <- lambda_g / norm_lst[j] * 
        (diag(1, M, M) - alpha_fit[j,] %*% t(alpha_fit[j,]) / norm_lst[j]^2)
    }
  }
  
  df <- sum(diag(solve(H_S + partial) %*% H_S))
  #print(df)
  if (type == 'BIC'){
    GIC <- - 2 * llh / sum(length_lst) + df * log(sum(length_lst)) / sum(length_lst)
  }
  if (type == 'mBIC'){
    GIC <- - 2 * llh / sum(length_lst) + df * log(max(exp(1), log(p))) * log(sum(length_lst)) / sum(length_lst)
  }
  if (type == 'AIC'){
    GIC <- - 2 * llh / sum(length_lst) + df * 2 / sum(length_lst)
  }
  
  return(list(GIC = GIC, beta = beta_fit))
}


# Function to evaluate the performance:

Evaluate_beta <- function(fit_coef, beta_fit, X_lst, beta_true_lst,
                          beta_center_true, beta_delta_true){
  M <- length(beta_true_lst)
  p <- length(beta_center_true) - 1
  
  homo_support <- which(beta_center_true[2:(p + 1)] != 0)
  heter_support <- which(beta_delta_true[[1]][2:(p + 1)] != 0)
  beta_support <- which(beta_true_lst[[1]][2:(p + 1)] != 0)
  
  true_num_homo <- length(intersect(homo_support, which(fit_coef[2:(p + 1)] != 0)))
  true_num_heter <- length(intersect(heter_support, which(fit_coef[(p + 3):(2 * p + 2)] != 0))) 
  true_num_beta <- length(intersect(beta_support, which(beta_fit[-1, 1] != 0)))
  
  
  beta_center_fit <- fit_coef[1:(p + 1)]
  beta_delta_fit <- c()
  for (m in 1:(M - 1)) {
    beta_delta_fit <- cbind(beta_delta_fit, fit_coef[((p + 1) * m + 1):((p + 1) * (m + 1))])
  }
  
  beta_center_abs <- abs(beta_center_fit)[-1]
  beta_delta_mean <- rowSums(beta_delta_fit^2)[-1]
 
  false_num_homo <- length(which(fit_coef[2:(p + 1)] != 0)) - true_num_homo
  false_num_heter <- length(which(fit_coef[(p + 3):(2 * p + 2)] != 0)) - true_num_heter
  false_num_beta <- length(which(beta_fit[-1, 1] != 0)) - true_num_beta
  
  pred_loss_lst <- c()
  est_loss_lst <- c()
  
  for (i in 1:M){
    pred_loss <- mean((cbind(rep(1, length(X_lst[[i]][,1])), X_lst[[i]]) %*% 
                         (beta_fit[,i] - beta_true_lst[[i]]))^2)
    est_loss <- sum(abs(beta_fit[,i] - beta_true_lst[[i]]))
    pred_loss_lst <- c(pred_loss_lst, pred_loss)
    est_loss_lst <- c(est_loss_lst, est_loss)
  }
  return(list(true.beta = true_num_beta, false.beta = false_num_beta,
              pred.loss = mean(pred_loss_lst), est.loss = sum(est_loss_lst)))
}


# Additional function for runing SMA + SIS


SIS <- function(Y_lst, X_lst, d, type = 'logit'){
  M <- length(Y_lst)
  S <- c()
  for (m in 1:M){
    Xm <- X_lst[[m]]
    Ym <- Y_lst[[m]]
    p <- length(Xm[1, ])
    n <- length(Xm[ ,1])
    Xm <- Xm - rep(1, n) %*% t(colMeans(Xm))
    Xm <- Xm / rep(1, n) %*% t(sqrt(colMeans(Xm^2)))
    if (type == 'logit'){
      import_m <- abs(colMeans(Xm[which(Ym == 1), ]) - 
                        colMeans(Xm[which(Ym == 0), ]))
      Sm <- which(rank(import_m) > p - d)
      S <- union(S, Sm)
    }
  }
  for (m in 1:M){
    X_lst[[m]] <- X_lst[[m]][ ,S]
  }
  return(list(S = S, X = X_lst, Y = Y_lst))
}


SIS_merge <- function(SMA_SIS_fit, support, p, intercept = T){
  M <- length(SMA_SIS_fit$min.beta[1, ])
  if (intercept == T){
    beta_sma <- matrix(0, p + 1, M)
    beta_sma[c(1, support + 1), ] <- SMA_SIS_fit$min.beta
    coef_sma <- matrix(0, (p + 1) * M, 1)
    q <- length(support)
    
    for (j in 1:M) {
      coef_sma[c(1 + (p + 1) * (j - 1), 
                 SIS_result$S + (p + 1) * (j - 1) + 1),] <- SMA_SIS_fit$min.coef[(q + 1) * (j - 1) + c(1:(q + 1)),]
    }
  }
  return(list(min.coef = coef_sma, min.beta = beta_sma))
}



