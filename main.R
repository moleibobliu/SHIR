
library(glmnet)

library(grplasso)
library(MASS)

source('function.R')


# Sample size at each site:

n <- 400

# Dimension:
p <- 800

# Number of sites:
M <- 4

r_lst <- c(0.15, 0.25, 0.35, 0.45)

# Generate data under Setting (i) in the SHIR paper:

simu_data <- Generate_data(n, p, M, magn_mu = 0.5, magn_alpha = 0.35, 
                           r_lst = r_lst, model = 'correct_sparse')


#### Other potential settings to try ####
"
# Generate data under Setting (ii) in the SHIR paper:

simu_data <- Generate_data(n, p, M, magn_mu = 0.2, magn_alpha = 0.15, 
                           r_lst = r_lst, model = 'correct_sparse')


# Generate data under Setting (iii) in the SHIR paper:

simu_data <- Generate_data(n, p, M, magn_mu = 0.5, magn_alpha = 0.35, 
                           r_lst = r_lst, model = 'correct_dense')


# Generate data under Setting (iv) in the SHIR paper:

simu_data <- Generate_data(n, p, M, magn_mu = 0.2, magn_alpha = 0.15, 
                           r_lst = r_lst, model = 'correct_dense')


# Generate data under Setting (v) in the SHIR paper:

simu_data <- Generate_data(n, p, M, magn_mu = 0.25, magn_alpha = 0.15, 
                           r_lst = r_lst, model = 'wrong')

"

# List of data for training:

X_lst <- simu_data$X
Y_lst <- simu_data$Y

beta_true_lst <- simu_data$beta_true
beta_center_true <- simu_data$beta_center_true
beta_delta_true <- simu_data$beta_delta_true


########## SHIR ##########

## Sample size of the sites

length_lst <- c()
for (m in 1:M){
  length_lst <- c(length_lst, length(X_lst[[m]][,1]))
}


## Local fit to derive the summary statistics

local_summary <- Local_fit(Y_lst, X_lst, family = 'binomial',
                           lambda_lst = 0.1 * c(4:20) * sqrt((log(p)) / n))

# U_lst: the derived gradient, I_lst: the derived hessian matrix

U_lst <- local_summary$U
I_lst <- local_summary$I

# Pre-specified candidat sets of tuning parameters:


lambda_lst = sqrt(n * log(p)) * 0.3 * c(7:30)
lambda_g_lst = c(0.75)


# Integrative regression

SHIR_train <- SHIR_fit(I_lst, U_lst, length_lst, lambda_lst = lambda_lst, 
                       lambda_g_lst = lambda_g_lst, tune = 'BIC')

SHIR_train$min.beta


# Evaluation:

SHIR_eval <- Evaluate_beta(SHIR_train$min.coef, SHIR_train$min.beta, X_lst, 
                           beta_true_lst, beta_center_true, beta_delta_true)

# absolute astimation error (AEE)

SHIR_eval$pred.loss

# prediction error (PE)

SHIR_eval$est.loss

# Number of true positive in discovering beta:

SHIR_eval$true.beta

# Number of false positive in discovering beta:

SHIR_eval$false.beta

# To obtain the fitted mean effect mu:

rowMeans(SHIR_train$min.beta)

# To obtain the fitted heterogeneous effect alpha (the m-th column loads alpha for the m-th site):

SHIR_train$min.beta - rowMeans(SHIR_train$min.beta)


###### Debias approach ######

tau_lst <- 0.1 * 5:15 * sqrt((log(p)) / n)
nu1_lst <- 0.2 * c(10:60) * sqrt(log(p) / sum(length_lst))
nu2_lst <- 0.2 * c(10:60) * sqrt((log(p) + M) / n)

# Integrative regression

debias_train <- debias_fit(I_lst, U_lst, local_summary$beta, n_lst = length_lst,
                           tau_lst = tau_lst, nu1_lst = nu1_lst, 
                           nu2_lst = nu2_lst, tune = 'BIC')

# Evaluation

debias_eval <- Evaluate_beta(debias_train$coef, debias_train$beta,
                             X_lst, beta_true_lst, beta_center_true, beta_delta_true)

############### IPD ###############


lambda_lst = sqrt(n * log(p)) * 0.2 * c(5:20)
lambda_g_lst =  c(0.75)

# Integrative regression

IPD_train <- IPD_fit(X_lst, Y_lst, length_lst, lambda_lst = lambda_lst, 
                     lambda_g_lst = lambda_g_lst, tune = 'BIC')

# Evaluation

IPD_eval <- Evaluate_beta(IPD_train$min.coef, IPD_train$min.beta, X_lst, 
                          beta_true_lst, beta_center_true, beta_delta_true)


############### SMA (+ SIS) ###############

# Screening at the local sites:

SIS_result <- SIS(Y_lst, X_lst, n / (3 * log(n)))

# Derive summary data with the selected set:

local_summary_sma <- Local_fit(SIS_result$Y, SIS_result$X, lambda_lst = 0)
Grad_lst_sma <- local_summary_sma$U
I_lst_sma <- local_summary_sma$I
p_ <- length(local_summary_sma$U[[1]])

# Integrative regression:


lambda_lst = sqrt(n * log(p_)) * 0.1 * c(3:20)
lambda_g_lst = c(0.75)


SMA_SIS_fit <- SHIR_fit(I_lst_sma, Grad_lst_sma, length_lst,
                        lambda_lst = lambda_lst, lambda_g_lst = lambda_g_lst, 
                        tune = 'BIC')

# Get back to obtain the 

SMA_train <- SIS_merge(SMA_SIS_fit, SIS_result$S, p)

SMA_SIS_eval <- Evaluate_beta(SMA_train$min.coef, SMA_train$min.beta, X_lst, 
                              beta_true_lst, beta_center_true, beta_delta_true)






