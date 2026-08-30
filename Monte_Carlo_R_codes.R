#rm(list = ls())


### Section 4.3
#######################

# Black-Scholes valuation of a European call
##############################################
S0 <- 80        # Initial stock price
K <- 90         # Strike price
r <- 0.05       # Risk-free rate
sigma <- 0.20   # Volatility
T <- 1          # Time to maturity
t <- 0          # Current time

tau<- T - t   # Time remaining to maturity

# Calculate d1 & d2
d1 <- (log(S0 / K) + (r + 0.5 * sigma^2) * tau) /
  (sigma * sqrt(tau))

# Calculate d2
d2 <- d1 - sigma * sqrt(tau)

# Calculate N(d1) and N(d2)
N_d1 <- pnorm(d1)
N_d2 <- pnorm(d2)

# Black-Scholes call price
call_BS <- S0 * N_d1 - K * exp(-r * tau) * N_d2

# Black-Scholes European put price
put_BS <- K * exp(-r * tau) * pnorm(-d2) - S0 * pnorm(-d1)


# Display results
call_BS   #  4.164392
put_BS    # 9.7751


#   MC Option pricing under the risk-neutral measure r
##################################################################
# Simulate stock prices under Geometric Brownian Motion (GBM) model. 
# The GBM model is a widely used model in finance to simulate stock prices over time.
##============================================================
# Set parameters
S0 <- 80        # Initial stock price
r <- 0.05       # Risk-free rate
sigma <- 0.20   # Annual volatility
T <- 1         # Time to maturity (years)
M <- 252        # Number of monitoring dates


dt <- T / M     # Time step

# Generate standard normal random variables
set.seed(123)   # For reproducibility
Z <- rnorm(M)

# Simulate stock prices using exact GBM discretisation
S <- numeric(M + 1)
S[1] <- S0

for (t in 1:M) {
  S[t + 1] <- S[t] * exp((r - 0.5 * sigma^2) * dt +  sigma * sqrt(dt) * Z[t]) }

##====================================================================
# Simulate 100,000 asset-price paths under the risk-neutral GBM model and compute the Monte Carlo estimates for European options.
##====================================================================

# Simulate N asset-price paths under the risk-neutral GBM model

S0 <- 80
r <- 0.05
sigma <- 0.20
T <- 1
N <- 100000        # Number of simulated paths
M <- 252          # Number of monitoring dates

dt <- T / M

set.seed(123)

# Generate standard normal random variables
Z <- matrix(rnorm(N * M), nrow = N, ncol = M)

# Initialise asset-price matrix
prices <- matrix(0, nrow = N, ncol = M + 1)
prices[, 1] <- S0

# Simulate asset-price paths using exact GBM discretisation
for (t in 1:M) {
  prices[, t + 1] <- prices[, t] * exp((r - 0.5 * sigma^2) * dt + sigma * sqrt(dt) * Z[, t])}


# Plot the first 50 simulated paths

matplot( 0:M, t(prices[1:50, ]), type = "l", lty = 1, xlab = "Monitoring Date",
         ylab = "Asset Price", main = "50 Simulated Asset-Price Paths under GBM")


# European option payoffs using the simulated asset-price paths
####==========================================================
K <- 90
# Terminal stock price from all 100,000 simulated paths
ST <- prices[, M + 1]

# European call payoff (payoff from each simulated path)
call_payoff <- pmax(ST - K, 0)

# European put payoff
put_payoff <- pmax(K - ST, 0)


# Monte Carlo option prices
call_MC <- exp(-r * T) * mean(call_payoff)
put_Mc  <- exp(-r * T) * mean(put_payoff)

call_MC  #  4.167621
put_Mc   #   9.763013

##=========================
# Monte Carlo standard error
call_se <- exp(-r * T) * sd(call_payoff) / sqrt(N)
put_se  <- exp(-r * T) * sd(put_payoff) / sqrt(N)

call_se    #  0.02750649
put_se     # 0.03224782

#================================================
# Absolute errors between Monte Carlo and Black-Scholes prices
#=======================================================
abs_error_call <- abs(call_MC - call_BS)
abs_error_put  <- abs(put_Mc - put_BS)

abs_error_call  # 0.003228792
abs_error_put   # 0.01202713

####################################################
# Section 4.5
###################################################
# Plot first 30 simulated asset-price paths
##########################################
matplot(0:M, t(prices[1:30, ]), type = "l", lty = 1, xlab = "Monitoring Date",ylab = "Asset Price",
        main = "Simulated Asset-Price Paths under GBM")

## Using Time instead of monitoring dates for x-axis
time <- seq(0, T, length.out = M + 1)

matplot(time, t(prices[1:20, ]), type = "l", lty = 1, xlab = "Time (Years)",
  ylab = "Asset Price", main = "Simulated Asset-Price Paths under GBM")



#######################################
## Section 4.4
########################################
# Monte Carlo Sensitivity to Number of Simulations (European Option)
# ============================================================

# Parameters
S0 <- 80
K <- 90
r <- 0.05
sigma <- 0.20
T <- 1

# Simulation sizes
N_values <- c(1000, 5000, 10000, 100000, 200000)

# Store results
call_MC <- numeric(length(N_values))
put_MC  <- numeric(length(N_values))
call_SE <- numeric(length(N_values))
put_SE  <- numeric(length(N_values))

set.seed(123)

# =====
# Monte Carlo simulation
# =========

for (k in seq_along(N_values)) {
  N <- N_values[k]
  
  # Generate terminal stock prices directly from exact GBM
  Z <- rnorm(N)
  
  ST <- S0 * exp((r - 0.5 * sigma^2) * T + sigma * sqrt(T) * Z )
  
  # European call payoff
  call_payoff <- pmax(ST - K, 0)
  
  # European put payoff
  put_payoff <- pmax(K - ST, 0)
  
  # Monte Carlo option prices
  call_MC[k] <- exp(-r * T) * mean(call_payoff)
  put_MC[k]  <- exp(-r * T) * mean(put_payoff)
  
  # Monte Carlo standard errors
  call_SE[k] <- exp(-r * T) * sd(call_payoff) / sqrt(N)
  
  put_SE[k] <- exp(-r * T) * sd(put_payoff) / sqrt(N)
}

# ======
# Black-Scholes benchmark
# =====

d1 <- (log(S0 / K) + (r + 0.5 * sigma^2) * T) / (sigma * sqrt(T))

d2 <- d1 - sigma * sqrt(T)

call_BS <- S0 * pnorm(d1) - K * exp(-r * T) * pnorm(d2)

put_BS <- K * exp(-r * T) * pnorm(-d2) - S0 * pnorm(-d1)

# ========
# Absolute errors
# ========

call_error <- abs(call_MC - call_BS)
put_error  <- abs(put_MC - put_BS)

# 95% Confidence Intervals for European Options
# ============================================================

# Critical value for 95% confidence interval
z <- 1.96

# European Call
call_lower <- call_MC - z * call_SE
call_upper <- call_MC + z * call_SE

# European Put
put_lower <- put_MC - z * put_SE
put_upper <- put_MC + z * put_SE

# Create results table  (Black- Scholes included for reference)
CI_results <- data.frame(
  N = N_values,
  
  Call_BS = round(call_BS, 4),
  Call_MC = round(call_MC, 4),
  Call_SE = round(call_SE, 4),
  Call_Lower_95 = round(call_lower, 4),
  Call_Upper_95 = round(call_upper, 4),
  
  Put_BS = round(put_BS, 4),
  Put_MC = round(put_MC, 4),
  Put_SE = round(put_SE, 4),
  Put_Lower_95 = round(put_lower, 4),
  Put_Upper_95 = round(put_upper, 4)
)

CI_results

#    N     Call_BS Call_MC Call_SE     Call_Lower_95  Call_Upper_95 Put_BS Put_MC Put_SE     Put_Lower_95
# 1 1e+03  4.1644  4.2217  0.2817        3.6696        4.7738       9.775 9.5956 0.3175       8.9733
# 2 5e+03  4.1644  4.1150  0.1223        3.8753        4.3547       9.775 9.6943 0.1442       9.4117
# 3 1e+04  4.1644  4.1822  0.0867        4.0123        4.3521       9.775 9.8052 0.1022       9.6050
# 4 1e+05  4.1644  4.1841  0.0275        4.1302        4.2381       9.775 9.7484 0.0323       9.6851
# 5 2e+05  4.1644  4.1870  0.0196        4.1486        4.2254       9.775 9.7758 0.0228       9.7310
# 6 5e+05  4.1644  4.1562  0.0123        4.1322        4.1802       9.775 9.7840 0.0144       9.7557       


#     Put_Upper_95
# 1      10.2178
# 2       9.9769
# 3      10.0054
# 4       9.8117
# 5       9.8205


##===================================================================
# Plot: Monte Carlo Standard Error versus Simulation Size N
#=====================================================================

plot(N_values, call_SE, type = "b", pch = 19, col= "red", log = "x",
     xlab = "Number of Simulations (N)", ylab = "Monte Carlo Standard Error",
     main = "Monte Carlo Standard Error versus Simulation Size")

lines(N_values, put_SE, type = "b",col="black", pch = 17, lty = 2)

legend( "topright", legend = c(  "European Call", "European Put" ),col=c("red","black"), pch = c(19, 17),
        lty = c(1, 2))






###########################################
 ## Section 4.6
#######################

# Lookback Options Valuation
###====================
S0 <- 80
K <- 90
r <- 0.05
q <- 0
sigma <- 0.20
T <- 1

N <- 100000       # number of paths
M <- 252           # monitoring dates

dt <- T / M

set.seed(123)

# Initial values
S <- rep(S0, N)

# Running maximum and minimum
S_max <- S
S_min <- S

# Simulate one time step at a time

for (j in 1:M) {
  Z <- rnorm(N)
  S <- S * exp( (r - q - 0.5 * sigma^2) * dt + sigma * sqrt(dt) * Z )
  
  # Update running maximum and minimum
  S_max <- pmax(S_max, S)
  S_min <- pmin(S_min, S)
}

# FLOATING-STRIKE LOOKBACK PAYOFFS
Floating_Call_Payoff <- S - S_min
Floating_Put_Payoff  <- S_max - S


# FIXED-STRIKE LOOKBACK PAYOFFS
Fixed_Call_Payoff <- pmax(S_max - K, 0)
Fixed_Put_Payoff  <- pmax(K - S_min, 0)

# MONTE CARLO PRICES
MC_Floating_Call <- exp(-r * T) * mean(Floating_Call_Payoff)

MC_Floating_Put <- exp(-r * T) * mean(Floating_Put_Payoff)

MC_Fixed_Call <- exp(-r * T) * mean(Fixed_Call_Payoff)

MC_Fixed_Put <- exp(-r * T) * mean(Fixed_Put_Payoff)

# STANDARD ERRORS
####==========================

SE_Floating_Call <- exp(-r * T) * sd(Floating_Call_Payoff) / sqrt(N)

SE_Floating_Put <- exp(-r * T) * sd(Floating_Put_Payoff) / sqrt(N)

SE_Fixed_Call <- exp(-r * T) * sd(Fixed_Call_Payoff) / sqrt(N)

SE_Fixed_Put <- exp(-r * T) * sd(Fixed_Put_Payoff) / sqrt(N)


# DISPLAY RESULTS
results <- data.frame(
  
  Option = c("Floating-Strike Call", "Floating-Strike Put", "Fixed-Strike Call", "Fixed-Strike Put" ),
  
  MC_Value = c( MC_Floating_Call,  MC_Floating_Put,  MC_Fixed_Call,  MC_Fixed_Put),
  
  Standard_Error = c( SE_Floating_Call, SE_Floating_Put, SE_Fixed_Call, SE_Fixed_Put) )

results$MC_Value <- round(results$MC_Value, 4)
results$Standard_Error <- round(results$Standard_Error, 4)

print(results)

#   Option                MC_Value     Standard_Error
# 1 Floating-Strike Call  13.3100         0.0368
# 2 Floating-Strike Put   10.7734         0.0249
# 3 Fixed-Strike Call     7.3506          0.0332
# 4 Fixed-Strike Put      18.9054         0.0227


# 95% Confidence Intervals for Lookback Options
# ================================================
# Floating-strike confidence intervals
float_call_CI <- c( LB_float_call - 1.96 * SE_float_call, LB_float_call + 1.96 * SE_float_call )
float_put_CI <- c( LB_float_put - 1.96 * SE_float_put,  LB_float_put + 1.96 * SE_float_put )

# Fixed-strike confidence intervals
fixed_call_CI <- c( LB_fixed_call - 1.96 * SE_fixed_call, LB_fixed_call + 1.96 * SE_fixed_call )

fixed_put_CI <- c( LB_fixed_put - 1.96 * SE_fixed_put, LB_fixed_put + 1.96 * SE_fixed_put )

###=============
CI_results <- data.frame( Option = c( "Floating Call", "Floating Put", "Fixed Call", "Fixed Put" ),
                          Lower_95 = c(  float_call_CI[1],  float_put_CI[1],  fixed_call_CI[1],  fixed_put_CI[1] ),
                          Upper_95 = c( float_call_CI[2], float_put_CI[2], fixed_call_CI[2], fixed_put_CI[2] ) )

CI_results$Lower_95 <- round(CI_results$Lower_95, 4)
CI_results$Upper_95 <- round(CI_results$Upper_95, 4)
CI_results

#   Option         Lower_95 Upper_95
# 1 Floating Call  15.2921  15.4746
# 2 Floating Put   14.5400  14.7195
# 3 Fixed Call     7.4984   7.6297
# 4 Fixed Put      19.0736  19.1622

############################################################
# LOOKBACK OPTIONS: EFFECT OF INCREASING TIME STEPS
# Discrete Monitoring vs Continuous-Monitoring Analytical Benchmarks
############################################################
# Model Parameters unchanged
n_paths <- 100000

# Increasing number of monitoring dates
M_values <- c( 252, 504, 1008, 2016)

set.seed(123)

# ANALYTICAL BENCHMARK VALUES
#Closed-form analytical values obtained from the Black-Scholes lookback formulas.

Analytical_Floating_Call <- 13.7734
Analytical_Floating_Put  <- 11.4325
Analytical_Fixed_Call    <- 7.7458
Analytical_Fixed_Put     <- 19.3841

#LOOKBACK MONTE CARLO FUNCTION

simulate_lookback <- function(M, n_paths) {
  
  dt <- T / M
  
  # Generate random numbers
  Z <- matrix(
    rnorm(n_paths * M),
    nrow = n_paths,
    ncol = M
  )
  
  # Initial price
  S <- rep(S0, n_paths)
  
  # Running maximum and minimum
  S_max <- rep(S0, n_paths)
  S_min <- rep(S0, n_paths)
  
  
  # SIMULATE GBM PATHS USING EXACT TRANSITION
  
  for (j in 1:M) {
    
    S <- S * exp(
      (r - q - 0.5 * sigma^2) * dt +
        sigma * sqrt(dt) * Z[, j]
    )
    
    # Update running maximum and minimum
    S_max <- pmax(S_max, S)
    S_min <- pmin(S_min, S)
  }
  
  # LOOKBACK PAYOFFS
  floating_call_payoff <- S - S_min
  
  floating_put_payoff <- S_max - S
  
  fixed_call_payoff <- pmax(S_max - K, 0)
  
  fixed_put_payoff <- pmax(K - S_min, 0)
  
  # DISCOUNT PAYOFFS
  FC <- exp(-r * T) * floating_call_payoff
  FP <- exp(-r * T) * floating_put_payoff
  XC <- exp(-r * T) * fixed_call_payoff
  XP <- exp(-r * T) * fixed_put_payoff
  
  # MONTE CARLO PRICES
  FC_MC <- mean(FC)
  FP_MC <- mean(FP)
  XC_MC <- mean(XC)
  XP_MC <- mean(XP)
  
  # STANDARD ERRORS
  FC_SE <- sd(FC) / sqrt(n_paths)
  FP_SE <- sd(FP) / sqrt(n_paths)
  XC_SE <- sd(XC) / sqrt(n_paths)
  XP_SE <- sd(XP) / sqrt(n_paths)
  
  # 95% CONFIDENCE INTERVALS
  FC_L <- FC_MC - 1.96 * FC_SE
  FC_U <- FC_MC + 1.96 * FC_SE
  
  FP_L <- FP_MC - 1.96 * FP_SE
  FP_U <- FP_MC + 1.96 * FP_SE
  
  XC_L <- XC_MC - 1.96 * XC_SE
  XC_U <- XC_MC + 1.96 * XC_SE
  
  XP_L <- XP_MC - 1.96 * XP_SE
  XP_U <- XP_MC + 1.96 * XP_SE
  
  # ERRORS RELATIVE TO ANALYTICAL VALUES
  FC_AE <- abs(FC_MC - Analytical_Floating_Call)
  FP_AE <- abs(FP_MC - Analytical_Floating_Put)
  XC_AE <- abs(XC_MC - Analytical_Fixed_Call)
  XP_AE <- abs(XP_MC - Analytical_Fixed_Put)
  
  FC_PE <- 100 * FC_AE / Analytical_Floating_Call
  FP_PE <- 100 * FP_AE / Analytical_Floating_Put
  XC_PE <- 100 * XC_AE / Analytical_Fixed_Call
  XP_PE <- 100 * XP_AE / Analytical_Fixed_Put
  
  # RETURN RESULTS
  data.frame(
    
    M = M,
    
    Floating_Call_MC = FC_MC,
    Floating_Call_SE = FC_SE,
    Floating_Call_Lower_95 = FC_L,
    Floating_Call_Upper_95 = FC_U,
    Floating_Call_Absolute_Error = FC_AE,
    Floating_Call_Percentage_Error = FC_PE,
    
    Floating_Put_MC = FP_MC,
    Floating_Put_SE = FP_SE,
    Floating_Put_Lower_95 = FP_L,
    Floating_Put_Upper_95 = FP_U,
    Floating_Put_Absolute_Error = FP_AE,
    Floating_Put_Percentage_Error = FP_PE,
    
    Fixed_Call_MC = XC_MC,
    Fixed_Call_SE = XC_SE,
    Fixed_Call_Lower_95 = XC_L,
    Fixed_Call_Upper_95 = XC_U,
    Fixed_Call_Absolute_Error = XC_AE,
    Fixed_Call_Percentage_Error = XC_PE,
    
    Fixed_Put_MC = XP_MC,
    Fixed_Put_SE = XP_SE,
    Fixed_Put_Lower_95 = XP_L,
    Fixed_Put_Upper_95 = XP_U,
    Fixed_Put_Absolute_Error = XP_AE,
    Fixed_Put_Percentage_Error = XP_PE
  )
}

# RUN FOR INCREASING M
results_list <- list()

for (i in seq_along(M_values)) {
  
  M <- M_values[i]
  
  cat("Running M =", M, "\n")
  
  results_list[[i]] <-
    simulate_lookback(M, n_paths)
}

# COMBINE RESULTS
results_M <- do.call(rbind, results_list)
rownames(results_M) <- NULL
results_M[, -1] <-round(results_M[, -1], 4)


print(results_M)

# ============================================================
# Continuous vs Discrete Monitoring of a Lookback Option
# ============================================================

t <- seq(0, 1, length.out = 1000)

# Illustrative continuous asset-price path
S <- 100 + 7 * sin(2 * pi * t) + 20 * exp(-((t - 0.48) / 0.07)^2)

# Discrete monitoring dates
# monitoring_times <- c(0.25, 0.50, 0.75)
monitoring_times <- c(0.1, 0.2, 0.3, 0.4, 0.5,0.6, 0.7, 0.8, 0.9, 1 )

# Prices observed only at monitoring dates
monitoring_prices <- approx(
  t, S,
  xout = monitoring_times
)$y

# Continuous maximum
i_max <- which.max(S)
t_max <- t[i_max]
S_max <- S[i_max]

# Discrete maximum
i_disc_max <- which.max(monitoring_prices)
t_disc_max <- monitoring_times[i_disc_max]
S_disc_max <- monitoring_prices[i_disc_max]

# ------------------------------------------------------------
# Plot
# ------------------------------------------------------------

plot(t, S,
     type = "l",
     lwd = 2,
     xlab = "Time",
     ylab = "Asset Price",
     main = "Continuous and Discrete Monitoring of a Lookback Option")

# Monitoring dates
abline(v = monitoring_times,
       lty = 2)

# Discrete observations
points(monitoring_times,
       monitoring_prices,
       pch = 19,
       cex = 1.2)

# True continuous maximum
points(t_max,
       S_max,
       pch = 19,
       cex = 1.4)

# Discrete maximum
points(t_disc_max,
       S_disc_max,
       pch = 17,
       cex = 1.4)

# Lines indicating the two maxima
abline(h = S_max, lty = 3)
abline(h = S_disc_max, lty = 3)

# Labels
text(t_max, S_max,
     labels = "Continuous maximum",
     pos = 3)

text(t_disc_max, S_disc_max,
     labels = "Discrete maximum",
     pos = 1)

# Monitoring labels
text(monitoring_times,
     monitoring_prices,
     labels = c(expression(t[1]),
                expression(t[2]),
                expression(t[3])),
     pos = 3)

################################
# Effect of M on LB option Prices
##############################################
n_paths_values <-100000

# Increasing number of monitoring dates
M_values <- c(1000, 2000, 3000, 5000, 10000, 20000)

set.seed(123)

# ANALYTICAL BENCHMARK VALUES
Analytical_Floating_Call <- 13.7734
Analytical_Floating_Put  <- 11.4325
Analytical_Fixed_Call    <- 7.7458
Analytical_Fixed_Put     <- 19.3841

#  LOOKBACK MONTE CARLO FUNCTION
simulate_lookback <- function(M, n_paths, batch_size = 5000) {
  
  dt <- T / M
  
  # Number of batches
  n_batches <- ceiling(n_paths / batch_size)
  
  # Storage for discounted payoffs
  FC_all <- numeric(n_paths)
  FP_all <- numeric(n_paths)
  XC_all <- numeric(n_paths)
  XP_all <- numeric(n_paths)
  
  start <- 1
  
  # SIMULATE IN BATCHES
  for (b in 1:n_batches) {
    
    # Number of paths in this batch
    end <- min(start + batch_size - 1, n_paths)
    n_batch <- end - start + 1
    
    # Generate random numbers only for this batch
    Z <- matrix(
      rnorm(n_batch * M),
      nrow = n_batch,
      ncol = M
    )
    
    # Initial price
    S <- rep(S0, n_batch)
    
    # Running maximum and minimum
    S_max <- rep(S0, n_batch)
    S_min <- rep(S0, n_batch)
    
    # EXACT GBM TRANSITION
    drift <- (r - q - 0.5 * sigma^2) * dt
    diffusion <- sigma * sqrt(dt)
    
    for (j in 1:M) {
      
      S <- S * exp(drift + diffusion * Z[, j] )
      
      S_max <- pmax(S_max, S)
      S_min <- pmin(S_min, S)
    }
    
    # LOOKBACK PAYOFFS
    floating_call_payoff <- S - S_min
    floating_put_payoff  <- S_max - S
    
    fixed_call_payoff <- pmax(S_max - K, 0)
    fixed_put_payoff  <- pmax(K - S_min, 0)
    
    # DISCOUNT
    discount <- exp(-r * T)
    
    FC_all[start:end] <- discount * floating_call_payoff
    FP_all[start:end] <- discount * floating_put_payoff
    XC_all[start:end] <- discount * fixed_call_payoff
    XP_all[start:end] <- discount * fixed_put_payoff
    
    start <- end + 1
    
    # Progress indicator
    cat(
      "\r  Batch", b, "of", n_batches,
      "for M =", M,
      "and paths =", n_paths
    )
  }
  
  cat("\n")
  
  # MONTE CARLO PRICES
  FC_MC <- mean(FC_all)
  FP_MC <- mean(FP_all)
  XC_MC <- mean(XC_all)
  XP_MC <- mean(XP_all)
  
  # STANDARD ERRORS
  FC_SE <- sd(FC_all) / sqrt(n_paths)
  FP_SE <- sd(FP_all) / sqrt(n_paths)
  XC_SE <- sd(XC_all) / sqrt(n_paths)
  XP_SE <- sd(XP_all) / sqrt(n_paths)
  
  # 95% CONFIDENCE INTERVALS
  
  FC_L <- FC_MC - 1.96 * FC_SE
  FC_U <- FC_MC + 1.96 * FC_SE
  
  FP_L <- FP_MC - 1.96 * FP_SE
  FP_U <- FP_MC + 1.96 * FP_SE
  
  XC_L <- XC_MC - 1.96 * XC_SE
  XC_U <- XC_MC + 1.96 * XC_SE
  
  XP_L <- XP_MC - 1.96 * XP_SE
  XP_U <- XP_MC + 1.96 * XP_SE
  
  # ABSOLUTE ERRORS
  FC_AE <- abs(FC_MC - Analytical_Floating_Call)
  FP_AE <- abs(FP_MC - Analytical_Floating_Put)
  XC_AE <- abs(XC_MC - Analytical_Fixed_Call)
  XP_AE <- abs(XP_MC - Analytical_Fixed_Put)
  
  # PERCENTAGE ERRORS
  FC_PE <- 100 * FC_AE / Analytical_Floating_Call
  FP_PE <- 100 * FP_AE / Analytical_Floating_Put
  XC_PE <- 100 * XC_AE / Analytical_Fixed_Call
  XP_PE <- 100 * XP_AE / Analytical_Fixed_Put
  
  # RETURN RESULTS
  data.frame(
    
    M = M,
    n_paths = n_paths,
    
    Floating_Call_MC = FC_MC,
    Floating_Call_SE = FC_SE,
    Floating_Call_Lower_95 = FC_L,
    Floating_Call_Upper_95 = FC_U,
    Floating_Call_Absolute_Error = FC_AE,
    Floating_Call_Percentage_Error = FC_PE,
    
    Floating_Put_MC = FP_MC,
    Floating_Put_SE = FP_SE,
    Floating_Put_Lower_95 = FP_L,
    Floating_Put_Upper_95 = FP_U,
    Floating_Put_Absolute_Error = FP_AE,
    Floating_Put_Percentage_Error = FP_PE,
    
    Fixed_Call_MC = XC_MC,
    Fixed_Call_SE = XC_SE,
    Fixed_Call_Lower_95 = XC_L,
    Fixed_Call_Upper_95 = XC_U,
    Fixed_Call_Absolute_Error = XC_AE,
    Fixed_Call_Percentage_Error = XC_PE,
    
    Fixed_Put_MC = XP_MC,
    Fixed_Put_SE = XP_SE,
    Fixed_Put_Lower_95 = XP_L,
    Fixed_Put_Upper_95 = XP_U,
    Fixed_Put_Absolute_Error = XP_AE,
    Fixed_Put_Percentage_Error = XP_PE
  )
}

# RUN ALL COMBINATIONS
results_list <- list()

counter <- 1

for (n in n_paths_values) {
  
  for (M in M_values) {
    
    cat("\n============================================\n")
    cat("Running n_paths =", n, "\n")
    cat("Running M       =", M, "\n")
    cat("============================================\n")
    
    results_list[[counter]] <-
      simulate_lookback(
        M = M,
        n_paths = n,
        batch_size = 5000
      )
    
    counter <- counter + 1
  }
}

# COMBINE RESULTS
results_M <- do.call(rbind, results_list)

rownames(results_M) <- NULL

# ROUND RESULTS
results_M[, -c(1, 2)] <-
  round(results_M[, -c(1, 2)], 4)

print(results_M)



###########################################################
###########################################################
# Section 4.7
#################

#####################################################################
# COMPARISON OF DISCRETISATION METHODS
# Exact GBM vs Euler-Maruyama vs Milstein
############################################################
#  MODEL PARAMETERS

S0    <- 80
K     <- 90
r     <- 0.05
q     <- 0
sigma <- 0.20
T     <- 1

n_paths <- 100000

# Number of time steps to investigate
M_values <- c(126, 252, 504, 1008, 2016)
set.seed(123)

# BLACK-SCHOLES ANALYTICAL VALUES
d1 <- (log(S0 / K) +(r - q + 0.5 * sigma^2) * T) /(sigma * sqrt(T))

d2 <- d1 - sigma * sqrt(T)

BS_Call <- S0 * exp(-q * T) * pnorm(d1) - K * exp(-r * T) * pnorm(d2)

BS_Put <- K * exp(-r * T) * pnorm(-d2) - S0 * exp(-q * T) * pnorm(-d1)

cat("Black-Scholes Call =", round(BS_Call, 4), "\n")
cat("Black-Scholes Put  =", round(BS_Put, 4), "\n")

# FUNCTION FOR ONE DISCRETISATION METHOD
############################################################

simulate_european <- function(method, M, Z) {
  
  dt <- T / M
  
  S <- rep(S0, n_paths)
  
  for (j in 1:M) {
    
    Zj <- Z[, j]
    
    if (method == "Exact") {S <- S * exp((r - q - 0.5 * sigma^2) * dt + sigma * sqrt(dt) * Zj )
    
    } else if (method == "Euler-Maruyama") { S <- S + (r - q) * S * dt + sigma * S * sqrt(dt) * Zj
    
    # Prevent negative stock prices
    S <- pmax(S, 0)} else if (method == "Milstein") { S <- S + (r - q) * S * dt + sigma * S * sqrt(dt) * Zj + 0.5 * sigma^2 * S * (Zj^2 - 1) * dt
    
    # Prevent negative stock prices
    S <- pmax(S, 0) }
  }
  
  
  # PAYOFFS
  call_payoff <- pmax(S - K, 0)
  put_payoff  <- pmax(K - S, 0)
  
  # DISCOUNTED PAYOFFS
  discounted_call <- exp(-r * T) * call_payoff
  discounted_put  <- exp(-r * T) * put_payoff
  
  # MONTE CARLO VALUES
  Call_MC <- mean(discounted_call)
  Put_MC  <- mean(discounted_put)
  
  # STANDARD ERRORS
  Call_SE <- sd(discounted_call) / sqrt(n_paths)
  Put_SE  <- sd(discounted_put) / sqrt(n_paths)
  
  
  # 95% CONFIDENCE INTERVALS
  Call_Lower <- Call_MC - 1.96 * Call_SE
  Call_Upper <- Call_MC + 1.96 * Call_SE
  
  Put_Lower <- Put_MC - 1.96 * Put_SE
  Put_Upper <- Put_MC + 1.96 * Put_SE
  
  # ERRORS
  Call_Abs_Error <- abs(Call_MC - BS_Call)
  Put_Abs_Error  <- abs(Put_MC - BS_Put)
  
  Call_Pct_Error <- 100 * Call_Abs_Error / BS_Call
  Put_Pct_Error  <- 100 * Put_Abs_Error / BS_Put
  
  # RETURN RESULTS
  return(data.frame( Method = method,
                     M = M,
                     Call_MC = Call_MC, Call_SE = Call_SE,
                     Call_Lower_95 = Call_Lower, Call_Upper_95 = Call_Upper,
                     Call_Absolute_Error = Call_Abs_Error, Call_Percentage_Error = Call_Pct_Error,
                     
                     Put_MC = Put_MC, Put_SE = Put_SE,
                     Put_Lower_95 = Put_Lower, Put_Upper_95 = Put_Upper,
                     Put_Absolute_Error = Put_Abs_Error, Put_Percentage_Error = Put_Pct_Error )) }

#  RUN ALL THREE METHODS
########################################
results_list <- list()

counter <- 1

for (M in M_values) { # Generate common random numbers
  Z <- matrix( rnorm(n_paths * M),  nrow = n_paths, ncol = M  )
  
  # Exact GBM
  results_list[[counter]] <- simulate_european("Exact", M, Z)
  counter <- counter + 1
  
  # Euler-Maruyama
  results_list[[counter]] <- simulate_european("Euler-Maruyama", M, Z)
  counter <- counter + 1
  
  # Milstein
  results_list[[counter]] <- simulate_european("Milstein", M, Z)
  counter <- counter + 1
}

# COMBINE RESULTS
###############
results <- do.call(rbind, results_list)
rownames(results) <- NULL

#ROUND RESULTS FOR PRESENTATION
results$Call_MC <- round(results$Call_MC, 4)
results$Call_SE <- round(results$Call_SE, 5)
results$Call_Lower_95 <- round(results$Call_Lower_95, 4)
results$Call_Upper_95 <- round(results$Call_Upper_95, 4)
results$Call_Absolute_Error <- round(
  results$Call_Absolute_Error, 4 )
results$Call_Percentage_Error <- round(
  results$Call_Percentage_Error, 4 )

results$Put_MC <- round(results$Put_MC, 4)
results$Put_SE <- round(results$Put_SE, 5)
results$Put_Lower_95 <- round(results$Put_Lower_95, 4)
results$Put_Upper_95 <- round(results$Put_Upper_95, 4)
results$Put_Absolute_Error <- round(
  results$Put_Absolute_Error, 4 )
results$Put_Percentage_Error <- round(
  results$Put_Percentage_Error, 4 )

print(results)

############################################################
# LOOKBACK OPTIONS: EXACT GBM vs EULER-MARUYAMA vs MILSTEIN
############################################################
n_paths <- 100000

M_values <- c(126, 252, 504, 1008)

set.seed(123)

# SIMULATION FUNCTION
####========================
simulate_lookback <- function(method, M, Z) {
  dt <- T / M
  
  # Initial stock price
  S <- rep(S0, n_paths)
  
  # Running maximum and minimum
  S_max <- rep(S0, n_paths)
  S_min <- rep(S0, n_paths)
  
  # SIMULATE PRICE PATHS
  
  for (j in 1:M) { Zj <- Z[, j]
  
  if (method == "Exact") {S <- S * exp( (r - q - 0.5 * sigma^2) * dt + sigma * sqrt(dt) * Zj )
  
  } else if (method == "Euler-Maruyama") { S <- S + (r - q) * S * dt +  sigma * S * sqrt(dt) * Zj
  
  # Avoid negative prices
  S <- pmax(S, 0)
  
  } else if (method == "Milstein") { S <- S + (r - q) * S * dt +  sigma * S * sqrt(dt) * Zj + 0.5 * sigma^2 * S * (Zj^2 - 1) * dt
  
  # Avoid negative prices
  S <- pmax(S, 0)
  }
  
  # Update running maximum and minimum
  S_max <- pmax(S_max, S)
  S_min <- pmin(S_min, S)
  }
  
  # LOOKBACK PAYOFFS
  ############
  
  # Floating-strike call
  floating_call_payoff <- S - S_min
  
  # Floating-strike put
  floating_put_payoff <- S_max - S
  
  # Fixed-strike call
  fixed_call_payoff <- pmax(S_max - K, 0)
  
  # Fixed-strike put
  fixed_put_payoff <- pmax(K - S_min, 0)
  
  # DISCOUNTED PAYOFFS
  ##############
  FC <- exp(-r * T) * floating_call_payoff
  FP <- exp(-r * T) * floating_put_payoff
  
  XC <- exp(-r * T) * fixed_call_payoff
  XP <- exp(-r * T) * fixed_put_payoff
  
  # MONTE CARLO PRICES
  #####################
  FC_price <- mean(FC)
  FP_price <- mean(FP)
  
  XC_price <- mean(XC)
  XP_price <- mean(XP)
  
  
  # STANDARD ERRORS
  ##########
  FC_SE <- sd(FC) / sqrt(n_paths)
  FP_SE <- sd(FP) / sqrt(n_paths)
  
  XC_SE <- sd(XC) / sqrt(n_paths)
  XP_SE <- sd(XP) / sqrt(n_paths)
  
  # 95% CONFIDENCE INTERVALS
  ############
  FC_Lower <- FC_price - 1.96 * FC_SE
  FC_Upper <- FC_price + 1.96 * FC_SE
  
  FP_Lower <- FP_price - 1.96 * FP_SE
  FP_Upper <- FP_price + 1.96 * FP_SE
  
  XC_Lower <- XC_price - 1.96 * XC_SE
  XC_Upper <- XC_price + 1.96 * XC_SE
  
  XP_Lower <- XP_price - 1.96 * XP_SE
  XP_Upper <- XP_price + 1.96 * XP_SE
  
  
  # RETURN RESULTS
  ########
  data.frame( Method = method, M = M,
              
              Floating_Call_MC = FC_price,
              Floating_Call_SE = FC_SE,
              Floating_Call_Lower_95 = FC_Lower,
              Floating_Call_Upper_95 = FC_Upper,
              
              Floating_Put_MC = FP_price,
              Floating_Put_SE = FP_SE,
              Floating_Put_Lower_95 = FP_Lower,
              Floating_Put_Upper_95 = FP_Upper,
              
              Fixed_Call_MC = XC_price,
              Fixed_Call_SE = XC_SE,
              Fixed_Call_Lower_95 = XC_Lower,
              Fixed_Call_Upper_95 = XC_Upper,
              
              Fixed_Put_MC = XP_price,
              Fixed_Put_SE = XP_SE,
              Fixed_Put_Lower_95 = XP_Lower,
              Fixed_Put_Upper_95 = XP_Upper
  )
}

# RUN ALL THREE METHODS
################
results_list <- list()

counter <- 1

for (M in M_values) {
  cat("Running M =", M, "\n")
  
  # Common random numbers
  Z <- matrix( rnorm(n_paths * M), nrow = n_paths, ncol = M )
  
  # Exact
  results_list[[counter]] <- simulate_lookback("Exact", M, Z)
  counter <- counter + 1
  
  # Euler-Maruyama
  results_list[[counter]] <- simulate_lookback("Euler-Maruyama", M, Z)
  counter <- counter + 1
  
  # Milstein
  results_list[[counter]] <- simulate_lookback("Milstein", M, Z)
  counter <- counter + 1
  
  rm(Z)
}

# COMBINE RESULTS
########
results_lookback <- do.call( rbind,  results_list )
rownames(results_lookback) <- NULL

# ROUND FOR PRESENTATION
##########
results_lookback[,-c(1,2)] <- round(results_lookback[,-c(1,2)], 4)

print(results_lookback)


####################################################################
###########################################################
# Section 4.8
#################
############################################
# Variance Reduction Techniques
#############################################
# Antithetic Variates Monte Carlo - European Options
# ===========================

S0    <- 80
K     <- 90
r     <- 0.05
q     <- 0
sigma <- 0.20
T     <- 1

set.seed(123)

N <- 100000

# Generate standard normal random variables
Z <- rnorm(N)

# Antithetic random variables
Z_anti <- -Z

# Simulate terminal stock prices
S_T <- S0 * exp((r - q - 0.5 * sigma^2) * T + sigma * sqrt(T) * Z)

S_T_anti <- S0 * exp((r - q - 0.5 * sigma^2) * T + sigma * sqrt(T) * Z_anti)

# ------------------------------------------------------------
# European Call Payoffs
# ------------------------------------------------------------

call_payoff <- pmax(S_T - K, 0)
call_payoff_anti <- pmax(S_T_anti - K, 0)

# Average antithetic payoffs
call_anti_payoff <- (call_payoff + call_payoff_anti) / 2

# Monte Carlo price
call_anti_MC <- exp(-r * T) * mean(call_anti_payoff)

# Standard error
call_anti_SE <- exp(-r * T) * sd(call_anti_payoff) / sqrt(N)

# ------------------------------------------------------------
# European Put Payoffs
# ------------------------------------------------------------

put_payoff <- pmax(K - S_T, 0)
put_payoff_anti <- pmax(K - S_T_anti, 0)

# Average antithetic payoffs
put_anti_payoff <- (put_payoff + put_payoff_anti) / 2

# Monte Carlo price
put_anti_MC <- exp(-r * T) * mean(put_anti_payoff)

# Standard error
put_anti_SE <- exp(-r * T) * sd(put_anti_payoff) / sqrt(N)

# ------------------------------------------------------------
# Results
# ------------------------------------------------------------

call_anti_MC   #  4.1619
call_anti_SE   #  0.01703826

put_anti_MC    # 9.773469
put_anti_SE    # 0.01012911
###########################################

# ============================================================
# Standard Monte Carlo vs Antithetic Variates
# European Options
# ============================================================

# Standard Monte Carlo prices
call_MC <- exp(-r * T) * mean(call_payoff)
put_MC  <- exp(-r * T) * mean(put_payoff)

# Standard Monte Carlo standard errors
call_SE <- exp(-r * T) * sd(call_payoff) / sqrt(N)
put_SE  <- exp(-r * T) * sd(put_payoff) / sqrt(N)

# Antithetic Monte Carlo prices
call_anti_MC <- exp(-r * T) * mean(call_anti_payoff)
put_anti_MC  <- exp(-r * T) * mean(put_anti_payoff)

# Antithetic standard errors
call_anti_SE <- exp(-r * T) * sd(call_anti_payoff) / sqrt(N)
put_anti_SE  <- exp(-r * T) * sd(put_anti_payoff) / sqrt(N)

# Results
results <- data.frame(
  Option = c("European Call", "European Put"),
  Standard_MC = c(call_MC, put_MC),
  Antithetic_MC = c(call_anti_MC, put_anti_MC),
  Standard_SE = c(call_SE, put_SE),
  Antithetic_SE = c(call_anti_SE, put_anti_SE)
)

print(results)

#   Option           Standard_MC  Antithetic_MC Standard_SE  Antithetic_SE
# 1 European Call   4.168328    4.161900        0.02743396    0.01703826
# 2 European Put    9.764681       9.773469     0.03228162    0.01012911
#############################################################################
#######################################################################

# ============================================================
# Control Variate: W_T -European Options
# ============================================================
# Model Parameters unchanged
set.seed(123)

N <- 100000

# Simulate W_T ~ N(0, T)
Z <- rnorm(N)
W_T <- sqrt(T) * Z


# Simulate terminal stock price
S_T <- S0 * exp( (r - q - 0.5 * sigma^2) * T + sigma * W_T )

# ------------------------------------------------------------
# European call and put payoffs
# ------------------------------------------------------------

call_payoff <- pmax(S_T - K, 0)
put_payoff  <- pmax(K - S_T, 0)

# ------------------------------------------------------------
# Estimate optimal control coefficients
# Covariance(payoff, W_T) / Var(W_T)
# ------------------------------------------------------------

b_call <- cov(call_payoff, W_T) / var(W_T)
b_put  <- cov(put_payoff, W_T) / var(W_T)

# ------------------------------------------------------------
# Control variate adjusted payoffs
# E[W_T] = 0
# ------------------------------------------------------------

call_cv_payoff <- call_payoff - b_call * W_T
put_cv_payoff  <- put_payoff - b_put * W_T

# ------------------------------------------------------------
# Control variate Monte Carlo prices
# ------------------------------------------------------------

call_CV <- exp(-r * T) * mean(call_cv_payoff)
put_CV  <- exp(-r * T) * mean(put_cv_payoff)

# ------------------------------------------------------------
# Standard errors
# ------------------------------------------------------------

call_CV_SE <- exp(-r * T) * sd(call_cv_payoff) / sqrt(N)

put_CV_SE <- exp(-r * T) * sd(put_cv_payoff) / sqrt(N)

# Results
call_CV      #   4.161998
call_CV_SE   #   0.01824446

put_CV       #  9.773975
put_CV_SE    # 0.01171456


# ============================================================
# Control Variate: S_T
# European Options
# ============================================================
# Parameters unchanged
set.seed(123)

N <- 100000

# Simulate terminal stock price
Z <- rnorm(N)

S_T <- S0 * exp((r - q - 0.5 * sigma^2) * T +sigma * sqrt(T) * Z )

# Expected value of S_T under the risk-neutral measure
E_ST <- S0 * exp((r - q) * T)

# ------------------------------------------------------------
# European option payoffs
# ------------------------------------------------------------

call_payoff <- pmax(S_T - K, 0)
put_payoff  <- pmax(K - S_T, 0)

# ------------------------------------------------------------
# Optimal control coefficients
# ------------------------------------------------------------

b_call <- cov(call_payoff, S_T) / var(S_T)
b_put  <- cov(put_payoff, S_T) / var(S_T)

# ------------------------------------------------------------
# Control variate adjusted payoffs
# ------------------------------------------------------------

call_cv_payoff <- call_payoff - b_call * (S_T - E_ST)

put_cv_payoff <- put_payoff - b_put * (S_T - E_ST)

# ------------------------------------------------------------
# Control variate Monte Carlo prices
# ------------------------------------------------------------

call_ST_CV <- exp(-r * T) * mean(call_cv_payoff)
put_ST_CV  <- exp(-r * T) * mean(put_cv_payoff)

# ------------------------------------------------------------
# Standard errors
# ------------------------------------------------------------

call_ST_CV_SE <- exp(-r * T) * sd(call_cv_payoff) / sqrt(N)

put_ST_CV_SE <- exp(-r * T) * sd(put_cv_payoff) / sqrt(N)

# Results
call_ST_CV     #   4.161974
call_ST_CV_SE  #  0.01539924

put_ST_CV      # 9.772622
put_ST_CV_SE  #  0.01539924

###########################################################
# ANALYTICAL VS MONTE CARLO COMPARISON of LB OPTIONS
############################################################

# ANALYTICAL VALUATION OF FLOATING-STRIKE LOOKBACK OPTIONS (Black-Scholes framework)
#======================================================

#  Model and option parameters 
# ---------------------------------
K=90
S0 <- 80          # Initial stock price
r <- 0.05         # Risk-free interest rate
q <- 0            # Dividend yield
sigma <- 0.20     # Volatility
T <- 1            # Time to maturity
t <- 0            # Current time

tau <- T - t      # Time remaining to maturity

# At inception:
# Running minimum = initial stock price
# Running maximum = initial stock price

m0 <- S0          # Running minimum
M0 <- S0          # Running maximum

N <- pnorm   # Standard normal cumulative distribution function


# ======================
#  FLOATING-STRIKE LOOKBACK CALL
# ==========================

floating_call <- function(S, m, r, q, sigma, tau) {
  
  # d1 and d2
  d1 <- (log(S / m) + (r - q + 0.5 * sigma^2) * tau) / (sigma * sqrt(tau))
  
  d2 <- d1 - sigma * sqrt(tau)
  

  # Special case: r = q

  if (abs(r - q) < 1e-10) {
    price <- S * exp(-r * tau) * N(d1) - m * exp(-r * tau) * N(d2) + S * exp(-r * tau) *
      sigma * sqrt(tau) *( dnorm(d1) - d1 * N(-d1) )
    
  } else {# General case: r != q
    
    a3 <- -d1 + (2 * (r - q))/sigma * sqrt(tau)
    
    price <- S * exp(-q * tau) * N(d1) - m * exp(-r * tau) * N(d2) + (S * sigma^2)/ (2 * (r - q)) * ( exp(-r * tau)*(S/m)^(-(2*(r-q))/sigma^2) * N(a3)-exp(-q*tau)*N(-d1))
                                                                                                   
  }

  return(price)
}

# =================
#FLOATING-STRIKE LOOKBACK PUT
# =================

floating_put <- function(S, M, r, q, sigma, tau) {
  # e1 and e2
  e1 <- (log(S / M) +(r - q + 0.5 * sigma^2) * tau) /  (sigma * sqrt(tau))
  
  e2 <- e1 - sigma * sqrt(tau)
  
  # ------------
  # Special case: r = q
  # ------------
  
  if (abs(r - q) < 1e-10) {
    price <- M * exp(-r * tau) * N(-e2) -  S * exp(-r * tau) * N(-e1) + S * exp(-r * tau) * sigma * sqrt(tau) * ( dnorm(e1) + e1 * N(e1) )
    
    } else {
    
    # General case: r != q
      
    e3 <- e1 - (2 * (r - q))/sigma * sqrt(tau) 
    
    price <- M * exp(-r * tau) * N(-e2) - S * exp(-q * tau) * N(-e1) + (S * sigma^2) / (2 * (r - q)) * (-exp(-r*tau)*(S/M)^(-2*(r-q)/sigma^2)*N(e1-(2*(r-q)/sigma)* sqrt(tau))+exp(-q*tau)*N(e1))
    

    }
 
  
  return(price)
}
#=======================
#  Calculate analytical prices
# ===============================

LB_Floating_Call <- floating_call(S = S0, m = m0, r = r, q = q, sigma = sigma, tau = tau )

LB_Floating_Put <- floating_put( S = S0, M = M0, r = r, q = q, sigma = sigma, tau = tau )

LB_Floating_Call  # 13.77344
LB_Floating_Put   # 11.4325

#############################################################

# ANALYTICAL VALUATION OF FIXED-STRIKE LB OPTIONS
############################################################

#  Model and option parameters remain unchanged

# =========
# FIXED-STRIKE LOOKBACK CALL (q !=r)
# ===========
K=90
fixed_call <- function(S, K, M, r, q, sigma, tau) {

  # Case 1: K > M
 #----------------------
  
  if (K > M) {
    
    d <- (log(S / K) + (r - q + 0.5 * sigma^2) * tau) / (sigma * sqrt(tau))
    
    d4 <- d - sigma * sqrt(tau)
    
    d5 <- d - (2 * (r - q))/ sigma * sqrt(tau) 
    
    price <- S * exp(-q * tau) * N(d) - K * exp(-r * tau) * N(d4) + sigma^2/(2 *(r - q)) *S* (-exp(-r * tau)*(S / K)^(-2 * (r - q) / sigma^2) * N(d5) + exp(-q * tau) * N(d) )
    
  } else {
    
    # ---------
    # Case 2: K <= M
    # -----
    
    d1 <- (log(S / m) + (r - q + 0.5 * sigma^2) * tau) / (sigma * sqrt(tau))
    
    d2 <- d1 - sigma * sqrt(tau)
    
    b <-  d1 - (2 * (r - q))/ sigma * sqrt(tau) 
    
    price <- exp(-r * tau) * (M - K) + S * exp(-q * tau) * N(d1) - M * exp(-r * tau) * N(d2)   + sigma^2/ (2 *(r - q))  * S * ( -exp(-r*tau)*(S / M)^(-2 * (r - q) / sigma^2) * N(b) +
          exp(- q * tau) * N(d1)  )
  }
  
  return(price)
}

# =================
# FIXED-STRIKE LOOKBACK PUT
# ==============

fixed_put <- function(S, K, m, r, q, sigma, tau) {
  
  # Case 1: K < m
  # --------------------------------------------------------
  
  if (K < m) {
    
    d<- (log(S / K) + (r - q + 0.5 * sigma^2) * tau) / (sigma * sqrt(tau))
    
    d4 <- d - sigma * sqrt(tau)
    
    d3<- d1 - (2 * (r - q))/ sigma * sqrt(tau) 
    
    price <- -exp(-q * tau)*S * N(-d) + exp(-r*tau)* K * N * (-d + sigma*sqrt(tau)) +
      sigma^2/(2 * (r - q)) * S * (exp(-r * tau)  * ((S / K)^(-2 * (r - q) / sigma^2) *  
      N(-d3) - exp(- q * tau) * N(-d) ))
  } else {
    
    # -------
    # Case 2: K >= m
    # -------
    
    d1 <- (log(S / m) + (r - q + 0.5 * sigma^2) * tau) / (sigma * sqrt(tau))
    
    d2 <- d1 - sigma * sqrt(tau)
    
    b<- d1 -( 2 * (r - q)) / sigma * sqrt(tau)
    
    price <- exp(-r * tau)* (K - m) - exp(-q * tau)* S * N(-d1) + exp(-r * tau) *m *N(-d1+ sigma*sqrt(tau)) + sigma^2/(2*(r-q)) * S*(exp(-r * tau)* (S / m)^((-2 * (r - q)) / sigma^2) * N(-d1 + (2*(r-q)/sigma)*sqrt(tau))- exp( - q * tau) * N(-d1) )
  }
  
  return(price)
}

# =============
# Calculate analytical fixed-strike prices
# ==============

LB_Fixed_Call <- fixed_call( S = S0, K = K, M = M0, r = r, q = q, sigma = sigma, tau = tau )
LB_Fixed_Put <- fixed_put( S = S0, K = K,  m = m0, r = r, q = q, sigma = sigma, tau = tau )

LB_Fixed_Call  # 7.7458
LB_Fixed_Put   # 19.3841

####################################################################################
############################################################
# LOOKBACK OPTIONS: ANALYTICAL VS MONTE CARLO COMPARISON
# ==========================================================
# 1. MONTE CARLO RESULTS
# =======================
MC_Floating_Call <- 13.3100 
SE_Floating_Call <- 0.0368

MC_Floating_Put <- 10.7734  
SE_Floating_Put <- 0.0249

MC_Fixed_Call <- 7.3506
SE_Fixed_Call <- 0.0332

MC_Fixed_Put <- 18.9054
SE_Fixed_Put <- 0.0227


# 2. ANALYTICAL VALUES
# ====================
# These should come from your closed-form functions.

Analytical_Floating_Call <- LB_Floating_Call
Analytical_Floating_Put  <- LB_Floating_Put

Analytical_Fixed_Call <- LB_Fixed_Call
Analytical_Fixed_Put  <- LB_Fixed_Put


# 3. 95% CONFIDENCE INTERVALS
# ==============================

Floating_Call_Lower_95 <- MC_Floating_Call - 1.96 * SE_Floating_Call
Floating_Call_Upper_95 <- MC_Floating_Call + 1.96 * SE_Floating_Call


Floating_Put_Lower_95 <- MC_Floating_Put - 1.96 * SE_Floating_Put
Floating_Put_Upper_95 <-  MC_Floating_Put + 1.96 * SE_Floating_Put


Fixed_Call_Lower_95 <- MC_Fixed_Call - 1.96 * SE_Fixed_Call
Fixed_Call_Upper_95 <- MC_Fixed_Call + 1.96 * SE_Fixed_Call


Fixed_Put_Lower_95 <-  MC_Fixed_Put - 1.96 * SE_Fixed_Put
Fixed_Put_Upper_95 <- MC_Fixed_Put + 1.96 * SE_Fixed_Put


# 4. ABSOLUTE ERRORS
# =====================

Absolute_Error_Floating_Call <- abs(MC_Floating_Call - Analytical_Floating_Call)

Absolute_Error_Floating_Put <- abs(MC_Floating_Put - Analytical_Floating_Put)

Absolute_Error_Fixed_Call <- abs(MC_Fixed_Call - Analytical_Fixed_Call)

Absolute_Error_Fixed_Put <-  abs(MC_Fixed_Put - Analytical_Fixed_Put)


# 5. PERCENTAGE ERRORS
# ================

Percentage_Error_Floating_Call <- Absolute_Error_Floating_Call /Analytical_Floating_Call * 100

Percentage_Error_Floating_Put <- Absolute_Error_Floating_Put / Analytical_Floating_Put * 100

Percentage_Error_Fixed_Call <- Absolute_Error_Fixed_Call / Analytical_Fixed_Call * 100

Percentage_Error_Fixed_Put <- Absolute_Error_Fixed_Put / Analytical_Fixed_Put * 100

#=================================
# 6. CREATE COMPARISON TABLE
# ===============================

Lookback_Comparison <- data.frame(
  
  Option = c( "Floating-Strike Call", "Floating-Strike Put", "Fixed-Strike Call", "Fixed-Strike Put" ),
  
  Analytical_Value = c( Analytical_Floating_Call, Analytical_Floating_Put,  Analytical_Fixed_Call, Analytical_Fixed_Put),
  
  Monte_Carlo_Value = c( MC_Floating_Call, MC_Floating_Put,  MC_Fixed_Call, MC_Fixed_Put ),
  
  Standard_Error = c( SE_Floating_Call, SE_Floating_Put, SE_Fixed_Call, SE_Fixed_Put ),
  
  Lower_95 = c( Floating_Call_Lower_95, Floating_Put_Lower_95, Fixed_Call_Lower_95, Fixed_Put_Lower_95 ),
  
  Upper_95 = c( Floating_Call_Upper_95, Floating_Put_Upper_95, Fixed_Call_Upper_95, Fixed_Put_Upper_95 ),
  
  Absolute_Error = c( Absolute_Error_Floating_Call, Absolute_Error_Floating_Put, Absolute_Error_Fixed_Call, Absolute_Error_Fixed_Put ),
  
  Percentage_Error = c( Percentage_Error_Floating_Call, Percentage_Error_Floating_Put, Percentage_Error_Fixed_Call,
 Percentage_Error_Fixed_Put ))


# ==========
# 7. ROUND RESULTS FOR PRESENTATION
# =============

Lookback_Comparison$Analytical_Value <- round(Lookback_Comparison$Analytical_Value, 4)

Lookback_Comparison$Monte_Carlo_Value <- round(Lookback_Comparison$Monte_Carlo_Value, 4)

Lookback_Comparison$Standard_Error <- round(Lookback_Comparison$Standard_Error, 4)

Lookback_Comparison$Lower_95 <- round(Lookback_Comparison$Lower_95, 4)

Lookback_Comparison$Upper_95 <- round(Lookback_Comparison$Upper_95, 4)

Lookback_Comparison$Absolute_Error <- round(Lookback_Comparison$Absolute_Error, 4)

Lookback_Comparison$Percentage_Error <- round(Lookback_Comparison$Percentage_Error, 4)


# =======
# 8. DISPLAY TABLE
# ========

print(Lookback_Comparison)

#      Option                  Analytical_Value  Monte_Carlo_Value Standard_Error
# 1 Floating-Strike Call          13.7734           13.3100         0.0368
# 2  Floating-Strike Put          11.4325           10.7734         0.0249
# 3    Fixed-Strike Call           7.7458            7.3506         0.0332
# 4     Fixed-Strike Put          19.3841           18.9054         0.0227
#    Lower_95 Upper_95       Absolute_Error Percentage_Error
# 1  13.2379  13.3821         0.4634           3.3647
# 2  10.7246  10.8222         0.6591           5.7648
# 3   7.2855   7.4157         0.3952           5.1020
# 4  18.8609  18.9499         0.4787           2.4695


# ============================================================
## VARIANCE, STANDARD ERROR, AND VARIANCE REDUCTION CALCULATIONS FOR EUROPEAN OPTIONS
#===================================================================

# Results from your Monte Carlo simulations
results_variance <- data.frame(
  Option = c("European Call", "European Put"),
  
  Standard_SE = c(0.02743, 0.03228),
  
  Antithetic_SE = c(0.01704, 0.01013),
  
  WT_CV_SE = c(0.01824, 0.01171),
  
  ST_CV_SE = c(0.01540, 0.01540) )

# Number of simulations
N <- 100000


# ------------------------------------------------------------
# CALCULATE VARIANCES
# ------------------------------------------------------------
#
# SE = sqrt(Variance / N)
# Therefore:  Variance = SE^2 * N
#

results_variance$Standard_Variance <- results_variance$Standard_SE^2 * N

results_variance$Antithetic_Variance <- results_variance$Antithetic_SE^2 * N

results_variance$WT_CV_Variance <-  results_variance$WT_CV_SE^2 * N

results_variance$ST_CV_Variance <- results_variance$ST_CV_SE^2 * N


# ------------------------------------------------------------
# ABSOLUTE VARIANCE REDUCTION
# ------------------------------------------------------------
#
# Variance Reduction = Standard Variance - Reduced Variance
#

results_variance$Antithetic_Variance_Reduction <-
  results_variance$Standard_Variance -
  results_variance$Antithetic_Variance

results_variance$WT_CV_Variance_Reduction <-
  results_variance$Standard_Variance -
  results_variance$WT_CV_Variance

results_variance$ST_CV_Variance_Reduction <-
  results_variance$Standard_Variance -
  results_variance$ST_CV_Variance


# ------------------------------------------------------------
# PERCENTAGE VARIANCE REDUCTION
# ------------------------------------------------------------
#
# % Variance Reduction = 100 * (1 - Reduced Variance / Standard Variance)
#
results_variance$Antithetic_VR_Percent <-
  100 * ( 1 -
            results_variance$Antithetic_Variance /
            results_variance$Standard_Variance
  )

results_variance$WT_CV_VR_Percent <-
  100 * ( 1 -
            results_variance$WT_CV_Variance /
            results_variance$Standard_Variance
  )

results_variance$ST_CV_VR_Percent <-
  100 * (1 -
           results_variance$ST_CV_Variance /
           results_variance$Standard_Variance
  )


# ------------------------------------------------------------
# STANDARD ERROR REDUCTION
# ------------------------------------------------------------
#
# % SE Reduction = 100 * (1 - Reduced SE / Standard SE)
#

results_variance$Antithetic_SE_Reduction_Percent <-
  100 * ( 1 -
            results_variance$Antithetic_SE /
            results_variance$Standard_SE
  )

results_variance$WT_CV_SE_Reduction_Percent <-
  100 * (  1 -
             results_variance$WT_CV_SE /
             results_variance$Standard_SE
  )

results_variance$ST_CV_SE_Reduction_Percent <-
  100 * ( 1 -
            results_variance$ST_CV_SE /
            results_variance$Standard_SE
  )
# ------------------------------------------------------------
# ROUND RESULTS
# ------------------------------------------------------------

results_variance[-1] <- round(results_variance[-1], 5)

print(results_variance)

#  Option            Standard_SE Antithetic_SE WT_CV_SE ST_CV_SE  Standard_Variance
# 1 European Call     0.02743       0.01704   0.01824   0.0154          75.24049
# 2  European Put     0.03228       0.01013   0.01171   0.0154         104.19984
#           Antithetic_Variance WT_CV_Variance ST_CV_Variance
# 1            29.03616       33.26976         23.716
# 2            10.26169       13.71241         23.716
#                 Antithetic_Variance_Reduction  WT_CV_Variance_Reduction
# 1                      46.20433                 41.97073
# 2                      93.93815                 90.48743
#               ST_CV_Variance_Reduction Antithetic_VR_Percent WT_CV_VR_Percent
# 1                 51.52449              61.40886             55.78211
# 2                 80.48384              90.15191             86.84028
#        ST_CV_VR_Percent           Antithetic_SE_Reduction_Percent   WT_CV_SE_Reduction_Percent
# 1         68.47974                        37.87824                   33.50346
# 2         77.23989                        68.61834                   63.72367
#       ST_CV_SE_Reduction_Percent
# 1                   43.85709
# 2                   52.29244
###############################################################################################
####################################################################################

# ============================================================
# ANALYTICAL VARIANCE, VARIANCE REDUCTION
# AND MONTE CARLO COMPARISON OF EUROPEAN OPTIONS
# ============================================================

# Model Parameters:
S0    <- 80
K     <- 90
r     <- 0.05
q     <- 0
sigma <- 0.20
T     <- 1
N     <- 100000

Phi <- pnorm
phi <- dnorm


# ============================================================
# BLACK-SCHOLES QUANTITIES
# ============================================================

d1 <- (
  log(S0 / K) +
    (r - q + 0.5 * sigma^2) * T
) / (sigma * sqrt(T))

d2 <- d1 - sigma * sqrt(T)


# Black-Scholes prices

Call_BS <-
  S0 * exp(-q * T) * Phi(d1) -
  K * exp(-r * T) * Phi(d2)

Put_BS <-
  K * exp(-r * T) * Phi(-d2) -
  S0 * exp(-q * T) * Phi(-d1)


# ============================================================
# ANALYTICAL STANDARD MONTE CARLO VARIANCE
# ============================================================

# ------------------------------------------------------------
# European Call
# ------------------------------------------------------------

A0 <- Phi(d2)

A1 <-
  S0 *
  exp((r - q) * T) *
  Phi(d1)

A2 <-
  S0^2 *
  exp(2 * (r - q) * T + sigma^2 * T) *
  Phi(d1 + sigma * sqrt(T))


Call_E2 <-
  exp(-2 * r * T) *
  (
    A2 -
      2 * K * A1 +
      K^2 * A0
  )

Call_Var <-
  Call_E2 - Call_BS^2


# ------------------------------------------------------------
# European Put
# ------------------------------------------------------------

B0 <- Phi(-d2)

B1 <-
  S0 *
  exp((r - q) * T) *
  Phi(-d1)

B2 <-
  S0^2 *
  exp(2 * (r - q) * T + sigma^2 * T) *
  Phi(-d1 - sigma * sqrt(T))


Put_E2 <-
  exp(-2 * r * T) *
  (
    K^2 * B0 -
      2 * K * B1 +
      B2
  )

Put_Var <-
  Put_E2 - Put_BS^2


# Standard errors

Call_SE <- sqrt(Call_Var / N)
Put_SE  <- sqrt(Put_Var / N)


# ============================================================
# ANTITHETIC VARIATES
# ============================================================

mu <- (r - q - 0.5 * sigma^2) * T
b  <- sigma * sqrt(T)


# ------------------------------------------------------------
# Call
# ------------------------------------------------------------

Call_AV_integrand <- function(z) {
  
  ST_z <- S0 * exp(mu + b * z)
  
  ST_minus_z <- S0 * exp(mu - b * z)
  
  payoff_z <-
    exp(-r * T) *
    pmax(ST_z - K, 0)
  
  payoff_minus_z <-
    exp(-r * T) *
    pmax(ST_minus_z - K, 0)
  
  payoff_z *
    payoff_minus_z *
    dnorm(z)
}


Call_AV_product <-
  integrate(
    Call_AV_integrand,
    lower = -10,
    upper = 10
  )$value


Call_AV_Cov <-
  Call_AV_product -
  Call_BS^2


Call_AV_Var <-
  (
    Call_Var +
      Call_AV_Cov
  ) / 2


Call_AV_SE <-
  sqrt(Call_AV_Var / N)


# ------------------------------------------------------------
# Put
# ------------------------------------------------------------

Put_AV_integrand <- function(z) {
  
  ST_z <- S0 * exp(mu + b * z)
  
  ST_minus_z <- S0 * exp(mu - b * z)
  
  payoff_z <-
    exp(-r * T) *
    pmax(K - ST_z, 0)
  
  payoff_minus_z <-
    exp(-r * T) *
    pmax(K - ST_minus_z, 0)
  
  payoff_z *
    payoff_minus_z *
    dnorm(z)
}


Put_AV_product <-
  integrate(
    Put_AV_integrand,
    lower = -10,
    upper = 10
  )$value


Put_AV_Cov <-
  Put_AV_product -
  Put_BS^2


Put_AV_Var <-
  (
    Put_Var +
      Put_AV_Cov
  ) / 2


Put_AV_SE <-
  sqrt(Put_AV_Var / N)


# ============================================================
# W_T CONTROL VARIATE
# ============================================================

# W_T ~ N(0,T)

Var_W <- T


# ------------------------------------------------------------
# Call covariance with W_T
# ------------------------------------------------------------

Call_Cov_W <-
  exp(-r * T) *
  sqrt(T) *
  (
    S0 *
      exp((r - q) * T) *
      (
        sigma * sqrt(T) * Phi(d1) +
          phi(d1)
      )
    -
      K * phi(d2)
  )


Call_W_Var <-
  Call_Var -
  Call_Cov_W^2 / Var_W


Call_W_SE <-
  sqrt(Call_W_Var / N)


Call_rho_W <-
  Call_Cov_W /
  sqrt(Call_Var * Var_W)


# ------------------------------------------------------------
# Put covariance with W_T
# ------------------------------------------------------------

Put_Cov_W <-
  exp(-r * T) *
  sqrt(T) *
  (
    -K * phi(d2)
    -
      S0 *
      exp((r - q) * T) *
      (
        sigma * sqrt(T) * Phi(-d1) -
          phi(d1)
      )
  )


Put_W_Var <-
  Put_Var -
  Put_Cov_W^2 / Var_W


Put_W_SE <-
  sqrt(Put_W_Var / N)


Put_rho_W <-
  Put_Cov_W /
  sqrt(Put_Var * Var_W)


# ============================================================
# S_T CONTROL VARIATE
# ============================================================

# Control variable:
# Y = exp(-rT) S_T
#
# E[Y] = S0 exp(-qT)

E_ST_Control <-
  S0 * exp(-q * T)


Var_ST_Control <-
  S0^2 *
  exp(-2 * q * T) *
  (
    exp(sigma^2 * T) - 1
  )


# ------------------------------------------------------------
# Call covariance with discounted S_T
# ------------------------------------------------------------

Call_Cov_ST <-
  exp(-2 * r * T) *
  (
    A2 -
      K * A1
  ) -
  Call_BS * E_ST_Control


Call_ST_Var <-
  Call_Var -
  Call_Cov_ST^2 /
  Var_ST_Control


Call_ST_SE <-
  sqrt(Call_ST_Var / N)


Call_rho_ST <-
  Call_Cov_ST /
  sqrt(
    Call_Var *
      Var_ST_Control
  )


# ------------------------------------------------------------
# Put covariance with discounted S_T
# ------------------------------------------------------------

Put_Cov_ST <-
  exp(-2 * r * T) *
  (
    K * B1 -
      B2
  ) -
  Put_BS * E_ST_Control


Put_ST_Var <-
  Put_Var -
  Put_Cov_ST^2 /
  Var_ST_Control


Put_ST_SE <-
  sqrt(Put_ST_Var / N)


Put_rho_ST <-
  Put_Cov_ST /
  sqrt(
    Put_Var *
      Var_ST_Control
  )


# ============================================================
# ANALYTICAL VARIANCE REDUCTION PERCENTAGES
# ============================================================

# Antithetic

Call_AV_VR <-
  100 *
  (
    1 -
      Call_AV_Var /
      Call_Var
  )

Put_AV_VR <-
  100 *
  (
    1 -
      Put_AV_Var /
      Put_Var
  )


# W_T control variate

Call_W_VR <-
  100 *
  (
    1 -
      Call_W_Var /
      Call_Var
  )

Put_W_VR <-
  100 *
  (
    1 -
      Put_W_Var /
      Put_Var
  )


# S_T control variate

Call_ST_VR <-
  100 *
  (
    1 -
      Call_ST_Var /
      Call_Var
  )

Put_ST_VR <-
  100 *
  (
    1 -
      Put_ST_Var /
      Put_Var
  )


# ============================================================
# MONTE CARLO RESULTS
# ============================================================

# Values obtained from your simulations

MC_Call_Var <- c(
  Standard   = 75.24049,
  Antithetic = 29.03616,
  WT_CV      = 33.26976,
  ST_CV      = 23.71600
)

MC_Put_Var <- c(
  Standard   = 104.19984,
  Antithetic = 10.26169,
  WT_CV      = 13.71241,
  ST_CV      = 23.71600
)


MC_Call_SE <- c(
  Standard   = 0.02743,
  Antithetic = 0.01704,
  WT_CV      = 0.01824,
  ST_CV      = 0.01540
)

MC_Put_SE <- c(
  Standard   = 0.03228,
  Antithetic = 0.01013,
  WT_CV      = 0.01171,
  ST_CV      = 0.01540
)


# ============================================================
# MONTE CARLO VARIANCE REDUCTION
# ============================================================

MC_Call_VR <- c(
  Standard   = 0,
  Antithetic =
    100 * (
      1 -
        MC_Call_Var["Antithetic"] /
        MC_Call_Var["Standard"]
    ),
  
  WT_CV =
    100 * (
      1 -
        MC_Call_Var["WT_CV"] /
        MC_Call_Var["Standard"]
    ),
  
  ST_CV =
    100 * (
      1 -
        MC_Call_Var["ST_CV"] /
        MC_Call_Var["Standard"]
    )
)


MC_Put_VR <- c(
  Standard   = 0,
  Antithetic =
    100 * (
      1 -
        MC_Put_Var["Antithetic"] /
        MC_Put_Var["Standard"]
    ),
  
  WT_CV =
    100 * (
      1 -
        MC_Put_Var["WT_CV"] /
        MC_Put_Var["Standard"]
    ),
  
  ST_CV =
    100 * (
      1 -
        MC_Put_Var["ST_CV"] /
        MC_Put_Var["Standard"]
    )
)


# ============================================================
# ANALYTICAL RESULTS
# ============================================================

Analytical_Call_Var <- c(
  Standard   = Call_Var,
  Antithetic = Call_AV_Var,
  WT_CV      = Call_W_Var,
  ST_CV      = Call_ST_Var
)


Analytical_Put_Var <- c(
  Standard   = Put_Var,
  Antithetic = Put_AV_Var,
  WT_CV      = Put_W_Var,
  ST_CV      = Put_ST_Var
)


Analytical_Call_SE <- c(
  Standard   = Call_SE,
  Antithetic = Call_AV_SE,
  WT_CV      = Call_W_SE,
  ST_CV      = Call_ST_SE
)


Analytical_Put_SE <- c(
  Standard   = Put_SE,
  Antithetic = Put_AV_SE,
  WT_CV      = Put_W_SE,
  ST_CV      = Put_ST_SE
)


Analytical_Call_VR <- c(
  Standard   = 0,
  Antithetic = Call_AV_VR,
  WT_CV      = Call_W_VR,
  ST_CV      = Call_ST_VR
)


Analytical_Put_VR <- c(
  Standard   = 0,
  Antithetic = Put_AV_VR,
  WT_CV      = Put_W_VR,
  ST_CV      = Put_ST_VR
)


# ============================================================
# ANALYTICAL VS MONTE CARLO COMPARISON
# ============================================================

comparison_table <- data.frame(
  
  Option = c(
    rep("European Call", 4),
    rep("European Put", 4)
  ),
  
  Method = rep(
    c(
      "Standard MC",
      "Antithetic",
      "W_T CV",
      "S_T CV"
    ),
    2
  ),
  
  Analytical_Variance = c(
    Analytical_Call_Var,
    Analytical_Put_Var
  ),
  
  MC_Variance = c(
    MC_Call_Var,
    MC_Put_Var
  ),
  
  Analytical_SE = c(
    Analytical_Call_SE,
    Analytical_Put_SE
  ),
  
  MC_SE = c(
    MC_Call_SE,
    MC_Put_SE
  ),
  
  Analytical_VR_Percent = c(
    Analytical_Call_VR,
    Analytical_Put_VR
  ),
  
  MC_VR_Percent = c(
    MC_Call_VR,
    MC_Put_VR
  )
)


# ============================================================
#  DIFFERENCES BETWEEN ANALYTICAL AND MC RESULTS
# ============================================================

comparison_table$Variance_Absolute_Difference <-
  abs(
    comparison_table$Analytical_Variance -
      comparison_table$MC_Variance
  )


comparison_table$Variance_Percent_Difference <-
  100 *
  comparison_table$Variance_Absolute_Difference /
  comparison_table$Analytical_Variance


comparison_table$SE_Absolute_Difference <-
  abs(
    comparison_table$Analytical_SE -
      comparison_table$MC_SE
  )


# ============================================================
# ROUND RESULTS
# ============================================================

comparison_table$Analytical_Variance <-
  round(
    comparison_table$Analytical_Variance,
    5
  )

comparison_table$MC_Variance <-
  round(
    comparison_table$MC_Variance,
    5
  )

comparison_table$Analytical_SE <-
  round(
    comparison_table$Analytical_SE,
    5
  )

comparison_table$MC_SE <-
  round(
    comparison_table$MC_SE,
    5
  )

comparison_table$Analytical_VR_Percent <-
  round(
    comparison_table$Analytical_VR_Percent,
    5
  )

comparison_table$MC_VR_Percent <-
  round(
    comparison_table$MC_VR_Percent,
    5
  )

comparison_table$Variance_Absolute_Difference <-
  round(
    comparison_table$Variance_Absolute_Difference,
    5
  )

comparison_table$Variance_Percent_Difference <-
  round(
    comparison_table$Variance_Percent_Difference,
    5
  )

comparison_table$SE_Absolute_Difference <-
  round(
    comparison_table$SE_Absolute_Difference,
    5
  )


# ============================================================
# DISPLAY FULL COMPARISON TABLE
# ============================================================

print(comparison_table)

#  Option      Method Analytical_Variance MC_Variance Analytical_SE
# 1 European Call Standard MC            75.55152    75.24049       0.02749
# 2 European Call  Antithetic            29.10468    29.03616       0.01706
# 3 European Call      W_T CV            33.43958    33.26976       0.01829
# 4 European Call      S_T CV            23.80329    23.71600       0.01543
# 5  European Put Standard MC           104.22323   104.19984       0.03228
# 6  European Put  Antithetic            10.26533    10.26169       0.01013
# 7  European Put      W_T CV            13.77117    13.71241       0.01174
# 8  European Put      S_T CV            23.80329    23.71600       0.01543
# MC_SE Analytical_VR_Percent MC_VR_Percent Variance_Absolute_Difference
# 1 0.02743               0.00000       0.00000                      0.31103
# 2 0.01704              61.47704      61.40886                      0.06852
# 3 0.01824              55.73936      55.78211                      0.16982
# 4 0.01540              68.49397      68.47974                      0.08729
# 5 0.03228               0.00000       0.00000                      0.02339
# 6 0.01013              90.15064      90.15191                      0.00364
# 7 0.01171              86.78685      86.84028                      0.05876
# 8 0.01540              77.16125      77.23989                      0.08729
# Variance_Percent_Difference SE_Absolute_Difference
# 1                     0.41168                  6e-05
# 2                     0.23542                  2e-05
# 3                     0.50785                  5e-05
# 4                     0.36670                  3e-05
# 5                     0.02245                  0e+00
# 6                     0.03542                  0e+00
# 7                     0.42667                  3e-05
# 8                     0.36670                  3e-05



# ============================================================
# COMPACT VARIANCE COMPARISON TABLE
# ============================================================

variance_comparison <- data.frame(
  
  Option = comparison_table$Option,
  
  Method = comparison_table$Method,
  
  Analytical_Variance =
    comparison_table$Analytical_Variance,
  
  MC_Variance =
    comparison_table$MC_Variance,
  
  Absolute_Difference =
    comparison_table$Variance_Absolute_Difference,
  
  Percent_Difference =
    comparison_table$Variance_Percent_Difference
)

print(variance_comparison)

#       Option      Method    Analytical_Variance   MC_Variance   Absolute_Difference
# 1 European Call Standard MC            75.55152    75.24049             0.31103
# 2 European Call  Antithetic            29.10468    29.03616             0.06852
# 3 European Call      W_T CV            33.43958    33.26976             0.16982
# 4 European Call      S_T CV            23.80329    23.71600             0.08729
# 5  European Put Standard MC           104.22323   104.19984             0.02339
# 6  European Put  Antithetic            10.26533    10.26169             0.00364
# 7  European Put      W_T CV            13.77117    13.71241             0.05876
# 8  European Put      S_T CV            23.80329    23.71600             0.08729
# Percent_Difference
# 1            0.41168
# 2            0.23542
# 3            0.50785
# 4            0.36670
# 5            0.02245
# 6            0.03542
# 7            0.42667
# 8            0.36670

# ============================================================
# COMPACT VARIANCE REDUCTION COMPARISON
# ============================================================

VR_comparison <- data.frame(
  
  Option = c(
    "European Call",
    "European Put"
  ),
  
  Analytical_Antithetic_VR = c(
    Call_AV_VR,
    Put_AV_VR
  ),
  
  MC_Antithetic_VR = c(
    MC_Call_VR["Antithetic"],
    MC_Put_VR["Antithetic"]
  ),
  
  Analytical_WT_CV_VR = c(
    Call_W_VR,
    Put_W_VR
  ),
  
  MC_WT_CV_VR = c(
    MC_Call_VR["WT_CV"],
    MC_Put_VR["WT_CV"]
  ),
  
  Analytical_ST_CV_VR = c(
    Call_ST_VR,
    Put_ST_VR
  ),
  
  MC_ST_CV_VR = c(
    MC_Call_VR["ST_CV"],
    MC_Put_VR["ST_CV"]
  )
)


VR_comparison[-1] <-
  round(
    VR_comparison[-1],
    5
  )

print(VR_comparison)


# ============================================================
#  BLACK-SCHOLES AND CORRELATION INFORMATION
# ============================================================

additional_results <- data.frame(
  
  Option = c(
    "European Call",
    "European Put"
  ),
  
  Black_Scholes_Price = c(
    Call_BS,
    Put_BS
  ),
  
  WT_Correlation = c(
    Call_rho_W,
    Put_rho_W
  ),
  
  ST_Correlation = c(
    Call_rho_ST,
    Put_rho_ST
  )
)


additional_results[-1] <-
  round(
    additional_results[-1],
    5
  )

print(additional_results)


#####################################################
# ============================================================
# ANALYTICAL CORRELATIONS BETWEEN EUROPEAN OPTION PAYOFFS
# AND CONTROL VARIATES W_T AND S_T
# ============================================================

# ------------------------------------------------------------
# PARAMETERS
# ------------------------------------------------------------

S0    <- 80
K     <- 90
r     <- 0.05
q     <- 0
sigma <- 0.20
T     <- 1

Phi <- pnorm
phi <- dnorm


# ------------------------------------------------------------
#  BLACK-SCHOLES QUANTITIES
# ------------------------------------------------------------

d1 <- (
  log(S0 / K) +
    (r - q + 0.5 * sigma^2) * T
) / (sigma * sqrt(T))

d2 <- d1 - sigma * sqrt(T)


# Black-Scholes prices

Call_BS <-
  S0 * exp(-q * T) * Phi(d1) -
  K * exp(-r * T) * Phi(d2)

Put_BS <-
  K * exp(-r * T) * Phi(-d2) -
  S0 * exp(-q * T) * Phi(-d1)


# ============================================================
# 3 ANALYTICAL VARIANCE OF OPTION PAYOFFS
# ============================================================

# ------------------------------------------------------------
# Call truncated moments
# ------------------------------------------------------------

A0 <- Phi(d2)

A1 <-
  S0 *
  exp((r - q) * T) *
  Phi(d1)

A2 <-
  S0^2 *
  exp(2 * (r - q) * T + sigma^2 * T) *
  Phi(d1 + sigma * sqrt(T))


Call_E2 <-
  exp(-2 * r * T) *
  (
    A2 -
      2 * K * A1 +
      K^2 * A0
  )

Call_Var <-
  Call_E2 - Call_BS^2


# ------------------------------------------------------------
# Put truncated moments
# ------------------------------------------------------------

B0 <- Phi(-d2)

B1 <-
  S0 *
  exp((r - q) * T) *
  Phi(-d1)

B2 <-
  S0^2 *
  exp(2 * (r - q) * T + sigma^2 * T) *
  Phi(-d1 - sigma * sqrt(T))


Put_E2 <-
  exp(-2 * r * T) *
  (
    K^2 * B0 -
      2 * K * B1 +
      B2
  )

Put_Var <-
  Put_E2 - Put_BS^2


# ============================================================
#  CORRELATION WITH W_T
# ============================================================

# W_T ~ N(0,T)
# Therefore Var(W_T) = T

Var_WT <- T


# ------------------------------------------------------------
# Call covariance with W_T
# ------------------------------------------------------------

Call_Cov_WT <-
  exp(-r * T) *
  sqrt(T) *
  (
    S0 *
      exp((r - q) * T) *
      (
        sigma * sqrt(T) * Phi(d1) +
          phi(d1)
      )
    -
      K * phi(d2)
  )


# Call correlation

Call_Corr_WT <-
  Call_Cov_WT /
  sqrt(Call_Var * Var_WT)


# ------------------------------------------------------------
# Put covariance with W_T
# ------------------------------------------------------------

Put_Cov_WT <-
  exp(-r * T) *
  sqrt(T) *
  (
    -K * phi(d2)
    -
      S0 *
      exp((r - q) * T) *
      (
        sigma * sqrt(T) * Phi(-d1) -
          phi(d1)
      )
  )


# Put correlation

Put_Corr_WT <-
  Put_Cov_WT /
  sqrt(Put_Var * Var_WT)


# ============================================================
# CORRELATION WITH DISCOUNTED S_T
# ============================================================

# Control variable:
# Y = exp(-rT) * S_T

# Expected value of control

E_ST_Control <-
  S0 * exp(-q * T)


# Variance of discounted S_T

Var_ST_Control <-
  S0^2 *
  exp(-2 * q * T) *
  (
    exp(sigma^2 * T) - 1
  )


# ------------------------------------------------------------
# Call covariance with discounted S_T
# ------------------------------------------------------------

Call_Cov_ST <-
  exp(-2 * r * T) *
  (
    A2 -
      K * A1
  ) -
  Call_BS * E_ST_Control


# Call correlation

Call_Corr_ST <-
  Call_Cov_ST /
  sqrt(
    Call_Var *
      Var_ST_Control
  )


# ------------------------------------------------------------
# Put covariance with discounted S_T
# ------------------------------------------------------------

Put_Cov_ST <-
  exp(-2 * r * T) *
  (
    K * B1 -
      B2
  ) -
  Put_BS * E_ST_Control


# Put correlation

Put_Corr_ST <-
  Put_Cov_ST /
  sqrt(
    Put_Var *
      Var_ST_Control
  )


# ============================================================
# RESULTS
# ============================================================

correlation_results <- data.frame(
  
  Option = c(
    "European Call",
    "European Put"
  ),
  
  WT_Covariance = c(
    Call_Cov_WT,
    Put_Cov_WT
  ),
  
  WT_Correlation = c(
    Call_Corr_WT,
    Put_Corr_WT
  ),
  
  ST_Covariance = c(
    Call_Cov_ST,
    Put_Cov_ST
  ),
  
  ST_Correlation = c(
    Call_Corr_ST,
    Put_Corr_ST
  )
)


# Round results

correlation_results[-1] <-
  round(
    correlation_results[-1],
    5
  )


print(correlation_results)

#    Option        WT_Covariance    WT_Correlation ST_Covariance ST_Correlation
# 1  European Call    6.48937        0.74659       116.2586        0.82761
# 2  European Put      -9.51063     -0.93159       -144.9303       -0.87841

# ============================================================
# CHECK CORRELATION AGAINST VARIANCE REDUCTION
# ============================================================

# For an optimal control variate:
#
# Var(X*) = Var(X)(1-rho^2)
#
# Therefore: VR (%) = 100*rho^2


correlation_VR <- data.frame(
  
  Option = c( "European Call", "European Put" ),
  
  WT_Correlation = c(Call_Corr_WT, Put_Corr_WT ),
  WT_Implied_VR_Percent = c( 100 * Call_Corr_WT^2, 100 * Put_Corr_WT^2 ),
  
  ST_Correlation = c( Call_Corr_ST, Put_Corr_ST ),
  ST_Implied_VR_Percent = c( 100 * Call_Corr_ST^2, 100 * Put_Corr_ST^2 ) )

correlation_VR[-1] <- round( correlation_VR[-1],5 )

print(correlation_VR)
























































