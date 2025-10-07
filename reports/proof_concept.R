set.seed(3423457) # set seed for reproducibility
nrep <- 10000   # number of replications
BASE <- matrix(nrow = nrep, ncol = 3)
BASEalt <- matrix(nrow = nrep, ncol = 3)
OLS <- matrix(nrow = nrep, ncol = 3)
OLSalt <- matrix(nrow = nrep, ncol = 3)
MINT <- matrix(nrow = nrep, ncol = 3)
MINTalt <- matrix(nrow = nrep, ncol = 3)
ICOMB <- matrix(nrow = nrep, ncol = 3)
ICOMBalt <- matrix(nrow = nrep, ncol = 3)
sigb1 <- 1   # the error st dev for the first bottom level
sigb2 <- 1   # the error st dev for the second bottom level
sigf <- 1    # the error st dev for the common factor
rho <- 0.6   # the AR parameter of the common factor
n <- 1000000 # lazy way of getting the theoretical values from simulation
e1 <- rnorm(n + 1)
e2 <- rnorm(n + 1)
ef <- rnorm(n + 1)
f <- rep(0, n + 1)
b1f <- rep(0, n + 1) # base forecast for the first bottom level
b2f <- rep(0, n + 1) # base forecast for the second bottom level
yf <- rep(0, n + 1) # base forecast for the aggregate

f[1] <- 1.25 * ef[1] # initial draw from the unconditional distribution
for (i in 2:(n + 1)) {
  f[i] <- 0.6 * f[i - 1] + ef[i]
}
b1 <- 1 + 0.8 * f + e1
b2 <- 1 + 0.8 * f + e2
y <- b1 + b2
Y <- cbind(y, b1, b2)
Yin <- Y[2:(n + 1), ] # estimation sample

# Conditioning on the first observation being known
b1f[1] <- b1[1] 
b2f[1] <- b2[1] 
yf[1] <- y[1] 

# Generating base forecasts from the implied univariate ARMA(1,1) models
for (i in 2:(n+1)) {
  b1f[i] <- 1 + 0.6 * (b1[i - 1] - 1) - (b1[i - 1] - b1f[i - 1])/3
  b2f[i] <- 1 + 0.6 * (b2[i - 1] - 1) - (b2[i - 1] - b2f[i - 1])/3
  yf[i] <- 2 + 0.6 * (y[i - 1] - 2) - 0.2404082 * (y[i - 1] - yf[i - 1])
}

b2falt <- 1 + 0.8 * f # forecast of the second bottom level assuming knowledge of f_t
# base forecasts for the estimation and evaluation samples
Yfb <- cbind(yf, b1f, b2f) 
Yfbin <- Yfb[2:(n + 1), ]
Yfbalt <- cbind(yf, b1f, b2falt)
Yfbaltin <- Yfbalt[2:(n + 1), ]

E <- Yin - Yfbin # errors in the base forecast in the estimation sample
Ealt <- Yin - Yfbaltin

W <- var(E)
Walt <- var(Ealt)
invW <- solve(W)
invWalt <- solve(Walt)

S <- matrix(c(1, 1, 0, 1, 0, 1), nrow = 3, ncol = 2)

Bmint <- invW %*% S %*% solve(t(S) %*% invW %*% S) %*% t(S)

Bmintalt <- invWalt %*% S %*% solve(t(S) %*% invWalt %*% S) %*% t(S)

Bicomb <- solve(t(Yfbin) %*% Yfbin) %*% t(Yfbin) %*% Yin

Bicombalt <- solve(t(Yfbaltin) %*% Yfbaltin) %*% t(Yfbaltin) %*% Yin

Bols <- S %*% solve(t(S) %*% S) %*% t(S)

n <- 100     # just to remove the effect of initial conditions
h <- 1       # forecast horizon
sigb1 <- 1   # the error st dev for the first bottom level
sigb2 <- 1   # the error st dev for the second bottom level
sigf <- 1    # the error st dev for the common factor
rho <- 0.6   # the AR parameter of the common factor
for (irep in 1:nrep) {
  # All analysis is conditional on the first observation being known
  # So, generate one more data point than necessary
  e1 <- rnorm(n + h + 1)
  e2 <- rnorm(n + h + 1)
  ef <- rnorm(n + h + 1)
  f <- rep(0, n + h + 1)
  b1f <- rep(0, n + h + 1) # base forecast for the first bottom level
  b2f <- rep(0, n + h + 1) # base forecast for the second bottom level
  yf <- rep(0, n + h + 1) # base forecast for the aggregate
  
  f[1] <- 1.25 * ef[1] # initial draw from the unconditional distribution
  for (i in 2:(n + h + 1)) {
    f[i] <- 0.6 * f[i - 1] + ef[i]
  }
  b1 <- 1 + 0.8 * f + e1
  b2 <- 1 + 0.8 * f + e2
  y <- b1 + b2
  Y <- cbind(y, b1, b2)
  
  Yout <- Y[n + h + 1, ]  # evaluation observation
  
  # Conditioning on the first observation being known
  b1f[1] <- b1[1] 
  b2f[1] <- b2[1] 
  yf[1] <- y[1] 
  
  # Generating base forecasts from the implied univariate ARMA(1,1) models
  for (i in 2:(n + h + 1)) {
    b1f[i] <- 1 + 0.6 * (b1[i - 1] - 1) - (b1[i - 1] - b1f[i - 1])/3
    b2f[i] <- 1 + 0.6 * (b2[i - 1] - 1) - (b2[i - 1] - b2f[i - 1])/3
    yf[i] <- 2 + 0.6 * (y[i - 1] - 2) - 0.2404082 * (y[i - 1] - yf[i - 1])
  }
  b2falt <- 1 + 0.8 * f # forecast of the second bottom level assuming knowledge of f_t
  # base forecasts for evaluation
  Yfb <- cbind(yf, b1f, b2f) 
  Yfbout <- Yfb[n + h + 1, ] # evaluation forecast
  Yfbalt <- cbind(yf, b1f, b2falt)
  Yfbaltout <- Yfbalt[n + h + 1, ]
  
  mintf <- Yfbout %*% Bmint
  mintfalt <- Yfbaltout %*% Bmintalt
  
  icombf <- Yfbout %*% Bicomb
  icombfalt <- Yfbaltout %*% Bicombalt
  
  olsf <- Yfbout %*% Bols
  olsfalt <- Yfbaltout %*% Bols
  
  BASE[irep, ] <- Yout - Yfbout
  BASEalt[irep, ] <- Yout - Yfbaltout
  
  OLS[irep, ] <- Yout - olsf
  OLSalt[irep, ] <- Yout - olsfalt
  
  MINT[irep, ] <- Yout - mintf
  MINTalt[irep, ] <- Yout - mintfalt
  
  ICOMB[irep, ] <- Yout - icombf
  ICOMBalt[irep, ] <- Yout - icombfalt
}

baseMSE <- sum(BASE^2)/nrep
olsMSE <- sum(OLS^2)/nrep
mintMSE <- sum(MINT^2)/nrep
icombMSE <- sum(ICOMB^2)/nrep
MSE <- data.frame(baseMSE, olsMSE, mintMSE, icombMSE)

basealtMSE <- sum(BASEalt^2)/nrep
olsaltMSE <- sum(OLSalt^2)/nrep
mintaltMSE <- sum(MINTalt^2)/nrep
icombaltMSE <- sum(ICOMBalt^2)/nrep
altMSE <- data.frame(basealtMSE, olsaltMSE, mintaltMSE, icombaltMSE)

print("MSE based on 3 univariate ARMA forecasts")

print.data.frame(MSE)

print("MSE based on 2 ARMA + oracle forecast for 2nd bottom level")

print.data.frame(altMSE)


