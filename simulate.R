#-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-#
#-*-*-*-*-*-*- routines: paper Maximum likelihood Mtb model (Wang) -*-*-*-*-*-#
#-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-#

# --- Draw sample:
# pj capture probabilities (initial)
gera_dados <- function(N,K,pj,c){
  u = m = n = M = c()
  etaj = log(pj/(1-pj))
  cj   = exp(etaj+c)/(1+exp(etaj+c))  # probability of recapture
  u[1] = rbinom(1,N,pj[1])            # number of animals captured
  m[1] = 0                            # number of animals tagged
  n[1] = u[1]+m[1]                    # number of animals selected
  M[1] = 0                            # number of animals tagged in the population
  # second sample ahead
  for(i in 2:K){
    M[i] = M[i-1] + u[i-1]
    u[i] = rbinom(1,N-M[i],pj[i])
    m[i] = rbinom(1,M[i],cj[i])
    n[i] = u[i] + m[i]
  }
  r = M[K] + u[K]
  return(c(u,m,M,n,r))
}

library(Rcpp)
sourceCpp("gera_dados.cpp")

# --- Likelihood function Model: Mtb
Likelihood_function_Mtb <- function(parameters,K,m,M,n,r){
  N   = parameters[1]
  eta = parameters[2:(K+1)]
  c   = parameters[K+2]
  l1  = lgamma(N+1)
  l2  = - lgamma(N - r+1)
  l3  = sum(eta*n)
  l4  = c*sum(m)
  l5  = - sum((N-M)*log(1 + exp(eta)))
  l6  = - sum(M*log(1 + exp(c + eta)))
  log.like = l1 + l2 + l3 + l4 + l5 +l6
  return(-log.like)
}

# --- Likelihood function Model: Mt
Likelihood_function_Mt <- function(parameters,K,m,M,n,r){
  N = parameters[1]
  eta = parameters[2:(K+1)]
  l1 = lgamma(N+1)
  l2 = - lgamma(N - r+1)
  l3 = sum(eta*n)
  l4 = - sum((N-M)*log(1 + exp(eta))) 
  l5 = -  sum(M*log(1 + exp(eta)))
  log.like = l1 + l2 + l3 + l4 + l5
  return(-log.like)
}

# --- função de estimação de maxima verossimilhanca Mtb
emv.Mtb <- function(PAR, k) {
  u = PAR[1:k]
  m = PAR[seq(1 + k, 2 * k, 1)]
  M = PAR[seq(1 + 2 * k, 3 * k, 1)]
  n = PAR[seq(1 + 3 * k, 4 * k, 1)]
  r = PAR[1 + 4 * k]
  N.estMtb <- NA
  c.estMtb <- NA
  L.maxMtb <- NA
  
  if (sum(m) != 0) {
    par = c(r + 100, rep(0, k + 1))
    maximization <- try(optim(par, Likelihood_function_Mtb,
                              method = "L-BFGS-B",
                              lower = c(r, rep(-Inf, k + 1)),
                              K = k, m = m, M = M, n = n, r = r),
                        silent = TRUE)
    
    if (!inherits(maximization, "try-error")) {
      N.estMtb <- maximization$par[1]
      c.estMtb <- maximization$par[k + 2]
      L.maxMtb <- maximization$value
    }
  }
  
  return(rbind('N' = N.estMtb, 'c' = c.estMtb, "loglik" = L.maxMtb, "R" = r))
}

# --- função de estimação de maxima verossimilhanca Mt
emv.Mt <- function(PAR, k) {
  u = PAR[1:k]
  m = PAR[seq(1 + k, 2 * k, 1)]
  M = PAR[seq(1 + 2 * k, 3 * k, 1)]
  n = PAR[seq(1 + 3 * k, 4 * k, 1)]
  r = PAR[1 + 4 * k]
  N.estMt <- NA
  L.maxMt <- NA  
  
  if (sum(m) != 0) {
    par = c(r + 100, rep(0, k))
    maximization <- try(optim(par, Likelihood_function_Mt,
                             method = "L-BFGS-B",
                             lower = c(r, rep(-Inf, k)),
                             K = k, m = m, M = M, n = n, r = r),
                       silent = TRUE)
    
    if (!inherits(maximization, "try-error")) { 
      N.estMt <- maximization$par[1]
      L.maxMt <- maximization$value
    }
  }
  
  return(rbind('N' = N.estMt, "loglik" = L.maxMt, "R" = r))
}

# --- Scenario:
N.true <- c(500,1000)
K.true <- seq(5,15,2)
c.true <- seq(-2,2,0.2)
no.simulation <- 10

AIC          <- BIC        <- AICc       <- HQIC <- matrix(NA,length(c.true),length(N.true))
TRV.90       <- TRV.95     <- TRV.99     <- matrix(NA,length(c.true),length(N.true))
Vies.Mt.N    <- Vies.Mtb.N <- Vies.Mtb.c <- matrix(NA,length(c.true),length(N.true))
EQM.Mt.N     <- EQM.Mtb.N  <- EQM.Mtb.c  <- matrix(NA,length(c.true),length(N.true))

resultados   <- matrix.N <- matrix.c <-  list()

library(parallel)
cl <- makeCluster(detectCores() - 3)

clusterEvalQ(cl, {
  library(Rcpp)
  sourceCpp("gera_dados.cpp")
})

# Exporte todas as dependências
clusterExport(cl, c("N.true", "K.true", "c.true", "gera_dados_cpp",
                    "emv.Mt", "emv.Mtb", 
                    "Likelihood_function_Mt", "Likelihood_function_Mtb"))


tictoc::tic()
for (k in 1:length(K.true)) {
  set.seed(k)
  pj <- runif(K.true[k] * no.simulation, 0.10, 0.30)
  X  <- matrix(pj, K.true[k], no.simulation)
  
  # Initialize fit.N.T and fit.c.T for each k with appropriate columns
  fit.N.T <- matrix(nrow = 0, ncol = 2)  
  fit.c.T <- matrix(nrow = 0, ncol = 2)
  
  for (j in 1:length(c.true)){
    # Initialize fit.N and fit.c as empty matrices for each j
    fit.N <- matrix(nrow = no.simulation * 2, ncol = 0)
    fit.c <- matrix(nrow = no.simulation, ncol = 0)
    
    for (i in 1:length(N.true)) {
      
      # Gerandor de dados 
      dados <- apply(X, 2, function(N, K, pj, c)gera_dados_cpp(N, K, pj, c),
                     N=N.true[i], K=K.true[k], c=c.true[j])
      
      # Estimativas Mt
      fit.Mt  <- parApply(cl,dados, 2, emv.Mt, k=K.true[k])
      
      AIC.Mt  <- 2 * fit.Mt[2,] + 2 * (K.true[k]+1)
      BIC.Mt  <- 2 * fit.Mt[2,] + (K.true[k]+1)*log(fit.Mt[3,])
      AICc.Mt <- AIC.Mt + 2 * ((K.true[k]+1)^2+K.true[k]+1)/(fit.Mt[3,]-(K.true[k]+1)+1)
      HQIC.Mt <- 2 * fit.Mt[2,] + 2  *(K.true[k]+1)  * log(log(fit.Mt[3,])) 
      
      # Estimativas Mtb
      fit.Mtb  <- apply(dados,2, emv.Mtb, k=K.true[k])
      
      AIC.Mtb  <- 2 * fit.Mtb[3,] + 2 * (K.true[k]+2)
      BIC.Mtb  <- 2 * fit.Mtb[3,] + (K.true[k]+2)*log(fit.Mtb[4,])
      AICc.Mtb <- AIC.Mtb + 2*((K.true[k]+2)^2+K.true[k]+2)/(fit.Mtb[4,]-(K.true[k]+2)+1)
      HQIC.Mtb <- + 2 *fit.Mtb[3,] + 2  * (K.true[k]+2)  * log(log(fit.Mtb[4,]))
      
      # Collect N estimates
      fit.N <- cbind(fit.N,c(fit.Mt[1,],fit.Mtb[1,]))
      
      # Collect c estimates
      fit.c <- cbind(fit.c,fit.Mtb[2,])
      
      # Critério
      AIC[j,i]  <- mean(1*(AIC.Mtb < AIC.Mt),na.rm=TRUE)
      BIC[j,i]  <- mean(1*(BIC.Mtb < BIC.Mt),na.rm=TRUE)
      AICc[j,i] <- mean(1*(AICc.Mtb < AICc.Mt),na.rm=TRUE)
      HQIC[j,i] <- mean(1*(HQIC.Mtb < HQIC.Mt),na.rm=TRUE)
      
      # Teste 
      TRV.90[j,i]    <- mean(1*(2*(fit.Mt[2,]-fit.Mtb[3,])>qchisq(0.90,1)),na.rm=TRUE)
      TRV.95[j,i]    <- mean(1*(2*(fit.Mt[2,]-fit.Mtb[3,])>qchisq(0.95,1)),na.rm=TRUE)
      TRV.99[j,i]    <- mean(1*(2*(fit.Mt[2,]-fit.Mtb[3,])>qchisq(0.99,1)),na.rm=TRUE)
      
      # Viés
      Vies.Mt.N[j,i]  <- mean((fit.Mt[1,]  - N.true[i])/N.true[i],na.rm = TRUE)
      Vies.Mtb.N[j,i] <- mean((fit.Mtb[1,] - N.true[i])/N.true[i],na.rm = TRUE)
      Vies.Mtb.c[j,i] <- mean(fit.Mtb[2,] - c.true[j],na.rm = TRUE)
      
      # EQM
      EQM.Mt.N[j,i]  <- sqrt(mean((fit.Mt[1,]  - N.true[i])^2,na.rm = TRUE))
      EQM.Mtb.N[j,i] <- sqrt(mean((fit.Mtb[1,] - N.true[i])^2,na.rm = TRUE))
      EQM.Mtb.c[j,i] <- sqrt(mean((fit.Mtb[2,] - c.true[j])^2,na.rm = TRUE))
      
      cat(K.true[k],N.true[i],c.true[j],"\n")
    }
    
    fit.N.T <- rbind(fit.N.T, fit.N )
    fit.c.T <- rbind(fit.c.T, fit.c )
    
  }
  
  
  matrix.N[[k]]   <- fit.N.T
  matrix.c[[k]]   <- fit.c.T
  
  resultados[[k]] <-  list("AIC"=AIC, "BIC"=BIC, "AICc"=AICc, "HQIC"=HQIC,
                           "TRV.90"=TRV.90, "TRV.95"=TRV.95, "TRV.99"=TRV.99,
                           "Vies.Mt.N"=Vies.Mt.N, "Vies.Mtb.N"=Vies.Mtb.N, "Vies.Mtb.c"=Vies.Mtb.c,
                           "EQM.Mt.N"=EQM.Mt.N, "EQM.Mtb.N"=EQM.Mtb.N, "EQM.Mtb.c"=EQM.Mtb.c)
}
tictoc::toc()

