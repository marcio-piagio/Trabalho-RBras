# -*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
# -*-*-*- Rotinas: Estimação por Máxima Verossimilhança para o Modelo Mtb -*-
# -*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*

# --- Função para gerar dados de captura e recaptura 
gera_dados <- function(N,K,pj,c){
  # Inicialização de vetores para armazenar os dados de captura e marcação
  u = m = n = M = c()  
  etaj = log(pj/(1-pj))               # Transformação logística para cálculo da probabilidade
  cj   = exp(etaj+c)/(1+exp(etaj+c))  # Probabilidade de recaptura ajustada
  u[1] = rbinom(1,N,pj[1])            # Número de animais capturados inicialmente
  m[1] = 0                            # Nenhum animal marcado na primeira amostragem
  n[1] = u[1]+m[1]                    # Total de animais capturados na primeira amostragem
  M[1] = 0                            # Nenhum animal marcado na população
  
  # Segunda amostragem em diante
  for(i in 2:K){
    M[i] = M[i-1] + u[i-1]            # Atualização do número de animais marcados na população
    u[i] = rbinom(1,N-M[i],pj[i])     # Número de novos animais capturados
    m[i] = rbinom(1,M[i],cj[i])       # Número de animais recapturados
    n[i] = u[i] + m[i]                # Total de animais capturados na amostragem i
  }
  r = M[K] + u[K]                     # Número total de animais capturados ao longo do estudo
  return(c(u,m,M,n,r))                # Retorno dos resultados em um vetor
}

library(Rcpp)
sourceCpp("gera_dados.cpp")  # Importa versão otimizada da função escrita em C++

# --- Função de Log-Verossimilhança do Modelo Mtb
Likelihood_function_Mtb <- function(parameters,K,m,M,n,r){
  N   = parameters[1]                  # Estimativa da população total
  eta = parameters[2:(K+1)]             # Parâmetros da probabilidade de captura
  c   = parameters[K+2]                 # Parâmetro da probabilidade de recaptura
  
  # Cálculo da função de log-verossimilhança baseada no modelo
  l1  = lgamma(N+1)
  l2  = - lgamma(N - r+1)
  l3  = sum(eta*n)
  l4  = c*sum(m)
  l5  = - sum((N-M)*log(1 + exp(eta)))
  l6  = - sum(M*log(1 + exp(c + eta)))
  
  log.like = l1 + l2 + l3 + l4 + l5 + l6   # Cálculo final
  return(-log.like)                        # Retorno do log-verossimilhança negativo (para otimização)
}

# --- Função de Log-Verossimilhança do Modelo Mt (sem ajuste de recaptura)
Likelihood_function_Mt <- function(parameters,K,m,M,n,r){
  N   = parameters[1]                  # Estimativa da população total
  eta = parameters[2:(K+1)]             # Parâmetros da probabilidade de captura
  
  # Cálculo da função de log-verossimilhança baseada no modelo simplificado
  l1 = lgamma(N+1)
  l2 = - lgamma(N - r+1)
  l3 = sum(eta*n)
  l4 = - sum((N-M)*log(1 + exp(eta))) 
  l5 = - sum(M*log(1 + exp(eta)))
  
  log.like = l1 + l2 + l3 + l4 + l5  # Cálculo final
  return(-log.like)                   # Retorno do log-verossimilhança negativo (para otimização)
}

# --- Estimação por Máxima Verossimilhança para Mtb
emv.Mtb <- function(PAR, k) {
  u = PAR[1:k]
  m = PAR[seq(1 + k, 2 * k, 1)]
  M = PAR[seq(1 + 2 * k, 3 * k, 1)]
  n = PAR[seq(1 + 3 * k, 4 * k, 1)]
  r = PAR[1 + 4 * k]
  
  # Inicialização das estimativas
  N.estMtb <- NA
  c.estMtb <- NA
  L.maxMtb <- NA
  
  if (sum(m) != 0) {  # Apenas executa se houve recapturas
    par = c(r + 100, rep(0, k + 1))   # Parâmetros iniciais
    maximization <- try(optim(par, Likelihood_function_Mtb,
                              method = "L-BFGS-B",
                              lower = c(r, rep(-Inf, k + 1)),
                              K = k, m = m, M = M, n = n, r = r),
                        silent = TRUE)  # Tentativa de otimização silenciosa
    
    if (!inherits(maximization, "try-error")) {  
      N.estMtb <- maximization$par[1]   # Estimativa do tamanho populacional
      c.estMtb <- maximization$par[k + 2]  # Estimativa do parâmetro de recaptura
      L.maxMtb <- maximization$value   # Valor máximo da função de verossimilhança
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

# --- Configuração dos parâmetros do estudo
N.true <- seq(50,500,50)          # Valores reais do tamanho populacional
K.true <- seq(5,15,1)             # Número de amostragens variando
c.true <- seq(-2,2,0.1)           # Valores reais do parâmetro de recaptura
no.simulation <- 10           # Número de simulações

# --- Scenario: teste 
# N.true <- c(500,1000)
# K.true <- seq(5,15,2)
# c.true <- seq(-2,2,0.2)
# no.simulation <- 100

# Matrizes para armazenar os critérios estatísticos das simulações
M <- matrix(nrow = length(c.true),ncol = length(N.true))
AIC        <- BIC        <- AICc <- HQIC <- CAIC <- M
TRV.90     <- TRV.95     <- TRV.99               <- M
Vies.Mt.N  <- Vies.Mtb.N <- Vies.Mtb.c           <- M
EQM.Mt.N   <- EQM.Mtb.N  <- EQM.Mtb.c            <- M


resultados   <- matrix.N <- matrix.c <-  list()  # Listas para armazenar os resultados da simulação

library(parallel)
# --- Configuração do processamento paralelo
cl <- makeCluster(detectCores() - 3)  # Cria um cluster, deixando 3 núcleos livres para outras tarefas

# Inicializa os workers do cluster com os pacotes e funções necessárias
clusterEvalQ(cl, {
  library(Rcpp)
  sourceCpp("gera_dados.cpp")  # Carrega a versão otimizada da geração de dados em C++
})

# Exportação de variáveis para o cluster paralelo
clusterExport(cl, c("N.true", "K.true", "c.true", "gera_dados_cpp",
                    "emv.Mt", "emv.Mtb", 
                    "Likelihood_function_Mt", "Likelihood_function_Mtb"))


# --- Início da simulação
tictoc::tic()  # Inicia a medição do tempo total de execução

for (k in 1:length(K.true)) {
  set.seed(k)  # Define a semente para replicabilidade
  pj <- runif(K.true[k] * no.simulation, 0.10, 0.30)  # Gera probabilidades de captura aleatórias
  X  <- matrix(pj, K.true[k], no.simulation)  # Organiza os valores em uma matriz
  
  # --- Inicializa matrizes para armazenar as estimativas
  fit.N.T <- matrix(nrow = 0, ncol = length(N.true) )  
  fit.c.T <- matrix(nrow = 0, ncol = length(N.true) )
  
  for (j in 1:length(c.true)){
    fit.N <- matrix(nrow = no.simulation * 2, ncol = 0)
    fit.c <- matrix(nrow = no.simulation, ncol = 0)
    
    for (i in 1:length(N.true)) {
      
      # --- Geração dos dados de captura e recaptura
      dados <- apply(X, 2, function(N, K, pj, c) gera_dados_cpp(N, K, pj, c),
                     N=N.true[i], K=K.true[k], c=c.true[j])
      
      # --- Estimação por Máxima Verossimilhança usando os modelos
      fit.Mt  <- parApply(cl, dados, 2, emv.Mt, k=K.true[k])        # Modelo Mt
      fit.Mtb <- parApply(cl, dados, 2, emv.Mtb, k=K.true[k])       # Modelo Mtb
      
      # --- Cálculo dos critérios estatísticos para avaliação dos modelos
      AIC.Mt  <- 2 * fit.Mt[2,] + 2 * (K.true[k]+1)
      BIC.Mt  <- 2 * fit.Mt[2,] + (K.true[k]+1) * log(fit.Mt[3,])
      AICc.Mt <- AIC.Mt + 2 * ((K.true[k]+1)^2 + K.true[k]+1) / (fit.Mt[3,] - (K.true[k]+1) + 1)
      HQIC.Mt <- 2 * fit.Mt[2,] + 2 * (K.true[k]+1) * log(log(fit.Mt[3,])) 
      CAIC.Mt <- 2 * fit.Mt[2,] + (K.true[k]+1) * (log(fit.Mt[3,]) + 1)
      
      AIC.Mtb  <- 2 * fit.Mtb[3,] + 2 * (K.true[k]+2)
      BIC.Mtb  <- 2 * fit.Mtb[3,] + (K.true[k]+2) * log(fit.Mtb[4,])
      AICc.Mtb <- AIC.Mtb + 2 * ((K.true[k]+2)^2 + K.true[k]+2) / (fit.Mtb[4,] - (K.true[k]+2) + 1)
      HQIC.Mtb <- 2 * fit.Mtb[3,] + 2 * (K.true[k]+2) * log(log(fit.Mtb[4,]))
      CAIC.Mtb <- 2 * fit.Mtb[3,] + (K.true[k]+1) * (log(fit.Mtb[4,]) + 1)
      
      AIC[j,i]  <- mean(1*(AIC.Mtb < AIC.Mt),na.rm=TRUE)
      BIC[j,i]  <- mean(1*(BIC.Mtb < BIC.Mt),na.rm=TRUE)
      AICc[j,i] <- mean(1*(AICc.Mtb < AICc.Mt),na.rm=TRUE)
      HQIC[j,i] <- mean(1*(HQIC.Mtb < HQIC.Mt),na.rm=TRUE)
      CAIC[j,i] <- mean(1*(CAIC.Mtb < CAIC.Mt),na.rm=TRUE)
      
      # --- Armazena as estimativas obtidas
      fit.N <- cbind(fit.N, c(fit.Mt[1,], fit.Mtb[1,]))
      fit.c <- cbind(fit.c, fit.Mtb[2,])
      
      # --- Teste de razão de verossimilhança (TRV)
      TRV.90[j,i] <- mean(1 * (2 * (fit.Mt[2,] - fit.Mtb[3,]) > qchisq(0.90,1)), na.rm=TRUE)
      TRV.95[j,i] <- mean(1 * (2 * (fit.Mt[2,] - fit.Mtb[3,]) > qchisq(0.95,1)), na.rm=TRUE)
      TRV.99[j,i] <- mean(1 * (2 * (fit.Mt[2,] - fit.Mtb[3,]) > qchisq(0.99,1)), na.rm=TRUE)
      
      # --- Cálculo do viés
      Vies.Mt.N[j,i]  <- mean((fit.Mt[1,]  - N.true[i]) / N.true[i], na.rm=TRUE)
      Vies.Mtb.N[j,i] <- mean((fit.Mtb[1,] - N.true[i]) / N.true[i], na.rm=TRUE)
      Vies.Mtb.c[j,i] <- mean(fit.Mtb[2,] - c.true[j], na.rm=TRUE)
      
      # --- Cálculo do Erro Quadrático Médio (EQM)
      EQM.Mt.N[j,i]  <- sqrt(mean((fit.Mt[1,]  - N.true[i])^2, na.rm=TRUE))
      EQM.Mtb.N[j,i] <- sqrt(mean((fit.Mtb[1,] - N.true[i])^2, na.rm=TRUE))
      EQM.Mtb.c[j,i] <- sqrt(mean((fit.Mtb[2,] - c.true[j])^2, na.rm=TRUE))
      
      cat(K.true[k], N.true[i], c.true[j], "\n")  # Impressão para acompanhar o progresso
    }
    
    # --- Adicione os resultados para o j atual de fit.N.T e fit.c.T
    fit.N.T <- rbind(fit.N.T, fit.N)
    fit.c.T <- rbind(fit.c.T, fit.c)
  }
  # --- Armazene os resultados para o k atual
  matrix.N[[k]] <- fit.N.T
  matrix.c[[k]] <- fit.c.T
  
  resultados[[k]] <-  list("AIC"=AIC, "BIC"=BIC, "AICc"=AICc, "HQIC"=HQIC, "CAIC"=CAIC,
                           "TRV.90"=TRV.90, "TRV.95"=TRV.95, "TRV.99"=TRV.99,
                           "Vies.Mt.N"=Vies.Mt.N, "Vies.Mtb.N"=Vies.Mtb.N, "Vies.Mtb.c"=Vies.Mtb.c,
                           "EQM.Mt.N"=EQM.Mt.N, "EQM.Mtb.N"=EQM.Mtb.N, "EQM.Mtb.c"=EQM.Mtb.c)
}
tictoc::toc()  # Finaliza a medição do tempo total de execução

stopCluster(cl)

matrix.N[[k]]
