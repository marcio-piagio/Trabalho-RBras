# -*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
# -*-*-*- Otimização com família apply e paralelismo eficiente -*-
# -*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
install_if_missing <- function(packages) {
  for (pkg in packages) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      install.packages(pkg)
    }
    library(pkg, character.only = TRUE)
  }
}

# Lista de pacotes que você deseja verificar e instalar
packages <- c("Rcpp", "parallel", "tictoc")

install_if_missing(packages)

Rcpp::sourceCpp("gera_dados.cpp")  # Importa versão otimizada da função escrita em C++

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

# Função principal de processamento
process_k <- function(k) {
  set.seed(k)
  K <- K.true[k]
  X <- matrix(
    runif(K * no.simulation, 0.10, 0.30),
    nrow = K,
    ncol = no.simulation
  )
  
  process_j <- function(j) {
    process_i <- function(i) {
      set.seed(k * i)
      dados <- apply(X, 2, function(col) {
        gera_dados_cpp(N.true[i], K, col, c.true[j])
      })
      
      fit.Mt <- apply(dados, 2, emv.Mt, k = K)
      fit.Mtb <- apply(dados, 2, emv.Mtb, k = K)
      
      # Cálculo de métricas
      calc_metrics(fit.Mt, fit.Mtb, K, i, j)
    }
    
    i_results <- lapply(seq_along(N.true), process_i)
    combine_i_results(i_results)
  }
  
  j_results <- lapply(seq_along(c.true), process_j)
  combine_j_results(j_results, K)
}

# Funções auxiliares
calc_metrics <- function(fit.Mt, fit.Mtb, K, i, j) {
  
  # --- Cálculo dos critérios estatísticos para avaliação dos modelos
  AIC.Mt  = 2 * fit.Mt[2,] + 2 * (K+1)
  BIC.Mt  = 2 * fit.Mt[2,] + (K+1) * log(fit.Mt[3,])
  AICc.Mt = AIC.Mt + 2 * ((K+1)^2 + K+1) / (fit.Mt[3,] - (K+1) + 1)
  HQIC.Mt = 2 * fit.Mt[2,] + 2 * (K+1) * log(log(fit.Mt[3,])) 
  CAIC.Mt = 2 * fit.Mt[2,] + (K+1) * (log(fit.Mt[3,]) + 1)
  
  AIC.Mtb  = 2 * fit.Mtb[3,] + 2 * (K+2)
  BIC.Mtb  = 2 * fit.Mtb[3,] + (K+2) * log(fit.Mtb[4,])
  AICc.Mtb = AIC.Mtb + 2 * ((K+2)^2 + K+2) / (fit.Mtb[4,] - (K+2) + 1)
  HQIC.Mtb = 2 * fit.Mtb[3,] + 2 * (K+2) * log(log(fit.Mtb[4,]))
  CAIC.Mtb = 2 * fit.Mtb[3,] + (K+2) * (log(fit.Mtb[4,]) + 1)

  
  list(
    AIC  = mean(1*(AIC.Mtb < AIC.Mt),na.rm=TRUE),
    BIC  = mean(1*(BIC.Mtb < BIC.Mt),na.rm=TRUE),
    AICc = mean(1*(AICc.Mtb < AICc.Mt),na.rm=TRUE),
    HQIC = mean(1*(HQIC.Mtb < HQIC.Mt),na.rm=TRUE),
    CAIC = mean(1*(CAIC.Mtb < CAIC.Mt),na.rm=TRUE),
    
    # --- Teste de razão de verossimilhança (TRV)
    TRV.90 = mean(1 * (2 * (fit.Mt[2,] - fit.Mtb[3,]) > qchisq(0.90,1)), na.rm=TRUE),
    TRV.95 = mean(1 * (2 * (fit.Mt[2,] - fit.Mtb[3,]) > qchisq(0.95,1)), na.rm=TRUE),
    TRV.99 = mean(1 * (2 * (fit.Mt[2,] - fit.Mtb[3,]) > qchisq(0.99,1)), na.rm=TRUE),
    
    # --- Cálculo do viés
    Vies.Mt.N  = mean((fit.Mt[1,]  - N.true[i]) / N.true[i], na.rm=TRUE),
    Vies.Mtb.N = mean((fit.Mtb[1,] - N.true[i]) / N.true[i], na.rm=TRUE),
    Vies.Mtb.c = mean(fit.Mtb[2,] - c.true[j], na.rm=TRUE),
    
    # --- Cálculo do Erro Quadrático Médio (EQM)
    EQM.Mt.N  = sqrt(mean((fit.Mt[1,]  - N.true[i])^2, na.rm=TRUE)),
    EQM.Mtb.N = sqrt(mean((fit.Mtb[1,] - N.true[i])^2, na.rm=TRUE)),
    EQM.Mtb.c = sqrt(mean((fit.Mtb[2,] - c.true[j])^2, na.rm=TRUE)),
    
    # --- Armazena as estimativas obtidas
    fit.N =c(fit.Mt[1,], fit.Mtb[1,]),
    fit.c = fit.Mtb[2,])
}

combine_i_results <- function(i_results) {
  list(
    AIC = sapply(i_results, `[[`, "AIC"),
    BIC = sapply(i_results, `[[`, "BIC"),
    AICc = sapply(i_results, `[[`, "AICc"),
    HQIC = sapply(i_results, `[[`, "HQIC"),
    CAIC = sapply(i_results, `[[`, "CAIC"),
    
    TRV.90 = sapply(i_results, `[[`, "TRV.90"),
    TRV.95 = sapply(i_results, `[[`, "TRV.95"),
    TRV.99 = sapply(i_results, `[[`, "TRV.99"),
    
    
    Vies.Mt.N = sapply(i_results, `[[`, "Vies.Mt.N"),
    Vies.Mtb.N = sapply(i_results, `[[`, "Vies.Mtb.N"),
    Vies.Mtb.c = sapply(i_results, `[[`, "Vies.Mtb.c"),
    
    EQM.Mt.N = sapply(i_results, `[[`, "EQM.Mt.N"),
    EQM.Mtb.N = sapply(i_results, `[[`, "EQM.Mtb.N"),
    EQM.Mtb.c = sapply(i_results, `[[`, "EQM.Mtb.c"),
    
    
    fit.N = do.call(cbind, lapply(i_results, function(x) x$fit.N)),
    fit.c = do.call(cbind, lapply(i_results, function(x) x$fit.c))
  )
}

combine_j_results <- function(j_results, K) {
  list(
    matrix.N = do.call(rbind, lapply(j_results, function(x) x$fit.N)),
    matrix.c = do.call(rbind, lapply(j_results, function(x) x$fit.c)),
    
    AIC = do.call(rbind, lapply(j_results, function(x) x$AIC)),
    BIC = do.call(rbind, lapply(j_results, function(x) x$BIC)),
    AICc = do.call(rbind, lapply(j_results, function(x) x$AICc)),
    HQIC = do.call(rbind, lapply(j_results, function(x) x$HQIC)),
    CAIC = do.call(rbind, lapply(j_results, function(x) x$CAIC)),
    
    TRV.90 = do.call(rbind, lapply(j_results, function(x) x$TRV.90)),
    TRV.95 = do.call(rbind, lapply(j_results, function(x) x$TRV.95)),
    TRV.99 = do.call(rbind, lapply(j_results, function(x) x$TRV.99)),
    
    Vies.Mt.N = do.call(rbind, lapply(j_results, function(x) x$Vies.Mt.N)),
    Vies.Mtb.N = do.call(rbind, lapply(j_results, function(x) x$Vies.Mtb.N)),
    Vies.Mtb.c = do.call(rbind, lapply(j_results, function(x) x$Vies.Mtb.c)),
    
    EQM.Mt.N = do.call(rbind, lapply(j_results, function(x) x$EQM.Mt.N)),
    EQM.Mtb.N = do.call(rbind, lapply(j_results, function(x) x$EQM.Mtb.N)),
    EQM.Mtb.c = do.call(rbind, lapply(j_results, function(x) x$EQM.Mtb.c))
    )
}

# --- Configuração dos parâmetros do estudo
N.true <- seq(50,500,50)          # Valores reais do tamanho populacional
K.true <- seq(5,15,1)             # Número de amostragens variando
c.true <- seq(-2,2,0.1)           # Valores reais do parâmetro de recaptura
no.simulation <- 100000           # Número de simulações

# Configuração paralela
cl <- parallel::makeCluster(parallel::detectCores()/2)
parallel::clusterEvalQ(cl, {
  library(Rcpp)
  sourceCpp("gera_dados.cpp")  # Garante que o C++ está carregado nos workers
})
parallel::clusterExport(cl, c("N.true", "K.true", "c.true", "no.simulation",
                              "emv.Mt", "emv.Mtb", "Likelihood_function_Mt", 
                              "Likelihood_function_Mtb", "calc_metrics", 
                              "combine_i_results", "combine_j_results"))

# Execução principal
tictoc::tic()
resultados_list <- parallel::parLapply(cl, seq_along(K.true), process_k)
tictoc::toc()

parallel::stopCluster(cl)


matrix.N <- lapply(resultados_list, `[[`, "matrix.N")
matrix.c <- lapply(resultados_list, `[[`, "matrix.c")
AIC <- lapply(resultados_list, `[[`, "AIC")
BIC <- lapply(resultados_list, `[[`, "BIC")
AICc <- lapply(resultados_list, `[[`, "AICc")
HQIC <- lapply(resultados_list, `[[`, "HQIC")
CAIC <- lapply(resultados_list, `[[`, "CAIC")
TRV.90 <- lapply(resultados_list, `[[`, "TRV.90")
TRV.95 <- lapply(resultados_list, `[[`, "TRV.95")
TRV.99 <- lapply(resultados_list, `[[`, "TRV.99")
Vies.Mt.N <- lapply(resultados_list, `[[`, "Vies.Mt.N")
Vies.Mtb.N <- lapply(resultados_list, `[[`, "Vies.Mtb.N")
Vies.Mtb.c <- lapply(resultados_list, `[[`, "Vies.Mtb.c")
EQM.Mt.N <- lapply(resultados_list, `[[`, "EQM.Mt.N")
EQM.Mtb.N <- lapply(resultados_list, `[[`, "EQM.Mtb.N")
EQM.Mtb.c <- lapply(resultados_list, `[[`, "EQM.Mtb.c")


save(matrix.N, file = "dados/matrix.N.RData")
save(matrix.c, file = "dados/matrix.c.RData")
save(AIC, file = "dados/AIC.RData")
save(BIC, file = "dados/BIC.RData")
save(AICc, file = "dados/AICc.RData")
save(HQIC, file = "dados/HQIC.RData")
save(CAIC, file = "dados/CAIC.RData")
save(TRV.90, file = "dados/TRV.90.RData")
save(TRV.95, file = "dados/TRV.95.RData")
save(TRV.99, file = "dados/TRV.99.RData")
save(Vies.Mt.N, file = "dados/Vies.Mt.N.RData")
save(Vies.Mtb.N, file = "dados/Vies.Mtb.N.RData")
save(Vies.Mtb.c, file = "dados/Vies.Mtb.c.RData")
save(EQM.Mt.N, file = "dados/EQM.Mt.N.RData")
save(EQM.Mtb.N, file = "dados/EQM.Mtb.N.RData")
save(EQM.Mtb.c, file = "dados/EQM.Mtb.c.RData")
