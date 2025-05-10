library(Rcpp)

# Função R
gera_dados <- function(N, K, pj, c) {
  u <- m <- n <- M <- numeric(K)
  etaj <- log(pj / (1 - pj))
  cj <- exp(etaj + c) / (1 + exp(etaj + c))
  u[1] <- rbinom(1, N, pj[1])
  m[1] <- 0
  n[1] <- u[1] + m[1]
  M[1] <- 0
  for (i in 2:K) {
    M[i] <- M[i-1] + u[i-1]
    u[i] <- rbinom(1, N - M[i], pj[i])
    m[i] <- rbinom(1, M[i], cj[i])
    n[i] <- u[i] + m[i]
  }
  r <- M[K] + u[K]
  return(c(u, m, M, n, r))
}

# Compilar função C++
Rcpp::sourceCpp(code = '
#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericVector gera_dados_cpp(int N, int K, NumericVector pj, double c) {
  NumericVector u(K), m(K), M(K), n(K);
  
  NumericVector etaj = log(pj / (1 - pj));
  NumericVector cj = exp(etaj + c) / (1 + exp(etaj + c));
  
  u[0] = R::rbinom(N, pj[0]);
  m[0] = 0;
  n[0] = u[0] + m[0];
  M[0] = 0;
  
  for (int i = 1; i < K; i++) {
    M[i] = M[i-1] + u[i-1];
    u[i] = R::rbinom(N - M[i], pj[i]);
    m[i] = R::rbinom(M[i], cj[i]);
    n[i] = u[i] + m[i];
  }
  
  double r = M[K-1] + u[K-1];
  
  // Concatenação manual dos vetores (u, m, M, n, r)
  int total_size = u.size() + m.size() + M.size() + n.size() + 1;
  NumericVector result(total_size);
  
  int pos = 0;
  // Copia u
  for (int i = 0; i < u.size(); i++) result[pos++] = u[i];
  // Copia m
  for (int i = 0; i < m.size(); i++) result[pos++] = m[i];
  // Copia M
  for (int i = 0; i < M.size(); i++) result[pos++] = M[i];
  // Copia n
  for (int i = 0; i < n.size(); i++) result[pos++] = n[i];
  // Adiciona r
  result[pos] = r;
  
  return result;
}')

# Verificar a Equivalência dos Resultados

# Definir parâmetros
N <- 100   # Tamanho da população
K <- 11     # Número de amostragens
pj <- runif(K, 0.1, 0.9)  # Probabilidades de captura
c <- -0.5    # Efeito de recaptura

# Gerar resultados com semente fixa
set.seed(666)
(result_r <- gera_dados(N, K, pj, c))

set.seed(666)
(result_cpp <- gera_dados_cpp(N, K, pj, c))

# Verificar igualdade
identical(result_r, result_cpp)  # Deve retornar TRUE

# Comparar o Tempo de Execução

library(microbenchmark)

# Definir parâmetros de teste
N <- 100  # População grande para destacar diferenças
K <- 20     # Amostragens múltiplas

# Benchmarking
bench <- microbenchmark(
  R = gera_dados(N, K, pj, c),
  Rcpp = gera_dados_cpp(N, K, pj, c),
  times = 10000  # Número de repetições
)

# Resultados
print(bench)

# Comparar o Uso de Memória

library(peakRAM)

# Medir uso de memória
mem_usage <- peakRAM(
  R = gera_dados(N, K, pj, c),
  Rcpp = gera_dados_cpp(N, K, pj, c)
)

print(mem_usage)


