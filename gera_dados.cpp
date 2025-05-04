#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericVector gera_dados_cpp(int N, int K, NumericVector pj, double c) {
  NumericVector u(K);
  NumericVector m(K);
  NumericVector M(K);
  NumericVector n(K);
  
  NumericVector etaj = log(pj / (1 - pj));
  NumericVector cj = exp(etaj+c)/(1+exp(etaj+c));
  
  u[0] = R::rbinom(N, pj[0]);
  m[0] = 0;
  n[0] = u[0] + m[0];
  M[0] = 0;
  
  for(int i = 1; i < K; i++) {
    M[i] = M[i-1] + u[i-1];
    u[i] = R::rbinom(N - M[i], pj[i]);
    m[i] = R::rbinom(M[i], cj[i]);
    n[i] = u[i] + m[i];
  }
  
  double r = M[K-1] + u[K-1];
  
  int total_size = u.size() + m.size() + M.size() + n.size() + 1;
  NumericVector result(total_size);
  
  int pos = 0;
  // Copia u
  for(int i = 0; i < u.size(); i++) {
    result[pos++] = u[i];
  }
  // Copia m
  for(int i = 0; i < m.size(); i++) {
    result[pos++] = m[i];
  }
  // Copia M
  for(int i = 0; i < M.size(); i++) {
    result[pos++] = M[i];
  }
  // Copia n
  for(int i = 0; i < n.size(); i++) {
    result[pos++] = n[i];
  }
  // Adiciona r no final
  result[pos] = r;
  
  return result;
}