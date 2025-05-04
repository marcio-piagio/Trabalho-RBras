#include <Rcpp.h>
#include <cuda_runtime.h>
#include <curand_kernel.h>
#include <Rcpp.h>
#include <random>

#ifdef __CUDACC__
#define RCPP_CUDA
#endif

// Modifique a inicialização do curand
curand_init((unsigned long long)clock() + idx, 0, 0, &state);  // Seed mais robusta

__global__ void gera_dados_kernel(double *dados, int N, int K, double *pj, double *cj, int num_sim) {
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (idx < num_sim) {
        curandState state;
        curand_init(clock64(), idx, 0, &state);

        // Variáveis locais para cada simulação
        int u[100], m[100], M[100], n[100]; // Ajuste o tamanho conforme K máximo
        double r;

        // Inicialização para i=0
        M[0] = 0;
        u[0] = 0;
        for(int j = 0; j < N; ++j) {
            if(curand_uniform_double(&state) < pj[0]) u[0]++;
        }
        m[0] = 0;
        n[0] = u[0] + m[0];

        // Loop para i >= 1
        for(int i = 1; i < K; ++i) {
            M[i] = M[i-1] + u[i-1];
            int remaining = N - M[i];
            
            u[i] = 0;
            for(int j = 0; j < remaining; ++j) {
                if(curand_uniform_double(&state) < pj[i]) u[i]++;
            }
            
            m[i] = 0;
            for(int j = 0; j < M[i]; ++j) {
                if(curand_uniform_double(&state) < cj[i]) m[i]++;
            }
            
            n[i] = u[i] + m[i];
        }

        r = M[K-1] + u[K-1];

        // Armazenamento dos resultados
        int base = idx * (4*K + 1);
        for(int i = 0; i < K; ++i) dados[base++] = u[i];
        for(int i = 0; i < K; ++i) dados[base++] = m[i];
        for(int i = 0; i < K; ++i) dados[base++] = M[i];
        for(int i = 0; i < K; ++i) dados[base++] = n[i];
        dados[base] = r;
    }
}



// [[Rcpp::export]]
Rcpp::NumericMatrix gera_dados_cpu(int N, int K, Rcpp::NumericVector pj, double c, int num_sim) {
    Rcpp::NumericMatrix resultados(num_sim, 4*K + 1);
    
    #pragma omp parallel for
    for(int idx = 0; idx < num_sim; ++idx) {
        std::mt19937 gen(idx);
        std::uniform_real_distribution<> dis(0.0, 1.0);
        
        int u[100], m[100], M[100], n[100];
        // ... (restante da lógica original)
    }
    return resultados;
}