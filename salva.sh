#!/bin/bash

# Definir variáveis
R_SCRIPT="simulate-melhorado.R"  # Nome do script R
OUTPUT_DIR="output_$(date +%Y-%m-%d_%H-%M-%S)"  # Diretório com timestamp

# Criar diretório de saída
mkdir -p "$OUTPUT_DIR"

# Executar o script R
Rscript "$R_SCRIPT"

# Mover arquivos específicos para o diretório de saída
mv matrix.N.RData matrix.c.RData resultados.RData "$OUTPUT_DIR"

echo "Arquivos salvos em $OUTPUT_DIR"
