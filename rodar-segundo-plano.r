# Inicia o processo em segundo plano e salva o PID:
system("Rscript teste.r & echo $! > pid.txt")

# Lê o PID do arquivo:
pid <- readLines("pid.txt")

# Mata o processo:
system(paste("kill", pid))

# Verifica se o processo foi encerrado:
system(paste("ps -p", pid))  # Deve retornar "Nenhum processo"