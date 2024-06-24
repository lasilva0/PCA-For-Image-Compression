library(jpeg)
library(gridExtra)
library(grid)
library(xtable)

imagem <- readJPEG("C:/Users/120895/Desktop/PCA Compression/rakan.jpeg", native = FALSE)
dir   <- "C:/Users/120895/Desktop/PCA Compression/"

# Função para calcular o EQM
calculate_eqm <- function(original, reconstructed) {
  return(mean((original - reconstructed)^2))
}


# Exiba as propriedades da imagem
cat("Dimensões da imagem:\n")
cat("Altura:", dim(imagem)[1], "\n")
cat("Largura:", dim(imagem)[2], "\n")
cat("Número de canais:", ifelse(length(dim(imagem)) == 3, dim(imagem)[3], 1), "\n")

# Visualize a imagem
plot(1:2, type='n')  # Cria uma área de plotagem vazia
rasterImage(imagem, 1, 1, 2, 2)  # Plota a imagem na área de plotagem

red   <- imagem[,,1]
green <- imagem[,,2]
blue  <- imagem[,,3]

r_svd <- svd(t(red)   %*% red)
g_svd <- svd(t(green) %*% green)
b_svd <- svd(t(blue)  %*% blue)

ncolunas <- c(1,10,50,100,500,1000)
rho   <- c()
gamma <- c()
eqm   <- c()
memory_vector <- c()

k = 1
for(i in ncolunas){
  r_projecao <- red %*% r_svd$u[,1:i]
  # reconstrução
  r_approx <- r_projecao %*% t(r_svd$u[,1:i])
  
  g_projecao <- green %*% g_svd$u[,1:i]
  # reconstrução
  g_approx <- g_projecao %*% t(g_svd$u[,1:i])
  
  b_projecao <- blue %*% b_svd$u[,1:i]
  # reconstrução
  b_approx <- b_projecao %*% t(b_svd$u[,1:i])
  
  imagemFinal <- array(NA, dim = c(nrow(red), ncol(red), 3))
  imagemFinal[,,1] <- r_approx
  imagemFinal[,,2] <- g_approx
  imagemFinal[,,3] <- b_approx
  
  imagemFinal <- ifelse(imagemFinal < 0, 0, ifelse(imagemFinal > 1,1,imagemFinal))
  
  # Cria um nome de arquivo único para cada imagem baseado no número de componentes
  nome_arquivo <- paste0("imagem_reconstruida_", i, ".png")
  
  # Caminho completo do arquivo de saída
  caminho_arquivo <- file.path(dir, nome_arquivo)
  
  # Salva a imagem reconstruída como um arquivo PNG
  png(file = caminho_arquivo)
  grid.newpage()  # Cria uma nova página de plotagem
  grid.raster(imagemFinal)
  dev.off()
  
  num_elements_original   <- prod(dim(red)) * 3 # 3 canais, por isso 3x
  num_elements_compressed <- prod(dim(r_projecao))    + 
                             prod(dim(r_svd$u[,1:i])) + 
                             prod(dim(r_svd$v[,1:i]))
  
  
  # Calcule o EQM para cada canal
  eqm_red   <- calculate_eqm(red,   r_approx)
  eqm_green <- calculate_eqm(green, g_approx)
  eqm_blue  <- calculate_eqm(blue,  b_approx)
  eqm_total <- (eqm_red + eqm_green + eqm_blue) / 3
  memory    <- i*dim(imagem)[2]
  
  compression_factor  <- (num_elements_compressed / num_elements_original)
  compression_rate    <- 1 - compression_factor 
  
  
  rho[k]   <- compression_factor
  gamma[k] <- 1 - rho[k]
  eqm[k]   <- eqm_total 
  memory_vector[k] <- memory
  
  k = k + 1
  
  cat("Número de componentes =", i, "\n", 
      "Fator de compressão = ", round(compression_factor,4), "\n",
      "EQM = ", round(eqm_total,4),"\n" ,
      "Taxa de compressão = ", round(compression_rate,4), "\n",
      "Memória necessária (dados finais) = ", memory, "\n")
}


dados <- data.frame(
  "Número de componentes" = ncolunas,
  "Fator de compressão" = rho,
  "Taxa de compressão" = gamma,
  "EQM" = eqm,
  "Memória necessária (MB)" = memory_vector
)

print(dados)


tabela_latex <- xtable(dados, caption = "Resultados da compressão de imagens utilizando PCA")

tabela_latex

# Imprimir a tabela LaTeX
print(tabela_latex, caption.placement = "top", digits = 4)


