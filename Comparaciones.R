### Agregar probabilidad minima a matrics de probabilidad. ya funciona. 
### Esto ya está. Solo falta tratar distintas thetas y ver que onda con 17. 


valor_sel <- c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "20", "21", "22", "23", "24", "25", "26")
rangos <- (1,2,3,4)

for (s in valor_sel)  {
  matriz_prob <- read.csv(paste("matriz_probabilidad_", s, "Ns_100runs_set2_rangos_1.csv", sep=""))	
  matriz_prob <- matriz_prob[,-1]
  matriz_cuenta <- matriz_prob + 0.00000001

  tabla_referencia <- matriz_cuenta[-10,] #Reemplazar 11 por ultima fila

  for (i in 1:8)  {
	 matriz_cuenta[10,i] <- 1 - sum(tabla_referencia[,i])
  }
  write.csv(matriz_cuenta, paste("matriz_probabilidad_", s, "Ns_100runs_rangos1_chido.csv", sep=""))
}



### #Adaptación a distinta agrupación de linajes (rangos)

replicas <- seq(1,49)
valor_sel <- c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "20", "21", "22", "23", "24", "25", "26")
rangos <- c(1,2,3,4)


inference_matrix <- matrix(nrow=49, ncol=26) #replicas * valors 2Ns 
colnames(inference_matrix) <- c("P(0Ns)", "P(0.001Ns)", "P(0.003Ns)", "P(0.005Ns)", "P(0.01Ns)", "P(0.017Ns)", "P(0.03Ns)", "P(0.05Ns)", "P(0.1Ns)", "P(0.17Ns)", "P(0.31Ns)", "P(0.56Ns)", "P(1Ns)", "P(1.77Ns)", "P(3.16Ns)", "P(5.62Ns)", "P(10Ns)", "P(17Ns)", "P(31Ns)", "P(56Ns)", "P(100Ns)", "P(177Ns)", "P(316Ns)", "P(562Ns)", "P(1000Ns)", "P(1778Ns)")

matriz_inferencia <- matrix(nrow=1274, ncol=2) #26x49x4 = 5200 
colnames(matriz_inferencia) <- c("loglikelihood_estimate", "seleccion")

for (l in rangos)  {  
  for (b in valor_sel)  {
    p <- 1
    for (a in valor_sel)  {
      for (f in replicas)  { 
        Prob_rangos_ <- read.csv(paste("matriz_probabilidad_", a, "Ns_100runs_rangos", l, "_chido.csv", sep="")) #Para cambiar valor Ns tengo que cambiarlo
        Prob_rangos_ <- Prob_rangos_[,-1] #Borra la primer columna bugeada
        matriz_conteo_rangos <- read.csv(paste("matriz_conteo_smooth_", b, "Ns_100runs_rango", l, "_replica", f, ".csv", sep="")) #10 de esto
        matriz_conteo_rangos <- matriz_conteo_rangos[,-1]          
        likelihoods <- vector()
        i <- 1
        for (r in 1:nrow(matriz_conteo_rangos))  {
          for (c in 1:ncol(matriz_conteo_rangos))  {
            if (matriz_conteo_rangos[r, c] != 0)  {
              likelihoods[i] <- matriz_conteo_rangos[r,c] * log(Prob_rangos_[r, c]) #Necesito que este archivo no cambie
              i <- i +1
            }
          }
        }
        loglikelihood_estimate <- sum(likelihoods) 
        inference_matrix[f, p] <- loglikelihood_estimate
      }
      p <- p +1 
    } 
    write.csv(inference_matrix, paste("inferencia_parametros_", b, "Ns_100runs_rango", l, ".csv", sep=""))
  }
}  



#Si funcionó, deberían de haber... 26 archivos de inferencia... y todo completo 104
#Ahora a la parte de abajo.

### Visualización parametros maximos ###

valor_sel <- c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "20", "21", "22", "23", "24", "25", "26")
rangos <- c(1,2,3,4)

matriz_maximos <- matrix(nrow=5096, ncol =3) #nrow = sel*replicas*rangos = 
colnames(matriz_maximos) <- c("parametro_estimado", "parametro_real", "agrupacion")

i <- 1
for (a in valor_sel)  {
  for (l in rangos)  {
    tabla_inferencia <- read.csv(paste("inferencia_parametros_", a, "Ns_100runs_rango", l, ".csv", sep=""))
    tabla_inferencia <- tabla_inferencia[,-1]
    colnames(tabla_inferencia) <- c("0Ns", "0.001Ns", "0.003Ns", "0.005Ns", "0.01Ns", "0.017Ns", "0.03Ns", "0.05Ns", "0.1Ns", "0.17Ns", "0.31Ns", "0.56Ns", "1Ns", "1.77Ns", "3.16Ns", "5.62Ns", "10Ns", "17Ns", "31Ns", "56Ns", "100Ns", "177Ns", "316Ns", "562Ns", "1000Ns", "1778Ns")

    for (f in 1:nrow(tabla_inferencia))  {
      valor_maximo <- max(tabla_inferencia[f,])
      parametro_maximo <- which(tabla_inferencia[f,] == valor_maximo)
      matriz_maximos[i,1] <- colnames(tabla_inferencia[parametro_maximo[1]])
      matriz_maximos[i,2] <- a
      matriz_maximos[i,3] <- l
      i <- i +1
    }   
  }
}

matriz_maximos <- as.data.frame(matriz_maximos)
write.csv(matriz_maximos, "matriz_maximos_rangos")

ggplot(data = matriz_maximos, mapping = aes(x=parametro_real, y=parametro_estimado, color=agrupacion)) + 
geom_jitter(width=0.2, height=0.2, size=2, alpha=0.5) + facet_wrap(~agrupacion, scales="free") + scale_x_discrete(limit = c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "20", "21", "22", "23", "24", "25", "26")
 + scale_y_discrete(limit=c("0Ns", "0.001Ns", "0.003Ns", "0.005Ns", "0.01Ns", "0.017Ns", "0.03Ns", "0.05Ns", "0.1Ns", "0.17Ns", "0.31Ns", "0.56Ns", "1Ns", "1.77Ns", "3.16Ns", "5.62Ns", "10Ns", "17Ns", "31Ns", "56Ns", "100Ns", "177Ns", "316Ns", "562Ns", "1000Ns", "1778Ns")
) + labs(title="Estimados de máxima verosimilitud por agrupacion")

### Version bubbleplot solo para valores relevantes! 
ggplot(data = matriz_maximos) + 
  geom_count(mapping = aes(x=parametro_real, y=parametro_estimado), alpha=0.5, color=6) + 
  scale_size_area(max_size=7) + 
  scale_y_discrete(limit=c("0", "0.001", "0.003", "0.005", "0.01", "0.017", "0.03", "0.05", "0.1", "0.17", "0.31", "0.56", "1", "1.77", "3.16", "5.62", "10", "17", "31", "56", "100", "177", "316", "562", "1000", "1778"), breaks=c("0", "0.01", "0.1", "1", "10", "100", "1000")) + 
  scale_x_discrete(limits=c("0", "0.01", "0.1", "1", "10", "100", "1000"), breaks=c("0", "0.01", "0.1", "1", "10", "100", "1000"))


