#Script para generar tabla de comparasiones de verosimilitud
#Contiene codigo para graficar.
#Hecho para comparar el efecto del numero de arboles en la calidad de las inferencias. 

replicas <- c(1,2,3)
valor_sel <- c("0.01", "0.1", "0", "1", "10", "100")
runs <- c(1, 10, 50, 100)

inference_matrix <- matrix(nrow=3, ncol=6)
colnames(inference_matrix) <- c("P(0.01Ns)", "P(0.1Ns)", "P(0Ns)", "P(1Ns)", "P(10Ns)", "P(100Ns)")

matriz_inferencia <- matrix(nrow=18, ncol=2)
colnames(matriz_inferencia) <- c("loglikelihood_estimate", "seleccion")

for (b in valor_sel)  {
  for (l in runs)  {
    p <- 1
    for (a in valor_sel)  { 
      for (f in replicas)  { 
        Prob_rangos_ <- read.csv(paste("matriz_probabilidad_", a, "Ns_100runs_set2_norm_smooth_1.csv", sep="")) #Para cambiar valor Ns tengo que cambiarlo
        Prob_rangos_ <- Prob_rangos_[,-1] #Borra la primer columna bugeada
        matriz_conteo_rangos <- read.csv(paste("matriz_conteo_rangos_", b, "Ns_", l, "runs_replica", f, ".csv", sep="")) #10 de esto
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

    write.csv(inference_matrix, paste("inferencia_parametros_smooth_norm_v2_", b, "Ns_", l, "runs.csv", sep=""))
  }  
  
}


### Visualización ###

valor_sel <- c("0.01", "0.1", "0", "1", "10", "100")
runs <- c(1, 10, 50, 100)

inferencia_Ns <- matrix(nrow=72, ncol=3)
colnames(inferencia_Ns) <- c("loglikelihood", "seleccion", "num_corridas")

for (a in valor_sel)  {
  i <- 1
  for (l in runs)  {
    tabla_inferencia <- read.csv(paste("inferencia_parametros_smooth_norm_v2_", a, "Ns_", l, "runs.csv", sep=""))
    tabla_inferencia <- tabla_inferencia[,-1]
    colnames(tabla_inferencia) <- c("0.01Ns", "0.1Ns", "0Ns", "1Ns", "10Ns", "100Ns")

    for (c in 1:ncol(tabla_inferencia))  {
      for (r in 1:nrow(tabla_inferencia))  {
        inferencia_Ns[i,1] <- tabla_inferencia[r,c]
        inferencia_Ns[i,2] <- colnames(tabla_inferencia[c])
        inferencia_Ns[i,3] <- l
        i <- i +1
      }
    }
  }

  write.csv(inferencia_Ns, paste("inferencia_", a, "Ns_.csv", sep=""))

}

inferencia_Ns <-read.csv("inferencia_0.1Ns.csv")
inferencia_Ns <- inferencia_Ns[,-1]

inferencia_Ns <- as.data.frame(inferencia_Ns)
inferencia_Ns$loglikelihood <- as.numeric(inferencia_Ns$loglikelihood)

grafica <- ggplot(data=inferencia_Ns, mapping = aes(x=seleccion, y=loglikelihood, colour=num_corridas))
         + geom_point(size=3, alpha=0.5) + facet_wrap(~num_corridas) + labs(title="Parametro real 100Ns")


grafica <- ggplot(data=inferencia_Ns, mapping = aes(x=seleccion, y=loglikelihood, colour=num_corridas)) + geom_point(size=3, alpha=0.5) + facet_wrap(~num_corridas) + labs(title="Parametro real 1Ns")
 
