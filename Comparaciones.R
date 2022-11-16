

valor_ns <- c(0, 0.028, 0.28, 2.8)
loglike_matrix <- matrix(nrow = 10, ncol=length(valor_ns))
colnames(loglike_matrix) <- c("0Ns", "0.028Ns", "0.28Ns", "2.8Ns")
x <- 0


for (a in valor_ns)  {
  x <- x + 1
  for (n in 1:10)  {
    Prob_rangos_ <- read.csv(paste("matriz_probabilidad_", a, "Ns_sim", n, ".csv", sep="")) #Para cambiar valor Ns tengo que cambiarlo
    Prob_rangos_ <- Prob_rangos_[,-1]
    Prob_rangos_ <- Prob_rangos_[,-c(8,9,10,11,12)]

    matriz_conteo_rangos <- read.csv(paste("matriz_conteo_Relate.csv", sep=""))
    matriz_conteo_rangos <- matriz_conteo_rangos[,-1]
    matriz_conteo_rangos <- matriz_conteo_rangos[,-c(8,9,10,11,12)]

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
    
    loglikelihood_estimate <- sum(likelihoods) #sum(likelihoods) ## El segundo nombre cambia
    #Ahora looglikelihood lo tengo que meter a la matriz likelihood estiamtes
    
    loglike_matrix[n, x] <- loglikelihood_estimate

    #Matriz donde cada fila es una repeticion y cada columna un vlaor 2Ns
  }
}


#Ahora para plottear.

colnames(loglike_matrix) <- c("0", "0.028", "0.28", "2.8")
rownames(loglike_matrix) <- 1:10 

loglike_matrix <- as.data.frame(loglike_matrix)

#Lo transformo a una columna de valores log y otra columna de valores 2Ns para graficar facilmente. 

loglike_matrix_plot <- matrix(nrow=40, ncol=2)
loglike_matrix_plot[1,n]
loglike_matrix_plot[2,n] <- colnames

i <- 1
for (c in 1:ncol(loglike_matrix))  { 
  for (r in 1:nrow(loglike_matrix))  {   
    z <- loglike_matrix[r, c]
    loglike_matrix_plot[i, 2] <- colnames(loglike_matrix[c])
    loglike_matrix_plot[i, 1] <- z
    i <- i + 1 
  }
}

loglike_matrix_plot <- as.data.frame(loglike_matrix_plot)

#Plot. Tengo que asegurarme que datos y grupos sean valores numericos y no se transformen a char. 
datos <- as.numeric(loglikelihood_estimate$V1)
grupos <- as.numeric(loglikelihood_estimate$V2) 


stripchart(datos ~ grupos, group.names = c("0Ns","0.028Ns","0.28Ns","2.8Ns"), vertical = T, pch = 20, method = "jitter", col = 2, las = 2)




#write.csv(loglikelihood_estimates, paste("loglikelihood_estimates_50Ns.csv"))
#jpeg("50Ns_plot", width=900, height=600)
#boxplot(loglikelihood_estimates, main="Valor real 50Ns", xlab= "Valores 2Ns", ylab="Estimados de verosimilitud (-log)", col=4)
#dev.off()

