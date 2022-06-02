#Este script se corre solamente para las replicas bootstrap. Pero es como hacer combinaciones entonces tendria que ir en un script diferente? 
#Yo creo que si. 

#Tengo que leer los inputs que en este caso son archivos cvs
# Son dos: matriz_conteo_rangos y Prob_rangos
#Las matrices de prob son son vlaores boot, las de conteo_rangos si son por valor boot

Matriz_Ns <- c(0, 10, 20, 30, 40, 50, 60)
Matriz_bootstrap <- c(seq(1, 10)) #No olvidar cambiar esto a 100 para todo el data set. 


Prob_rangos_ <- read.csv(paste("matriz_probabilidad_50Ns.csv", sep="")) #Para cambiar valor Ns tengo que cambiarlo
Prob_rangos_ <- Prob_rangos_[,-1]

loglikelihood_estimates <- matrix(ncol=7, nrow=10) #Aqui voy a pegar los resultados para graficarlos facilmente.
colnames(loglikelihood_estimates) <- c("0Ns", "10Ns", "20Ns", "30Ns", "40Ns", "50Ns", "60Ns") 

for (b in 1:length(Matriz_Ns))  {
  for (a in 1:length(Matriz_bootstrap))  {
    matriz_conteo_rangos <- read.csv(paste("matriz_", Matriz_Ns[b], "Ns_conteo_b", Matriz_bootstrap[a], ".csv", sep = ""))
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
    
    
    
    
    loglikelihood_estimates[a, b] <- sum(likelihoods) #sum(likelihoods) ## El segundo nombre cambia
    #Ahora looglikelihood lo tengo que meter a la matriz likelihood estiamtes
    
    print(sum(likelihoods)) #Hasta aqui todo funciona de maravilla
    
  }

}

write.csv(loglikelihood_estimates, paste("loglikelihood_estimates_50Ns.csv"))
jpeg("50Ns_plot", width=900, height=600)
boxplot(loglikelihood_estimates, main="Valor real 50Ns", xlab= "Valores 2Ns", ylab="Estimados de verosimilitud (-log)", col=4)
dev.off()

#El tesultado de esto es toda la matriz llena pero solo comparando contra un valor de matriz real. 
#Y en total son 7 de estas 