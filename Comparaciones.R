#Script comparaciones maxima verosimilitud para una sola comparacion a la vez. 

Prob_rangos_ <- read.csv(paste("matriz_probabilidad_100kNe_ajuste4N.csv", sep="")) #Para cambiar valor Ns tengo que cambiarlo
Prob_rangos_ <- Prob_rangos_[,-1]

#loglikelihood_estimates <- matrix(ncol=1, nrow=10) #Aqui voy a pegar los resultados para graficarlos facilmente.
#colnames(loglikelihood_estimates) <- c("0Ns") 


matriz_conteo_rangos <- read.csv(paste("matriz_conteo_Relate.csv", sep=""))
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
       
loglikelihood_estimate <- sum(likelihoods) #sum(likelihoods) ## El segundo nombre cambia
#Ahora looglikelihood lo tengo que meter a la matriz likelihood estiamtes
    
print(sum(likelihoods)) #Hasta aqui todo funciona de maravilla
    





#write.csv(loglikelihood_estimates, paste("loglikelihood_estimates_50Ns.csv"))
#jpeg("50Ns_plot", width=900, height=600)
#boxplot(loglikelihood_estimates, main="Valor real 50Ns", xlab= "Valores 2Ns", ylab="Estimados de verosimilitud (-log)", col=4)
#dev.off()


