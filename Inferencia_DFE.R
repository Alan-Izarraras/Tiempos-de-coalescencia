
### Codigo para generar estimados de verosimilitud de parametros gamma a partir de arboles simulados con una DFE conocida ###
#Generar tabla de parametros gamma.
#Si quieres usar más corridas tienes que modificar algunas cosas... luego las comento o corrijo para que se haga automatico#

#Reemplazo '2.5' -> 1
#Reemplazo 5 -> 2

UpperThreshold = 60
Limit = UpperThreshold + 2.5
Runs = 60/5
RunsPlusTwo = Runs + 3

#### These two lists define the grid of the alpha and gamma parameters
AlphaGrid <- 0.02*1:10
GammaGrid <- 2*1:10

TwoNsValues <- 0:Runs * 5 + 2.5

Table <- matrix(ncol=RunsPlusTwo,nrow=0)

for (j in AlphaGrid){
    for (k in GammaGrid){
		Probability <- 0
		Row <- c(j,k)
		for (i in TwoNsValues){
			
			if ( i == 2.5){
				Probability <- pgamma(2.5,j,1/k)
				Row <- c(Row,Probability)
			}else if (i==Limit){
				Probability <- ( 1 - pgamma(i-5,j,1/k) )
				Row <- c(Row,Probability)
			}else{
				Probability <-(pgamma(i,j,1/k) - pgamma(i-5,j,1/k))
				Row <- c(Row,Probability)
			}
			
		}
		Table <- rbind(Table,Row)
#print (Probability)
    }
}

colnames(Table) <- c("alpha", "gamma", "P(2Ns < 2.5)", "P(2.5 < 2Ns < 7.5)", "P(7.5 < 2Ns < 12.5)", "P(12.5 < 2Ns < 17.5)", "P(17.5 < 2Ns < 22.5)", "P(22.5 < 2Ns < 27.5)", "P(27.5 < 2Ns < 32.5)", "P(32.5 < 2Ns < 37.5)", "P(37.5 < 2Ns < 42.5)", "P(42.5 < 2Ns < 47.5)", "P(47.5 < 2Ns < 52.5)", "P(52.5 < 2Ns < 57.5)", "P( 2Ns > 57.5)")

write.table(Table,file="TableOfProbabilities.txt",row.names=FALSE,col.names=FALSE,sep="\t")


#Leo las matrices de probabilidad en formato csv y las paso a una lista

Lista_matrices <- list()
Vector_Ns <- seq(0, 60, 5)

for (v in 1:length(Vector_Ns))  {
  uno <- read.csv(paste("matriz_probabilidad_", Vector_Ns[v], "Ns.csv", sep=""))
  uno <- uno[,-1]
  Lista_matrices[[v]] <- uno
}

#Agarro esa lista de matrices y las multiplico por el primer parametro gamma. 
#Donde cada intervalo corresponde aun valor Ns fijo. 
#Agrego un pasito para pasar de probabilidad a conteo. Le pido el valor de 20,4 y divido la matriz entre eso. 
#Si para pasar de conteo a probabilidad tengo que multiplicar por el valor maximo entonces para pasar de probabilidad a conteo
#Tengo que dividir por el valor maximo. Pero aqui no tengo conocimiento del valor maximo. 

Lista_Suma <- list()

for (r in 1:nrow(Table))  {
  i <- 3
  for (v in 1:length(Vector_Ns))  {
    Lista_Suma[[v]] <- Lista_matrices[[v]]*Table[r,i]
    i <- i+1
  }
  matriz_gamma_suma <- Lista_Suma[[1]] + Lista_Suma[[2]] + Lista_Suma[[3]] + Lista_Suma[[4]] + Lista_Suma[[5]] + Lista_Suma[[6]] + Lista_Suma[[7]] + Lista_Suma[[8]] + Lista_Suma[[9]] + Lista_Suma[[10]] + Lista_Suma[[11]] + Lista_Suma[[12]] + Lista_Suma[[13]]
  write.csv(matriz_gamma_suma, paste("matriz_gamma_", r, ".csv", sep=""))
}


#Ah ok me equivoque tantito porque para ahcer las comparaciones necesito los conteos y no las probabilidades. 
#Como paso de probabilidad a conteo? solamente divido la matriz entre en n max. 

#La idea es guardar las matrices de conteo en una lista de matrices
#Lo mismo con las likelihoods y likelihood_estimates

loglikelihoods <- vector()
likelihood_estimates <- matrix(nrow=100, ncol=3)
lista_DFE_conteo <- list()
lista_loglikelihoods <- list()
lista_likelihood_estimates <- list()

for (x in 1:10)  { #Por ahora estoy usando 10 bootstraps
  DFE_conteo <- read.csv(paste("Matriz_DFEboot", x, "_conteo.csv", sep=""))
  DFE_conteo <- DFE_conteo[,-1]
  lista_DFE_conteo[[x]] <- DFE_conteo

  for (a in 1:100)  {
    matriz_gamma <- read.csv(paste("matriz_gamma_", a, ".csv", sep=""))
    matriz_gamma <- matriz_gamma[,-1]
    i<-1
    for (r in 1:nrow(DFE_conteo))  {
      for (c in 1:ncol(DFE_conteo))  {
        #DFE_conteo <- lista_DFE_conteo[x]  #Testear esta linea
        if (DFE_conteo[r, c] != 0)  {
          loglikelihoods[i] <- DFE_conteo[r,c] * log(matriz_gamma[r, c]) #Necesito que este archivo no cambie
          lista_loglikelihoods[[x]] <- loglikelihoods
          i <- i +1
        }
      }
    }
    likelihood_estimates[a,3] <- sum(loglikelihoods)
    lista_likelihood_estimates[[x]] <- likelihood_estimates
  }
}

#Lo de abaho es meramente para graficar.

vector_maximos <- vector()

#Codigo para pasar alfa y theta de cada estimado.
for (q in 1:10)  {
  likelihood_estimates <- lista_likelihood_estimates[[q]]
  for (x in 1:100) {
    likelihood_estimates[x,1] <- Table[x,1] 
    likelihood_estimates[x,2] <- Table[x,2]
    #lista_likelihood_estimates[[q]] <- likelihood_estimates
  }
  likelihood_estimates <- cbind(likelihood_estimates, seq(1:100))
  likelihood_estimates <- as.data.frame(likelihood_estimates)
  #lista_likelihood_estimates[[q]] <- likelihood_estimates
  estimado_maximo <- order(likelihood_estimates$V3, decreasing = T)
  estimado_maximo[1]
  vector_maximos[q] <- estimado_maximo[1]
  jpeg(filename = paste("DFE_bootstrap", q, ".jpeg", sep=""), width=900)
  stripchart(likelihood_estimates$V3~likelihood_estimates$V4)
  stripchart(likelihood_estimates$V3~likelihood_estimates$V4, vertical=T, main="Curva de verosimilitud", ylab = "loglikelihood values", xlab = "Gamma distribution parameter combinations", col=1 , pch = 20)
  abline(v=95, col=2, pch = 7, cex = 2, lwd = 2, lty = 5) #Valor real Rojo = real
  abline(v=estimado_maximo[1], col = 3, pch = 7, cex = 2, lwd = 2, lty = 5) #Valor estimado Verde = Estimado
  dev.off()  
}

#Bien, con eso ya tengo todas las graficas individuales. Ahora falta la grafica que incorpora todos los bootstraps.

#2D scatterplot
#Falta agregar codigo para construir la matriz maximos

maximos <- matrix(nrow=10, ncol=2)

for (i in 1:10)  {
  indice <- vector_maximos[i] 
  maximos[i,1] <- likelihood_estimates[indice,1] 
  maximos[i,2] <- likelihood_estimates[indice,2]  
}

jpeg(filename = "Inferencia de parametros gamma.jpeg", width=500)
plot(maximos, main = "Inferencia de parametros de una distribucion gamma", xlim = c(0.02, 0.20), ylim = c(2,20), xlab = "Valores Alpha", ylab = "Valores Gamma", pch = 5)
abline(v=0.20, h=10, col=2, pch = 7, cex = 2, lwd = 2, lty = 5)
points(y = 10, x= 0.20, pch = 20, col=2)
dev.off()


#Metodo alternativo de plottear, Stripchart fallido
#stripchart(likelihood_estimates$V1~likelihood_estimates$V2, main="Estimacion de parametros de una distribución gamma",ylab = "Parámetro alfa", xlab = "Parametro Beta", vertical=T, method="jitter", pch=19, col=F)
#abline(v=5, col=2, pch = 7, cex = 2, lwd = 2, lty = 5)
#abline(h=0.20, col=2, pch = 7, cex = 2, lwd = 2, lty = 5)
#points(y = 0.20, x= 6)

