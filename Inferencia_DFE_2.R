#Inferencia DFE para valores Ns = seq(0,60,5)

UpperThreshold = 60 
Limit = UpperThreshold + 2.5 
Runs = 60/5 
RunsPlusTwo = Runs + 3 

#### These two lists define the grid of the alpha and gamma parameters
AlphaGrid <- 0.02*1:10 
GammaGrid <- 2*1:10 

TwoNsValues <- 0:Runs * 5 + 2.5   

Table <- matrix(ncol=RunsPlusTwo,nrow=0)  

for (j in AlphaGrid)  {
  for (k in GammaGrid)  {
    Probability <- 0
    Row <- c(j,k)
    for (i in TwoNsValues)  {
      if ( i == 2.5)  { #Cambio
        Probability <- pgamma(2.5,j,1/k) 
        Row <- c(Row,Probability)
      }else if (i==Limit)  {
        Probability <- (1 - pgamma(i-5,j,1/k)) 
        Row <- c(Row,Probability)
      }else  {
        Probability <- (pgamma(i,j,1/k) - pgamma(i-5,j,1/k)) 
        Row <- c(Row,Probability)
      }
      
    }
    Table <- rbind(Table,Row)
#print (Probability)
  }
}

write.table(Table,file="TableOfProbabilities.txt",row.names=FALSE,col.names=FALSE,sep="\t")

Lista_matrices <- list()
Vector_Ns <- seq(0, 60, 5) #Aqui tambien cambia

for (v in 1:length(Vector_Ns))  {
  uno <- read.csv(paste("matriz_probabilidad_", Vector_Ns[v], "Ns.csv", sep=""))
  uno <- uno[,-1]
  Lista_matrices[[v]] <- uno
}

Lista_Suma <- list()

for (r in 1:nrow(Table))  {
  i <- 3
  for (v in 1:length(Vector_Ns))  {
    Lista_Suma[[v]] <- Lista_matrices[[v]]*Table[r,i]
    i <- i+1
  } #Aqui cambia
  matriz_gamma_suma <- Lista_Suma[[1]] + Lista_Suma[[2]] + Lista_Suma[[3]] + Lista_Suma[[4]] + Lista_Suma[[5]] + Lista_Suma[[6]] + Lista_Suma[[7]] + Lista_Suma[[8]] + Lista_Suma[[9]] + Lista_Suma[[10]] + Lista_Suma[[11]] + Lista_Suma[[12]] + Lista_Suma[[13]] 
  write.csv(matriz_gamma_suma, paste("matriz_gamma_", r, ".csv", sep=""))
}

loglikelihoods <- vector()
likelihood_estimates <- matrix(nrow=100, ncol=3)
lista_DFE_conteo <- list()
lista_loglikelihoods <- list()
lista_likelihood_estimates <- list()

for (x in 1:10)  { #Por ahora estoy usando 10 bootstraps
  DFE_conteo <- read.csv(paste("matriz_DFEboot", x, "_conteo.csv", sep=""))
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



#vector_maximos <- vector_maximos[!duplicated(vector_maximos)]
intervalo <- matrix(nrow=length(vector_maximos), ncol=2)
i <- 1
#vector_maximos <- vector_maximos -1 #Tengo que hcer esto porque el vector maximos esta basando en la tabla de likelihoods que tiene fila 1, Pero estoy oprenaod sobre la tabla de probabilidades que tiene fila 0. 
for (n in vector_maximos)  {
  intervalo[i,1] <- Table[n,3] 
  intervalo[i,2] <- sum(Table[n,4], Table[n,5], Table[n,6], Table[n,7], Table[n,8], Table[n,9], Table[n,10], Table[n,11], Table[n,12], Table[n,13], Table[n,14], Table[n,15])
  i <- i + 1
}

vector_intervalos <- c("(0.20, 12)", "(0.20, 14)")
jpeg(filename = "Comparaciones_Intervalo_1.jpeg", width=600)
barplot(height = intervalo[,2], main= "Comparison of P(2Ns >= 1)", names=vector_intervalos, ylim=c(0,1), ylab = "Probability", col= 3)
dev.off()
jpeg(filename = "Comparaciones_Intervalo_2.jpeg", width=600)
barplot(height = intervalo[,1], main= "Comparison of P(1 < 2Ns < 61)", names=vector_intervalos, ylim=c(0,1), ylab = "Probability", col= 3)
dev.off()

rownames(intervalo) <- vector_intervalos

write.table(intervalo,file="ComparacionDistribuciones.txt",row.names=T,col.names=FALSE,sep="\t")

#Agregar el valor real a la matriz intervalo, nombrar las columnas y pasarlo a un cvs. 





