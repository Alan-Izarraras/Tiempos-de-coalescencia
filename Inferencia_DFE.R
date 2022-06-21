
### Codigo para generar estimados de verosimilitud de parametros gamma a partir de arboles simulados con una DFE conocida ###
#Generar tabla de parametros gamma.

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

DFE_conteo <- read.csv("matriz_DFE1_conteo.csv", sep=",")
DFE_conteo <- DFE_conteo[,-1]
loglikelihoods <- vector()
likelihood_estimates <- matrix(nrow=100, ncol=3)

for (a in 1:100)  {
  matriz_gamma <- read.csv(paste("matriz_gamma_", a, ".csv", sep=""))
  matriz_gamma <- matriz_gamma[,-1]
  i<-1
  for (r in 1:nrow(DFE_conteo))  {
    for (c in 1:ncol(DFE_conteo))  {
      if (DFE_conteo[r, c] != 0)  {
        loglikelihoods[i] <- DFE_conteo[r,c] * log(matriz_gamma[r, c]) #Necesito que este archivo no cambie
        i <- i +1
      }
    }
  }
  likelihood_estimates[a,3] <- sum(loglikelihoods)
}


#Codigo para pasar alfa y theta de cada estimado.
for (x in 1:100) {
  likelihood_estimates[x,1] <- Table[x,1] 
  likelihood_estimates[x,2]	<- Table[x,2]
}

#Que no cunda el panico. 

likelihood_estimates <- cbind(likelihood_estimates, seq(1:100))
likelihood_estimates <- as.data.frame(likelihood_estimates)

stripchart(likelihood_estimates$V3~likelihood_estimates$V4)
stripchart(likelihood_estimates$V3~likelihood_estimates$V4, vertical=T, main="Curva de verosimilitud", ylab = "loglikelihood values", xlab = "Gamma distribution parameter combinations", col=2 , pch = 16,)

#Necesito eso pero sin los pountos, solo poninendo el punto del estimado y despues el punto real.

stripchart(likelihood_estimates$V1~likelihood_estimates$V2, main="Parametros de una distribución gamma",ylab = "Parámetro alfa", xlab = "Parametro Beta", vertical=T, method="jitter", pch=19, col=3)


#Falta hacer grafica de barras con las porporciones de la DFE real y la DFE estimada.

Table_finalistas <- matrix()
Table_finalistas <- rbind(Table[95,], Table[96,])
Table_finalistas <- Table_finalistas[,-1]
Table_finalistas <- Table_finalistas[,-1]
#Ahora esto tengo que... quitarle las priemras dos columnas y luego voltear las dimensiones.

Data_table <- matrix(nrow=13, ncol=3)
categorias <- colnames(Table_finalistas)

i <- 1
for (x in 1:13)  {
	Data_table[x,1] <- Table_finalistas[1, x]
	Data_table[x,2] <- Table_finalistas[2, x]
	Data_table[x,3] <- categorias[x]
	i <- i+1
}

Data_table <- as.data.frame(Data_table)
Data_table$V1 <- as.numeric(Data_table$V1)
Data_table$V2 <- as.numeric(Data_table$V2)

barplot(height = Data_table$V1, names=Data_table$V3)

 c("2Ns < 2.5", "2.5 < 2Ns < 7.5", "7.5 < 2Ns < 12.5", "12.5 < 2Ns < 17.5", "17.5 < 2Ns < 22.5", "22.5 < 2Ns < 27.5", "27.5 < 2Ns < 32.5", "32.5 < 2Ns < 37.5", "37.5 < 2Ns < 42.5", "42.5 < 2Ns < 47.5", "47.5 < 2Ns < 52.5", "52.5 < 2Ns < 57.5", "2Ns > 57.5")

barplot(height = Data_table$V1, names=Data_table$V3, col = 4, ylab="Probabilidad", xlab = "Rango 2Ns", main = "Distribución gamma", ylim = seq(0,1), width = 0.1, las=2, names.arg = c("2Ns < 2.5", "2.5 < 2Ns < 7.5", "7.5 < 2Ns < 12.5", "12.5 < 2Ns < 17.5", "17.5 < 2Ns < 22.5", "22.5 < 2Ns < 27.5", "27.5 < 2Ns < 32.5", "32.5 < 2Ns < 37.5", "37.5 < 2Ns < 42.5", "42.5 < 2Ns < 47.5", "47.5 < 2Ns < 52.5", "52.5 < 2Ns < 57.5", "2Ns > 57.5"))
