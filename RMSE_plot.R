
### RMSE ver1 usando valores 2Ns ###

matriz_maximos <- read.csv("matriz_maximos_100runs.csv")
matriz_maximos <- matriz_maximos[,-1]

valor_ns <- c(0,0.001, 0.003 ,0.005 ,0.01 ,0.017 ,0.03 ,0.05 ,0.1 ,0.17 ,0.31 ,0.56 ,1,1.77 ,3.16 ,5.62 ,10 ,17 ,31,56,100 ,177 ,316,562,1000,1778)

i <- 1
v <- 1
for (r in 1:nrow(matriz_maximos_grupo4))  {
	matriz_maximos_grupo4[r,2] <- valor_ns[v]
	i <- i +1
	matriz_maximos_grupo4[r,1] <- gsub("Ns", "", matriz_maximos_grupo4[r,1])
	if (i==50)  {
	  v <- v + 1
	  i <- 0
	}
}

#matriz_maximos_grupo2[1300,2] <- v
matriz_maximos_grupo4[,1] <- as.numeric(matriz_maximos_grupo4[,1])

sqrt(mean((matriz_maximos_grupo4$parametro_real - matriz_maximos_grupo4$parametro_estimado )^2))

### RMSE plot ###
library(ggplot)

p <- ggplot(data=RMSE, aes(x=theta, y=RMSE, color=theta, fill=theta)) +
+     geom_bar(stat="identity")


#Para el caso en el que tengo varias agrupaciones que quiero probar dentro de una misma matriz. 

matriz_maximos_grupo1 <- matriz_maximos[which(matriz_maximos[,"agrupacion"]==1), ] 
matriz_maximos_grupo2 <- matriz_maximos[which(matriz_maximos[,"agrupacion"]==2), ]
matriz_maximos_grupo3 <- matriz_maximos[which(matriz_maximos[,"agrupacion"]==3), ]
matriz_maximos_grupo4 <- matriz_maximos[which(matriz_maximos[,"agrupacion"]==4), ]   
