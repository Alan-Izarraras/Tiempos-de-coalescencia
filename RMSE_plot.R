
### RMSE ver1 usando valores 2Ns ###

Necesito los valores 2 los tengo en el script de alado

matriz_maximos <- read.csv("matriz_maximos_100runs.csv")
matriz_maximos <- matriz_maximos[,-1]

valor_ns <- c(0,0.001, 0.003 ,0.005 ,0.01 ,0.017 ,0.03 ,0.05 ,0.1 ,0.17 ,0.31 ,0.56 ,1,1.77 ,3.16 ,5.62 ,10 ,17 ,31,56,100 ,177 ,316,562,1000,1778)

i <- 1
v <- 1
for (r in 1:nrow(matriz_maximos))  {
	matriz_maximos[r,2] <- valor_ns[v]
	i <- i +1
	matriz_maximos[r,1] <- gsub("Ns", "", matriz_maximos[r,1])
	if (i==50)  {
	  v <- v + 1
	  i <- 0
	}
}

matriz_maximos[1300,2] <- v
matriz_maximos[,1] <- as.numeric(matriz_maximos[,1])

sqrt(mean((matriz_maximos$parametro_real - matriz_maximos$parametro_estimado )^2))

### RMSE plot ###
library(ggplot)

p <- ggplot(data=RMSE, aes(x=theta, y=RMSE, color=theta, fill=theta)) +
+     geom_bar(stat="identity")







