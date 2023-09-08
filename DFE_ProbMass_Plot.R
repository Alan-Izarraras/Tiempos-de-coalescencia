### Preparar matriz para hacer barplots para distintos valores de s. DFE real es columna 1 y DFE estimada es columna dos  
#Okey pero de donde sale vector_maximos...??? 
#Kim et al. (a=0.186, b=706) --> posicion 1747
#Boyko et al. (a=0.184, b=6397) --> posicion 1750

vector_maximos <- read.csv("vector_maximos_boykoDFE.csv")
vector_maximos <-  vector_maximos[,-1]

Table <- read.csv("Table_of_probabilities_boykoDFE.csv", sep=",")
Table <- Table[,-1]

#Genero la masa de probabilidad para los intervalos (columnas) de las 50 corridas independientes (filas)
intervalo <- matrix(nrow=length(vector_maximos), ncol=5)
i <- 1
vector_maximos <- vector_maximos -1 #Tengo que hcer esto porque el vector maximos esta basando en la tabla de likelihoods que tiene fila 1, Pero estoy oprenaod sobre la tabla de probabilidades que tiene fila 0. 
for (n in vector_maximos)  {
  intervalo[i,1] <- sum(Table[n,3], Table[n,4], Table[n,5], Table[n,6], Table[n,7], Table[n,8], Table[n,9], Table[n,10], Table[n,11]) #de 0 a 10^5 Ns
  intervalo[i,2] <- sum(Table[n,12], Table[n, 13], Table[n,14], Table[n,15]) #de 10^5 a 10^4
  intervalo[i,3] <- sum(Table[n,16], Table[n,17], Table[n,18], Table[n,19]) # de 10^4 a 10^3
  intervalo[i,4] <- sum(Table[n,20], Table[n,21], Table[n,22], Table[n,23]) # de 10^3 a 10^2
  intervalo[i,5] <- sum(Table[n,24], Table[n,25], Table[n,26], Table[n,27], Table[n,28]) # > a 10^2 
  i <- i + 1
}

#Ahora, busco partirla en dos: Una matriz con el valor promedio y otra con los minimos y maximos. 
intervalo_mean <- matrix(nrow=2, ncol=5)
row.names(intervalo_mean) <- c("estimado", "boyko")
for (i in 1:5)  {
  intervalo_mean[1,i] <- mean(intervalo[,i])
}

#Fila 1 es mi estimación 
#Doble check que fila 1747 corresponda a distribucion real de boyko
#n Cambia para kim et al. 

n <- 1750 #Indice de gammas verdaderas. #Llenar segunda fila correspondiente bines de la dist real.
intervalo_mean[2,1] <- sum(Table[n,3], Table[n,4], Table[n,5], Table[n,6], Table[n,7], Table[n,8], Table[n,9], Table[n,10], Table[n,11])
intervalo_mean[2,2] <- sum(Table[n,12], Table[n, 13], Table[n,14], Table[n,15])
intervalo_mean[2,3] <- sum(Table[n,16], Table[n,17], Table[n,18], Table[n,19])
intervalo_mean[2,4] <- sum(Table[n,20], Table[n,21], Table[n,22], Table[n,23])
intervalo_mean[2,5] <- sum(Table[n,24], Table[n,25], Table[n,26], Table[n,27], Table[n,28])

#Vector de minimos y maximos para colocar el rango de error. 
intervalo_min_max <- matrix(nrow=2, ncol=5)
row.names(intervalo_min_max) <- c("min", "max")
for (i in 1:5)  {
  intervalo_min_max[1,i] <- min(intervalo[,i])
  intervalo_min_max[2,i] <- max(intervalo[,i])
}  

##Bar plots para los 5 intervalos
#A los graficos de R les gusta que estne ordenados en dos columnas: Una para el valor y otra para la etiqueta del valor. 
#Quiero plottear intervalo_mean 

categories <- c("1", "2", "3", "4", "5")
subgroups <- c("estimado", "boyko", "estimado", "boyko", "estimado", "boyko", "estimado", "boyko", "estimado", "boyko")
#barplot(intervalo_mean, main="Probability mass", beside = TRUE, names.arg = categories, ylim=c(0.0, 0.8))

min_values <- intervalo_min_max[1,]
max_values <- intervalo_min_max[2,]

#Con ggplot

data_estimated <- data.frame(
  Category = categories,
  Subgroup = "estimated",
  Mean = intervalo_mean[1, ],
  Min = min_values,
  Max = max_values
)

data_boyko <- data.frame(
  Category = categories,
  Subgroup = "boyko",
  Mean = intervalo_mean[2, ],
  Min = NA,
  Max = NA
)

data <- rbind(data_estimated, data_boyko)

data$Subgroup <- factor(data$Subgroup, levels = c("estimated", "boyko"))

ggplot(data, aes(x = Category, y = Mean, fill = factor(Category), alpha = Subgroup)) +
  geom_bar(stat = "identity", position = "dodge") +
  geom_errorbar(aes(ymin = Min, ymax = Max), width = 0.2, position = position_dodge(width = 0.9)) +
  labs(title = "Grouped Bar Plot Example", x = "Categories", y = "Values") +
  ylim(0, 0.6) +
  scale_fill_discrete(name = "Categories") +
  scale_alpha_manual(name = "Subgroups", values = c("estimated" = 1, "boyko" = 0.6))
  scale_x_discrete(labels = c("1" = "0 < s < 10^-5", "2" = "10^-5 < s < 10^-4", "3" = "10^-4 < s < 10^-3", "4" = "10^-3 < s < 10^-2", "5" = "s > 10^-2"))
  theme_minimal()

#Quedó pendiente subir la version de este codigo a Github.
#Doble chequear las coordenadas de la distr real porque fueron las mismos
#Aplicarlo una y otra vez para los datos de demografias. 



