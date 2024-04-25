#Para ingresar mis arboles de stdpopsim a R. 
#Este es el equivalente a tener la matriz de intervalos...? 

datos <- readLines("prueba_distancias_2.csv")

# Step 2: Split each row into a list of values
lista_intervalos <- lapply(datos, function(row) strsplit(row, ",")[[1]])
#If length the la lista es 3, entonces borrar el primer numero --> caso para los dobletones. 

num.obs <- sapply(lista_intervalos, length)
seq.max <- seq_len(max(num.obs))
matriz_intervalos <- t(sapply(lista_intervalos, "[", i = seq.max))

#Tengo una matriz donde las columnas son el numero de linaje. Excepto el ultimo valor, que es cuando voy de 1 a 0 linajes. 
#Osease, si la fila tiene 2 valores, es singleton, si tiene 3 es dobleton y así sucesivamente. 
#Este es el equivalente a matriz de tiempo acumulado
#Requiero ahora convertir a matriz de linajes, debería ser casi igual salvo que tengo tiempos "duplicados"

#1) Relleno todos los NAs con 0s 
#2) 

#Transformar las matrices a numerico, estan como string.

matriz_intervalos <- apply(matriz_intervalos, 2, as.numeric)

#Borrar renglones que inicien en 0. 
rows_to_keep <- matriz_intervalos[,1] != 0
matriz_intervalos <- matriz_intervalos[rows_to_keep, ]

matriz_intervalos[is.na(matriz_intervalos)] <- 0
tiempo_acumulado <- matriz_intervalos

max_linajes <- ncol(matriz_intervalos) -1 #Aqui otro cambio ncol(matriz_intervalos)
num_linajes <- 1

#Lo que hace es imprimir el numero de columna donde ocurre el 0 para tener reflejado el numero de linajes.
for (i in 1:nrow(tiempo_acumulado))  {
  num_linajes[i] <- match(0, tiempo_acumulado[i,])
  if (0 %in% num_linajes[i])  {
    num_linajes[i] <- max_linajes
  }
}

#Potencial bug en este codigo pero veamos

for (i in 1:length(num_linajes)) {
  if (num_linajes[i] == 2) {
    num_linajes[i] <- 1
  } else if (num_linajes[i] != 2) {
    num_linajes[i] <- num_linajes[i] - 2
  }
}

#Sale error... pero output parece estar bien. Tomemoslo en mente. 

 
num_linajes[is.na(num_linajes)] <- max_linajes

#rangos_tiempo <- c(0.000002, 0.000020, 0.000200, 0.002000, 0.020000, 0.200000, 2.000000, 16)
rangos_tiempo <- c(0.000002, 0.000500, 0.001000, 0.0020000, 0.020000, 0.200000, 2.000000, 16)
rangos_tiempo <- rangos_tiempo * 60000

matriz_linajes <- matrix(nrow=nrow(tiempo_acumulado), ncol=length(rangos_tiempo)) #24 columnas porque es la division de tiempo usada.

#matriz_linajes <- matriz_linajes^0 #Requerimos que todo esté en 1s para empezar el proceso.
#Creo que es mejor inicialziar en -1 
matriz_linajes <- matriz_linajes * -1

#Esto transforma tiempo acumulado a matriz linajes 
#Me esta identificando muy bien los tiempos en los que todo coalesce, pero yo buscaba contar las coalescencias? 
#Por lo menos no para singletones, para los otros parece hacerlo. 
#Creo funciona bien, solo es que los 0s deben ser 1 y los negativo deben ser 0 y atras de los negativos todos son 0 

for (r in 1:nrow(tiempo_acumulado))  {
  num_linajes[r]
  i <- 1
  for (c in 1:ncol(tiempo_acumulado))  {
    while (tiempo_acumulado[r,c] > rangos_tiempo[i])  {
      matriz_linajes[r, i] <- num_linajes[r]
      i <- i +1
    }
    num_linajes[r] <- num_linajes[r] - 1
    matriz_linajes[r,i] <- num_linajes[r]
  }
}

#Arreglo de la amtriz linajes,. 
#Marcar el numero negativo como coalescencia total aka 0 
matriz_linajes[matriz_linajes==0] <- 1
matriz_linajes[matriz_linajes<0] <- 0

matriz_conteo <- matrix(nrow=ncol(matriz_intervalos), ncol=length(rangos_tiempo)) # Para esto creamos matriz de conteo, con rows = num linajes y ncol = rangos de tiempo.
matriz_conteo <- matriz_conteo^0*0 #Queremos la matriz en 0s.

#Esto rellena la matriz conteo. 
#Hay bugs. Encontrar los bugs. 

for (c in 1:ncol(matriz_linajes))  {
  freq_col <- as.data.frame(table(matriz_linajes[,c]))
  for (x in nrow(matriz_conteo):1)  { #de 200 a 1 
    #Si el valor X se encuentra en la columna de valores, entonces pegar su freq en matriz_conteo[x,c]
    if (x %in% freq_col$Var1)  { # True si se encuentra.
      indice <- match(x, freq_col$Var1)
      fila <- (max_linajes + 1) -x
      matriz_conteo[fila, c] <- freq_col$Freq[indice] #Poner el 215 en una variable.
    }
  }
}

matriz_conteo[199,7:8] <- nrow(matriz_intervalos)
write.csv(matriz_conteo, paste("matriz_conteo_stdpopsim.csv", sep=""))


matriz_conteo_rangos <- matrix(nrow= 20, ncol=length(rangos_tiempo))


for (c in 1:ncol(matriz_conteo))  {
  suma <- 0
  i <- 1
  contador = 0
  for (r in 1:nrow(matriz_conteo))  {
    contador = contador + 1
    suma <- suma + matriz_conteo[r,c]
    if (contador==10)  { #Cambiar aqui
      matriz_conteo_rangos[i,c] <- suma
      suma <- 0
      i <- i+1
      contador = 0
    }
    else  {
      matriz_conteo_rangos[i,c] <- suma
    }
  }
}

write.csv(matriz_conteo_rangos, paste("matriz_conteo_rangos_yriDFE.csv", sep=""))


#Agregar fila invariables

vector_1 <- vector()
max_sitios <- 100000 #Checar cuantos son los arboles maximos ... Depende del tamaño de la muestra? Creo que no.

for (i in 1:ncol(matriz_conteo))  {
  vector_1[i] <- max_sitios - matriz_conteo_rangos[20,8] #numero maximo de arboles siempre esta en esa coordenada
}

matriz_conteo_invariables <- rbind(matriz_conteo_rangos, vector_1) #Hasta aqui tengo matriz con num de sitios invariables.
rownames(matriz_conteo_invariables) <- c(1:21)

write.csv(matriz_conteo_invariables, paste("matriz_conteo_yriDFE_invariables.csv", sep=""))




