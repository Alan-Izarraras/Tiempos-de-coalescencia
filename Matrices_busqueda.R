#Este script corrige arboles vacios Y manda todos los valores.
#Esto ya debería estar listo para ejecutarse en el cluster pero uno nunca sabe que bugs pasen

library(ape)

lista_arboles <- read.tree(paste('trees_22Ns_10runs_replica30.txt', sep=''))
lista_intervalos <- list()

if (length(lista_arboles) == 0)  {

  matriz_intervalos <- t(sapply(lista_intervalos, "[", i = seq.max))
  singletones <- read.table(paste('22Ns_singleton_10runs_replica30_ready.txt', sep=''))
  colnames(singletones) <- c(1,2)
  matriz_intervalos <- singletones

  tiempo_acumulado <- matrix(nrow=nrow(matriz_intervalos), ncol=ncol(matriz_intervalos))

  tiempo_acumulado <- singletones

  max_linajes <- ncol(matriz_intervalos) #Aqui otro cambio ncol(matriz_intervalos)
  num_linajes <- 1

  for (i in 1:nrow(tiempo_acumulado))  {
    num_linajes[i] <- match(0, tiempo_acumulado[i,])
    if (0 %in% num_linajes[i])  {
      num_linajes[i] <- max_linajes
    }
  }

  num_linajes[is.na(num_linajes)] <- max_linajes

  rangos_tiempo <- c(0.000002, 0.000020, 0.000200, 0.002000, 0.020000, 0.200000, 2.000000, 16)

  matriz_linajes <- matrix(nrow=nrow(tiempo_acumulado), ncol=length(rangos_tiempo)) #24 columnas porque es la division de tiempo usada.

  matriz_linajes <- matriz_linajes^0 #Requerimos que todo esté en 1s para empezar el proceso.

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


  matriz_linajes[matriz_linajes<=0] <- 1

  matriz_conteo <- matrix(nrow=ncol(matriz_intervalos), ncol=length(rangos_tiempo)) # Para esto creamos matriz de conteo, con rows = num linajes y ncol = rangos de tiempo.
  matriz_conteo <- matriz_conteo^0*0 #Queremos la matriz en 0s.

  for (c in 1:ncol(matriz_linajes))  {
    freq_col <- as.data.frame(table(matriz_linajes[,c]))
    for (x in nrow(matriz_conteo):1)  {
      #Si el valor X se encuentra en la columna de valores, entonces pegar su freq en matriz_conteo[x,c]
      if (x %in% freq_col$Var1)  { # True si se encuentra.
        indice <- match(x, freq_col$Var1)
        fila <- (max_linajes + 1) -x
        matriz_conteo[fila, c] <- freq_col$Freq[indice] #Poner el 215 en una variable.
      }
    }
  }

  filas_extra <- 99 - ncol(matriz_intervalos)
  matriz_extra <- matrix(nrow = filas_extra, ncol = ncol(matriz_conteo))
  matriz_extra[is.na(matriz_extra)] <- 0

  matriz_conteo <- rbind(matriz_extra, matriz_conteo)
  matriz_conteo_rangos <- matrix(nrow= 10, ncol=length(rangos_tiempo))

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

  write.csv(matriz_conteo_rangos, paste("matriz_conteo_smooth_22Ns_10runs_replica30.csv", sep=""))

} else {

  cuenta <- read.table("trees_22Ns_10runs_replica30.txt")
  if (nrow(cuenta) == 1)  {
    ci <- coalescent.intervals(lista_arboles)
    lista_intervalos[1] <- c(ci[2])
  }
  else  {
    for (i in 1:length(lista_arboles))  {
      ci <- coalescent.intervals(lista_arboles[[i]])
      lista_intervalos[i] <- c(ci[2])
    }
  }

  num.obs <- sapply(lista_intervalos, length)
  seq.max <- seq_len(max(num.obs))
  matriz_intervalos <- t(sapply(lista_intervalos, "[", i = seq.max))


  ###Incorporar longitudes de singletones y dobletones a la matriz de intervalos. Usar para arboles Relate

  singletones <- read.table(paste('22Ns_singleton_10runs_replica30_ready.txt', sep=''))
  colnames(singletones) <- c(1,2)
  matriz_temp <- matrix(nrow=nrow(singletones), ncol=ncol(matriz_intervalos)-2)
  matriz_temp <- cbind(singletones, matriz_temp)
  colnames(matriz_temp) <- seq(1:ncol(matriz_intervalos))
  colnames(matriz_intervalos) <- seq(1:ncol(matriz_intervalos))
  matriz_intervalos <- rbind(matriz_temp, matriz_intervalos)


  tiempo_acumulado <- matrix(nrow=nrow(matriz_intervalos), ncol=ncol(matriz_intervalos))

  for (r in 1:nrow(matriz_intervalos))  {
    i<-1
    suma <- matriz_intervalos[r,i]
    tiempo_acumulado[r,i] <- suma
    for (c in 1:ncol(matriz_intervalos))  {
      suma <- suma + matriz_intervalos[r, i+1]
      tiempo_acumulado[r, i+1] <- suma
      i <- i+1
      if (tiempo_acumulado[r,1] - tiempo_acumulado[r,2] == 0)  {
        tiempo_acumulado[r,2] = 0
      }
      if (i==ncol(matriz_intervalos))  { #Aqui el cambio 99 por ncol(matriz_intervalos)
        break
      }
    }
  }


  matriz_intervalos[is.na(matriz_intervalos)] <- 0
  tiempo_acumulado[is.na(tiempo_acumulado)] <- 0

  max_linajes <- ncol(matriz_intervalos) #Aqui otro cambio ncol(matriz_intervalos)
  num_linajes <- 1

  for (i in 1:nrow(tiempo_acumulado))  {
    num_linajes[i] <- match(0, tiempo_acumulado[i,])
    if (0 %in% num_linajes[i])  {
      num_linajes[i] <- max_linajes
    }
  }

  num_linajes[is.na(num_linajes)] <- max_linajes

  rangos_tiempo <- c(0.000002, 0.000020, 0.000200, 0.002000, 0.020000, 0.200000, 2.000000, 16)

  matriz_linajes <- matrix(nrow=nrow(tiempo_acumulado), ncol=length(rangos_tiempo)) #24 columnas porque es la division de tiempo usada.

  matriz_linajes <- matriz_linajes^0 #Requerimos que todo esté en 1s para empezar el proceso.


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

  matriz_linajes[matriz_linajes<=0] <- 1

  matriz_conteo <- matrix(nrow=ncol(matriz_intervalos), ncol=length(rangos_tiempo)) # Para esto creamos matriz de conteo, con rows = num linajes y ncol = rangos de tiempo.
  matriz_conteo <- matriz_conteo^0*0 #Queremos la matriz en 0s.


  for (c in 1:ncol(matriz_linajes))  {
    freq_col <- as.data.frame(table(matriz_linajes[,c]))
    for (x in nrow(matriz_conteo):1)  {
      #Si el valor X se encuentra en la columna de valores, entonces pegar su freq en matriz_conteo[x,c]
      if (x %in% freq_col$Var1)  { # True si se encuentra.
        indice <- match(x, freq_col$Var1)
        fila <- (max_linajes + 1) -x
        matriz_conteo[fila, c] <- freq_col$Freq[indice] #Poner el 215 en una variable.
      }
    }
  }

  filas_extra <- 99 - ncol(matriz_intervalos)
  matriz_extra <- matrix(nrow = filas_extra, ncol = ncol(matriz_conteo))
  matriz_extra[is.na(matriz_extra)] <- 0

  matriz_conteo <- rbind(matriz_extra, matriz_conteo)
  matriz_conteo_rangos <- matrix(nrow= 10, ncol=length(rangos_tiempo))

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

  write.csv(matriz_conteo_rangos, paste("matriz_conteo_smooth_22Ns_10runs_replica30.csv", sep=""))

}

