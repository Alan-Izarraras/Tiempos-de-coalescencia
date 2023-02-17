#Script para generar matrices de probabilidad/referencia correspondiendo al set2 de datos.
#Aun utilizando datos empiricos estas matrices provienen de datos de simulación porque son las referencias
#Para usarlo junto con datos empiricos de Relate, requiere un paso de ajuste demografico que está comentado en este script

library(ape)
valor_ns <- c(0.01, 0.1, 0, 1, 10, 100)

for (a in valor_ns)  {

  lista_arboles <- read.tree(paste(a, 'Ns_treespace_100runs_set2.txt', sep=''))
  lista_intervalos <- list()


  for (i in 1:length(lista_arboles))  {
    ci <- coalescent.intervals(lista_arboles[[i]])
    lista_intervalos[i] <- c(ci[2])
  }

  num.obs <- sapply(lista_intervalos, length)
  seq.max <- seq_len(max(num.obs))
  matriz_intervalos <- t(sapply(lista_intervalos, "[", i = seq.max))

    #Aqui es donde tengo que hacer el ajuste de multiplicar por 4N_diploide que tambien podria hacer un expeirmento de 4N_haploide
    #Entonces tengo que multiplicar por 4(5,000)

    #matriz_intervalos <- matriz_intervalos * 57896 #Podria ser 14,474 o 28,948

    #Incorporar longitudes de singletones y dobletones a la matriz de intervalos. Usar para arboles Relate

  singletones <- read.table(paste(a, 'Ns_singleton_100runs_ready_set2.txt', sep=''))
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

    
  #rangos_tiempo <- c(0.000002, 0.00002, 0.0002, 0.002, 0.02, 0.2, 2, 20, 200, 2000, 20000, 400000) #Relate
  #rangos_tiempo_Relate <- c(0.000002, 0.000016, 0.000128, 0.001024, 0.008192, 0.065536, 0.524288, 4.194304, 33, 268, 2147, 17179, 137438) #x8
  #rangos_tiempo <- c(0.000002, 0.000016, 0.000128, 16) #x8 #t4
  rangos_tiempo <- c(0.000002, 0.000020, 0.000200, 0.002000, 0.020000, 0.200000, 2.000000, 16)
  #rangos_tiempo <- c(0.000002, 0.000016, 0.000128, 0.001024, 0.008192, 0.065536,  0.524288, 4.194304, 16)
  #rangos_tiempo <- c(0.000002, 0.000004, 0.000008, 0.000016, 0.000032, 0.000064, 0.000128, 0.000256, 0.000512, 0.001024, 0.002048, 0.004096, 0.008192, 0.016392, 0.032784, 0.065568, 0.131136, 0.264472, 0.528944, 1.057888, 2.115776, 4.231552, 8.643104, 16)
  #Falta un rango de tiempo alternativo que no llegue hasta el tiempo maximo.


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
  matriz_conteo <- matriz_conteo + 0.1  #Laplace Creo que hacerlo desde aqui es lo mas correcto. 

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


  Prob_rangos_ <- matriz_conteo_rangos / matriz_conteo_rangos[10,8]

  write.csv(Prob_rangos_, paste("matriz_probabilidad_", a, "Ns_100runs_set2_smooth.csv", sep=""))
}
