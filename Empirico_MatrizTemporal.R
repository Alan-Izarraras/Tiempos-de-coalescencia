library(ape)


  lista_arboles <- read.tree('Crom22_subset_limpio.txt') 
  lista_intervalos <- list()


  for (i in 1:length(lista_arboles))  {
    ci <- coalescent.intervals(lista_arboles[[i]])
    lista_intervalos[i] <- c(ci[2])
  }
  
  num.obs <- sapply(lista_intervalos, length)
  seq.max <- seq_len(max(num.obs))
  matriz_intervalos <- t(sapply(lista_intervalos, "[", i = seq.max))

  #El objetivo ahora es appendear la matriz de singletones a la matriz de intervalos. 
  #rbind() soloi funciona cuando el numero de columnas son iguales, entonces lo que se puede hacer es insertar ese numero de columnas
  #A mi matriz de singletones. 
  #En este caso el maximo son 214 columnasentonces insertaria 212
  colnames(datos) <- c(1,2)
  vacio <- matrix(nrow=nrow(datos), ncol=212)
  datos_v2 <- cbind(datos, vacio)
  colnames(matriz_intervalos) <- seq(1:214)
  fusion <- rbind(datos_v2, matriz_intervalos)
  #Fusion es la matriz de interalos incluyendo singletones y dobletones! 

  ###Recuerda cambiar nombres de variables de esta parte nueva###


  tiempo_acumulado <- matrix(nrow=nrow(matriz_intervalos), ncol=ncol(matriz_intervalos))

  for (r in 1:nrow(matriz_intervalos))  {  
    i<-1                                                    
    suma <- matriz_intervalos[r,i]  
    tiempo_acumulado[r,i] <- suma   
    for (c in 1:ncol(matriz_intervalos))  { 
      suma <- suma + matriz_intervalos[r, i+1]  
      tiempo_acumulado[r, i+1] <- suma    
      i <- i+1                                                               
      if (i==ncol(matriz_intervalos))  { #Aqui el cambio 99 por ncol(matriz_intervalos)
        break  
      }
    }
  }

  #Se tarda mucho la mac aqui, entonces lo mejor seria... subsetear y trabajar con poquitos arboles
  #Estandarizarlo de esa manera
  #Despues ya desarrollo una pipeline/scripts formales y se los doy al cluster. 

  max_linajes <- ncol(matriz_intervalos) #Aqui otro cambio ncol(matriz_intervalos)
  num_linajes <- 1

  for (i in 1:nrow(tiempo_acumulado))  {
    num_linajes[i] <- match(NA, tiempo_acumulado[i,])
    if (NA %in% num_linajes[i])  {
      num_linajes[i] <- max_linajes
    }
  }
  
  rangos_tiempo <- c(0.000002, 0.00002, 0.0002, 0.002, 0.02, 0.2, 2, 20, 200, 2000, 20000, 200000) #Relate
  #rangos_tiempo_Relate <- c(0.000002, 0.000016, 0.000128, 0.001024, 0.008192, 0.065536, 0.524288, 4.194304, 33, 268, 2147, 17179, 137438) #x8
  #rangos_tiempo <- c(0.000002, 0.000016, 0.000128, 16) #x8 #t4 
  #rangos_tiempo <- c(0.000002, 0.000020, 0.000200, 0.002000, 0.020000, 0.200000, 2.000000, 16)
  #rangos_tiempo <- c(0.000002, 0.000016, 0.000128, 0.001024, 0.008192, 0.065536,  0.524288, 4.194304, 16)
  #rangos_tiempo <- c(0.000002, 0.000004, 0.000008, 0.000016, 0.000032, 0.000064, 0.000128, 0.000256, 0.000512, 0.001024, 0.002048, 0.004096, 0.008192, 0.016392, 0.032784, 0.065568, 0.131136, 0.264472, 0.528944, 1.057888, 2.115776, 4.231552, 8.643104, 16)
  #Falta un rango de tiempo alternativo que no llegue hasta el tiempo maximo. 

  matriz_intervalos[is.na(matriz_intervalos)] <- 0
  tiempo_acumulado[is.na(tiempo_acumulado)] <- 0

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

##Hmmmm no estoy seguro del valor de nrow. Ahhh ya me acorde, si tiene que ser el numero maximo de linajes. 
  matriz_conteo <- matrix(nrow=ncol(matriz_intervalos), ncol=length(rangos_tiempo)) # Para esto creamos matriz de conteo, con rows = num linajes y ncol = rangos de tiempo.
  matriz_conteo <- matriz_conteo^0*0 #Queremos la matriz en 0s.

#Ugh. Tengo que transformar esto para 214 linajes. 
#Seria buen momento para averiguar un algoritmo más pro que no involucre 200+ lineas xd
#Listo, me ahorre muchas lineas. 

#Para cada columna de la matriz de linajes te da la frecuencia de cada valor observado en una matriz nueva
#En esta nueva matriz con la frecuencia de cada valor, preguntamos si tiene cada valor en orden descendente (mayor a menor)
#Si si lo tiene entonces le preguntamos el indice para poder sacar su frecuencia y lo imprimimos en la matriz conteo en el 
#lugar que correspindiente que se saca haciendo una resta del numero maximo de linajes con el valor a buscar 
#Ya que los linajes en matriz conteo van de mayor a menor. 

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

#Aqui la idea es disminuir el numero de rangos, e.g. agrupar de 10 en 10 linajes + el pilon.
  matriz_conteo_rangos <- matrix(nrow=(max_linajes/10 + 1), ncol=length(rangos_tiempo)) #Agrupando de 10 en 10. 
  #matriz_conteo_rangos <- matrix(nrow=20, ncol=length(rangos_tiempo)) #Cambiar aqui
     
#Ejemplo para divir en rangos de 10 + el pilon 
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

  Prob_rangos_ <- matriz_conteo_rangos / nrow(matriz_linajes)

  write.csv(matriz_conteo_rangos, paste("matriz_DFE1_conteo.csv", sep=""))
  write.csv(Prob_rangos_, paste("matriz_probabilidad_DFE1.csv", sep=""))
 }
