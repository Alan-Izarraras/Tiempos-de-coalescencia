library(ape)


valor_ns <- c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26)
replica <- seq(1,50)
corridas <- c(100)

for (l in corridas)  {
  for (a in valor_ns)  {
    for (n in replica)  {

      lista_arboles <- read.tree(paste('trees_', a, 'Ns_', l, 'runs_replica', n, '.txt', sep=''))
      lista_intervalos <- list()

      if (length(lista_arboles) == 0)  {

        matriz_intervalos <- t(sapply(lista_intervalos, "[", i = seq.max))
        singletones <- read.table(paste(a, 'Ns_singleton_', l, 'runs_replica', n,  '_ready.txt', sep=''))
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
        matriz_conteo_rangos_1 <- matrix(nrow= 20, ncol=length(rangos_tiempo))
        matriz_conteo_rangos_2 <- matrix(nrow= 33, ncol=length(rangos_tiempo))
        matriz_conteo_rangos_3 <- matrix(nrow= 10, ncol=length(rangos_tiempo))

        write.csv(matriz_conteo, paste("matriz_conteo_full_", a, "Ns_", l, "runs_replica", n, ".csv", sep=""))

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

        ## Version 1
        for (c in 1:ncol(matriz_conteo))  {
          suma <- 0
          i <- 1
          contador = 0
          for (r in 1:nrow(matriz_conteo))  {
            contador = contador + 1
            suma <- suma + matriz_conteo[r,c]
            if (contador==5)  { #Cambiar aqui
              matriz_conteo_rangos_1[i,c] <- suma
              suma <- 0
              i <- i+1
              contador = 0
            }
            else  {
              matriz_conteo_rangos_1[i,c] <- suma
            }
          }
        }

        ## Version 2

        for (c in 1:ncol(matriz_conteo))  {
          suma <- 0
          i <- 1
          contador = 0
          for (r in 1:nrow(matriz_conteo))  {
            contador = contador + 1
            suma <- suma + matriz_conteo[r,c]
            if (contador==3)  { #Cambiar aqui
              matriz_conteo_rangos_2[i,c] <- suma
              suma <- 0
              i <- i+1
              contador = 0
            }
            else  {
              matriz_conteo_rangos_2[i,c] <- suma
            }
          }
        }

        #Version 3

        for (c in 1:ncol(matriz_conteo))  {
          suma <- 0
          i <- 1
          contador = 0
          for (r in 1:nrow(matriz_conteo))  {
            contador = contador + 1
            suma <- suma + matriz_conteo[r,c]
            if (contador==10)  { #Cambiar aqui
              matriz_conteo_rangos_3[i,c] <- suma
              suma <- 0
              i <- i+1
              contador = 0
            }
            else  {
              matriz_conteo_rangos_3[i,c] <- suma
            }
          }
        }

        matriz_corta <- matriz_conteo[91:99,] #lol si funcionó.
        num_arboles <- matriz_conteo_rangos[10,8]
        matriz_conteo_rangos <- matriz_conteo_rangos[-10,]
        matriz_conteo_rangos_3 <- rbind(matriz_conteo_rangos, matriz_corta)

        rownames(matriz_conteo_rangos_3) <- seq(1, 18)

        ###Parte para... agregar el renglon de sitios invariables###

        vector_1 <- 1
        max_sitios_1 <- 100000

        for (i in 1:ncol(matriz_conteo_rangos))  {
          vector_1[i] <- max_sitios_1 - num_arboles
        }

        matriz_theta_1 <- rbind(matriz_conteo_rangos, vector_1)
        matriz_theta_2 <- rbind(matriz_conteo_rangos_1, vector_1)
        matriz_theta_3 <- rbind(matriz_conteo_rangos_2, vector_1)
        matriz_theta_4 <- rbind(matriz_conteo_rangos_3, vector_1)

        write.csv(matriz_theta_1, paste("matriz_conteo_smooth_", a, "Ns_", l, "runs_rango1_replica", n, ".csv", sep=""))
        write.csv(matriz_theta_2, paste("matriz_conteo_smooth_", a, "Ns_", l, "runs_rango2_replica", n, ".csv", sep=""))
        write.csv(matriz_theta_3, paste("matriz_conteo_smooth_", a, "Ns_", l, "runs_rango3_replica", n, ".csv", sep=""))
        write.csv(matriz_theta_4, paste("matriz_conteo_smooth_", a, "Ns_", l, "runs_rango4_replica", n, ".csv", sep=""))


      } else {

        for (i in 1:length(lista_arboles))  {
          ci <- coalescent.intervals(lista_arboles[[i]])
          lista_intervalos[i] <- c(ci[2])
        }

        num.obs <- sapply(lista_intervalos, length)
        seq.max <- seq_len(max(num.obs))
        matriz_intervalos <- t(sapply(lista_intervalos, "[", i = seq.max))


        ###Incorporar longitudes de singletones y dobletones a la matriz de intervalos. Usar para arboles Relate

        singletones <- read.table(paste(a, 'Ns_singleton_', l, 'runs_replica', n,  '_ready.txt', sep=''))
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

        matriz_conteo_rangos <- matrix(nrow= 10, ncol=length(rangos_tiempo))
        matriz_conteo_rangos_1 <- matrix(nrow= 20, ncol=length(rangos_tiempo))
        matriz_conteo_rangos_2 <- matrix(nrow= 33, ncol=length(rangos_tiempo))
        matriz_conteo_rangos_3 <- matrix(nrow= 10, ncol=length(rangos_tiempo))

        write.csv(matriz_conteo, paste("matriz_conteo_full_", a, "Ns_", l, "runs_replica", n, ".csv", sep=""))

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

        ## Version 1
        for (c in 1:ncol(matriz_conteo))  {
          suma <- 0
          i <- 1
          contador = 0
          for (r in 1:nrow(matriz_conteo))  {
            contador = contador + 1
            suma <- suma + matriz_conteo[r,c]
            if (contador==5)  { #Cambiar aqui
              matriz_conteo_rangos_1[i,c] <- suma
              suma <- 0
              i <- i+1
              contador = 0
            }
            else  {
              matriz_conteo_rangos_1[i,c] <- suma
            }
          }
        }

        ## Version 2

        for (c in 1:ncol(matriz_conteo))  {
          suma <- 0
          i <- 1
          contador = 0
          for (r in 1:nrow(matriz_conteo))  {
            contador = contador + 1
            suma <- suma + matriz_conteo[r,c]
            if (contador==3)  { #Cambiar aqui
              matriz_conteo_rangos_2[i,c] <- suma
              suma <- 0
              i <- i+1
              contador = 0
            }
            else  {
              matriz_conteo_rangos_2[i,c] <- suma
            }
          }
        }

        #Version 3

        for (c in 1:ncol(matriz_conteo))  {
          suma <- 0
          i <- 1
          contador = 0
          for (r in 1:nrow(matriz_conteo))  {
            contador = contador + 1
            suma <- suma + matriz_conteo[r,c]
            if (contador==10)  { #Cambiar aqui
              matriz_conteo_rangos_3[i,c] <- suma
              suma <- 0
              i <- i+1
              contador = 0
            }
            else  {
              matriz_conteo_rangos_3[i,c] <- suma
            }
          }
        }

        matriz_corta <- matriz_conteo[91:99,] #lol si funcionó.
        num_arboles <- matriz_conteo_rangos[10,8]
        matriz_conteo_rangos <- matriz_conteo_rangos[-10,]
        matriz_conteo_rangos_3 <- rbind(matriz_conteo_rangos, matriz_corta)

        rownames(matriz_conteo_rangos_3) <- seq(1, 18)

        ###Parte para... agregar el renglon de sitios invariables###

        vector_1 <- 1
        max_sitios_1 <- 100000


        for (i in 1:ncol(matriz_conteo_rangos))  {
          vector_1[i] <- max_sitios_1 - num_arboles
        }

        matriz_theta_1 <- rbind(matriz_conteo_rangos, vector_1)
        matriz_theta_2 <- rbind(matriz_conteo_rangos_1, vector_1)
        matriz_theta_3 <- rbind(matriz_conteo_rangos_2, vector_1)
        matriz_theta_4 <- rbind(matriz_conteo_rangos_3, vector_1)

        write.csv(matriz_theta_1, paste("matriz_conteo_smooth_", a, "Ns_", l, "runs_rango1_replica", n, ".csv", sep=""))
        write.csv(matriz_theta_2, paste("matriz_conteo_smooth_", a, "Ns_", l, "runs_rango2_replica", n, ".csv", sep=""))
        write.csv(matriz_theta_3, paste("matriz_conteo_smooth_", a, "Ns_", l, "runs_rango3_replica", n, ".csv", sep=""))
        write.csv(matriz_theta_4, paste("matriz_conteo_smooth_", a, "Ns_", l, "runs_rango4_replica", n, ".csv", sep=""))


      }

    }

  }

}
