
library(ape)

Valores_2Ns <- seq(2, 60, 2)

for (a in 1:length(Valores_2Ns))  {
  lista_arboles <- read.tree(paste('merged_', Valores_2Ns[a], 'Ns.txt', sep="")) 
  lista_intervalos <- list()


  for (i in 1:length(lista_arboles))  {
    ci <- coalescent.intervals(lista_arboles[[i]])
    lista_intervalos[i] <- c(ci[2])
  }
  
  num.obs <- sapply(lista_intervalos, length)
  seq.max <- seq_len(max(num.obs))
  matriz_intervalos <- t(sapply(lista_intervalos, "[", i = seq.max))

  tiempo_acumulado <- matrix(nrow=nrow(matriz_intervalos), ncol=99)

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

  max_linajes <- 98 #Aqui otro cambio
  num_linajes <- 1

  for (i in 1:nrow(tiempo_acumulado))  {
    num_linajes[i] <- match(NA, tiempo_acumulado[i,])
    if (NA %in% num_linajes[i])  {
      num_linajes[i] <- max_linajes
    }
  }

  rangos_tiempo <- c(0.000002, 0.000016, 0.000128, 16) #x8 #t4 
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

  matriz_conteo <- matrix(nrow=99, ncol=length(rangos_tiempo)) # Para esto creamos matriz de conteo, con rows = num linajes y ncol = rangos de tiempo.
  matriz_conteo <- matriz_conteo^0*0 #Queremos la matriz en 0s.


  for (c in 1:ncol(matriz_linajes))  {
    for (r in 1:nrow(matriz_linajes))  {
      if (matriz_linajes[r, c]==99)
        matriz_conteo[1,c] <- matriz_conteo[1,c] +1
      else if (matriz_linajes[r, c]==98)
        matriz_conteo[2,c] <- matriz_conteo[2,c] +1
      else if (matriz_linajes[r, c]==97)
        matriz_conteo[3,c] <- matriz_conteo[3,c] +1
      else if (matriz_linajes[r, c]==96)
        matriz_conteo[4,c] <- matriz_conteo[4,c] +1
      else if (matriz_linajes[r, c]==95)
        matriz_conteo[5,c] <- matriz_conteo[5,c] +1
      else if (matriz_linajes[r, c]==94)
        matriz_conteo[6,c] <- matriz_conteo[6,c] +1
      else if (matriz_linajes[r, c]==93)
        matriz_conteo[7,c] <- matriz_conteo[7,c] +1
      else if (matriz_linajes[r, c]==92)
        matriz_conteo[8,c] <- matriz_conteo[8,c] +1
      else if (matriz_linajes[r, c]==91)
        matriz_conteo[9,c] <- matriz_conteo[9,c] +1
      else if (matriz_linajes[r, c]==90)
        matriz_conteo[10,c] <- matriz_conteo[10,c] +1
      else if (matriz_linajes[r, c]==89)
        matriz_conteo[11,c] <- matriz_conteo[11,c] +1
      else if (matriz_linajes[r, c]==88)
        matriz_conteo[12,c] <- matriz_conteo[12,c] +1
      else if (matriz_linajes[r, c]==87)
        matriz_conteo[13,c] <- matriz_conteo[13,c] +1
      else if (matriz_linajes[r, c]==86)
        matriz_conteo[14,c] <- matriz_conteo[14,c] +1
      else if (matriz_linajes[r, c]==85)
        matriz_conteo[15,c] <- matriz_conteo[15,c] +1
      else if (matriz_linajes[r, c]==84)
        matriz_conteo[16,c] <- matriz_conteo[16,c] +1
      else if (matriz_linajes[r, c]==83)
        matriz_conteo[17,c] <- matriz_conteo[17,c] +1
      else if (matriz_linajes[r, c]==82)
        matriz_conteo[18,c] <- matriz_conteo[18,c] +1
      else if (matriz_linajes[r, c]==81)
        matriz_conteo[19,c] <- matriz_conteo[19,c] +1
      else if (matriz_linajes[r, c]==80)
        matriz_conteo[20,c] <- matriz_conteo[20,c] +1
      else if (matriz_linajes[r, c]==79)
        matriz_conteo[21,c] <- matriz_conteo[21,c] +1
      else if (matriz_linajes[r, c]==78)
        matriz_conteo[22,c] <- matriz_conteo[22,c] +1
      else if (matriz_linajes[r, c]==77)
        matriz_conteo[23,c] <- matriz_conteo[23,c] +1
      else if (matriz_linajes[r, c]==76)
        matriz_conteo[24,c] <- matriz_conteo[24,c] +1
      else if (matriz_linajes[r, c]==75)
        matriz_conteo[25,c] <- matriz_conteo[25,c] +1
      else if (matriz_linajes[r, c]==74)
        matriz_conteo[26,c] <- matriz_conteo[26,c] +1
      else if (matriz_linajes[r, c]==73)
        matriz_conteo[27,c] <- matriz_conteo[27,c] +1
      else if (matriz_linajes[r, c]==72)
        matriz_conteo[28,c] <- matriz_conteo[28,c] +1
      else if (matriz_linajes[r, c]==71)
        matriz_conteo[29,c] <- matriz_conteo[29,c] +1
      else if (matriz_linajes[r, c]==70)
        matriz_conteo[30,c] <- matriz_conteo[30,c] +1
      else if (matriz_linajes[r, c]==69)
        matriz_conteo[31,c] <- matriz_conteo[31,c] +1
      else if (matriz_linajes[r, c]==68)
        matriz_conteo[32,c] <- matriz_conteo[32,c] +1
      else if (matriz_linajes[r, c]==67)
        matriz_conteo[33,c] <- matriz_conteo[33,c] +1
      else if (matriz_linajes[r, c]==66)
        matriz_conteo[34,c] <- matriz_conteo[34,c] +1
      else if (matriz_linajes[r, c]==65)
        matriz_conteo[35,c] <- matriz_conteo[35,c] +1
      else if (matriz_linajes[r, c]==64)
        matriz_conteo[36,c] <- matriz_conteo[36,c] +1
      else if (matriz_linajes[r, c]==63)
        matriz_conteo[37,c] <- matriz_conteo[37,c] +1
      else if (matriz_linajes[r, c]==62)
        matriz_conteo[38,c] <- matriz_conteo[38,c] +1
      else if (matriz_linajes[r, c]==61)
        matriz_conteo[39,c] <- matriz_conteo[39,c] +1
      else if (matriz_linajes[r, c]==60)
        matriz_conteo[40,c] <- matriz_conteo[40,c] +1
      else if (matriz_linajes[r, c]==59)
        matriz_conteo[41,c] <- matriz_conteo[41,c] +1
      else if (matriz_linajes[r, c]==58)
        matriz_conteo[42,c] <- matriz_conteo[42,c] +1
      else if (matriz_linajes[r, c]==57)
        matriz_conteo[43,c] <- matriz_conteo[43,c] +1
      else if (matriz_linajes[r, c]==56)
        matriz_conteo[44,c] <- matriz_conteo[44,c] +1
      else if (matriz_linajes[r, c]==55)
        matriz_conteo[45,c] <- matriz_conteo[45,c] +1
      else if (matriz_linajes[r, c]==54)
        matriz_conteo[46,c] <- matriz_conteo[46,c] +1
      else if (matriz_linajes[r, c]==53)
        matriz_conteo[47,c] <- matriz_conteo[47,c] +1
      else if (matriz_linajes[r, c]==52)
        matriz_conteo[48,c] <- matriz_conteo[48,c] +1
      else if (matriz_linajes[r, c]==51)
        matriz_conteo[49,c] <- matriz_conteo[49,c] +1
      else if (matriz_linajes[r, c]==50)
        matriz_conteo[50,c] <- matriz_conteo[50,c] +1
      else if (matriz_linajes[r, c]==49)
        matriz_conteo[51,c] <- matriz_conteo[51,c] +1
      else if (matriz_linajes[r, c]==48)
        matriz_conteo[52,c] <- matriz_conteo[52,c] +1
      else if (matriz_linajes[r, c]==47)
        matriz_conteo[53,c] <- matriz_conteo[53,c] +1
      else if (matriz_linajes[r, c]==46)
        matriz_conteo[54,c] <- matriz_conteo[54,c] +1
      else if (matriz_linajes[r, c]==45)
        matriz_conteo[55,c] <- matriz_conteo[55,c] +1
      else if (matriz_linajes[r, c]==44)
        matriz_conteo[56,c] <- matriz_conteo[56,c] +1
      else if (matriz_linajes[r, c]==43)
        matriz_conteo[57,c] <- matriz_conteo[57,c] +1
      else if (matriz_linajes[r, c]==42)
        matriz_conteo[58,c] <- matriz_conteo[58,c] +1
      else if (matriz_linajes[r, c]==41)
        matriz_conteo[59,c] <- matriz_conteo[59,c] +1
      else if (matriz_linajes[r, c]==40)
        matriz_conteo[60,c] <- matriz_conteo[60,c] +1
      else if (matriz_linajes[r, c]==39)
        matriz_conteo[61,c] <- matriz_conteo[61,c] +1
      else if (matriz_linajes[r, c]==38)
        matriz_conteo[62,c] <- matriz_conteo[62,c] +1
      else if (matriz_linajes[r, c]==37)
        matriz_conteo[63,c] <- matriz_conteo[63,c] +1
      else if (matriz_linajes[r, c]==36)
        matriz_conteo[64,c] <- matriz_conteo[64,c] +1
      else if (matriz_linajes[r, c]==35)
        matriz_conteo[65,c] <- matriz_conteo[65,c] +1
      else if (matriz_linajes[r, c]==34)
        matriz_conteo[66,c] <- matriz_conteo[66,c] +1
      else if (matriz_linajes[r, c]==33)
        matriz_conteo[67,c] <- matriz_conteo[67,c] +1
      else if (matriz_linajes[r, c]==32)
        matriz_conteo[68,c] <- matriz_conteo[68,c] +1
      else if (matriz_linajes[r, c]==31)
        matriz_conteo[69,c] <- matriz_conteo[69,c] +1
      else if (matriz_linajes[r, c]==30)
        matriz_conteo[70,c] <- matriz_conteo[70,c] +1
      else if (matriz_linajes[r, c]==29)
        matriz_conteo[71,c] <- matriz_conteo[71,c] +1
      else if (matriz_linajes[r, c]==28)
        matriz_conteo[72,c] <- matriz_conteo[72,c] +1
      else if (matriz_linajes[r, c]==27)
        matriz_conteo[73,c] <- matriz_conteo[73,c] +1
      else if (matriz_linajes[r, c]==26)
        matriz_conteo[74,c] <- matriz_conteo[74,c] +1
      else if (matriz_linajes[r, c]==25)
        matriz_conteo[75,c] <- matriz_conteo[75,c] +1
      else if (matriz_linajes[r, c]==24)
        matriz_conteo[76,c] <- matriz_conteo[76,c] +1
      else if (matriz_linajes[r, c]==23)
        matriz_conteo[77,c] <- matriz_conteo[77,c] +1
      else if (matriz_linajes[r, c]==22)
        matriz_conteo[78,c] <- matriz_conteo[78,c] +1
      else if (matriz_linajes[r, c]==21)
        matriz_conteo[79,c] <- matriz_conteo[79,c] +1
      else if (matriz_linajes[r, c]==20)
        matriz_conteo[80,c] <- matriz_conteo[80,c] +1
      else if (matriz_linajes[r, c]==19)
        matriz_conteo[81,c] <- matriz_conteo[81,c] +1
      else if (matriz_linajes[r, c]==18)
        matriz_conteo[82,c] <- matriz_conteo[82,c] +1
      else if (matriz_linajes[r, c]==17)
        matriz_conteo[83,c] <- matriz_conteo[83,c] +1
      else if (matriz_linajes[r, c]==16)
        matriz_conteo[84,c] <- matriz_conteo[84,c] +1
      else if (matriz_linajes[r, c]==15)
        matriz_conteo[85,c] <- matriz_conteo[85,c] +1
      else if (matriz_linajes[r, c]==14)
        matriz_conteo[86,c] <- matriz_conteo[86,c] +1
      else if (matriz_linajes[r, c]==13)
        matriz_conteo[87,c] <- matriz_conteo[87,c] +1
      else if (matriz_linajes[r, c]==12)
        matriz_conteo[88,c] <- matriz_conteo[88,c] +1
      else if (matriz_linajes[r, c]==11)
        matriz_conteo[89,c] <- matriz_conteo[89,c] +1
      else if (matriz_linajes[r, c]==10)
        matriz_conteo[90,c] <- matriz_conteo[90,c] +1
      else if (matriz_linajes[r, c]==9)
        matriz_conteo[91,c] <- matriz_conteo[91,c] +1
      else if (matriz_linajes[r, c]==8)
        matriz_conteo[92,c] <- matriz_conteo[92,c] +1
      else if (matriz_linajes[r, c]==7)
        matriz_conteo[93,c] <- matriz_conteo[93,c] +1
      else if (matriz_linajes[r, c]==6)
        matriz_conteo[94,c] <- matriz_conteo[94,c] +1
      else if (matriz_linajes[r, c]==5)
        matriz_conteo[95,c] <- matriz_conteo[95,c] +1
      else if (matriz_linajes[r, c]==4)
        matriz_conteo[96,c] <- matriz_conteo[96,c] +1
      else if (matriz_linajes[r, c]==3)
        matriz_conteo[97,c] <- matriz_conteo[97,c] +1
      else if (matriz_linajes[r, c]==2)
        matriz_conteo[98,c] <- matriz_conteo[98,c] +1
      else if (matriz_linajes[r, c]==1)
        matriz_conteo[99,c] <- matriz_conteo[99,c] +1
    }
  }  

  matriz_conteo_rangos <- matrix(nrow=20, ncol=length(rangos_tiempo)) #Cambiar aqui
      
  for (c in 1:ncol(matriz_conteo))  {
    suma <- 0 
    i <- 1
    for (r in 1:nrow(matriz_conteo))  {
      suma <- suma + matriz_conteo[r,c]
      if (r%%5==0)  { #Cambiar aqui
        matriz_conteo_rangos[i,c] <- suma
        suma <- 0
        i <- i+1
      } 
      else if (r%%9==0)  {
        matriz_conteo_rangos[i,c] <- suma
      }
    }       
  }

  Prob_rangos_ <- matriz_conteo_rangos / nrow(matriz_linajes)

  write.csv(matriz_conteo_rangos, paste("matriz_DFE1_conteo.csv", sep=""))
  write.csv(Prob_rangos_, paste("matriz_probabilidad_DFE1.csv", sep=""))
 }
