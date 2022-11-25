import tskit
import numpy as np
import pandas as pd

chr22 = tskit.load("yri_200_tskit.trees")

vector_coordenadas = np.loadtxt('lista_missense.txt', delimiter=' ')

with open("missense_distancias_arbolitos.txt", 'w') as f:
    with open("missense_newick_arboles.txt", 'w') as w: #Genero y escribo archivo de texto
        for arbol in chr22.trees(): #Entro a la topologia
            for sitio in arbol.sites(): #La toplología tiene varios Site-IDs los cuales mapean a una coordenada genomica
                tabla_mutacion = getattr(sitio, 'mutations') #Cada sitio tiene una tabla de mutacion. Guardo esa info 
                i = 0 #Que se resetee para cada sitio
                if len(tabla_mutacion) == 1: #Pero hay sitios sin mutaciones. Con el if muestro que tengan mutacion
                    sub_arbol = getattr(tabla_mutacion[0], 'node') #Si si hay mutacion, entonces dime el nodo donde ocurre
                    hijo = arbol.children(sub_arbol) #Dime informacion sobre los hijos, dos opciones = 0 o 2
                    coordenada = getattr(sitio, 'position') #Dime tambien la coordenada del sitio
                    x = vector_coordenadas[i] #Dame la primer coordenada de mi vector de vep
                    while coordenada != x: #Valor coordenada se queda estatico y valor x va cambiando.
                        i = i+1
                        x = vector_coordenadas[i]
                        if i == 248:
                            break
                    if coordenada == x:
                        if (arbol.is_leaf(sub_arbol) == 1): #Si es True que sea leaf entonces = 1 y pasa. 
                            singleton = arbol.branch_length(sub_arbol) #Guardo su distancia en una varaible
                            f.write(str(singleton)) #Lo escribo convirtiendo de float a string
                            f.write(' 0')
                            f.write('\n') #Inserto un salto de linea
                        elif ((len(hijo) == 2) and (arbol.is_leaf(hijo[0]) == 1)): #Doble condicion para identificar dobleton. nodo debe tener hijos y estos ser hojas.
                            distancia_hijo = (arbol.branch_length(hijo[1]))
                            distancia_padre = (arbol.branch_length(sub_arbol))
                            #print(distancia_padre, distancia_hijo)
                            f.write(str(distancia_hijo))
                            f.write(' ') #Espacio para que archivo se lea como space delimited
                            f.write(str(distancia_padre))
                            f.write('\n')        
                        elif (arbol.is_leaf(sub_arbol) == 0): #Solo proceder si el arbol NO es hoja. --> Eliminar arboles singletones.
                            arbolito = arbol.as_newick(root=sub_arbol) #Dame el arbol en el que estamos (topologia) pero a partir del nodo mutado
                            w.write(arbolito) #Escribelo en el archivo
                            w.write('\n') #Ponle un salto de linea


