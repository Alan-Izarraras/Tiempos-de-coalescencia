chr1 = tskit.load("Complete_chr1_stdpopsim_replica1.trees")

lista_nodos = []
lista_tiempos = []
a = 0
b = 0
c = 0
temp_ = []

with open("chr1_stdpopsim_distancias.txt", 'w') as w: 
    for arbol in chr1.trees(): #Cada arbol de la secuencia de arboles
        for mutations in arbol.mutations(): #Cada mutación del arbol
            mut = getattr(mutations, 'metadata') #Me interesa extraer los metadatos donde encuentro el tipo de mutacion
            mut_type = mut["mutation_list"][0]["mutation_type"] #Hay 2 tipos de mutación, mutaciones dos son deletereas y de interes
            if (mut_type == 2): #Si es mutación 2 entonces
                sub_arbol = getattr(mutations, "node") #Nodo donde ocurre la mutación
                tiempo = getattr(mutations, "time") #Tiempo exacto de la mutación
                hijo = arbol.children(sub_arbol) #Para saber si es singleton o no
                arbolito = arbol.as_newick(root=sub_arbol) #El arbol en formato Newick
                distancia = arbol.branch_length(sub_arbol) #La distancia total del subarbol
                num_muestras = arbol.num_samples(sub_arbol) #total de linajes, podría ser util
                #print("separacion")
                lista_sub_arbol = [] #Contiene las distancias de cada nodo de un subarbol
                #Por que necesito que se reescriba en cada loop de mutations?? 
                lista_subarbol_ordenada = []
                lista_sub_arbol.append(tiempo) #Aqui estoy poniendo el tiempo de toodas las mutaciones m2
             
                for muestra in arbol.samples(sub_arbol): #Entro a los nodos del subarbol. Cada loop es la info de un nodo individual                     
                    if (num_muestras != 1): #Si no es singleton      
                        distancia_del_nodo = arbol.branch_length(muestra) #obtengo la distancia a cada nodo de todo el subarbol
                        lista_sub_arbol.append(distancia_del_nodo) 
                        lista_subarbol_ordenada = sorted(lista_sub_arbol) #Se tiene que sortear para ir en orden de tiempo
                        c = c + 1
                    else: #Si es singleton
                        lista_subarbol_ordenada = sorted(lista_sub_arbol)
                        #No estoy imprimiendo los singletones..? 
                if (len(lista_subarbol_ordenada) <= 200): #Esta condición excluye mutaciones que se fijaron.
                    temp_.append(lista_subarbol_ordenada)
                    b = b + 1
                else:
                    a = a + 1
                    
    ##Sortea por frecuencia    
    temp = sorted(temp_, key=len)
    
    for i in temp:
        w.write(str(i))
        w.write("\n")

### Quita los square brackets        
with open("chr1_stdpopsim_distancias.txt", "r") as input_file:
    contenido = input_file.read()
    brackets = contenido.replace("[", "").replace("]", "")
with open("chr1_stdpopsim_distancias.txt", "w") as output_file:
    output_file.write(brackets)
