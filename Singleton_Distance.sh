###Codigo para preparar extraer las distancias de los arboles singletones en formato Newick que provienen de PReFerSim.

cut -d "," -f 1 singletones_prueba.txt | gcut -b 1-3 --complement | sed 's/$/ 0/' > singletones_prueba_2.txt

#Eso ya funciona, btw el signo de dinero en el comando sed refiere al final de linea.  
