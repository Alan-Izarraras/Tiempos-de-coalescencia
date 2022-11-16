###Codigo para preparar extraer las distancias de los arboles singletones en formato Newick que provienen de PReFerSim.

#cut -d "," -f 1 singletones_prueba.txt | gcut -b 1-3 --complement | sed 's/$/ 0/' > singletones_prueba_2.txt

#Eso ya funciona, btw el signo de dinero en el comando sed refiere al final de linea.  
#Algo pasa con la version de sed en el cluster no es la misma que mi mac. 
cut -d "," -f 1 trees_2.8Ns_singleton_replica$i.txt | cut -b 1-3 --complement | sed 's/$ /0/' 
#Un for para valor Ns
#Otro for para numero de replica.


for i in {1..10}
do
  cut -d "," -f 1 trees_2.8Ns_singleton_replica$i.txt | cut -b 1-3 --complement | sed 's/$ /0/' > trees_2.8Ns_singleton_replica$i_ready.txt  
  echo $i
done


cut -d "," -f 1 trees_2.8Ns_singleton_replica1.txt | cut -b 1-3 --complement | sed 's/$/ 0/' > trees_2.8Ns_singleton_replica1_ready.txt
cut -d "," -f 1 trees_2.8Ns_singleton_replica2.txt | cut -b 1-3 --complement | sed 's/$/ 0/' > trees_2.8Ns_singleton_replica2_ready.txt
cut -d "," -f 1 trees_2.8Ns_singleton_replica3.txt | cut -b 1-3 --complement | sed 's/$/ 0/' > trees_2.8Ns_singleton_replica3_ready.txt
cut -d "," -f 1 trees_2.8Ns_singleton_replica4.txt | cut -b 1-3 --complement | sed 's/$/ 0/' > trees_2.8Ns_singleton_replica4_ready.txt
cut -d "," -f 1 trees_2.8Ns_singleton_replica5.txt | cut -b 1-3 --complement | sed 's/$/ 0/' > trees_2.8Ns_singleton_replica5_ready.txt
cut -d "," -f 1 trees_2.8Ns_singleton_replica6.txt | cut -b 1-3 --complement | sed 's/$/ 0/' > trees_2.8Ns_singleton_replica6_ready.txt
cut -d "," -f 1 trees_2.8Ns_singleton_replica7.txt | cut -b 1-3 --complement | sed 's/$/ 0/' > trees_2.8Ns_singleton_replica7_ready.txt
cut -d "," -f 1 trees_2.8Ns_singleton_replica8.txt | cut -b 1-3 --complement | sed 's/$/ 0/' > trees_2.8Ns_singleton_replica8_ready.txt
cut -d "," -f 1 trees_2.8Ns_singleton_replica9.txt | cut -b 1-3 --complement | sed 's/$/ 0/' > trees_2.8Ns_singleton_replica9_ready.txt
cut -d "," -f 1 trees_2.8Ns_singleton_replica10.txt | cut -b 1-3 --complement | sed 's/$/ 0/' > trees_2.8Ns_singleton_replica10_ready.txt  
