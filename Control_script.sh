### Script de bash mara mandar a correr muchos qsub
### Ejecuta MasterScript_AllFreqs.sge con distintos valores de qsub. 


Sel=()
Sel=({0..1200000..100000})
num_corridas=100


#echo ${Freq[97]} #Pido el indice 37 del array. Indice maximo es 97
#echo ${Sel[12]} #Indice maximo es 12


for x in {0..12} #0 a 12
do
  	qsub1=$(expr ${Sel[x]} + 1)
        qsub2=$(expr ${Sel[x]} + $num_corridas)
        qsub -t $qsub1-$qsub2 MasterScript_AllFreqs.sge
        echo se esta corriendo $qsub1 a $qsub2
done

#echo se corrieron las siguientes IDs: $qsub1 a $qsub2