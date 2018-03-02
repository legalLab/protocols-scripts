#!/bin/bash


ARRAY=(Aotus_nancymaae Callithrix_jacchus  Carlito_syrichta Cebus_capucinus Cercocebus_atys Chlorocebus_sabaeus Colobus_angolensis Daubentonia_madagascariensis Eulemur_flavifrons Eulemur_macaco Gorilla_gorilla Homo_sapiens Macaca_fascicularis Macaca_mulatta Macaca_nemestrina Mandrillus_leucophaeus Microcebus_murinus Nasalis_larvatus Nomascus_leucogenys Otolemur_garnettii Pan_paniscus Pan_troglodytes Papio_anubis Piliocolobus_tephrosceles Pongo_abelii Propithecus_coquereli Rhinopithecus_bieti Rhinopithecus_roxellana Saimiri_boliviensis)

#echo ${ARRAY[*]}

for i in ${ARRAY[*]}; do
	echo "********"
	echo "***creating blast database***"
	echo "*******$i********"
	cd ../genomes/$i
	makeblastdb -in $i.fna -out ./blast/$i -dbtype nucl -title $i -parse_seqids
	cd ..
done

