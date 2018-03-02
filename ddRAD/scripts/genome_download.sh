#!/usr/bin/env bash
# script to download all primate genomes used in this study
# summary can be found at: http://www.ncbi.nlm.nih.gov/genome/[insert GenBank No.]

# 451 Otolemur garnettii (small-eared galago)
mkdir -p ../genomes/primata/Otolemur_garnettii
wget "ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF_000181295.1_OtoGar3/GCF_000181295.1_OtoGar3_genomic.fna.gz" -O ../genomes/primata/Otolemur_garnettii/Otolemur_garnettii.fna.gz
gzip -d ../genomes/primata/Otolemur_garnettii/Otolemur_garnettii.fna.gz

# 24390 Propithecus coquereli (Coquerel's sifaka)
mkdir -p ../genomes/primata/Propithecus_coquereli
wget "ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF_000956105.1_Pcoq_1.0/GCF_000956105.1_Pcoq_1.0_genomic.fna.gz" -O ../genomes/primata/Propithecus_coquereli/Propithecus_coquereli.fna.gz
gzip -d ../genomes/primata/Propithecus_coquereli/Propithecus_coquereli.fna.gz

# 777 Microcebus murinus (gray mouse lemur)
mkdir -p ../genomes/primata/Microcebus_murinus
wget "ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF_000165445.1_Mmur_2.0/GCF_000165445.1_Mmur_2.0_genomic.fna.gz" -O ../genomes/primata/Microcebus_murinus/Microcebus_murinus.fna.gz
gzip -d ../genomes/primata/Microcebus_murinus/Microcebus_murinus.fna.gz

# 766 Carlito syrichta (Philippine tarsier)
mkdir -p ../genomes/primata/Carlito_syrichta
wget "ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF_000164805.1_Tarsius_syrichta-2.0.1/GCF_000164805.1_Tarsius_syrichta-2.0.1_genomic.fna.gz" -O ../genomes/primata/Carlito_syrichta/Carlito_syrichta.fna.gz
gzip -d ../genomes/primata/Carlito_syrichta/Carlito_syrichta.fna.gz

# 14430 Aotus nancymaae (Ma's night monkey)
mkdir -p ../genomes/primata/Aotus_nancymaae
wget "ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF_000952055.1_Anan_1.0/GCF_000952055.1_Anan_1.0_genomic.fna.gz" -O ../genomes/primata/Aotus_nancymaae/Aotus_nancymaae.fna.gz
gzip -d ../genomes/primata/Aotus_nancymaae/Aotus_nancymaae.fna.gz

# 442 Callithrix jacchus (white-tufted-ear marmoset)
mkdir -p ../genomes/primata/Callithrix_jacchus
wget "ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF_000004665.1_Callithrix_jacchus-3.2/GCF_000004665.1_Callithrix_jacchus-3.2_genomic.fna.gz" -O ../genomes/primata/Callithrix_jacchus/Callithrix_jacchus.fna.gz
gzip -d ../genomes/primata/Callithrix_jacchus/Callithrix_jacchus.fna.gz

# 6907 Saimiri boliviensis boliviensis (Bolivian squirrel monkey)
mkdir -p ../genomes/primata/Saimiri_boliviensis
wget "ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF_000235385.1_SaiBol1.0/GCF_000235385.1_SaiBol1.0_genomic.fna.gz" -O ../genomes/primata/Saimiri_boliviensis/Saimiri_boliviensis.fna.gz
gzip -d ../genomes/primata/Saimiri_boliviensis/Saimiri_boliviensis.fna.gz

# 480 Nomascus leucogenys (northern white-cheeked gibbon)
mkdir -p ../genomes/primata/Nomascus_leucogenys
wget "ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF_000146795.2_Nleu_3.0/GCF_000146795.2_Nleu_3.0_genomic.fna.gz" -O ../genomes/primata/Nomascus_leucogenys/Nomascus_leucogenys.fna.gz
gzip -d ../genomes/primata/Nomascus_leucogenys/Nomascus_leucogenys.fna.gz

# 325 Pongo abelii (Sumatran orangutan)
mkdir -p ../genomes/primata/Pongo_abelii
wget "ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF_000001545.4_P_pygmaeus_2.0.2/GCF_000001545.4_P_pygmaeus_2.0.2_genomic.fna.gz" -O ../genomes/primata/Pongo_abelii/Pongo_abelii.fna.gz
gzip -d ../genomes/primata/Pongo_abelii/Pongo_abelii.fna.gz

# 2156 Gorilla gorilla (western gorilla)
mkdir -p ../genomes/primata/Gorilla_gorilla
wget "ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF_000151905.1_gorGor3.1/GCF_000151905.1_gorGor3.1_genomic.fna.gz" -O ../genomes/primata/Gorilla_gorilla/Gorilla_gorilla.fna.gz
gzip -d ../genomes/primata/Gorilla_gorilla/Gorilla_gorilla.fna.gz

# 51 Homo sapiens (human)
mkdir -p ../genomes/primata/Homo_sapiens
wget "ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF_000001405.31_GRCh38.p5/GCF_000001405.31_GRCh38.p5_genomic.fna.gz" -O ../genomes/primata/Homo_sapiens/Homo_sapiens.fna.gz
gzip -d ../genomes/primata/Homo_sapiens/Homo_sapiens.fna.gz

# 202 Pan troglodytes (chimpanzee)
mkdir -p ../genomes/primata/Pan_troglodytes
wget "ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF_000001515.6_Pan_troglodytes-2.1.4/GCF_000001515.6_Pan_troglodytes-2.1.4_genomic.fna.gz" -O ../genomes/primata/Pan_troglodytes/Pan_troglodytes.fna.gz
gzip -d ../genomes/primata/Pan_troglodytes/Pan_troglodytes.fna.gz

# 10729 Pan paniscus (pygmy chimpanzee)
mkdir -p ../genomes/primata/Pan_paniscus
wget "ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF_000258655.2_panpan1.1/GCF_000258655.2_panpan1.1_genomic.fna.gz" -O ../genomes/primata/Pan_paniscus/Pan_paniscus.fna.gz
gzip -d ../genomes/primata/Pan_paniscus/Pan_paniscus.fna.gz

# 36539 Colobus angolensis (Angolan colobus)
mkdir -p ../genomes/primata/Colobus_angolensis
wget "ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF_000951035.1_Cang.pa_1.0/GCF_000951035.1_Cang.pa_1.0_genomic.fna.gz" -O ../genomes/primata/Colobus_angolensis/Colobus_angolensis.fna.gz
gzip -d ../genomes/primata/Colobus_angolensis/Colobus_angolensis.fna.gz

# 7996 Rhinopithecus roxellana (golden snub-nosed monkey)
mkdir -p ../genomes/primata/Rhinopithecus_roxellana
wget "ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF_000769185.1_Rrox_v1/GCF_000769185.1_Rrox_v1_genomic.fna.gz" -O ../genomes/primata/Rhinopithecus_roxellana/Rhinopithecus_roxellana.fna.gz
gzip -d ../genomes/primata/Rhinopithecus_roxellana/Rhinopithecus_roxellana.fna.gz

# 7994 Nasalis larvatus (proboscis monkey)
mkdir -p ../genomes/primata/Nasalis_larvatus
wget "ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA_000772465.1_Charlie1.0/GCA_000772465.1_Charlie1.0_genomic.fna.gz" -O ../genomes/primata/Nasalis_larvatus/Nasalis_larvatus.fna.gz
gzip -d ../genomes/primata/Nasalis_larvatus/Nasalis_larvatus.fna.gz

# 13136 Chlorocebus sabaeus (green monkey)
mkdir -p ../genomes/primata/Chlorocebus_sabaeus
wget "ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF_000409795.2_Chlorocebus_sabeus_1.1/GCF_000409795.2_Chlorocebus_sabeus_1.1_genomic.fna.gz" -O ../genomes/primata/Chlorocebus_sabaeus/Chlorocebus_sabaeus.fna.gz
gzip -d ../genomes/primata/Chlorocebus_sabaeus/Chlorocebus_sabaeus.fna.gz

# 394 Papio anubis (olive baboon)
mkdir -p ../genomes/primata/Papio_anubis
wget "ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF_000264685.2_Panu_2.0/GCF_000264685.2_Panu_2.0_genomic.fna.gz" -O ../genomes/primata/Papio_anubis/Papio_anubis.fna.gz
gzip -d ../genomes/primata/Papio_anubis/Papio_anubis.fna.gz

# 13303 Cercocebus atys (sooty mangabey)
mkdir -p ../genomes/primata/Cercocebus_atys
wget "ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF_000955945.1_Caty_1.0/GCF_000955945.1_Caty_1.0_genomic.fna.gz" -O ../genomes/primata/Cercocebus_atys/Cercocebus_atys.fna.gz
gzip -d ../genomes/primata/Cercocebus_atys/Cercocebus_atys.fna.gz

# 36538 Mandrillus leucophaeus (drill)
mkdir -p ../genomes/primata/Mandrillus_leucophaeus
wget "ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF_000951045.1_Mleu.le_1.0/GCF_000951045.1_Mleu.le_1.0_genomic.fna.gz" -O ../genomes/primata/Mandrillus_leucophaeus/Mandrillus_leucophaeus.fna.gz
gzip -d ../genomes/primata/Mandrillus_leucophaeus/Mandrillus_leucophaeus.fna.gz

# 13267 Macaca nemestrina (pig-tailed macaque)
mkdir -p ../genomes/primata/Macaca_nemestrina
wget "ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF_000956065.1_Mnem_1.0/GCF_000956065.1_Mnem_1.0_genomic.fna.gz" -O ../genomes/primata/Macaca_nemestrina/Macaca_nemestrina.fna.gz
gzip -d ../genomes/primata/Macaca_nemestrina/Macaca_nemestrina.fna.gz

# 215 Macaca mulatta (Rhesus monkey)
mkdir -p ../genomes/primata/Macaca_mulatta
wget "ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF_000002255.3_Mmul_051212/GCF_000002255.3_Mmul_051212_genomic.fna.gz" -O ../genomes/primata/Macaca_mulatta/Macaca_mulatta.fna.gz
gzip -d ../genomes/primata/Macaca_mulatta/Macaca_mulatta.fna.gz

# 776 Macaca fascicularis (crab-eating macaque)
mkdir -p ../genomes/primata/Macaca_fascicularis
wget "ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF_000364345.1_Macaca_fascicularis_5.0/GCF_000364345.1_Macaca_fascicularis_5.0_genomic.fna.gz" -O ../genomes/primata/Macaca_fascicularis/Macaca_fascicularis.fna.gz
gzip -d ../genomes/primata/Macaca_fascicularis/Macaca_fascicularis.fna.gz

# 11864 Daubentonia madagascariensis (aye-aye)
mkdir -p ../genomes/primata/Daubentonia_madagascariensis
wget "ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/241/425/GCA_000241425.1_DauMad_1.0/GCA_000241425.1_DauMad_1.0_genomic.fna.gz" -O ../genomes/primata/Daubentonia_madagascariensis/Daubentonia_madagascariensis.fna.gz
gzip -d ../genomes/primata/Daubentonia_madagascariensis/Daubentonia_madagascariensis.fna.gz

# 44416 Cebus capucinus (white-faced sapajou)
mkdir -p ../genomes/primata/Cebus_capucinus
wget "ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/001/604/975/GCF_001604975.1_Cebus_imitator-1.0/GCF_001604975.1_Cebus_imitator-1.0_genomic.fna.gz" -O ../genomes/primata/Cebus_capucinus/Cebus_capucinus.fna.gz
gzip -d ../genomes/primata/Cebus_capucinus/Cebus_capucinus.fna.gz

# 39851 Eulemur flavifrons (Sclater's lemur)
mkdir -p ../genomes/primata/Eulemur_flavifrons 
wget "ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/001/262/665/GCA_001262665.1_Eflavifronsk33QCA/GCA_001262665.1_Eflavifronsk33QCA_genomic.fna.gz" -O ../genomes/primata/Eulemur_flavifrons/Eulemur_flavifrons.fna.gz
gzip -d ../genomes/primata/Eulemur_flavifrons/Eulemur_flavifrons.fna.gz

# 39850 Eulemur macaco (black lemur)
mkdir -p ../genomes/primata/Eulemur_macaco
wget "ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/001/262/655/GCA_001262655.1_Emacaco_refEf_BWA_oneround/GCA_001262655.1_Emacaco_refEf_BWA_oneround_genomic.fna.gz" -O ../genomes/primata/Eulemur_macaco/Eulemur_macaco.fna.gz
gzip -d ../genomes/primata/Eulemur_macaco/Eulemur_macaco.fna.gz

#  64575 Piliocolobus tephrosceles (Ugandan red Colobus)
mkdir -p ../genomes/primata/Piliocolobus_tephrosceles
wget "ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/002/776/525/GCF_002776525.1_ASM277652v1/GCF_002776525.1_ASM277652v1_genomic.fna.gz" -O ../genomes/primata/Piliocolobus_tephrosceles/Piliocolobus_tephrosceles.fna.gz
gzip -d ../genomes/primata/Piliocolobus_tephrosceles/Piliocolobus_tephrosceles.fna.gz

