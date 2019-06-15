#!/bin/bash

#obj: generating braken reports. taking each kraken output, making kraken report, and then estimating abundance with braken. ending w/abundances

#cd into proper directory
cd ~/rapidevo/hardy/kraken_outs/

#mkdir for abund outs
mkdir ~/rapidevo/hardy/kraken_outs/abund

#load kraken
module load kraken

#isolate the name of files. 

for i in $(ls *.kraken | rev | cut -c 13- | rev )

do

kraken-report --db ~/minikraken/minikraken_8gb/ ~/rapidevo/hardy/kraken_outs/${i}_seqs.kraken > ~/rapidevo/hardy/kraken_outs/abund/${i}_report.txt

done

#move to where estimating abundances and call each one

cd ~/rapidevo/hardy/kraken_outs/abund

#isolate new names in abund dir

for i in $(ls *report.txt | rev | cut -c 12- | rev)

do

python ~/kraken_db/est_abundance.py -i ${i}_report.txt -k ~/minikraken/minikraken_8GB_75mers_distrib.txt -l 'F' -t 5 -o ${i}_fam_abund.txt

done
