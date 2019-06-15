#!/bin/bash

#obj: read through each sra, convert and zip to fastq.gz, then pass to kraken. prepare kraken report

#cd into proper directory
cd ~/rapidevo/hardy/cleanseqs/

#make directory to move
mkdir ~/rapidevo/hardy/kraken_outs

#load kraken
module load kraken

#isolate the name of files. chop names to preserve order

for i in $(ls *.gz | rev | cut -c 18- | rev )

do
#to echo SRR ID 
>&2 echo ${i}

kraken --db ~/minikraken/minikraken_8gb/ --threads 10 ${i}_cleaned.fastq.gz > ~/rapidevo/hardy/kraken_outs/${i}_seqs.kraken 

done
