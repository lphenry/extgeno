#!/bin/bash

#obj: clean up fastq.gz files from unpacked sras

#cd into proper directory
cd ~/ncbi/public/sra/

#make appropriate directory to save the cleanseqs
mkdir ~/rapidevo/hardy/cleanseqs

#isolate the name of files. chop names to preserve order
#this pattern depends on if there are single or paired ends. 
#for SE: cut -c 10- for file names within directory
#for PE: cut -c 12- for file names within directory 

#standard trimming to remove trim adapters, low quality sequences using trimmomatic

for i in $(ls *.gz | rev | cut -c 10- | rev )

do

java -jar ~/programs/trimmomatic-0.36.jar SE -threads 10 ${i}.fastq.gz ~/rapidevo/hardy/cleanseqs/${i}_cleaned.fastq.gz ILLUMINACLIP:AdaptersTrim.fasta:1:30:7 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:20

done
