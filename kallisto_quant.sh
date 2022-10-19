#!/usr/bin/env bash

set -euo pipefail  # useful option for bash scripts. More here https://gist.github.com/mohanpedala/1e2ff5661761d3abd0385e8223e16425

reads_base=/home/iuliiasolomennikova/Desktop/trans_practice/HW2_trans/sra_data
  # directory with trimmed reads
index=/home/iuliiasolomennikova/Desktop/trans_practice/HW2_trans/mus_musculus/transcriptome.idx  # index created with kallisto index

for i in `ls $reads_base/*.fastq | sed s/.fastq//`;  # iteration over all fastq prefixes
do
  echo 'quantifying' $i;
  kallisto quant -i $index -o "$i".kallisto --single -l 200 -s 30 "$i".fastq  2> $i.kallisto.log ;  # kallisto quant for all read pairs, consult kallisto quant --help for more options
done
