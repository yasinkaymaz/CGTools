#!/bin/bash
dir=`pwd`

toolDir='/Users/yasinkaymaz/Dropbox/codes/CGTools/'

#Take input dna multi fasta file. this file has to end with .fasta
InputFasta=$1

#Convert input fasta to tab file. returns "$InputFasta".tab
python $toolDir/fasta_to_tab.py "$InputFasta"


while read line;
do
  #Transform lowercase bases to UPPERCASEs only sequences
  sample_name=`echo $line|cut -d " " -f1`;
  seq=`echo $line|cut -d " " -f2`;
  newseq=`echo $seq | tr "atgcn" "ATGCN"`;
  #echo -e "$sample_name $newseq";
  #Revert tab to fasta for a given line.
  echo -e ">$sample_name\n$newseq" > dna.fa
  #translate DNA to protein.
  perl $toolDir/translateDna.pl dna.fa protein.fa
  #write the output
  cat protein.fa >> "${InputFasta%.fasta}".protein.fasta

done < "$InputFasta".tab;

rm "$InputFasta".tab protein.fa dna.fa;

echo "DNA to protein sequence conversion is complete."

#~/Documents/Tools/weblogo/seqlogo -k 0 -M  -c -n -Y  -F PDF -f epitopes.bed.aln.fasta > test.pdf
#weblogo -A protein -c hydrophobicity -l 10 -u 20 -D fasta -F pdf -f my_ICed_EBNA3C_protein.fasta -o test2.pdf
