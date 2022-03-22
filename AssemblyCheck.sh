#!/bin/bash
InFasta=$1
RefFasta=$2
Prefix=$3
mummerplot=$4

if [ "$mummerplot" = '1' ];
then
    nucmer --mum $InFasta $RefFasta -p $Prefix && \
    delta-filter -l 1000 -q $Prefix.delta > $Prefix.filtered.delta && \
    mummerplot -png -f $Prefix.filtered.delta -p $Prefix;
else
    echo "no mummerplot"
fi


#Reads to _assembly:
# mkdir Bowtie2map
# cd Bowtie2map
# bowtie2-build $InFasta "${Prefix}"_idx
#
# bowtie2 -x "${Prefix}"_idx -1 s7_1.fq.gz -2 s7_2.fq.gz -S Lens_ervo_readsback.sam
#
# bash ~/Dropbox/codes/CGTools/bamprocess.sh Lens_ervo_readsback.sam
#
# #https://www.biostars.org/p/5165/
# samtools depth  Lens_ervo_readsback.sorted.bam |  awk '{sum+=$3} END { print "Average = ", sum/NR}'
# #Average =  6800.19
