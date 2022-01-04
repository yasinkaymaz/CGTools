
#bowtie2 -x Cari.idx -1 s4_1.fq.gz s4_2.fq.gz -S echi2ari.sam
#bowtie2 -x Cechi.idx -1 s4_1.fq.gz s4_2.fq.gz -S echi2echi.sam

InSam=$1

#sam prosessing

#SAM to BAM
samtools view -S -b "$InSam" > "${InSam%.sam}".bam

samtools view -b -F 4 "${InSam%.sam}".bam > "${InSam%.sam}".F4.bam

#Sorting
samtools sort "${InSam%.sam}".F4.bam -o "${InSam%.sam}".sorted.bam

#Indexing
samtools index "${InSam%.sam}".sorted.bam

#Compute the depth of coverage
samtools depth  "${InSam%.sam}".sorted.bam |  awk '{sum+=$3} END { print "Average = ", sum/NR}'

#Subsampling 5%
samtools view -s 0.05 -b "${InSam%.sam}".sorted.bam > "${InSam%.sam}".sorted.0_05.bam

samtools index "${InSam%.sam}".sorted.0_05.bam
