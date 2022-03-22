#!/bin/bash

#This input genome MSA fasta file has to have a Reference genome.
GenomeMSAfasta=$1
# This script also takes a bed file for "CDS regions" of each gene.
#This bed should have 6 columns of which 4th is gene names and 6the is the strand of the transcript.
CDSbed=$2
#Specify the Reference genome name that matches the fasta sequence
RefName=$3
outDir=$4
#cd $outDir/


toolDir='/Users/yasinkaymaz/Dropbox/codes/CGTools'
DistanceCalc=1

use-java () {
  export JAVA_HOME=`/usr/libexec/java_home -v 1.$1`
}


if [ "$DistanceCalc" = "1" ]
then
  #do some stuff
  #Create an array of genomic positions of repetitive and non-coding regions.
  #module load bedtools/2.17.0

  #echo "" > pairwiseMeandNdS.txt
  echo -e "#Gene\tNum_Tot_Vars\tNum_N\tNum_S\N-Changes" > $outDir/Genes.Ns.Ss.txt
  #rm *pairwisedistances*
  #rm *pairwiseMeandNdS*
  #for each coding gene
  for geneOfvariant in $(cut -f4 $CDSbed|sort|uniq);
  do
    echo $geneOfvariant;
    grep -w $geneOfvariant $CDSbed > $outDir/"${geneOfvariant}".bed;#Change with awk $4 == $geneOfvariant
    strand=`cut -f6 $outDir/"$geneOfvariant".bed|uniq`;
    echo -e "Strand is $strand";

    #Extract gene sequence from Ref sequence.
    python $toolDir/MSA_gene_extractor.py $GenomeMSAfasta $RefName $outDir/"$geneOfvariant".bed

    #<--->
    echo "!--- gene extracted"
    ## this returns "$geneOfvariant".bed.aln.fasta.tab
    python $toolDir/fasta_to_tab.py $outDir/"$geneOfvariant".bed.aln.fasta
    echo "!--- file converted to tab"
    #is gene expressed from minus transcript? if yes, take reverse transcript of multi-fasta file.
    if [ "$strand" = "-" ]
    then
      #if gene is from reverse strand:
      while read ll;
      do
        sample_name=`echo $ll|cut -d " " -f1`;
        seq=`echo $ll|cut -d " " -f2`;
        rcseq=`echo $seq|rev | tr "ATGCNatgcn" "TACGNTACGN"`;
        echo -e "$sample_name\t$rcseq";
      done < $outDir/"$geneOfvariant".bed.aln.fasta.tab > $outDir/"$geneOfvariant".aln.rc.tab;

      awk '{print ">"$1"\n"$2 }' $outDir/"$geneOfvariant".aln.rc.tab > $outDir/"$geneOfvariant".aln.rc.fasta

      genemultifasta=$outDir/"$geneOfvariant".aln.rc.fasta

    #if no, keep multi-fasta as it is. just convert to upper case
    else

      while read ll;
      do
        sample_name=`echo $ll|cut -d " " -f1`;
        seq=`echo $ll|cut -d " " -f2`;
        rcseq=`echo $seq| tr "ATGCNatgcn" "ATGCNATGCN"`;
        echo -e "$sample_name\t$rcseq";
      done < $outDir/"$geneOfvariant".bed.aln.fasta.tab > $outDir/"$geneOfvariant".aln.c.tab

      awk '{print ">"$1"\n"$2 }' $outDir/"$geneOfvariant".aln.c.tab > $outDir/"$geneOfvariant".aln.fasta

      genemultifasta=$outDir/"$geneOfvariant".aln.fasta

    fi
    #<--->

    bash $toolDir/DNA2Protein.sh $genemultifasta

    #Calculate d-Kimura two parameter for All genomes
#    python $toolDir/bin/MSA_distanceCalc.py MeanPairwiseKimuraDist -af $genemultifasta -rn $geneOfvariant >> pairwisedistances.txt
#    python $toolDir/MSA_distanceCalc.py MeandNdS -af $genemultifasta -rn $geneOfvariant >> pairwiseMeandNdS.txt

    python $toolDir/VarLocate.py "${genemultifasta%.fasta}".protein.fasta $RefName
    Ns=`cat "${genemultifasta%.fasta}".protein.fasta_SNV_positions.bed|wc -l|sed 's/ //g'`
    list_Ns=`cat "${genemultifasta%.fasta}".protein.fasta_SNV_positions.bed|cut -f4 |tr '\n' ','|sed 's/.$//'`

    python $toolDir/VarLocate.py "$genemultifasta" $RefName
    SNVs=`cat "$genemultifasta"_SNV_positions.bed|wc -l|sed 's/ //g'`

    let Ss=SNVs-Ns;

    echo -e "$geneOfvariant\t$SNVs\t$Ns\t$Ss\t$list_Ns" >> $outDir/Genes.Ns.Ss.txt

    #Clean up intermediate files.
    rm $outDir/"$geneOfvariant".bed $outDir/"$geneOfvariant".aln.*c.tab $outDir/"$geneOfvariant".aln*.fasta.tmp.file.aln $outDir/"$geneOfvariant".bed.aln.fasta.tab $outDir/"$geneOfvariant".bed.aln.fasta

  done

else
  echo "skip gene extracting..."
fi
