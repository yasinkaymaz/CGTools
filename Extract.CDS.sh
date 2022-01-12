
FASTA=$1
GFF=$2

gffread -C $GFF -g $FASTA -x CDS.fa



for i in `grep ">" CDS.fa|cut -b 2-`;
do
    GENEID=$(grep -w $i $GFF|awk '($3 == "gene")'|cut -f9|tr ';' '\n'|grep "Name="|cut -b 6- |uniq)
    echo $GENEID
    sed "s/$i/$GENEID/" CDS.fa > genes.CDS.fa
    mv genes.CDS.fa CDS.fa

done

#bash ~/Dropbox/codes/CGTools/Extract.CDS.sh NC_028171.1.Astragalus_nakaianus.fasta Astragalus_nakaianus.gff3
