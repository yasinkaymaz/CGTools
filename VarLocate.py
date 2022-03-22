#!/usr/bin/python
"""
Determining single nucleotide variant locations from multiple sequence alignment files. Pairwise! Max two sequences.
"""
import os
import pandas  as pd
import sys
from Bio import AlignIO
import Bio.Align
dir = os.path.dirname(__file__)
print(os.getcwd())
print(dir)
from collections import defaultdict, Counter
#from collections import namedtuple
import collections

alndata = []

#sys.argv[1] -> Alignment File
#sys.argv[2] -> Refname
#sys.argv[3] -> Ref repeats bed file
#sys.argv[4] ->

if len(sys.argv) < 2:
    print("Welcome to new era of MSA.")
    print("Please provide required arguments in proper order:")
    print("VarLocate.py AlignmentFile Refname")
    print("Feed with an alignment file in fasta format.")
    sys.exit(1)
Ref = sys.argv[2]
#Make sure that repeat files are in the same directory
#if sys.argv[2] == '1':
#    Ref='NC_007605'
#    Repeat_inputFile = '%(directory)s/../resources/Annotation/Type1/NC_007605_miropeat_default_run_repeats.bed' %{"directory":dir}
#else:
#    Ref='NC_009334'
#    Repeat_inputFile = '%(directory)s/../resources/Annotation/Type2/NC_009334_miropeat_default_run_repeats.bed' %{"directory":dir}



outfile1 = open(sys.argv[1]+"_SNV_positions.bed", "w")

tmpoutfile = open(sys.argv[1]+".tmp.file.aln", "w")
alignment = AlignIO.read(open(sys.argv[1]), "fasta")
recordslist = []
for record in alignment:
#    print(record.id, record.seq)
    recordslist.append(record.id)
    tmpoutfile.write(str(record.id)+"\t"+str(record.seq)+"\n")
tmpoutfile.close()

for id in recordslist:
    if id != Ref:
        Alt = id
#print(recordslist, Alt))

#with open(sys.argv[3], "r") as alnfile:
with open(sys.argv[1]+".tmp.file.aln", "r") as alnfile:
    df = pd.read_csv(alnfile,sep="\t", header=None)
    dfseq = []
    # take the names of sequences to index
    index = df[0]
    #str_seq1 = df.ix[1][1]
    str_seq1 = df.iloc[1][1]
    print("len : " +  str(len(str_seq1)))
    # create an empty list of length of sequence letters
    columns = list(range(len(str_seq1)))
    #new data frame
    new_df = pd.DataFrame(index=index, columns= columns)
    #print(new_df)
    #create new dataframe
    for i in range(len(df)):
        # change a string of sequence to a list of sequence
        #dfseq = list(df.ix[i][1])
        dfseq = list(df.iloc[i][1])
        # You can edit a subset of a dataframe by using loc:
        # df.loc[<row selection>, <column selection>]
        # place the list of sequence to the row that it belongs to in the new data frame
        new_df.loc[i:,:] = dfseq

    # for each row:
    Ref_pos = 0
    MatchCount = 0
    MissMatchCount = 0
    Insertion_l = 0
    Inserted_seq = ''
    Deletion_l = 0
    Deleted_seq = ''
    SNV_s=None
    SeqCount=len(df)
    for i in range(len(str_seq1)):
        # place the the sequence in a position to a set
        #print(collections.Counter(new_df.loc[:,i]), collections.Counter(new_df.loc[:,i]).most_common(4))
        nuccompos = collections.Counter(new_df.loc[:,i])
    #    print("n_count:", nuccompos['n'], "gap_count:", nuccompos['-'])
        del nuccompos['n']
###        del nuccompos['-']
        #print(nuccompos.most_common(4))

        set_seq = set(new_df.loc[:,i])

        if new_df.loc[Ref,i] != '-':
            Ref_pos = Ref_pos +1
        else:
            pass
        #put the bases at the location of the query sequence 1 and 2 into a set
        #set_pair = set(new_df.loc[ [sys.argv[1],sys.argv[2] ],i  ])
    #    set_pair = set(new_df.loc[ :,i  ])

        #skip indels and Ns in the alignment
    #    if 'n' in set_pair or '-' in set_pair:
    #        pass
        #print(new_df.loc[:,i].isin(['-']))
        if new_df.loc[:,i].isin(['-']).any() and Ref_pos:

            #print("A")
            if new_df.loc[Ref,i] == '-' and new_df.loc[Ref,i] == new_df.loc[Ref,i-1]:
                Insertion_l = Insertion_l + 1
                Insertion_start = i - Insertion_l +1
                Inserted_seq = Inserted_seq+str(new_df.loc[Alt,i])
                print("Insertion", Insertion_start, Insertion_l, Inserted_seq)
                SNV_s = Insertion_start
                Ref_nuc = '-'
                Alt_nuc = Inserted_seq
                continue
            else:
                Inserted_seq = ''
                Insertion_l = 0

            if new_df.loc[Alt,i] == '-' and new_df.loc[Alt,i] == new_df.loc[Alt,i-1]:
                Deletion_l = Deletion_l + 1
                Deletion_start = i - Deletion_l
                Deleted_seq = Deleted_seq+str(new_df.loc[Ref,i])
                print("Deletion", Deletion_start, Deletion_l, Deleted_seq)
                SNV_s = Deletion_start
                Ref_nuc = Deleted_seq
                Alt_nuc = '-'
            else:
                Deleted_seq = ''
                Deletion_l = 0
            continue

#        if Deleted_seq != '' or Inserted_seq != '':
#            outfile1.write(str(Ref)+"\t"+str(SNV_s)+"\t"+str(SNV_s+1)+"\t"+str(Ref_nuc)+str(Ref_pos)+str(Alt_nuc)+"\n")
        if Deleted_seq != '' and '-' not in nuccompos.most_common(4) :
            outfile1.write(str(Ref)+"\t"+str(Deletion_start)+"\t"+str(Deletion_start+Deletion_l)+"\t"+str(Ref_nuc)+str(Ref_pos)+str(Alt_nuc)+"\n")
            Deleted_seq = ''

        elif Inserted_seq != '' and '-' not in nuccompos.most_common(4) :
            print("insertion write")
            outfile1.write(str(Ref)+"\t"+str(Insertion_start)+"\t"+str(Insertion_start+Insertion_l)+"\t"+str(Ref_nuc)+str(Ref_pos)+str(Alt_nuc)+"\n")
            Inserted_seq = ''


        #if there is a mismatch error and this position is not in repeat regions, count as sequencing error.
        if len(nuccompos.most_common(4))>1 and '-' not in nuccompos.most_common(4) and Ref_pos and float(100*nuccompos.most_common(4)[1][1]/SeqCount) > 10.00:
#           print(nuccompos.most_common(4), Ref_pos, SeqCount, 100*nuccompos.most_common(4)[0][1]/SeqCount,100*nuccompos.most_common(4)[1][1]/SeqCount)

            SNV_s = Ref_pos
            Ref_nuc = new_df.loc[Ref,i]
            Alt_nuc = new_df.loc[Alt,i]

            MissMatchCount = MissMatchCount +1
            outfile1.write(str(Ref)+"\t"+str(SNV_s)+"\t"+str(SNV_s+1)+"\t"+str(Ref_nuc)+str(Ref_pos)+str(Alt_nuc)+"\n")

        #If there is a perfect match and the position is not in the repeat regions, count as match.
        elif len(nuccompos.most_common(4)) == 1 and Ref_pos:
            MatchCount = MatchCount +1






#outfile2.write(sys.argv[1]+"\t"+sys.argv[2]+"\t"+str(MissMatchCount)+"\t"+str(MatchCount)+"\t"+str(float(MissMatchCount)/(MissMatchCount+MatchCount)))
print("Match Count:", MatchCount, Ref_pos)
outfile1.close()
