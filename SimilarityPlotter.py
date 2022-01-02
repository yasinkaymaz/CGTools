#!/usr/bin/python

import os
import pandas  as pd
import sys
from Bio import AlignIO
import Bio.Align
dir = os.path.dirname(__file__)
print os.getcwd()
print dir
alndata = []
#sys.argv[1] -> Ref name
#sys.argv[2] -> Sample genome
#sys.argv[3] -> Alignment File
#sys.argv[4] -> Repeat file of reference
#
if len(sys.argv) < 4:
	print "This code computes the similarities between a given genome and its reference using their multiple sequence alignments."
	print "Please provide required arguments in proper order:"
	print "python SimilarityPlotter.py RefName GenomeName AlignmentFile [Repeat file]"
	print "Feed with an alignment file in fasta format."
	sys.exit(1)

RefName=sys.argv[1]
GenomeName = sys.argv[2]

# Ref=''
# #Make sure that repeat files are in the same directory
# if sys.argv[4] == '1':
InputAlign=sys.argv[3]
Repeat_inputFile = sys.argv[4]
# else:
#Ref2='NC_009334'
# 	Repeat_inputFile = '%(directory)s/../resources/Annotation/Type2/NC_009334_miropeat_default_run_repeats.bed' %{"directory":dir}
#
# #Find the genomic locations fall into miropeats repeat regions
RepeatList = []
with open(Repeat_inputFile) as repeatFile:
	for line in repeatFile:
		repstart = int(line.strip().split("\t")[1])
		repend = int(line.strip().split("\t")[2])
		for i in range(repstart,repend):
			RepeatList.append(i)
RepeatList = set(RepeatList)
print len(set(RepeatList))
#
#
outfile1 = open(GenomeName+"_Mismatch_positions.bed", "w")

# outfile2 = open(sys.argv[1]+sys.argv[2]+"_Mismatch_error_rates.txt","w")
#
tmpoutfile = open(InputAlign+".tmp.file.aln","w")
alignment = AlignIO.read(open(InputAlign), "fasta")
for record in alignment:
	#print record.id, record.seq
	tmpoutfile.write(str(record.id)+"\t"+str(record.seq)+"\n")
tmpoutfile.close()
#
#
UncoveredPos_1=[]
#UncoveredPos_2=[]
with open(InputAlign+".tmp.file.aln", "r") as alnfile:
	df = pd.read_csv(alnfile,sep="\t",header=None)
	dfseq = []
#     # take the names of sequences to index
	index = df[0]
	str_seq1 = df.ix[1][1]
	print "len : " +  str(len(str_seq1))
#     # create an empty list of length of sequence letters
	columns = list(range(len(str_seq1)))
#     #new data frame
	new_df = pd.DataFrame(index=index, columns= columns)
    #create new dataframe
	for i in range(len(df)):
#         # change a string of sequence to a list of sequence
		dfseq = list(df.ix[i][1])
#         # You can edit a subset of a dataframe by using loc:
#         # df.loc[<row selection>, <column selection>]
#         # place the list of sequence to the row that it belongs to in the new data frame
		new_df.loc[i:,:] = dfseq

#
#     # for each row:
	Ref_pos = 0
#	MatchCount = 0
	MissMatchCount = 0
	for i in range(len(str_seq1)):
#         # place the the sequence in a position to a set
		set_seq = set(new_df.loc[:,i])

		if new_df.loc[RefName, i] != '-':
			Ref_pos = Ref_pos +1
		else:
			pass
# 		#put the bases at the location of the query sequence 1 and 2 into a set
		set_pair = set(new_df.loc[ [GenomeName, RefName ],i  ])

		if 'N' in set_pair or 'n' in set_pair:
			UncoveredPos_1.append(Ref_pos)
			print set_pair
		#if there is a mismatch error and this position is not in repeat regions, count as sequencing error.
		elif len(set_pair)>1 and Ref_pos not in RepeatList and '-' not in set_pair:
			print set_pair, Ref_pos
			MissMatchCount = MissMatchCount +1
			outfile1.write(str(RefName)+"\t"+str(Ref_pos)+"\t"+str(Ref_pos+1)+"\n")

		######## Do the same things for Ref2, which is type 2 reference genome.
		# set_pair = set(new_df.loc[ [sys.argv[1],Ref2 ],i  ])
		# if 'N' in set_pair or 'n' in set_pair:
		# 	UncoveredPos_2.append(Ref_pos)
		# 	print set_pair
		# #if there is a mismatch error and this position is not in repeat regions, count as sequencing error.
		# elif len(set_pair)>1 and Ref_pos not in RepeatList and '-' not in set_pair:
		# 	print set_pair, Ref_pos
		# 	MissMatchCount = MissMatchCount +1
		# 	outfile2.write(str(Ref1)+"\t"+str(Ref_pos)+"\t"+str(Ref_pos+1)+"\n")
#

outfile1.close()
#outfile2.close()
print "Now smoothing"
outfile1 = open(GenomeName+"_Mismatch_positions_smoothed.bed","w")
#outfile2 = open(sys.argv[1]+"_Mismatch_positions_with_type2_smooth.bed","w")

#smoothing:
#Define a window, typically 100 bases.
win=600
#Define a smoothing interval.
sm=200
genomeLen=Ref_pos
#genomeLen=171725

#Read the list of genomic positions that are mismatch between the aligned two genomes. Relative to Reference.
MismatchPositions =[]
with open(GenomeName+"_Mismatch_positions.bed","r") as snpfile:
    for line in snpfile:
        MismatchPositions.append(int(line.strip().split("\t")[1]))

outfile1.write("Chrom"+"\t"+"start"+"\t"+"end"+"\t"+"Sim2Ref"+"\n")
for i in range(0,genomeLen-win+sm,sm):
    windowPoss = range(i, i+win)
    #print windowPoss
    dissim=len(set(windowPoss)&set(MismatchPositions))
    CoveredLen=win-len(set(windowPoss)&set(UncoveredPos_1))+1
    PercentDis = 100*float(dissim)/CoveredLen

    #print "CoveredLen: ", CoveredLen, "dissim: ", dissim, "PercentDis: ", PercentDis
    #print i, dissim
#    outfile1.write("NC_007605"+"\t"+str(i)+"\t"+str(i+sm)+"\t"+str(100-dissim/10)+"\n")
    outfile1.write(RefName+"\t"+str(i)+"\t"+str(i+sm)+"\t"+str(PercentDis)+"\n")

###### ----- ######
#Do the same things for Ref2, which is type 2 reference genome.
#
# MismatchPositions =[]
# with open(sys.argv[1]+"_Mismatch_positions_with_type2.bed","r") as snpfile:
#     for line in snpfile:
#         MismatchPositions.append(int(line.strip().split("\t")[1]))
#
# outfile2.write("Chrom"+"\t"+"start"+"\t"+"end"+"\t"+"Sim2Type2"+"\n")
# for i in range(0,genomeLen-win+sm,sm):
#     windowPoss = range(i, i+win)
#     #print windowPoss
#     dissim=len(set(windowPoss)&set(MismatchPositions))
#     CoveredLen=win-len(set(windowPoss)&set(UncoveredPos_1))+1
#     PercentSim=100-(100*dissim/CoveredLen)
#     #print i, dissim
#     outfile2.write("NC_007605"+"\t"+str(i)+"\t"+str(i+sm)+"\t"+str(PercentSim)+"\n")


outfile1.close()
#outfile2.close()
