#!/usr/bin/python

import os
import pandas as pd
import sys
import collections
from Bio import AlignIO
import Bio.Align

alndata = []
#sys.argv[1] -> Alignment File
#sys.argv[2] -> type of virus

if len(sys.argv) < 3:
	print("Please provide required arguments in proper order:")
	print("MSA_parser_cleaner.py AlignmentFile (fasta format) RefName Gene_bed_file (bed3 format)")
	sys.exit(1)

GeneBedfile = sys.argv[3]

Ref=sys.argv[2]

#if sys.argv[3] == "RepeatFilter":

#Find the genomic locations fall into miropeats repeat regions
#	RepeatList = []
#	with open(Repeat_inputFile) as repeatFile:
#		for line in repeatFile:
#			repstart = int(line.strip().split("\t")[1])
#			repend = int(line.strip().split("\t")[2])
#			for i in range(repstart,repend):
#				RepeatList.append(i)
#	RepeatList = set(RepeatList)
#	print len(set(RepeatList))
#else:
#	RepeatList = []
#	pass

#compare = lambda x, y: collections.Counter(x) == collections.Counter(y)

def PullLoci(BedFile):
	"""
	this function takes a bed for multiple loci relative to a reference genome
	and dumps a list of genomic positions.
	"""
	PositionList = []
	with open(BedFile) as bfile:
		for line in bfile:
			start = int(line.strip().split("\t")[1])
			end = int(line.strip().split("\t")[2])
			for i in range(start,end+1):
				PositionList.append(i)
	PositionList = set(PositionList)

	return PositionList



alignment = AlignIO.read(open(sys.argv[1]), "fasta")
print("Alignment length %i" % alignment.get_alignment_length())

#for record in alignment :
#    print(record.seq + " " + record.id)
#    print(record.id), record.seq.count('-'), record.seq.count('n')

#good alignment positions = gap
#gap = []
#bad alignment positions = bap
#bap = []

edited = alignment[:,1:2]
#print edited

TargetGenesLoci = PullLoci(GeneBedfile)

Ref_pos = 0

#List for alignment positions for target gene
targetGeneAlignment=[]
for record in alignment:
	if record.id == Ref:
		for i in range(len(record.seq)):
			if record.seq[i] == '-':
				Ref_pos = Ref_pos
			else:
				Ref_pos = Ref_pos + 1
			if Ref_pos in TargetGenesLoci:
				targetGeneAlignment.append(i)
			else:
				pass


#print "GAP:", len(gap), "BAP:", len(bap)
print("Target gene loci length:", len(targetGeneAlignment))

for i in range(alignment.get_alignment_length()):

	#alignment Position Base Composition = apbc
#	apbc = []

#	for record in alignment:
#		apbc.append(record.seq[i])

	if i in set(targetGeneAlignment):
		edited = edited + alignment[:,i:i+1]
	else:
		pass

#	if compare(set(apbc), set(['-','n']) ) or compare(set(apbc), set(['-']) ) or compare(set(apbc), set(['n']) ) or (i in set(bap)):
		#print record.seq[i], record.id
#		pass
#	else:
#		edited = edited + alignment[:,i:i+1]


edited = edited[:,1:]
#print edited

#>>> align.add_sequence("Gamma", "ACTGCTAGATAG")

#AlignIO.write(edited, "my_TargetGeneAlignment_"+sys.argv[1], "fasta")
#AlignIO.write(edited, os.path.basename(sys.argv[3])+".aln.fasta", "fasta")
AlignIO.write(edited, sys.argv[3]+".aln.fasta", "fasta")
#print len(bap)
#print len(gap)
#outfile1 = open(sys.argv[1]+sys.argv[2]+"_Mismatch_error_positions.bed","w")
#outfile2 = open(sys.argv[1]+sys.argv[2]+"_Mismatch_error_rates.txt","w")





#outfile1.close()
#outfile2.close()
