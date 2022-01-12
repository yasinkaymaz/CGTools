#!/usr/bin/python
#mainly from: https://bioinformatics.stackexchange.com/questions/13656/relaxed-and-sequential-phylip-format-conversion
from Bio import AlignIO
from Bio.AlignIO.PhylipIO import PhylipWriter
import sys

infile = sys.argv[1]
outfile = sys.argv[2]
idwidth = int(sys.argv[3])

with open(infile, 'r') as input_handle:
    align = AlignIO.read(input_handle, "fasta")
    with open(outfile, 'w') as output_handle:
        PhylipWriter(output_handle).write_alignment(align, id_width=idwidth)
