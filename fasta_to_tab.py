from Bio import SeqIO
import sys

records = SeqIO.parse(sys.argv[1], "fasta")
count = SeqIO.write(records, sys.argv[1]+".tab", "tab")
print("Converted %i fasta records to tab" % count)
