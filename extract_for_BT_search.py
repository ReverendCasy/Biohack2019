import sys
from Bio import SeqIO


recs =list()
for record in SeqIO.parse(open(sys.argv[1]),"fasta"):
    if 'ACQ99547.1' not in record.id and 'OTW51788.1' not in record.id and 'OTX69600.1' not in record.id and 'OTY54786.1' not in record.id and 'WP_060703593.1' not in record.id and 'WP_103570866.1' not in record.id:
        recs.append(record)

SeqIO.write(recs,'seq_iter.faa', "fasta")
