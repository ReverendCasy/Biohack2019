import sys
from Bio import SeqIO
import subprocess

print(sys.argv[1][:-6]+'_orf.fasta')

cmd = subprocess.call('./ORFfinder -in {0} -out {1}'.format(sys.argv[1], sys.argv[1][:-6]+'_orf.fasta'), shell=True)

max_rec = list(SeqIO.parse(open(sys.argv[1][:-6]+'_orf.fasta'),"fasta"))[0]
for record in SeqIO.parse(open(sys.argv[1][:-6]+'_orf.fasta'),"fasta"):
    if len(record.seq)>len(max_rec.seq):
        max_rec=record
print(max_rec)
SeqIO.write(max_rec, "{}".format(sys.argv[1][:-6]+'_l.fasta'), "fasta")
