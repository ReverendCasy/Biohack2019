import sys
from Bio import Entrez,SeqIO
import csv
import time 
import math

with open(sys.argv[1], 'r') as csv_file:
    my_reader = csv.reader(csv_file, delimiter='\t') 
    ids = [row[2] for row in my_reader]
    for step in range(int(math.ceil(len(ids)/3.0))):
        if step == 0:
            start=0
            stop=3
        elif step < 3:
            start=step*3 
            stop=(step+1)*3  
        else:
            start=step*3 
            stop=len(ids)+1 

        print(ids[start:stop])
        uidList = ','.join(ids[start:stop]) 
        handle = Entrez.efetch(db="nucleotide", id=uidList, rettype="fasta", retmode="text")
        print(SeqIO.parse(handle, 'fasta')	)
        time.sleep(3)
        SeqIO.write(list(SeqIO.parse(handle, 'fasta')),'etest{}.fasta'.format(step), "fasta")
