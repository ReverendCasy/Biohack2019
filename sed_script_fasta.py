from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import subprocess
import os
import re

os.chdir('/home/reverend_casy/Biohack/Hs_APN')
default_file = 'uniprot.fasta'


write_list = []
for entry in SeqIO.parse(default_file, 'fasta'):
    descr = entry.description
    if re.search(f'OS=Homo', descr) and len(entry.seq) > 900:
        print(descr, len(entry.seq))
        write_list.append(entry)

print(len(write_list))
subprocess.run(f"cat /dev/null > {default_file}", shell=True)
SeqIO.write(write_list, default_file, 'fasta')

