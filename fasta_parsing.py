import sys
from Bio import Entrez,SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna
from Bio.SeqRecord import SeqRecord
import csv
from collections import defaultdict

dat_dict = defaultdict(list)
req_dict = dict()
for req in list(SeqIO.parse(open(sys.argv[2]),"fasta")):
    req_dict[req.id]= req.seq

with open(sys.argv[1], 'r') as csv_file:
    my_reader = csv.reader(csv_file, delimiter='\t') 
    for row in my_reader:
        dat_dict[row[7]+row[8].replace(' ', '_')]=[row[2], row[3], row[4], row[5] ,row[6],row[7]]
wrongids=list()
for nucl in dat_dict:
    try:
    #print(dat_dict[nucl][0], dat_dict[nucl][1], dat_dict[nucl][2])
    #print(req_dict[dat_dict[nucl][0]][int(dat_dict[nucl][1]):int(dat_dict[nucl][2])+1])
        if dat_dict[nucl][3]=='+':
            obj = Seq(str(req_dict[dat_dict[nucl][0]][int(dat_dict[nucl][1])-1:int(dat_dict[nucl][2])]), generic_dna)
        else:
            obj = Seq(str(req_dict[dat_dict[nucl][0]][int(dat_dict[nucl][1])-1:int(dat_dict[nucl][2])]), generic_dna).reverse_complement()
    #record = SeqRecord(obj,id=('_').join(nucl.split()).replace('.','_').replace('/','_'))
    #prot_rec=SeqRecord(obj.translate(),id=('_').join(nucl.split()).replace('.','_').replace('/','_'))
        record = SeqRecord(obj,id=dat_dict[nucl][0],description = nucl)
        prot_rec=SeqRecord(obj.translate(),id=dat_dict[nucl][4],description = nucl)
        SeqIO.write(record,'{}.fna'.format(('_').join(nucl.split()).replace('.','_').replace('/','_')), "fasta")
        SeqIO.write(prot_rec,'{}.faa'.format(('_').join(nucl.split()).replace('.','_').replace('/','_')), "fasta")
    except:
        wrongids.append(nucl)
with open('wrong_ids.txt', 'w') as f:
    for item in wrongids:
        f.write("%s\n" % item)
