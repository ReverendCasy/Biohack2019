import sys
from Bio import Entrez,SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna, generic_protein
from Bio.SeqRecord import SeqRecord
import csv
from collections import defaultdict

req_dict = defaultdict(list)
coord_dict= defaultdict(list)


for req in list(open(sys.argv[1]))[1:]:
        req_dict[req.strip()]= list([None]*5)
for req in list(SeqIO.parse(open(sys.argv[2]),"fasta")):
    if "lcl|" in req.id or '|' not in req.id :
           idr = req.id.split('/')[0]
    else:
           idr = req.id.split('/')[0].split('|')[1]
    if idr in req_dict:
        req_dict[idr][0]= str(req.seq)
        req_dict[idr][4]= 0
new_dict = defaultdict(list)
for req in req_dict:
    if req_dict[req][0] is not None:
        new_dict[req]=req_dict[req]


for record in list(SeqIO.parse(open(sys.argv[3]),"fasta")):
    try:
        if "lcl|" in record.id or '|' not in record.id :
            idr = record.id.split('/')[0]
            new_dict[idr][1]=str(record.seq)
            new_dict[idr][4]+=1
            #coord_dict[idr].extend(record.id.split('/')[1].split('-'))
        else:
            idr = record.id.split('/')[0].split('|')[1]
            new_dict[idr][1]=str(record.seq)
            new_dict[idr][4]+=1
            #coord_dict[idr].extend(record.id.split('/')[1].split('-'))
        #print(new_dict[idr])
    except:
        pass

for record in list(SeqIO.parse(open(sys.argv[4]),"fasta")):
    try:
        if "lcl|" in record.id or '|' not in record.id :
            idr = record.id.split('/')[0]
            new_dict[idr][2]=str(record.seq)
            new_dict[idr][4]+=1
            coord_dict[idr].extend(record.id.split('/')[1].split('-'))
        else:
            idr = record.id.split('/')[0].split('|')[1]
            new_dict[idr][2]=str(record.seq)
            new_dict[idr][4]+=1
            coord_dict[idr].extend(record.id.split('/')[1].split('-'))
       # print(new_dict[idr])
    except:
        pass

for record in list(SeqIO.parse(open(sys.argv[5]),"fasta")):
    try:
        if "lcl|" in record.id or '|' not in record.id :
            idr = record.id.split('/')[0]
            new_dict[idr][3]=str(record.seq)
            new_dict[idr][4]+=1
            coord_dict[idr].extend(record.id.split('/')[1].split('-'))
        else:
            idr = record.id.split('/')[0].split('|')[1]
            new_dict[idr][3]=str(record.seq)
            new_dict[idr][4]+=1
            coord_dict[idr].extend(record.id.split('/')[1].split('-'))
    except:
        pass

substr = list()
for req in new_dict:
	if len(req_dict[req])>=5 and req_dict[req][4]>=3:
            try:
                start = req_dict[req][0].index(req_dict[req][1])
                stop = req_dict[req][0].index(req_dict[req][3])+len(req_dict[req][3])
                substr.append(SeqRecord(Seq(req_dict[req][0][start:stop],generic_protein),id=req))
            except:
                start = min([int(x) for x in coord_dict[req]])-1
                stop = max([int(x) for x in coord_dict[req]])-1
                substr.append(SeqRecord(Seq(req_dict[req][0][start:stop],generic_protein),id=req))
                pass
print(substr)
SeqIO.write(substr,'C2+C3_domains_extr.faa', "fasta")
