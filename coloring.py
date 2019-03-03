import sys
from Bio import Entrez,SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna, generic_protein
from Bio.SeqRecord import SeqRecord
import csv
from collections import defaultdict
from Bio import Align


coord_dict = defaultdict(dict)
req_dict = dict()

for req in list(SeqIO.parse(open(sys.argv[1]),"fasta")):
    req_dict[req.id] = req.seq
    for ind in range(len(req.seq)):
        if ind not in coord_dict:
            coord_dict[ind]=dict()
       # print(ind)
        if req.seq[ind] not in coord_dict[ind]:
            coord_dict[ind][req.seq[ind]]=1 
        else:
            coord_dict[ind][req.seq[ind]]+=1
for pos in coord_dict:
    for letter in coord_dict[pos]:
        coord_dict[pos][letter]=coord_dict[pos][letter]/float(len(coord_dict))
#print(coord_dict)

#init_list=list()
#for req in list(SeqIO.parse(open(sys.argv[2]),"fasta")):
#    init_list.append(req)

query = list(SeqIO.parse(open(sys.argv[2]),"fasta"))[0].seq
#print(query)


aln_query=req_dict['Cry1Aa2-domain_3']
print(aln_query)

flag=251
with open('{}.txt'.format('Cry1Aa2-domain_3'),'w') as my_file:
    for ind in range(len(aln_query)):
        if aln_query[ind] != '-':
            my_file.write('set_color {}= ['.format('color' +str(flag)) + str(round(coord_dict[ind][aln_query[ind]]*1.5,4)) +',' + str((1- round(coord_dict[ind][aln_query[ind]],4))/4.0)+ ','+str((1-round(coord_dict[ind][aln_query[ind]],4))/4.0)  + ']'+ '\n')
            my_file.write('color color' +str(flag) +', resi '+ str(flag) + '\n')
            flag+=1
       
my_file.close()

#color $color, resi $pos

