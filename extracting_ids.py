import sys
from Bio import SeqIO

recs_id =set()
out_recs =dict()
for record in SeqIO.parse(open(sys.argv[1]),"fasta"):
    if "lcl|" in record.id or '|' not in record.id :
        recs_id.add(record.id.split('/')[0])
    else:
        recs_id.add(record.id.split('/')[0].split('|')[1])
#print(recs_id)

with open('hmm_ids.txt', 'w') as f:
    for item in list(recs_id):
        f.write("%s\n" % item)

for record in SeqIO.parse(open(sys.argv[2]),"fasta"):
    if record.id in recs_id:
        out_recs.append(record)
SeqIO.write(out_recs,'extr_D1.faa', "fasta")
print(recs_recs_idid)

for req in list(open(sys.argv[1])):
	out_recs[req]=1
for req in list(open(sys.argv[2])):
	if req in out_recs:
            out_recs[req]+=1 
        else:
            out_recs[req]=1
for req in list(open(sys.argv[3])):
	if req in out_recs:
            out_recs[req]+=1 
        else:
            out_recs[req]=1
s=0
for rec in out_recs:
    if out_recs[rec]==3:
       print(rec)
 
       s+=1
