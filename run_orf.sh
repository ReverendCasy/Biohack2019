#!/bin/bash

DATENOW=$( date )
echo "Orfs finding started at: $DATENOW"

for FQ in *.fasta
do
	echo "scanning $FQ"
        TMPTAG=${FQ%%.fasta}
	python largest_orf.py $FQ
        rm $TMPTAG"_orf.fasta"
	mv $FQ fasta/
done
DATENOW=$( date )
echo "Orfs finding finished at: $DATENOW"

