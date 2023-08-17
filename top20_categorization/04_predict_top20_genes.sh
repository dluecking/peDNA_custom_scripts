#!/bin/bash -l
# needs bioinf env

for i in `cat helper_files/top20_contigs_unique.txt`;
do
    if [ ! -f top20/genes/${i}.genes.fasta ]
    then
        prodigal -i contigs/${i} -d top20/genes/${i}.genes.fasta -a top20/proteins/${i}.proteins.faa -p meta
    else
        echo "Already done!"
    fi
done
