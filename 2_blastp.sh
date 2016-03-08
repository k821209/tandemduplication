makeblastdb -in ${1} -dbtype prot
blastp -query ${1} -db ${1} -evalue 1e-5 -outfmt 7 -out ${1}.bp.ev1e5.out7 -num_threads 5 

