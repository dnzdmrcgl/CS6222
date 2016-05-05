#!/bin/bash
#$ -N sql7
#$ -S /bin/bash
#$ -pe OpenMP 1
#$ -cwd
#$ -o /home/ddemircioglu/CS6222/sql7.out
#$ -e /home/ddemircioglu/CS6222/sql7.err
#$ -M ddemircioglu@gis.a-star.edu.sg
#$ -m bea
#$ -q medium.q
#$ -l mem_free=40G
#$ -l h_rt=40:00:00
#$ -V

mysql -h genome-mysql.cse.ucsc.edu -u genome -A -D hg19 --skip-column-names < real2_snps.sql > real2_snps.txt
