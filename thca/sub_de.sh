#!/bin/bash

#PBS -N czysz_tcga
#PBS -S /bin/bash
#PBS -l walltime=30:00:00
#PBS -l nodes=1:ppn=4
#PBS -l mem=32gb

#PBS -o $HOME/thca.out
#PBS -e $HOME/thca.err

module load R/3.2.0

#export TEMPDIR=/home/t.cri.cczysz/tmp
time Rscript /home/t.cri.cczysz/gdc/thca/count_DE.R
