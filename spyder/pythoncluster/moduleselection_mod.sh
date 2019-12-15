#!/bin/bash

#MSUB -l nodes=1:ppn=1,mem=100gb,walltime=48:00:00
#MSUB -M s091d230@ku.edu
#MSUB -m abe
#MSUB -N LNrate7132018

module load  python/3.6 
python LNrate7132018.py