#!/bin/bash

#MSUB -l nodes=1:ppn=1,mem=16gb,walltime=980:00:00
#MSUB -M s091d230@ku.edu
#MSUB -m ae
#MSUB -d /home/s091d230/LangevinPython
#MSUB -e /home/s091d230/LangevinPython/${PBS_JOBNAME}-${PBS_JOBID}.err
#MSUB -o /home/s091d230/LangevinPython/${PBS_JOBNAME}-${PBS_JOBID}.out

module load  python/3.6 
python Langevindemon.py