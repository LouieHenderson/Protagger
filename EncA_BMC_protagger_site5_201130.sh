#!/bin/sh
#SBATCH --ntasks=18 --partition=heartofgoldCPU

./bin/runProtagger.py -b ./Databases/database_200625 -o EncA_BMC_site5_201130 -p 4pt2_sub -n 130 -c 163 -r 3 -u 17 -t 2.5 -v 3 -l > EncA_BMC_site5_17core_201130.log
