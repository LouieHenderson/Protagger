#!/bin/sh
#SBATCH --ntasks=18 --partition=heartofgoldCPU

./bin/runProtagger.py -b ./Databases/database_200625 -o EncA_BMC_site3_201130 -p 4pt2_sub -n 218 -c 228 -r 3 -u 17 -t 2.5 -v 3 -l > EncA_BMC_site3_17core_201130.log
