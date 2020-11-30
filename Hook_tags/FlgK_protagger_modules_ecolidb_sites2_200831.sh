#!/bin/sh
#SBATCH --ntasks=15 --partition=heartofgoldCPU

./bin/runProtagger.py -b ./database_200625 -o FlgK_multimer_site2_ecolidb_200831 -p FlgEKL_partial_multimer -n 377 -c 384 -r 3 -u 15 -m 80 -t 2.5 -v 3 -x A,B,C,D,E -l > FlgK_multimer_site2_ecolidb_15core_o3_200831.log
