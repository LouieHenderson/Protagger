#!/bin/sh
#SBATCH --ntasks=15

./bin/runProtagger.py -b ./thermdata_200707 -o FlgK_multimer_site1_200831 -p FlgEKL_partial_multimer -n 340 -c 350 -r 3 -u 15 -m 80 -t 2.5 -v 3 -x A,B,C,D,E -l > FlgK_multimer_site1_15core_o3_200831.log
