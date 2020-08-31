#!/bin/sh
#SBATCH --ntasks=15

./bin/runProtagger.py -b ./thermdata_200707 -o FlgL_multimer_site2_200831 -p FlgEKL_partial_multimer -n 464 -c 493 -r 3 -u 15 -m 80 -t 2.5 -v 3 -x L,M,N,O,P,Q -l > FlgL_multimer_site2_15core_o3_200831.log
