#!/bin/sh
#SBATCH --ntasks=15

./bin/runProtagger.py -b ./database_200625 -o FliG_multimer_site1_200625 -p C-ring_portion_15mer -n 255 -c 262 -r 3 -u 15 -m 80 -t 2.5 -v 3 -x I,J,K,L,M -l > FliG_multimer_site1_15core_o3_200625.log
