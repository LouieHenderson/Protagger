#!/bin/sh
#SBATCH --ntasks=15

./bin/runProtagger.py -b ./database_200625 -o FliM_multimer_site4_200628 -p C-ring_portion_15mer -n 94 -c 101 -r 3 -u 15 -m 80 -t 2.5 -v 3 -x W,X,Y,Z,a -l > FliM_multimer_site4_15core_o3_200628.log
