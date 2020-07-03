#!/bin/sh
#SBATCH --ntasks=50

./bin/runProtagger.py -b ./database_200606 -o FliI_multimer_site2_200616 -p HIC6_model -n 378 -c 383 -r 3 -u 50 -t 2.5 -v 3 -x A,D,G,J,M,P -l > FliI_multimer_site2_50core_03_200616.log
