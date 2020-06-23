#!/bin/sh
#SBATCH --ntasks=20

./bin/runProtagger.py -b ./database_200606 -o FliI_multimer_site3_attempt2_200620 -p HIC6_model -n 196 -c 204 -r 3 -u 20 -m 100 -t 2.5 -v 2 -x A,D,G,J,M,P -l > FliI_multimer_site3_20core_o2_200620.log
