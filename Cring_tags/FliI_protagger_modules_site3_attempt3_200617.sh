#!/bin/sh
#SBATCH --ntasks=20

./bin/runProtagger.py -b ./database_200606 -o FliI_multimer_site3_deep_200620 -p HIC6_model -n 188 -c 208 -r 3 -u 20 -m 100 -t 2.5 -v 3 -x A,D,G,J,M,P -l -d > FliI_multimer_site3_20core_o3_deep_200620.log
