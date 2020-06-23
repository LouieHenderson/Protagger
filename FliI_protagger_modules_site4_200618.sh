#!/bin/sh
#SBATCH --ntasks=20

./bin/runProtagger.py -b ./database_200606 -o FliI_multimer_site4_200618 -p HIC6_model -n 431 -c 443 -r 3 -u 20 -m 100 -t 2.5 -v 3 -x A,D,G,J,M,P -l > FliI_multimer_site4_20core_o3_200618.log
