#!/bin/sh
#SBATCH --ntasks=20

./bin/runProtagger.py -b ./database_200606 -o FliH_multimer_site1_200622 -p HIC6_model -n 138 -c 149 -r 3 -u 20 -m 100 -t 2.5 -v 3 -x B,C,E,F,H,I,K,L,N,O,Q,R -l > FliH_multimer_site1_20core_o3_200622.log
