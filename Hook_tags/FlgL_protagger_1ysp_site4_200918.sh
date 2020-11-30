#!/bin/sh
#SBATCH --ntasks=15 --partition=heartofgoldCPU

./bin/runProtagger.py -b ./1ysp_db -o FlgLonly_multimer_1ysp_site4_200918 -p FlgL_Lstart -n 569 -c 590 -r 3 -u 15 -m 80 -t 2.5 -v 3 -x L,M,N,O,P,Q -l > FlgLonly_multimer_site4_1ysp_15core_o3_200918.log
