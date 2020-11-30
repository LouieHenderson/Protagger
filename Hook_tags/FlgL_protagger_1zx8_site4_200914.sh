#!/bin/sh
#SBATCH --ntasks=15 --partition=heartofgoldCPU

./bin/runProtagger.py -b ./1zx8_db -o FlgL_multimer_1zx8_site4_200914 -p FlgEKL_partial_multimer -n 569 -c 590 -r 3 -u 15 -m 80 -t 2.5 -v 3 -x L,M,N,O,P,Q -l > FlgL_multimer_site4_1zx8_15core_o3_200914.log
