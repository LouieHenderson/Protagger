#!/bin/sh
#SBATCH --ntasks=50

./bin/runProtagger.py -b ../ecoli_tags_small_220218 -o Campy_FliG_Ndomain_tester_280519 -p Campy_FliG -n 13 -c 99 -r 3 -u 50 -t 2.5 -l -d > protagger_sbatch_Campy_FliG_Ndomain_tester_280519.log
