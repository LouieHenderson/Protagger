#!/bin/sh
#SBATCH --ntasks=50

./bin/runProtagger.py -b ./database_200606 -o 5jxl_multimer_test_site1_200614 -p 5jxl_multimer -n 542 -c 554 -r 3 -u 50 -t 2.5 -v 2 -l > sbatchmodules_50core_200614_site1_o2.log
