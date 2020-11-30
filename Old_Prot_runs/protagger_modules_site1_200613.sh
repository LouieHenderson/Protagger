#!/bin/sh
#SBATCH --ntasks=50

./bin/runProtagger.py -b ./database_200606 -o 5jxl_multimer_test_site3_200607 -p 5jxl_multimer -n 542 -c 554 -r 3 -u 50 -t 2.5 -v 4 -l > sbatchmodules_50core_200613_site1_o4.log
