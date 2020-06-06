#!/bin/sh
#SBATCH --ntasks=50

./bin/runProtagger.py -b ./database_200605 -o 5jxl_multimer_test_site2_200605 -p 5jxl_multimer -n 543 -c 554 -r 3 -u 50 -t 2.5 -l > sbatchmodules_50core_200605_site2_o3.log
