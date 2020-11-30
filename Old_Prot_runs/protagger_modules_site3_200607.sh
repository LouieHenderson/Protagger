#!/bin/sh
#SBATCH --ntasks=50

./bin/runProtagger.py -b ./database_200606 -o 5jxl_multimer_test_site3_200607 -p 5jxl_multimer -n 428 -c 460 -r 3 -u 50 -t 2.5 -l > sbatchmodules_50core_200607_site3_o3.log
