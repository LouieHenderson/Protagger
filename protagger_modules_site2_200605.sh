#!/bin/sh
#SBATCH --ntasks=50

./bin/runProtagger.py -b ./database_200606 -o 5jxl_multimer_test_site2_200606 -p 5jxl_multimer -n 543 -c 554 -r 3 -u 50 -t 2.5 -l > sbatchmodules_50core_200606_site2_o3.log
