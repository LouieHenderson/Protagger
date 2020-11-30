#!/bin/sh
#SBATCH --ntasks=50

./bin/runProtagger.py -b ../../ecoli_tags_small_220218 -o 5jxl_multimer_test_200603 -p 5jxl_multimer -n 543 -c 554 -r 3 -u 50 -t 2.5 -l > sbatchmodules_50core_200603_o3.log
