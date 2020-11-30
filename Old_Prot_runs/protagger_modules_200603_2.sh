#!/bin/sh
#SBATCH --ntasks=50

./bin/runProtagger.py -b ../../ecoli_tags_small_220218 -o 5jxl_multimer_test_200604 -p 5jxl_multimer -n 264 -c 282 -r 3 -u 50 -t 2.5 -l > sbatchmodules_50core_200604_o3.log
