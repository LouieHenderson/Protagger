#!/bin/sh
#SBATCH --ntasks=50

./bin/runProtagger.py -d ../ecoli_tags_small_220218 -p YYYY -n 191 -c 283 -o 3 -r 3 -u 50 -l > sbatchmodules_50core_130319_o3.log
