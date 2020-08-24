#!/bin/sh
#SBATCH --ntasks=15

./bin/runProtagger.py -b ./1tqg_db -o Campy_FliGn_site2_1tqg_240920 -p Campy_FliMGnc_cring -n 37 -c 43 -r 3 -u 4 -m 80 -t 2.5 -v 3 -x I,J,K,L,M -l > Campy_FliGn_site2_15core_o3_1tqg_240820.log
