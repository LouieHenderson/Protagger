#!/bin/sh
#SBATCH --ntasks=15

./bin/runProtagger.py -b ./thermdata_no3idu_200707 -o Campy_FliGn_site4_no3idu_200924 -p Campy_FliMGnc_cring -n 76 -c 86 -r 3 -u 15 -m 80 -t 2.5 -v 3 -x I,J,K,L,M -l > Campy_FliGn_site4_15core_o3_no3idu_200824.log
