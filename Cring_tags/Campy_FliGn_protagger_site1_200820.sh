#!/bin/sh
#SBATCH --ntasks=15

./bin/runProtagger.py -b ./thermdata_200707 -o Campy_FliGn_site1_therm_200920 -p Campy_FliMGnc_cring -n 26 -c 32 -r 3 -u 15 -m 80 -t 2.5 -v 3 -x I,J,K,L,M -l > Campy_FliGn_site1_15core_o3_therm_200820.log
