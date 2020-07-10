#!/bin/sh
#SBATCH --ntasks=10

./bin/runProtagger.py -b ./thermdata_200707 -o Campy_FliG_site1_therm_200707 -p Campy_FliGM_cring -n 263 -c 269 -r 3 -u 10 -m 80 -t 2.5 -v 3 -x I,J,K,L,M -l > Campy_FliG_site1_therm_10core_o3_200707.log
