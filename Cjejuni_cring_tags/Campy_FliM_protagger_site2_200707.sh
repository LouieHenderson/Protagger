#!/bin/sh
#SBATCH --ntasks=10

./bin/runProtagger.py -b ./thermdata_200707 -o Campy_FliM_site2_therm_200707 -p Campy_FliGM_cring -n 119 -c 124 -r 3 -u 10 -m 80 -t 2.5 -v 3 -x W,X,Y,Z,a -l > Campy_FliM_site2_10core_o3_therm_200707.log
