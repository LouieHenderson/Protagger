#!/bin/sh
#SBATCH --ntasks=10

./bin/runProtagger.py -b ./thermdata_200707 -o Campy_FliM_site3_therm_200709 -p Campy_FliGM_cring -n 82 -c 88 -r 3 -u 10 -m 80 -t 2.5 -v 3 -x W,X,Y,Z,a -l > Campy_FliM_site3_10core_o3_therm_200709.log
