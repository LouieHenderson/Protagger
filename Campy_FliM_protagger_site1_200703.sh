#!/bin/sh
#SBATCH --ntasks=10

./bin/runProtagger.py -b ./database_200625 -o Campy_FliM_site1_200706 -p Campy_FliGM_cring -n 201 -c 218 -r 3 -u 10 -m 80 -t 2.5 -v 3 -x W,X,Y,Z,a -l -d > Campy_FliM_site1_10core_o3_200706.log
