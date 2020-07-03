#!/bin/sh
#SBATCH --ntasks=10

./bin/runProtagger.py -b ./FliM_1tag -o FliM_multimer_site3_200703 -p Campy_FliGM_cring -n 86 -c 91 -r 4 -u 10 -m 80 -t 2.5 -v 3 -x W,X,Y,Z,a -l > Campy_FliM_site3_10core_o3_200703.log

mv Campy_FliM_site3_10core_o3_200703.log ./FliM_multimer_site3_200703
