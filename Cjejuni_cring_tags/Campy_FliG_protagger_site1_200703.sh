#!/bin/sh
#SBATCH --ntasks=10 --nodelist=heartofgold2

./bin/runProtagger.py -b ./FliG_1tag -o Campy_FliG_site1_200703 -p Campy_FliGM_cring -n 251 -c 279 -r 3 -u 10 -m 80 -t 2.5 -v 3 -x I,J,K,L,M -l -d > Campy_FliG_site1_10core_o3_200703.log
