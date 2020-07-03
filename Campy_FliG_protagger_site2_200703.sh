#!/bin/sh
#SBATCH --ntasks=10

./bin/runProtagger.py -b ./FliG_1tag -o Campy_FliG_site2_200703 -p Campy_FliGM_cring.pdb -n 221 -c 229 -r 4 -u 10 -m 80 -t 2.5 -v 3 -x I,J,K,L,M -l > Campy_FliG_site2_10core_o3_200703.log

mv Campy_FliG_site2_10core_o3_200703.log ./Campy_FliG_site2_200703
