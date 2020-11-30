#!/bin/sh
#SBATCH --ntasks=20

./bin/runProtagger.py -b ./thermdata_200707 -o Campy_FliI_site2_therm_200919 -p HIC6_model -n 378 -c 383 -r 3 -u 20 -m 80 -t 2.5 -v 3 -x A,D,G,J,M,P -l > Campy_FliI_site2_20core_o3_therm_200819.log
