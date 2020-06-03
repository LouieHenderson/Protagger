#!/bin/sh
#SBATCH --ntasks=50

./bin/runProtagger.py -b ../../tagdatabase_2fwg -o 1yab_2fwg -p 1yab -n 107 -c 146 -r 3 -u 4 -t 2.5 -v 3 -e 1.5 -l > test_1yab_290819.log
