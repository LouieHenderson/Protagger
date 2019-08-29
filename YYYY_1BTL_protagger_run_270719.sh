#!/bin/sh
#SBATCH --ntasks=50

./bin/runProtagger.py -b ../../1BTL_ecoli_tag/ -o YYYY_1BTL_db_190_283_270819 -p YYYY -n 190 -c 283 -r 3 -u 50 -t 2.5 -v 3 -e 1.5 -l > protagger_YYYY_1BTL_270819.log

