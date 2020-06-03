#!/bin/sh
#SBATCH --ntasks=50

./bin/runProtagger.py -b ../../tagdatabase_1tag/ -o YYYY_1610_db_205_210_270819 -p YYYY -n 205 -c 210 -r 3 -u 50 -t 2.5 -v 4 -m 150 -e 1.5 -l > protagger_YYYY_1610_db_1m40_270819.log

