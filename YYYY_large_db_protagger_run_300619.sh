#!/bin/sh
#SBATCH --ntasks=50

./bin/runProtagger.py -b ../../tagdatabase_250619/ -o YYYY_1610_db_190819 -p YYYY -n 209 -c 243 -r 3 -u 50 -t 2.5 -v 4 -l > protagger_YYYY_1610_db_209_243_190819.log

