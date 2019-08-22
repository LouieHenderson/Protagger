#!/bin/sh
#SBATCH --ntasks=50

./bin/runProtagger.py -b ../../tagdatabase_250619/ -o YYYY_1610_db_specific_200_230_210819 -p YYYY -n 200 -c 230 -r 3 -u 50 -t 2.5 -v 4 -m 150 -e 1.5 -d -l > protagger_YYYY_1610_db_spec_200_230_210819.log

