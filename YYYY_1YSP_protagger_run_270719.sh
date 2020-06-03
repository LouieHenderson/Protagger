#!/bin/sh
#SBATCH --ntasks=50

./bin/runProtagger.py -b ../../1YSP_only/ -o YYYY_1YSP_db_200_220_spec_270819 -p YYYY -n 200 -c 220 -r 3 -u 50 -t 2.5 -v 3 -e 1.5 -l -d > protagger_YYYY_1YSP_spec_280819.log

