#!/bin/sh
#SBATCH --ntasks=50 --partition=heartofgoldCPU

../bin/runProtagger.py -b ../1m40_db/ -o YYYY_1m40_200918 -p YYYY -n 205 -c 212 -r 3 -u 50 -t 2.5 -v 4 -l > protagger_YYYY_1m40_db_200918.log

