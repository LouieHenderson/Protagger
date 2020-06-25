#!/bin/sh

set -e

# Test 1: the NOMERGE option
# Merge original rnase with a file of symmetry mates

#pdbcur xyzin $CEXAM/rnase/rnase.pdb xyzout $CCP4_SCR/rnase_symop.pdb <<EOF
#symop -X+1/2,-Y,Z+1/2
#symcommit
#EOF

pdb_merge xyzin1 ./total.cif \
          xyzout ./tester.cif <<eof
nomerge
end
eof

