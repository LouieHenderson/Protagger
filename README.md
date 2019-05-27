# Protagger
Primary repository for Protagger scripts and updates

```
         **********
       **         **
      **          **
     **         **
    ************
   **
  ** ******  ***
 ** **     **  **
** **       ***

This is the main repository
```

Test code:

```
./bin/runProtagger.py -d ./tagdatabase_1tag -p 1yab -n 107 -c 146 -o 3 -r 3 -u 3 -t 3 -l -s

./bin/runProtagger.py -d ./tagdatabase_1tag/ -p YYYY -n 191 -c 283 -o 3 -r 2.34 -u 2 -l

Parameters
~~~~~~~~~~

-d/--tagdatabase: File directory location of tag database to be checked
-p/--pdb: Name of pdb file to be tagged (not including file extension)
-n/--nterm: Enter N-terminal residue position number of tag site on acceptor
-c/--cterm: Enter C-terminal residue position number of tag site on acceptor
-o/--overlap: Enter overlap length used for N/C terminal alignment
-r/--RMSDthreshold: Enter RMSD threshold under which to check. Recommended value = 3, larger values increase runtime
-u/--CPU: Enter the number of CPUs to use during parallel processing steps. Multiple cores recommended for greater numbers of chains to tag.
-t/--atomthreshold: Enter the threshold (in angstrom) of detection for clashes between neighbouring atoms. Default/minimum is 2.5, user higher values for stricter clash check (e.g. 4)
-l/--local: Acceptor pdb is stored locally
-x/--chain: Enter chain identifiers of those you wish to tag. E.g. A-C,F,X
-s/--specific: Find the ideal tagsite for a specific protein within range N/C terminals. Use larger distances between N/C insert site for better results. Recommended to use alongside multiple cores OR fewer potential tags.
```
