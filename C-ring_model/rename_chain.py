#!/usr/bin/python
import Bio.PDB
from Bio import PDB
from argparse import ArgumentParser

parser = ArgumentParser(description='Rename protein chains')

#Adds parser arguments which the user will input
parser.add_argument('-p', '--pdb', type=str, help='Enter pdb to rename chains (lower case only)', required = True)
parser.add_argument('-o', '--out', type=str, help='Output name of altered pdb', required = True)
parser.add_argument('-c', '--chain', type=str, help='Enter chain identifier to start new chain ID')

#Associates arguments to objects
args = parser.parse_args()
pdb_input    = args.pdb
pdb_output   = str(args.out)
chain_id     = args.chain

#Read pdb file into structure object
def read_pdb(pdb_input):

    #Creates structure object of protein to be tagged
    parser = PDB.PDBParser(PERMISSIVE=1)
    pdb_change = parser.get_structure(pdb_input,"./"+pdb_input+".pdb")

    return pdb_change

#Run program
Init_PDB = read_pdb(pdb_input)

id_len = len(Init_PDB.get_list()[0].get_list())

# Create alphabet list of lowercase letters
alphabet = []

for letter in range(65, 91):
    alphabet.append(chain_id + chr(letter))

for letter in range(97,123):
    alphabet.append(chain_id + chr(letter))

#for model in Init_PDB:
#    print model
#    for chain in model:
#        print chain
#        for residue in chain:
#            #print residue.id[1]
#            if (residue.id[1] == 5) == True:
#                print residue

#Rename chain ids
for chain, new_id in zip(Init_PDB.get_chains(), alphabet):
    print type(new_id)
    print type(str(new_id))
    chain.id = new_id

#for model in Init_PDB:
#    print model
#    for chain in model:
#        print chain
#        for residue in chain:
#            #print residue.id[1]
#            if (residue.id[1] == 5) == True:
#                print residue


#print pdb_output
print "Complete!"

#Saves the PDB of the optimal final acceptor
io=Bio.PDB.PDBIO()
io.set_structure(Init_PDB)
io.save(".pdb")
