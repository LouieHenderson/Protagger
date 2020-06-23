#!/usr/bin/python
import Bio.PDB
from Bio import PDB
from argparse import ArgumentParser
import sys, os
import time

aminoAcids = {'ALA':1, 'ARG':1, 'ASN':1, 'ASP':1, 'CYS':1, 'GLN':1, 'GLU':1, 'GLY':1, 'HIS':1, 'ILE':1, 'LEU':1, 'LYS':1, 'MET':1, 'PHE':1, 'PRO':1, 'SER':1, 'THR':1, 'TRP':1, 'TYR':1, 'VAL':1}

def HETscrubber(PDBstructure):
    for chain in PDBstructure[0]:
        for residue in chain.get_list():
            resid = residue.id
            if resid[0] != ' ' or residue.__contains__('CA') == False or residue.get_resname() in aminoAcids == False:
                chain.detach_child(resid)

    return PDBstructure

parser = ArgumentParser(description='Rename protein chains')

#Adds parser arguments which the user will input
parser.add_argument('-p', '--pdb', type=str, help='Enter pdb to rename chains (lower case only)')
parser.add_argument('-o', '--out', type=str, help='Output name of altered pdb')

#Associates arguments to objects
args         = parser.parse_args()
pdb_input    = args.pdb
pdb_output   = str(args.out)

parser = PDB.MMCIFParser(QUIET = True)
structure = parser.get_structure(pdb_input, pdb_input)

scrubbed_cif = HETscrubber(structure)

for x in scrubbed_cif:
    print "yoy", x
    for y in x:
        print "yay", y
    #    for z in y:
    #        print "yey", z

print scrubbed_cif

#Saves the PDB of the optimal final acceptor
io=Bio.PDB.PDBIO()
io.set_structure(scrubbed_cif)
io.save(pdb_output + ".cif", scrubbed_cif)

print "Complete!"

# pymol launching: quiet (-q), without GUI (-c) and with arguments from command line
#import pymol
#pymol.pymol_argv = ['pymol','-qc',pdb_input + '.cif']
#pymol.finish_launching()
#cmd = pymol.cmd
#time.sleep(1)

#chainlist = []
#for id in cmd.get_chains():#
#    chainlist.append(id)

#print "Here's the initial list", chainlist

# Create alphabet list of lowercase letters
#alphabet        = []

#for letter in range(65, 91):
#    alphabet.append(chain_id + chr(letter))

#for letter in range(97,123):
#    alphabet.append(chain_id + chr(letter))

#print "Before", cmd.get_chains()

#for chain in chainlist:
#    print chain

#Rename chain ids
#for chain_init, new_id in zip(chainlist, alphabet):
#    if chain_init.isalpha() == True or chain_init.isdigit() == True:
#        #print chain, new_id
#        #print "chain "+ chain, "chain=\""+new_id+"\""
#        cmd.alter("chain "+ chain_init, "chain=\""+new_id+"\"")
#        cmd.sort("chain "+new_id)
#    else:
#        #print chain, new_id
#        #print "chain " + "\""+chain+"\"", "chain=\""+new_id+"\""
#        cmd.alter("chain " + "\""+chain_init+"\"", "chain=\""+new_id+"\"")
#        cmd.alter("chain " + "\\"+chain_init, "chain=\""+new_id+"\"")
#        cmd.sort("chain "+new_id)

#cmd.sort(pdb_input)

#print "After", cmd.get_chains()

#cmd.save(pdb_output + ".cif", pdb_input, 0)
