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
parser.add_argument('-c1', '--cif1', type=str, help='Enter cif 1')
parser.add_argument('-c2', '--cif2', type=str, help='Enter cif 2')
parser.add_argument('-c3', '--cif3', type=str, help='Enter cif 3')
parser.add_argument('-o', '--out', type=str, help='Output name of combined cif')

#Associates arguments to objects
args            = parser.parse_args()
cif1_input      = args.cif1
cif2_input      = args.cif2
cif3_input      = args.cif3
cif_output      = str(args.out)

parser = PDB.MMCIFParser(QUIET = True)
structure1 = parser.get_structure(cif1_input, cif1_input)
structure2 = parser.get_structure(cif2_input, cif2_input)
structure3 = parser.get_structure(cif3_input, cif3_input)

scrubbed_cif1 = HETscrubber(structure1)
scrubbed_cif2 = HETscrubber(structure2)
scrubbed_cif3 = HETscrubber(structure3)

print "First cif"
for x in scrubbed_cif1:
    print "yoy", x
    for y in x:
        print "yay", y

print "Second cif"
for x in scrubbed_cif2:
    print "yoy", x
    for y in x:
        print "yay", y

print "Third cif"
for x in scrubbed_cif3:
    print "yoy", x
    for y in x:
        print "yay", y
    #    for z in y:
    #        print "yey", z

print scrubbed_cif1

#Saves the PDB of the optimal final acceptor
io=Bio.PDB.PDBIO()
io.set_structure(structure1)
io.save(cif_output + ".cif", structure1)

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
