#!/usr/bin/python
#import Bio.PDB
#from Bio import PDB
from argparse import ArgumentParser
import sys, os
import time

parser = ArgumentParser(description='Rename protein chains')

#Adds parser arguments which the user will input
parser.add_argument('-p', '--pdb', type=str, help='Enter pdb to rename chains (lower case only)')
parser.add_argument('-o', '--out', type=str, help='Output name of altered pdb')
parser.add_argument('-c', '--chain', type=str, help='Enter chain identifier to start new chain ID')

#Associates arguments to objects
args = parser.parse_args()
pdb_input    = args.pdb
pdb_output   = str(args.out)
chain_id     = args.chain

# pymol launching: quiet (-q), without GUI (-c) and with arguments from command line
import pymol
pymol.pymol_argv = ['pymol','-qc',pdb_input + '.cif']
pymol.finish_launching()
cmd = pymol.cmd
time.sleep(1)

chainlist = []
for id in cmd.get_chains():
    chainlist.append(id)

#print "Here's the initial list", chainlist

# Create alphabet list of lowercase letters
alphabet        = []

for letter in range(65, 91):
    alphabet.append(chain_id + chr(letter))

for letter in range(97,123):
    alphabet.append(chain_id + chr(letter))

#print "Before", cmd.get_chains()

#for chain in chainlist:
#    print chain

#Rename chain ids
for chain_init, new_id in zip(chainlist, alphabet):
    if chain_init.isalpha() == True or chain_init.isdigit() == True:
        #print chain, new_id
        #print "chain "+ chain, "chain=\""+new_id+"\""
        cmd.alter("chain "+ chain_init, "chain=\""+new_id+"\"")
        cmd.sort("chain "+new_id)
    else:
        #print chain, new_id
        #print "chain " + "\""+chain+"\"", "chain=\""+new_id+"\""
        cmd.alter("chain " + "\""+chain_init+"\"", "chain=\""+new_id+"\"")
        cmd.alter("chain " + "\\"+chain_init, "chain=\""+new_id+"\"")
        cmd.sort("chain "+new_id)

cmd.sort(pdb_input)

print "After", cmd.get_chains()

cmd.save(pdb_output + ".cif", pdb_input)

print "Complete!"
