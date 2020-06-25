#!/usr/bin/python
#import Bio.PDB
#from Bio import PDB
from argparse import ArgumentParser
import sys, os
import time

parser = ArgumentParser(description='Rename protein chains')

#Adds parser arguments which the user will input
parser.add_argument('-p1', '--pdb1', type=str, help='Enter pdb 1')
parser.add_argument('-p2', '--pdb2', type=str, help='Enter pdb 2')
parser.add_argument('-p3', '--pdb3', type=str, help='Enter pdb 3')
parser.add_argument('-o', '--out', type=str, help='Output name of altered pdb')

#Associates arguments to objects
args = parser.parse_args()
pdb_input1    = args.pdb1
pdb_input2    = args.pdb2
pdb_input3    = args.pdb3
pdb_output   = str(args.out)


# pymol launching: quiet (-q), without GUI (-c) and with arguments from command line
import pymol
pymol.pymol_argv = ['pymol','-qc',pdb_input1 + '.pdb']#,pdb_input2 + '.pdb',pdb_input3 + '.pdb',]
pymol.finish_launching()
cmd = pymol.cmd
time.sleep(1)

chainlist1 = []
for id in cmd.get_chains(pdb_input1):
    chainlist1.append(id)

print "Here's the initial list", cmd.get_chains(pdb_input1)

# Create alphabet list of lowercase letters
alphabet        = []

for letter in range(65, 91):
    alphabet.append(chr(letter))

for letter in range(97,123):
    alphabet.append(chr(letter))

print "Before 1", cmd.get_chains(pdb_input1)
#for chain in chainlist:
#    print chain

#Rename chain ids
for chain_init, new_id in zip(cmd.get_chains(pdb_input1), alphabet):
    if chain_init.isalpha() == True or chain_init.isdigit() == True:
        #print chain, new_id
        #print "chain "+ chain, "chain=\""+new_id+"\""
        cmd.alter("chain "+ chain_init, "chain=\""+new_id+"\"")
        cmd.sort("chain "+new_id)
#        print new_id
        alphabet.remove(new_id)
    else:
        #print chain, new_id
        #print "chain " + "\""+chain+"\"", "chain=\""+new_id+"\""
        cmd.alter("chain " + "\""+chain_init+"\"", "chain=\""+new_id+"\"")
        cmd.alter("chain " + "\\"+chain_init, "chain=\""+new_id+"\"")
        cmd.sort("chain "+new_id)
#        print new_id
        alphabet.remove(new_id)

#print "yoyoyoy", alphabet

chainlist1 = []
cmd.sort(pdb_input1)

print "After pdb 1", cmd.get_chains(pdb_input1)

cmd.load(pdb_input2+".pdb")
time.sleep(1)

chainlist2 = []
for id in cmd.get_chains(pdb_input2):
    chainlist2.append(id)

print "Before 2", cmd.get_chains(pdb_input2)

for chain_init, new_id in zip(cmd.get_chains(pdb_input2), alphabet):
    if chain_init.isalpha() == True or chain_init.isdigit() == True:
        #print chain, new_id
        #print "chain "+ chain, "chain=\""+new_id+"\""
        cmd.alter("chain "+ chain_init, "chain=\""+new_id+"\"")
        cmd.sort("chain "+new_id)
        alphabet.remove(new_id)
    else:
        #print chain, new_id
        #print "chain " + "\""+chain+"\"", "chain=\""+new_id+"\""
        cmd.alter("chain " + "\""+chain_init+"\"", "chain=\""+new_id+"\"")
        cmd.alter("chain " + "\\"+chain_init, "chain=\""+new_id+"\"")
        cmd.sort("chain "+new_id)
        alphabet.remove(new_id)

#print "yoyoyoy", alphabet

chainlist2 = []
cmd.sort(pdb_input2)

print "After pdb 2", cmd.get_chains(pdb_input2)

cmd.load(pdb_input3+".pdb")
time.sleep(1)

chainlist3 = []
for id in cmd.get_chains(pdb_input3):
    chainlist3.append(id)

print "Before 3", cmd.get_chains(pdb_input3)

for chain_init, new_id in zip(cmd.get_chains(pdb_input3), alphabet):
    if chain_init.isalpha() == True or chain_init.isdigit() == True:
        #print chain, new_id
        #print "chain "+ chain, "chain=\""+new_id+"\""
        cmd.alter("chain "+ chain_init, "chain=\""+new_id+"\"")
        cmd.sort("chain "+new_id)
        alphabet.remove(new_id)
    else:
        #print chain, new_id
        #print "chain " + "\""+chain+"\"", "chain=\""+new_id+"\""
        cmd.alter("chain " + "\""+chain_init+"\"", "chain=\""+new_id+"\"")
        cmd.alter("chain " + "\\"+chain_init, "chain=\""+new_id+"\"")
        cmd.sort("chain "+new_id)
        alphabet.remove(new_id)

#print "yoyoyoy", alphabet

cmd.sort(pdb_input3)

print "After pdb 3", cmd.get_chains(pdb_input3)

cmd.save(pdb_output + ".pdb")

print "Complete!"
