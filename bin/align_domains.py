#!/usr/bin/env python
import re
import copy
import argparse
from Bio import PDB
import Bio.PDB
from Bio.PDB import PDBList
from Bio.PDB import PDBParser
from Bio.PDB import StructureBuilder
from Bio.PDB import Selection
from Bio.PDB import Entity
from Bio.PDB import PDBIO

aminoAcids = {'ALA':1, 'ARG':1, 'ASN':1, 'ASP':1, 'CYS':1, 'GLN':1, 'GLU':1, 'GLY':1, 'HIS':1, 'ILE':1, 'LEU':1, 'LYS':1, 'MET':1, 'PHE':1, 'PRO':1, 'SER':1, 'THR':1, 'TRP':1, 'TYR':1, 'VAL':1}

#Iterates through list of amino acids, if heteroatom is found, residue is removed
def HETscrubber(PDBstructure):

    for chain in PDBstructure[0]:
        for residue in chain.get_list():
            resid = residue.id
            if resid[0] != ' ' or residue.__contains__('CA') == False or residue.get_resname() in aminoAcids == False:
                chain.detach_child(resid)

    return PDBstructure

def align_regions(FixedPDB, MovingPDB, FixedResListN, FixedResListC, MovingResListN, MovingResListC, chainregx):

    parser = PDB.PDBParser(PERMISSIVE=1)
    FixedInit = parser.get_structure(FixedPDB,"./"+FixedPDB+".pdb")

    parser = PDB.PDBParser(PERMISSIVE=1)
    MovingInit = parser.get_structure(MovingPDB,"./"+MovingPDB+".pdb")

    Fixed   = HETscrubber(FixedInit)
    Moving  = HETscrubber(MovingInit)

    FixedCA = []
    MovingCA = []

    #Removes residues not in range of region to align
    for model in Fixed:
        for chain in model:
            if (chain.id == chainregx[0]) == True:
                for residue in chain.get_list():
                    if residue.get_id()[1] in range(FixedResListN, FixedResListC + 1):
                        FixedCA.append(residue['CA'])

    for model in Moving:
        for chain in model:
            chain.id = chainregx[0]
            for residue in chain.get_list():
                if residue.get_id()[1] in range(MovingResListN, MovingResListC + 1):
                    MovingCA.append(residue['CA'])

    super_imposer = Bio.PDB.Superimposer()
    super_imposer.set_atoms(FixedCA, MovingCA)
    super_imposer.apply(Moving.get_atoms())

    return Moving

def main():
    #Defines arguments to be included when running the program
    parser= argparse.ArgumentParser(description='Align regions of interest')

    #Adds parser arguments which the user will input
    parser.add_argument('-a', '--fixedpdb', type=str, help='Enter the pdb file to use as fixed list', required = True)
    parser.add_argument('-b', '--movepdb', type=str, help='Enter the pdb file to use as moving list', required = True)
    parser.add_argument('-c', '--fixedlistN', type=int, help='Enter the N residue to align in fixed list', required = True)
    parser.add_argument('-d', '--fixedlistC', type=int, help='Enter the C residue to align in fixed list', required = True)
    parser.add_argument('-e', '--movinglistN', type=int, help='Enter the N residue to align in moving list', required = True)
    parser.add_argument('-f', '--movinglistC', type=int, help='Enter the C residue to align in moving list', required = True)
    parser.add_argument('-x', '--chain', type=str, help='Enter chain identifiers of those you wish to tag. E.g. A-C,F,X')

    args = parser.parse_args()

    FixedPDB        =   args.fixedpdb
    MovingPDB       =   args.movepdb
    FixedListN      =   args.fixedlistN
    FixedListC      =   args.fixedlistC
    MovingListN     =   args.movinglistN
    MovingListC     =   args.movinglistC
    chainlabel      =   args.chain

    #Creates the list containing parameters for the regular expression for chain IDs
    chainregx   = []

    if chainlabel != None and ((',' in chainlabel) == True):
        chainregx = str(chainlabel).split(",")
    elif chainlabel != None and ((',' in chainlabel) != True):
        chainregx.append(str(chainlabel))
    else:
        chainregx.append('\D')

    OptMoving = align_regions(FixedPDB, MovingPDB, FixedListN, FixedListC, MovingListN, MovingListC, chainregx)

    io=Bio.PDB.PDBIO()
    io.set_structure(OptMoving)
    io.save("./"+MovingPDB+'_aligned_to_'+FixedPDB+'_chain_'+chainregx[0]+'.pdb')

    print "Regions successfully aligned!"

main()
