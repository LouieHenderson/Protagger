#!/usr/bin/python

import numpy
import sys
import os
import argparse
import re
import glob
import itertools
import copy
import time
import multiprocessing as mp
from multiprocessing import Pool, Pipe, Queue, Process
import cPickle as pickle
from argparse import ArgumentParser
from itertools import chain
from functools import partial
from Bio import PDB
import Bio.PDB
from Bio.PDB import PDBList
from Bio.PDB import PDBParser
from Bio.PDB import StructureBuilder
from Bio.PDB import Selection
from Bio.PDB import Entity
from Bio.PDB import PDBIO
from CDROM import DVD_CGMS_UNRESTRICTED

aminoAcids = {'ALA':1, 'ARG':1, 'ASN':1, 'ASP':1, 'CYS':1, 'GLN':1, 'GLU':1, 'GLY':1, 'HIS':1, 'ILE':1, 'LEU':1, 'LYS':1, 'MET':1, 'PHE':1, 'PRO':1, 'SER':1, 'THR':1, 'TRP':1, 'TYR':1, 'VAL':1}

#Global objects used in optimal RMSD determination and tag structure formation
Checkpoint     = 0
Startime       = 0
Optimal        = []
AccepTup       = []
AcceptorIter   = None
DonorIter      = None
sphereclash1   = []
donorchains    = []
TrueCheck      = [True, True]
threshold_A    = 0
CPUs           = 0
chainregx      = []
Queueregulator = False

def downloadpdb(pdb_acceptor):
    # Downloads pdb file to local folder
    folder = pdb_acceptor
    pdbl = PDBList()
    data_set = pdbl.retrieve_pdb_file(pdb_acceptor, pdir=pdb_acceptor)
  
    #Parses PDB file to return a structure object
    parser = PDB.PDBParser(PERMISSIVE=1)

    #Creates structure object of protein to be tagged
    acceptor = parser.get_structure(pdb_acceptor,folder+'/'+pdb_acceptor+".pdb")
   
    return acceptor

def localpdb(pdb_acceptor):    
    #Parses PDB file to return a structure object
    parser = PDB.PDBParser(PERMISSIVE=1)

    #Creates structure object of protein to be tagged
    acceptor = parser.get_structure(pdb_acceptor,"./"+pdb_acceptor+".pdb")

    return acceptor

def findClash(listofchains, start_tagged_residue, end_tagged_residue, spherelist, start_potentialTag, end_potentialTag, overlap, sphereclash):
   #Provides the threshold at which atomic distances for clashes are rejected
    global threshold_A #= 3

    acceptorchains = []
    donorchains = []
    spherex = []

    #Separates acceptor chains from donor chains into respective objects
    for chain in listofchains.get_list()[0].get_list():
        if (re.findall('(?<=.)tag', str(chain.get_id()))) == ['tag']:
            donorchains.append(chain)
        else:
            acceptorchains.append(chain)    
     
    #Provides a checkpoint object which acts as gatekeeper, if any clash is found, all subsequent checks abandoned
    #Compares every atom within every residue of the N/C tag sites with all atoms in the tags minus NC truncation
    #Excludes clash check of region between N/C terminal of Acceptor that will be excised
    global TrueCheck
    TrueCheck = [True, True] 
    
    #Iterates through acceptor/donor pairs, provides rough spherical approximations which indicate overlap to identify which donors to clash check in finer grain
    for chain in acceptorchains:#listofchains.get_list()[0].get_list():#acceptorchains:
        #print 'chain', chain
        if TrueCheck[0] == True and TrueCheck[1] == True:
            spherex = []
            checkspheres = []
            checkdonors = []
            spectag = [] 
            
            for spheres in sphereclash:
                if spheres[0][2] == chain.get_id():
                    spherex.append(spheres)
            for spheres in spherex:
                if spheres[2] == False or spheres[1][2] == (str(chain.get_id())+'tag'):
                    checkspheres.append(spheres)
            for donors in donorchains:
                for sphere in checkspheres:
                    if donors.get_id() == sphere[1][2]:
                        checkdonors.append(donors)              
            
            #Chunk object created to start multiprocessing child processes
            reschunks = list(chunks(chain.get_list(), (len(chain.get_list()) / CPUs)+1))
            
            #print 'reschunks', reschunks
            #Multiprocessing block which generates fine grain clash check processes.
            #Uses TrueCheck object booleans to track the identification of a False return of clash check in any of the child processes 
            if __name__ == '__main__': 
                queue1 = Queue()
                for chunk in reschunks:                                   
                    process = Process(target=restruecheck, args=(start_tagged_residue, end_tagged_residue, start_potentialTag, end_potentialTag, overlap, threshold_A, checkdonors, queue1, chunk))
                    process.Daemon = True
                    process.start()
                for chunk in reschunks:
                    queue1.put(TrueCheck)
                    while len(mp.active_children()) > 1 and TrueCheck[0] != False and TrueCheck[1] != False:          
                        TrueCheck = queue1.get()
                        if TrueCheck[0] == False or TrueCheck[1] == False:
                            process.join()
                            break
                        queue1.put(TrueCheck)                 
            
            #Control steps to prevent runaway process generation and to ensure correct process/queue closure  
            while len(mp.active_children()) > 1:
                time.sleep(0.1)
                  
            process.join() 
            time.sleep(0.1)
            queue1.close()
            
            
    #print 'Final TrueCheck being returned = ', TrueCheck
        
    return TrueCheck
       
def chunks(list, n):
    for i in range(0, len(list), n):
        yield list[i:i+n]  

def restruecheck(start_tagged_residue, end_tagged_residue, start_potentialTag, end_potentialTag, overlap, threshold_A, checkdonors, queue1, residlist):     
     global TrueCheck
     
     #Iterate through each atom in residues of compared acceptor/donor pairs and checks for potential clash
     for resid1 in residlist:
         Checker = queue1.get()
         queue1.put(Checker)
         if Checker[0] == False or Checker[1] == False:
             TrueCheck = Checker
             queue1.put(TrueCheck)
         CheckpointA = resid1.get_id()[1] <= (start_tagged_residue - 1 - overlap) or resid1.get_id()[1] >= (end_tagged_residue + 1 + overlap)
         if CheckpointA == True: 
             if TrueCheck[0] == True:
                 for atom1 in resid1:
                     for chain2 in checkdonors:                                                   
                         for resid2 in chain2:
                             CheckpointB = resid2.get_id()[1] >= (start_potentialTag + 1 + overlap) and resid2.get_id()[1] <= (end_potentialTag - 1 - overlap)
                             if CheckpointB == True and TrueCheck[0] == True:
                                 for atom2 in resid2:
                                     Atomdist = atom1 - atom2
                                     Boolean1 = Atomdist >= threshold_A
                                     if Boolean1 == True:
                                         TrueCheck[0] = True
                                     elif Boolean1 == False:
                                         TrueCheck[0] = False
                                         queue1.put(TrueCheck) 
                                         break                           
                                     else:
                                         pass 
                             else:
                                 pass                           
             elif TrueCheck[0] == False:
                 queue1.put(TrueCheck)
                 break
         else:
             pass         

def parsetag(donor):   
    #Regular expression search for .ent files in Tag folder, produces file path
    #Parses PDB file to return a structure object
    parser  = PDB.PDBParser(PERMISSIVE=1, QUIET=True)
    tags    = []
    tags2   = []
    tags3   = []
    tags4   = []
    
    #Iterates through list of tags and returns structure objects
    #Turns file path list objects to string
    tag_id = ''.join(donor)
    
    #Searches string, pulls out pdb 4 character code
    tag_id1 = re.findall('....(?=\....)', donor)
    
    #Creates structure object for tag
    tags2 = parser.get_structure(tag_id1,tag_id)
    
    #Removes heteroatoms from structure
    tags3 = HETscrubber(tags2)
    
    #Removes flexible N/C terminal residues 
    Flextimer = time.time()
    tags4 = Flexscrubber(tags3)

    tags.append(tags4)
 
    return tags
 
def getFirstAndLastResiduesInChain(theChain):
    firstResid = -1
    lastResid  = -1
    firstResidObj = None
    lastResidObj = None 

    #Find the first and last residues in the PDB file
    for model in theChain[0]:  
	for chain in model:  
            for tag_resid in chain:
                if tag_resid.get_resname() in aminoAcids:
                    if ((firstResid == -1) or (tag_resid.get_id()[1] < firstResid)):
                        firstResid = tag_resid.get_id()[1]
                        firstResidObj = tag_resid
                    if ((lastResid == -1) or (tag_resid.get_id()[1] > lastResid)):
                        lastResid = tag_resid.get_id()[1]
                        lastResidObj = tag_resid

    return (firstResid,lastResid,firstResidObj,lastResidObj) 

def tagallchain(taggedProtein, CompositePDB1, chunk, AccDonIter, queue):
    TagChain1output = []

    for chains in chunk:
        AcceptorDonorIter2 = pickle.loads(AccDonIter)
        AcceptorIter2 = AcceptorDonorIter2[0]
        DonorIter2 = AcceptorDonorIter2[1]
        TagChain1 = AcceptorDonorIter2[2]
        
        #Creates a superimposer, aligns atoms in residue of tag chain A with current chain to get translational matrix, applies this to Acceptor list  
        Donors = []
        Acceptors = [] 
    
        #Pulls out the specific N/Cterminal Acceptor/Donor CA atom lists for use with alignment
        for residue in TagChain1.get_list():
            for atom in residue:
                for donor in DonorIter2:
                    if atom.get_full_id()[3][1] == donor.get_full_id()[3][1] and atom.get_id() == 'CA':
                        Donors.append(residue['CA'])
                       
        for residue in chains.get_list():
            for atom in residue:
                for acceptor in AcceptorIter2:
                    if atom.get_full_id()[3][1] == acceptor.get_full_id()[3][1] and atom.get_full_id()[4][0] == 'CA':
                        Acceptors.append(residue['CA'])
    
        #Create superimposer object for transformations
        super_imposerx = Bio.PDB.Superimposer()
    
        #Using newly transformed Acceptor list, aligns Donor list and applies this transformation to instance of TagChain
        super_imposerx.set_atoms(Acceptors, Donors)
        super_imposerx.apply(TagChain1)
        TagChain1.id = str(chains.id)+str('tag')
        TagChain1.detach_parent()
        
        TagChain1output.append(TagChain1)
        
    queue.put(TagChain1output)

def spheredefine(Structure, AccN, AccC, DonN, DonC, overlap):    
    #Center of sphere
    atomcoords = [0, 0, 0]

    #Radial maximum 
    atomradmax = [0, 0, 0]

    #Vector magnitude
    vectormax = 0
    
    #Total number of atoms iterated through
    atomcount = 0

    #This code block iterates through atoms, adds up x, y, z coords and divides them by max number to find the mean 
    for residue in Structure:
        CheckpointA = residue.get_id()[1] <= (AccN - 1 - overlap) or residue.get_id()[1] >= (AccC + 1 + overlap)
        if CheckpointA == True:
            for atom in residue:
                atomcoords[0] = atomcoords[0] + atom.get_coord()[0]
                atomcoords[1] = atomcoords[1] + atom.get_coord()[1]
                atomcoords[2] = atomcoords[2] + atom.get_coord()[2]
                atomcount = atomcount + 1 
    #'print 'heyhye', atomcoords[0], atomcount

    atomcoords[0] = atomcoords[0] / atomcount
    atomcoords[1] = atomcoords[1] / atomcount
    atomcoords[2] = atomcoords[2] / atomcount

    #This code block finds the magnitude of the vectors for each atom from the centre point, compares and returns the coords of the vector with the maximum magnitude 
    for residue in Structure:
        CheckpointC = residue.get_id()[1] <= (AccN - 1 - overlap) or residue.get_id()[1] >= (AccC + 1 + overlap)
        if CheckpointC ==True:
            for atom in residue:
                if (((atom.get_coord()[0] - atomcoords[0])**2 + (atom.get_coord()[1] - atomcoords[1])**2 + (atom.get_coord()[2] - atomcoords[2])**2)**0.5) > vectormax:
                    vectormax = (((atom.get_coord()[0] - atomcoords[0])**2 + (atom.get_coord()[1] - atomcoords[1])**2 + (atom.get_coord()[2] - atomcoords[2])**2)**0.5)
                    atomradmax[0] = atom.get_coord()[0]
                    atomradmax[1] = atom.get_coord()[1]
                    atomradmax[2] = atom.get_coord()[2]             

    return vectormax, atomcoords, Structure.get_id()

def sphereclash(Spherelist):
    
    #Iterates through list of spheres and performs a pairwise comparison to determine if two given spheres have the potential to clash in space
    noclash = None 
    whichclash = []
    for sphere in Spherelist:
        for sphere2 in Spherelist:
            if (sphere[1] != sphere2[1]) == True and noclash != False:
                spherecenter = (((sphere2[1][0] - sphere[1][0])**2 + (sphere2[1][1] - sphere[1][1])**2 + (sphere2[1][2] - sphere[1][2])**2)**0.5)
                rmax = sphere[0] + sphere2[0]                
                noclash = spherecenter > rmax
                listobject = []
                listobject.append(sphere)
                listobject.append(sphere2)
                listobject.append(noclash)
                whichclash.append(listobject)
            if (sphere[1] != sphere2[1]) == True and noclash == False:
                spherecenter = (((sphere2[1][0] - sphere[1][0])**2 + (sphere2[1][1] - sphere[1][1])**2 + (sphere2[1][2] - sphere[1][2])**2)**0.5)
                rmax = sphere[0] + sphere2[0]
                currentclash = spherecenter > rmax
                listobject = []
                listobject.append(sphere)
                listobject.append(sphere2)
                listobject.append(currentclash)
                whichclash.append(listobject)

    return noclash, whichclash 

def alignCoordinates(taggedProtein, potentialTag, start_tagged_residue, end_tagged_residue, start_potentialTag, end_potentialTag, tag_overlap_length):
    #Provides series of objects used for RMSD/clash checking
    OptimalName     = []
    OptimalDonor    = []
    OptimalAcceptor = []
    ClashCheck      = []
    courseclash     = []
    Bestfalse       = 10
    NoTrue          = 0
    global Checkpoint
    global CompositePDB 
    global Chainnumber
    CompositePDB = taggedProtein
    Chainnumber = 1


    #Creates a series of objects to append CA atoms too
    #First number refers to Nend(_,2,3), second refers to Cend(_,2,3)
    tagged_atoms11 = []
    tagged_atoms12 = []
    tagged_atoms13 = []
    tagged_atoms21 = []
    tagged_atoms22 = []
    tagged_atoms23 = []
    tagged_atoms31 = []
    tagged_atoms32 = []
    tagged_atoms33 = []

    tag_atoms11    = []
    tag_atoms12    = []
    tag_atoms13    = []
    tag_atoms21    = []
    tag_atoms22    = []
    tag_atoms23    = []
    tag_atoms31    = []
    tag_atoms32    = []
    tag_atoms33    = []

    #Creates a series of objects containing model/chain of tag protein(donor)
    for model in potentialTag:
        firstTagModel = model
        break
    for chain in firstTagModel:
        firstTagChain = chain
        break

    #Creates object containing residues of tag(donor), appends atoms from these residues to another object
    for residue in firstTagChain:
        tag_res = residue
        Nend    = range(start_potentialTag, start_potentialTag + tag_overlap_length)
        Cend    = range(end_potentialTag + 1 - tag_overlap_length, end_potentialTag + 1)
        Nend2   = range(start_potentialTag + 1, start_potentialTag + tag_overlap_length + 1)
        Cend2   = range(end_potentialTag - tag_overlap_length, end_potentialTag)
        Nend3   = range(start_potentialTag + 2, start_potentialTag + tag_overlap_length + 2)
        Cend3   = range(end_potentialTag - 1 - tag_overlap_length, end_potentialTag - 1)
        #print '1', residue, residue.id, residue.__contains__('CA')#, Nend, Cend, Nend2, Cend2, Nend3, Cend3
        id = tag_res.id
        #print id[1], id
        if id[1] in Nend and residue.__contains__('CA') == True:              
            tag_atoms11.append(tag_res['CA'])
            tag_atoms12.append(tag_res['CA'])
            tag_atoms13.append(tag_res['CA'])
        if id[1] in Cend and residue.__contains__('CA') == True:
            tag_atoms11.append(tag_res['CA'])
            tag_atoms21.append(tag_res['CA'])
            tag_atoms31.append(tag_res['CA'])
        if id[1] in Nend2 and residue.__contains__('CA') == True:
            tag_atoms21.append(tag_res['CA'])
            tag_atoms22.append(tag_res['CA'])
            tag_atoms23.append(tag_res['CA'])
        if id[1] in Cend2 and residue.__contains__('CA') == True:
            tag_atoms12.append(tag_res['CA'])
            tag_atoms22.append(tag_res['CA'])
            tag_atoms32.append(tag_res['CA'])
        if id[1] in Nend3 and residue.__contains__('CA') == True:
            tag_atoms31.append(tag_res['CA'])
            tag_atoms32.append(tag_res['CA'])
            tag_atoms33.append(tag_res['CA'])
        if id[1] in Cend3 and residue.__contains__('CA') == True:
            tag_atoms13.append(tag_res['CA'])
            tag_atoms23.append(tag_res['CA'])
            tag_atoms33.append(tag_res['CA'])
   
    #print 'Don', len(tag_atoms11), len(tag_atoms12), len(tag_atoms13), len(tag_atoms21), len(tag_atoms22), len(tag_atoms23), len(tag_atoms31), len(tag_atoms32), len(tag_atoms33)
    #print Nend, Cend, Nend2, Cend2, Nend3, Cend3

    allchains = []
    
    #Creates a series of objects containing model/chain/atom of tagged protein(acceptor)
    for model in taggedProtein:
        firstModel = model
        break
    for chain in firstModel:
        if chainregx != []:
            for regex in chainregx:
                if re.search('['+regex+']', chain.id):
                    firstChain = chain
                    break
        else:
            firstChain = chain 
            break
    
    #Creates an object containing each chain
    for chain in firstModel:
        if chainregx != []:
            for regex in chainregx:
                if re.search('['+regex+']', chain.id):
    	            allchains.append(chain)
                    break
        else:
            allchains.append(chain)
            break

    #Iterates through residues in chain, appends CA atoms of residues to object if they match the region of interest  
    for residue in firstChain:
        tagged_res  = residue
        Nend        = range(start_tagged_residue + 1 - tag_overlap_length, start_tagged_residue + 1)
        Cend        = range(end_tagged_residue, end_tagged_residue + tag_overlap_length)
        Nend2       = range(start_tagged_residue - tag_overlap_length, start_tagged_residue)
        Cend2       = range(end_tagged_residue + 1, end_tagged_residue + 1 + tag_overlap_length)
        Nend3       = range(start_tagged_residue - 1 - tag_overlap_length, start_tagged_residue - 1)
        Cend3       = range(end_tagged_residue + 2, end_tagged_residue + 2 + tag_overlap_length)
        
        id = tagged_res.id
        
        if id[1] in Nend:
            tagged_atoms11.append(tagged_res['CA'])
            tagged_atoms12.append(tagged_res['CA'])
            tagged_atoms13.append(tagged_res['CA']) 
        if id[1] in Cend:
            tagged_atoms11.append(tagged_res['CA'])
            tagged_atoms21.append(tagged_res['CA'])
            tagged_atoms31.append(tagged_res['CA'])
        if id[1] in Nend2:
            tagged_atoms21.append(tagged_res['CA'])
            tagged_atoms22.append(tagged_res['CA'])
            tagged_atoms23.append(tagged_res['CA'])
        if id[1] in Cend2:
            tagged_atoms12.append(tagged_res['CA'])
            tagged_atoms22.append(tagged_res['CA'])
            tagged_atoms32.append(tagged_res['CA'])
        if id[1] in Nend3:
            tagged_atoms31.append(tagged_res['CA'])
            tagged_atoms32.append(tagged_res['CA'])
            tagged_atoms33.append(tagged_res['CA'])
        if id[1] in Cend3:
            tagged_atoms13.append(tagged_res['CA'])
            tagged_atoms23.append(tagged_res['CA'])
            tagged_atoms33.append(tagged_res['CA'])
                    
    #print 'Acc', len(tagged_atoms11), len(tagged_atoms12), len(tagged_atoms13), len(tagged_atoms21), len(tagged_atoms22), len(tagged_atoms23), len(tagged_atoms31), len(tagged_atoms32), len(tagged_atoms33)
    #print Nend, Cend, Nend2, Cend2, Nend3, Cend3
    
    print ' '        
    print '~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
    print 'RMSD of aligned acceptor/donor pairs'
    print '~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
    
    #Iterates through all permutations of Acceptor triplets and Donor triplets, outputs RMSD for each
    for Acceptor in tagged_atoms11, tagged_atoms12, tagged_atoms13, tagged_atoms21, tagged_atoms22, tagged_atoms23, tagged_atoms31, tagged_atoms32, tagged_atoms33:           
        
        # Creates superimposer and imposes atoms from residues of tagged(acceptor)/tag regions(donor) 
        super_imposer = Bio.PDB.Superimposer()
        for Donor in tag_atoms11, tag_atoms12, tag_atoms13, tag_atoms21, tag_atoms22, tag_atoms23, tag_atoms31, tag_atoms32, tag_atoms33:
                  
            if len(Acceptor) == (tag_overlap_length * 2) and len(Donor) == (tag_overlap_length * 2):               
                #Superimposed the Acceptor/Donor of the permutations of chewed back ends and produces RMSDs   
                super_imposer.set_atoms(Acceptor, Donor)
                super_imposer.apply(potentialTag.get_atoms())
                RMSD = super_imposer.rms
                    
                #Produces a series of objects with N and C terminals of Acceptor/Donor for clash checking
                ANterm1 = Acceptor[tag_overlap_length - (tag_overlap_length - 1)]
                ANterm2 = ANterm1.get_parent()
                ANterm3 = ANterm2.id            
        
                ACterm1 = Acceptor[tag_overlap_length + (tag_overlap_length -1)]
                ACterm2 = ACterm1.get_parent()
                ACterm3 = ACterm2.id
        
                DNterm1 = Donor[tag_overlap_length - 1]
                DNterm2 = DNterm1.get_parent()
                DNterm3 = DNterm2.id 
        
                DCterm1 = Donor[tag_overlap_length]
                DCterm2 = DCterm1.get_parent()
                DCterm3 = DCterm2.id
       
                print RMSD
     
                #Compares RMSD of subsequent tags, checks for clashes of acceptor/donor 
                if RMSD < Checkpoint:
                    CompositePDB1 = None
                    CompositePDB1 = copy.deepcopy(CompositePDB)
                   
                    AcceptorDonor = []
                    AcceptorDonor.append(Acceptor)
                    AcceptorDonor.append(Donor)
                    AcceptorDonor.append(firstTagChain)
                    AcceptorDonorIter = pickle.dumps(AcceptorDonor)
                    
                    #Iterates through each chain of Donor and aligns Acceptor regions to the equivalent positions
                    Begin = time.time()
                    
                    if len(allchains) <= CPUs:
                        chainchunks = list(chunks(allchains, 1))
                                           
                    if len(allchains) > CPUs: 
                        chainchunks = list(chunks(allchains, (len(allchains) / CPUs)+1))   
                          
                    if __name__ == '__main__': 
                        queue = Queue()
                        for chunk in chainchunks:                                   
                            process = Process(target=tagallchain, args=(firstChain, CompositePDB1, chunk, AcceptorDonorIter, queue))
                            process.Daemon = True
                            process.start()
                        for chunk in chainchunks:
                            Tagchainadd = queue.get()
                            for chainx in Tagchainadd:
                                CompositePDB1.get_list()[0].add(chainx)
                                                        
                        process.join()                          
                        queue.close()
    
                    print 'Time taken (tagallchain) =', time.time() - Begin
                    
                    #Control statement to prevent runaway processes
                    while len(mp.active_children()) > 1:
                        time.sleep(0.1)
                        pass
    
                    allsphereslist = []          
     
                    Begin = time.time()
                    
                    for chain in CompositePDB1.get_list()[0].get_list():
                        sphereget = spheredefine(chain, ANterm3[1], ACterm3[1], DNterm3[1], DCterm3[1], tag_overlap_length)
                        allsphereslist.append(sphereget)
                        
                   # print 'Contents of allsphereslist', allsphereslist  
                      
                    print 'Time taken (sphereget) =', time.time() - Begin
                     
                    Begin = time.time()
                    courseclash = sphereclash(allsphereslist)
                    #print ' '
                    #print 'Contents of courseclash', courseclash 
                    
                    print 'Time taken (courseclash) =', time.time() - Begin
    
                    if courseclash[0] == False: 
                        Begin = time.time()
                        clash = findClash(CompositePDB1, ANterm3[1], ACterm3[1], allsphereslist, DNterm3[1], DCterm3[1], tag_overlap_length, courseclash[1])
                        if clash[0] == True and clash[1] == True:
                            clash = True
                        else:
                            clash = False
                        print 'Time taken (fineclash) =', time.time() - Begin 
                        print 'Fine clash passed?', clash
                        print ' '   
     
                    #If the RMSD is lower or equal to the current best and no clashes are detected, new return is created
                    if RMSD <= Checkpoint and (clash == True or courseclash == True):
                        OptimalName = potentialTag
                        print 'Donor =', OptimalName
                        Checkpoint = RMSD 
                        OptimalDonor = Donor
                        OptimalAcceptor = Acceptor                
                        ClashCheck = 'No clashes?', clash
                        print 'RMSD of current optimal Donor without clashes', Checkpoint
                        print ' '
                    else:
                        break
            else:
                print 'Residue missing during alignment: Acceptor = ', len(Acceptor), 'Donor = ', len(Donor)
                pass
            
    #Returns OptimalName = Optimal donor, OptimalRMSD = Optimal RMSD between OptimalDonor/Acceptor
    return OptimalName, Checkpoint, ClashCheck, OptimalDonor, OptimalAcceptor      

def HETscrubber(PDBstructure):
    #Iterates through list of amino acids, if heteroatom is found, residue is removed
    for chain in PDBstructure[0]:
        for residue in chain.get_list():
            resid = residue.id
            if resid[0] != ' ' or residue.__contains__('CA') == False or residue.get_resname() in aminoAcids == False:
                chain.detach_child(resid)
            
    return PDBstructure  

def Flexscrubber(PDBstructure):
    #Create an object to increase count to detect secondary structure
    Secstructure = 0
    Rescheck     = 0
    nreslist     = []
    creslist     = []
    
    #Base on fact adjacent Calpha distance = ~3.8A and Ca to Ca + 2 trans-peptide distance = 7.6A, if more than 2 within 6A, secondary structure detected
    for model in PDBstructure:
        for chain in model:
            Secstructure = 0
            for residue in chain:
                if residue.__contains__('CA') == True:
                    Rescheck = 0              
                    for residue2 in chain:
                        if residue2.__contains__('CA') == True and (residue != residue2) == True and Secstructure < 2 and residue2.get_id()[1] in range(residue.get_id()[1] + 1, residue.get_id()[1] + 3, 1):# or (residue.get_id()[1] + 2) or (residue.get_id()[1] + 3):
                            Rescheck = Rescheck + 1
                            if (residue['CA'] != residue2['CA']) == True and (residue['CA'] - residue2['CA']) > 4 and (residue['CA'] - residue2['CA']) < 6 and Secstructure < 2:
                                Secstructure = Secstructure + 1 
                            elif Secstructure < 2:
                                if residue not in nreslist:
                                    nreslist.append(residue)                       
                            elif Secstructure >= 2:
                                break
                        else:
                            pass
        
        #Performs the same algorithm on the inversed list of residues to search from C - terminus 
        for chain in model:
            Secstructure = 0
            for residue3 in list(reversed(chain.get_list())):
                if residue3.__contains__('CA') == True:
                    Rescheck = 0              
                    for residue4 in list(reversed(chain.get_list())):
                        if residue4.__contains__('CA') == True and (residue3 != residue4) == True and Secstructure < 2 and residue4.get_id()[1] in range(residue3.get_id()[1] - 3, residue3.get_id()[1] - 1, 1):
                            Rescheck = Rescheck + 1
                            if (residue3['CA'] != residue4['CA']) == True and (residue3['CA'] - residue4['CA']) > 4 and (residue3['CA'] - residue4['CA']) < 6 and Secstructure < 2:
                                Secstructure = Secstructure + 1
                            elif Secstructure < 2: 
                                if residue3 not in creslist:
                                    creslist.append(residue3)
                            elif Secstructure >= 2:
                                break
                        else:
                            pass
    
    #Outputs the degree to which the sequences were pruned
    print 'Tag pruned:'
    if nreslist != []:
        print 'N-terminal chewed by', len(nreslist), 'residues to', nreslist[len(nreslist)-1]
    if creslist != []:
        print 'C-terminal chewed by', len(creslist), 'residues to', creslist[len(creslist)-1]
    if nreslist == [] and creslist == []:
        print 'No change'
    
    
    for residue in nreslist:   
        residue.get_parent().detach_child(residue.id)
    for residue in creslist:   
        residue.get_parent().detach_child(residue.id)
               
    return PDBstructure                 
                                
def prunePDB(PDBstructure, Nprune, Cprune):
    #Chews back N/C terminal ends of PDB based on N/Cprune input, returns new structure
    for model in PDBstructure:
        for chain in model:
            for residue in chain:
                x = residue.id
                for x in Nprune:
                    if (' ', x, ' ') in chain:
                        chain.detach_child((' ', x, ' '))             
                for x in Cprune:
                    if (' ', x, ' ') in chain:
                        chain.detach_child((' ', x, ' '))

    return PDBstructure   

def AccDonsuperimposer(AccDon, AccAlign, DonAlign):
    #Superimposes the optimal Acceptor/Donor
    PDBinterest = AccDon   
 
    super_imposer = Bio.PDB.Superimposer()
    super_imposer.set_atoms(AccAlign, DonAlign)
    super_imposer.apply(PDBinterest.get_atoms())
    
    return PDBinterest

def AccPairGenerator(Acceptor, distance, nter, cter):
    
    NCtagsites = []
    
    #Generates a series of potential tag sites within the Acceptor delimited by the input N/C terminal residues
    #Uses the N/C-terminal distance of the Donor and pairs residues which match the same distance +- a wiggle factor
    for Model in Acceptor:
        for Chain in Model:
            for Residue1 in Chain:
                if (Residue1.id[1] > nter) == True and (Residue1.id[1] < cter) == True:
                    for Residue2 in Chain:
                        if (Residue2.id[1] > nter) == True and (Residue2.id[1] < cter) == True:
                            DistCheck = Residue1['CA'] - Residue2['CA']
                            if (DistCheck > (int(distance-0.5))) == True and (DistCheck < (int(distance+0.5))) == True and (Residue2,Residue1) not in NCtagsites and (Residue1,Residue2) not in NCtagsites:
                                NCtagsites.append((Residue1,Residue2))
                        else:
                            pass
                else:
                    pass
            break  
        
    return NCtagsites

def main():
    #Creates time object to compare runtimes
    Startime = time.time()
    global CPUs, threshold_A

    #Defines arguments to be included when running the program
    #Example - ./protagger_v0_17.py -d ./1BTL_ecoli_tag -p YYYY -n 190 -c 283 -o 3  -r 3 -u 3 -l -x A-C

    parser= argparse.ArgumentParser(description='Find the optimal tag for your protein')

    #Adds parser arguments which the user will input
    parser.add_argument('-d', '--tagdatabase', type=str, help='Enter the folder containing the tags to be tested', required = True)
    parser.add_argument('-p', '--pdb', type=str, help='Enter pdb to tag (lower case only)', required = True)
    parser.add_argument('-n', '--nterm', type=int, help='Enter N-terminal residue number of tag site', required = True)
    parser.add_argument('-c', '--cterm', type=int, help='Enter C-terminal residue number of tag site', required = True)
    parser.add_argument('-o', '--overlap', type=int, help='Enter overlap length used for N/C terminal alignment', required = True)
    parser.add_argument('-r', '--RMSDthreshold', type=float, help='Enter RMSD threshold under which to check. Recommended value = 3, larger values increase runtime', required = True)
    parser.add_argument('-u', '--CPU', type=int, help='Enter the number of CPUs to use during parallel processing steps. Multiple cores recommended for greater numbers of chains to tag.', required = True)
    parser.add_argument('-t', '--atomthreshold', type=float, help='Enter the threshold (in angstrom) of detection for clashes between neighboring atoms. Default is 2.5, user higher values for stricter clash check (e.g. 4)')
    parser.add_argument('-l', '--local', action='store_true', help='The PDB to be tagged is stored locally')
    parser.add_argument('-s', '--specific', action='store_true', help='Find the ideal tagsite for a specific protein within range N/C terminals')
    parser.add_argument('-x', '--chain', type=str, help='Enter chain identifiers of those you wish to tag. E.g. A-C,F,X')

    #Creates an object containing parsed argument names
    args = parser.parse_args()

    #Creates objects with values of parsed user input
    local1          = args.local
    specsearch      = args.specific     
    pdb_acceptor    = args.pdb
    chainlabel      = args.chain
    nterm           = args.nterm
    cterm           = args.cterm
    overlap         = args.overlap
    CPUs            = args.CPU
    RMSDthreshold   = args.RMSDthreshold
    threshold_A_in  = args.atomthreshold
    Databaseloc     = args.tagdatabase + '/*'

    if threshold_A_in != None and threshold_A_in < 2.5:
        print "ERROR: atomthreshold must be 2.5 or greater"
        sys.exit()
    elif threshold_A_in != None and threshold_A_in >= 2.5:
        threshold_A = threshold_A_in
    else:
        threshold_A = 2.5

    #Returns user input
    print "Database of tags =", Databaseloc
    print "Target pdb =", pdb_acceptor
    print "PDB stored locally? =", local1 == True  
    print "Run specific search? =", specsearch == True
    print "N-terminal of insertion =", nterm
    print "C-terminal of insertion =", cterm
    print "Overlap alignment length =", overlap 
    print "RMSD threshold =", RMSDthreshold 
    print "Atomic clash threshold =", threshold_A
    print "Number of CPUs used for processing=", CPUs  
    if chainlabel != None:
        print "ID's of chains to tag =", chainlabel 
    else:
        print "ID's of chains to tag = all chains"

    #Creates an object used to generate a run specific log file
    logname = 'Protagger_'+str(pdb_acceptor)+"_N_"+str(nterm)+"_C_"+str(cterm)+"_overlap_"+str(overlap)+"_RMSD_"+str(RMSDthreshold)+"_atom_threshold_"+str(threshold_A)+".log"
    log = open(logname,"w")
    log.write(logname+"\n")
    
    #Creates empty objects for use in parsetag
    firstlast = []

    #Returns structure object of protein to be tagged (acceptor) and removes hetero atoms 
    if local1 == False:
        Acceptor1 = downloadpdb(pdb_acceptor)
    else:
        Acceptor1 = localpdb(pdb_acceptor)

    Acceptor = HETscrubber(Acceptor1)
    print 'Acceptor acquired'

    #Returns list of structure object tags (donors)
    taglist = glob.glob(Databaseloc)
    print "Donor list made"
    
    #Creates objects to monitor remaining tags
    tagsrem = len(taglist)
    tagstot = len(taglist)
    
    #Create a series of objects for subsequent steps
    global Checkpoint
    global Optimal    
    global chainregx
    Checkpoint  = RMSDthreshold
    Besttag     = []
    Notagsfound = 0
    
    #Creates the list containing parameters for the regular expression for chain IDs
    if chainlabel != None and ((',' in chainlabel) == True):
        chainregx = str(chainlabel).split(",")
    elif chainlabel != None and ((',' in chainlabel) != True):
        chainregx.append(str(chainlabel))
    else:
        chainregx.append('\D')
    
    #Iterate through donors and apply tagging algorithms
    for donor in taglist:
        Begin = time.time()
        tagsrem = tagsrem - 1
        print " "
        print "################################################################"
        print " "
        print "Tags remaining =", tagsrem, '/', tagstot
        print "Elapsed time (minutes) =",  ((time.time() - Startime)/60)
        print "Current tag being checked =", donor
        if Optimal != []:
            print "Current optimal =", Optimal[0]
        else:
            print "Current optimal = None"
        print "RMSD threshold =", Checkpoint
        print "\n"
        
        #Finds first and last residues in donor chains 
        Donor = parsetag(donor)
        firstlast = []
        firstlast.append(getFirstAndLastResiduesInChain(Donor))
          
        #Makes object for storing Acceptor tag sites and generates threshold criteria from Donor N/C terminals  
        NCtagsites = []
        firlasdist = firstlast[0][2]['CA']-firstlast[0][3]['CA']
        
        #Appends standard input N/C terminal from user input of Acceptor
        for Model in Acceptor:
            for Chain in Model:
                Nres = Chain.__getitem__(nterm)
                Cres = Chain.__getitem__(cterm)
                break
        
        NCtagsites.append((Nres,Cres))
        
        #If deep search algorithm is to be used, various Acceptor sites for tagging are generated
        if specsearch == True:
            NCtagsites = AccPairGenerator(Acceptor, firlasdist, nterm, cterm)
            
        #Object to monitor number of insert sites checked
        tupnum = None
        tuptop = len(NCtagsites)
        tuprem = len(NCtagsites)
        
        #Iterates through Acceptor tag sites and applies     
        for tup in NCtagsites:
            print ' '
            print 'Insert sites remaining =', tuprem, '/', tuptop 
            print 'Current insert site N/C tested =', tup 
            #print tup, tup[0].id[1], tup[1].id[1]
            tuprem = tuprem - 1
            
            Donorfl = []
            Donorfl = zip(Donor, firstlast)
            
            global Optimal, AccepTup
            Begintup = time.time()
            #Iterates through tag tuples and aligns the N/C terminals of each with the corresponding tag site
            #Outputs an object containing information of optimal tag
            
            output = alignCoordinates(Acceptor, Donorfl[0][0], tup[0].id[1], tup[1].id[1], Donorfl[0][1][0], Donorfl[0][1][1], overlap)
            print ' '
            print '~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
            print 'Time taken (insert site tested) =', time.time() - Begintup, 'Donor =', donor
            print 'Tag site N =', tup[0].id[1], 'Tag site C =', tup[1].id[1], 'Current checkpoint =', Checkpoint
            print '~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
            if output[1] == Checkpoint and output[0] != []:
                global Optimal
                Optimal     = []
                Optimal     = output 
                Checkpoint  = output[1]    
                Besttag     = Donorfl
                Notagsfound = Notagsfound + 1
                AccepTup    = []
                AccepTup    = tup
                print 'New optimal found! =', Optimal[0], Optimal[1], Optimal[2]
                
                #Output every potential tag found to file
                #log = open(logname,"w")
                log.write(" "+"\n\n")
                log.write("# " + str(Notagsfound) + " Name = " + str(Optimal[0]) + ", RMSD = " + str(Optimal[1]) + "\n")
                log.write("Donor N term = " + str(Optimal[3][0].get_parent()) + "\n")
                log.write("Donor C term = " + str(Optimal[3][overlap * 2 - 1].get_parent()) + "\n")
                log.write("Acceptor N term = " + str(Optimal[4][overlap - 1].get_parent()) + "\n")
                log.write("Acceptor C term = " + str(Optimal[4][overlap].get_parent()))
        print ' '
        print '========================================================================================================='
        print 'Time taken (tag tested) =', time.time() - Begin, 'Donor =', donor, 'Current checkpoint =', Checkpoint
        print '========================================================================================================='
    #If no tags are found, quits program and returns error message 
    if Optimal == []:
        print 'Time taken (failed) =', time.time() - Startime
        sys.exit("~~~~~~~~~~ No optimal tag found :( ~~~~~~~~~")

    #Wait 120 seconds to allow other processes to finish up
    time.sleep(60)
    
    print 'Optimal tag =', Optimal[0], Optimal[1], Optimal[2]
    
    #Creates a series of objects where the residue id is pulled out for the NCterms of Acceptor/Donor
    NAccterm1   = Optimal[4][overlap - 1]
    NAccterm2   = NAccterm1.get_parent() 
    NAccterm    = NAccterm2.id

    CAccterm1   = Optimal[4][overlap]
    CAccterm2   = CAccterm1.get_parent()
    CAccterm    = CAccterm2.id

    NDonterm1   = Optimal[3][0]
    NDonterm2   = NDonterm1.get_parent()
    NDonterm    = NDonterm2.id

    CDonterm1   = Optimal[3][overlap * 2 - 1]
    CDonterm2   = CDonterm1.get_parent()
    CDonterm    = CDonterm2.id

    TruncDonN = (NDonterm[1]+overlap)
    TruncDonC = (CDonterm[1]-overlap)
    
    #Prints the N/C terminals of Acceptor/Donor for insertion
    print " "
    print "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~Best Tag~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
    print "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
    print " RMSD match of N/C overlap = ", Optimal[1], "                                          "   
    print "                                                                                       "                                                
    print " Acceptor = ", Acceptor, "                                                             " 
    print " Acceptor optimal N terminal start residue is", NAccterm2, "                           " 
    print " Acceptor optimal C terminal end residue is", CAccterm2, "                             "
    print "                                                                                       "
    print " Donor = ", Optimal[0], "                                                              "
    print " Donor optimal N terminal start residue is", NDonterm2, "                              "
    print " Donor optimal C terminal end residue is", CDonterm2, "                                " 
    print "                                                                                       "
    print " Protagger suggests truncating Donor to residue numbers:                               "                                                                                                             
    print " N terminal = ", TruncDonN, "                                                          "
    print " C terminal = ", TruncDonC, "                                                          "
    print "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
    print "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~" 

    #Writes out best tag results to file
    log.write("\n")
    log.write("\n")
    log.write("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~Best Tag~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"+"\n")
    log.write("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"+"\n")
    log.write(" RMSD match of N/C overlap = " + str(Optimal[1]) + "                             "+"\n")   
    log.write("                                                                                 "+"\n")                                                
    log.write(" Acceptor = " + str(Acceptor) + "                                                "+"\n") 
    log.write(" Acceptor optimal N terminal start residue is" + str(NAccterm2) + "              "+"\n") 
    log.write(" Acceptor optimal C terminal end residue is"+ str(CAccterm2) + "                 "+"\n")
    log.write("                                                                                 "+"\n")
    log.write(" Donor = " + str(Optimal[0]) + "                                                 "+"\n")
    log.write(" Donor optimal N terminal start residue is" + str(NDonterm2) + "                 "+"\n")
    log.write(" Donor optimal C terminal end residue is" + str(CDonterm2) + "                   "+"\n") 
    log.write("                                                                                 "+"\n")
    log.write(" Protagger suggests truncating Donor to residue numbers:                         "+"\n")
    log.write(" N terminal = "+  str((NDonterm[1]+overlap)) + "                                 "+"\n")
    log.write(" C terminal = "+ str((CDonterm[1]-overlap)) + "                                  "+"\n")
    log.write("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"+"\n")
    log.write("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"+"\n") 
    log.write("\n")
    log.write('Total time taken (seconds) = '+ str(time.time() - Startime) + "\n")
    log.close() 
     
    #Uses the above objects with unwanted residue ids to pull out all residue ids to be deleted
    NdeleteAcc = range(NAccterm[1] + 1, CAccterm[1])
    CdeleteAcc = range(AccepTup[1].id[1] -1, CAccterm[1])
    NdeleteDon = range(Besttag[0][1][0] - 1, NDonterm[1])
    CdeleteDon = range(CDonterm[1] + 1, Besttag[0][1][1] + 1) 
  
    #Chews back residues of the Acceptor/Donor to provide structure files with optimal residues
    PrunedDonor = prunePDB(Optimal[0], NdeleteDon, CdeleteDon)
    PrunedAcceptor = prunePDB(Acceptor, NdeleteAcc, CdeleteAcc)    
    print 'Optimal acceptor/donor pruned'

    OptimalAcceptorfinal = AccDonsuperimposer(PrunedAcceptor, Optimal[4], Optimal[3])
    OptimalDonorfinal = AccDonsuperimposer(PrunedDonor, Optimal[4], Optimal[3])
 
    AccDonIter = []
    AccDonIter.append(Optimal[4])
    AccDonIter.append(Optimal[3])
    AccDonIter.append(OptimalDonorfinal)
    AccDonIter = pickle.dumps(AccDonIter)
    
    Begin = time.time()

    #Pulls name of optimal donor tag
    Donametmp = re.search("'....'", str(Optimal[0]))
    Doname = (Donametmp.group(0)).strip("'")

    #This section appends a tag to each of the monomers in the multimer, saves the output pdb
    for chains in OptimalAcceptorfinal.get_list()[0].get_list():
        if chainregx != []:
            for regex in chainregx:
                if re.search('['+regex+']', chains.id):
                    AcceptorDonorIter2 = pickle.loads(AccDonIter)
                    AcceptorIter2 = AcceptorDonorIter2[0]
                    DonorIter2 = AcceptorDonorIter2[1]
                    TagModel = AcceptorDonorIter2[2][0].get_list()[0]
            
                    Donors = []
                    Acceptors = []
            
                    #Pulls out the specific N/Cterminal Acceptor/Donor CA atom lists for use with alignment
                    for residue in TagModel:
                        for atom in residue: 
                            for donor in DonorIter2:
                                if atom.get_full_id()[3][1] == donor.get_full_id()[3][1] and atom.get_id() == 'CA':
                                    Donors.append(atom)
                   
                    for residue in chains.get_list():
                        for atom in residue:
                            for acceptor in AcceptorIter2:
                                if atom.get_full_id()[3][1] == acceptor.get_full_id()[3][1] and atom.get_id() == 'CA':
                                    Acceptors.append(residue['CA'])
            
                    #Using the pulled out lists, superimposer aligns and applies transformation to atoms in tag chain
                    super_imposerx = None
                    super_imposerx = Bio.PDB.Superimposer()
                    super_imposerx.set_atoms(Acceptors, Donors)
                    super_imposerx.apply(TagModel.get_atoms())
            
                    #Saves the new tag chain
                    io=Bio.PDB.PDBIO()
                    io.set_structure(TagModel)
                    io.save((str(pdb_acceptor)+'_'+Doname+'_'+str(chains.get_id())+'tag.pdb').replace(" ", ""))
    
                    print str(chains.get_id()),'tag.pdb saved'
                    
        elif chainregx == []:
            AcceptorDonorIter2 = pickle.loads(AccDonIter)
            AcceptorIter2 = AcceptorDonorIter2[0]
            DonorIter2 = AcceptorDonorIter2[1]
            TagModel = AcceptorDonorIter2[2][0].get_list()[0]
            
            Donors = []
            Acceptors = []
            
            #Pulls out the specific N/Cterminal Acceptor/Donor CA atom lists for use with alignment
            for residue in TagModel:
                for atom in residue: 
                    for donor in DonorIter2:
                        if atom.get_full_id()[3][1] == donor.get_full_id()[3][1] and atom.get_id() == 'CA':
                            Donors.append(atom)
                   
            for residue in chains.get_list():
                for atom in residue:
                    for acceptor in AcceptorIter2:
                        if atom.get_full_id()[3][1] == acceptor.get_full_id()[3][1] and atom.get_id() == 'CA':
                            Acceptors.append(residue['CA'])
            
            #Using the pulled out lists, superimposer aligns and applies transformation to atoms in tag chain
            super_imposerx = None
            super_imposerx = Bio.PDB.Superimposer()
            super_imposerx.set_atoms(Acceptors, Donors)
            super_imposerx.apply(TagModel.get_atoms())
            
            #Saves the new tag chain
            io=Bio.PDB.PDBIO()
            io.set_structure(TagModel)
            io.save((str(pdb_acceptor)+'_'+Doname+'_'+str(chains.get_id())+'tag.pdb').replace(" ", ""))
    
            print str(chains.get_id()),'tag.pdb saved'

    #Saves the PDB of the optimal final acceptor 
    io=Bio.PDB.PDBIO()
    io.set_structure(OptimalAcceptorfinal)
    io.save(pdb_acceptor+'_optimised_for_'+Doname+'.pdb')
    
    print 'Optimised', pdb_acceptor, '.pdb saved'
    print ' '
    print 'Optimal acceptor/donor pairs generated and saved'
    print 'Protagger is done!'
    print 'Time taken (total) =', time.time() - Startime

main()
