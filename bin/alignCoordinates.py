import globalvars
import re
import copy
import cPickle as pickle
import time
import multiprocessing as mp
from multiprocessing import Pool, Pipe, Queue, Process
from Bio import PDB
import Bio.PDB
from Bio.PDB import PDBList
from Bio.PDB import PDBParser
from Bio.PDB import StructureBuilder
from Bio.PDB import Selection
from Bio.PDB import Entity
from Bio.PDB import PDBIO

from sphereclash import sphereclash
from spheredefine import spheredefine
from findClash import findClash
from tagallchain import tagallchain

def chunks(list, n):
    for i in range(0, len(list), n):
        yield list[i:i+n]

def alignCoordinates(taggedProtein, potentialTag, start_tagged_residue, end_tagged_residue, start_potentialTag, end_potentialTag, tag_overlap_length):

    #Provides series of objects used for RMSD/clash checking
    OptimalName     = []
    OptimalDonor    = []
    OptimalAcceptor = []
    ClashCheck      = []
    courseclash     = []
    Bestfalse       = 10
    NoTrue          = 0

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
        if globalvars.chainregx != []:
            for regex in globalvars.chainregx:
                if re.search('['+regex+']', chain.id):
                    firstChain = chain
                    break
        else:
            firstChain = chain
            break

    #Creates an object containing each chain
    for chain in firstModel:
        if globalvars.chainregx != []:
            for regex in globalvars.chainregx:
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
                if RMSD < globalvars.Checkpoint:
                    CompositePDB1 = None
                    CompositePDB1 = copy.deepcopy(CompositePDB)

                    AcceptorDonor = []
                    AcceptorDonor.append(Acceptor)
                    AcceptorDonor.append(Donor)
                    AcceptorDonor.append(firstTagChain)
                    AcceptorDonorIter = pickle.dumps(AcceptorDonor)

                    #Iterates through each chain of Donor and aligns Acceptor regions to the equivalent positions
                    Begin = time.time()

                    if len(allchains) <= globalvars.CPUs:
                        chainchunks = list(chunks(allchains, 1))

                    if len(allchains) > globalvars.CPUs:
                        chainchunks = list(chunks(allchains, (len(allchains) / globalvars.CPUs)+1))

                    if __name__ == 'alignCoordinates':
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
                    if RMSD <= globalvars.Checkpoint and (clash == True or courseclash == True):
                        OptimalName = potentialTag
                        print 'Donor =', OptimalName
                        globalvars.Checkpoint = RMSD
                        OptimalDonor = Donor
                        OptimalAcceptor = Acceptor
                        ClashCheck = 'No clashes?', clash
                        print 'RMSD of current optimal Donor without clashes', globalvars.Checkpoint
                        print ' '
                    else:
                        break
            else:
                print 'Residue missing during alignment: Acceptor = ', len(Acceptor), 'Donor = ', len(Donor)
                pass

    #Returns OptimalName = Optimal donor, OptimalRMSD = Optimal RMSD between OptimalDonor/Acceptor
    return OptimalName, globalvars.Checkpoint, ClashCheck, OptimalDonor, OptimalAcceptor
