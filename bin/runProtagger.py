#!/usr/bin/python

#Import dependencies used by main and support scripts
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
import Bio.PDB
from Bio.PDB import *
import subprocess

#Import Protagger functions from support scripts
from alignCoordinates import alignCoordinates, chunks
from downloadpdb import downloadpdb
from findClash import findClash
from flexscrubber import Flexscrubber
from getFirstLast import getFirstAndLastResiduesInChain
from HETscrubber import HETscrubber
from localpdb import localpdb
from pairgenerator import AccPairGenerator
from parsetag import parsetag
from prunePDB import prunePDB
from restruecheck import restruecheck
from sphereclash import sphereclash
from spheredefine import spheredefine
from superimposer import AccDonsuperimposer
from tagallchain import tagallchain
from databaseorder import databaseorder
import globalvars

#Initialise global values used by modules
globalvars.initialise()

#Creates time object to compare runtimes
globalvars.Startime = time.time()

#Defines arguments to be included when running the program
#Example - ./protagger_v0_17.py -d ./1BTL_ecoli_tag -p YYYY -n 190 -c 283 -o 3  -r 3 -u 3 -l -x A-C

parser= argparse.ArgumentParser(description='Find the optimal tag for your protein')

#Adds parser arguments which the user will input
parser.add_argument('-b', '--tagdatabase', type=str, help='Enter the folder containing the tags to be tested', required = True)
parser.add_argument('-o', '--outfolder', type=str, help='Enter the name of the output folder')
parser.add_argument('-p', '--pdb', type=str, help='Enter pdb to tag (lower case only)', required = True)
parser.add_argument('-n', '--nterm', type=int, help='Enter N-terminal residue number of tag site', required = True)
parser.add_argument('-c', '--cterm', type=int, help='Enter C-terminal residue number of tag site', required = True)
parser.add_argument('-v', '--overlap', type=int, help='Enter overlap length used for N/C terminal Ca alignment, default value is 3')
parser.add_argument('-r', '--RMSDthreshold', type=float, help='Enter RMSD threshold under which to check. Recommended value = 3, larger values increase runtime', required = True)
parser.add_argument('-u', '--CPU', type=int, help='Enter the number of CPUs to use during parallel processing steps. Multiple cores recommended for greater numbers of chains to tag.', required = True)
parser.add_argument('-t', '--atomthreshold', type=float, help='Enter the threshold (in angstrom) of detection for clashes between neighboring atoms. Default/minimum is 2.5, use higher values for stricter clash check (e.g. 4)')
parser.add_argument('-l', '--local', action='store_true', help='The PDB to be tagged is stored locally')
parser.add_argument('-d', '--deepsearch', action='store_true', help='Find the ideal tagsite for a specific protein within range N/C terminals')
parser.add_argument('-x', '--chain', type=str, help='Enter chain identifiers of those you wish to tag. E.g. A-C,F,X')
parser.add_argument('-e', '--RMSDexit', type=float, help='Enter RMSD threshold under which terminate search. Recommended value')
parser.add_argument('-m', '--MinLengthDoner', type=int, help='Enter minimum residue length of pruned Donor chains')



#Creates an object containing parsed argument names
args = parser.parse_args()

#Creates objects with values of parsed user input
local1          = args.local
outfolder       = "./" + args.outfolder + "/"
specsearch      = args.deepsearch
pdb_acceptor    = args.pdb
chainlabel      = args.chain
nterm           = args.nterm
cterm           = args.cterm
overlap         = args.overlap
globalvars.CPUs = args.CPU
RMSDthreshold   = args.RMSDthreshold
threshold_A_in  = args.atomthreshold
Databaseloc     = args.tagdatabase + '/*'
RMSDexit        = args.RMSDexit
Minlength	= args.MinLengthDoner

#Primes optional user input variables with default values
if threshold_A_in != None and threshold_A_in < 2.5:
    print "ERROR: atomthreshold must be 2.5 or greater"
    sys.exit()
elif threshold_A_in != None and threshold_A_in >= 2.5:
    globalvars.threshold_A = threshold_A_in
else:
    globalvars.threshold_A = 2.5

if overlap == None:
    overlap = 3
elif overlap < 1:
    print "Overlap must be > 1, defaulting to 3"
    overlap = 3

if RMSDexit == None:
    RMSDexit = 0
elif RMSDexit > RMSDthreshold:
    print "RMSDexit must be lower than RMSDthreshold, setting to 0..."
    RMSDexit = 0

if Minlength == None:
    Minlength = 0

#Returns user input variables
print "Database of tags =", Databaseloc
print "Target pdb =", pdb_acceptor
print "PDB stored locally? = ", local1 == True
print "Run specific search? = ", specsearch == True
print "N-terminal of insertion =", nterm
print "C-terminal of insertion =", cterm
print "Overlap alignment length =", overlap
print "RMSD search threshold =", RMSDthreshold
print "RMSD exit threshold =", RMSDexit
print "Atomic clash threshold =", globalvars.threshold_A
print "Number of CPUs used for processing=", globalvars.CPUs
print "Minimum length of pruned Donor chain =", Minlength
if chainlabel != None:
    print "ID's of chains to tag =", chainlabel
else:
    print "ID's of chains to tag = all chains"


#Creates an object used to generate a run specific log file
subprocess.call(["mkdir", outfolder])
logname = 'Protagger_'+str(pdb_acceptor)+"_N_"+str(nterm)+"_C_"+str(cterm)+"_RMSD_"+str(int(RMSDthreshold))+".log"
log = open(outfolder + logname,"w")
log.write(logname + "\n" +
"Database of tags = " + str(Databaseloc) + "\n" +
"Target pdb = " + str(pdb_acceptor) + "\n" +
"PDB stored locally? =" + str(local1 == True) + "\n" +
"Run specific search? =" + str(specsearch == True) + "\n" +
"N-terminal of insertion = " + str(nterm) + "\n" +
"C-terminal of insertion = " + str(cterm) + "\n" +
"Overlap alignment length = " + str(overlap) + "\n" +
"RMSD threshold = " + str(RMSDthreshold) + "\n" +
"Atomic clash threshold = " + str(globalvars.threshold_A) + "\n" +
"Number of CPUs used for processing = " + str(globalvars.CPUs) + "\n"
"Minimum length of pruned Donor chain =" + str(Minlength) + "\n")


#Creates empty objects for use in parsetag
firstlast = []

#Returns structure object of protein to be tagged (acceptor) and removes hetero atoms
if local1 == False:
    Acceptor1 = downloadpdb(pdb_acceptor)
else:
    Acceptor1 = localpdb(pdb_acceptor)

Acceptor = HETscrubber(Acceptor1)
print 'Acceptor acquired'

#Creates the list containing parameters for the regular expression for chain IDs
if chainlabel != None and ((',' in chainlabel) == True):
    globalvars.chainregx = str(chainlabel).split(",")
elif chainlabel != None and ((',' in chainlabel) != True):
    globalvars.chainregx.append(str(chainlabel))
else:
    globalvars.chainregx.append('\D')

#Returns list of structure object tags (donors)
#Ranks list of potential tags based on pairwise distance of N/C terminals of donors
Acchains = Acceptor.get_chains()
accfirstlast = []
for chain in Acchains:
    if globalvars.chainregx != []:
        for regex in globalvars.chainregx:
            if re.search('['+regex+']', chain.id):
                print chain.id
                #print chain.id
                for residue in chain:
                #    print residue.id
                    accfirstlast.append(chain.__getitem__(nterm))
                    accfirstlast.append(chain.__getitem__(cterm))
                    break
            else:
                break
    else:
        #print chain.id
        for residue in chain:
        #    print residue.id
            accfirstlast.append(chain.__getitem__(nterm))
            accfirstlast.append(chain.__getitem__(cterm))
            break

print "yo", accfirstlast

print accfirstlast[0].get_full_id()
print accfirstlast[1].get_full_id()

accfirlasdist = accfirstlast[0]['CA']-accfirstlast[1]['CA']
print "ya", accfirlasdist
taglist = databaseorder(Databaseloc, accfirlasdist)
print taglist
print "Donor list made"

#Creates objects to monitor remaining tags
tagsrem = len(taglist)
tagstot = len(taglist)

#Create a series of objects for subsequent steps
globalvars.Checkpoint  = RMSDthreshold
Besttag     = []
Notagsfound = 0

#Iterate through donors and apply tagging algorithms
for donor in taglist:
    Begin = time.time()
    tagsrem = tagsrem - 1
    print " "
    print "################################################################"
    print " "
    print "Tags remaining =", tagsrem, '/', tagstot
    print "Elapsed time (minutes) =",  ((time.time() - globalvars.Startime)/60)
    print "Current tag being checked =", donor
    if globalvars.Optimal != []:
        print "Current optimal =", globalvars.Optimal[0]
    else:
        print "Current optimal = None"
    print "RMSD threshold =", globalvars.Checkpoint
    print "\n"

    if globalvars.Checkpoint > RMSDexit:
        try:
            #Attempts to parse tag, skips tag if error found
            Donor = parsetag(donor)
            #print "This is a Donor", Donor
	   # print "run test length", Donor[0].get_list()[0].get_list()[0].get_list()

            if len(Donor[0].get_list()[0].get_list()[0].get_list()) < Minlength:
		        print "Donor length below minimum threshold, aborting loop..."
		        continue

            #Identifies first/last residues in Donors
            firstlast = []
            firstlast.append(getFirstAndLastResiduesInChain(Donor))

            #Makes object for storing Acceptor tag sites and generates threshold criteria from Donor N/C terminals
            NCtagsites = []
            firlasdist = firstlast[0][2]['CA']-firstlast[0][3]['CA']

            #Appends standard input N/C terminal from user input of Acceptor
            #for Model in Acceptor:
            #    for Chain in Model:
            #        Nres = Chain.__getitem__(nterm)
            #        Cres = Chain.__getitem__(cterm)
            #        break

            for Model in Acceptor:
                for chain in Model:
                    if globalvars.chainregx != []:
                        for regex in globalvars.chainregx:
                            if re.search('['+regex+']', chain.id):
                                Nres = chain.__getitem__(nterm)
                                Cres = chain.__getitem__(cterm)
                                break
                    else:
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

                #global globalvars.Optimal, AccepTup
                Begintup = time.time()
                #Iterates through tag tuples and aligns the N/C terminals of each with the corresponding tag site
                #Outputs an object containing information of optimal tag

                output = alignCoordinates(Acceptor, Donorfl[0][0], tup[0].id[1], tup[1].id[1], Donorfl[0][1][0], Donorfl[0][1][1], overlap)
                print ' '
                print '~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
                print 'Time taken (insert site tested) =', time.time() - Begintup, 'Donor =', donor
                print 'Tag site N =', tup[0].id[1], 'Tag site C =', tup[1].id[1], 'Current checkpoint =', globalvars.Checkpoint
                print '~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
                if output[1] == globalvars.Checkpoint and output[0] != []:
                    #global Optimal
                    globalvars.Optimal     = []
                    globalvars.Optimal     = output
                    globalvars.Checkpoint  = output[1]
                    Besttag     = Donorfl
                    Notagsfound = Notagsfound + 1
                    globalvars.AccepTup    = []
                    globalvars.AccepTup    = tup
                    print 'New optimal found! =', globalvars.Optimal[0], globalvars.Optimal[1], globalvars.Optimal[2]

                    #Output every potential tag found to file
                    #log = open(logname,"w")
                    log.write(" "+"\n\n")
                    log.write("# " + str(Notagsfound) + " Name = " + str(globalvars.Optimal[0]) + ", RMSD = " + str(globalvars.Optimal[1]) + "\n")
                    log.write("Donor N term = " + str(globalvars.Optimal[3][0].get_parent()) + "\n")
                    log.write("Donor C term = " + str(globalvars.Optimal[3][overlap * 2 - 1].get_parent()) + "\n")
                    log.write("Acceptor N term = " + str(globalvars.Optimal[4][overlap - 1].get_parent()) + "\n")
                    log.write("Acceptor C term = " + str(globalvars.Optimal[4][overlap].get_parent()))
            print ' '
            print '========================================================================================================='
            print 'Time taken (tag tested) =', time.time() - Begin, 'Donor =', donor, 'Current checkpoint =', globalvars.Checkpoint
            print '========================================================================================================='
        except:
            print " "
            print "Tag", donor, "parsed incorrectly, skipping..."
    else:
        print "Exit threshold met, aborting further tests"
        break
#If no tags are found, quits program and returns error message
if globalvars.Optimal == []:
    print 'Time taken (failed) =', time.time() - globalvars.Startime
    sys.exit("~~~~~~~~~~ No optimal tag found :( ~~~~~~~~~")

#Wait 60 seconds to allow other processes to finish up
#time.sleep(60)

print 'Optimal tag =', globalvars.Optimal[0], globalvars.Optimal[1], globalvars.Optimal[2]

#Creates a series of objects where the residue id is pulled out for the NCterms of Acceptor/Donor
Nfuse = False

#for i in globalvars.Optimal[4]:
#    print i.get_full_id()

if Nfuse != True:
    NAccterm1   = globalvars.Optimal[4][overlap - 1]
    NAccterm2   = NAccterm1.get_parent()
    NAccterm    = NAccterm2.id

    CAccterm1   = globalvars.Optimal[4][overlap]
    CAccterm2   = CAccterm1.get_parent()
    CAccterm    = CAccterm2.id

    NDonterm1   = globalvars.Optimal[3][0]
    NDonterm2   = NDonterm1.get_parent()
    NDonterm    = NDonterm2.id

    CDonterm1   = globalvars.Optimal[3][overlap * 2 - 1]
    CDonterm2   = CDonterm1.get_parent()
    CDonterm    = CDonterm2.id

if Nfuse == True:
    ######################
    #######################
    #Find the last residue for donor and first for acceptor in list
    CAccterm1 = globalvars.Optimal[4][overlap - 1]
    CAccterm2 = CAccterm1.get_parent()
    CAccterm  = CAccterm2.id

    NAccterm2    = CAccterm2
    NAccterm     = CAccterm

    NDonterm1   = globalvars.Optimal[3][0]
    NDonterm2   = NDonterm1.get_parent()
    NDonterm    = NDonterm2.id

    print "yayayaya", globalvars.Optimal[0].get_list()[0].get_list()[0].get_list()[-1], globalvars.Optimal[0].get_list()[0].get_list()[0].get_list()[-1].get_full_id()

    CDonterm2   = globalvars.Optimal[0].get_list()[0].get_list()[0].get_list()[-1]
    CDonterm    = CDonterm2.get_full_id()[3]

    print "blah", CDonterm, globalvars.Optimal[3]

TruncDonN = (NDonterm[1]+overlap)
TruncDonC = (CDonterm[1]-overlap)

#Prints the N/C terminals of Acceptor/Donor for insertion
print " "
print "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~Best Tag~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
print "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
print " RMSD match of N/C overlap = ", globalvars.Optimal[1], "                                          "
print "                                                                                       "
print " Acceptor = ", Acceptor, "                                                             "
print " Acceptor optimal N terminal start residue is", NAccterm2, "                           "
print " Acceptor optimal C terminal end residue is", CAccterm2, "                             "
print "                                                                                       "
print " Donor = ", globalvars.Optimal[0], "                                                              "
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
log.write(" RMSD match of N/C overlap = " + str(globalvars.Optimal[1]) + "                             "+"\n")
log.write("                                                                                 "+"\n")
log.write(" Acceptor = " + str(Acceptor) + "                                                "+"\n")
log.write(" Acceptor optimal N terminal start residue is" + str(NAccterm2) + "              "+"\n")
log.write(" Acceptor optimal C terminal end residue is"+ str(CAccterm2) + "                 "+"\n")
log.write("                                                                                 "+"\n")
log.write(" Donor = " + str(globalvars.Optimal[0]) + "                                                 "+"\n")
log.write(" Donor optimal N terminal start residue is" + str(NDonterm2) + "                 "+"\n")
log.write(" Donor optimal C terminal end residue is" + str(CDonterm2) + "                   "+"\n")
log.write("                                                                                 "+"\n")
log.write(" Protagger suggests truncating Donor to residue numbers:                         "+"\n")
log.write(" N terminal = "+  str((NDonterm[1]+overlap)) + "                                 "+"\n")
log.write(" C terminal = "+ str((CDonterm[1]-overlap)) + "                                  "+"\n")
log.write("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"+"\n")
log.write("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"+"\n")
log.write("\n")
log.write('Total time taken (seconds) = '+ str(time.time() - globalvars.Startime) + "\n")
log.close()



#Uses the above objects with unwanted residue ids to pull out all residue ids to be deleted
NdeleteAcc = range(NAccterm[1] + 1, CAccterm[1])
CdeleteAcc = range(globalvars.AccepTup[1].id[1] -1, CAccterm[1])
NdeleteDon = range(Besttag[0][1][0] - 1, NDonterm[1])
CdeleteDon = range(CDonterm[1] + 1, Besttag[0][1][1] + 1)

#for model in Acceptor:
#    for chain in model:
#        for residue in chain:
#            print residue

#print "tester 5", NdeleteAcc, CdeleteAcc, NAccterm[1] + 1, CAccterm[1]
#print "This is a test", NdeleteDon, CdeleteDon
#Chews back residues of the Acceptor/Donor to provide structure files with optimal residues
PrunedDonor = prunePDB(globalvars.Optimal[0], NdeleteDon, CdeleteDon)
PrunedAcceptor = prunePDB(Acceptor, NdeleteAcc, CdeleteAcc)
print 'Optimal acceptor/donor pruned'

#for model in Acceptor:
#    for chain in model:
#        for residue in chain:
#            print residue

OptimalAcceptorfinal = AccDonsuperimposer(PrunedAcceptor, globalvars.Optimal[4], globalvars.Optimal[3])
OptimalDonorfinal = AccDonsuperimposer(PrunedDonor, globalvars.Optimal[4], globalvars.Optimal[3])

AccDonIter = []
AccDonIter.append(globalvars.Optimal[4])
AccDonIter.append(globalvars.Optimal[3])
AccDonIter.append(OptimalDonorfinal)
AccDonIter = pickle.dumps(AccDonIter)

Begin = time.time()

#Pulls name of optimal donor tag
Donametmp = re.search("'....'", str(globalvars.Optimal[0]))
Doname = (Donametmp.group(0)).strip("'")

#print "testoroony", OptimalAcceptorfinal.get_list()[0].get_list()

#This section appends a tag to each of the monomers in the multimer, saves the output pdb
for chains in OptimalAcceptorfinal.get_list()[0].get_list():
    if globalvars.chainregx != []:
        for regex in globalvars.chainregx:
            if re.search('['+regex+']', chains.id):
                AcceptorDonorIter2 = pickle.loads(AccDonIter)
                AcceptorIter2 = AcceptorDonorIter2[0]
                DonorIter2 = AcceptorDonorIter2[1]
                TagModel = AcceptorDonorIter2[2][0].get_list()[0]

                Donors = []
                Acceptors = []
                #print "yet another test", AcceptorIter2
                #Pulls out the specific N/Cterminal Acceptor/Donor CA atom lists for use with alignment
                for residue in TagModel:
                    for atom in residue:
                        for donor in DonorIter2:
                            if atom.get_full_id()[3][1] == donor.get_full_id()[3][1] and atom.get_id() == 'CA':
                                Donors.append(atom)

                for residue in chains.get_list():
                    for atom in residue:
                        for acceptor in AcceptorIter2:
                            #print "woah", acceptor
                            #print "yeye", atom.get_full_id()[3][1], acceptor.get_full_id()[3][1]
                            if atom.get_full_id()[3][1] == acceptor.get_full_id()[3][1] and atom.get_id() == 'CA':
                                Acceptors.append(residue['CA'])

                #Using the pulled out lists, superimposer aligns and applies transformation to atoms in tag chain
               # print "another test", Acceptors, Donors
                super_imposerx = None
                super_imposerx = Bio.PDB.Superimposer()
                super_imposerx.set_atoms(Acceptors, Donors)
                super_imposerx.apply(TagModel.get_atoms())

                #Saves the new tag chain
                io=Bio.PDB.PDBIO()
                io.set_structure(TagModel)
                io.save((outfolder+str(pdb_acceptor)+'_'+Doname+'_'+str(chains.get_id())+'tag.pdb').replace(" ", ""))

                print str(chains.get_id()),'tag.pdb saved'

    elif globalvars.chainregx == []:
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
        io.save((outfolder+str(pdb_acceptor)+'_'+Doname+'_'+str(chains.get_id())+'tag.pdb').replace(" ", ""))

        print str(chains.get_id()),'tag.pdb saved'

#Saves the PDB of the optimal final acceptor
io=Bio.PDB.PDBIO()
io.set_structure(OptimalAcceptorfinal)
io.save(outfolder+pdb_acceptor+'_optimised_for_'+Doname+'.pdb')

print 'Optimised', pdb_acceptor, '.pdb saved'
print ' '
print 'Optimal acceptor/donor pairs generated and saved'
print 'Protagger is done!'
print 'Time taken (total) =', time.time() - globalvars.Startime
