#!/usr/bin/python

import glob
import argparse
import time
import globalvars
import Bio.PDB
from parsetag import parsetag
from getFirstLast import getFirstAndLastResiduesInChain

def databaseorder(Databaseloc, accdist):
    #Load up the database of tags to check
    taglist = glob.glob(Databaseloc)

    #Creates objects to monitor remaining tags
    tagsrem = len(taglist)
    tagstot = len(taglist)

    #Create list to identify distances
    distlist = []

    Begin = time.time()

    #Iterate through database of tags and identifies N/C terminal distances
    for donor in taglist:
        tagsrem = tagsrem - 1
        try:
            #Generates threshold criteria from Donor N/C terminals
            Donor = parsetag(donor)
            firstlast = []
            firstlast.append(getFirstAndLastResiduesInChain(Donor))
            firlasdist = firstlast[0][2]['CA']-firstlast[0][3]['CA']
            #Converts to positive int and appends tup of tag + dist to chekc list
            truedist = abs(accdist - firlasdist)
            distlist.append((donor, truedist))
        except:
            print "This tag parsed incorrectly = ", donor

    #Sorts list based on distance compared to Acceptor tag site
    sortlist = sorted((distlist), key=lambda x: x[1])
    donorlist = []

    for tup in sortlist:
        donorlist.append(tup[0])

    print "Time taken to order all tags = ", time.time() - Begin

    return donorlist
