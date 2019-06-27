#!/usr/bin/python

import glob
import argparse
import time
import globalvars
import Bio.PDB
from parsetag import parsetag

globalvars.initialise()

#Adds parser arguments to gain database folder
parser= argparse.ArgumentParser(description='Find the optimal tag for your protein')
parser.add_argument('-b', '--tagdatabase', type=str, help='Enter the folder containing the tags to be tested', required = True)

#Creates an object containing parsed argument names
args = parser.parse_args()
Databaseloc     = args.tagdatabase + '/*'

#Load up the database of tags to check
print "Database of tags = " + str(Databaseloc)
taglist = glob.glob(Databaseloc)

#Creates objects to monitor remaining tags
tagsrem = len(taglist)
tagstot = len(taglist)

#Create list to identify problematic tags
rmlist = []

Begin = time.time()

#Iterate through database of tags and check for errors in iteration
for donor in taglist:
    tagsrem = tagsrem - 1
    print " "
    print "Tags remaining =", tagsrem, '/', tagstot
    print "Current tag being checked =", donor
    try:
        parsetag(donor)
    except:
        print "This tag should be removed = ", donor
        rmlist.append(donor)
    print "Tag parsed successfully"

print " "
print "Time taken to parse all tags= ", time.time() - Begin
print "List of tags to remove", rmlist
