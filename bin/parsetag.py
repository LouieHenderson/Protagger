from Bio import PDB
import re
import time

from flexscrubber import Flexscrubber
from HETscrubber import HETscrubber

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
