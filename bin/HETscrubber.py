import globalvars

#print globalvars.CPUs

def HETscrubber(PDBstructure):
    #Iterates through list of amino acids, if heteroatom is found, residue is removed
    print globalvars.CPUs
    for chain in PDBstructure[0]:
        for residue in chain.get_list():
            resid = residue.id
            if resid[0] != ' ' or residue.__contains__('CA') == False or residue.get_resname() in globalvars.aminoAcids == False:
                chain.detach_child(resid)

    return PDBstructure
