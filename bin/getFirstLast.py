import globalvars

def getFirstAndLastResiduesInChain(theChain):
    firstResid = -1
    lastResid  = -1
    firstResidObj = None
    lastResidObj = None

    #Find the first and last residues in the PDB file
    for model in theChain[0]:
    	for chain in model:
                for tag_resid in chain:
                    if tag_resid.get_resname() in globalvars.aminoAcids:
                        if ((firstResid == -1) or (tag_resid.get_id()[1] < firstResid)):
                            firstResid = tag_resid.get_id()[1]
                            firstResidObj = tag_resid
                        if ((lastResid == -1) or (tag_resid.get_id()[1] > lastResid)):
                            lastResid = tag_resid.get_id()[1]
                            lastResidObj = tag_resid

    return (firstResid,lastResid,firstResidObj,lastResidObj)
