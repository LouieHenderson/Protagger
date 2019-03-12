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
