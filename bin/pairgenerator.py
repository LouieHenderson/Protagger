import re

def AccPairGenerator(Acceptor, distance, nter, cter, regexid):

    NCtagsites = []

    #Generates a series of potential tag sites within the Acceptor delimited by the input N/C terminal residues
    #Uses the N/C-terminal distance of the Donor and pairs residues which match the same distance +- a wiggle factor
    for Model in Acceptor:
        for Chain in Model:
            if regexid != []:
                for regex in regexid:
                    if re.search('['+regex+']', Chain.id):
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
            else:
                for Residue1 in Chain:
                    if (Residue1.id[1] > nter) == True and (Residue1.id[1] < cter) == True:
                        for Residue2 in Chain:
                            if (Residue2.id[1] > nter) == True and (Residue2.id[1] < cter) == True:
                                DistCheck = Residue1['CA'] - Residue2['CA']
                                print "yoyo", distance, int(distance-0.5), ((DistCheck > (int(distance-0.5))) == True), ((DistCheck < (int(distance+0.5))) == True)
                                if (DistCheck > (int(distance-0.5))) == True and (DistCheck < (int(distance+0.5))) == True and (Residue2,Residue1) not in NCtagsites and (Residue1,Residue2) not in NCtagsites:
                                    NCtagsites.append((Residue1,Residue2))
                            else:
                                pass
                    else:
                        pass
                break

    return NCtagsites
