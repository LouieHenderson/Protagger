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
    #print 'Tag pruned:'
    #if nreslist != []:
    #    print 'N-terminal chewed by', len(nreslist), 'residues to', nreslist[len(nreslist)-1]
    #if creslist != []:
    #    print 'C-terminal chewed by', len(creslist), 'residues to', creslist[len(creslist)-1]
    #if nreslist == [] and creslist == []:
    #    print 'No change'


    for residue in nreslist:
        residue.get_parent().detach_child(residue.id)
    for residue in creslist:
        residue.get_parent().detach_child(residue.id)

    return PDBstructure
