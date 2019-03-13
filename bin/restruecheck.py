import globalvars

def restruecheck(start_tagged_residue, end_tagged_residue, start_potentialTag, end_potentialTag, overlap, threshold_A, checkdonors, queue1, residlist):

    #Iterate through each atom in residues of compared acceptor/donor pairs and checks for potential clash
    for resid1 in residlist:
        #print "this is resid parent", resid1.get_parent().get_parent().get_parent()
        Checker = queue1.get()
        queue1.put(Checker)
        if Checker[0] == False or Checker[1] == False:
            globalvars.TrueCheck = Checker
            queue1.put(globalvars.TrueCheck)
        CheckpointA = resid1.get_id()[1] <= (start_tagged_residue - 1 - overlap) or resid1.get_id()[1] >= (end_tagged_residue + 1 + overlap)
        if CheckpointA == True:
            if globalvars.TrueCheck[0] == True:
                for atom1 in resid1:
                    for chain2 in checkdonors:
                        print "this is being checked", chain2, resid1.get_parent().get_parent().get_parent()
                        for resid2 in chain2:
                            CheckpointB = resid2.get_id()[1] >= (start_potentialTag + 1 + overlap) and resid2.get_id()[1] <= (end_potentialTag - 1 - overlap)
                            if CheckpointB == True and globalvars.TrueCheck[0] == True:
                                for atom2 in resid2:
                                    Atomdist = atom1 - atom2
                                    Boolean1 = Atomdist >= globalvars.threshold_A
                                #    print Atomdist, Boolean1
                                    if Boolean1 == True:
                                        globalvars.TrueCheck[0] = True
                                    elif Boolean1 == False:
                                        globalvars.TrueCheck[0] = False
                                        print "False was found"
                                        queue1.put(globalvars.TrueCheck)
                                        break
                                    else:
                                        pass
                            else:
                                pass
            elif globalvars.TrueCheck[0] == False:
                queue1.put(globalvars.TrueCheck)
                break
        else:
            pass
