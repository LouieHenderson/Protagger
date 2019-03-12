def spheredefine(Structure, AccN, AccC, DonN, DonC, overlap):
    #Center of sphere
    atomcoords = [0, 0, 0]

    #Radial maximum
    atomradmax = [0, 0, 0]

    #Vector magnitude
    vectormax = 0

    #Total number of atoms iterated through
    atomcount = 0

    #This code block iterates through atoms, adds up x, y, z coords and divides them by max number to find the mean
    for residue in Structure:
        CheckpointA = residue.get_id()[1] <= (AccN - 1 - overlap) or residue.get_id()[1] >= (AccC + 1 + overlap)
        if CheckpointA == True:
            for atom in residue:
                atomcoords[0] = atomcoords[0] + atom.get_coord()[0]
                atomcoords[1] = atomcoords[1] + atom.get_coord()[1]
                atomcoords[2] = atomcoords[2] + atom.get_coord()[2]
                atomcount = atomcount + 1
    #'print 'heyhye', atomcoords[0], atomcount

    atomcoords[0] = atomcoords[0] / atomcount
    atomcoords[1] = atomcoords[1] / atomcount
    atomcoords[2] = atomcoords[2] / atomcount

    #This code block finds the magnitude of the vectors for each atom from the centre point, compares and returns the coords of the vector with the maximum magnitude
    for residue in Structure:
        CheckpointC = residue.get_id()[1] <= (AccN - 1 - overlap) or residue.get_id()[1] >= (AccC + 1 + overlap)
        if CheckpointC == True:
            for atom in residue:
                if (((atom.get_coord()[0] - atomcoords[0])**2 + (atom.get_coord()[1] - atomcoords[1])**2 + (atom.get_coord()[2] - atomcoords[2])**2)**0.5) > vectormax:
                    vectormax = (((atom.get_coord()[0] - atomcoords[0])**2 + (atom.get_coord()[1] - atomcoords[1])**2 + (atom.get_coord()[2] - atomcoords[2])**2)**0.5)
                    atomradmax[0] = atom.get_coord()[0]
                    atomradmax[1] = atom.get_coord()[1]
                    atomradmax[2] = atom.get_coord()[2]

    return vectormax, atomcoords, Structure.get_id()
