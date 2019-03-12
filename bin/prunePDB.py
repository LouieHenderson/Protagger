def prunePDB(PDBstructure, Nprune, Cprune):

    #Chews back N/C terminal ends of PDB based on N/Cprune input, returns new structure
    for model in PDBstructure:
        for chain in model:
            for residue in chain:
                x = residue.id
                for x in Nprune:
                    if (' ', x, ' ') in chain:
                        chain.detach_child((' ', x, ' '))
                for x in Cprune:
                    if (' ', x, ' ') in chain:
                        chain.detach_child((' ', x, ' '))

    return PDBstructure
