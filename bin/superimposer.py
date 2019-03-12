import Bio.PDB

def AccDonsuperimposer(AccDon, AccAlign, DonAlign):
    #Superimposes the optimal Acceptor/Donor
    PDBinterest = AccDon

    super_imposer = Bio.PDB.Superimposer()
    super_imposer.set_atoms(AccAlign, DonAlign)
    super_imposer.apply(PDBinterest.get_atoms())

    return PDBinterest
