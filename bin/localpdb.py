from Bio import PDB

def localpdb(pdb_acceptor):

    #Creates structure object of protein to be tagged
    parser = PDB.PDBParser(PERMISSIVE=1)
    acceptor = parser.get_structure(pdb_acceptor,"./"+pdb_acceptor+".pdb")

    return acceptor
