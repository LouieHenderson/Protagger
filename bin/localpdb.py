from Bio import PDB

def localpdb(pdb_acceptor):

    parser = PDB.PDBParser(PERMISSIVE=1)

    #Creates structure object of protein to be tagged
    acceptor = parser.get_structure(pdb_acceptor,"./"+pdb_acceptor+".pdb")

    return acceptor
