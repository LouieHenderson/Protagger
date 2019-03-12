def downloadpdb(pdb_acceptor):
    # Downloads pdb file to local folder
    folder = pdb_acceptor
    pdbl = PDBList()
    data_set = pdbl.retrieve_pdb_file(pdb_acceptor, pdir=pdb_acceptor)

    #Parses PDB file to return a structure object
    parser = PDB.PDBParser(PERMISSIVE=1)

    #Creates structure object of protein to be tagged
    acceptor = parser.get_structure(pdb_acceptor,folder+'/'+pdb_acceptor+".pdb")

    return acceptor
