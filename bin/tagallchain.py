import Bio.PDB
import cPickle as pickle

def tagallchain(taggedProtein, CompositePDB1, chunk, AccDonIter, queue):
    TagChain1output = []

    for chains in chunk:
        AcceptorDonorIter2 = pickle.loads(AccDonIter)
        AcceptorIter2 = AcceptorDonorIter2[0]
        DonorIter2 = AcceptorDonorIter2[1]
        TagChain1 = AcceptorDonorIter2[2]

        #Creates a superimposer, aligns atoms in residue of tag chain A with current chain to get translational matrix, applies this to Acceptor list
        Donors = []
        Acceptors = []

        #Pulls out the specific N/Cterminal Acceptor/Donor CA atom lists for use with alignment
        for residue in TagChain1.get_list():
            for atom in residue:
                for donor in DonorIter2:
                    if atom.get_full_id()[3][1] == donor.get_full_id()[3][1] and atom.get_id() == 'CA':
                        Donors.append(residue['CA'])

        for residue in chains.get_list():
            for atom in residue:
                for acceptor in AcceptorIter2:
                    if atom.get_full_id()[3][1] == acceptor.get_full_id()[3][1] and atom.get_full_id()[4][0] == 'CA':
                        Acceptors.append(residue['CA'])

        #Create superimposer object for transformations
        super_imposerx = Bio.PDB.Superimposer()

        #Using newly transformed Acceptor list, aligns Donor list and applies this transformation to instance of TagChain
        super_imposerx.set_atoms(Acceptors, Donors)
        super_imposerx.apply(TagChain1)
        TagChain1.id = str(chains.id)+str('tag')
        TagChain1.detach_parent()

        TagChain1output.append(TagChain1)

    queue.put(TagChain1output)
