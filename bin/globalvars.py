#Global objects used in optimal RMSD determination and tag structure formation
def initialise():
    global Checkpoint
    global Startime
    global Optimal
    global AccepTup
    global AcceptorIter
    global DonorIter
    global sphereclash1
    global donorchains
    global TrueCheck
    global threshold_A
    global CPUs
    global chainregx
    #global Queueregulator
    global aminoAcids


    Checkpoint     = 0
    Startime       = 0
    Optimal        = []
    AccepTup       = []
    AcceptorIter   = None
    DonorIter      = None
    sphereclash1   = []
    donorchains    = []
    TrueCheck      = [True, True]
    threshold_A    = 0
    CPUs           = 0
    chainregx      = []
    #Queueregulator = False
    aminoAcids = {'ALA':1, 'ARG':1, 'ASN':1, 'ASP':1, 'CYS':1, 'GLN':1, 'GLU':1, 'GLY':1, 'HIS':1, 'ILE':1, 'LEU':1, 'LYS':1, 'MET':1, 'PHE':1, 'PRO':1, 'SER':1, 'THR':1, 'TRP':1, 'TYR':1, 'VAL':1}
