import re
import globalvars
import multiprocessing as mp
from multiprocessing import Pool, Pipe, Queue, Process
import time
from restruecheck import restruecheck

def chunks(list, n):
    for i in range(0, len(list), n):
        yield list[i:i+n]

def findClash(listofchains, start_tagged_residue, end_tagged_residue, spherelist, start_potentialTag, end_potentialTag, overlap, sphereclash):

    #Provides list objects which are appended to during the main clash checking script
    acceptorchains = []
    globalvars.donorchains = []
    spherex = []

    #Separates acceptor chains from donor chains into respective objects
    for chain in listofchains.get_list()[0].get_list():
        if (re.findall('(?<=.)tag', str(chain.get_id()))) == ['tag']:
            globalvars.donorchains.append(chain)
        else:
            acceptorchains.append(chain)

    #Provides a checkpoint object which acts as gatekeeper, if any clash is found, all subsequent checks abandoned
    #Compares every atom within every residue of the N/C tag sites with all atoms in the tags minus NC truncation
    #Excludes clash check of region between N/C terminal of Acceptor that will be excised

    globalvars.TrueCheck = [True, True]

    #Iterates through acceptor/donor pairs, provides rough spherical approximations
    #which indicate overlap to identify which donors to clash check in finer grain
    for chain in (acceptorchains + globalvars.donorchains):
        #print "chain", chain.get_id()
        if globalvars.TrueCheck[0] == True and globalvars.TrueCheck[1] == True:
            spherex = []
            checkspheres = []
            checkdonors = []
            spectag = []
            for spheres in sphereclash:
                if spheres[0][2] == chain.get_id():
                    spherex.append(spheres)
            for spheres in spherex:
                if spheres[2] == False or spheres[1][2] == (str(chain.get_id())+'tag') and spheres[1][2] != (str(chain.get_id())):
                    checkspheres.append(spheres)
            for donors in globalvars.donorchains:
                #print "this is the donor", donors
                for sphere in checkspheres:
                    if donors.get_id() == sphere[1][2]:
                    #    print "donor", donors.get_id(), "chain", chain.get_id()
                        checkdonors.append(donors)

            #print "These are the lists \n"
            #print "sphereclash", sphereclash, "\n"
            #print "checkspheres", checkspheres, "\n"
            #print "spherex", spherex, "\n"
            #print "globalvars.donorchains", globalvars.donorchains, "\n"

            #Chunk object created to start multiprocessing child processes
            reschunks = list(chunks(chain.get_list(), (len(chain.get_list()) / globalvars.CPUs)+1))

            #Multiprocessing block which generates fine grain clash check processes.
            #Uses TrueCheck object booleans to track the identification of a False
            #return of clash check in any of the child processes
            if __name__ == 'findClash':
                queue1 = Queue()
                for chunk in reschunks:
                    process = Process(target=restruecheck, args=(start_tagged_residue, end_tagged_residue, start_potentialTag, end_potentialTag, overlap, globalvars.threshold_A, checkdonors, queue1, chunk))
                    process.Daemon = True
                    process.start()
                for chunk in reschunks:
                    queue1.put(globalvars.TrueCheck)
                    while len(mp.active_children()) > 1 and globalvars.TrueCheck[0] != False and globalvars.TrueCheck[1] != False:
                        globalvars.TrueCheck = queue1.get()
                        if globalvars.TrueCheck[0] == False or globalvars.TrueCheck[1] == False:
                            process.join()
                            break
                        queue1.put(globalvars.TrueCheck)

            #Control steps to prevent runaway process generation and to ensure correct process/queue closure
            while len(mp.active_children()) > 1:
                time.sleep(0.1)

            process.join()
            time.sleep(0.1)
            queue1.close()


    #print 'Final TrueCheck being returned = ', TrueCheck

    return globalvars.TrueCheck
