//
//     PDB_MERGE merge two PDB files
//     Copyright (C) 2005 Martyn Winn
//
//     This code is distributed under the terms and conditions of the
//     CCP4 Program Suite Licence Agreement as a CCP4 Application.
//     A copy of the CCP4 licence can be obtained by writing to the
//     CCP4 Secretary, Daresbury Laboratory, Warrington WA4 4AD, UK.
//
// =========================================================================

// TODO
// * Output combined files as separate models e.g. to be used as
//   NMR ensemble in Molrep

#include <string.h>
#include <stdlib.h>

#include "mmdb/mmdb_manager.h"
#include "ccp4/ccp4_parser.h"
#include "ccp4/ccp4_general.h"
#include "ccp4/ccp4_program.h"

using namespace CCP4;

void CheckFileInput(CMMDBManager & MMDBManager, const char *logname, int & RC) {
int       lcount;
char      S[500], errstr[500];

  // Check for possible errors:
  if (RC) {
    //  An error was encountered. MMDB provides an error messenger
    //  function for easy error message printing.
    printf ( " ***** ERROR #%i READ:\n\n %s\n\n",RC,GetErrorDescription(RC) );
    //  Location of the error may be identified as precise as line
    //  number and the line itself (PDB only. Errors in mmCIF are
    //  located by category/item name. Errors of reading BINary files
    //  are not locatable and the files are not editable). This
    //  information is now retrieved from MMDB input buffer:
    MMDBManager.GetInputBuffer ( S,lcount );
    if (lcount>=0) 
      printf ( "       LINE #%i:\n%s\n\n",lcount,S );
    else if (lcount==-1)
      printf ( "       CIF ITEM: %s\n\n",S );
    //  now quit
    sprintf(errstr,"Error reading %s",logname);
    ccperror ( 1, errstr );
  } else  {
    //  MMDB allows to identify the type of file that has been just
    //  read:
    switch (MMDBManager.GetFileType())  {
      case MMDB_FILE_PDB    : printf ( " PDB"         );  break;
      case MMDB_FILE_CIF    : printf ( " mmCIF"       );  break;
      case MMDB_FILE_Binary : printf ( " MMDB binary" );  break;
      default : printf ( " Unknown format (report as a bug!)" );
    }
    printf ( " file %s has been read in.\n",getenv(logname) );
  }
}

char GetChainID(const PCChain chain) {
  return chain->GetChainID()[0] ? chain->GetChainID()[0] : ' ';
}

int MakeNewChainIDs(CMMDBManager & MMDB, int & nch1, char *chainIDs) {

// adjust the chain names in MMDB so that they are different from those
// listed in chainIDs

int            ich,i,natom,nmod,imod,nch=0,nchain,ichain;
int            got_duplicate, changed_chain_names = 0;
PPCAtom        atomTable;
char           *allChainIDs;
char           chID,new_chain_id,new_chain_id_str[2];
PCChain        chain,new_chain;
PCResidue      res;

  // get list of all chain IDs in MMDB - needed for generating new
  // chain IDs
  nmod = MMDB.GetNumberOfModels();
  for (imod=1;imod<=nmod;imod++)
    nch += MMDB.GetNumberOfChains(imod);
  allChainIDs = new char[nch+nch1+1];
  ich = 0;
  for (imod=1;imod<=nmod;imod++) {
    nchain = MMDB.GetNumberOfChains(imod);
    for (ichain=0;ichain<nchain;ichain++) {
      chain = MMDB.GetChain(imod,ichain);
      allChainIDs[ich++] = GetChainID(chain);
    }
  }
  allChainIDs[nch] = '\0';

  // loop over input chain IDs
  for (ich=0;ich<nch1;ich++) {
    chID = chainIDs[ich];
    printf(" Checking for duplication of chain %c .... \n", chID);

    got_duplicate = 0;

    for (imod=1;imod<=nmod;imod++) {
     nchain = MMDB.GetNumberOfChains(imod);
     for (ichain=0;ichain<nchain;ichain++) {
      chain = MMDB.GetChain(imod,ichain);
      if (GetChainID(chain) == chID) {

       if (!got_duplicate) {
        // get first unused chainID
        new_chain_id = 'A';
        while (strchr(allChainIDs,new_chain_id) || 
               strchr(chainIDs,new_chain_id)) new_chain_id = char(int(new_chain_id)+1);
        allChainIDs[nch++] = new_chain_id;
        allChainIDs[nch] = '\0';

        got_duplicate = 1;
        changed_chain_names = 1;
       }

       printf ( " ... replacing duplicated chain ID %c by new_chain_id %c \n", 
                      chID, new_chain_id );
       new_chain_id_str[0] = new_chain_id;
       new_chain_id_str[1] = '\0';
       chain->SetChainID(new_chain_id_str);
      }
     }
    }
  }
  delete[] allChainIDs;

  return (changed_chain_names);
}

int main ( int argc, char ** argv, char ** env )  {
CMMDBManager   MMDBManager,MMDBManager2;
int            i,kout,RC,natoms1,nmod1,nmod2,nch1=0,natoms2,lnewIDs;
int            imod,ich,nchain,ichain;
int            merge_chains=1;
char           *chainIDs;
PCModel        model1;
PCChain        chain;

// input parser parameters
int           ntok=0;
char          line[201],*key;
CCP4PARSERTOKEN * token=NULL;
CCP4PARSERARRAY * parser;

  //  1.  General CCP4 initializations
  ccp4fyp         ( argc,argv );
  ccp4ProgramName ( "PDB_MERGE"   );
  ccp4_banner();

  //  2.  Make routine initializations, which must always be done
  //      before working with MMDB
  InitMatType();
  kout = MMDB_FILE_PDB;

  // 3. Keyworded input

  /* Initialise a parser array used by cparser
     This is used to return the tokens and associated info
     Set maximum number of tokens per line to 20 */
  parser = (CCP4PARSERARRAY *) ccp4_parse_start(20);

  if (parser == NULL) ccperror ( 1,"Couldn't create parser array" );
  
  /* Set some convenient pointers to members of the parser array */
  key   = parser->keyword;
  token = parser->token;

  /* Read lines from stdin until END/end keyword is entered or
     EOF is reached */
  RC   = 0;

  while (!RC) {

    /* Always blank the line before calling cparser
       to force reading from stdin */
    line[0] = '\0';

    /* Call cparser to read input line and break into tokens
       Returns the number of tokens, or zero for eof */
    ntok = ccp4_parser(line,200,parser,1);

    if (ntok < 1) {

      /* End of file encountered */
      RC = 111;

    } else {      

      /* Keyword interpretation starts here */
      if (ccp4_keymatch("END",key))  {

        /* End of keyworded input */
        RC = 111;

      } else if (ccp4_keymatch("OUTPUT",key))  {
	if (ntok >= 2) {
          if (!strncmp(UpperCase(token[1].fullstring),"PDB",3)) kout = MMDB_FILE_PDB;
          if (!strncmp(UpperCase(token[1].fullstring),"CIF",3)) kout = MMDB_FILE_CIF;
          if (!strncmp(UpperCase(token[1].fullstring),"BIN",3)) kout = MMDB_FILE_Binary;
	}

      } else if (ccp4_keymatch("NOMERGE",key))  {
        merge_chains=0;
      }
    }
  }

  if (RC==111)  RC = 0;  // normal return from the parser loop

  /* Clean up parser array */
  ccp4_parse_end ( parser );

  // Read first coordinate file by its logical name
  RC = MMDBManager.ReadCoorFile1 ( "XYZIN1" );
  CheckFileInput(MMDBManager, "XYZIN1", RC);

  natoms1 = MMDBManager.GetNumberOfAtoms();
  printf ( " Read in %d atoms from XYZIN1\n", natoms1 );

  // If we are going to keep chains from the 2 files separate, then
  // we need to collect information on chains in file 1
  if (!merge_chains) { 
    // get list of chain IDs in XYZIN1
    nmod1 = MMDBManager.GetNumberOfModels();
    for (imod=1;imod<=nmod1;imod++)
      nch1 += MMDBManager.GetNumberOfChains(imod);
    chainIDs = new char[nch1+1];
    ich = 0;
    for (imod=1;imod<=nmod1;imod++) {
      nchain = MMDBManager.GetNumberOfChains(imod);
      for (ichain=0;ichain<nchain;ichain++) {
        chain = MMDBManager.GetChain(imod,ichain);
        chainIDs[ich++] = GetChainID(chain);
      }
    }
    chainIDs[nch1] = '\0';
  }

  // Read and add second coordinate file by its logical name
  // N.B. can't use ReadCoorFile1 again because it resets the manager
  if (merge_chains) {
    RC = MMDBManager.AddPDBASCII1 ( "XYZIN2" );
    CheckFileInput(MMDBManager, "XYZIN2", RC);

    natoms2 = MMDBManager.GetNumberOfAtoms();
    printf ( " Read in %d atoms from XYZIN2\n", natoms2 - natoms1 );
 
    printf ( "\n XYZIN1 and XYZIN2 will be merged together preserving residue and chain labels.\n\n" );

 } else {
    RC = MMDBManager2.ReadCoorFile1 ( "XYZIN2" );
    CheckFileInput(MMDBManager2, "XYZIN2", RC);

    natoms2 = MMDBManager2.GetNumberOfAtoms();
    printf ( " Read in %d atoms from XYZIN2\n", natoms2 );
 
    printf ( "\n XYZIN2 will be appended to XYZIN1 keeping constituent chains separate.\n\n" );

    // change change IDs if necessary so don't clash with XYZIN1
    lnewIDs = MakeNewChainIDs(MMDBManager2, nch1, chainIDs );
    if (lnewIDs)
      printf ( " Some chains from XYZIN2 renamed.\n");

    // now copy chains to first model of MMDBManager
    // (this requires more work when multiple models)
    model1 = MMDBManager.GetModel ( 1 ); 

    nmod2 = MMDBManager2.GetNumberOfModels();
    for (imod=1;imod<=nmod2;imod++) {
      nchain = MMDBManager2.GetNumberOfChains(imod);
      for (ichain=0;ichain<nchain;ichain++) {
        chain = MMDBManager2.GetChain(imod,ichain);

        // now add chain to first model of first structure
        model1->AddChain ( chain );
      }
    }
  }

  // Tidy up output
  MMDBManager.PDBCleanup ( PDBCLEAN_INDEX );
  MMDBManager.PDBCleanup ( PDBCLEAN_SERIAL );

  //  5.  Now output the whole MMDB in the required format:
  switch (kout)  {
    default               :
    case MMDB_FILE_PDB    : 
        printf ( " Writing PDB file %s\n",getenv("XYZOUT") );
        RC = MMDBManager.WritePDBASCII1 ( "XYZOUT" );
      break;
    case MMDB_FILE_CIF    :
        printf ( " Writing CIF file %s\n",getenv("XYZOUT") );
        RC = MMDBManager.WriteCIFASCII1 ( "XYZOUT" );
      break;
    case MMDB_FILE_Binary :
        printf ( " Writing MMDB binary file %s\n",getenv("XYZOUT") );
        RC = MMDBManager.WriteMMDBF1 ("XYZOUT"  );
  }

  if (!merge_chains) delete[] chainIDs;
  ccperror ( 0,"Normal termination" );

  return 0;

}
