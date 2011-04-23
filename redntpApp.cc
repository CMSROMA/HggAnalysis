#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <string>

#include <TTree.h>
#include <TChain.h>
#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TH2.h>
#include <TH1F.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TVector2.h>
#include <TVector3.h>
#include <TAxis.h>
#include <TString.h>

#include "RedNtpTree.h"
#include "IsGJet.h"
#include "CrossSection.h"

using namespace std;

int main(int argc, char* argv[]) {

      //================ Parameters 
      if(argc < 3 ) {
        cout << "Usage:  ./tmp/redntpApp  listfile   outputfile   selection \n" 
             << "    listfile:    list of root files incusing protocol eg dcap:/// .....\n"
             << "    outputfile:  name of output root file  eg output.root\n"
             << "    selection:   selection for preselecting events"  
             << endl;
        exit(-1);
      }


      //  1st option: nome del file contenete lista di root file
     
      // Input list
      char listName[500];
      sprintf(listName,argv[1]); 

      // Output filename (.root)  
      TString OutputFileName(argv[2]);
      
      // Name of input tree objects in (.root) files 
      char treeName[100] = "myanalysis/pippo";
      //sprintf(treeName,argv[2]);

      // fai TChain
      TChain *chain = new TChain(treeName);
      char pName[500];
      ifstream is(listName);
      if(! is.good()) {
         cout << "int main() >> ERROR : file " << listName << " not read" << endl;
         is.close();
         exit(-1);
      }
      cout << "Reading list : " << listName << " ......." << endl;
  
      while( is.getline(pName, 500, '\n') ) {
	 if (pName[0] == '#') continue;
	   //cout << "   Add: " << pName << endl;
	   chain->Add(pName); 
      }
      is.close();


      //4th option:  name of flat file with cuts
      char  selection[100];
      sprintf(selection,argv[3]);
      string finder(selection);
      if(finder == "") sprintf(selection,"looseeg");
      cout << "Photon selection is : " << selection << endl;

       // find cross section for this list
       float myxsec = CrossSection(listName);

       // filter for 2gam + jets. this is included in GJets samples but we use dedicated DiPhotonjets-madgraph
       int isGJet = IsGJet(listName);

       // compute equivalent luminosity
       //Long64_t  ntot = chain->GetEntries();
       Long64_t  ntot = 1;
       double lumi = ntot/myxsec;
       cout << "#events: " << ntot << "\t xsec: " << myxsec << " pb\t equiv. lumi: " 
            << lumi/1000. << " fb-1"
            << endl;

       // run analysis code
       RedNtpTree tool(chain, OutputFileName);
       tool.SetNtotXsection( ntot, myxsec );
       tool.Loop(isGJet, selection);
}
