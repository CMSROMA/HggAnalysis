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

#include "AnalysisTool.h"
#include "IsGJet.h"

using namespace std;

int main(int argc, char* argv[]) {

      //================ Parameters 
      if(argc < 3 ) {
        cout << "Usage:  ./tmp/analysis  listfile   outputfile  [cutfile]\n" 
             << "    listfile:    list of root files incusing protocol eg dcap:/// .....\n"
             << "    outputfile:  name of output root file  eg output.root\n"
             //<< "    histo_color: fill color for histograms\n"
             << "optional:\n"
             << "    cutfile:     flat file with list of cuts in correct order"  
             << endl;
        exit(-1);
      }


      //  1st option: nome del file contenete lista di root file
     
      // Input list
      char listName[500];
      sprintf(listName,argv[1]); 

      // Output filename (.root)  
      TString OutputFileName(argv[2]);
      //TString OutputFileName("out_analysis.root");
      
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

      //3rd option:  name of flat file with cuts
      char  cutfile[200];
      sprintf(cutfile,"cuts.txt");
      if( argc >3 ) {
        sprintf(cutfile,argv[3]);
        cout << "cuts to be read from file: " << cutfile << endl;
      }

      //4th option: hfill color   OBSOLETE
      int color = 1;
      //color = atoi( argv[3] );
      //cout << "Hai scelto color: " << color << endl;

       //================ Run analysis
 
       AnalysisTool tool(chain,OutputFileName);
       tool.SetHFillColor( color );

       // filter out 2gamma+jets events when running on GJets sample
       if( IsGJet(listName) > 0 ) tool.Filter2GammaJet();

       tool.ReadCutsFromFile(cutfile);
       tool.Loop();
}
