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

using namespace std;

int main(int argc, char* argv[]) {

      //================ Parameters 
      if(argc < 5 ) {
        cout << "Usage:  ./tmp/analysis  listfile   outputfile   histo_color cutfile\n" 
             << "    listfile:    list of root files incusing protocol eg dcap:/// .....\n"
             << "    outputfile:  name of output root file  eg output.root\n"
             << "    histo_color: fill color for histograms\n"
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

      //3rd option: hfill color
      int color = 1;
      color = atoi( argv[3] );
      cout << "Hai scelto color: " << color << endl;

      //4th option:  name of flat file with cuts
      char  cutfile[200];
      sprintf(cutfile,argv[4]);
      cout << "cuts to be read from file: " << cutfile << endl;

       //================ Run analysis
 
       AnalysisTool tool(chain,OutputFileName);
       tool.SetHFillColor( color );

       tool.ReadCutsFromFile(cutfile);
       tool.Loop();
}
