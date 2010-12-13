#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <string>
#include <map>

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
#include <TRegexp.h>

#include "RedNtpTree.h"

using namespace std;

float CrossSection(const char*);

int main(int argc, char* argv[]) {

      //================ Parameters 
      if(argc < 3 ) {
        cout << "Usage:  ./tmp/redntpApp  listfile   outputfile   histo_color cutfile\n" 
             << "    listfile:    list of root files incusing protocol eg dcap:/// .....\n"
             << "    outputfile:  name of output root file  eg output.root\n"
             //<< "    cutfile:     flat file with list of cuts in correct order"  
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

/*
      //4th option:  name of flat file with cuts
      char  cutfile[200];
      sprintf(cutfile,argv[3]);
      cout << "cuts to be read from file: " << cutfile << endl;
*/

       // find cross section for this list
       float myxsec = CrossSection(listName);

       TString alist(listName);
       // filter for 2gam + jets. this is included in GJets samples but we use dedicated DiPhotonjets-madgraph
       int isGJet = 0;
       int pos = alist.Index( TRegexp("GJet") );
       if(pos>=0) {
          isGJet = 1;
          cout << "GJet* samples. will  filter out 2g+jet events included in dedicated DiPhotonJets-madgraph" << endl;
       }

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
       tool.Loop(isGJet);
}

float CrossSection(const char* sample) {
       TString alist(sample);
       //  cross sections
       std::map<TString, double> xsec; // in pb as in https://twiki.cern.ch/twiki/bin/view/CMS/ProductionFall2010
       xsec["GJet_Pt-20_doubleEMEnriched_TuneZ2_7TeV-pythia6"]       = 493.44;
       xsec["DiPhotonJets_7TeV-madgraph"] 			     = 134.;
       xsec["DiPhotonBox_Pt10to25_TrackingParticles_7TeV-pythia6"]   = 358.2;
       xsec["DiPhotonBox_Pt25to250_TrackingParticles_7TeV-pythia6"]  = 12.37;
       xsec["DiPhotonBox_Pt250toinf_TrackingParticles_7TeV-pythia6"] = 2.08e-4;
       xsec["DYJetsToLL_TuneZ2_M-50_7TeV-madgraph"]                  = 2321.;
       xsec["QCD_Pt-30to40_doubleEMEnriched_TuneZ2_7TeV-pythia6"]    = 9.61e3;
       xsec["QCD_Pt-40_doubleEMEnriched_TuneZ2_7TeV-pythia6"]        = 4.04e4;
       xsec["run2010"]                                               = 1.;

       xsec["GluGluToHToGG_M-90_7TeV-powheg-pythia6"] = 1.;
       xsec["GluGluToHToGG_M-100_7TeV-powheg-pythia6"] = 1.;
       xsec["GluGluToHToGG_M-105_7TeV-powheg-pythia6"] = 1.;
       xsec["GluGluToHToGG_M-110_7TeV-powheg-pythia6"] = 1.;
       xsec["GluGluToHToGG_M-115_7TeV-powheg-pythia6"] = 1.;
       xsec["GluGluToHToGG_M-120_7TeV-powheg-pythia6"] = 1.;
       xsec["GluGluToHToGG_M-130_7TeV-powheg-pythia6"] = 1.;
       xsec["GluGluToHToGG_M-140_7TeV-powheg-pythia6"] = 1.;
       xsec["VBF_HToGG_M-90_7TeV-powheg-pythia6"] = 1.;
       xsec["VBF_HToGG_M-95_7TeV-powheg-pythia6"] = 1.;
       xsec["VBF_HToGG_M-100_7TeV-powheg-pythia6"] = 1.;
       xsec["VBF_HToGG_M-105_7TeV-powheg-pythia6"] = 1.;
       xsec["VBF_HToGG_M-110_7TeV-powheg-pythia6"] = 1.;
       xsec["VBF_HToGG_M-115_7TeV-powheg-pythia6"] = 1.;
       xsec["VBF_HToGG_M-120_7TeV-powheg-pythia6"] = 1.;
       xsec["VBF_HToGG_M-130_7TeV-powheg-pythia6"] = 1.;
       xsec["VBF_HToGG_M-140_7TeV-powheg-pythia6"] = 1.;
       xsec["WH_ZH_TTH_HToGG_M-90_7TeV-pythia6"] = 1.;
       xsec["WH_ZH_TTH_HToGG_M-95_7TeV-pythia6"] = 1.;
       xsec["WH_ZH_TTH_HToGG_M-100_7TeV-pythia6"] = 1.;
       xsec["WH_ZH_TTH_HToGG_M-105_7TeV-pythia6"] = 1.;
       xsec["WH_ZH_TTH_HToGG_M-110_7TeV-pythia6"] = 1.;
       xsec["WH_ZH_TTH_HToGG_M-115_7TeV-pythia6"] = 1.;
       xsec["WH_ZH_TTH_HToGG_M-120_7TeV-pythia6"] = 1.;
       xsec["WH_ZH_TTH_HToGG_M-130_7TeV-pythia6"] = 1.;
       xsec["WH_ZH_TTH_HToGG_M-140_7TeV-pythia6"] = 1.;

       //xsec[""] = ;
       //xsec[""] = ;

       vector<TString> keys;
       for(std::map<TString, double>::const_iterator k = xsec.begin(); k != xsec.end(); ++k) {
          keys.push_back( k->first );
       }
       //  cross sections end


       //cout << "input list: <" << alist << ">" << endl;
       double myxsec = -1;
       for(int i=0; i< keys.size() && myxsec<0.; ++i) {
          //cout << "key: " << keys[i] << endl;
          int pos = alist.Index( TRegexp(keys[i]) );
          if(pos>=0) {
             myxsec = xsec[keys[i]];
             cout << "xsec: " << myxsec << "\t for " << alist << endl;
          }
       } 
       if(myxsec<0) {
         cout << "No xsection found for " << alist << endl;
         cout << "exiting..." << endl; 
         exit(-1);
       } else {
         return myxsec;
       }
}
