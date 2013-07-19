#ifndef CrossSection_H
#define CrossSection_H

#include "TString.h"
#include <iostream>
#include <map>

float CrossSection(const char* sample) {
       using namespace std;
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
       xsec["QCD_Pt-20to30-EMEnriched_TuneZ2_7TeV-pythia6"]    = 2.45e6; 
       xsec["QCD_Pt-30to80-EMEnriched_TuneZ2_7TeV-pythia6"]    = 3.87e6; 
       xsec["QCD_Pt-80to170-EMEnriched_TuneZ2_7TeV-pythia6"]    = 1.40e5; 
       xsec["run2010"]                                               = 1.;
       xsec["run2011"]                                               = 1.;

       xsec["GluGluToHToGG_M-90_7TeV-pythia6"] = 1.;
       xsec["GluGluToHToGG_M-100_7TeV-pythia6"] = 1.;
       xsec["GluGluToHToGG_M-105_7TeV-pythia6"] = 1.;
       xsec["GluGluToHToGG_M-110_7TeV-pythia6"] = 1.;
       xsec["GluGluToHToGG_M-115_7TeV-pythia6"] = 1.;
       xsec["GluGluToHToGG_M-120_7TeV-pythia6"] = 1.;
       xsec["GluGluToHToGG_M-130_7TeV-pythia6"] = 1.;
       xsec["GluGluToHToGG_M-140_7TeV-pythia6"] = 1.;
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
         cout << "settinx xsection = 1. !" << endl;
         return 1.;
         //cout << "exiting..." << endl; 
         //exit(-1);
       } else {
         return myxsec;
       }
}
#endif
