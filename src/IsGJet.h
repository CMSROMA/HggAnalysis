#ifndef IsGJet_HH
#define IsGJet_HH

#include "TRegexp.h"
#include "TString.h"
#include <iostream>

// filter for 2gam + jets. this is included in GJets samples but we use dedicated DiPhotonjets-madgraph
int IsGJet(const char* sample) {
   using namespace std;
   int isGJet = 0;
   TString alist(sample);
   int posgjet = alist.Index( TRegexp("GJet") );
   int posqcd = alist.Index( TRegexp("QCD") );
   if(posgjet>=0) {
      isGJet = 1;
      cout << "GJet* samples. will  filter out 2g+jet events included in dedicated DiPhotonJets-madgraph" << endl;
   }
   if(posqcd>=0) {
      isGJet = 1;
      cout << "QCD* samples. will  filter out 2g+jet events included in dedicated DiPhotonJets-madgraph" << endl;
   }
   return isGJet;
}
#endif
