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
   int pos = alist.Index( TRegexp("GJet") );
   if(pos>=0) {
      isGJet = 1;
      cout << "GJet* samples. will  filter out 2g+jet events included in dedicated DiPhotonJets-madgraph" << endl;
   }
   return isGJet;
}
#endif
