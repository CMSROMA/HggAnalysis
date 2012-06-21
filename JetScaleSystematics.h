#ifndef __JETSCALESYSTEMATICS__
#define __JETSCALESYSTEMATICS__

#include <string>
#include <TString.h>
#include <vector>
#include <iostream>
#include <map>


// ------------------------------------------------------------------------------------
class JetScaleSystematics {
public:
  
  JetScaleSystematics(TString inputfile);
  float getJESUncertainty(double eta, double pt);

private:
  float mineta[38];
  float minpt[38][39];
  float unc[38][39];

};


#endif
