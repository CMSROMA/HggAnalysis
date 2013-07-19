#ifndef AnalysisTool_H
#define AnalysisTool_H

//#include "tree_reader_V1.h"
//#include "tree_reader_V2.h"
//#include "tree_reader_V3.h"
#include "tree_reader_V6.h"
#include "PhotonIdCuts.h"

#include "TString.h"

class AnalysisTool : public tree_reader_V6 {

private:
  TString outname_; // output file
  int hfillcolor; // colore istogrammi

public:

  //   "Digitare 0: no selezione sui fotoni reco asso,"
  //   "Digitare 1: selezione Pt>40 Pt>20 sui fotoni reco asso,"
  //   "Digitare 2: selezione ed isolamento (non ancora implementata (3 NOV) )"

  void SetHFillColor(int col) { hfillcolor = col; }

  AnalysisTool(TTree *tree, const TString& outname);
  virtual ~AnalysisTool() { };

  void ordinamento(float *,int *, int );  

  virtual void     Loop();
  virtual std::vector <int>  firsttwo(Float_t *vec, std::vector<bool> *asso);
  bool cutID(int i, photonidcuts const& pid, std::vector<bool> *vpass = 0);

  void ReadCutsFromFile(const char* fname);

  // VBF jet cuts
  double ptJet1_cut, ptJet2_cut;
  double detaJets_cut;
  double zepp_cut;
  double mjj_cut; 

  // photon ID cuts
  photonidcuts mediumid;
  photonidcuts looseid;


  bool filter2GammaJet;
  void Filter2GammaJet() { filter2GammaJet = true; }
  


};
#endif
