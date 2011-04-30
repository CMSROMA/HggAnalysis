#ifndef RedNtpTree_h
#define RedNtpTree_h

//#include "higgsanal_tree_V1.h"
//#include "tree_reader_V2.h"
//#include "tree_reader_V3.h"
#include "tree_reader_V6.h"
#include "PhotonIdCuts.h"


#include <TFile.h>
#include <TString.h>
#include<vector>
#include<string>
using std::string;
using std::vector;


class RedNtpTree : public tree_reader_V6 {

public:

   RedNtpTree(TTree *tree=0, const TString& outname="redntp.root");
   virtual ~RedNtpTree();
   virtual void     Loop(int isgjet=0, char* selection = "loose");
   void SetNtotXsection(int ntot, float xsec) {
      NtotEvents = ntot;
      xsection = xsec;
      EquivLumi = ntot/xsec;
   }


private:
   TFile* hOutputFile ;
   TTree * ana_tree ;

   Int_t SampleID;
   Int_t  NtotEvents;
   float xsection;
   float EquivLumi;

   virtual vector<int>    firsttwo(Float_t * vec, vector<bool> *asso);
   bool cutID(int i, photonidcuts const& pid, std::vector<bool> *vpass = 0);
   bool cutIDEG(int i, photonidegcuts const& pid, std::vector<bool> *vpass = 0);
   bool cutIDele(int i, photonidelecuts const& pid, std::vector<bool> *vpass = 0);
   bool cutIDpresel(int i, photonidcuts const& pid, std::vector<bool> *vpass = 0);
   bool cutIDcs(int i, photonidcuts const& pid, std::vector<bool> *vpass = 0); 

   Float_t massgg;
   Float_t ptgg;
   Float_t ptphot1;
   Float_t ptphot2;
   Float_t etaphot1;
   Float_t etaphot2;
   Float_t phiphot1;
   Float_t phiphot2;
   Float_t timephot1; 
   Float_t timephot2; 
   Float_t E1phot1;
   Float_t E1phot2;
   Float_t E9phot1;
   Float_t E9phot2;
   Float_t ptjet1;
   Float_t ptjet2;
   Float_t ptcorrjet1;
   Float_t ptcorrjet2;
   Float_t etajet1;
   Float_t etajet2;
   Float_t phijet1;
   Float_t phijet2;
   Float_t deltaeta;
   Float_t zeppenjet;
   Float_t invmassjet;
   Float_t invmass2g1j;
   Float_t invmass2g2j;
   Float_t nvtx;
   Float_t met;
   Int_t isemEGphot1;
   Int_t isemEGphot2;
   Int_t idloosenewEGphot1;
   Int_t idloosenewEGphot2;
   Int_t idloose006newEGphot1;
   Int_t idloose006newEGphot2;
   Int_t idtightnewEGphot1;
   Int_t idtightnewEGphot2;
   Int_t idhggtightnewEGphot1;
   Int_t idhggtightnewEGphot2;
   Int_t idlooseEGphot1;
   Int_t idlooseEGphot2;
   Int_t idtightEGphot1;
   Int_t idtightEGphot2;
   Int_t idloosephot1; 
   Int_t idloosephot2; 
   Int_t idmediumphot1; 
   Int_t idmediumphot2; 
   Int_t idloosecsphot1;
   Int_t idloosecsphot2;
   Int_t idmediumcsphot1;
   Int_t idmediumcsphot2;
   Int_t idelephot1;
   Int_t idelephot2;
   Int_t     pid_haspixelseedphot1; 
   Int_t     pid_haspixelseedphot2; 
   Int_t     pid_isEMphot1;
   Int_t     pid_isEMphot2;
   Float_t   pid_jurECALphot1;
   Float_t   pid_jurECALphot2;
   Float_t   pid_twrHCALphot1;
   Float_t   pid_twrHCALphot2;
   Float_t   pid_HoverEphot1;
   Float_t   pid_HoverEphot2;
   Float_t   pid_hlwTrackphot1;
   Float_t   pid_hlwTrackphot2;
   Float_t   pid_etawidphot1;
   Float_t   pid_etawidphot2;
   Float_t   pid_sminphot1;
   Float_t   pid_sminphot2;
   Float_t   pid_smajphot1;
   Float_t   pid_smajphot2;
   Int_t     pid_ntrkphot1;
   Int_t     pid_ntrkphot2;
   Float_t   pid_ptisophot1;
   Float_t   pid_ptisophot2;
   Int_t     pid_ntrkcsphot1; 
   Int_t     pid_ntrkcsphot2; 
   Float_t   pid_ptisocsphot1; 
   Float_t   pid_ptisocsphot2; 
   Float_t   pid_ecalisophot1;
   Float_t   pid_ecalisophot2;
   Float_t   pid_hcalisophot1;
   Float_t   pid_hcalisophot2;
   Int_t runRN;
   Int_t eventRN;
   Int_t lumi;

};
#endif
