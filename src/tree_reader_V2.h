//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Tue Dec  7 12:08:12 2010 by ROOT version 5.22/00d
// from TTree pippo/Analysis tree
// found on file: dcap://cmsrm-se01.roma1.infn.it//pnfs/roma1.infn.it/data/cms/store/user/rahatlou/data/Cert_132440-149442_7TeV_StreamExpress_Collisions10_JSON_v3_NOVRERECO_2010B_smaller_new2/output_100_1_4fg.root
//////////////////////////////////////////////////////////

#ifndef tree_reader_V2_h
#define tree_reader_V2_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include<vector>
#include<string>
using std::vector;
using std::string;

/*

//Photon ID DelRe                                                            
struct photonidcuts {
  int tracknb;
  float trackiso_rel;
  float ecaliso_rel;
  float ecaliso_abs;
  float hcaliso_rel;
  float hcaliso_abs;
  float sminmin;
  float sminmin_min;
  float smajmaj;
};

struct photonidegcuts {
  float trackiso_rel;
  float trackiso_abs;
  float ecaliso_rel;
  float ecaliso_abs;
  float hovereiso;
  float hcaliso_rel;
  float hcaliso_abs;
  float setaetaEB;
  float setaetaEE;
};

struct photonidelecuts {
  float trackiso_relEB;
  float ecaliso_relEB;
  float hovereisoEB;
  float hcaliso_relEB;
  float setaetaEB;
  float detaEB;
  float dphiEB;
  int minhitsEB;
  float dcotEB;
  float distEB;
  float trackiso_relEE;
  float ecaliso_relEE;
  float hovereisoEE;
  float hcaliso_relEE;
  float setaetaEE;
  float detaEE;
  float dphiEE;
  int minhitsEE;
  float dcotEE;
  float distEE;
};
*/

class tree_reader_V2 {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

   // Declaration of leaf types
   Float_t         genpt;
   Bool_t          isMC;
   Int_t           store;
   Int_t           lbn;
   Int_t           bx;
   Int_t           orbit;
   Int_t           run;
   Int_t           event;
   Int_t           nMC;
   Int_t           pdgIdMC[150];   //[nMC]
   Int_t           statusMC[150];   //[nMC]
   Int_t           motherIDMC[150];   //[nMC]
   Float_t         ptMC [150];   //[nMC]
   Float_t         eMC  [150];   //[nMC]
   Float_t         etaMC[150];   //[nMC]
   Float_t         phiMC[150];   //[nMC]
   Int_t           nPhot;
   Float_t         ptPhot [40];   //[nPhot]
   Float_t         ePhot  [40];   //[nPhot]
   Float_t         escPhot  [40];   //[nPhot]
   Float_t         eseedPhot  [40];   //[nPhot]
   Float_t         etaPhot[40];   //[nPhot]
   Float_t         phiPhot[40];   //[nPhot]
   Float_t         timePhot[40];   //[nPhot]
   Float_t         e4SwissCrossPhot[40];   //[nPhot]
   Int_t           hasPixelSeedPhot[40];   //[nPhot]
   Bool_t          pid_isEM[40];   //[nPhot]
   Bool_t          pid_isLoose[40];   //[nPhot]
   Bool_t          pid_isTight[40];   //[nPhot]
   Float_t         pid_jurECAL[40];   //[nPhot]
   Float_t         pid_twrHCAL[40];   //[nPhot]
   Float_t         pid_HoverE[40];   //[nPhot]
   Float_t         pid_hlwTrack[40];   //[nPhot]
   Float_t         pid_etawid[40];   //[nPhot]
   Float_t         ptiso0015Phot[40];   //[nPhot]
   Int_t           ntrkiso0015Phot[40];   //[nPhot]
   Float_t         ptiso035Phot[40];   //[nPhot]
   Int_t           ntrkiso035Phot[40];   //[nPhot]
   Float_t         ptiso04Phot[40];   //[nPhot]
   Int_t           ntrkiso04Phot[40];   //[nPhot]
   Float_t         hcalovecal04Phot[40];   //[nPhot]
   Float_t         ecaliso04Phot[40];   //[nPhot]
   Float_t         sMajMajPhot[40];   //[nPhot]
   Float_t         sMinMinPhot[40];   //[nPhot]
   Float_t         alphaPhot[40];   //[nPhot]
   Float_t         sEtaEtaPhot[40];   //[nPhot]
   Float_t         sEtaPhiPhot[40];   //[nPhot]
   Float_t         sPhiPhiPhot[40];   //[nPhot]
   Float_t         E1Phot[40];   //[nPhot]
   Float_t         E9Phot[40];   //[nPhot]
   Float_t         E25Phot[40];   //[nPhot]
   Int_t           ieleassocPhot[40];   //[nPhot]
   Int_t           nElePhot;
   Float_t         pid_jurECALElePhot [40];   //[nElePhot]
   Float_t         pid_twrHCALElePhot [40];   //[nElePhot]
   Float_t         pid_HoverEElePhot [40];   //[nElePhot]
   Float_t         pid_hlwTrackElePhot [40];   //[nElePhot]
   Float_t         pid_etawidElePhot [40];   //[nElePhot]
   Float_t         pid_dphivtxElePhot [40];   //[nElePhot]
   Float_t         pid_detavtxElePhot [40];   //[nElePhot]
   Float_t         pid_dcotElePhot [40];   //[nElePhot]
   Float_t         pid_distElePhot [40];   //[nElePhot]
   Float_t         pid_ptElePhot [40];   //[nElePhot]
   Int_t           pid_mishitsElePhot [40];   //[nElePhot]
   Int_t           nJet_akt5;
   Float_t         ptJet_akt5 [100];   //[nJet_akt5]
   Float_t         ptCorrJet_akt5 [100];   //[nJet_akt5]
   Float_t         eJet_akt5  [100];   //[nJet_akt5]
   Float_t         etaJet_akt5[100];   //[nJet_akt5]
   Float_t         phiJet_akt5[100];   //[nJet_akt5]
   Float_t         emfJet_akt5[100];   //[nJet_akt5]
   Float_t         n90Jet_akt5[100];   //[nJet_akt5]
   Float_t         n90HitsJet_akt5[100];   //[nJet_akt5]
   Float_t         fHPDJet_akt5[100];   //[nJet_akt5]
   Float_t         fRBXJet_akt5[100];   //[nJet_akt5]
   Int_t           nJet_akt7;
   Float_t         ptJet_akt7 [100];   //[nJet_akt7]
   Float_t         ptCorrJet_akt7 [100];   //[nJet_akt5]
   Float_t         eJet_akt7  [100];   //[nJet_akt7]
   Float_t         etaJet_akt7[100];   //[nJet_akt7]
   Float_t         phiJet_akt7[100];   //[nJet_akt7]
   Float_t         emfJet_akt7[100];   //[nJet_akt7]
   Float_t         n90Jet_akt7[100];   //[nJet_akt7]
   Float_t         n90HitsJet_akt7[100];   //[nJet_akt7]
   Float_t         fHPDJet_akt7[100];   //[nJet_akt7]
   Float_t         fRBXJet_akt7[100];   //[nJet_akt7]
   Int_t           nJet_jptak5;
   Float_t         ptJet_jptak5 [100];   //[nJet_jptak5]
   Float_t         eJet_jptak5  [100];   //[nJet_jptak5]
   Float_t         etaJet_jptak5[100];   //[nJet_jptak5]
   Float_t         phiJet_jptak5[100];   //[nJet_jptak5]
   Float_t         emfJet_jptak5[100];   //[nJet_jptak5]
   Int_t           nJet_pfkt4;
   Float_t         ptJet_pfkt4 [100];   //[nJet_pfkt4]
   Float_t         eJet_pfkt4  [100];   //[nJet_pfkt4]
   Float_t         etaJet_pfkt4[100];   //[nJet_pfkt4]
   Float_t         phiJet_pfkt4[100];   //[nJet_pfkt4]
   Int_t           nJet_pfakt5;
   Float_t         ptJet_pfakt5 [100];   //[nJet_pfakt5]
   Float_t         ptCorrJet_pfakt5 [100];   //[nJet_pfakt5]
   Float_t         eJet_pfakt5  [100];   //[nJet_pfakt5]
   Float_t         etaJet_pfakt5[100];   //[nJet_pfakt5]
   Float_t         phiJet_pfakt5[100];   //[nJet_pfakt5]
   Int_t           nChargedHadrons_pfakt5[100];   //[nJet_pfakt5]
   Int_t           nPhotons_pfakt5[100];   //[nJet_pfakt5]
   Int_t           nMuons_pfakt5[100];   //[nJet_pfakt5]
   Int_t           nElectrons_pfakt5[100];   //[nJet_pfakt5]
   Int_t           nNeutralHadrons_pfakt5[100];   //[nJet_pfakt5]
   Int_t           nHFHadrons_pfakt5[100];   //[nJet_pfakt5]
   Int_t           nHFEM_pfakt5[100];   //[nJet_pfakt5]
   Float_t         eChargedHadrons_pfakt5[100];   //[nJet_pfakt5]
   Float_t         ePhotons_pfakt5[100];   //[nJet_pfakt5]
   Float_t         eMuons_pfakt5[100];   //[nJet_pfakt5]
   Float_t         eElectrons_pfakt5[100];   //[nJet_pfakt5]
   Float_t         eNeutralHadrons_pfakt5[100];   //[nJet_pfakt5]
   Float_t         eHFHadrons_pfakt5[100];   //[nJet_pfakt5]
   Float_t         eHFEM_pfakt5[100];   //[nJet_pfakt5]
   Int_t           nJet_pfakt7;
   Float_t         ptJet_pfakt7 [100];   //[nJet_pfakt7]
   Float_t         ptCorrJet_pfakt7 [100];   //[nJet_pfakt7]
   Float_t         eJet_pfakt7  [100];   //[nJet_pfakt7]
   Float_t         etaJet_pfakt7[100];   //[nJet_pfakt7]
   Float_t         phiJet_pfakt7[100];   //[nJet_pfakt7]
   Int_t           nJet_pfkt6;
   Float_t         ptJet_pfkt6 [100];   //[nJet_pfkt6]
   Float_t         eJet_pfkt6  [100];   //[nJet_pfkt6]
   Float_t         etaJet_pfkt6[100];   //[nJet_pfkt6]
   Float_t         phiJet_pfkt6[100];   //[nJet_pfkt6]
   Int_t           nJetGen_akt5;
   Float_t         ptJetGen_akt5 [100];   //[nJetGen_akt5]
   Float_t         eJetGen_akt5  [100];   //[nJetGen_akt5]
   Float_t         etaJetGen_akt5[100];   //[nJetGen_akt5]
   Float_t         phiJetGen_akt5[100];   //[nJetGen_akt5]
   Int_t           nMuonsGen_akt5[100];   //[nJetGen_akt5]
   Int_t           nElectronsGen_akt5[100];   //[nJetGen_akt5]
   Int_t           nPhotonsGen_akt5[100];   //[nJetGen_akt5]
   Int_t           nTracksGen_akt5[100];   //[nJetGen_akt5]
   Int_t           nNeutralHadronsGen_akt5[100];   //[nJetGen_akt5]
   Int_t           nHFHadronsGen_akt5[100];   //[nJetGen_akt5]
   Int_t           nHFEMGen_akt5[100];   //[nJetGen_akt5]
   Int_t           nNeutronsGen_akt5[100];   //[nJetGen_akt5]
   Int_t           nK0LGen_akt5[100];   //[nJetGen_akt5]
   Int_t           nK0SGen_akt5[100];   //[nJetGen_akt5]
   Int_t           nLambdasGen_akt5[100];   //[nJetGen_akt5]
   Int_t           nCsiGen_akt5[100];   //[nJetGen_akt5]
   Int_t           nOtherNeutralHadronsGen_akt5[100];   //[nJetGen_akt5]
   Float_t         eMuonsGen_akt5[100];   //[nJetGen_akt5]
   Float_t         eElectronsGen_akt5[100];   //[nJetGen_akt5]
   Float_t         ePhotonsGen_akt5[100];   //[nJetGen_akt5]
   Float_t         eTracksGen_akt5[100];   //[nJetGen_akt5]
   Float_t         eNeutralHadronsGen_akt5[100];   //[nJetGen_akt5]
   Float_t         eHFHadronsGen_akt5[100];   //[nJetGen_akt5]
   Float_t         eHFEMGen_akt5[100];   //[nJetGen_akt5]
   Float_t         eNeutronsGen_akt5[100];   //[nJetGen_akt5]
   Float_t         eK0LGen_akt5[100];   //[nJetGen_akt5]
   Float_t         eK0SGen_akt5[100];   //[nJetGen_akt5]
   Float_t         eLambdasGen_akt5[100];   //[nJetGen_akt5]
   Float_t         eCsiGen_akt5[100];   //[nJetGen_akt5]
   Float_t         eOtherNeutralHadronsGen_akt5[100];   //[nJetGen_akt5]
   Int_t           nJetGen_akt7;
   Float_t         ptJetGen_akt7 [100];   //[nJetGen_akt7]
   Float_t         eJetGen_akt7  [100];   //[nJetGen_akt7]
   Float_t         etaJetGen_akt7[100];   //[nJetGen_akt7]
   Float_t         phiJetGen_akt7[100];   //[nJetGen_akt7]
   Int_t           nJetGen_kt4;
   Float_t         ptJetGen_kt4 [100];   //[nJetGen_kt4]
   Float_t         eJetGen_kt4  [100];   //[nJetGen_kt4]
   Float_t         etaJetGen_kt4[100];   //[nJetGen_kt4]
   Float_t         phiJetGen_kt4[100];   //[nJetGen_kt4]
   Int_t           nJetGen_kt6;
   Float_t         ptJetGen_kt6 [100];   //[nJetGen_kt6]
   Float_t         eJetGen_kt6  [100];   //[nJetGen_kt6]
   Float_t         etaJetGen_kt6[100];   //[nJetGen_kt6]
   Float_t         phiJetGen_kt6[100];   //[nJetGen_kt6]
   Float_t         sMet  ;
   Float_t         eMet  ;
   Float_t         phiMet;
   Float_t         signifMet;
   Float_t         sCorrMet  ;
   Float_t         eCorrMet  ;
   Float_t         phiCorrMet;
   Float_t         signifCorrMet;
   Float_t         smuCorrMet  ;
   Float_t         emuCorrMet  ;
   Float_t         phimuCorrMet;
   Float_t         signifmuCorrMet;
   Float_t         sNoHFMet  ;
   Float_t         eNoHFMet  ;
   Float_t         phiNoHFMet;
   Float_t         signifNoHFMet;
   Float_t         stcMet  ;
   Float_t         etcMet  ;
   Float_t         phitcMet;
   Float_t         signiftcMet;
   Float_t         spfMet  ;
   Float_t         epfMet  ;
   Float_t         phipfMet;
   Float_t         signifpfMet;
   Float_t         sMetGen  ;
   Float_t         eMetGen  ;
   Float_t         phiMetGen;
   Float_t         signifMetGen;
   Float_t         sMetGen2  ;
   Float_t         eMetGen2  ;
   Float_t         phiMetGen2;
   Int_t           nvertex;
   Float_t         vxMC;
   Float_t         vyMC;
   Float_t         vzMC;
   Float_t         vx[100];   //[nvertex]
   Float_t         vy[100];   //[nvertex]
   Float_t         vz[100];   //[nvertex]
   Float_t         vntracks[100];   //[nvertex]
   Float_t         vchi2[100];   //[nvertex]
   Float_t         vndof[100];   //[nvertex]
   Int_t           nHLT;
   Int_t           hltNamesLen;
   vector<string>  *HLTNames;
   vector<bool>    *HLTResults;
   Double_t        Xsec;

   // List of branches
   TBranch        *b_genpt;   //!
   TBranch        *b_isMC;   //!
   TBranch        *b_store;   //!
   TBranch        *b_lbn;   //!
   TBranch        *b_bx;   //!
   TBranch        *b_orbit;   //!
   TBranch        *b_run;   //!
   TBranch        *b_event;   //!
   TBranch        *b_nMC;   //!
   TBranch        *b_pdgIdMC;   //!
   TBranch        *b_statusMC;   //!
   TBranch        *b_motherIDMC;   //!
   TBranch        *b_ptMC ;   //!
   TBranch        *b_eMC  ;   //!
   TBranch        *b_etaMC;   //!
   TBranch        *b_phiMC;   //!
   TBranch        *b_nPhot;   //!
   TBranch        *b_ptPhot ;   //!
   TBranch        *b_ePhot  ;   //!
   TBranch        *b_escPhot  ;   //!
   TBranch        *b_eseedPhot  ;   //!
   TBranch        *b_etaPhot;   //!
   TBranch        *b_phiPhot;   //!
   TBranch        *b_timePhot;   //!
   TBranch        *b_e4SwissCrossPhot;   //!
   TBranch        *b_hasPixelSeedPhot;   //!
   TBranch        *b_pid_isEM;   //!
   TBranch        *b_pid_isLoose;   //!
   TBranch        *b_pid_isTight;   //!
   TBranch        *b_pid_jurECAL;   //!
   TBranch        *b_pid_twrHCAL;   //!
   TBranch        *b_pid_HoverE;   //!
   TBranch        *b_pid_hlwTrack;   //!
   TBranch        *b_pid_etawid;   //!
   TBranch        *b_ptiso0015Phot;   //!
   TBranch        *b_ntrkiso0015Phot;   //!
   TBranch        *b_ptiso035Phot;   //!
   TBranch        *b_ntrkiso035Phot;   //!
   TBranch        *b_ptiso04Phot;   //!
   TBranch        *b_ntrkiso04Phot;   //!
   TBranch        *b_hcalovecal04Phot;   //!
   TBranch        *b_ecaliso04Phot;   //!
   TBranch        *b_sMajMajPhot;   //!
   TBranch        *b_sMinMinPhot;   //!
   TBranch        *b_alphaPhot;   //!
   TBranch        *b_sEtaEtaPhot;   //!
   TBranch        *b_sEtaPhiPhot;   //!
   TBranch        *b_sPhiPhiPhot;   //!
   TBranch        *b_E1Phot;   //!
   TBranch        *b_E9Phot;   //!
   TBranch        *b_E25Phot;   //!
   TBranch        *b_ieleassocPhot;   //!
   TBranch        *b_nElePhot;   //!
   TBranch        *b_pid_jurECALElePhot ;   //!
   TBranch        *b_pid_twrHCALElePhot ;   //!
   TBranch        *b_pid_HoverEElePhot ;   //!
   TBranch        *b_pid_hlwTrackElePhot ;   //!
   TBranch        *b_pid_etawidElePhot ;   //!
   TBranch        *b_pid_dphivtxElePhot ;   //!
   TBranch        *b_pid_dcotElePhot ;   //!
   TBranch        *b_pid_distElePhot ;   //!
   TBranch        *b_pid_detavtxElePhot ;   //!
   TBranch        *b_pid_ptElePhot ;   //!
   TBranch        *b_pid_mishitsElePhot ;   //!
   TBranch        *b_nJet_akt5;   //!
   TBranch        *b_ptJet_akt5 ;   //!
   TBranch        *b_ptCorrJet_akt5 ;   //!
   TBranch        *b_eJet_akt5  ;   //!
   TBranch        *b_etaJet_akt5;   //!
   TBranch        *b_phiJet_akt5;   //!
   TBranch        *b_emfJet_akt5;   //!
   TBranch        *b_n90Jet_akt5;   //!
   TBranch        *b_n90HitsJet_akt5;   //!
   TBranch        *b_fHPDJet_akt5;   //!
   TBranch        *b_fRBXJet_akt5;   //!
   TBranch        *b_nJet_akt7;   //!
   TBranch        *b_ptJet_akt7 ;   //!
   TBranch        *b_ptCorrJet_akt7 ;   //!
   TBranch        *b_eJet_akt7  ;   //!
   TBranch        *b_etaJet_akt7;   //!
   TBranch        *b_phiJet_akt7;   //!
   TBranch        *b_emfJet_akt7;   //!
   TBranch        *b_n90Jet_akt7;   //!
   TBranch        *b_n90HitsJet_akt7;   //!
   TBranch        *b_fHPDJet_akt7;   //!
   TBranch        *b_fRBXJet_akt7;   //!
   TBranch        *b_nJet_jptak5;   //!
   TBranch        *b_ptJet_jptak5 ;   //!
   TBranch        *b_eJet_jptak5  ;   //!
   TBranch        *b_etaJet_jptak5;   //!
   TBranch        *b_phiJet_jptak5;   //!
   TBranch        *b_emfJet_jptak5;   //!
   TBranch        *b_nJet_pfkt4;   //!
   TBranch        *b_ptJet_pfkt4 ;   //!
   TBranch        *b_eJet_pfkt4  ;   //!
   TBranch        *b_etaJet_pfkt4;   //!
   TBranch        *b_phiJet_pfkt4;   //!
   TBranch        *b_nJet_pfakt5;   //!
   TBranch        *b_ptJet_pfakt5 ;   //!
   TBranch        *b_ptCorrJet_pfakt5 ;   //!
   TBranch        *b_eJet_pfakt5  ;   //!
   TBranch        *b_etaJet_pfakt5;   //!
   TBranch        *b_phiJet_pfakt5;   //!
   TBranch        *b_nChargedHadrons_pfakt5;   //!
   TBranch        *b_nPhotons_pfakt5;   //!
   TBranch        *b_nMuons_pfakt5;   //!
   TBranch        *b_nElectrons_pfakt5;   //!
   TBranch        *b_nNeutralHadrons_pfakt5;   //!
   TBranch        *b_nHFHadrons_pfakt5;   //!
   TBranch        *b_nHFEM_pfakt5;   //!
   TBranch        *b_eChargedHadrons_pfakt5;   //!
   TBranch        *b_ePhotons_pfakt5;   //!
   TBranch        *b_eMuons_pfakt5;   //!
   TBranch        *b_eElectrons_pfakt5;   //!
   TBranch        *b_eNeutralHadrons_pfakt5;   //!
   TBranch        *b_eHFHadrons_pfakt5;   //!
   TBranch        *b_eHFEM_pfakt5;   //!
   TBranch        *b_nJet_pfakt7;   //!
   TBranch        *b_ptJet_pfakt7 ;   //!
   TBranch        *b_ptCorrJet_pfakt7 ;   //!
   TBranch        *b_eJet_pfakt7  ;   //!
   TBranch        *b_etaJet_pfakt7;   //!
   TBranch        *b_phiJet_pfakt7;   //!
   TBranch        *b_nJet_pfkt6;   //!
   TBranch        *b_ptJet_pfkt6 ;   //!
   TBranch        *b_eJet_pfkt6  ;   //!
   TBranch        *b_etaJet_pfkt6;   //!
   TBranch        *b_phiJet_pfkt6;   //!
   TBranch        *b_nJetGen_akt5;   //!
   TBranch        *b_ptJetGen_akt5 ;   //!
   TBranch        *b_eJetGen_akt5  ;   //!
   TBranch        *b_etaJetGen_akt5;   //!
   TBranch        *b_phiJetGen_akt5;   //!
   TBranch        *b_nMuonsGen_akt5;   //!
   TBranch        *b_nElectronsGen_akt5;   //!
   TBranch        *b_nPhotonsGen_akt5;   //!
   TBranch        *b_nTracksGen_akt5;   //!
   TBranch        *b_nNeutralHadronsGen_akt5;   //!
   TBranch        *b_nHFHadronsGen_akt5;   //!
   TBranch        *b_nHFEMGen_akt5;   //!
   TBranch        *b_nNeutronsGen_akt5;   //!
   TBranch        *b_nK0LGen_akt5;   //!
   TBranch        *b_nK0SGen_akt5;   //!
   TBranch        *b_nLambdasGen_akt5;   //!
   TBranch        *b_nCsiGen_akt5;   //!
   TBranch        *b_nOtherNeutralHadronsGen_akt5;   //!
   TBranch        *b_eMuonsGen_akt5;   //!
   TBranch        *b_eElectronsGen_akt5;   //!
   TBranch        *b_ePhotonsGen_akt5;   //!
   TBranch        *b_eTracksGen_akt5;   //!
   TBranch        *b_eNeutralHadronsGen_akt5;   //!
   TBranch        *b_eHFHadronsGen_akt5;   //!
   TBranch        *b_eHFEMGen_akt5;   //!
   TBranch        *b_eNeutronsGen_akt5;   //!
   TBranch        *b_eK0LGen_akt5;   //!
   TBranch        *b_eK0SGen_akt5;   //!
   TBranch        *b_eLambdasGen_akt5;   //!
   TBranch        *b_eCsiGen_akt5;   //!
   TBranch        *b_eOtherNeutralHadronsGen_akt5;   //!
   TBranch        *b_nJetGen_akt7;   //!
   TBranch        *b_ptJetGen_akt7 ;   //!
   TBranch        *b_eJetGen_akt7  ;   //!
   TBranch        *b_etaJetGen_akt7;   //!
   TBranch        *b_phiJetGen_akt7;   //!
   TBranch        *b_nJetGen_kt4;   //!
   TBranch        *b_ptJetGen_kt4 ;   //!
   TBranch        *b_eJetGen_kt4  ;   //!
   TBranch        *b_etaJetGen_kt4;   //!
   TBranch        *b_phiJetGen_kt4;   //!
   TBranch        *b_nJetGen_kt6;   //!
   TBranch        *b_ptJetGen_kt6 ;   //!
   TBranch        *b_eJetGen_kt6  ;   //!
   TBranch        *b_etaJetGen_kt6;   //!
   TBranch        *b_phiJetGen_kt6;   //!
   TBranch        *b_sMet;   //!
   TBranch        *b_eMet;   //!
   TBranch        *b_phiMet;   //!
   TBranch        *b_signifMet;   //!
   TBranch        *b_sCorrMet;   //!
   TBranch        *b_eCorrMet;   //!
   TBranch        *b_phiCorrMet;   //!
   TBranch        *b_signifCorrMet;   //!
   TBranch        *b_smuCorrMet;   //!
   TBranch        *b_emuCorrMet;   //!
   TBranch        *b_phimuCorrMet;   //!
   TBranch        *b_signifmuCorrMet;   //!
   TBranch        *b_sNoHFMet;   //!
   TBranch        *b_eNoHFMet;   //!
   TBranch        *b_phiNoHFMet;   //!
   TBranch        *b_signifNoHFMet;   //!
   TBranch        *b_stcMet;   //!
   TBranch        *b_etcMet;   //!
   TBranch        *b_phitcMet;   //!
   TBranch        *b_signiftcMet;   //!
   TBranch        *b_spfMet;   //!
   TBranch        *b_epfMet;   //!
   TBranch        *b_phipfMet;   //!
   TBranch        *b_signifpfMet;   //!
   TBranch        *b_sMetGen;   //!
   TBranch        *b_eMetGen;   //!
   TBranch        *b_phiMetGen;   //!
   TBranch        *b_signifMetGen;   //!
   TBranch        *b_sMetGen2;   //!
   TBranch        *b_eMetGen2;   //!
   TBranch        *b_phiMetGen2;   //!
   TBranch        *b_nvertex;   //!
   TBranch        *b_vxMC;   //!
   TBranch        *b_vyMC;   //!
   TBranch        *b_vzMC;   //!
   TBranch        *b_vx;   //!
   TBranch        *b_vy;   //!
   TBranch        *b_vz;   //!
   TBranch        *b_vntracks;   //!
   TBranch        *b_vchi2;   //!
   TBranch        *b_vndof;   //!
   TBranch        *b_nHLT;   //!
   TBranch        *b_hltNamesLen;   //!
   TBranch        *b_HLTNames;   //!
   TBranch        *b_HLTResults;   //!
   TBranch        *b_Xsec;   //!

   tree_reader_V2(TTree *tree=0);
   virtual ~tree_reader_V2();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif
