//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Wed Nov  3 11:17:06 2010 by ROOT version 5.22/00d
// from TTree pippo/Analysis tree
// found on file: /cmshome/fanellic/CMSSW_3_6_3/src/rootfiles/crab_VBF/output_1_1_NnK.root
//////////////////////////////////////////////////////////

#ifndef tree_reader_V1_h
#define tree_reader_V1_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

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

class tree_reader_V1 {

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
   Float_t         ptMC[150];   //[nMC]
   Float_t         eMC[150];   //[nMC]
   Float_t         etaMC[150];   //[nMC]
   Float_t         phiMC[150];   //[nMC]
   Int_t           nSIM;
   Int_t           pdgIdSIM[150];   //[nSIM]
   Int_t           statusSIM[150];   //[nSIM]
   Float_t         ptSIM[150];   //[nSIM]
   Float_t         eSIM[150];   //[nSIM]
   Float_t         etaSIM[150];   //[nSIM]
   Float_t         phiSIM[150];   //[nSIM]
   Float_t         rSIM[150];   //[nSIM]
   Float_t         zSIM[150];   //[nSIM]
   Int_t           nPF;
   Int_t           pdgIdPF[150];   //[nPF]
   Float_t         ptPF[150];   //[nPF]
   Float_t         ePF[150];   //[nPF]
   Float_t         etaPF[150];   //[nPF]
   Float_t         phiPF[150];   //[nPF]
   Int_t           nPhot;
   Float_t         ptPhot[40];   //[nPhot]
   Float_t         ePhot[40];   //[nPhot]
   Float_t         escPhot[40];   //[nPhot]
   Float_t         eseedPhot[40];   //[nPhot]
   Float_t         etaPhot[40];   //[nPhot]
   Float_t         phiPhot[40];   //[nPhot]
   Float_t         timePhot[40];   //[nPhot]
   Float_t         e4SwissCrossPhot[40];   //[nPhot]
   Int_t           nconvPhot;
   Float_t         chi2convPhot[10];   //[nconvPhot]
   Float_t         ndofconvPhot[10];   //[nconvPhot]
   Float_t         rconvPhot[10];   //[nconvPhot]
   Float_t         phiconvPhot[10];   //[nconvPhot]
   Float_t         zconvPhot[10];   //[nconvPhot]
   Int_t           ntrkconvPhot[10];   //[nconvPhot]
   Float_t         eovpconvPhot[10];   //[nconvPhot]
   Float_t         etaecalconvPhot[10];   //[nconvPhot]
   Float_t         phiecalconvPhot[10];   //[nconvPhot]
   Float_t         eecalconvPhot[10];   //[nconvPhot]
   Int_t           algoconvPhot[10];   //[nconvPhot]
   Float_t         d0convPhot[10];   //[nconvPhot]
   Float_t         detaecalconvPhot[10];   //[nconvPhot]
   Float_t         dphiecalconvPhot[10];   //[nconvPhot]
   Float_t         dphivtxconvPhot[10];   //[nconvPhot]
   Float_t         pairsepconvPhot[10];   //[nconvPhot]
   Float_t         pairmassconvPhot[10];   //[nconvPhot]
   Float_t         trchi21convPhot[10];   //[nconvPhot]
   Float_t         trndof1convPhot[10];   //[nconvPhot]
   Int_t           trqual1convPhot[10];   //[nconvPhot]
   Float_t         trpt1convPhot[10];   //[nconvPhot]
   Float_t         trerr1convPhot[10];   //[nconvPhot]
   Float_t         trchi22convPhot[10];   //[nconvPhot]
   Float_t         trndof2convPhot[10];   //[nconvPhot]
   Int_t           trqual2convPhot[10];   //[nconvPhot]
   Float_t         trpt2convPhot[10];   //[nconvPhot]
   Float_t         trerr2convPhot[10];   //[nconvPhot]
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
   Float_t         ptiso05Phot[40];   //[nPhot]
   Int_t           ntrkiso05Phot[40];   //[nPhot]
   Float_t         ptiso07Phot[40];   //[nPhot]
   Int_t           ntrkiso07Phot[40];   //[nPhot]
   Float_t         ptiso1Phot[40];   //[nPhot]
   Int_t           ntrkiso1Phot[40];   //[nPhot]
   Float_t         hcalovecal01Phot[40];   //[nPhot]
   Float_t         hcalovecal015Phot[40];   //[nPhot]
   Float_t         hcalovecal04Phot[40];   //[nPhot]
   Float_t         hcalovecal05Phot[40];   //[nPhot]
   Float_t         hcalovecal07Phot[40];   //[nPhot]
   Float_t         hcalovecal1Phot[40];   //[nPhot]
   Float_t         ecaliso01Phot[40];   //[nPhot]
   Float_t         ecaliso015Phot[40];   //[nPhot]
   Float_t         ecaliso04Phot[40];   //[nPhot]
   Float_t         ecaliso05Phot[40];   //[nPhot]
   Float_t         ecaliso07Phot[40];   //[nPhot]
   Float_t         ecaliso1Phot[40];   //[nPhot]
   Float_t         LATPhot[40];   //[nPhot]
   Float_t         sMajMajPhot[40];   //[nPhot]
   Float_t         sMinMinPhot[40];   //[nPhot]
   Float_t         alphaPhot[40];   //[nPhot]
   Float_t         sEtaEtaPhot[40];   //[nPhot]
   Float_t         sEtaPhiPhot[40];   //[nPhot]
   Float_t         sPhiPhiPhot[40];   //[nPhot]
   Float_t         E1Phot[40];   //[nPhot]
   Float_t         E9Phot[40];   //[nPhot]
   Float_t         E25Phot[40];   //[nPhot]
   Float_t         FisherPhot[40];   //[nPhot]
   Int_t           nJet_ite;
   Float_t         ptJet_ite[100];   //[nJet_ite]
   Float_t         eJet_ite [100];   //[nJet_ite]
   Float_t         etaJet_ite[100];   //[nJet_ite]
   Float_t         phiJet_ite[100];   //[nJet_ite]
   Float_t         emfJet_ite[100];   //[nJet_ite]
   Int_t           nJet_kt4;
   Float_t         ptJet_kt4[100];   //[nJet_kt4]
   Float_t         eJet_kt4[100];   //[nJet_kt4]
   Float_t         etaJet_kt4[100];   //[nJet_kt4]
   Float_t         phiJet_kt4[100];   //[nJet_kt4]
   Float_t         emfJet_kt4[100];   //[nJet_kt4]
   Int_t           nJet_kt6;
   Float_t         ptJet_kt6[100];   //[nJet_kt6]
   Float_t         eJet_kt6[100];   //[nJet_kt6]
   Float_t         etaJet_kt6[100];   //[nJet_kt6]
   Float_t         phiJet_kt6[100];   //[nJet_kt6]
   Float_t         emfJet_kt6[100];   //[nJet_kt6]
   Int_t           nJet_akt5;
   Float_t         ptJet_akt5[100];   //[nJet_akt5]
   Float_t         ptCorrJet_akt5[100];   //[nJet_akt5]
   Float_t         eJet_akt5[100];   //[nJet_akt5]
   Float_t         etaJet_akt5[100];   //[nJet_akt5]
   Float_t         phiJet_akt5[100];   //[nJet_akt5]
   Float_t         emfJet_akt5[100];   //[nJet_akt5]
   Int_t           nJet_sis5;
   Float_t         ptJet_sis5[100];   //[nJet_sis5]
   Float_t         eJet_sis5[100];   //[nJet_sis5]
   Float_t         etaJet_sis5[100];   //[nJet_sis5]
   Float_t         phiJet_sis5[100];   //[nJet_sis5]
   Float_t         emfJet_sis5[100];   //[nJet_sis5]
   Int_t           nJet_sis7;
   Float_t         ptJet_sis7[100];   //[nJet_sis7]
   Float_t         eJet_sis7[100];   //[nJet_sis7]
   Float_t         etaJet_sis7[100];   //[nJet_sis7]
   Float_t         phiJet_sis7[100];   //[nJet_sis7]
   Float_t         emfJet_sis7[100];   //[nJet_sis7]
   Int_t           nJet_jptak5;
   Float_t         ptJet_jptak5[100];   //[nJet_jptak5]
   Float_t         eJet_jptak5[100];   //[nJet_jptak5]
   Float_t         etaJet_jptak5[100];   //[nJet_jptak5]
   Float_t         phiJet_jptak5[100];   //[nJet_jptak5]
   Float_t         emfJet_jptak5[100];   //[nJet_jptak5]
   Int_t           nJet_pfite;
   Float_t         ptJet_pfite[100];   //[nJet_pfite]
   Float_t         eJet_pfite[100];   //[nJet_pfite]
   Float_t         etaJet_pfite[100];   //[nJet_pfite]
   Float_t         phiJet_pfite[100];   //[nJet_pfite]
   Int_t           nJet_pfkt4;
   Float_t         ptJet_pfkt4[100];   //[nJet_pfkt4]
   Float_t         eJet_pfkt4[100];   //[nJet_pfkt4]
   Float_t         etaJet_pfkt4[100];   //[nJet_pfkt4]
   Float_t         phiJet_pfkt4[100];   //[nJet_pfkt4]
   Int_t           nJet_pfakt5;
   Float_t         ptJet_pfakt5[100];   //[nJet_pfakt5]
   Float_t         ptCorrJet_pfakt5[100];   //[nJet_pfakt5]
   Float_t         eJet_pfakt5[100];   //[nJet_pfakt5]
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
   Float_t         ptChargedHadrons_pfakt5[100];   //[nJet_pfakt5]
   Float_t         ptPhotons_pfakt5[100];   //[nJet_pfakt5]
   Float_t         ptMuons_pfakt5[100];   //[nJet_pfakt5]
   Float_t         ptElectrons_pfakt5[100];   //[nJet_pfakt5]
   Float_t         ptNeutralHadrons_pfakt5[100];   //[nJet_pfakt5]
   Float_t         ptHFHadrons_pfakt5[100];   //[nJet_pfakt5]
   Float_t         ptHFEM_pfakt5[100];   //[nJet_pfakt5]
   Float_t         phiChargedHadrons_pfakt5[100];   //[nJet_pfakt5]
   Float_t         phiPhotons_pfakt5[100];   //[nJet_pfakt5]
   Float_t         phiMuons_pfakt5[100];   //[nJet_pfakt5]
   Float_t         phiElectrons_pfakt5[100];   //[nJet_pfakt5]
   Float_t         phiNeutralHadrons_pfakt5[100];   //[nJet_pfakt5]
   Float_t         phiHFHadrons_pfakt5[100];   //[nJet_pfakt5]
   Float_t         phiHFEM_pfakt5[100];   //[nJet_pfakt5]
   Float_t         etaChargedHadrons_pfakt5[100];   //[nJet_pfakt5]
   Float_t         etaPhotons_pfakt5[100];   //[nJet_pfakt5]
   Float_t         etaMuons_pfakt5[100];   //[nJet_pfakt5]
   Float_t         etaElectrons_pfakt5[100];   //[nJet_pfakt5]
   Float_t         etaNeutralHadrons_pfakt5[100];   //[nJet_pfakt5]
   Float_t         etaHFHadrons_pfakt5[100];   //[nJet_pfakt5]
   Float_t         etaHFEM_pfakt5[100];   //[nJet_pfakt5]
   Int_t           nJet_pfakt7;
   Float_t         ptJet_pfakt7[100];   //[nJet_pfakt7]
   Float_t         ptCorrJet_pfakt7[100];   //[nJet_pfakt7]
   Float_t         eJet_pfakt7[100];   //[nJet_pfakt7]
   Float_t         etaJet_pfakt7[100];   //[nJet_pfakt7]
   Float_t         phiJet_pfakt7[100];   //[nJet_pfakt7]
   Int_t           nChargedHadrons_pfakt7[100];   //[nJet_pfakt7]
   Int_t           nPhotons_pfakt7[100];   //[nJet_pfakt7]
   Int_t           nMuons_pfakt7[100];   //[nJet_pfakt7]
   Int_t           nElectrons_pfakt7[100];   //[nJet_pfakt7]
   Int_t           nNeutralHadrons_pfakt7[100];   //[nJet_pfakt7]
   Int_t           nHFHadrons_pfakt7[100];   //[nJet_pfakt7]
   Int_t           nHFEM_pfakt7[100];   //[nJet_pfakt7]
   Float_t         eChargedHadrons_pfakt7[100];   //[nJet_pfakt7]
   Float_t         ePhotons_pfakt7[100];   //[nJet_pfakt7]
   Float_t         eMuons_pfakt7[100];   //[nJet_pfakt7]
   Float_t         eElectrons_pfakt7[100];   //[nJet_pfakt7]
   Float_t         eNeutralHadrons_pfakt7[100];   //[nJet_pfakt7]
   Float_t         eHFHadrons_pfakt7[100];   //[nJet_pfakt7]
   Float_t         eHFEM_pfakt7[100];   //[nJet_pfakt7]
   Float_t         ptChargedHadrons_pfakt7[100];   //[nJet_pfakt7]
   Float_t         ptPhotons_pfakt7[100];   //[nJet_pfakt7]
   Float_t         ptMuons_pfakt7[100];   //[nJet_pfakt7]
   Float_t         ptElectrons_pfakt7[100];   //[nJet_pfakt7]
   Float_t         ptNeutralHadrons_pfakt7[100];   //[nJet_pfakt7]
   Float_t         ptHFHadrons_pfakt7[100];   //[nJet_pfakt7]
   Float_t         ptHFEM_pfakt7[100];   //[nJet_pfakt7]
   Float_t         phiChargedHadrons_pfakt7[100];   //[nJet_pfakt7]
   Float_t         phiPhotons_pfakt7[100];   //[nJet_pfakt7]
   Float_t         phiMuons_pfakt7[100];   //[nJet_pfakt7]
   Float_t         phiElectrons_pfakt7[100];   //[nJet_pfakt7]
   Float_t         phiNeutralHadrons_pfakt7[100];   //[nJet_pfakt7]
   Float_t         phiHFHadrons_pfakt7[100];   //[nJet_pfakt7]
   Float_t         phiHFEM_pfakt7[100];   //[nJet_pfakt7]
   Float_t         etaChargedHadrons_pfakt7[100];   //[nJet_pfakt7]
   Float_t         etaPhotons_pfakt7[100];   //[nJet_pfakt7]
   Float_t         etaMuons_pfakt7[100];   //[nJet_pfakt7]
   Float_t         etaElectrons_pfakt7[100];   //[nJet_pfakt7]
   Float_t         etaNeutralHadrons_pfakt7[100];   //[nJet_pfakt7]
   Float_t         etaHFHadrons_pfakt7[100];   //[nJet_pfakt7]
   Float_t         etaHFEM_pfakt7[100];   //[nJet_pfakt7]
   Int_t           nJet_pfsis5;
   Float_t         ptJet_pfsis5[100];   //[nJet_pfsis5]
   Float_t         eJet_pfsis5[100];   //[nJet_pfsis5]
   Float_t         etaJet_pfsis5[100];   //[nJet_pfsis5]
   Float_t         phiJet_pfsis5[100];   //[nJet_pfsis5]
   Int_t           nJet_pfkt6;
   Float_t         ptJet_pfkt6[100];   //[nJet_pfkt6]
   Float_t         eJet_pfkt6[100];   //[nJet_pfkt6]
   Float_t         etaJet_pfkt6[100];   //[nJet_pfkt6]
   Float_t         phiJet_pfkt6[100];   //[nJet_pfkt6]
   Int_t           nJet_pfsis7;
   Float_t         ptJet_pfsis7[100];   //[nJet_pfsis7]
   Float_t         eJet_pfsis7[100];   //[nJet_pfsis7]
   Float_t         etaJet_pfsis7[100];   //[nJet_pfsis7]
   Float_t         phiJet_pfsis7[100];   //[nJet_pfsis7]
   Int_t           nJetGen_ite;
   Float_t         ptJetGen_ite[100];   //[nJetGen_ite]
   Float_t         eJetGen_ite[100];   //[nJetGen_ite]
   Float_t         etaJetGen_ite[100];   //[nJetGen_ite]
   Float_t         phiJetGen_ite[100];   //[nJetGen_ite]
   Int_t           nJetGen_akt5;
   Float_t         ptJetGen_akt5[100];   //[nJetGen_akt5]
   Float_t         eJetGen_akt5[100];   //[nJetGen_akt5]
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
   Float_t         ptMuonsGen_akt5[100];   //[nJetGen_akt5]
   Float_t         ptElectronsGen_akt5[100];   //[nJetGen_akt5]
   Float_t         ptPhotonsGen_akt5[100];   //[nJetGen_akt5]
   Float_t         ptTracksGen_akt5[100];   //[nJetGen_akt5]
   Float_t         ptNeutralHadronsGen_akt5[100];   //[nJetGen_akt5]
   Float_t         ptHFHadronsGen_akt5[100];   //[nJetGen_akt5]
   Float_t         ptHFEMGen_akt5[100];   //[nJetGen_akt5]
   Float_t         phiMuonsGen_akt5[100];   //[nJetGen_akt5]
   Float_t         phiElectronsGen_akt5[100];   //[nJetGen_akt5]
   Float_t         phiPhotonsGen_akt5[100];   //[nJetGen_akt5]
   Float_t         phiTracksGen_akt5[100];   //[nJetGen_akt5]
   Float_t         phiNeutralHadronsGen_akt5[100];   //[nJetGen_akt5]
   Float_t         phiHFHadronsGen_akt5[100];   //[nJetGen_akt5]
   Float_t         phiHFEMGen_akt5[100];   //[nJetGen_akt5]
   Float_t         etaMuonsGen_akt5[100];   //[nJetGen_akt5]
   Float_t         etaElectronsGen_akt5[100];   //[nJetGen_akt5]
   Float_t         etaPhotonsGen_akt5[100];   //[nJetGen_akt5]
   Float_t         etaTracksGen_akt5[100];   //[nJetGen_akt5]
   Float_t         etaNeutralHadronsGen_akt5[100];   //[nJetGen_akt5]
   Float_t         etaHFHadronsGen_akt5[100];   //[nJetGen_akt5]
   Float_t         etaHFEMGen_akt5[100];   //[nJetGen_akt5]
   Int_t           nJetGen_akt7;
   Float_t         ptJetGen_akt7[100];   //[nJetGen_akt7]
   Float_t         eJetGen_akt7[100];   //[nJetGen_akt7]
   Float_t         etaJetGen_akt7[100];   //[nJetGen_akt7]
   Float_t         phiJetGen_akt7[100];   //[nJetGen_akt7]
   Int_t           nJetGen_kt4;
   Float_t         ptJetGen_kt4[100];   //[nJetGen_kt4]
   Float_t         eJetGen_kt4[100];   //[nJetGen_kt4]
   Float_t         etaJetGen_kt4[100];   //[nJetGen_kt4]
   Float_t         phiJetGen_kt4[100];   //[nJetGen_kt4]
   Int_t           nJetGen_kt6;
   Float_t         ptJetGen_kt6[100];   //[nJetGen_kt6]
   Float_t         eJetGen_kt6[100];   //[nJetGen_kt6]
   Float_t         etaJetGen_kt6[100];   //[nJetGen_kt6]
   Float_t         phiJetGen_kt6[100];   //[nJetGen_kt6]
   Int_t           nJetGen_sis5;
   Float_t         ptJetGen_sis5[100];   //[nJetGen_sis5]
   Float_t         eJetGen_sis5[100];   //[nJetGen_sis5]
   Float_t         etaJetGen_sis5[100];   //[nJetGen_sis5]
   Float_t         phiJetGen_sis5[100];   //[nJetGen_sis5]
   Int_t           nJetGen_sis7;
   Float_t         ptJetGen_sis7[100];   //[nJetGen_sis7]
   Float_t         eJetGen_sis7[100];   //[nJetGen_sis7]
   Float_t         etaJetGen_sis7[100];   //[nJetGen_sis7]
   Float_t         phiJetGen_sis7[100];   //[nJetGen_sis7]
   Float_t         sMet  ;
   Float_t         eMet  ;
   Float_t         phiMet;
   Float_t         stcMet  ;
   Float_t         etcMet  ;
   Float_t         phitcMet;
   Float_t         spfMet  ;
   Float_t         epfMet  ;
   Float_t         phipfMet;
   Float_t         sMetGen  ;
   Float_t         eMetGen  ;
   Float_t         phiMetGen;
   Float_t         sMetGen2  ;
   Float_t         eMetGen2  ;
   Float_t         phiMetGen2;
   Int_t           nvertex;
   Float_t         vxMC;
   Float_t         vyMC;
   Float_t         vzMC;
   Float_t         vx;
   Float_t         vy;
   Float_t         vz;
   Float_t         vntracks;
   Float_t         vchi2;
   Float_t         vndof;
   Bool_t          hltPass;
   Int_t           nHLT;
   Int_t           hltNamesLen;
   Char_t          HLTNames[6000];   //[hltNamesLen]
   Bool_t          HLTResults[200];   //[nHLT]

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
   TBranch        *b_nSIM;   //!
   TBranch        *b_pdgIdSIM;   //!
   TBranch        *b_statusSIM;   //!
   TBranch        *b_ptSIM ;   //!
   TBranch        *b_eSIM  ;   //!
   TBranch        *b_etaSIM;   //!
   TBranch        *b_phiSIM;   //!
   TBranch        *b_rSIM;   //!
   TBranch        *b_zSIM;   //!
   TBranch        *b_nPF;   //!
   TBranch        *b_pdgIdPF;   //!
   TBranch        *b_ptPF ;   //!
   TBranch        *b_ePF  ;   //!
   TBranch        *b_etaPF;   //!
   TBranch        *b_phiPF;   //!
   TBranch        *b_nPhot;   //!
   TBranch        *b_ptPhot ;   //!
   TBranch        *b_ePhot  ;   //!
   TBranch        *b_escPhot  ;   //!
   TBranch        *b_eseedPhot  ;   //!
   TBranch        *b_etaPhot;   //!
   TBranch        *b_phiPhot;   //!
   TBranch        *b_timePhot;   //!
   TBranch        *b_e4SwissCrossPhot;   //!
   TBranch        *b_nconvPhot;   //!
   TBranch        *b_chi2convPhot;   //!
   TBranch        *b_ndofconvPhot;   //!
   TBranch        *b_rconvPhot;   //!
   TBranch        *b_phiconvPhot;   //!
   TBranch        *b_zconvPhot;   //!
   TBranch        *b_ntrkconvPhot;   //!
   TBranch        *b_eovpconvPhot;   //!
   TBranch        *b_etaecalconvPhot;   //!
   TBranch        *b_phiecalconvPhot;   //!
   TBranch        *b_eecalconvPhot;   //!
   TBranch        *b_algoconvPhot;   //!
   TBranch        *b_d0convPhot;   //!
   TBranch        *b_detaecalconvPhot;   //!
   TBranch        *b_dphiecalconvPhot;   //!
   TBranch        *b_dphivtxconvPhot;   //!
   TBranch        *b_pairsepconvPhot;   //!
   TBranch        *b_pairmassconvPhot;   //!
   TBranch        *b_trchi21convPhot;   //!
   TBranch        *b_trndof1convPhot;   //!
   TBranch        *b_trqual1convPhot;   //!
   TBranch        *b_trpt1convPhot;   //!
   TBranch        *b_trerr1convPhot;   //!
   TBranch        *b_trchi22convPhot;   //!
   TBranch        *b_trndof2convPhot;   //!
   TBranch        *b_trqual2convPhot;   //!
   TBranch        *b_trpt2convPhot;   //!
   TBranch        *b_trerr2convPhot;   //!
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
   TBranch        *b_ptiso05Phot;   //!
   TBranch        *b_ntrkiso05Phot;   //!
   TBranch        *b_ptiso07Phot;   //!
   TBranch        *b_ntrkiso07Phot;   //!
   TBranch        *b_ptiso1Phot;   //!
   TBranch        *b_ntrkiso1Phot;   //!
   TBranch        *b_hcalovecal01Phot;   //!
   TBranch        *b_hcalovecal015Phot;   //!
   TBranch        *b_hcalovecal04Phot;   //!
   TBranch        *b_hcalovecal05Phot;   //!
   TBranch        *b_hcalovecal07Phot;   //!
   TBranch        *b_hcalovecal1Phot;   //!
   TBranch        *b_ecaliso01Phot;   //!
   TBranch        *b_ecaliso015Phot;   //!
   TBranch        *b_ecaliso04Phot;   //!
   TBranch        *b_ecaliso05Phot;   //!
   TBranch        *b_ecaliso07Phot;   //!
   TBranch        *b_ecaliso1Phot;   //!
   TBranch        *b_LATPhot;   //!
   TBranch        *b_sMajMajPhot;   //!
   TBranch        *b_sMinMinPhot;   //!
   TBranch        *b_alphaPhot;   //!
   TBranch        *b_sEtaEtaPhot;   //!
   TBranch        *b_sEtaPhiPhot;   //!
   TBranch        *b_sPhiPhiPhot;   //!
   TBranch        *b_E1Phot;   //!
   TBranch        *b_E9Phot;   //!
   TBranch        *b_E25Phot;   //!
   TBranch        *b_FisherPhot;   //!
   TBranch        *b_nJet_ite;   //!
   TBranch        *b_ptJet_ite ;   //!
   TBranch        *b_eJet_ite  ;   //!
   TBranch        *b_etaJet_ite;   //!
   TBranch        *b_phiJet_ite;   //!
   TBranch        *b_emfJet_ite;   //!
   TBranch        *b_nJet_kt4;   //!
   TBranch        *b_ptJet_kt4 ;   //!
   TBranch        *b_eJet_kt4  ;   //!
   TBranch        *b_etaJet_kt4;   //!
   TBranch        *b_phiJet_kt4;   //!
   TBranch        *b_emfJet_kt4;   //!
   TBranch        *b_nJet_kt6;   //!
   TBranch        *b_ptJet_kt6 ;   //!
   TBranch        *b_eJet_kt6  ;   //!
   TBranch        *b_etaJet_kt6;   //!
   TBranch        *b_phiJet_kt6;   //!
   TBranch        *b_emfJet_kt6;   //!
   TBranch        *b_nJet_akt5;   //!
   TBranch        *b_ptJet_akt5 ;   //!
   TBranch        *b_ptCorrJet_akt5 ;   //!
   TBranch        *b_eJet_akt5  ;   //!
   TBranch        *b_etaJet_akt5;   //!
   TBranch        *b_phiJet_akt5;   //!
   TBranch        *b_emfJet_akt5;   //!
   TBranch        *b_nJet_sis5;   //!
   TBranch        *b_ptJet_sis5 ;   //!
   TBranch        *b_eJet_sis5  ;   //!
   TBranch        *b_etaJet_sis5;   //!
   TBranch        *b_phiJet_sis5;   //!
   TBranch        *b_emfJet_sis5;   //!
   TBranch        *b_nJet_sis7;   //!
   TBranch        *b_ptJet_sis7 ;   //!
   TBranch        *b_eJet_sis7  ;   //!
   TBranch        *b_etaJet_sis7;   //!
   TBranch        *b_phiJet_sis7;   //!
   TBranch        *b_emfJet_sis7;   //!
   TBranch        *b_nJet_jptak5;   //!
   TBranch        *b_ptJet_jptak5 ;   //!
   TBranch        *b_eJet_jptak5  ;   //!
   TBranch        *b_etaJet_jptak5;   //!
   TBranch        *b_phiJet_jptak5;   //!
   TBranch        *b_emfJet_jptak5;   //!
   TBranch        *b_nJet_pfite;   //!
   TBranch        *b_ptJet_pfite ;   //!
   TBranch        *b_eJet_pfite  ;   //!
   TBranch        *b_etaJet_pfite;   //!
   TBranch        *b_phiJet_pfite;   //!
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
   TBranch        *b_ptChargedHadrons_pfakt5;   //!
   TBranch        *b_ptPhotons_pfakt5;   //!
   TBranch        *b_ptMuons_pfakt5;   //!
   TBranch        *b_ptElectrons_pfakt5;   //!
   TBranch        *b_ptNeutralHadrons_pfakt5;   //!
   TBranch        *b_ptHFHadrons_pfakt5;   //!
   TBranch        *b_ptHFEM_pfakt5;   //!
   TBranch        *b_phiChargedHadrons_pfakt5;   //!
   TBranch        *b_phiPhotons_pfakt5;   //!
   TBranch        *b_phiMuons_pfakt5;   //!
   TBranch        *b_phiElectrons_pfakt5;   //!
   TBranch        *b_phiNeutralHadrons_pfakt5;   //!
   TBranch        *b_phiHFHadrons_pfakt5;   //!
   TBranch        *b_phiHFEM_pfakt5;   //!
   TBranch        *b_etaChargedHadrons_pfakt5;   //!
   TBranch        *b_etaPhotons_pfakt5;   //!
   TBranch        *b_etaMuons_pfakt5;   //!
   TBranch        *b_etaElectrons_pfakt5;   //!
   TBranch        *b_etaNeutralHadrons_pfakt5;   //!
   TBranch        *b_etaHFHadrons_pfakt5;   //!
   TBranch        *b_etaHFEM_pfakt5;   //!
   TBranch        *b_nJet_pfakt7;   //!
   TBranch        *b_ptJet_pfakt7 ;   //!
   TBranch        *b_ptCorrJet_pfakt7 ;   //!
   TBranch        *b_eJet_pfakt7  ;   //!
   TBranch        *b_etaJet_pfakt7;   //!
   TBranch        *b_phiJet_pfakt7;   //!
   TBranch        *b_nChargedHadrons_pfakt7;   //!
   TBranch        *b_nPhotons_pfakt7;   //!
   TBranch        *b_nMuons_pfakt7;   //!
   TBranch        *b_nElectrons_pfakt7;   //!
   TBranch        *b_nNeutralHadrons_pfakt7;   //!
   TBranch        *b_nHFHadrons_pfakt7;   //!
   TBranch        *b_nHFEM_pfakt7;   //!
   TBranch        *b_eChargedHadrons_pfakt7;   //!
   TBranch        *b_ePhotons_pfakt7;   //!
   TBranch        *b_eMuons_pfakt7;   //!
   TBranch        *b_eElectrons_pfakt7;   //!
   TBranch        *b_eNeutralHadrons_pfakt7;   //!
   TBranch        *b_eHFHadrons_pfakt7;   //!
   TBranch        *b_eHFEM_pfakt7;   //!
   TBranch        *b_ptChargedHadrons_pfakt7;   //!
   TBranch        *b_ptPhotons_pfakt7;   //!
   TBranch        *b_ptMuons_pfakt7;   //!
   TBranch        *b_ptElectrons_pfakt7;   //!
   TBranch        *b_ptNeutralHadrons_pfakt7;   //!
   TBranch        *b_ptHFHadrons_pfakt7;   //!
   TBranch        *b_ptHFEM_pfakt7;   //!
   TBranch        *b_phiChargedHadrons_pfakt7;   //!
   TBranch        *b_phiPhotons_pfakt7;   //!
   TBranch        *b_phiMuons_pfakt7;   //!
   TBranch        *b_phiElectrons_pfakt7;   //!
   TBranch        *b_phiNeutralHadrons_pfakt7;   //!
   TBranch        *b_phiHFHadrons_pfakt7;   //!
   TBranch        *b_phiHFEM_pfakt7;   //!
   TBranch        *b_etaChargedHadrons_pfakt7;   //!
   TBranch        *b_etaPhotons_pfakt7;   //!
   TBranch        *b_etaMuons_pfakt7;   //!
   TBranch        *b_etaElectrons_pfakt7;   //!
   TBranch        *b_etaNeutralHadrons_pfakt7;   //!
   TBranch        *b_etaHFHadrons_pfakt7;   //!
   TBranch        *b_etaHFEM_pfakt7;   //!
   TBranch        *b_nJet_pfsis5;   //!
   TBranch        *b_ptJet_pfsis5 ;   //!
   TBranch        *b_eJet_pfsis5  ;   //!
   TBranch        *b_etaJet_pfsis5;   //!
   TBranch        *b_phiJet_pfsis5;   //!
   TBranch        *b_nJet_pfkt6;   //!
   TBranch        *b_ptJet_pfkt6 ;   //!
   TBranch        *b_eJet_pfkt6  ;   //!
   TBranch        *b_etaJet_pfkt6;   //!
   TBranch        *b_phiJet_pfkt6;   //!
   TBranch        *b_nJet_pfsis7;   //!
   TBranch        *b_ptJet_pfsis7 ;   //!
   TBranch        *b_eJet_pfsis7  ;   //!
   TBranch        *b_etaJet_pfsis7;   //!
   TBranch        *b_phiJet_pfsis7;   //!
   TBranch        *b_nJetGen_ite;   //!
   TBranch        *b_ptJetGen_ite ;   //!
   TBranch        *b_eJetGen_ite  ;   //!
   TBranch        *b_etaJetGen_ite;   //!
   TBranch        *b_phiJetGen_ite;   //!
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
   TBranch        *b_ptMuonsGen_akt5;   //!
   TBranch        *b_ptElectronsGen_akt5;   //!
   TBranch        *b_ptPhotonsGen_akt5;   //!
   TBranch        *b_ptTracksGen_akt5;   //!
   TBranch        *b_ptNeutralHadronsGen_akt5;   //!
   TBranch        *b_ptHFHadronsGen_akt5;   //!
   TBranch        *b_ptHFEMGen_akt5;   //!
   TBranch        *b_phiMuonsGen_akt5;   //!
   TBranch        *b_phiElectronsGen_akt5;   //!
   TBranch        *b_phiPhotonsGen_akt5;   //!
   TBranch        *b_phiTracksGen_akt5;   //!
   TBranch        *b_phiNeutralHadronsGen_akt5;   //!
   TBranch        *b_phiHFHadronsGen_akt5;   //!
   TBranch        *b_phiHFEMGen_akt5;   //!
   TBranch        *b_etaMuonsGen_akt5;   //!
   TBranch        *b_etaElectronsGen_akt5;   //!
   TBranch        *b_etaPhotonsGen_akt5;   //!
   TBranch        *b_etaTracksGen_akt5;   //!
   TBranch        *b_etaNeutralHadronsGen_akt5;   //!
   TBranch        *b_etaHFHadronsGen_akt5;   //!
   TBranch        *b_etaHFEMGen_akt5;   //!
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
   TBranch        *b_nJetGen_sis5;   //!
   TBranch        *b_ptJetGen_sis5;   //!
   TBranch        *b_eJetGen_sis5  ;   //!
   TBranch        *b_etaJetGen_sis5;   //!
   TBranch        *b_phiJetGen_sis5;   //!
   TBranch        *b_nJetGen_sis7;   //!
   TBranch        *b_ptJetGen_sis7 ;   //!
   TBranch        *b_eJetGen_sis7  ;   //!
   TBranch        *b_etaJetGen_sis7;   //!
   TBranch        *b_phiJetGen_sis7;   //!
   TBranch        *b_sMet;   //!
   TBranch        *b_eMet;   //!
   TBranch        *b_phiMet;   //!
   TBranch        *b_stcMet;   //!
   TBranch        *b_etcMet;   //!
   TBranch        *b_phitcMet;   //!
   TBranch        *b_spfMet;   //!
   TBranch        *b_epfMet;   //!
   TBranch        *b_phipfMet;   //!
   TBranch        *b_sMetGen;   //!
   TBranch        *b_eMetGen;   //!
   TBranch        *b_phiMetGen;   //!
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
   TBranch        *b_hltPass;   //!
   TBranch        *b_nHLT;   //!
   TBranch        *b_hltNamesLen;   //!
   TBranch        *b_HLTNames;   //!
   TBranch        *b_HLTResults;   //!

   tree_reader_V1(TTree *tree=0);
   virtual ~tree_reader_V1();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
   //   virtual vector<int>    firsttwo(Float_t * vec, vector<bool> *asso);

};

#endif
