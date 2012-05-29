//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Fri Apr 29 21:17:22 2011 by ROOT version 5.27/06b
// from TTree pippo/Analysis tree
// found on file: dcap:///pnfs/roma1.infn.it/data/cms/store/user/rahatlou/MC/41xV7/QCD_Pt-30to40_doubleEMEnriched_TuneZ2_7TeV-pythia6-41x_ntpv1/output_89_1_z0g.root
//////////////////////////////////////////////////////////

#ifndef tree_reader_V7_h
#define tree_reader_V7_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

#include <string>
#include <vector>
using std::vector;
using std::string;

//#define SMALL_VERTEX_VECTOR

#ifdef SMALL_VERTEX_VECTOR
#define VECTOR_ARRAY_SIZE 60
#else
#define VECTOR_ARRAY_SIZE 100
#endif

class tree_reader_V7 {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

   // Declaration of leaf types
   Float_t         genpt;
   Int_t           genProcessId;
   Float_t         genQScale;
   Double_t        qPDF;
   Double_t        x1PDF;
   Double_t        x2PDF;
   Double_t        id1PDF;
   Double_t        id2PDF;
   Int_t           nWeightsPDF[10];
   Double_t        pdfWeight[10][150];
   Bool_t          isMC;
   Int_t           store;
   Int_t           lbn;
   Int_t           bx;
   Int_t           orbit;
   Int_t           run;
   Int_t           event;
   Float_t         rhoPF;
   Float_t         rhoCalo;
   Float_t         rhoAllJets;
   Int_t           nMC;
   Int_t           pdgIdMC[150];   //[nMC]
   Int_t           statusMC[150];   //[nMC]
   Int_t           motherIDMC[150];   //[nMC]
   Float_t         ptMC [150];   //[nMC]
   Float_t         eMC  [150];   //[nMC]
   Float_t         etaMC[150];   //[nMC]
   Float_t         phiMC[150];   //[nMC]
   Int_t           pu_n;
   Float_t         pu_zpos[100];   //[pu_n]
   Float_t         pu_sumpt_lowpt[100];   //[pu_n]
   Float_t         pu_sumpt_highpt[100];   //[pu_n]
   Float_t         pu_ntrks_lowpt[100];   //[pu_n]
   Float_t         pu_ntrks_highpt[100];   //[pu_n]
   Int_t           nPhot;
   Float_t         ptPhot [40];   //[nPhot]
   Float_t         ePhot  [40];   //[nPhot]
   Float_t         escPhot  [40];   //[nPhot]
   Float_t         escRawPhot  [40];   //[nPhot]
   Float_t         escRegrPhot  [40];   //[nPhot]
   Float_t         escRegrPhotError  [40];   //[nPhot]
   Float_t         escPhFixPhot  [40];   //[nPhot]
   Float_t         escPhFixPhotError  [40];   //[nPhot]
   Float_t         etascPhot  [40];   //[nPhot]
   Float_t         phiscPhot  [40];   //[nPhot]
   Float_t         xscPhot  [40];   //[nPhot]
   Float_t         yscPhot  [40];   //[nPhot]
   Float_t         zscPhot  [40];   //[nPhot]
   Float_t         xcaloPhot  [40];   //[nPhot]
   Float_t         ycaloPhot  [40];   //[nPhot]
   Float_t         zcaloPhot  [40];   //[nPhot]
   Float_t         eseedPhot  [40];   //[nPhot]
   Float_t         etaPhot[40];   //[nPhot]
   Float_t         phiPhot[40];   //[nPhot]
   Float_t         timePhot[40];   //[nPhot]
   Float_t         e4SwissCrossPhot[40];   //[nPhot]
   Int_t           hasPixelSeedPhot[40];   //[nPhot]
   Int_t           hasMatchedPromptElePhot[40];   //[nPhot]
   Int_t           hasMatchedConvPhot[40];   //[nPhot]
   Bool_t          isEBPhot[40];
   Bool_t          isEEPhot[40];
   Bool_t          isEBEEGapPhot[40];
   Int_t           ntracksConvPhot[40];   //[nPhot]
   Bool_t          isValidVtxConvPhot[40];   //[nPhot]
   Float_t         pairInvmassConvPhot[40];   //[nPhot]
   Float_t         pairCotThetaSeperationConvPhot[40];   //[nPhot]
   Float_t         pairmomentumXConvPhot[40];   //[nPhot]
   Float_t         pairmomentumYConvPhot[40];   //[nPhot]
   Float_t         pairmomentumZConvPhot[40];   //[nPhot]
   Float_t         chi2ConvPhot[40];   //[nPhot]
   Float_t         nDofConvPhot[40];   //[nPhot]
   Float_t         eOverPConvPhot[40];   //[nPhot]
   Float_t         convVxConvPhot[40];   //[nPhot]
   Float_t         convVyConvPhot[40];   //[nPhot]
   Float_t         convVzConvPhot[40];   //[nPhot]
   Float_t         distOfMinimumApproachConvPhot[40];   //[nPhot]
   Float_t         dPhiTracksAtVtxConvPhot[40];   //[nPhot]
   Bool_t          pid_isEM[40];   //[nPhot]
   Bool_t          pid_isLoose[40];   //[nPhot]
   Bool_t          pid_isTight[40];   //[nPhot]
   Float_t         pid_jurECAL[40];   //[nPhot]
   Float_t         pid_twrHCAL[40];   //[nPhot]
   Float_t         pid_HoverE[40];   //[nPhot]
   Float_t         pid_hlwTrack[40];   //[nPhot]
   Float_t         pid_hlwTrackForCiC[40][VECTOR_ARRAY_SIZE];
   Float_t         pid_hlwTrackNoDz[40];   //[nPhot]
   Float_t         pid_etawid[40];   //[nPhot]
   Float_t         pid_jurECAL03[40];   //[nPhot]
   Float_t         pid_twrHCAL03[40];   //[nPhot]
   Float_t         pid_hlwTrack03[40];   //[nPhot]
   Float_t         pid_hlwTrack03NoDz[40];   //[nPhot]
   Float_t         pid_deltaRToTrackPhot[40];   //[nPhot]
   Float_t         pid_hlwTrack03ForCiC[40][VECTOR_ARRAY_SIZE];
   Float_t         ptiso004Phot[40];   //[nPhot]
   Int_t           ntrkiso004Phot[40];   //[nPhot]
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
   Float_t         E2OverE9Phot[40];   //[nPhot]
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
   Int_t           pid_mishitsElePhot [40];   //[nElePhot]
   Float_t         pid_distElePhot [40];   //[nElePhot]
   Float_t         pid_dcotElePhot [40];   //[nElePhot]
   Float_t         pid_ptElePhot [40];   //[nElePhot]
   Int_t           nJet_akt5;
   Float_t         ptJet_akt5 [200];   //[nJet_akt5]
   Float_t         ptCorrJet_akt5 [200];   //[nJet_akt5]
   Float_t         eJet_akt5  [200];   //[nJet_akt5]
   Float_t         etaJet_akt5[200];   //[nJet_akt5]
   Float_t         phiJet_akt5[200];   //[nJet_akt5]
   Float_t         emfJet_akt5[200];   //[nJet_akt5]
   Float_t         n90Jet_akt5[200];   //[nJet_akt5]
   Float_t         n90HitsJet_akt5[200];   //[nJet_akt5]
   Float_t         fHPDJet_akt5[200];   //[nJet_akt5]
   Float_t         fRBXJet_akt5[200];   //[nJet_akt5]
   Int_t           nJet_akt7;
   Float_t         ptJet_akt7 [200];   //[nJet_akt7]
   Float_t         ptCorrJet_akt7 [200];   //[nJet_akt5]
   Float_t         eJet_akt7  [200];   //[nJet_akt7]
   Float_t         etaJet_akt7[200];   //[nJet_akt7]
   Float_t         phiJet_akt7[200];   //[nJet_akt7]
   Float_t         emfJet_akt7[200];   //[nJet_akt7]
   Float_t         n90Jet_akt7[200];   //[nJet_akt7]
   Float_t         n90HitsJet_akt7[200];   //[nJet_akt7]
   Float_t         fHPDJet_akt7[200];   //[nJet_akt7]
   Float_t         fRBXJet_akt7[200];   //[nJet_akt7]
   Int_t           nJet_pfakt5;
   Float_t         ptJet_pfakt5 [200];   //[nJet_pfakt5]
   Float_t         ptCorrJet_pfakt5 [200];   //[nJet_pfakt5]
   Float_t         eJet_pfakt5  [200];   //[nJet_pfakt5]
   Float_t         etaJet_pfakt5[200];   //[nJet_pfakt5]
   Float_t         phiJet_pfakt5[200];   //[nJet_pfakt5]
   Float_t         ptDJet_pfakt5[200];   //[nJet_pfakt5]
   Float_t         rmsCandJet_pfakt5[200];   //[nJet_pfakt5]
   Float_t         jetId_dRMean_pfakt5[100];   //[nJet_pfakt5]
   Float_t         jetId_frac01_pfakt5[100];   //[nJet_pfakt5]
   Float_t         jetId_frac02_pfakt5[100];   //[nJet_pfakt5]
   Float_t         jetId_frac03_pfakt5[100];   //[nJet_pfakt5]
   Float_t         jetId_frac04_pfakt5[100];   //[nJet_pfakt5]
   Float_t         jetId_frac05_pfakt5[100];   //[nJet_pfakt5]
   Float_t         jetId_nNeutrals_pfakt5[100];   //[nJet_pfakt5]
   Float_t         jetId_beta_pfakt5[100];   //[nJet_pfakt5]
   Float_t         jetId_betaStar_pfakt5[100];   //[nJet_pfakt5]
   Float_t         jetId_dZ_pfakt5[100];   //[nJet_pfakt5]
   Float_t         jetId_nCharged_pfakt5[100];   //[nJet_pfakt5]
   Float_t         jetId_dR2Mean_pfakt5[100];   //[nJet_pfakt5]
   Float_t         jetId_betaStarClassic_pfakt5[100];   //[nJet_pfakt5]
   Float_t         jetIdSimple_mva_pfakt5[100];   //[nJet_pfakt5]
   Int_t           jetIdSimple_wp_pfakt5[100];   //[nJet_pfakt5]
   Float_t         jetIdFull_mva_pfakt5[100];   //[nJet_pfakt5]
   Int_t           jetIdFull_wp_pfakt5[100];   //[nJet_pfakt5]
   Float_t         jetIdCutBased_mva_pfakt5[100];   //[nJet_pfakt5]
   Int_t           jetIdCutBased_wp_pfakt5[100];   //[nJet_pfakt5]
   Float_t         beta_pfakt5[100][100];
   Float_t         betaStar_pfakt5[100][100];
   Float_t         combinedSecondaryVertexBJetTags[200];   //[nJet_pfakt5]
   Float_t         combinedSecondaryVertexMVABJetTags[200];   //[nJet_pfakt5]
   Float_t         jetBProbabilityBJetTags[200];   //[nJet_pfakt5]
   Float_t         jetProbabilityBJetTags[200];   //[nJet_pfakt5]
   Float_t         simpleSecondaryVertexHighEffBJetTags[200];   //[nJet_pfakt5]
   Float_t         simpleSecondaryVertexHighPurBJetTags[200];   //[nJet_pfakt5]
   Float_t         softMuonBJetTags[200];   //[nJet_pfakt5]
   Float_t         softMuonByIP3dBJetTags[200];   //[nJet_pfakt5]
   Float_t         softMuonByPtBJetTags[200];   //[nJet_pfakt5]
   Float_t         softElectronByIP3dBJetTags[200];   //[nJet_pfakt5]
   Float_t         softElectronByPtBJetTags[200];   //[nJet_pfakt5]
   Float_t         trackCountingHighPurBJetTags[200];   //[nJet_pfakt5]
   Float_t         trackCountingHighEffBJetTags[200];   //[nJet_pfakt5]
   Int_t           nChargedHadrons_pfakt5[200];   //[nJet_pfakt5]
   Int_t           nPhotons_pfakt5[200];   //[nJet_pfakt5]
   Int_t           nMuons_pfakt5[200];   //[nJet_pfakt5]
   Int_t           nElectrons_pfakt5[200];   //[nJet_pfakt5]
   Int_t           nNeutralHadrons_pfakt5[200];   //[nJet_pfakt5]
   Int_t           nHFHadrons_pfakt5[200];   //[nJet_pfakt5]
   Int_t           nHFEM_pfakt5[200];   //[nJet_pfakt5]
   Float_t         eChargedHadrons_pfakt5[200];   //[nJet_pfakt5]
   Float_t         ePhotons_pfakt5[200];   //[nJet_pfakt5]
   Float_t         eMuons_pfakt5[200];   //[nJet_pfakt5]
   Float_t         eElectrons_pfakt5[200];   //[nJet_pfakt5]
   Float_t         eNeutralHadrons_pfakt5[200];   //[nJet_pfakt5]
   Float_t         eHFHadrons_pfakt5[200];   //[nJet_pfakt5]
   Float_t         eHFEM_pfakt5[200];   //[nJet_pfakt5]
   Int_t           nJet_pfakt7;
   Float_t         ptJet_pfakt7 [200];   //[nJet_pfakt7]
   Float_t         ptCorrJet_pfakt7 [200];   //[nJet_pfakt7]
   Float_t         eJet_pfakt7  [200];   //[nJet_pfakt7]
   Float_t         etaJet_pfakt7[200];   //[nJet_pfakt7]
   Float_t         phiJet_pfakt7[200];   //[nJet_pfakt7]
   Int_t           nJetGen_akt5;
   Float_t         ptJetGen_akt5 [200];   //[nJetGen_akt5]
   Float_t         eJetGen_akt5  [200];   //[nJetGen_akt5]
   Float_t         etaJetGen_akt5[200];   //[nJetGen_akt5]
   Float_t         phiJetGen_akt5[200];   //[nJetGen_akt5]
   Int_t           nMuonsGen_akt5[200];   //[nJetGen_akt5]
   Int_t           nElectronsGen_akt5[200];   //[nJetGen_akt5]
   Int_t           nPhotonsGen_akt5[200];   //[nJetGen_akt5]
   Int_t           nTracksGen_akt5[200];   //[nJetGen_akt5]
   Int_t           nNeutralHadronsGen_akt5[200];   //[nJetGen_akt5]
   Int_t           nHFHadronsGen_akt5[200];   //[nJetGen_akt5]
   Int_t           nHFEMGen_akt5[200];   //[nJetGen_akt5]
   Int_t           nNeutronsGen_akt5[200];   //[nJetGen_akt5]
   Int_t           nK0LGen_akt5[200];   //[nJetGen_akt5]
   Int_t           nK0SGen_akt5[200];   //[nJetGen_akt5]
   Int_t           nLambdasGen_akt5[200];   //[nJetGen_akt5]
   Int_t           nCsiGen_akt5[200];   //[nJetGen_akt5]
   Int_t           nOtherNeutralHadronsGen_akt5[200];   //[nJetGen_akt5]
   Float_t         eMuonsGen_akt5[200];   //[nJetGen_akt5]
   Float_t         eElectronsGen_akt5[200];   //[nJetGen_akt5]
   Float_t         ePhotonsGen_akt5[200];   //[nJetGen_akt5]
   Float_t         eTracksGen_akt5[200];   //[nJetGen_akt5]
   Float_t         eNeutralHadronsGen_akt5[200];   //[nJetGen_akt5]
   Float_t         eHFHadronsGen_akt5[200];   //[nJetGen_akt5]
   Float_t         eHFEMGen_akt5[200];   //[nJetGen_akt5]
   Float_t         eNeutronsGen_akt5[200];   //[nJetGen_akt5]
   Float_t         eK0LGen_akt5[200];   //[nJetGen_akt5]
   Float_t         eK0SGen_akt5[200];   //[nJetGen_akt5]
   Float_t         eLambdasGen_akt5[200];   //[nJetGen_akt5]
   Float_t         eCsiGen_akt5[200];   //[nJetGen_akt5]
   Float_t         eOtherNeutralHadronsGen_akt5[200];   //[nJetGen_akt5]
   Int_t           nJetGen_akt7;
   Float_t         ptJetGen_akt7 [200];   //[nJetGen_akt7]
   Float_t         eJetGen_akt7  [200];   //[nJetGen_akt7]
   Float_t         etaJetGen_akt7[200];   //[nJetGen_akt7]
   Float_t         phiJetGen_akt7[200];   //[nJetGen_akt7]
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
   Float_t         sglobalPfMet;
   Float_t         eglobalPfMet;
   Float_t         phiglobalPfMet;
   Float_t         signifglobalPfMet;
   Float_t         scentralPfMet;
   Float_t         ecentralPfMet;
   Float_t         phicentralPfMet;
   Float_t         signifcentralPfMet;
   Float_t         eassocPfMet[100];   //[nvertex]
   Float_t         phiassocPfMet[100];   //[nvertex]
   Float_t         signifassocPfMet[100];   //[nvertex]
   Float_t         eassocOtherVtxPfMet[100];   //[nvertex]
   Float_t         phiassocOtherVtxPfMet[100];   //[nvertex]
   Float_t         signifassocOtherVtxPfMet[100];   //[nvertex]
   Float_t         etrkPfMet[100];   //[nvertex]
   Float_t         phitrkPfMet[100];   //[nvertex]
   Float_t         signiftrkPfMet[100];   //[nvertex]
   Float_t         ecleanPfMet[100];   //[nvertex]
   Float_t         phicleanPfMet[100];   //[nvertex]
   Float_t         signifcleanPfMet[100];   //[nvertex]
   Float_t         ecleanedSaclayPfMet[100];   //[nvertex]
   Float_t         phicleanedSaclayPfMet[100];   //[nvertex]
   Float_t         signifcleanedSaclayPfMet[100];   //[nvertex]
   Float_t         eminTypeICleanSaclayPfMet[100];   //[nvertex]
   Float_t         phiminTypeICleanSaclayPfMet[100];   //[nvertex]
   Float_t         signifminTypeICleanSaclayPfMet[100];   //[nvertex]
   Float_t         globalPfSums[100];
   Float_t         spfMet  ;
   Float_t         epfMet  ;
   Float_t         phipfMet;
   Float_t         signifpfMet;
   Float_t         spfMetType1;
   Float_t         epfMetType1;
   Float_t         phipfMetType1;
   Float_t         signifpfMetType1;
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
   Float_t         vlogsumpt2[100];   //[nvertex]
   Int_t           nPreselPhotonPairs;
   Float_t         indexPreselPhot1[100];   //[nPreselPhotonPairs]
   Float_t         indexPreselPhot2[100];   //[nPreselPhotonPairs]
   Int_t           vrankPhotonPairs[100];   //[nPreselPhotonPairs]
   Float_t         vevtMvaPhotonPairs[100];   //[nPreselPhotonPairs]
   Float_t         vevtProbPhotonPairs[100];   //[nPreselPhotonPairs]
   Float_t         vptbalPhotonPairs[100];   //[nPreselPhotonPairs]
   Float_t         vptasymPhotonPairs[100];   //[nPreselPhotonPairs]
   Int_t           nHLT;
   Int_t           hltNamesLen;
   vector<string>  *HLTNames;
   vector<bool>    *HLTResults;
   Int_t           nEle;   
   Float_t         electron_px[40];   //[nEle]
   Float_t         electron_py[40];   //[nEle]
   Float_t         electron_pz[40];   //[nEle]
   Float_t         electron_vx[40];   //[nEle]
   Float_t         electron_vy[40];   //[nEle]
   Float_t         electron_vz[40];   //[nEle]
   Float_t         electron_pt[40];   //[nEle]
   Float_t         electron_eta[40];   //[nEle]
   Float_t         electron_phi[40];   //[nEle]
   Float_t         electron_energy[40];   //[nEle]
   Float_t         electron_charge[40];   //[nEle]
   Float_t         electron_fBrem[40];   //[nEle]
   Float_t         electron_dist[40];   //[nEle]
   Float_t         electron_dcot[40];   //[nEle]
   Int_t           electron_misHits[40];   //[nEle]
   Int_t           electron_seedType[40];   //[nEle]
   Float_t         electron_EoP[40];   //[nEle]
   Float_t         electron_OneOverEMinusOneOverP[40];   //[nEle]
   Float_t         electron_r9[40];   //[nEle]
   Int_t           electron_nSubClusters[40];   //[nEle]
   Float_t         electron_trkIso[40];   //[nEle]
   Float_t         electron_ecalIso[40];   //[nEle]
   Float_t         electron_hcalIso[40];   //[nEle]
   Float_t         electron_trkIso03[40];   //[nEle]
   Float_t         electron_ecalIso03[40];   //[nEle]
   Float_t         electron_hcalIso03[40];   //[nEle]
   Float_t         electron_SigmaIetaIeta[40];   //[nEle]
   Float_t         electron_SigmaIphiIphi[40];   //[nEle]
   Float_t         electron_dEtaIn[40];   //[nEle]
   Float_t         electron_dPhiIn[40];   //[nEle]
   Float_t         electron_HoE[40];   //[nEle]
   Float_t         electron_pFlowMVA[40];   //[nEle]
   Float_t         electron_sc_energy[40];   //[nEle]
   Float_t         electron_sc_eta[40];   //[nEle]
   Float_t         electron_sc_phi[40];   //[nEle]
   Bool_t          isBeamHaloGlobalLoosePass;
   Bool_t          isBeamHaloGlobalTightPass;
   Bool_t          isBeamHaloHcalLoosePass;
   Bool_t          isBeamHaloHcalTightPass;
   Bool_t          isBeamHaloCSCLoosePass;
   Bool_t          isBeamHaloCSCTightPass;
   Bool_t          isBeamHaloEcalLoosePass;
   Bool_t          isBeamHaloEcalTightPass;
   Bool_t          isBeamHaloIDTightPass;
   Bool_t          isBeamHaloIDLoosePass;
   Bool_t          isSmellsLikeHalo_Tag;
   Bool_t          isLooseHalo_Tag;
   Bool_t          isTightHalo_Tag;
   Bool_t          isExtremeTightHalo_Tag;
   Int_t           nMuons;
   Float_t         Muon_px[200];   //[nMuons]
   Float_t         Muon_py[200];   //[nMuons]
   Float_t         Muon_pz[200];   //[nMuons]
   Float_t         Muon_vx[200];   //[nMuons]
   Float_t         Muon_vy[200];   //[nMuons]
   Float_t         Muon_vz[200];   //[nMuons]
   Float_t         Muon_pt[200];   //[nMuons]
   Float_t         Muon_eta[200];   //[nMuons]
   Float_t         Muon_phi[200];   //[nMuons]
   Float_t         Muon_energy[200];   //[nMuons]
   Float_t         Muon_charge[200];   //[nMuons]
   Bool_t          Muon_isGlobalMuon[200];   //[nMuons]
   Bool_t          Muon_isTrackerMuon[200];   //[nMuons]
   Bool_t          Muon_isStandAloneMuon[200];   //[nMuons]
   Bool_t          Muon_InnerTrack_isNonnull[200];   //[nMuons]
   Bool_t          Muon_OuterTrack_isNonnull[200];   //[nMuons]
   Float_t         Muon_OuterPoint_x[200];   //[nMuons]
   Float_t         Muon_OuterPoint_y[200];   //[nMuons]
   Float_t         Muon_OuterPoint_z[200];   //[nMuons]
   Float_t         Muon_InnerPoint_x[200];   //[nMuons]
   Float_t         Muon_InnerPoint_y[200];   //[nMuons]
   Float_t         Muon_InnerPoint_z[200];   //[nMuons]
   Float_t         Muon_trackIso[200];   //[nMuons]
   Float_t         Muon_ecalIso[200];   //[nMuons]
   Float_t         Muon_hcalIso[200];   //[nMuons]
   Float_t         Muon_relIso[200];   //[nMuons]
   Int_t           Muon_normChi2[200];   //[nMuons]
   Int_t           Muon_validHits[200];   //[nMuons]
   Int_t           Muon_tkHits[200];   //[nMuons]
   Int_t           Muon_pixHits[200];   //[nMuons]
   Int_t           Muon_numberOfMatches[200];   //[nMuons]
   Int_t           nCosmicMuons;
   Float_t         CosmicMuon_px[200];   //[nCosmicMuons]
   Float_t         CosmicMuon_py[200];   //[nCosmicMuons]
   Float_t         CosmicMuon_pz[200];   //[nCosmicMuons]
   Float_t         CosmicMuon_pt[200];   //[nCosmicMuons]
   Float_t         CosmicMuon_eta[200];   //[nCosmicMuons]
   Float_t         CosmicMuon_phi[200];   //[nCosmicMuons]
   Float_t         CosmicMuon_energy[200];   //[nCosmicMuons]
   Float_t         CosmicMuon_charge[200];   //[nCosmicMuons]
   Bool_t          CosmicMuon_isGlobalMuon[200];   //[nCosmicMuons]
   Bool_t          CosmicMuon_isTrackerMuon[200];   //[nCosmicMuons]
   Bool_t          CosmicMuon_isStandAloneMuon[200];   //[nCosmicMuons]
   Bool_t          CosmicMuon_InnerTrack_isNonnull[200];   //[nCosmicMuons]
   Bool_t          CosmicMuon_OuterTrack_isNonnull[200];   //[nCosmicMuons]
   Float_t         CosmicMuon_OuterPoint_x[200];   //[nCosmicMuons]
   Float_t         CosmicMuon_OuterPoint_y[200];   //[nCosmicMuons]
   Float_t         CosmicMuon_OuterPoint_z[200];   //[nCosmicMuons]
   Double_t        Xsec;

   // List of branches
   TBranch        *b_genpt;   //!
   TBranch        *b_genProcessId;   //!
   TBranch        *b_genQScale;   //!
   TBranch        *b_qPDF;   //!
   TBranch        *b_x1PDF;   //!
   TBranch        *b_x2PDF;   //!
   TBranch        *b_id1PDF;   //!
   TBranch        *b_id2PDF;   //!
   TBranch        *b_nWeightsPDF;   //!
   TBranch        *b_pdfWeight;   //!
   TBranch        *b_isMC;   //!
   TBranch        *b_store;   //!
   TBranch        *b_lbn;   //!
   TBranch        *b_bx;   //!
   TBranch        *b_orbit;   //!
   TBranch        *b_run;   //!
   TBranch        *b_event;   //!
   TBranch        *b_rhoPF;   //!
   TBranch        *b_rhoCalo;   //!
   TBranch        *b_rhoAllJets;   //!
   TBranch        *b_nMC;   //!
   TBranch        *b_pdgIdMC;   //!
   TBranch        *b_statusMC;   //!
   TBranch        *b_motherIDMC;   //!
   TBranch        *b_ptMC ;   //!
   TBranch        *b_eMC  ;   //!
   TBranch        *b_etaMC;   //!
   TBranch        *b_phiMC;   //!
   TBranch        *b_pu_n;   //!
   TBranch        *b_pu_zpos;   //!
   TBranch        *b_pu_sumpt_lowpt;   //!
   TBranch        *b_pu_sumpt_highpt;   //!
   TBranch        *b_pu_ntrks_lowpt;   //!
   TBranch        *b_pu_ntrks_highpt;   //!
   TBranch        *b_nPhot;   //!
   TBranch        *b_ptPhot ;   //!
   TBranch        *b_ePhot  ;   //!
   TBranch        *b_escPhot  ;   //!
   TBranch        *b_escRawPhot  ;   //!
   TBranch        *b_escRegrPhot  ;   //!
   TBranch        *b_escRegrPhotError  ;   //!
   TBranch        *b_escPhFixPhot  ;   //!
   TBranch        *b_escPhFixPhotError  ;   //!
   TBranch        *b_etascPhot  ;   //!
   TBranch        *b_phiscPhot  ;   //!
   TBranch        *b_xscPhot  ;   //!
   TBranch        *b_yscPhot  ;   //!
   TBranch        *b_zscPhot  ;   //!
   TBranch        *b_xcaloPhot  ;   //!
   TBranch        *b_ycaloPhot  ;   //!
   TBranch        *b_zcaloPhot  ;   //!
   TBranch        *b_eseedPhot  ;   //!
   TBranch        *b_etaPhot;   //!
   TBranch        *b_phiPhot;   //!
   TBranch        *b_timePhot;   //!
   TBranch        *b_e4SwissCrossPhot;   //!
   TBranch        *b_hasPixelSeedPhot;   //!
   TBranch        *b_hasMatchedPromptElePhot;   //!
   TBranch        *b_hasMatchedConvPhot;   //!
   TBranch        *b_isEBPhot;   //!
   TBranch        *b_isEEPhot;   //!
   TBranch        *b_isEBEEGapPhot;   //!
   TBranch        *b_ntracksConvPhot;   //!
   TBranch        *b_isValidVtxConvPhot;   //!
   TBranch        *b_pairInvmassConvPhot;   //!
   TBranch        *b_pairCotThetaSeperationConvPhot;   //!
   TBranch        *b_pairmomentumXConvPhot;   //!
   TBranch        *b_pairmomentumYConvPhot;   //!
   TBranch        *b_pairmomentumZConvPhot;   //!
   TBranch        *b_chi2ConvPhot;   //!
   TBranch        *b_nDofConvPhot;   //!
   TBranch        *b_eOverPConvPhot;   //!
   TBranch        *b_convVxConvPhot;   //!
   TBranch        *b_convVyConvPhot;   //!
   TBranch        *b_convVzConvPhot;   //!
   TBranch        *b_distOfMinimumApproachConvPhot;   //!
   TBranch        *b_dPhiTracksAtVtxConvPhot;   //!
   TBranch        *b_pid_isEM;   //!
   TBranch        *b_pid_isLoose;   //!
   TBranch        *b_pid_isTight;   //!
   TBranch        *b_pid_jurECAL;   //!
   TBranch        *b_pid_twrHCAL;   //!
   TBranch        *b_pid_HoverE;   //!
   TBranch        *b_pid_hlwTrack;   //!
   TBranch        *b_pid_hlwTrackBestRank;   //!
   TBranch        *b_pid_hlwTrackNoDz;   //!
   TBranch        *b_pid_etawid;   //!
   TBranch        *b_pid_jurECAL03;   //!
   TBranch        *b_pid_twrHCAL03;   //!
   TBranch        *b_pid_hlwTrack03;   //!
   TBranch        *b_pid_hlwTrack03ForCiC;   //!
   TBranch        *b_pid_hlwTrack03NoDz;   //!
   TBranch        *b_ptiso004Phot;   //!
   TBranch        *b_ntrkiso004Phot;   //!
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
   TBranch        *b_E2OverE9Phot;   //!
   TBranch        *b_E9Phot;   //!
   TBranch        *b_E25Phot;   //!
   TBranch        *b_ieleassocPhot;   //!
   TBranch        *b_pid_deltaRToTrackPhot;   //!
   TBranch        *b_nElePhot;   //!
   TBranch        *b_pid_jurECALElePhot ;   //!
   TBranch        *b_pid_twrHCALElePhot ;   //!
   TBranch        *b_pid_HoverEElePhot ;   //!
   TBranch        *b_pid_hlwTrackElePhot ;   //!
   TBranch        *b_pid_etawidElePhot ;   //!
   TBranch        *b_pid_dphivtxElePhot ;   //!
   TBranch        *b_pid_detavtxElePhot ;   //!
   TBranch        *b_pid_mishitsElePhot ;   //!
   TBranch        *b_pid_distElePhot ;   //!
   TBranch        *b_pid_dcotElePhot ;   //!
   TBranch        *b_pid_ptElePhot ;   //!
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
   TBranch        *b_ptDJet_pfakt5;   //!
   TBranch        *b_rmsCandJet_pfakt5;   //!
   TBranch        *b_jetId_dRMean_pfakt5;   //!
   TBranch        *b_jetId_frac01_pfakt5;   //!
   TBranch        *b_jetId_frac02_pfakt5;   //!
   TBranch        *b_jetId_frac03_pfakt5;   //!
   TBranch        *b_jetId_frac04_pfakt5;   //!
   TBranch        *b_jetId_frac05_pfakt5;   //!
   TBranch        *b_jetId_nNeutrals_pfakt5;   //!
   TBranch        *b_jetId_beta_pfakt5;   //!
   TBranch        *b_jetId_betaStar_pfakt5;   //!
   TBranch        *b_jetId_dZ_pfakt5;   //!
   TBranch        *b_jetId_nCharged_pfakt5;   //!
   TBranch        *b_jetId_dR2Mean_pfakt5;   //!
   TBranch        *b_jetId_betaStarClassic_pfakt5;   //!
   TBranch        *b_jetIdSimple_mva_pfakt5;   //!
   TBranch        *b_jetIdSimple_wp_pfakt5;   //!
   TBranch        *b_jetIdFull_mva_pfakt5;   //!
   TBranch        *b_jetIdFull_wp_pfakt5;   //!
   TBranch        *b_jetIdCutBased_mva_pfakt5;   //!
   TBranch        *b_jetIdCutBased_wp_pfakt5;   //!
   TBranch        *b_beta_pfakt5;   //!
   TBranch        *b_betaStar_pfakt5;   //!
   TBranch        *b_combinedSecondaryVertexBJetTags;   //!
   TBranch        *b_combinedSecondaryVertexMVABJetTags;   //!
   TBranch        *b_jetBProbabilityBJetTags;   //!
   TBranch        *b_jetProbabilityBJetTags;   //!
   TBranch        *b_simpleSecondaryVertexHighEffBJetTags;   //!
   TBranch        *b_simpleSecondaryVertexHighPurBJetTags;   //!
   TBranch        *b_softMuonBJetTags;   //!
   TBranch        *b_softMuonByIP3dBJetTags;   //!
   TBranch        *b_softMuonByPtBJetTags;   //!
   TBranch        *b_softElectronByIP3dBJetTags;   //!
   TBranch        *b_softElectronByPtBJetTags;   //!
   TBranch        *b_trackCountingHighPurBJetTags;   //!
   TBranch        *b_trackCountingHighEffBJetTags;   //!
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
   TBranch        *b_sglobalPfMet;   //!
   TBranch        *b_eglobalPfMet;   //!
   TBranch        *b_phiglobalPfMet;   //!
   TBranch        *b_signifglobalPfMet;   //!
   TBranch        *b_scentralPfMet;   //!
   TBranch        *b_ecentralPfMet;   //!
   TBranch        *b_phicentralPfMet;   //!
   TBranch        *b_signifcentralPfMet;   //!
   TBranch        *b_eassocPfMet;   //!
   TBranch        *b_phiassocPfMet;   //!
   TBranch        *b_signifassocPfMet;   //!
   TBranch        *b_eassocOtherVtxPfMet;   //!
   TBranch        *b_phiassocOtherVtxPfMet;   //!
   TBranch        *b_signifassocOtherVtxPfMet;   //!
   TBranch        *b_etrkPfMet;   //!
   TBranch        *b_phitrkPfMet;   //!
   TBranch        *b_signiftrkPfMet;   //!
   TBranch        *b_ecleanPfMet;   //!
   TBranch        *b_phicleanPfMet;   //!
   TBranch        *b_signifcleanPfMet;   //!
   TBranch        *b_ecleanedSaclayPfMet;   //!
   TBranch        *b_phicleanedSaclayPfMet;   //!
   TBranch        *b_signifcleanedSaclayPfMet;   //!
   TBranch        *b_eminTypeICleanSaclayPfMet;   //!
   TBranch        *b_phiminTypeICleanSaclayPfMet;   //!
   TBranch        *b_signifminTypeICleanSaclayPfMet;   //!
   TBranch        *b_globalPfSums;   //!
   TBranch        *b_spfMet;   //!
   TBranch        *b_epfMet;   //!
   TBranch        *b_phipfMet;   //!
   TBranch        *b_signifpfMet;   //!
   TBranch        *b_spfMetType1;   //!
   TBranch        *b_epfMetType1;   //!
   TBranch        *b_phipfMetType1;   //!
   TBranch        *b_signifpfMetType1;   //!
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
   TBranch        *b_vlogsumpt2;   //!
   TBranch        *b_nPreselPhotonPairs;   //!
   TBranch        *b_indexPreselPhot1;   //!
   TBranch        *b_indexPreselPhot2;   //!
   TBranch        *b_vrankPhotonPairs;   //!
   TBranch        *b_vevtMvaPhotonPairs;   //!
   TBranch        *b_vevtProbPhotonPairs;   //!
   TBranch        *b_vptbalPhotonPairs;   //!
   TBranch        *b_vptasymPhotonPairs;   //!
   TBranch        *b_nHLT;   //!
   TBranch        *b_hltNamesLen;   //!
   TBranch        *b_HLTNames;   //!
   TBranch        *b_HLTResults;   //!
  TBranch        *b_nEle;   //!
   TBranch        *b_electron_px;   //!
   TBranch        *b_electron_py;   //!
   TBranch        *b_electron_pz;   //!
   TBranch        *b_electron_vx;   //!
   TBranch        *b_electron_vy;   //!
   TBranch        *b_electron_vz;   //!
   TBranch        *b_electron_pt;   //!
   TBranch        *b_electron_eta;   //!
   TBranch        *b_electron_phi;   //!
   TBranch        *b_electron_energy;   //!
   TBranch        *b_electron_charge;   //!
   TBranch        *b_electron_fBrem;   //!
   TBranch        *b_electron_dist;   //!
   TBranch        *b_electron_dcot;   //!
   TBranch        *b_electron_misHits;   //!
   TBranch        *b_electron_seedType;   //!
   TBranch        *b_electron_EoP;   //!
   TBranch        *b_electron_OneOverEMinusOneOverP;   //!
   TBranch        *b_electron_r9;   //!
   TBranch        *b_electron_nSubClusters;   //!
   TBranch        *b_electron_trkIso;   //!
   TBranch        *b_electron_ecalIso;   //!
   TBranch        *b_electron_hcalIso;   //!
   TBranch        *b_electron_trkIso03;   //!
   TBranch        *b_electron_ecalIso03;   //!
   TBranch        *b_electron_hcalIso03;   //!
   TBranch        *b_electron_SigmaIetaIeta;   //!
   TBranch        *b_electron_SigmaIphiIphi;   //!
   TBranch        *b_electron_dEtaIn;   //!
   TBranch        *b_electron_dPhiIn;   //!
   TBranch        *b_electron_HoE;   //!
   TBranch        *b_electron_pFlowMVA;   //!
   TBranch        *b_electron_sc_energy;   //!
   TBranch        *b_electron_sc_eta;   //!
   TBranch        *b_electron_sc_phi;   //!
   TBranch        *b_isBeamHaloGlobalLoosePass;   //!
   TBranch        *b_isBeamHaloGloablTightPass;   //!
   TBranch        *b_isBeamHaloHcalLoosePass;   //!
   TBranch        *b_isBeamHaloHcalTightPass;   //!
   TBranch        *b_isBeamHaloCSCLoosePass;   //!
   TBranch        *b_isBeamHaloCSCTightPass;   //!
   TBranch        *b_isBeamHaloEcalLoosePass;   //!
   TBranch        *b_isBeamHaloEcalTightPass;   //!
   TBranch        *b_isBeamHaloIDTightPass;   //!
   TBranch        *b_isBeamHaloIDLoosePass;   //!
   TBranch        *b_isSmellsLikeHalo_Tag;   //!
   TBranch        *b_isLooseHalo_Tag;   //!
   TBranch        *b_isTightHalo_Tag;   //!
   TBranch        *b_isExtremeTightHalo_Tag;   //!
   TBranch        *b_nMuons;   //!
   TBranch        *b_Muon_px;   //!
   TBranch        *b_Muon_py;   //!
   TBranch        *b_Muon_pz;   //!
   TBranch        *b_Muon_vx;   //!
   TBranch        *b_Muon_vy;   //!
   TBranch        *b_Muon_vz;   //!
   TBranch        *b_Muon_pt;   //!
   TBranch        *b_Muon_eta;   //!
   TBranch        *b_Muon_phi;   //!
   TBranch        *b_Muon_energy;   //!
   TBranch        *b_Muon_charge;   //!
   TBranch        *b_Muon_isGlobalMuon;   //!
   TBranch        *b_Muon_isTrackerMuon;   //!
   TBranch        *b_Muon_isStandAloneMuon;   //!
   TBranch        *b_Muon_InnerTrack_isNonnull;   //!
   TBranch        *b_Muon_OuterTrack_isNonnull;   //!
   TBranch        *b_Muon_OuterPoint_x;   //!
   TBranch        *b_Muon_OuterPoint_y;   //!
   TBranch        *b_Muon_OuterPoint_z;   //!
   TBranch        *b_Muon_InnerPoint_x;   //!
   TBranch        *b_Muon_InnerPoint_y;   //!
   TBranch        *b_Muon_InnerPoint_z;   //!
   TBranch        *b_Muon_trackIso;   //!
   TBranch        *b_Muon_ecalIso;   //!
   TBranch        *b_Muon_hcalIso;   //!
   TBranch        *b_Muon_relIso;   //!
   TBranch        *b_Muon_normChi2;   //!
   TBranch        *b_Muon_validHits;   //!
   TBranch        *b_Muon_tkHits;   //!
   TBranch        *b_Muon_pixHits;   //!
   TBranch        *b_Muon_numberOfMatches;   //!
   TBranch        *b_nCosmicMuons;   //!
   TBranch        *b_CosmicMuon_px;   //!
   TBranch        *b_CosmicMuon_py;   //!
   TBranch        *b_CosmicMuon_pz;   //!
   TBranch        *b_CosmicMuon_pt;   //!
   TBranch        *b_CosmicMuon_eta;   //!
   TBranch        *b_CosmicMuon_phi;   //!
   TBranch        *b_CosmicMuon_energy;   //!
   TBranch        *b_CosmicMuon_charge;   //!
   TBranch        *b_CosmicMuon_isGlobalMuon;   //!
   TBranch        *b_CosmicMuon_isTrackerMuon;   //!
   TBranch        *b_CosmicMuon_isStandAloneMuon;   //!
   TBranch        *b_CosmicMuon_InnerTrack_isNonnull;   //!
   TBranch        *b_CosmicMuon_OuterTrack_isNonnull;   //!
   TBranch        *b_CosmicMuon_OuterPoint_x;   //!
   TBranch        *b_CosmicMuon_OuterPoint_y;   //!
   TBranch        *b_CosmicMuon_OuterPoint_z;   //!
   TBranch        *b_Xsec;   //!

   tree_reader_V7(TTree *tree=0);
   virtual ~tree_reader_V7();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif
