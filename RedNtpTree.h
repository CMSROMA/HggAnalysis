#ifndef RedNtpTree_h
#define RedNtpTree_h

//#include "higgsanal_tree_V1.h"
//#include "tree_reader_V2.h"
//#include "tree_reader_V3.h"
//#include "tree_reader_V7.h"
#include "tree_reader_V8.h"
#include "PhotonIdCuts.h"
#include "LeptonIdCuts.h"
#include "EnergyScaleCorrection.h"
#include "JetScaleSystematics.h"
#include "ElectronEffectiveArea.h"
#include "TLorentzVector.h"

#include <TFile.h>
#include <TH1F.h>
#include <TH2D.h>
#include <TMath.h>
#include <TRandom3.h>
#include <TString.h>
#include<vector>
#include<string>
#include <TChain.h>

#include "TMVA/Reader.h"

using std::string;
using std::vector;


#define NGENJETS 200
#define NMC 150

class RedNtpTree : public tree_reader_V8 {

public:
  
  RedNtpTree(TTree *tree=0, const TString& outname="redntp.root");
    virtual ~RedNtpTree();
    virtual void     Loop(int isgjetqcd=0, char* selection = "loose");
    void SetJsonFile(const char* json) { jsonFile = json; };
    void SetPuWeights(std::string puWeightFile);
    void SetPtWeights(std::string ptWeightFile);
    void DoPDFWeighting();
    double ErrEt(double Et, double Eta);
    void SetNtotXsection(int ntot, float xsec) {      
      NtotEvents = ntot;
      xsection = xsec;
      EquivLumi = ntot/xsec;
   }
    void setEnergyScaleCorrections(TString correctionFile, TString correctionType)
   {
     std::cout << "Constructing new Scale Corrections Of Type " << correctionType<< std::endl;
     std::cout << "Constructing new Scale Corrections from file " << correctionFile << std::endl;
     scaleCorrections_=new EnergyScaleCorrection(correctionFile,correctionType);
   }
    void setJetSystematics(TString correctionFile, float typesyst)
   {
     std::cout << "Constructing JEC systematics from file " << correctionFile << std::endl;
     std::cout << "Type of JEC systematics " << typesyst << std::endl;
     jetsyst_=new JetScaleSystematics(correctionFile);
     typejetsyst_=typesyst;
   }
    
    /*
    TTree* myTree; 
    struct cicTree_structure_ {
      int runCIC;
      int eventCIC;
      float isosumoet;
      float isoecalet;
      float isohcalet;
      float isotrackeret;
      float isosumoetbad;
      float isoecaletbad;
      float isohcaletbad;
      float isotrackeretbad;
      float sieie;
      float hoe;
      float r9;
      float drtotk_25_99;
      float pixel;
    };

    cicTree_structure_ tree_;
    */

    std::string photonLevelNewIDMVA_EB;
    std::string photonLevelNewIDMVA_EE;
    std::string cicVersion;

private:
   TFile* hOutputFile ;
   TTree * ana_tree ;

   TRandom3* gen_;

   vector<bool> jetnoisophot;

   const char* jsonFile;
   
   Int_t SampleID;
   Int_t  NtotEvents;
   float xsection;
   float EquivLumi;
   bool doPDFweight;


   void SetAllRecoVarToMinus999();
   void SetAllGenVarToMinus999();

   virtual vector<int>    firstones(Float_t * vec, vector<bool> *asso, int number=4);
   bool cutID(int i, photonidcuts const& pid, std::vector<bool> *vpass = 0);
   bool cutIDEG(int i, photonidegcuts const& pid, std::vector<bool> *vpass = 0,  bool pu = 0);
   bool cutIDele(int i, photonidelecuts const& pid, std::vector<bool> *vpass = 0);
   bool cutIDpresel(int i, photonidcuts const& pid, std::vector<bool> *vpass = 0);
   bool cutIDcs(int i, photonidcuts const& pid, std::vector<bool> *vpass = 0); 
   bool mcID(int i); 
   bool assoJet(int i);
   void correctPhotons(bool energyRegression);
   void correctJets(int scale, float smear);
   TLorentzVector correctMet(TLorentzVector uncormet, bool smearing = 1, bool scale = 0, bool PUremoval = 0);
   TLorentzVector shiftMet(TLorentzVector uncormet);

   // defines photon CiC ID cuts for all cut levels
   enum phoCiCIDLevel { phoNOCUTS=0, phoLOOSE, phoMEDIUM, phoTIGHT, phoSUPERTIGHT, phoHYPERTIGHT1, phoHYPERTIGHT2, phoHYPERTIGHT3, phoHYPERTIGHT4, phoHYPERTIGHT5, phoHYPERTIGHT6, phoHYPERTIGHT7, phoHYPERTIGHT8, phoHYPERTIGHT9, phoNCUTLEVELS };
   enum phoCiCCuts { phoISOSUMOET=0,  phoISOSUMOETBAD,   phoTRKISOOETOM,   phoSIEIE,   phoHOVERE,   phoR9,   phoDRTOTK_25_99,   phoPIXEL, phoNCUTS };
   enum phoCiC6Categories { phoCiC6EBhighR9=0, phoCiC6EBmidR9, phoCiC6EBlowR9, phoCiC6EEhighR9, phoCiC6EEmidR9, phoCiC6EElowR9, phoCiC6NCATEGORIES };
   enum phoCiC4Categories { phoCiC4EBhighR9=0, phoCiC4EBlowR9, phoCiC4EEhighR9, phoCiC4EElowR9, phoCiC4NCATEGORIES };

   void SetPhotonCutsInCategories(phoCiCIDLevel cutlevel, float * cic6_cuts_lead, float * cic6_cuts_sublead, float * cic4_cuts_lead, float * cic4_cuts_sublead, float*, float*);
   void FillPhotonCiCSelectionVariable(int photon_index, int vtx_index);

   float cic6_cut_lead_isosumoet[phoNCUTLEVELS][6];
   float cic6_cut_lead_isosumoetbad[phoNCUTLEVELS][6];
   float cic6_cut_lead_trkisooet[phoNCUTLEVELS][6];
   float cic6_cut_lead_sieie[phoNCUTLEVELS][6];
   float cic6_cut_lead_hovere[phoNCUTLEVELS][6];
   float cic6_cut_lead_r9[phoNCUTLEVELS][6];
   float cic6_cut_lead_drtotk_25_99[phoNCUTLEVELS][6];
   float cic6_cut_lead_pixel[phoNCUTLEVELS][6];
   float cic6_cut_sublead_isosumoet[phoNCUTLEVELS][6];
   float cic6_cut_sublead_isosumoetbad[phoNCUTLEVELS][6];
   float cic6_cut_sublead_trkisooet[phoNCUTLEVELS][6];
   float cic6_cut_sublead_sieie[phoNCUTLEVELS][6];
   float cic6_cut_sublead_hovere[phoNCUTLEVELS][6];
   float cic6_cut_sublead_r9[phoNCUTLEVELS][6];
   float cic6_cut_sublead_drtotk_25_99[phoNCUTLEVELS][6];
   float cic6_cut_sublead_pixel[phoNCUTLEVELS][6];
   
   float cic4_cut_lead_isosumoet[phoNCUTLEVELS][4];
   float cic4_cut_lead_isosumoetbad[phoNCUTLEVELS][4];
   float cic4_cut_lead_trkisooet[phoNCUTLEVELS][4];
   float cic4_cut_lead_sieie[phoNCUTLEVELS][4];
   float cic4_cut_lead_hovere[phoNCUTLEVELS][4];
   float cic4_cut_lead_r9[phoNCUTLEVELS][4];
   float cic4_cut_lead_drtotk_25_99[phoNCUTLEVELS][4];
   float cic4_cut_lead_pixel[phoNCUTLEVELS][4];
   float cic4_cut_sublead_isosumoet[phoNCUTLEVELS][4];
   float cic4_cut_sublead_isosumoetbad[phoNCUTLEVELS][4];
   float cic4_cut_sublead_trkisooet[phoNCUTLEVELS][4];
   float cic4_cut_sublead_sieie[phoNCUTLEVELS][4];
   float cic4_cut_sublead_hovere[phoNCUTLEVELS][4];
   float cic4_cut_sublead_r9[phoNCUTLEVELS][4];
   float cic4_cut_sublead_drtotk_25_99[phoNCUTLEVELS][4];
   float cic4_cut_sublead_pixel[phoNCUTLEVELS][4];

   float cic4pf_cut_lead_isosumoet[phoNCUTLEVELS][4];
   float cic4pf_cut_lead_isosumoetbad[phoNCUTLEVELS][4];
   float cic4pf_cut_lead_trkisooet[phoNCUTLEVELS][4];
   float cic4pf_cut_lead_sieie[phoNCUTLEVELS][4];
   float cic4pf_cut_lead_hovere[phoNCUTLEVELS][4];
   float cic4pf_cut_lead_r9[phoNCUTLEVELS][4];
   float cic4pf_cut_lead_drtotk_25_99[phoNCUTLEVELS][4];
   float cic4pf_cut_lead_pixel[phoNCUTLEVELS][4];
   
   float cic4pf_cut_sublead_isosumoet[phoNCUTLEVELS][4];
   float cic4pf_cut_sublead_isosumoetbad[phoNCUTLEVELS][4];
   float cic4pf_cut_sublead_trkisooet[phoNCUTLEVELS][4];
   float cic4pf_cut_sublead_sieie[phoNCUTLEVELS][4];
   float cic4pf_cut_sublead_hovere[phoNCUTLEVELS][4];
   float cic4pf_cut_sublead_r9[phoNCUTLEVELS][4];
   float cic4pf_cut_sublead_drtotk_25_99[phoNCUTLEVELS][4];
   float cic4pf_cut_sublead_pixel[phoNCUTLEVELS][4];

   TH1F* cic4_cut_isosumoet[phoCiC4NCATEGORIES];
   TH1F* cic4_cut_isosumoetbad[phoCiC4NCATEGORIES];
   TH1F* cic4_cut_trkisooet[phoCiC4NCATEGORIES];
   TH1F* cic4_cut_sieie[phoCiC4NCATEGORIES];
   TH1F* cic4_cut_hovere[phoCiC4NCATEGORIES];
   TH1F* cic4_cut_r9[phoCiC4NCATEGORIES];
   TH1F* cic4_cut_drtotk_25_99[phoCiC4NCATEGORIES];
   TH1F* cic4_cut_pixel[phoCiC4NCATEGORIES];

   void SetAllMVA(); 
   int   PhotonCiCSelectionLevel( int photon_index, bool electronVeto, int vertex_index, bool usePF);
   bool  PhotonMITPreSelection( int photon_index, int vertex_index, bool electronVeto) ;

   //  TMVA::Reader *tmvaReader_dipho_MIT;
  TMVA::Reader *tmvaReaderID_Single_Barrel, *tmvaReaderID_Single_Endcap;
  Float_t tmva_photonid_pfchargedisogood03;
  Float_t tmva_photonid_pfchargedisobad03;
  Float_t tmva_photonid_pfphotoniso03;
  Float_t tmva_photonid_pfneutraliso03;
  Float_t tmva_photonid_sieie;
  Float_t tmva_photonid_sieip;
  Float_t tmva_photonid_etawidth;
  Float_t tmva_photonid_phiwidth;
  Float_t tmva_photonid_r9;
  Float_t tmva_photonid_s4ratio;
  Float_t tmva_photonid_lambdaratio;
  Float_t tmva_photonid_sceta;
  Float_t tmva_photonid_eventrho;
  Float_t tmva_photonid_ESEffSigmaRR;
  Float_t PhotonIDMVANew(Int_t iPhoton, Int_t vtx);  
  
  //photon category functions (r9 and eta)
  int PhotonCategory(int photonindex) { 
    return PhotonR9Category(photonindex) + 2*PhotonEtaCategory(photonindex);
  }
   Int_t PhotonR9Category(int photonindex) { 
     if(photonindex < 0) return -1;
     int r9cat = (Int_t)(E9Phot[photonindex]/escRawPhot[photonindex]<0.94);// 0, 1(high r9 --> low r9)
     return r9cat;
   }
   int PhotonEtaCategory(int photonindex) {
     if(photonindex < 0) return -1;
     //int etacat = (Int_t)(!isEBPhot[photonindex]);   // 0, 1 (barrel --> endcap)
     int etacat = (Int_t)(TMath::Abs(etascPhot[photonindex])>1.479);   // 0, 1 (barrel --> endcap)
     return  etacat;
   }


   // lepton tag
   bool leptonCutsEle2011(int iEle, electronidcuts const& pid, vector<bool> *vpass);
   bool leptonCutsEle2012(int iEle, electronidcuts2012 const& pid, vector<bool> *vpass);
   bool leptonCutsMu2011(int iMu, muonidcuts const& pid, vector<bool> *vpass);
   bool leptonCutsMu2012(int iMu, muonidcuts2012 const& pid, vector<bool> *vpass);
   double eleDzPV(int iele, int iPV);
   double eleDxyPV(int iele, int iPV);
   double muonDzPV(int imu, int iPV);
   double muonDxyPV(int imu, int iPV);
   double trackDzPV(TVector3 PVPos, TVector3 trackVPos, TVector3 trackMom);
   double trackDxyPV(TVector3 PVPos, TVector3 trackVPos, TVector3 trackMom);

   // gen level info
   int countLOGenGamma();
   int countISRGenGamma();
   int countFSRGenGamma();

   // vector of pu weights
   std::vector<Double_t> puweights_;
   TH1D* ptweights_;
   TH2D* jetDR;
   TH2D* jetresp_vs_pt;
   TH2D* jetresp_vs_eta;
   TH2D* jetresp_vs_npu;
   TH2D* jetresp_vs_eta_50;
   TH2D* jetresp_vs_npu_50;
   TH2D* jetresp_vs_eta_150;
   TH2D* jetresp_vs_eta_50_abs;
   TH2D* jetresp_vs_npu_150;
   TH2D* jetresp_vs_pt_forward;
   TH2D* jetresp_vs_npu_forward;
  
   EnergyScaleCorrection* scaleCorrections_;
   JetScaleSystematics* jetsyst_;
   Float_t typejetsyst_;
 
   Float_t massgg;
   Float_t ptgg;
   Float_t phigg;
   Float_t etagg;
   Float_t massggnewvtx;
   Float_t ptggnewvtx;
   Float_t ptphot1;
   Float_t ptphot2;
   Float_t deltaRToTrackphot1;
   Float_t deltaRToTrackphot2;
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
   Float_t energyErrphot1;
   Float_t energyErrphot2;
   Float_t energySmearingphot1;
   Float_t energySmearingphot2;
   Int_t njets;
   Float_t ptjet[10];
   Float_t ptcorrjet[10];
   Float_t ecorrjet[10];
   Float_t etajet[10];
   Float_t phijet[10];
   Float_t betajet[10];
   Float_t betastarjet[10];
   Float_t btagvtxjet[10];
   Float_t btagtrkjet[10];
   Float_t btagjprobjet[10];
   Float_t ptDjet[10];
   Float_t rmsjet[10];
   Int_t ntrkjet[10];
   Int_t nneutjet[10];
   Float_t jetIdSimple_mvajet[10];
   Float_t jetIdFull_mvajet[10];
   Float_t jetId_dR2Meanjet[10];
   Float_t jetId_betaStarClassicjet[10];
   Float_t jetId_frac01jet[10];
   Float_t jetId_frac02jet[10];
   Float_t jetId_frac03jet[10];
   Float_t jetId_frac04jet[10];
   Float_t jetId_frac05jet[10];
   Float_t jetId_betajet[10];
   Float_t jetId_betaStarjet[10];
   Int_t jetIdCutBased_wpjet[10];
   Int_t jetIdSimple_wpjet[10];
   Int_t jetIdFull_wpjet[10];
   Int_t assjet[10];
   Int_t partPdgIDjet[10];
   Int_t partMomPdgIDjet[10];

   Float_t deltaeta;
   Float_t zeppenjet;
   Float_t deltaphi;
   Float_t deltaphinewvtx;
   Float_t deltaphigg;
   Float_t eta2j;
   Float_t phi2j;
   Float_t pt2j;
   Float_t invmassjet;
   Float_t invmass2g1j;
   Float_t invmass2g2j;
   Float_t pt2g2j;           
   Float_t nvtx;
   Int_t vtxId;
   Float_t vtxPos_x;
   Float_t vtxPos_y;
   Float_t vtxPos_z;
   Float_t vtxIdMVA;
   Float_t vtxIdEvtProb;

   //////////////////////////////////////
   Float_t         sMet_  ;
   Float_t         eMet_  ;
   Float_t         phiMet_;
   Float_t         signifMet_;
   Float_t         eSmearedMet_;   
   Float_t         phiSmearedMet_;
   Float_t         eShiftedMet_;   
   Float_t         phiShiftedMet_;
   Float_t         eShiftedScaledMet_;   
   Float_t         phiShiftedScaledMet_;
   Float_t         eSmearedShiftedMet_;   
   Float_t         phiSmearedShiftedMet_;
   Float_t         eShiftedScaledMetPUcorr_;   
   Float_t         phiShiftedScaledMetPUcorr_;
   Float_t         eSmearedShiftedMetPUcorr_;   
   Float_t         phiSmearedShiftedMetPUcorr_;
   Float_t         sCorrMet_  ;
   Float_t         eCorrMet_  ;
   Float_t         phiCorrMet_;
   Float_t         signifCorrMet_;
   Float_t         smuCorrMet_  ;
   Float_t         emuCorrMet_  ;
   Float_t         phimuCorrMet_;
   Float_t         signifmuCorrMet_;
   Float_t         sNoHFMet_  ;
   Float_t         eNoHFMet_  ;
   Float_t         phiNoHFMet_;
   Float_t         signifNoHFMet_;
   Float_t         stcMet_  ;
   Float_t         etcMet_  ;
   Float_t         phitcMet_;
   Float_t         signiftcMet_;
   Float_t         sglobalPfMet_;
   Float_t         eglobalPfMet_;
   Float_t         phiglobalPfMet_;
   Float_t         signifglobalPfMet_;
   Float_t         scentralPfMet_;
   Float_t         ecentralPfMet_;
   Float_t         phicentralPfMet_;
   Float_t         signifcentralPfMet_;
   Float_t         eassocPfMet_;   //[nvertex]
   Float_t         phiassocPfMet_;   //[nvertex]
   Float_t         signifassocPfMet_;   //[nvertex]
   Float_t         eassocOtherVtxPfMet_;   //[nvertex]
   Float_t         phiassocOtherVtxPfMet_;   //[nvertex]
   Float_t         signifassocOtherVtxPfMet_;   //[nvertex]
   Float_t         etrkPfMet_;   //[nvertex]
   Float_t         phitrkPfMet_;   //[nvertex]
   Float_t         signiftrkPfMet_;   //[nvertex]
   Float_t         ecleanPfMet_;   //[nvertex]
   Float_t         phicleanPfMet_;   //[nvertex]
   Float_t         signifcleanPfMet_;   //[nvertex]
   Float_t         ecleanedSaclayPfMet_;   //[nvertex] 
   Float_t         phicleanedSaclayPfMet_;   //[nvertex] 
   Float_t         signifcleanedSaclayPfMet_;   //[nvertex] 
   Float_t         eminTypeICleanSaclayPfMet_;   //[nvertex] 
   Float_t         phiminTypeICleanSaclayPfMet_;   //[nvertex] 
   Float_t         signifminTypeICleanSaclayPfMet_;   //[nvertex]
   Float_t         globalPfSums_;
   Float_t         spfMet_  ;
   Float_t         epfMet_  ;
   Float_t         phipfMet_;
   Float_t         signifpfMet_;
   Float_t         spfMetType1_;
   Float_t         epfMetType1_;
   Float_t         phipfMetType1_;
   Float_t         signifpfMetType1_;
   Float_t         sMetGen_  ;
   Float_t         eMetGen_  ;
   Float_t         phiMetGen_;
   Float_t         signifMetGen_;
   Float_t         sMetGen2_  ;
   Float_t         eMetGen2_  ;
   Float_t         phiMetGen2_;
   //////////////////////////////////////

   // gen variables
   ///////////////////////
    
   Int_t gen_custom_processId;

   Float_t gen_pt_gamma1;
   Float_t gen_pt_gamma2;
   Float_t gen_eta_gamma1;
   Float_t gen_eta_gamma2;
   Float_t gen_phi_gamma1;
   Float_t gen_phi_gamma2;

   Float_t gen_pt_genjet1;
   Float_t gen_pt_genjet2;
   Float_t gen_eta_genjet1;
   Float_t gen_eta_genjet2;
   Float_t gen_phi_genjet1;
   Float_t gen_phi_genjet2;

  //  Float_t gen_pt_VectorBoson;
  //  Float_t gen_phi_VectorBoson;
  //  Float_t gen_eta_VectorBoson;

   Float_t gen_mass_diphoton;
   Float_t gen_pt_diphoton;
   Float_t gen_eta_diphoton;
   Float_t gen_phi_diphoton;

   Float_t gen_mass_dijet;
   Float_t gen_pt_dijet;
   Float_t gen_eta_dijet;
   Float_t gen_phi_dijet;

   Float_t gen_zeppenfeld;

   Float_t gen_pt_lep1,  gen_pt_lep2;
   Float_t gen_eta_lep1, gen_eta_lep2;
   Float_t gen_phi_lep1, gen_phi_lep2;
   Int_t gen_pid_lep1, gen_pid_lep2;
   ////////////////////////

   Int_t npu;
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
   Int_t idloosenewpuEGphot1;
   Int_t idloosenewpuEGphot2;
   Int_t idtightnewpuEGphot1;
   Int_t idtightnewpuEGphot2;
   Int_t idhggtightnewpuEGphot1;
   Int_t idhggtightnewpuEGphot2;
   Float_t idmvaphot1;
   Float_t idmvaphot2;
   Int_t idcicphot1;
   Int_t idcicphot2;
   Int_t idcicnoelvetophot1;
   Int_t idcicnoelvetophot2;
   Int_t idcicpfphot1;
   Int_t idcicpfphot2;
   Int_t idcicpfnoelvetophot1;
   Int_t idcicpfnoelvetophot2;
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
   Int_t promptGamma;
   Int_t LOGamma;
   Int_t ISRGamma;
   Int_t FSRGamma;

   Bool_t H_event;
   Bool_t V_event;
   Bool_t Zbb_event;
   Bool_t Vqq_event;

   Float_t   rhoPFRN;
   Float_t   pid_hlwTrackNoDzphot1;
   Float_t   pid_hlwTrackNoDzphot2;
   Int_t     pid_hasMatchedConvphot1;
   Int_t     pid_hasMatchedConvphot2;
   Int_t     pid_hasMatchedPromptElephot1;
   Int_t     pid_hasMatchedPromptElephot2;
   Float_t   r9phot1;
   Float_t   r9phot2;
   Float_t   etascphot1;
   Float_t   etascphot2;
   Float_t   phiscphot1;
   Float_t   phiscphot2;
   Float_t   pu_weight;
   Float_t   pt_weight;

   Int_t nWeightsPDF1;
   Int_t nWeightsPDF2;
   Int_t nWeightsPDF3;
   Int_t nWeightsPDF4;
   Int_t nWeightsPDF5;
   Int_t nWeightsPDF6;
   Int_t nWeightsPDF7;
   Int_t nWeightsPDF8;
   Int_t nWeightsPDF9;
   Int_t nWeightsPDF10;
   Float_t PDFweight1[150];
   Float_t PDFweight2[150];
   Float_t PDFweight3[150];
   Float_t PDFweight4[150];
   Float_t PDFweight5[150];
   Float_t PDFweight6[150];
   Float_t PDFweight7[150];
   Float_t PDFweight8[150];
   Float_t PDFweight9[150];
   Float_t PDFweight10[150];

   // lepton tag
   Float_t ptele1, ptele2;
   Float_t etaele1, etaele2;
   Float_t phiele1, phiele2;
   Float_t eneele1, eneele2;
   Float_t sIeIeele1, sIeIeele2;
   Float_t dphiele1, dphiele2;
   Float_t detaele1, detaele2;
   Int_t mhitsele1, mhitsele2;
   Float_t dcotele1, dcotele2;
   Float_t distele1, distele2;
   Float_t d0ele1, d0ele2;
   Float_t dzele1, dzele2;
   Float_t isoele1, isoele2;
   Float_t fullisoele1, fullisoele2;
   Float_t invMassele1g1, invMassele1g2;
   Float_t invMassele2g1, invMassele2g2;
   Float_t oEmoPele1, mvanotrigele1, mvatrigele1; 
   Int_t matchconvele1;
   Float_t chHadIso03ele1, nHadIso03ele1, photIso03ele1;
   Float_t oEmoPele2, mvanotrigele2, mvatrigele2; 
   Int_t matchconvele2;
   Float_t chHadIso03ele2, nHadIso03ele2, photIso03ele2;

   // loose electrons 
   Float_t pteleloose1, pteleloose2;
   Float_t etaeleloose1, etaeleloose2;
   Float_t phieleloose1, phieleloose2;
   Float_t eneeleloose1, eneeleloose2;
   Float_t sIeIeeleloose1, sIeIeeleloose2;
   Float_t dphieleloose1, dphieleloose2;
   Float_t detaeleloose1, detaeleloose2;
   Int_t mhitseleloose1, mhitseleloose2;
   Float_t dcoteleloose1, dcoteleloose2;
   Float_t disteleloose1, disteleloose2;
   Float_t d0eleloose1, d0eleloose2;
   Float_t dzeleloose1, dzeleloose2;
   Float_t isoeleloose1, isoeleloose2;
   Float_t fullisoeleloose1, fullisoeleloose2;
   Float_t invMasseleloose1g1, invMasseleloose1g2;
   Float_t invMasseleloose2g1, invMasseleloose2g2;
   Float_t oEmoPeleloose1, mvanotrigeleloose1, mvatrigeleloose1; 
   Int_t matchconveleloose1;
   Float_t chHadIso03eleloose1, nHadIso03eleloose1, photIso03eleloose1;
   Float_t oEmoPeleloose2, mvanotrigeleloose2, mvatrigeleloose2; 
   Int_t matchconveleloose2;
   Float_t chHadIso03eleloose2, nHadIso03eleloose2, photIso03eleloose2;

   Float_t ptmu1, ptmu2;
   Float_t etamu1, etamu2;
   Float_t phimu1, phimu2;
   Float_t enemu1, enemu2;
   Int_t pixhitsmu1, pixhitsmu2;
   Int_t trkhitsmu1, trkhitsmu2;
   Int_t hitsmu1, hitsmu2;
   Float_t chi2mu1, chi2mu2;
   Int_t matchmu1, matchmu2;
   Float_t d0mu1, d0mu2;
   Float_t dzmu1, dzmu2;
   Float_t isomu1,isomu2;
   Float_t chHadmu1, nHadmu1, photmu1, puptmu1;
   Float_t chHadmu2, nHadmu2, photmu2, puptmu2;

   // loose muons
   Float_t ptmuloose1, ptmuloose2;
   Float_t etamuloose1, etamuloose2;
   Float_t phimuloose1, phimuloose2;
   Float_t enemuloose1, enemuloose2;
   Int_t pixhitsmuloose1, pixhitsmuloose2;
   Int_t trkhitsmuloose1, trkhitsmuloose2;
   Int_t hitsmuloose1, hitsmuloose2;
   Float_t chi2muloose1, chi2muloose2;
   Int_t matchmuloose1, matchmuloose2;
   Float_t d0muloose1, d0muloose2;
   Float_t dzmuloose1, dzmuloose2;
   Float_t isomuloose1,isomuloose2;
   Float_t chHadmuloose1, nHadmuloose1, photmuloose1, puptmuloose1;
   Float_t chHadmuloose2, nHadmuloose2, photmuloose2, puptmuloose2;

   
   float weight;
};
#endif
