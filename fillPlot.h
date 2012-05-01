//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Wed Apr 11 21:12:57 2012 by ROOT version 5.30/02
// from TTree AnaTree/Reduced tree for final analysis
// found on file: ../../dati/Hgg/redntp.42xv6b_data.cicloose.regrPho_eCorr_30Nov.v2_jetid/merged/redntp_Photon-Run2011-30Nov2011-v1-DiPhotonSkimOnFly.root
//////////////////////////////////////////////////////////

#ifndef fillPlot_h
#define fillPlot_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TH1D.h>

class fillPlot {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

   struct diPhotonTree_structure_ {
     int run;
     int lumi;
     int event;
     float ptgg;
     int ebeb;
     float massggnewvtx;
     float weight;
   };

   diPhotonTree_structure_ tree_;

   //Cuts values                                                                
   double ptphot1cut;
   double ptphot2cut;
   double pthiggsmincut;
   double pthiggsmaxcut;
   double ptjet1cut;
   double ptjet2cut;
   double metcut;
   double phijetmetcut;
   double phiphot1metcut;
   double phiphot2metcut;
   double phihiggsmetcut;
   double phipho1pho2cut;
   double deltaetacut;
   double deltaphicut;
   double zeppencut;
   double invmassjetcut;
   int ebcat;
   int r9cat;
   int cicselection;
   bool thirdcat;
   bool leptontag;
   bool leptonveto;
   double ptphot1cut_2;
   double ptphot2cut_2;
   double pthiggsmincut_2;
   double pthiggsmaxcut_2;

   // bool to decide if we want to write output txt file 
   std::string writetxt;
   std::string writeRoot;

   // bool for switching on smearing and smearing parameters                    
   bool dosmear,dopureweight,doptreweight;
   double meansmear, spreadsmear;

   // vector of pu weights                                                      
   std::vector<Double_t> puweights_;
   double weights_[15][15];

   // Declaration of leaf types
   Int_t           run;
   Int_t           event;
   Int_t           lumi;
   Float_t         rhoPF;
   Float_t         massgg;
   Float_t         ptgg;
   Float_t         ptggnewvtx;
   Float_t         phigg;
   Float_t         etagg;
   Float_t         massggnewvtx;
   Float_t         ptphot1;
   Float_t         ptphot2;
   Float_t         deltaRToTrackphot1;
   Float_t         deltaRToTrackphot2;
   Float_t         timephot1;
   Float_t         timephot2;
   Float_t         etaphot1;
   Float_t         etaphot2;
   Float_t         phiphot1;
   Float_t         phiphot2;
   Float_t         etascphot1;
   Float_t         etascphot2;
   Float_t         phiscphot1;
   Float_t         phiscphot2;
   Float_t         E1phot1;
   Float_t         E1phot2;
   Float_t         E9phot1;
   Float_t         E9phot2;
   Float_t         r9phot1;
   Float_t         r9phot2;
   Int_t           isemEGphot1;
   Int_t           isemEGphot2;
   Int_t           promptGamma;
   Int_t           LOGamma;
   Int_t           ISRGamma;
   Int_t           FSRGamma;
   Int_t           idcicphot1;
   Int_t           idcicphot2;
   Int_t           idcicnoelvetophot1;
   Int_t           idcicnoelvetophot2;
   Int_t           idelephot1;
   Int_t           idelephot2;
   Int_t           pid_isEMphot1;
   Int_t           pid_isEMphot2;
   Int_t           pid_haspixelseedphot1;
   Int_t           pid_haspixelseedphot2;
   Float_t         pid_jurECALphot1;
   Float_t         pid_jurECALphot2;
   Float_t         pid_twrHCALphot1;
   Float_t         pid_twrHCALphot2;
   Float_t         pid_HoverEphot1;
   Float_t         pid_HoverEphot2;
   Float_t         pid_hlwTrackphot1;
   Float_t         pid_hlwTrackphot2;
   Float_t         pid_etawidphot1;
   Float_t         pid_etawidphot2;
   Float_t         pid_hlwTrackNoDzphot1;
   Float_t         pid_hlwTrackNoDzphot2;
   Int_t           pid_hasMatchedConvphot1;
   Int_t           pid_hasMatchedConvphot2;
   Int_t           pid_hasMatchedPromptElephot1;
   Int_t           pid_hasMatchedPromptElephot2;
   Float_t         pid_sminphot1;
   Float_t         pid_sminphot2;
   Float_t         pid_smajphot1;
   Float_t         pid_smajphot2;
   Int_t           pid_ntrkphot1;
   Int_t           pid_ntrkphot2;
   Float_t         pid_ptisophot1;
   Float_t         pid_ptisophot2;
   Int_t           pid_ntrkcsphot1;
   Int_t           pid_ntrkcsphot2;
   Float_t         pid_ptisocsphot1;
   Float_t         pid_ptisocsphot2;
   Float_t         pid_ecalisophot1;
   Float_t         pid_ecalisophot2;
   Float_t         pid_hcalisophot1;
   Float_t         pid_hcalisophot2;
   Float_t         ptjet1;
   Float_t         ptjet2;
   Float_t         ptjet3;
   Float_t         ptjet4;
   Float_t         ptcorrjet1;
   Float_t         ptcorrjet2;
   Float_t         ptcorrjet3;
   Float_t         ptcorrjet4;
   Float_t         etajet1;
   Float_t         etajet2;
   Float_t         etajet3;
   Float_t         etajet4;
   Float_t         phijet1;
   Float_t         phijet2;
   Float_t         phijet3;
   Float_t         phijet4;
   Float_t         betajet1;
   Float_t         betajet2;
   Float_t         betastarjet1;
   Float_t         betastarjet2;
   Float_t         btagvtxjet1;
   Float_t         btagtrkjet1;
   Float_t         btagvtxjet2;
   Float_t         btagtrkjet2;
   Float_t         ptDjet1;
   Float_t         rmsjet1;
   Int_t           ntrkjet1;
   Int_t           nneutjet1;
   Float_t         ptDjet2;
   Float_t         rmsjet2;
   Int_t           ntrkjet2;
   Int_t           nneutjet2;
   Int_t           assjet1;
   Int_t           assjet2;
   Float_t         deltaeta;
   Float_t         zeppenjet;
   Float_t         deltaphi;
   Float_t         deltaphinewvtx;
   Float_t         deltaphigg;
   Float_t         invmassjet;
   Float_t         invmass2g1j;
   Float_t         invmass2g2j;
   Float_t         pt2g2j;
   Float_t         eta2j;
   Float_t         phi2j;
   Float_t         pt2j;
   Float_t         nvtx;
   Float_t         sMet;
   Float_t         eMet;
   Float_t         phiMet;
   Float_t         signifMet;
   Float_t         eSmearedMet;
   Float_t         phiSmearedMet;
   Float_t         eShiftedMet;
   Float_t         phiShiftedMet;
   Float_t         eShiftedScaledMet;
   Float_t         phiShiftedScaledMet;
   Float_t         eSmearedShiftedMet;
   Float_t         phiSmearedShiftedMet;
   Float_t         sCorrMet;
   Float_t         eCorrMet;
   Float_t         phiCorrMet;
   Float_t         signifCorrMet;
   Float_t         smuCorrMet;
   Float_t         emuCorrMet;
   Float_t         phimuCorrMet;
   Float_t         signifmuCorrMet;
   Float_t         sNoHFMet;
   Float_t         eNoHFMet;
   Float_t         phiNoHFMet;
   Float_t         signifNoHFMet;
   Float_t         stcMet;
   Float_t         etcMet;
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
   Float_t         eassocPfMet;
   Float_t         phiassocPfMet;
   Float_t         signifassocPfMet;
   Float_t         eassocOtherVtxPfMet;
   Float_t         phiassocOtherVtxPfMet;
   Float_t         signifassocOtherVtxPfMet;
   Float_t         etrkPfMet;
   Float_t         phitrkPfMet;
   Float_t         signiftrkPfMet;
   Float_t         ecleanPfMet;
   Float_t         phicleanPfMet;
   Float_t         signifcleanPfMet;
   Float_t         ecleanedSaclayPfMet;
   Float_t         phicleanedSaclayPfMet;
   Float_t         signifcleanedSaclayPfMet;
   Float_t         eminTypeICleanSaclayPfMet;
   Float_t         phiminTypeICleanSaclayPfMet;
   Float_t         signifminTypeICleanSaclayPfMet;
   Float_t         globalPfSums;
   Float_t         spfMet;
   Float_t         epfMet;
   Float_t         phipfMet;
   Float_t         signifpfMet;
   Float_t         spfMetType1;
   Float_t         epfMetType1;
   Float_t         phipfMetType1;
   Float_t         signifpfMetType1;
   Float_t         sMetGen;
   Float_t         eMetGen;
   Float_t         phiMetGen;
   Float_t         signifMetGen;
   Float_t         sMetGen2;
   Float_t         eMetGen2;
   Float_t         phiMetGen2;
   Int_t           npu;
   Int_t           NtotEvents;
   Float_t         xsection;
   Float_t         EquivLumi;
   Int_t           SampleID;
   Float_t         pu_weight;
   Float_t         pt_weight;
   Int_t           gen_custom_processId;
   Float_t         gen_pt_gamma1;
   Float_t         gen_pt_gamma2;
   Float_t         gen_eta_gamma1;
   Float_t         gen_eta_gamma2;
   Float_t         gen_phi_gamma1;
   Float_t         gen_phi_gamma2;
   Float_t         gen_pt_genjet1;
   Float_t         gen_pt_genjet2;
   Float_t         gen_eta_genjet1;
   Float_t         gen_eta_genjet2;
   Float_t         gen_phi_genjet1;
   Float_t         gen_phi_genjet2;
   Float_t         gen_mass_diphoton;
   Float_t         gen_pt_diphoton;
   Float_t         gen_eta_diphoton;
   Float_t         gen_phi_diphoton;
   Float_t         gen_mass_dijet;
   Float_t         gen_pt_dijet;
   Float_t         gen_eta_dijet;
   Float_t         gen_phi_dijet;
   Float_t         gen_zeppenfeld;
   Float_t         gen_pt_lep1;
   Float_t         gen_pt_lep2;
   Float_t         gen_eta_lep1;
   Float_t         gen_eta_lep2;
   Float_t         gen_phi_lep1;
   Float_t         gen_phi_lep2;
   Int_t           gen_pid_lep1;
   Int_t           gen_pid_lep2;
   Float_t         ptele1;
   Float_t         ptele2;
   Float_t         etaele1;
   Float_t         etaele2;
   Float_t         phiele1;
   Float_t         phiele2;
   Float_t         eneele1;
   Float_t         eneele2;
   Float_t         sIeIeele1;
   Float_t         sIeIeele2;
   Float_t         dphiele1;
   Float_t         dphiele2;
   Float_t         detaele1;
   Float_t         detaele2;
   Int_t           mhitsele1;
   Int_t           mhitsele2;
   Float_t         dcotele1;
   Float_t         dcotele2;
   Float_t         distele1;
   Float_t         distele2;
   Float_t         d0ele1;
   Float_t         d0ele2;
   Float_t         dzele1;
   Float_t         dzele2;
   Float_t         isoele1;
   Float_t         isoele2;
   Float_t         fullisoele1;
   Float_t         fullisoele2;
   Float_t         invMassele1g1;
   Float_t         invMassele1g2;
   Float_t         invMassele2g1;
   Float_t         invMassele2g2;
   Float_t         ptmu1;
   Float_t         ptmu2;
   Float_t         etamu1;
   Float_t         etamu2;
   Float_t         phimu1;
   Float_t         phimu2;
   Float_t         enemu1;
   Float_t         enemu2;
   Int_t           pixhitsmu1;
   Int_t           pixhitsmu2;
   Int_t           trkhitsmu1;
   Int_t           trkhitsmu2;
   Int_t           hitsmu1;
   Int_t           hitsmu2;
   Float_t         chi2mu1;
   Float_t         chi2mu2;
   Int_t           matchmu1;
   Int_t           matchmu2;
   Float_t         d0mu1;
   Float_t         d0mu2;
   Float_t         dzmu1;
   Float_t         dzmu2;
   Float_t         isomu1;
   Float_t         isomu2;
   Int_t           nWeightsPDF1;
   Int_t           nWeightsPDF2;
   Int_t           nWeightsPDF3;
   Int_t           nWeightsPDF4;
   Int_t           nWeightsPDF5;
   Int_t           nWeightsPDF6;
   Int_t           nWeightsPDF7;
   Int_t           nWeightsPDF8;
   Int_t           nWeightsPDF9;
   Int_t           nWeightsPDF10;
   Float_t         PDFweight1[1];   //[nWeightsPDF1]
   Float_t         PDFweight2[1];   //[nWeightsPDF2]
   Float_t         PDFweight3[1];   //[nWeightsPDF3]
   Float_t         PDFweight4[1];   //[nWeightsPDF4]
   Float_t         PDFweight5[1];   //[nWeightsPDF5]
   Float_t         PDFweight6[1];   //[nWeightsPDF6]
   Float_t         PDFweight7[1];   //[nWeightsPDF7]
   Float_t         PDFweight8[1];   //[nWeightsPDF8]
   Float_t         PDFweight9[1];   //[nWeightsPDF9]
   Float_t         PDFweight10[1];   //[nWeightsPDF10]

   // List of branches
   TBranch        *b_run;   //!
   TBranch        *b_event;   //!
   TBranch        *b_lumi;   //!
   TBranch        *b_rhoPF;   //!
   TBranch        *b_massgg;   //!
   TBranch        *b_ptgg;   //!
   TBranch        *b_ptggnewvtx;   //!
   TBranch        *b_phigg;   //!
   TBranch        *b_etagg;   //!
   TBranch        *b_massggnewvtx;   //!
   TBranch        *b_ptphot1;   //!
   TBranch        *b_ptphot2;   //!
   TBranch        *b_deltaRToTrackphot1;   //!
   TBranch        *b_deltaRToTrackphot2;   //!
   TBranch        *b_timephot1;   //!
   TBranch        *b_timephot2;   //!
   TBranch        *b_etaphot1;   //!
   TBranch        *b_etaphot2;   //!
   TBranch        *b_phiphot1;   //!
   TBranch        *b_phiphot2;   //!
   TBranch        *b_etascphot1;   //!
   TBranch        *b_etascphot2;   //!
   TBranch        *b_phiscphot1;   //!
   TBranch        *b_phiscphot2;   //!
   TBranch        *b_E1phot1;   //!
   TBranch        *b_E1phot2;   //!
   TBranch        *b_E9phot1;   //!
   TBranch        *b_E9phot2;   //!
   TBranch        *b_r9phot1;   //!
   TBranch        *b_r9phot2;   //!
   TBranch        *b_isemEGphot1;   //!
   TBranch        *b_isemEGphot2;   //!
   TBranch        *b_promptGamma;   //!
   TBranch        *b_LOGamma;   //!
   TBranch        *b_ISRGamma;   //!
   TBranch        *b_FSRGamma;   //!
   TBranch        *b_idcicphot1;   //!
   TBranch        *b_idcicphot2;   //!
   TBranch        *b_idcicnoelvetophot1;   //!
   TBranch        *b_idcicnoelvetophot2;   //!
   TBranch        *b_idelephot1;   //!
   TBranch        *b_idelephot2;   //!
   TBranch        *b_pid_isEMphot1;   //!
   TBranch        *b_pid_isEMphot2;   //!
   TBranch        *b_pid_haspixelseedphot1;   //!
   TBranch        *b_pid_haspixelseedphot2;   //!
   TBranch        *b_pid_jurECALphot1;   //!
   TBranch        *b_pid_jurECALphot2;   //!
   TBranch        *b_pid_twrHCALphot1;   //!
   TBranch        *b_pid_twrHCALphot2;   //!
   TBranch        *b_pid_HoverEphot1;   //!
   TBranch        *b_pid_HoverEphot2;   //!
   TBranch        *b_pid_hlwTrackphot1;   //!
   TBranch        *b_pid_hlwTrackphot2;   //!
   TBranch        *b_pid_etawidphot1;   //!
   TBranch        *b_pid_etawidphot2;   //!
   TBranch        *b_pid_hlwTrackNoDzphot1;   //!
   TBranch        *b_pid_hlwTrackNoDzphot2;   //!
   TBranch        *b_pid_hasMatchedConvphot1;   //!
   TBranch        *b_pid_hasMatchedConvphot2;   //!
   TBranch        *b_pid_hasMatchedPromptElephot1;   //!
   TBranch        *b_pid_hasMatchedPromptElephot2;   //!
   TBranch        *b_pid_sminphot1;   //!
   TBranch        *b_pid_sminphot2;   //!
   TBranch        *b_pid_smajphot1;   //!
   TBranch        *b_pid_smajphot2;   //!
   TBranch        *b_pid_ntrkphot1;   //!
   TBranch        *b_pid_ntrkphot2;   //!
   TBranch        *b_pid_ptisophot1;   //!
   TBranch        *b_pid_ptisophot2;   //!
   TBranch        *b_pid_ntrkcsphot1;   //!
   TBranch        *b_pid_ntrkcsphot2;   //!
   TBranch        *b_pid_ptisocsphot1;   //!
   TBranch        *b_pid_ptisocsphot2;   //!
   TBranch        *b_pid_ecalisophot1;   //!
   TBranch        *b_pid_ecalisophot2;   //!
   TBranch        *b_pid_hcalisophot1;   //!
   TBranch        *b_pid_hcalisophot2;   //!
   TBranch        *b_ptjet1;   //!
   TBranch        *b_ptjet2;   //!
   TBranch        *b_ptjet3;   //!
   TBranch        *b_ptjet4;   //!
   TBranch        *b_ptcorrjet1;   //!
   TBranch        *b_ptcorrjet2;   //!
   TBranch        *b_ptcorrjet3;   //!
   TBranch        *b_ptcorrjet4;   //!
   TBranch        *b_etajet1;   //!
   TBranch        *b_etajet2;   //!
   TBranch        *b_etajet3;   //!
   TBranch        *b_etajet4;   //!
   TBranch        *b_phijet1;   //!
   TBranch        *b_phijet2;   //!
   TBranch        *b_phijet3;   //!
   TBranch        *b_phijet4;   //!
   TBranch        *b_betajet1;   //!
   TBranch        *b_betajet2;   //!
   TBranch        *b_betastarjet1;   //!
   TBranch        *b_betastarjet2;   //!
   TBranch        *b_btagvtxjet1;   //!
   TBranch        *b_btagtrkjet1;   //!
   TBranch        *b_btagvtxjet2;   //!
   TBranch        *b_btagtrkjet2;   //!
   TBranch        *b_ptDjet1;   //!
   TBranch        *b_rmsjet1;   //!
   TBranch        *b_ntrkjet1;   //!
   TBranch        *b_nneutjet1;   //!
   TBranch        *b_ptDjet2;   //!
   TBranch        *b_rmsjet2;   //!
   TBranch        *b_ntrkjet2;   //!
   TBranch        *b_nneutjet2;   //!
   TBranch        *b_assjet1;   //!
   TBranch        *b_assjet2;   //!
   TBranch        *b_deltaeta;   //!
   TBranch        *b_zeppenjet;   //!
   TBranch        *b_deltaphi;   //!
   TBranch        *b_deltaphinewvtx;   //!
   TBranch        *b_deltaphigg;   //!
   TBranch        *b_invmassjet;   //!
   TBranch        *b_invmass2g1j;   //!
   TBranch        *b_invmass2g2j;   //!
   TBranch        *b_pt2g2j;   //!
   TBranch        *b_eta2j;   //!
   TBranch        *b_phi2j;   //!
   TBranch        *b_pt2j;   //!
   TBranch        *b_nvtx;   //!
   TBranch        *b_sMet;   //!
   TBranch        *b_eMet;   //!
   TBranch        *b_phiMet;   //!
   TBranch        *b_signifMet;   //!
   TBranch        *b_eSmearedMet;   //!
   TBranch        *b_phiSmearedMet;   //!
   TBranch        *b_eShiftedMet;   //!
   TBranch        *b_phiShiftedMet;   //!
   TBranch        *b_eShiftedScaledMet;   //!
   TBranch        *b_phiShiftedScaledMet;   //!
   TBranch        *b_eSmearedShiftedMet;   //!
   TBranch        *b_phiSmearedShiftedMet;   //!
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
   TBranch        *b_npu;   //!
   TBranch        *b_NtotEvents;   //!
   TBranch        *b_xsection;   //!
   TBranch        *b_EquivLumi;   //!
   TBranch        *b_SampleID;   //!
   TBranch        *b_pu_weight;   //!
   TBranch        *b_pt_weight;   //!
   TBranch        *b_gen_custom_processId;   //!
   TBranch        *b_gen_pt_gamma1;   //!
   TBranch        *b_gen_pt_gamma2;   //!
   TBranch        *b_gen_eta_gamma1;   //!
   TBranch        *b_gen_eta_gamma2;   //!
   TBranch        *b_gen_phi_gamma1;   //!
   TBranch        *b_gen_phi_gamma2;   //!
   TBranch        *b_gen_pt_genjet1;   //!
   TBranch        *b_gen_pt_genjet2;   //!
   TBranch        *b_gen_eta_genjet1;   //!
   TBranch        *b_gen_eta_genjet2;   //!
   TBranch        *b_gen_phi_genjet1;   //!
   TBranch        *b_gen_phi_genjet2;   //!
   TBranch        *b_gen_mass_diphoton;   //!
   TBranch        *b_gen_pt_diphoton;   //!
   TBranch        *b_gen_eta_diphoton;   //!
   TBranch        *b_gen_phi_diphoton;   //!
   TBranch        *b_gen_mass_dijet;   //!
   TBranch        *b_gen_pt_dijet;   //!
   TBranch        *b_gen_eta_dijet;   //!
   TBranch        *b_gen_phi_dijet;   //!
   TBranch        *b_gen_zeppenfeld;   //!
   TBranch        *b_gen_pt_lep1;   //!
   TBranch        *b_gen_pt_lep2;   //!
   TBranch        *b_gen_eta_lep1;   //!
   TBranch        *b_gen_eta_lep2;   //!
   TBranch        *b_gen_phi_lep1;   //!
   TBranch        *b_gen_phi_lep2;   //!
   TBranch        *b_gen_pid_lep1;   //!
   TBranch        *b_gen_pid_lep2;   //!
   TBranch        *b_ptele1;   //!
   TBranch        *b_ptele2;   //!
   TBranch        *b_etaele1;   //!
   TBranch        *b_etaele2;   //!
   TBranch        *b_phiele1;   //!
   TBranch        *b_phiele2;   //!
   TBranch        *b_eneele1;   //!
   TBranch        *b_eneele2;   //!
   TBranch        *b_sIeIeele1;   //!
   TBranch        *b_sIeIeele2;   //!
   TBranch        *b_dphiele1;   //!
   TBranch        *b_dphiele2;   //!
   TBranch        *b_detaele1;   //!
   TBranch        *b_detaele2;   //!
   TBranch        *b_mhitsele1;   //!
   TBranch        *b_mhitsele2;   //!
   TBranch        *b_dcotele1;   //!
   TBranch        *b_dcotele2;   //!
   TBranch        *b_distele1;   //!
   TBranch        *b_distele2;   //!
   TBranch        *b_d0ele1;   //!
   TBranch        *b_d0ele2;   //!
   TBranch        *b_dzele1;   //!
   TBranch        *b_dzele2;   //!
   TBranch        *b_isoele1;   //!
   TBranch        *b_isoele2;   //!
   TBranch        *b_fullisoele1;   //!
   TBranch        *b_fullisoele2;   //!
   TBranch        *b_invMassele1g1;   //!
   TBranch        *b_invMassele1g2;   //!
   TBranch        *b_invMassele2g1;   //!
   TBranch        *b_invMassele2g2;   //!
   TBranch        *b_ptmu1;   //!
   TBranch        *b_ptmu2;   //!
   TBranch        *b_etamu1;   //!
   TBranch        *b_etamu2;   //!
   TBranch        *b_phimu1;   //!
   TBranch        *b_phimu2;   //!
   TBranch        *b_enemu1;   //!
   TBranch        *b_enemu2;   //!
   TBranch        *b_pixhitsmu1;   //!
   TBranch        *b_pixhitsmu2;   //!
   TBranch        *b_trkhitsmu1;   //!
   TBranch        *b_trkhitsmu2;   //!
   TBranch        *b_hitsmu1;   //!
   TBranch        *b_hitsmu2;   //!
   TBranch        *b_chi2mu1;   //!
   TBranch        *b_chi2mu2;   //!
   TBranch        *b_matchmu1;   //!
   TBranch        *b_matchmu2;   //!
   TBranch        *b_d0mu1;   //!
   TBranch        *b_d0mu2;   //!
   TBranch        *b_dzmu1;   //!
   TBranch        *b_dzmu2;   //!
   TBranch        *b_isomu1;   //!
   TBranch        *b_isomu2;   //!
   TBranch        *b_nWeightsPDF1;   //!
   TBranch        *b_nWeightsPDF2;   //!
   TBranch        *b_nWeightsPDF3;   //!
   TBranch        *b_nWeightsPDF4;   //!
   TBranch        *b_nWeightsPDF5;   //!
   TBranch        *b_nWeightsPDF6;   //!
   TBranch        *b_nWeightsPDF7;   //!
   TBranch        *b_nWeightsPDF8;   //!
   TBranch        *b_nWeightsPDF9;   //!
   TBranch        *b_nWeightsPDF10;   //!
   TBranch        *b_PDFweight1;   //!
   TBranch        *b_PDFweight2;   //!
   TBranch        *b_PDFweight3;   //!
   TBranch        *b_PDFweight4;   //!
   TBranch        *b_PDFweight5;   //!
   TBranch        *b_PDFweight6;   //!
   TBranch        *b_PDFweight7;   //!
   TBranch        *b_PDFweight8;   //!
   TBranch        *b_PDFweight9;   //!
   TBranch        *b_PDFweight10;   //!

   fillPlot(TTree *tree=0,bool isData=0);
   virtual ~fillPlot();
   virtual void     getweights();
   // virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);

   virtual void     Setcuts(double pt1=50, double pt2=30, double pthiggsmin=-100, double pthiggsmax=-100, double ptj1=20, double ptj2=15, double misset=70, double deltae=2.5, double zep=2.5, double mjj=300, double deltap=2.6, double jetmet=0., double p1met=0., double p2met=0., double hmet=0., double phigg=300., int eb = 1, int r9 = 1, bool thirdcat = 0, bool leptag = 0, bool lepveto = 0);
   virtual TH1D*    Plot(string var, string name, int nbin=200, double min=90, double max=190, bool cs=0, int signal=50);
   virtual bool     cutIDEG(double ptPhot, double etaPhot, double pid_hlwTrackNoDz, double pid_jurECAL, double pid_twrHCAL, double pid_HoverE, double pid_etawid, int scaletrk=100, int scaleecal=100, int scalehcal=100, int scalehove=100);
   virtual bool     exclSel();
   virtual void     setCic(int cic=5);
   virtual void     Writetxt(char * filename);
   virtual void     WriteRoot(char * filename);
   virtual void     SetPuWeights(bool isData = 0,std::string file = "");
   virtual void     DoSmearing(double mean, double spread);
   virtual void     DoPuReweight();
   virtual void     DoPtReweight();
   bool electronTagSelection();
   bool muonTagSelection();

   // virtual void     Loop();
   virtual Bool_t   Notify();
   // virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef fillPlot_cxx
fillPlot::fillPlot(TTree *tree,bool isData)
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("../../dati/Hgg/redntp.42xv6b_data.cicloose.regrPho_eCorr_30Nov.v2_jetid/merged/redntp_Photon-Run2011-30Nov2011-v1-DiPhotonSkimOnFly.root");
      if (!f || !f->IsOpen()) {
         f = new TFile("../../dati/Hgg/redntp.42xv6b_data.cicloose.regrPho_eCorr_30Nov.v2_jetid/merged/redntp_Photon-Run2011-30Nov2011-v1-DiPhotonSkimOnFly.root");
      }
      f->GetObject("AnaTree",tree);

   }
   Init(tree);
   writetxt = "";
   writeRoot = "";
   dosmear = 0;
   dopureweight = 0;
   doptreweight = 0;
   cicselection = -1;
   //   SetPuWeights(isData); 
}

fillPlot::~fillPlot()
{
  if (!fChain) return;
  // delete fChain->GetCurrentFile();
}

Int_t fillPlot::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t fillPlot::LoadTree(Long64_t entry)
{
// Set the environment to read one entry
   if (!fChain) return -5;
   Long64_t centry = fChain->LoadTree(entry);
   if (centry < 0) return centry;
   if (fChain->GetTreeNumber() != fCurrent) {
      fCurrent = fChain->GetTreeNumber();
      Notify();
   }
   return centry;
}

void fillPlot::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("run", &run, &b_run);
   fChain->SetBranchAddress("event", &event, &b_event);
   fChain->SetBranchAddress("lumi", &lumi, &b_lumi);
   fChain->SetBranchAddress("rhoPF", &rhoPF, &b_rhoPF);
   fChain->SetBranchAddress("massgg", &massgg, &b_massgg);
   fChain->SetBranchAddress("ptgg", &ptgg, &b_ptgg);
   fChain->SetBranchAddress("ptggnewvtx", &ptggnewvtx, &b_ptggnewvtx);
   fChain->SetBranchAddress("phigg", &phigg, &b_phigg);
   fChain->SetBranchAddress("etagg", &etagg, &b_etagg);
   fChain->SetBranchAddress("massggnewvtx", &massggnewvtx, &b_massggnewvtx);
   fChain->SetBranchAddress("ptphot1", &ptphot1, &b_ptphot1);
   fChain->SetBranchAddress("ptphot2", &ptphot2, &b_ptphot2);
   fChain->SetBranchAddress("deltaRToTrackphot1", &deltaRToTrackphot1, &b_deltaRToTrackphot1);
   fChain->SetBranchAddress("deltaRToTrackphot2", &deltaRToTrackphot2, &b_deltaRToTrackphot2);
   fChain->SetBranchAddress("timephot1", &timephot1, &b_timephot1);
   fChain->SetBranchAddress("timephot2", &timephot2, &b_timephot2);
   fChain->SetBranchAddress("etaphot1", &etaphot1, &b_etaphot1);
   fChain->SetBranchAddress("etaphot2", &etaphot2, &b_etaphot2);
   fChain->SetBranchAddress("phiphot1", &phiphot1, &b_phiphot1);
   fChain->SetBranchAddress("phiphot2", &phiphot2, &b_phiphot2);
   fChain->SetBranchAddress("etascphot1", &etascphot1, &b_etascphot1);
   fChain->SetBranchAddress("etascphot2", &etascphot2, &b_etascphot2);
   fChain->SetBranchAddress("phiscphot1", &phiscphot1, &b_phiscphot1);
   fChain->SetBranchAddress("phiscphot2", &phiscphot2, &b_phiscphot2);
   fChain->SetBranchAddress("E1phot1", &E1phot1, &b_E1phot1);
   fChain->SetBranchAddress("E1phot2", &E1phot2, &b_E1phot2);
   fChain->SetBranchAddress("E9phot1", &E9phot1, &b_E9phot1);
   fChain->SetBranchAddress("E9phot2", &E9phot2, &b_E9phot2);
   fChain->SetBranchAddress("r9phot1", &r9phot1, &b_r9phot1);
   fChain->SetBranchAddress("r9phot2", &r9phot2, &b_r9phot2);
   fChain->SetBranchAddress("isemEGphot1", &isemEGphot1, &b_isemEGphot1);
   fChain->SetBranchAddress("isemEGphot2", &isemEGphot2, &b_isemEGphot2);
   fChain->SetBranchAddress("promptGamma", &promptGamma, &b_promptGamma);
   fChain->SetBranchAddress("LOGamma", &LOGamma, &b_LOGamma);
   fChain->SetBranchAddress("ISRGamma", &ISRGamma, &b_ISRGamma);
   fChain->SetBranchAddress("FSRGamma", &FSRGamma, &b_FSRGamma);
   fChain->SetBranchAddress("idcicphot1", &idcicphot1, &b_idcicphot1);
   fChain->SetBranchAddress("idcicphot2", &idcicphot2, &b_idcicphot2);
   fChain->SetBranchAddress("idcicnoelvetophot1", &idcicnoelvetophot1, &b_idcicnoelvetophot1);
   fChain->SetBranchAddress("idcicnoelvetophot2", &idcicnoelvetophot2, &b_idcicnoelvetophot2);
   fChain->SetBranchAddress("idelephot1", &idelephot1, &b_idelephot1);
   fChain->SetBranchAddress("idelephot2", &idelephot2, &b_idelephot2);
   fChain->SetBranchAddress("pid_isEMphot1", &pid_isEMphot1, &b_pid_isEMphot1);
   fChain->SetBranchAddress("pid_isEMphot2", &pid_isEMphot2, &b_pid_isEMphot2);
   fChain->SetBranchAddress("pid_haspixelseedphot1", &pid_haspixelseedphot1, &b_pid_haspixelseedphot1);
   fChain->SetBranchAddress("pid_haspixelseedphot2", &pid_haspixelseedphot2, &b_pid_haspixelseedphot2);
   fChain->SetBranchAddress("pid_jurECALphot1", &pid_jurECALphot1, &b_pid_jurECALphot1);
   fChain->SetBranchAddress("pid_jurECALphot2", &pid_jurECALphot2, &b_pid_jurECALphot2);
   fChain->SetBranchAddress("pid_twrHCALphot1", &pid_twrHCALphot1, &b_pid_twrHCALphot1);
   fChain->SetBranchAddress("pid_twrHCALphot2", &pid_twrHCALphot2, &b_pid_twrHCALphot2);
   fChain->SetBranchAddress("pid_HoverEphot1", &pid_HoverEphot1, &b_pid_HoverEphot1);
   fChain->SetBranchAddress("pid_HoverEphot2", &pid_HoverEphot2, &b_pid_HoverEphot2);
   fChain->SetBranchAddress("pid_hlwTrackphot1", &pid_hlwTrackphot1, &b_pid_hlwTrackphot1);
   fChain->SetBranchAddress("pid_hlwTrackphot2", &pid_hlwTrackphot2, &b_pid_hlwTrackphot2);
   fChain->SetBranchAddress("pid_etawidphot1", &pid_etawidphot1, &b_pid_etawidphot1);
   fChain->SetBranchAddress("pid_etawidphot2", &pid_etawidphot2, &b_pid_etawidphot2);
   fChain->SetBranchAddress("pid_hlwTrackNoDzphot1", &pid_hlwTrackNoDzphot1, &b_pid_hlwTrackNoDzphot1);
   fChain->SetBranchAddress("pid_hlwTrackNoDzphot2", &pid_hlwTrackNoDzphot2, &b_pid_hlwTrackNoDzphot2);
   fChain->SetBranchAddress("pid_hasMatchedConvphot1", &pid_hasMatchedConvphot1, &b_pid_hasMatchedConvphot1);
   fChain->SetBranchAddress("pid_hasMatchedConvphot2", &pid_hasMatchedConvphot2, &b_pid_hasMatchedConvphot2);
   fChain->SetBranchAddress("pid_hasMatchedPromptElephot1", &pid_hasMatchedPromptElephot1, &b_pid_hasMatchedPromptElephot1);
   fChain->SetBranchAddress("pid_hasMatchedPromptElephot2", &pid_hasMatchedPromptElephot2, &b_pid_hasMatchedPromptElephot2);
   fChain->SetBranchAddress("pid_sminphot1", &pid_sminphot1, &b_pid_sminphot1);
   fChain->SetBranchAddress("pid_sminphot2", &pid_sminphot2, &b_pid_sminphot2);
   fChain->SetBranchAddress("pid_smajphot1", &pid_smajphot1, &b_pid_smajphot1);
   fChain->SetBranchAddress("pid_smajphot2", &pid_smajphot2, &b_pid_smajphot2);
   fChain->SetBranchAddress("pid_ntrkphot1", &pid_ntrkphot1, &b_pid_ntrkphot1);
   fChain->SetBranchAddress("pid_ntrkphot2", &pid_ntrkphot2, &b_pid_ntrkphot2);
   fChain->SetBranchAddress("pid_ptisophot1", &pid_ptisophot1, &b_pid_ptisophot1);
   fChain->SetBranchAddress("pid_ptisophot2", &pid_ptisophot2, &b_pid_ptisophot2);
   fChain->SetBranchAddress("pid_ntrkcsphot1", &pid_ntrkcsphot1, &b_pid_ntrkcsphot1);
   fChain->SetBranchAddress("pid_ntrkcsphot2", &pid_ntrkcsphot2, &b_pid_ntrkcsphot2);
   fChain->SetBranchAddress("pid_ptisocsphot1", &pid_ptisocsphot1, &b_pid_ptisocsphot1);
   fChain->SetBranchAddress("pid_ptisocsphot2", &pid_ptisocsphot2, &b_pid_ptisocsphot2);
   fChain->SetBranchAddress("pid_ecalisophot1", &pid_ecalisophot1, &b_pid_ecalisophot1);
   fChain->SetBranchAddress("pid_ecalisophot2", &pid_ecalisophot2, &b_pid_ecalisophot2);
   fChain->SetBranchAddress("pid_hcalisophot1", &pid_hcalisophot1, &b_pid_hcalisophot1);
   fChain->SetBranchAddress("pid_hcalisophot2", &pid_hcalisophot2, &b_pid_hcalisophot2);
   fChain->SetBranchAddress("ptjet1", &ptjet1, &b_ptjet1);
   fChain->SetBranchAddress("ptjet2", &ptjet2, &b_ptjet2);
   fChain->SetBranchAddress("ptjet3", &ptjet3, &b_ptjet3);
   fChain->SetBranchAddress("ptjet4", &ptjet4, &b_ptjet4);
   fChain->SetBranchAddress("ptcorrjet1", &ptcorrjet1, &b_ptcorrjet1);
   fChain->SetBranchAddress("ptcorrjet2", &ptcorrjet2, &b_ptcorrjet2);
   fChain->SetBranchAddress("ptcorrjet3", &ptcorrjet3, &b_ptcorrjet3);
   fChain->SetBranchAddress("ptcorrjet4", &ptcorrjet4, &b_ptcorrjet4);
   fChain->SetBranchAddress("etajet1", &etajet1, &b_etajet1);
   fChain->SetBranchAddress("etajet2", &etajet2, &b_etajet2);
   fChain->SetBranchAddress("etajet3", &etajet3, &b_etajet3);
   fChain->SetBranchAddress("etajet4", &etajet4, &b_etajet4);
   fChain->SetBranchAddress("phijet1", &phijet1, &b_phijet1);
   fChain->SetBranchAddress("phijet2", &phijet2, &b_phijet2);
   fChain->SetBranchAddress("phijet3", &phijet3, &b_phijet3);
   fChain->SetBranchAddress("phijet4", &phijet4, &b_phijet4);
   fChain->SetBranchAddress("betajet1", &betajet1, &b_betajet1);
   fChain->SetBranchAddress("betajet2", &betajet2, &b_betajet2);
   fChain->SetBranchAddress("betastarjet1", &betastarjet1, &b_betastarjet1);
   fChain->SetBranchAddress("betastarjet2", &betastarjet2, &b_betastarjet2);
   fChain->SetBranchAddress("btagvtxjet1", &btagvtxjet1, &b_btagvtxjet1);
   fChain->SetBranchAddress("btagtrkjet1", &btagtrkjet1, &b_btagtrkjet1);
   fChain->SetBranchAddress("btagvtxjet2", &btagvtxjet2, &b_btagvtxjet2);
   fChain->SetBranchAddress("btagtrkjet2", &btagtrkjet2, &b_btagtrkjet2);
   fChain->SetBranchAddress("ptDjet1", &ptDjet1, &b_ptDjet1);
   fChain->SetBranchAddress("rmsjet1", &rmsjet1, &b_rmsjet1);
   fChain->SetBranchAddress("ntrkjet1", &ntrkjet1, &b_ntrkjet1);
   fChain->SetBranchAddress("nneutjet1", &nneutjet1, &b_nneutjet1);
   fChain->SetBranchAddress("ptDjet2", &ptDjet2, &b_ptDjet2);
   fChain->SetBranchAddress("rmsjet2", &rmsjet2, &b_rmsjet2);
   fChain->SetBranchAddress("ntrkjet2", &ntrkjet2, &b_ntrkjet2);
   fChain->SetBranchAddress("nneutjet2", &nneutjet2, &b_nneutjet2);
   fChain->SetBranchAddress("assjet1", &assjet1, &b_assjet1);
   fChain->SetBranchAddress("assjet2", &assjet2, &b_assjet2);
   fChain->SetBranchAddress("deltaeta", &deltaeta, &b_deltaeta);
   fChain->SetBranchAddress("zeppenjet", &zeppenjet, &b_zeppenjet);
   fChain->SetBranchAddress("deltaphi", &deltaphi, &b_deltaphi);
   fChain->SetBranchAddress("deltaphinewvtx", &deltaphinewvtx, &b_deltaphinewvtx);
   fChain->SetBranchAddress("deltaphigg", &deltaphigg, &b_deltaphigg);
   fChain->SetBranchAddress("invmassjet", &invmassjet, &b_invmassjet);
   fChain->SetBranchAddress("invmass2g1j", &invmass2g1j, &b_invmass2g1j);
   fChain->SetBranchAddress("invmass2g2j", &invmass2g2j, &b_invmass2g2j);
   fChain->SetBranchAddress("pt2g2j", &pt2g2j, &b_pt2g2j);
   fChain->SetBranchAddress("eta2j", &eta2j, &b_eta2j);
   fChain->SetBranchAddress("phi2j", &phi2j, &b_phi2j);
   fChain->SetBranchAddress("pt2j", &pt2j, &b_pt2j);
   fChain->SetBranchAddress("nvtx", &nvtx, &b_nvtx);
   fChain->SetBranchAddress("sMet", &sMet, &b_sMet);
   fChain->SetBranchAddress("eMet", &eMet, &b_eMet);
   fChain->SetBranchAddress("phiMet", &phiMet, &b_phiMet);
   fChain->SetBranchAddress("signifMet", &signifMet, &b_signifMet);
   fChain->SetBranchAddress("eSmearedMet", &eSmearedMet, &b_eSmearedMet);
   fChain->SetBranchAddress("phiSmearedMet", &phiSmearedMet, &b_phiSmearedMet);
   fChain->SetBranchAddress("eShiftedMet", &eShiftedMet, &b_eShiftedMet);
   fChain->SetBranchAddress("phiShiftedMet", &phiShiftedMet, &b_phiShiftedMet);
   fChain->SetBranchAddress("eShiftedScaledMet", &eShiftedScaledMet, &b_eShiftedScaledMet);
   fChain->SetBranchAddress("phiShiftedScaledMet", &phiShiftedScaledMet, &b_phiShiftedScaledMet);
   fChain->SetBranchAddress("eSmearedShiftedMet", &eSmearedShiftedMet, &b_eSmearedShiftedMet);
   fChain->SetBranchAddress("phiSmearedShiftedMet", &phiSmearedShiftedMet, &b_phiSmearedShiftedMet);
   fChain->SetBranchAddress("sCorrMet", &sCorrMet, &b_sCorrMet);
   fChain->SetBranchAddress("eCorrMet", &eCorrMet, &b_eCorrMet);
   fChain->SetBranchAddress("phiCorrMet", &phiCorrMet, &b_phiCorrMet);
   fChain->SetBranchAddress("signifCorrMet", &signifCorrMet, &b_signifCorrMet);
   fChain->SetBranchAddress("smuCorrMet", &smuCorrMet, &b_smuCorrMet);
   fChain->SetBranchAddress("emuCorrMet", &emuCorrMet, &b_emuCorrMet);
   fChain->SetBranchAddress("phimuCorrMet", &phimuCorrMet, &b_phimuCorrMet);
   fChain->SetBranchAddress("signifmuCorrMet", &signifmuCorrMet, &b_signifmuCorrMet);
   fChain->SetBranchAddress("sNoHFMet", &sNoHFMet, &b_sNoHFMet);
   fChain->SetBranchAddress("eNoHFMet", &eNoHFMet, &b_eNoHFMet);
   fChain->SetBranchAddress("phiNoHFMet", &phiNoHFMet, &b_phiNoHFMet);
   fChain->SetBranchAddress("signifNoHFMet", &signifNoHFMet, &b_signifNoHFMet);
   fChain->SetBranchAddress("stcMet", &stcMet, &b_stcMet);
   fChain->SetBranchAddress("etcMet", &etcMet, &b_etcMet);
   fChain->SetBranchAddress("phitcMet", &phitcMet, &b_phitcMet);
   fChain->SetBranchAddress("signiftcMet", &signiftcMet, &b_signiftcMet);
   fChain->SetBranchAddress("sglobalPfMet", &sglobalPfMet, &b_sglobalPfMet);
   fChain->SetBranchAddress("eglobalPfMet", &eglobalPfMet, &b_eglobalPfMet);
   fChain->SetBranchAddress("phiglobalPfMet", &phiglobalPfMet, &b_phiglobalPfMet);
   fChain->SetBranchAddress("signifglobalPfMet", &signifglobalPfMet, &b_signifglobalPfMet);
   fChain->SetBranchAddress("scentralPfMet", &scentralPfMet, &b_scentralPfMet);
   fChain->SetBranchAddress("ecentralPfMet", &ecentralPfMet, &b_ecentralPfMet);
   fChain->SetBranchAddress("phicentralPfMet", &phicentralPfMet, &b_phicentralPfMet);
   fChain->SetBranchAddress("signifcentralPfMet", &signifcentralPfMet, &b_signifcentralPfMet);
   fChain->SetBranchAddress("eassocPfMet", &eassocPfMet, &b_eassocPfMet);
   fChain->SetBranchAddress("phiassocPfMet", &phiassocPfMet, &b_phiassocPfMet);
   fChain->SetBranchAddress("signifassocPfMet", &signifassocPfMet, &b_signifassocPfMet);
   fChain->SetBranchAddress("eassocOtherVtxPfMet", &eassocOtherVtxPfMet, &b_eassocOtherVtxPfMet);
   fChain->SetBranchAddress("phiassocOtherVtxPfMet", &phiassocOtherVtxPfMet, &b_phiassocOtherVtxPfMet);
   fChain->SetBranchAddress("signifassocOtherVtxPfMet", &signifassocOtherVtxPfMet, &b_signifassocOtherVtxPfMet);
   fChain->SetBranchAddress("etrkPfMet", &etrkPfMet, &b_etrkPfMet);
   fChain->SetBranchAddress("phitrkPfMet", &phitrkPfMet, &b_phitrkPfMet);
   fChain->SetBranchAddress("signiftrkPfMet", &signiftrkPfMet, &b_signiftrkPfMet);
   fChain->SetBranchAddress("ecleanPfMet", &ecleanPfMet, &b_ecleanPfMet);
   fChain->SetBranchAddress("phicleanPfMet", &phicleanPfMet, &b_phicleanPfMet);
   fChain->SetBranchAddress("signifcleanPfMet", &signifcleanPfMet, &b_signifcleanPfMet);
   fChain->SetBranchAddress("ecleanedSaclayPfMet", &ecleanedSaclayPfMet, &b_ecleanedSaclayPfMet);
   fChain->SetBranchAddress("phicleanedSaclayPfMet", &phicleanedSaclayPfMet, &b_phicleanedSaclayPfMet);
   fChain->SetBranchAddress("signifcleanedSaclayPfMet", &signifcleanedSaclayPfMet, &b_signifcleanedSaclayPfMet);
   fChain->SetBranchAddress("eminTypeICleanSaclayPfMet", &eminTypeICleanSaclayPfMet, &b_eminTypeICleanSaclayPfMet);
   fChain->SetBranchAddress("phiminTypeICleanSaclayPfMet", &phiminTypeICleanSaclayPfMet, &b_phiminTypeICleanSaclayPfMet);
   fChain->SetBranchAddress("signifminTypeICleanSaclayPfMet", &signifminTypeICleanSaclayPfMet, &b_signifminTypeICleanSaclayPfMet);
   fChain->SetBranchAddress("globalPfSums", &globalPfSums, &b_globalPfSums);
   fChain->SetBranchAddress("spfMet", &spfMet, &b_spfMet);
   fChain->SetBranchAddress("epfMet", &epfMet, &b_epfMet);
   fChain->SetBranchAddress("phipfMet", &phipfMet, &b_phipfMet);
   fChain->SetBranchAddress("signifpfMet", &signifpfMet, &b_signifpfMet);
   fChain->SetBranchAddress("spfMetType1", &spfMetType1, &b_spfMetType1);
   fChain->SetBranchAddress("epfMetType1", &epfMetType1, &b_epfMetType1);
   fChain->SetBranchAddress("phipfMetType1", &phipfMetType1, &b_phipfMetType1);
   fChain->SetBranchAddress("signifpfMetType1", &signifpfMetType1, &b_signifpfMetType1);
   fChain->SetBranchAddress("sMetGen", &sMetGen, &b_sMetGen);
   fChain->SetBranchAddress("eMetGen", &eMetGen, &b_eMetGen);
   fChain->SetBranchAddress("phiMetGen", &phiMetGen, &b_phiMetGen);
   fChain->SetBranchAddress("signifMetGen", &signifMetGen, &b_signifMetGen);
   fChain->SetBranchAddress("sMetGen2", &sMetGen2, &b_sMetGen2);
   fChain->SetBranchAddress("eMetGen2", &eMetGen2, &b_eMetGen2);
   fChain->SetBranchAddress("phiMetGen2", &phiMetGen2, &b_phiMetGen2);
   fChain->SetBranchAddress("npu", &npu, &b_npu);
   fChain->SetBranchAddress("NtotEvents", &NtotEvents, &b_NtotEvents);
   fChain->SetBranchAddress("xsection", &xsection, &b_xsection);
   fChain->SetBranchAddress("EquivLumi", &EquivLumi, &b_EquivLumi);
   fChain->SetBranchAddress("SampleID", &SampleID, &b_SampleID);
   fChain->SetBranchAddress("pu_weight", &pu_weight, &b_pu_weight);
   fChain->SetBranchAddress("pt_weight", &pt_weight, &b_pt_weight);
   fChain->SetBranchAddress("gen_custom_processId", &gen_custom_processId, &b_gen_custom_processId);
   fChain->SetBranchAddress("gen_pt_gamma1", &gen_pt_gamma1, &b_gen_pt_gamma1);
   fChain->SetBranchAddress("gen_pt_gamma2", &gen_pt_gamma2, &b_gen_pt_gamma2);
   fChain->SetBranchAddress("gen_eta_gamma1", &gen_eta_gamma1, &b_gen_eta_gamma1);
   fChain->SetBranchAddress("gen_eta_gamma2", &gen_eta_gamma2, &b_gen_eta_gamma2);
   fChain->SetBranchAddress("gen_phi_gamma1", &gen_phi_gamma1, &b_gen_phi_gamma1);
   fChain->SetBranchAddress("gen_phi_gamma2", &gen_phi_gamma2, &b_gen_phi_gamma2);
   fChain->SetBranchAddress("gen_pt_genjet1", &gen_pt_genjet1, &b_gen_pt_genjet1);
   fChain->SetBranchAddress("gen_pt_genjet2", &gen_pt_genjet2, &b_gen_pt_genjet2);
   fChain->SetBranchAddress("gen_eta_genjet1", &gen_eta_genjet1, &b_gen_eta_genjet1);
   fChain->SetBranchAddress("gen_eta_genjet2", &gen_eta_genjet2, &b_gen_eta_genjet2);
   fChain->SetBranchAddress("gen_phi_genjet1", &gen_phi_genjet1, &b_gen_phi_genjet1);
   fChain->SetBranchAddress("gen_phi_genjet2", &gen_phi_genjet2, &b_gen_phi_genjet2);
   fChain->SetBranchAddress("gen_mass_diphoton", &gen_mass_diphoton, &b_gen_mass_diphoton);
   fChain->SetBranchAddress("gen_pt_diphoton", &gen_pt_diphoton, &b_gen_pt_diphoton);
   fChain->SetBranchAddress("gen_eta_diphoton", &gen_eta_diphoton, &b_gen_eta_diphoton);
   fChain->SetBranchAddress("gen_phi_diphoton", &gen_phi_diphoton, &b_gen_phi_diphoton);
   fChain->SetBranchAddress("gen_mass_dijet", &gen_mass_dijet, &b_gen_mass_dijet);
   fChain->SetBranchAddress("gen_pt_dijet", &gen_pt_dijet, &b_gen_pt_dijet);
   fChain->SetBranchAddress("gen_eta_dijet", &gen_eta_dijet, &b_gen_eta_dijet);
   fChain->SetBranchAddress("gen_phi_dijet", &gen_phi_dijet, &b_gen_phi_dijet);
   fChain->SetBranchAddress("gen_zeppenfeld", &gen_zeppenfeld, &b_gen_zeppenfeld);
   fChain->SetBranchAddress("gen_pt_lep1", &gen_pt_lep1, &b_gen_pt_lep1);
   fChain->SetBranchAddress("gen_pt_lep2", &gen_pt_lep2, &b_gen_pt_lep2);
   fChain->SetBranchAddress("gen_eta_lep1", &gen_eta_lep1, &b_gen_eta_lep1);
   fChain->SetBranchAddress("gen_eta_lep2", &gen_eta_lep2, &b_gen_eta_lep2);
   fChain->SetBranchAddress("gen_phi_lep1", &gen_phi_lep1, &b_gen_phi_lep1);
   fChain->SetBranchAddress("gen_phi_lep2", &gen_phi_lep2, &b_gen_phi_lep2);
   fChain->SetBranchAddress("gen_pid_lep1", &gen_pid_lep1, &b_gen_pid_lep1);
   fChain->SetBranchAddress("gen_pid_lep2", &gen_pid_lep2, &b_gen_pid_lep2);
   fChain->SetBranchAddress("ptele1", &ptele1, &b_ptele1);
   fChain->SetBranchAddress("ptele2", &ptele2, &b_ptele2);
   fChain->SetBranchAddress("etaele1", &etaele1, &b_etaele1);
   fChain->SetBranchAddress("etaele2", &etaele2, &b_etaele2);
   fChain->SetBranchAddress("phiele1", &phiele1, &b_phiele1);
   fChain->SetBranchAddress("phiele2", &phiele2, &b_phiele2);
   fChain->SetBranchAddress("eneele1", &eneele1, &b_eneele1);
   fChain->SetBranchAddress("eneele2", &eneele2, &b_eneele2);
   fChain->SetBranchAddress("sIeIeele1", &sIeIeele1, &b_sIeIeele1);
   fChain->SetBranchAddress("sIeIeele2", &sIeIeele2, &b_sIeIeele2);
   fChain->SetBranchAddress("dphiele1", &dphiele1, &b_dphiele1);
   fChain->SetBranchAddress("dphiele2", &dphiele2, &b_dphiele2);
   fChain->SetBranchAddress("detaele1", &detaele1, &b_detaele1);
   fChain->SetBranchAddress("detaele2", &detaele2, &b_detaele2);
   fChain->SetBranchAddress("mhitsele1", &mhitsele1, &b_mhitsele1);
   fChain->SetBranchAddress("mhitsele2", &mhitsele2, &b_mhitsele2);
   fChain->SetBranchAddress("dcotele1", &dcotele1, &b_dcotele1);
   fChain->SetBranchAddress("dcotele2", &dcotele2, &b_dcotele2);
   fChain->SetBranchAddress("distele1", &distele1, &b_distele1);
   fChain->SetBranchAddress("distele2", &distele2, &b_distele2);
   fChain->SetBranchAddress("d0ele1", &d0ele1, &b_d0ele1);
   fChain->SetBranchAddress("d0ele2", &d0ele2, &b_d0ele2);
   fChain->SetBranchAddress("dzele1", &dzele1, &b_dzele1);
   fChain->SetBranchAddress("dzele2", &dzele2, &b_dzele2);
   fChain->SetBranchAddress("isoele1", &isoele1, &b_isoele1);
   fChain->SetBranchAddress("isoele2", &isoele2, &b_isoele2);
   fChain->SetBranchAddress("fullisoele1", &fullisoele1, &b_fullisoele1);
   fChain->SetBranchAddress("fullisoele2", &fullisoele2, &b_fullisoele2);
   fChain->SetBranchAddress("invMassele1g1", &invMassele1g1, &b_invMassele1g1);
   fChain->SetBranchAddress("invMassele1g2", &invMassele1g2, &b_invMassele1g2);
   fChain->SetBranchAddress("invMassele2g1", &invMassele2g1, &b_invMassele2g1);
   fChain->SetBranchAddress("invMassele2g2", &invMassele2g2, &b_invMassele2g2);
   fChain->SetBranchAddress("ptmu1", &ptmu1, &b_ptmu1);
   fChain->SetBranchAddress("ptmu2", &ptmu2, &b_ptmu2);
   fChain->SetBranchAddress("etamu1", &etamu1, &b_etamu1);
   fChain->SetBranchAddress("etamu2", &etamu2, &b_etamu2);
   fChain->SetBranchAddress("phimu1", &phimu1, &b_phimu1);
   fChain->SetBranchAddress("phimu2", &phimu2, &b_phimu2);
   fChain->SetBranchAddress("enemu1", &enemu1, &b_enemu1);
   fChain->SetBranchAddress("enemu2", &enemu2, &b_enemu2);
   fChain->SetBranchAddress("pixhitsmu1", &pixhitsmu1, &b_pixhitsmu1);
   fChain->SetBranchAddress("pixhitsmu2", &pixhitsmu2, &b_pixhitsmu2);
   fChain->SetBranchAddress("trkhitsmu1", &trkhitsmu1, &b_trkhitsmu1);
   fChain->SetBranchAddress("trkhitsmu2", &trkhitsmu2, &b_trkhitsmu2);
   fChain->SetBranchAddress("hitsmu1", &hitsmu1, &b_hitsmu1);
   fChain->SetBranchAddress("hitsmu2", &hitsmu2, &b_hitsmu2);
   fChain->SetBranchAddress("chi2mu1", &chi2mu1, &b_chi2mu1);
   fChain->SetBranchAddress("chi2mu2", &chi2mu2, &b_chi2mu2);
   fChain->SetBranchAddress("matchmu1", &matchmu1, &b_matchmu1);
   fChain->SetBranchAddress("matchmu2", &matchmu2, &b_matchmu2);
   fChain->SetBranchAddress("d0mu1", &d0mu1, &b_d0mu1);
   fChain->SetBranchAddress("d0mu2", &d0mu2, &b_d0mu2);
   fChain->SetBranchAddress("dzmu1", &dzmu1, &b_dzmu1);
   fChain->SetBranchAddress("dzmu2", &dzmu2, &b_dzmu2);
   fChain->SetBranchAddress("isomu1", &isomu1, &b_isomu1);
   fChain->SetBranchAddress("isomu2", &isomu2, &b_isomu2);
   fChain->SetBranchAddress("nWeightsPDF1", &nWeightsPDF1, &b_nWeightsPDF1);
   fChain->SetBranchAddress("nWeightsPDF2", &nWeightsPDF2, &b_nWeightsPDF2);
   fChain->SetBranchAddress("nWeightsPDF3", &nWeightsPDF3, &b_nWeightsPDF3);
   fChain->SetBranchAddress("nWeightsPDF4", &nWeightsPDF4, &b_nWeightsPDF4);
   fChain->SetBranchAddress("nWeightsPDF5", &nWeightsPDF5, &b_nWeightsPDF5);
   fChain->SetBranchAddress("nWeightsPDF6", &nWeightsPDF6, &b_nWeightsPDF6);
   fChain->SetBranchAddress("nWeightsPDF7", &nWeightsPDF7, &b_nWeightsPDF7);
   fChain->SetBranchAddress("nWeightsPDF8", &nWeightsPDF8, &b_nWeightsPDF8);
   fChain->SetBranchAddress("nWeightsPDF9", &nWeightsPDF9, &b_nWeightsPDF9);
   fChain->SetBranchAddress("nWeightsPDF10", &nWeightsPDF10, &b_nWeightsPDF10);
   fChain->SetBranchAddress("PDFweight1", &PDFweight1, &b_PDFweight1);
   fChain->SetBranchAddress("PDFweight2", &PDFweight2, &b_PDFweight2);
   fChain->SetBranchAddress("PDFweight3", &PDFweight3, &b_PDFweight3);
   fChain->SetBranchAddress("PDFweight4", &PDFweight4, &b_PDFweight4);
   fChain->SetBranchAddress("PDFweight5", &PDFweight5, &b_PDFweight5);
   fChain->SetBranchAddress("PDFweight6", &PDFweight6, &b_PDFweight6);
   fChain->SetBranchAddress("PDFweight7", &PDFweight7, &b_PDFweight7);
   fChain->SetBranchAddress("PDFweight8", &PDFweight8, &b_PDFweight8);
   fChain->SetBranchAddress("PDFweight9", &PDFweight9, &b_PDFweight9);
   fChain->SetBranchAddress("PDFweight10", &PDFweight10, &b_PDFweight10);
   Notify();
}

Bool_t fillPlot::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

#endif // #ifdef fillPlot_cxx
