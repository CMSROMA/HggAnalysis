#include <TTree.h>
#include <TFile.h>
#include <TLegend.h>
#include <TH2.h>
#include <TStyle.h>
#include <TROOT.h>
#include <TCanvas.h>
#include <TMath.h>
#include <iostream>
#include <fstream>
#include <vector>
#include <TLorentzVector.h>
#include <fillPlot.h>

inline double delta_phi(double phi1, double phi2) {

  double dphi = TMath::Abs(phi1 - phi2);
  return (dphi <= TMath::Pi())? dphi : TMath::TwoPi() - dphi;

}

// USE:
//
// .x finalize(data int lumi 2010, data int lumi 2011, pt1 cut, pt2 cut, ptj1 cut, ptj2 cut, deltae cut, zep cut, mjj cut, eb cat, r9 cat, cic selection)
// 
// example:
//
// .x finalize.C(33,230,40,25,20,15,2.5,2.5,300,1,1,4,1) 

vector <double> finalize(double int_exp_2010, double int_exp_2011, double pt1=50, double pt2=30, double pthiggsmin = -100, double pthiggsmax = -100, double ptj1=20, double ptj2=15, double misset=70, double deltae=2.5, double zep=2.5, double mjj=300, double deltap = 2.6, double jetmet=0., double p1met=0., double p2met=0., double hmet=0., double phigg=300., int eb = 1, int r9 = 1, int cic = 4, bool thirdcat = 0, bool leptontag = 0, bool leptonveto = 0, string variableMC = "massgg", string variableData = "massgg", int nbin = 200, double min = 90, double max = 190, string axis = "m(#gamma#gamma)[GeV]"){

  gROOT->SetStyle("Plain");
   
  gStyle->SetPalette(1);
  gStyle->SetOptStat(1111111);  // Show overflow, underflow + SumOfWeights 
  gStyle->SetOptFit(111110); 
  gStyle->SetOptFile(1); 
  
  gStyle->SetMarkerStyle(20);
  gStyle->SetMarkerSize(.3);
  gStyle->SetMarkerColor(1);

  TCanvas* c0 = new TCanvas("c0"," ",200,10,500,500);
  c0->Clear();

  //input files
  string mcnames[28];
  mcnames[0] = "box";
  mcnames[1] = "diphoton";
  mcnames[2] = "gjet";
  mcnames[3] = "qcdpt>40";
  mcnames[4] = "qcd30<pt<40";
  mcnames[5] = "dy";
  mcnames[6]  = "WenuG";
  mcnames[7]  = "WmnuG";
  mcnames[8]  = "WtnuG";
  mcnames[9]  = "ZeeG";
  mcnames[10] = "ZmmG";
  mcnames[11] = "ZttG";
  mcnames[12] = "ggWminus";
  mcnames[13] = "ggWplus";
  mcnames[14] = "ggZ";
  mcnames[15] = "ggtt";
  mcnames[16] = "ttjets";
  mcnames[17] = "Wjets";
  mcnames[18] = "WW";
  mcnames[19] = "WZ";
  mcnames[20] = "ZZ";
  mcnames[21] = "higgsgluglu";
  mcnames[22] = "higgsVBF";
  mcnames[23] = "higgsWH, Wleptonic";
  mcnames[24] = "higgsWH, Whadronic";
  mcnames[25] = "higgsZH, Zleptonic";
  mcnames[26] = "higgsZH, Zhadronic";
  mcnames[27] = "higgsZH, Zneutrinos";

  TFile* mc_2010[28];
  TFile* mc_2011[28];

  TFile* mc_gluglu_2011[7];
  TFile* mc_vbf_2011[7];
  TFile* mc_wzh_2011[7];
  TFile* mc_tth_2011[7];

  int h_masses[7] = {100,105,110,115,120,130,140};

  // TString redntpDir= "/xrootdfs/u2/xrootd/delre/Higgs/reduced/";

  // TString redntpDir= "root://pccmsrm27.cern.ch:1094//u2/xrootd/delre/Higgs/reduced/";
  TString redntpDir= "/xrootdfs/cms/local/delre/Higgs/reduced/";
//  TString redntpDir= "root://pccmsrm23.cern.ch:1094//u2/xrootd/meridian/Higgs/reduced/";
  TString preselectionLevel;


//    if (cic>0)
  preselectionLevel="cicloose";
//   else
//    preselectionLevel="preselectionCS";

  TString preselectionLevelCS="preselectionCS";
  // total data sample
  TFile* data = TFile::Open(redntpDir+"/redntp.42xv6b_data."+preselectionLevel+".regrPho_eCorr_30Nov.v2_jetid/merged/redntp_Photon-Run2011-30Nov2011-v1-DiPhotonSkimOnFly.root");
  TFile* datacs = TFile::Open(redntpDir+"/redntp.42xv6b_data.preselectionCS.regrPho_eCorr_30Nov.v2_jetid/merged/redntp_Photon-Run2011-30Nov2011-v1-DiPhotonSkimOnFly.root");

  if(int_exp_2010>0){
    // box samples
    mc_2010[0] = TFile::Open(redntpDir+"/redntp.42xv6b."+preselectionLevel+".regrPho_eCorr_30Nov.v2_jetid/merged/redntp_DiPhotonBox_Pt-25To250_7TeV-pythia6-Fall11-PU_S6_START42_V14B-v1.root");
    // diphoton jets samples                                                                                                                            
    mc_2010[1] = TFile::Open(redntpDir+"/redntp.42xv6b."+preselectionLevel+".regrPho_eCorr_30Nov.v2_jetid/merged/redntp_DiPhotonJets_7TeV-madgraph-Fall11-PU_S6_START42_V14B-v1.root");
    // gjet samples                                                                                                                                     
    mc_2010[2] = TFile::Open(redntpDir+"/redntp.42xv6b."+preselectionLevel+".regrPho_eCorr_30Nov.v2_jetid/merged/redntp_GJet_Pt-20_doubleEMEnriched_TuneZ2_7TeV-pythia6-Fall11-PU_S6_START42_V14B-v1.root");
    // qcd pt>40 samples                                                                                                                                
    mc_2010[3] = TFile::Open(redntpDir+"/redntp.42xv6b."+preselectionLevel+".regrPho_eCorr_30Nov.v2_jetid/merged/redntp_QCD_Pt-40_doubleEMEnriched_TuneZ2_7TeV-pythia6-Fall11-PU_S6_START42_V14B-v1.root");
    // qcd 30<pt<40 samples                                                                                                                            
    mc_2010[4] = TFile::Open(redntpDir+"/redntp.42xv6b."+preselectionLevel+".regrPho_eCorr_30Nov.v2_jetid/merged/redntp_QCD_Pt-30to40_doubleEMEnriched_TuneZ2_7TeV-pythia6-Fall11-PU_S6_START42_V14B-v1.root");
    // drell yan samples                                                                                                                                
    mc_2010[5] = TFile::Open(redntpDir+"/redntp.42xv6b_list2."+preselectionLevel+".regrPho_eCorr_30Nov.v2_jetid/merged/redntp_DYJetsToLL_TuneZ2_M-50_7TeV-madgraph-tauola-Fall11-PU_S6_START42_V14B-v1.root");
    // W+gamma, W->enu                                                                                                                                  
    mc_2010[6] = TFile::Open(redntpDir+"/redntp.42xv6b."+preselectionLevel+".regrPho_eCorr_30Nov.v2_jetid/merged/redntp_WGToENuG_TuneZ2_7TeV-madgraph-Fall11-PU_S6_START42_V14B-v1.root");
    // W+gamma, W->munu                                                                                                                                 
    mc_2010[7] = TFile::Open(redntpDir+"/redntp.42xv6b."+preselectionLevel+".regrPho_eCorr_30Nov.v2_jetid/merged/redntp_WGToMuNuG_TuneZ2_7TeV-madgraph-Fall11-PU_S6_START42_V14B-v1.root");
    // W+gamma, W->taunu                                                                                                                                
    mc_2010[8] = TFile::Open(redntpDir+"/redntp.42xv6b."+preselectionLevel+".regrPho_eCorr_30Nov.v2_jetid/merged/redntp_WGToTauNuG_TuneZ2_7TeV-madgraph-tauola-Fall11-PU_S6_START42_V14B-v1.root");
    // Z+gamma, Z->ee                                                                                                                                   
    mc_2010[9]  = TFile::Open(redntpDir+"/redntp.42xv6b."+preselectionLevel+".regrPho_eCorr_30Nov.v2_jetid/merged/redntp_ZGToEEG_TuneZ2_7TeV-madgraph-Fall11-PU_S6_START42_V14B-v1.root");
    // Z+gamma, Z->mumu                                                                                                                                 
    mc_2010[10] = TFile::Open(redntpDir+"/redntp.42xv6b."+preselectionLevel+".regrPho_eCorr_30Nov.v2_jetid/merged/redntp_ZGToMuMuG_TuneZ2_7TeV-madgraph-Fall11-PU_S6_START42_V14B-v1.root");
    // Z+gamma, Z->tautau                                                                                                                               
    mc_2010[11] = TFile::Open(redntpDir+"/redntp.42xv6b."+preselectionLevel+".regrPho_eCorr_30Nov.v2_jetid/merged/redntp_ZGToTauTauG_TuneZ2_7TeV-madgraph-tauola-Fall11-PU_S6_START42_V14B-v1.root");
    // W-gg                                                                                                                                             
    mc_2010[12] = TFile::Open(redntpDir+"/redntp.42xv6b."+preselectionLevel+".regrPho_eCorr_30Nov.v2_jetid/merged/redntp_Wminusgg_madgraph_Fall11_private_v2.root");
    // W+ gg                                                                                                                                            
    mc_2010[13] = TFile::Open(redntpDir+"/redntp.42xv6b."+preselectionLevel+".regrPho_eCorr_30Nov.v2_jetid/merged/redntp_Wplusgg_madgraph_Fall11_private_v2.root");
    // Zgg                                                                                                                                              
    mc_2010[14] = TFile::Open(redntpDir+"/redntp.42xv6b."+preselectionLevel+".regrPho_eCorr_30Nov.v2_jetid/merged/redntp_Zgg_madgraph_Fall11_private_v2.root");
    // ttgg                                                                                                                                             
    mc_2010[15] = TFile::Open(redntpDir+"/redntp.42xv6b."+preselectionLevel+".regrPho_eCorr_30Nov.v2_jetid/merged/redntp_ttgg_madgraph_Fall11_private_v2.root");
    // TTjets                                                                                                                                           
    mc_2010[16] = TFile::Open(redntpDir+"/redntp.42xv6b."+preselectionLevel+".regrPho_eCorr_30Nov.v2_jetid/merged/redntp_TTJets_TuneZ2_7TeV-madgraph-tauola-Fall11-PU_S6_START42_V14B-v1.root");
    // Wjets                                                                                                                                            
    mc_2010[17] = TFile::Open(redntpDir+"/redntp.42xv6b_list2."+preselectionLevel+".regrPho_eCorr_30Nov.v2_jetid/merged/redntp_WJetsToLNu_TuneZ2_7TeV-madgraph-tauola-Fall11-PU_S6_START42_V14B-v1.root");
    // WW                                                                                                                                               
    mc_2010[18] = TFile::Open(redntpDir+"/redntp.42xv6b."+preselectionLevel+".regrPho_eCorr_30Nov.v2_jetid/merged/redntp_WWTo2L2Nu_TuneZ2_7TeV_pythia6_tauola-Fall11-PU_S6_START42_V14B-v1.root");
    // WZ                                                                                                                                               
    mc_2010[19] = TFile::Open(redntpDir+"/redntp.42xv6b."+preselectionLevel+".regrPho_eCorr_30Nov.v2_jetid/merged/redntp_WZTo3LNu_TuneZ2_7TeV_pythia6_tauola-Fall11-PU_S6_START42_V14B-v1.root");
    // ZZ                                                                                                                                               
    mc_2010[20] = TFile::Open(redntpDir+"/redntp.42xv6b."+preselectionLevel+".regrPho_eCorr_30Nov.v2_jetid/merged/redntp_ZZTo2L2Nu_TuneZ2_7TeV_pythia6_tauola-Fall11-PU_S6_START42_V14B-v1.root");
    // gluglu higgs samples                                                                                                                             
    mc_2010[21] = TFile::Open(redntpDir+"/redntp.42xv6b."+preselectionLevel+".regrPho_eCorr_30Nov.v2_jetid/merged/redntp_GluGluToHToGG_M-120_7TeV-powheg-pythia6-Fall11-PU_S6_START42_V14B-v1.root");
    // vbf higgs samples                                                                                                                                
    mc_2010[22] = TFile::Open(redntpDir+"/redntp.42xv6b."+preselectionLevel+".regrPho_eCorr_30Nov.v2_jetid/merged/redntp_VBF_HToGG_M-120_7TeV-powheg-pythia6-Fall11-PU_S6_START42_V14B-v1.root");
    // W/Z/TT H higgs samples: only W->lnu                                                                                                              
    mc_2010[23] = TFile::Open(redntpDir+"/redntp.42xv6b."+preselectionLevel+".regrPho_eCorr_30Nov.v2_jetid/merged/redntp_WH_ZH_HToGG_M-120_7TeV-pythia6-Fall11-PU_S6_START42_V14B-v1.root");
    // W/Z/TT H higgs samples: only W->qq                                                                                                               
    mc_2010[24] = TFile::Open(redntpDir+"/redntp.42xv6b."+preselectionLevel+".regrPho_eCorr_30Nov.v2_jetid/merged/redntp_WH_ZH_HToGG_M-120_7TeV-pythia6-Fall11-PU_S6_START42_V14B-v1.root");
    // W/Z/TT H higgs samples: only Z->ll                                                                                                               
    mc_2010[25] = TFile::Open(redntpDir+"/redntp.42xv6b."+preselectionLevel+".regrPho_eCorr_30Nov.v2_jetid/merged/redntp_WH_ZH_HToGG_M-120_7TeV-pythia6-Fall11-PU_S6_START42_V14B-v1.root");
    // W/Z/TT H higgs samples: only Z->qq                                                                                                               
    mc_2010[26] = TFile::Open(redntpDir+"/redntp.42xv6b."+preselectionLevel+".regrPho_eCorr_30Nov.v2_jetid/merged/redntp_WH_ZH_HToGG_M-120_7TeV-pythia6-Fall11-PU_S6_START42_V14B-v1.root");
    // W/Z/TT H higgs samples: only Z->nunu                                                                                                             
    mc_2010[27] = TFile::Open(redntpDir+"/redntp.42xv6b."+preselectionLevel+".regrPho_eCorr_30Nov.v2_jetid/merged/redntp_WH_ZH_HToGG_M-120_7TeV-pythia6-Fall11-PU_S6_START42_V14B-v1.root");
  }

  if(int_exp_2011>0){
    // box samples
    mc_2011[0] = TFile::Open(redntpDir+"/redntp.42xv6b."+preselectionLevel+".regrPho_eCorr_30Nov.v2_jetid/merged/redntp_DiPhotonBox_Pt-25To250_7TeV-pythia6-Fall11-PU_S6_START42_V14B-v1.root");
    // diphoton jets samples                                                                                                                            
    mc_2011[1] = TFile::Open(redntpDir+"/redntp.42xv6b."+preselectionLevel+".regrPho_eCorr_30Nov.v2_jetid/merged/redntp_DiPhotonJets_7TeV-madgraph-Fall11-PU_S6_START42_V14B-v1.root");
    // gjet samples                                                                                                                                     
    mc_2011[2] = TFile::Open(redntpDir+"/redntp.42xv6b."+preselectionLevel+".regrPho_eCorr_30Nov.v2_jetid/merged/redntp_GJet_Pt-20_doubleEMEnriched_TuneZ2_7TeV-pythia6-Fall11-PU_S6_START42_V14B-v1.root");
    // qcd pt>40 samples                                                                                                                                
    mc_2011[3] = TFile::Open(redntpDir+"/redntp.42xv6b."+preselectionLevel+".regrPho_eCorr_30Nov.v2_jetid/merged/redntp_QCD_Pt-40_doubleEMEnriched_TuneZ2_7TeV-pythia6-Fall11-PU_S6_START42_V14B-v1.root");
    // qcd 30<pt<40 samples                                                                                                                            
    mc_2011[4] = TFile::Open(redntpDir+"/redntp.42xv6b."+preselectionLevel+".regrPho_eCorr_30Nov.v2_jetid/merged/redntp_QCD_Pt-30to40_doubleEMEnriched_TuneZ2_7TeV-pythia6-Fall11-PU_S6_START42_V14B-v1.root");
    // drell yan samples                                                                                                                                
    mc_2011[5] = TFile::Open(redntpDir+"/redntp.42xv6b_list2."+preselectionLevel+".regrPho_eCorr_30Nov.v2_jetid/merged/redntp_DYJetsToLL_TuneZ2_M-50_7TeV-madgraph-tauola-Fall11-PU_S6_START42_V14B-v1.root");
    // W+gamma, W->enu                                                                                                                                  
    mc_2011[6] = TFile::Open(redntpDir+"/redntp.42xv6b."+preselectionLevel+".regrPho_eCorr_30Nov.v2_jetid/merged/redntp_WGToENuG_TuneZ2_7TeV-madgraph-Fall11-PU_S6_START42_V14B-v1.root");
    // W+gamma, W->munu                                                                                                                                 
    mc_2011[7] = TFile::Open(redntpDir+"/redntp.42xv6b."+preselectionLevel+".regrPho_eCorr_30Nov.v2_jetid/merged/redntp_WGToMuNuG_TuneZ2_7TeV-madgraph-Fall11-PU_S6_START42_V14B-v1.root");
    // W+gamma, W->taunu                                                                                                                                
    mc_2011[8] = TFile::Open(redntpDir+"/redntp.42xv6b."+preselectionLevel+".regrPho_eCorr_30Nov.v2_jetid/merged/redntp_WGToTauNuG_TuneZ2_7TeV-madgraph-tauola-Fall11-PU_S6_START42_V14B-v1.root");
    // Z+gamma, Z->ee                                                                                                                                   
    mc_2011[9]  = TFile::Open(redntpDir+"/redntp.42xv6b."+preselectionLevel+".regrPho_eCorr_30Nov.v2_jetid/merged/redntp_ZGToEEG_TuneZ2_7TeV-madgraph-Fall11-PU_S6_START42_V14B-v1.root");
    // Z+gamma, Z->mumu                                                                                                                                 
    mc_2011[10] = TFile::Open(redntpDir+"/redntp.42xv6b."+preselectionLevel+".regrPho_eCorr_30Nov.v2_jetid/merged/redntp_ZGToMuMuG_TuneZ2_7TeV-madgraph-Fall11-PU_S6_START42_V14B-v1.root");
    // Z+gamma, Z->tautau                                                                                                                               
    mc_2011[11] = TFile::Open(redntpDir+"/redntp.42xv6b."+preselectionLevel+".regrPho_eCorr_30Nov.v2_jetid/merged/redntp_ZGToTauTauG_TuneZ2_7TeV-madgraph-tauola-Fall11-PU_S6_START42_V14B-v1.root");
    // W-gg                                                                                                                                             
    mc_2011[12] = TFile::Open(redntpDir+"/redntp.42xv6b."+preselectionLevel+".regrPho_eCorr_30Nov.v2_jetid/merged/redntp_Wminusgg_madgraph_Fall11_private_v2.root");
    // W+ gg                                                                                                                                            
    mc_2011[13] = TFile::Open(redntpDir+"/redntp.42xv6b."+preselectionLevel+".regrPho_eCorr_30Nov.v2_jetid/merged/redntp_Wplusgg_madgraph_Fall11_private_v2.root");
    // Zgg                                                                                                                                              
    mc_2011[14] = TFile::Open(redntpDir+"/redntp.42xv6b."+preselectionLevel+".regrPho_eCorr_30Nov.v2_jetid/merged/redntp_Zgg_madgraph_Fall11_private_v2.root");
    // ttgg                                                                                                                                             
    mc_2011[15] = TFile::Open(redntpDir+"/redntp.42xv6b."+preselectionLevel+".regrPho_eCorr_30Nov.v2_jetid/merged/redntp_ttgg_madgraph_Fall11_private_v2.root");
    // TTjets                                                                                                                                           
    mc_2011[16] = TFile::Open(redntpDir+"/redntp.42xv6b."+preselectionLevel+".regrPho_eCorr_30Nov.v2_jetid/merged/redntp_TTJets_TuneZ2_7TeV-madgraph-tauola-Fall11-PU_S6_START42_V14B-v1.root");
    // Wjets                                                                                                                                            
    mc_2011[17] = TFile::Open(redntpDir+"/redntp.42xv6b_list2."+preselectionLevel+".regrPho_eCorr_30Nov.v2_jetid/merged/redntp_WJetsToLNu_TuneZ2_7TeV-madgraph-tauola-Fall11-PU_S6_START42_V14B-v1.root");
    // WW                                                                                                                                               
    mc_2011[18] = TFile::Open(redntpDir+"/redntp.42xv6b."+preselectionLevel+".regrPho_eCorr_30Nov.v2_jetid/merged/redntp_WWTo2L2Nu_TuneZ2_7TeV_pythia6_tauola-Fall11-PU_S6_START42_V14B-v1.root");
    // WZ                                                                                                                                               
    mc_2011[19] = TFile::Open(redntpDir+"/redntp.42xv6b."+preselectionLevel+".regrPho_eCorr_30Nov.v2_jetid/merged/redntp_WZTo3LNu_TuneZ2_7TeV_pythia6_tauola-Fall11-PU_S6_START42_V14B-v1.root");
    // ZZ                                                                                                                                               
    mc_2011[20] = TFile::Open(redntpDir+"/redntp.42xv6b."+preselectionLevel+".regrPho_eCorr_30Nov.v2_jetid/merged/redntp_ZZTo2L2Nu_TuneZ2_7TeV_pythia6_tauola-Fall11-PU_S6_START42_V14B-v1.root");
    // gluglu higgs samples                                                                                                                             
    mc_2011[21] = TFile::Open(redntpDir+"/redntp.42xv6b."+preselectionLevel+".regrPho_eCorr_30Nov.v2_jetid/merged/redntp_GluGluToHToGG_M-120_7TeV-powheg-pythia6-Fall11-PU_S6_START42_V14B-v1.root");
    // vbf higgs samples                                                                                                                                
    mc_2011[22] = TFile::Open(redntpDir+"/redntp.42xv6b."+preselectionLevel+".regrPho_eCorr_30Nov.v2_jetid/merged/redntp_VBF_HToGG_M-120_7TeV-powheg-pythia6-Fall11-PU_S6_START42_V14B-v1.root");
    // W/Z/TT H higgs samples: only W->lnu                                                                                                              
    mc_2011[23] = TFile::Open(redntpDir+"/redntp.42xv6b."+preselectionLevel+".regrPho_eCorr_30Nov.v2_jetid/merged/redntp_WH_ZH_HToGG_M-120_7TeV-pythia6-Fall11-PU_S6_START42_V14B-v1.root");
    // W/Z/TT H higgs samples: only W->qq                                                                                                               
    mc_2011[24] = TFile::Open(redntpDir+"/redntp.42xv6b."+preselectionLevel+".regrPho_eCorr_30Nov.v2_jetid/merged/redntp_WH_ZH_HToGG_M-120_7TeV-pythia6-Fall11-PU_S6_START42_V14B-v1.root");
    // W/Z/TT H higgs samples: only Z->ll                                                                                                               
    mc_2011[25] = TFile::Open(redntpDir+"/redntp.42xv6b."+preselectionLevel+".regrPho_eCorr_30Nov.v2_jetid/merged/redntp_WH_ZH_HToGG_M-120_7TeV-pythia6-Fall11-PU_S6_START42_V14B-v1.root");
    // W/Z/TT H higgs samples: only Z->qq                                                                                                               
    mc_2011[26] = TFile::Open(redntpDir+"/redntp.42xv6b."+preselectionLevel+".regrPho_eCorr_30Nov.v2_jetid/merged/redntp_WH_ZH_HToGG_M-120_7TeV-pythia6-Fall11-PU_S6_START42_V14B-v1.root");
    // W/Z/TT H higgs samples: only Z->nunu                                                                                                             
    mc_2011[27] = TFile::Open(redntpDir+"/redntp.42xv6b."+preselectionLevel+".regrPho_eCorr_30Nov.v2_jetid/merged/redntp_WH_ZH_HToGG_M-120_7TeV-pythia6-Fall11-PU_S6_START42_V14B-v1.root");
  }
 
  for (int i=0;i<7;i++){
    char hmass[100]; sprintf(hmass,"%d",h_masses[i]);
    // gluglu samples
    mc_gluglu_2011[i] = TFile::Open(redntpDir+"/redntp.42xv6b."+preselectionLevel+".regrPho_eCorr_30Nov.v2_jetid/merged/redntp_GluGluToHToGG_M-"+ hmass+ "_7TeV-powheg-pythia6-Fall11-PU_S6_START42_V14B-v1.root");
    // vbf higgs samples  
    mc_vbf_2011[i] = TFile::Open(redntpDir+"/redntp.42xv6b."+preselectionLevel+".regrPho_eCorr_30Nov.v2_jetid/merged/redntp_VBF_HToGG_M-"+ hmass +"_7TeV-powheg-pythia6-Fall11-PU_S6_START42_V14B-v1.root");
    // W/Z H higgs samples  - all together 
    mc_wzh_2011[i] = TFile::Open(redntpDir+"/redntp.42xv6b."+preselectionLevel+".regrPho_eCorr_30Nov.v2_jetid/merged/redntp_WH_ZH_HToGG_M-"+ hmass +"_7TeV-pythia6-Fall11-PU_S6_START42_V14B-v1.root");
  }

 // k factors  
  double kfactordiphot = 1.3;
  double kfactordiphotmadgraph = 1.15;
  double kfactorgamjet = 1.3;
  double kfactorqcd = 1;
  double kfactordy = 1.15;
  double kfactorwg = 1;       // not 1, but I put the k-factor in the x-sec...
  double kfactorzg = 1;       // not 1, but I put the k-factor in the x-sec...
  double kfactorwmgg = 1;     // not 1, but I put the k-factor in the x-sec...
  double kfactorwpgg = 1;     // not 1, but I put the k-factor in the x-sec...
  double kfactorzgg = 1;      // not 1, but I put the k-factor in the x-sec...
  double kfactorttgg = 1;     // not 1, but I put the k-factor in the x-sec...
  double kfactorttjets = 1;   // not 1, but I put the k-factor in the x-sec...
  double kfactorwjets = 1;    // not 1, but I put the k-factor in the x-sec...
  double kfactorww = 1;       // not 1, but I put the k-factor in the x-sec...
  double kfactorwz = 1;       // not 1, but I put the k-factor in the x-sec...
  double kfactorzz = 1;       // not 1, but I put the k-factor in the x-sec...
  
 // cross sections and scaling
  double boosthiggs(1);
  double cross_mc[28];
  cross_mc[0] = 12.37 * kfactordiphot; // box
  cross_mc[1] = 134 * kfactordiphotmadgraph; // diphoton jets
  cross_mc[2] = 493.44 * kfactorgamjet; // gjet
  cross_mc[3] = 40392 * kfactorqcd; // qcd pt>40
  cross_mc[4] = 9610 * kfactorqcd; // qcd 30<pt<40 
  cross_mc[5] = 2321 * kfactordy; // drell yan
  cross_mc[6]  = 137.5 * kfactorwg;            // Wgamma, W->enu
  cross_mc[7]  = 137.5 * kfactorwg;            // Wgamma, W->munu 
  cross_mc[8]  = 137.5 * kfactorwg;            // Wgamma, W->taunu 
  cross_mc[9]  = 45.2  * kfactorzg;            // Zgamma, Z->ee 
  cross_mc[10] = 45.2  * kfactorzg;            // Zgamma, Z->mm
  cross_mc[11] = 45.2  * kfactorzg;            // Zgamma, Z->tt 
  // still to be double checked - init
  cross_mc[12] = 0.00581 * kfactorwmgg;        // private sample, taken from AN2012_035_v4 
  cross_mc[13] = 0.00888 * kfactorwpgg;        // private sample, taken from AN2012_035_v4
  cross_mc[14] = 0.00169 * kfactorzgg;         // private sample, taken from AN2012_035_v4
  cross_mc[15] = 0.00998 * kfactorttgg;        // private sample, taken from AN2012_035_v4  
  cross_mc[16] = 157.5 * kfactorttjets;        // taken from AN2012_035_v4
  cross_mc[17] = 31314. * kfactorwjets;        // taken from AN2012_035_v4
  cross_mc[18] = 4.51 * kfactorww;             // taken from AN2012_035_v4
  cross_mc[19] = 0.595 * kfactorwz;            // taken from AN2012_035_v4
  cross_mc[20] = 0.119 * kfactorzz;            // taken from AN2012_035_v4
  // still to be double checked - end
  /* standard model Higgs cross sections at 120 GeV
  cross_mc[21]  = 16.63  * 2.13e-03 * boosthiggs; // glu glu higgs SM 120 
  cross_mc[22]  = 1.269  * 2.13e-03 * boosthiggs; // vbf higgs SM 120
  cross_mc[23]  = 0.6561 * 2.13e-03 * boosthiggs * 0.3257;  // WH, W->lnu, higgs SM 120 
  cross_mc[24]  = 0.6561 * 2.13e-03 * boosthiggs * 0.6760;  // WH, W->qq,  higgs SM 120 
  cross_mc[25] = 0.3598 * 2.13e-03 * boosthiggs * 0.10096; // ZH, Z->ll,  higgs SM 120 
  cross_mc[26] = 0.3598 * 2.13e-03 * boosthiggs * 0.6991;  // ZH, Z->qq,  higgs SM 120 
  cross_mc[27] = 0.3598 * 2.13e-03 * boosthiggs * 0.20;    // ZH, Z->nn,  higgs SM 120   
  standard model Higgs cross sections at 120 GeV  */
  // fermiophobic Higgs cross sections at 120 GeV
  cross_mc[21]  = 0; // glu glu FP higgs 
  cross_mc[22]  = 1.269  * 0.0231 * boosthiggs; // vbf FP higgs
  cross_mc[23]  = 0.6561 * 0.0231 * boosthiggs * 0.3257;  // WH, W->lnu, higgs FP 120
  cross_mc[24]  = 0.6561 * 0.0231 * boosthiggs * 0.6760;  // WH, W->qq,  higgs FP 120 
  cross_mc[25]  = 0.3598 * 0.0231 * boosthiggs * 0.10096; // ZH, Z->ll,  higgs FP 120
  cross_mc[26]  = 0.3598 * 0.0231 * boosthiggs * 0.6991;  // ZH, Z->qq,  higgs FP 120
  cross_mc[27]  = 0.3598 * 0.0231 * boosthiggs * 0.20;    // ZH, Z->nn,  higgs FP 120 
 
  // getting the number of original events in each sample (processed with CMSSW)
  int n_mc_2010[28], n_mc_2011[28],n_gluglu_2011[7],n_vbf_2011[7],n_wzh_2011[7],n_tth_2011[7];
  for(int i=0; i<23; i++){
    n_mc_2010[i] = n_mc_2011[i] = 0;
    if(int_exp_2010>0) n_mc_2010[i] = ((TH1D*)mc_2010[i]->Get("ptphotgen1"))->GetEntries();
    if(int_exp_2011>0) n_mc_2011[i] = ((TH1D*)mc_2011[i]->Get("ptphotgen1"))->GetEntries();
  }
  if(int_exp_2010>0) {
    n_mc_2010[23]  = ((TH1D*)mc_2010[23]->Get("ptphotgen1wl"))->GetEntries();
    n_mc_2010[24]  = ((TH1D*)mc_2010[24]->Get("ptphotgen1wh"))->GetEntries();
    n_mc_2010[25]  = ((TH1D*)mc_2010[25]->Get("ptphotgen1zl"))->GetEntries();
    n_mc_2010[26]  = ((TH1D*)mc_2010[26]->Get("ptphotgen1zh"))->GetEntries();
    n_mc_2010[27]  = ((TH1D*)mc_2010[27]->Get("ptphotgen1zn"))->GetEntries();
  }
  if(int_exp_2011>0) {
    n_mc_2011[23]  = ((TH1D*)mc_2011[23]->Get("ptphotgen1wl"))->GetEntries();
    n_mc_2011[24]  = ((TH1D*)mc_2011[24]->Get("ptphotgen1wh"))->GetEntries();
    n_mc_2011[25]  = ((TH1D*)mc_2011[25]->Get("ptphotgen1zl"))->GetEntries();
    n_mc_2011[26]  = ((TH1D*)mc_2011[26]->Get("ptphotgen1zh"))->GetEntries();
    n_mc_2011[27]  = ((TH1D*)mc_2011[27]->Get("ptphotgen1zn"))->GetEntries();
  }

  for(int i=0; i<7; i++){
    n_gluglu_2011[i] = n_vbf_2011[i] = n_wzh_2011[i] = n_tth_2011[i] = 0;
    n_gluglu_2011[i] = ((TH1D*)mc_gluglu_2011[i]->Get("ptphotgen1"))->GetEntries();
    n_vbf_2011[i] = ((TH1D*)mc_vbf_2011[i]->Get("ptphotgen1"))->GetEntries();
    n_wzh_2011[i] = ((TH1D*)mc_wzh_2011[i]->Get("ptphotgen1"))->GetEntries();
//     n_tth_2011[i] = ((TH1D*)mc_tth_2011[i]->Get("ptphotgen1"))->GetEntries();
  }

  // setting the scaling factor to actual lumi
  double scale_mc_2010[28], scale_mc_2011[28];
  for(int i=0; i<28; i++){
    scale_mc_2010[i] = scale_mc_2011[i] = 0; 
    if(int_exp_2010>0) scale_mc_2010[i] = cross_mc[i] * int_exp_2010 / n_mc_2010[i];
    if(int_exp_2011>0) scale_mc_2011[i] = cross_mc[i] * int_exp_2011 / n_mc_2011[i];
  }

  // char for output name
  char name[1000];
  char allcut[3000];
  sprintf(allcut,"%3.1f%s%3.1f%s%3.1f%s%3.1f%s%3.1f%s%3.1f%s%3.1f%s%3.1f%s%3.1f%s%3.1f%s%3.1f%s%3.1f%s%3.1f%s%3.1f%s%3.1f%s%3.1f%s%3.1f%s%d%s%d%s%d%s%d%s%d",pt1,"_",pt2,"_",pthiggsmin,"_",pthiggsmax,"_",ptj1,"_",ptj2,"_",misset,"_",deltae,"_",zep,"_",mjj,"_",deltap,"_",jetmet,"_",p1met,"_",p2met,"_",hmet,"_",phigg,"_",eb,"_",r9,"_",thirdcat,"_",leptontag,"_",leptonveto,"_",cic);

  // output root file
  sprintf(name,"%s%s%s%s%s","results_gg/histo_",variableData.c_str(),"_",allcut,".root");
  TFile * hOutputFile   = new TFile(name, "RECREATE" ) ;

  // histograms needed by the machinery
  TH1D* vardata = new TH1D("vardata","vardata",nbin,min,max);
  TH1D* vardatacs = new TH1D("vardatacs","vardatacs",nbin,min,max);
  TH1D* var_mc_2010[28];
  TH1D* var_mc_2011[28];
  TH1D* var_gluglu_2011[7];
  TH1D* var_vbf_2011[7];
  TH1D* var_wzh_2011[7];
  TH1D* var_tth_2011[7];
  TH1D * var[14];
  for (int i=0; i<14; i++) {
    sprintf(name,"%s%d","var",i);
    var[i] = new TH1D(name,name,nbin,min,max);
  }
  for (int i=0; i<28; i++) {
    sprintf(name,"%s%d","var_mc_2010_",i);
    var_mc_2010[i] = new TH1D(name,name,nbin,min,max);
    sprintf(name,"%s%d","var_mc_2011_",i);
    var_mc_2011[i] = new TH1D(name,name,nbin,min,max);
  }
  for (int i=0; i<7; i++) {
    sprintf(name,"%s%d","var_gluglu_2011_",h_masses[i]);
    var_gluglu_2011[i] = new TH1D(name,name,nbin,min,max);
    sprintf(name,"%s%d","var_vbf_2011_",h_masses[i]);
    var_vbf_2011[i] = new TH1D(name,name,nbin,min,max);
    sprintf(name,"%s%d","var_wzh_2011_",h_masses[i]);
    var_wzh_2011[i] = new TH1D(name,name,nbin,min,max);
//     sprintf(name,"%s%d","var_tth_2011_",h_masses[i]);
//     var_tth_2011[i] = new TH1D(name,name,nbin,min,max);
  }


  // creating the fillers and setting cuts
  fillPlot data_fill((TTree*)data->Get("AnaTree"), 1);
  fillPlot datacs_fill((TTree*)datacs->Get("AnaTree"), 1);
  data_fill.Setcuts(pt1,pt2,pthiggsmin,pthiggsmax,ptj1,ptj2,misset,deltae,zep,mjj,deltap,jetmet,p1met,p2met,hmet,phigg,eb,r9,thirdcat,leptontag,leptonveto);
  datacs_fill.Setcuts(pt1,pt2,pthiggsmin,pthiggsmax,ptj1,ptj2,misset,deltae,zep,mjj,deltap,jetmet,p1met,p2met,hmet,phigg,eb,r9,thirdcat,leptontag,leptonveto);
  if(cic>0)
    {
      data_fill.setCic(cic);
      datacs_fill.setCic(cic);
    }

  fillPlot* mc_2010_fill[28];
  fillPlot* mc_2011_fill[28];
  fillPlot* mc_gluglu_2011_fill[7];  
  fillPlot* mc_vbf_2011_fill[7];  
  fillPlot* mc_wzh_2011_fill[7];  
  fillPlot* mc_tth_2011_fill[7];  

  for (int i=0; i<28; i++){
    if(int_exp_2010>0) mc_2010_fill[i] = new fillPlot((TTree*)mc_2010[i]->Get("AnaTree"), 1);
    if(int_exp_2011>0) mc_2011_fill[i] = new fillPlot((TTree*)mc_2011[i]->Get("AnaTree"), 1);
    if(int_exp_2010>0) mc_2010_fill[i]->Setcuts(pt1,pt2,pthiggsmin,pthiggsmax,ptj1,ptj2,misset,deltae,zep,mjj,deltap,jetmet,p1met,p2met,hmet,phigg,eb,r9,thirdcat,leptontag,leptonveto);
    if(int_exp_2011>0) mc_2011_fill[i]->Setcuts(pt1,pt2,pthiggsmin,pthiggsmax,ptj1,ptj2,misset,deltae,zep,mjj,deltap,jetmet,p1met,p2met,hmet,phigg,eb,r9,thirdcat,leptontag,leptonveto);
    
    mc_2011_fill[i]->DoPuReweight();
    mc_2011_fill[i]->DoPtReweight();
    
    if(cic>0){
      if(int_exp_2010>0) mc_2010_fill[i]->setCic(cic);
      if(int_exp_2011>0) mc_2011_fill[i]->setCic(cic);
    }
  }
 
  for (int i=0; i<7; i++){
    mc_gluglu_2011_fill[i] = new fillPlot((TTree*)mc_gluglu_2011[i]->Get("AnaTree"), 1);
    mc_vbf_2011_fill[i] = new fillPlot((TTree*)mc_vbf_2011[i]->Get("AnaTree"), 1);
    mc_wzh_2011_fill[i] = new fillPlot((TTree*)mc_wzh_2011[i]->Get("AnaTree"), 1);
//     mc_tth_2011_fill[i] = new fillPlot((TTree*)mc_tth_2011[i]->Get("AnaTree"), 1);
    mc_gluglu_2011_fill[i]->Setcuts(pt1,pt2,pthiggsmin,pthiggsmax,ptj1,ptj2,misset,deltae,zep,mjj,deltap,jetmet,p1met,p2met,hmet,phigg,eb,r9,thirdcat,leptontag,leptonveto);
    mc_vbf_2011_fill[i]->Setcuts(pt1,pt2,pthiggsmin,pthiggsmax,ptj1,ptj2,misset,deltae,zep,mjj,deltap,jetmet,p1met,p2met,hmet,phigg,eb,r9,thirdcat,leptontag,leptonveto);
    mc_wzh_2011_fill[i]->Setcuts(pt1,pt2,pthiggsmin,pthiggsmax,ptj1,ptj2,misset,deltae,zep,mjj,deltap,jetmet,p1met,p2met,hmet,phigg,eb,r9,thirdcat,leptontag,leptonveto);
//     mc_tth_2011_fill[i]->Setcuts(pt1,pt2,pthiggsmin,pthiggsmax,ptj1,ptj2,deltae,zep,mjj,deltap,eb,r9,thirdcat);
    mc_gluglu_2011_fill[i]->DoPuReweight();
    mc_vbf_2011_fill[i]->DoPuReweight();
    mc_wzh_2011_fill[i]->DoPuReweight();
//     mc_tth_2011_fill[i]->DoPuReweight();
    mc_gluglu_2011_fill[i]->DoPtReweight();
    mc_vbf_2011_fill[i]->DoPtReweight();
    mc_wzh_2011_fill[i]->DoPtReweight();
//     mc_tth_2011_fill[i]->DoPtReweight();
    mc_gluglu_2011_fill[i]->setCic(cic);
    mc_vbf_2011_fill[i]->setCic(cic);
    mc_wzh_2011_fill[i]->setCic(cic);
//     mc_tth_2011_fill[i]->setCic(cic);
 }
  // smear mc
//   for (int i=0; i<28; i++){
//     if(int_exp_2010>0) mc_2010_fill[i]->DoSmearing(1.,0.0001);
//     if(int_exp_2011>0) mc_2011_fill[i]->DoSmearing(1.,0.0001);
//   }

  // filling histograms
  std::cout << " ++++++++++++++ DATA ++++++++++++++++" << std::endl;
  cout << "running over " << ((TTree*)data->Get("AnaTree"))->GetEntries("") << " data events" <<  endl;
  if (variableData == "massgg") {  
    sprintf(name,"%s%s%s","results_gg/events_",allcut,".txt");
    data_fill.Writetxt(name);
    sprintf(name,"%s%s%s","results_gg/events_",allcut,".root");
    data_fill.WriteRoot(name);
  }
  vardata->Add(data_fill.Plot(variableData,"data", nbin, min, max,0,100)); 
  std::cout << "Selected events on data " << vardata->GetEntries() << std::endl;
  cout << "running over " << ((TTree*)datacs->Get("AnaTree"))->GetEntries("") << " data events (for cs)" <<  endl; 

  if (variableData == "massgg") {  
    sprintf(name,"%s%s%s","results_gg/events_",allcut,"_cs.txt");
    datacs_fill.Writetxt(name);
    sprintf(name,"%s%s%s","results_gg/events_",allcut,"_cs.root");
    datacs_fill.WriteRoot(name);
  }
  vardatacs->Add(datacs_fill.Plot(variableData,"datacs", nbin, min, max, 1,100)); 
  std::cout << "Selected events on data cs " << vardatacs->GetEntries() << std::endl;

  std::cout << " ++++++++++++++ MC ++++++++++++++++" << std::endl;
  for (int i=0; i<28; i++){ 
    sprintf(name,"%s%s",mcnames[i].c_str()," 2010");
    if(int_exp_2010>0) {
      cout << "running over " << ((TTree*)mc_2010[i]->Get("AnaTree"))->GetEntries("") << " " << name << " events" <<  endl; 
      sprintf(name,"%s%s%s%s%s","results_gg/events_",mcnames[i].c_str(),"_2010_",allcut,".root");
      //      if (variableMC == "massgg") mc_2010_fill[i]->WriteRoot(name);
      var_mc_2010[i]->Add(mc_2010_fill[i]->Plot(variableMC, name, nbin, min, max,0,i));
      std::cout << "Selected events on mc2010 " << name << " " << var_mc_2010[i]->GetEntries() << std::endl;
    }
    sprintf(name,"%s%s",mcnames[i].c_str()," 2011");
    if(int_exp_2011>0) {
      cout << "running over " << ((TTree*)mc_2011[i]->Get("AnaTree"))->GetEntries("") << " " << name << " events" <<  endl; 
      sprintf(name,"%s%s%s%s%s","results_gg/events_",mcnames[i].c_str(),"_2011_",allcut,".root");
      if (variableMC == "massgg") {
	//	mc_2011_fill[i]->WriteRoot(name);
// 	sprintf(name,"%s%s%s%s%s","results_gg/events_",mcnames[i].c_str(),"_2011_",allcut,".txt");
// 	if(i==7) mc_2011_fill[i]->Writetxt(name);
      }
      var_mc_2011[i]->Add(mc_2011_fill[i]->Plot(variableMC, name, nbin, min, max,0,i));
      std::cout << "Selected events on mc2011 " << name << " " << var_mc_2011[i]->GetEntries() << std::endl;
    }

  }

  std::cout << " ++++++++++++++ signal MC ++++++++++++++++" << std::endl;
  for (int i=0; i<7; i++){ 
    cout << "running over " << ((TTree*)mc_gluglu_2011[i]->Get("AnaTree"))->GetEntries("") << " gluglu M=" << h_masses[i] << " events" <<  endl; 
    sprintf(name,"%s%d%s%s%s","results_gg/events_gluglu",h_masses[i],"_2011_",allcut,".root");
    //    if (variableMC == "massgg") mc_gluglu_2011_fill[i]->WriteRoot(name);
    var_gluglu_2011[i]->Add(mc_gluglu_2011_fill[i]->Plot(variableMC, name, nbin, min, max,0,50));
    std::cout << "Selected events on mc2011 gluglu " << h_masses[i] << " " << var_gluglu_2011[i]->GetEntries() << std::endl;
 
    cout << "running over " << ((TTree*)mc_vbf_2011[i]->Get("AnaTree"))->GetEntries("") << " vbf M=" << h_masses[i] << " events" <<  endl; 
    sprintf(name,"%s%d%s%s%s","results_gg/events_vbf",h_masses[i],"_2011_",allcut,".root");
    //    if (variableMC == "massgg") mc_vbf_2011_fill[i]->WriteRoot(name);
    var_vbf_2011[i]->Add(mc_vbf_2011_fill[i]->Plot(variableMC, name, nbin, min, max, 0, 50));
    std::cout << "Selected events on mc2011 vbf " << h_masses[i] << " " << var_vbf_2011[i]->GetEntries() << std::endl;

    cout << "running over " << ((TTree*)mc_wzh_2011[i]->Get("AnaTree"))->GetEntries("") << " wzh M=" << h_masses[i] << " events" <<  endl; 
    sprintf(name,"%s%d%s%s%s","results_gg/events_wzh",h_masses[i],"_2011_",allcut,".root");
    //    if (variableMC == "massgg") mc_wzh_2011_fill[i]->WriteRoot(name);
    var_wzh_2011[i]->Add(mc_wzh_2011_fill[i]->Plot(variableMC, name, nbin, min, max, 0, 50));
    std::cout << "Selected events on mc2011 wzh " << h_masses[i] << " " << var_wzh_2011[i]->GetEntries() << std::endl;

//     cout << "running over " << ((TTree*)mc_tth_2011[i]->Get("AnaTree"))->GetEntries("") << " tth M=" << h_masses[i] << " events" <<  endl; 
//     sprintf(name,"%s%d%s%s%s","results_gg/events_tth",h_masses[i],"_2011_",allcut,".root");
//     if (variableMC == "massgg") mc_tth_2011_fill[i]->WriteRoot(name);
//     var_tth_2011[i]->Add(mc_tth_2011_fill[i]->Plot(variableMC, name, nbin, min, max, 0, 50));
//     std::cout << "Selected events on mc2011 tth " << h_masses[i] << " " << var_tth_2011[i]->GetEntries() << std::endl;
 }

  // scale mc to equivalent lumi
  for (int i=0; i<28; i++){ 
    if(int_exp_2010>0) var_mc_2010[i]->Scale(scale_mc_2010[i]);  
    if(int_exp_2011>0) var_mc_2011[i]->Scale(scale_mc_2011[i]);  
  }

  // counting number of events passing selection (scaled)
  double num_mc_2010[28],num_mc_2011[28],num_gluglu_2011[7],num_vbf_2011[7],num_wzh_2011[7],num_tth_2011[7]; 
  double num_uns_mc_2010[28],num_uns_mc_2011[28],num_uns_gluglu_2011[7],num_uns_vbf_2011[7],num_uns_wzh_2011[7],num_uns_tth_2011[7]; 

  for (int i=0; i<28; i++){ 
    num_mc_2010[i] = num_mc_2011[i] = 0;
    if(int_exp_2010>0) num_mc_2010[i] = var_mc_2010[i]->Integral();  
    if(int_exp_2011>0) num_mc_2011[i] = var_mc_2011[i]->Integral();  
    num_uns_mc_2010[i] = num_uns_mc_2011[i] = 0;
    if(int_exp_2010>0) num_uns_mc_2010[i] = var_mc_2010[i]->GetEntries();  
    if(int_exp_2011>0) num_uns_mc_2011[i] = var_mc_2011[i]->GetEntries();  
  }
  for(int i=0; i<7; i++){
    num_gluglu_2011[i] = num_vbf_2011[i] = num_wzh_2011[i] = num_tth_2011[i] = 0;
    num_gluglu_2011[i] = num_gluglu_2011[i] = var_gluglu_2011[i]->Integral(); 
    num_vbf_2011[i] = num_vbf_2011[i] = var_vbf_2011[i]->Integral(); 
    num_wzh_2011[i] = num_wzh_2011[i] = var_wzh_2011[i]->Integral(); 
    num_uns_gluglu_2011[i] = num_uns_vbf_2011[i] = num_uns_wzh_2011[i] = num_uns_tth_2011[i] = 0;
    num_uns_gluglu_2011[i] = num_uns_gluglu_2011[i] = var_gluglu_2011[i]->GetEntries(); 
    num_uns_vbf_2011[i] = num_uns_vbf_2011[i] = var_vbf_2011[i]->GetEntries(); 
    num_uns_wzh_2011[i] = num_uns_wzh_2011[i] = var_wzh_2011[i]->GetEntries(); 
  }

  // add two QCD bins
  // if(int_exp_2010>0) var_mc_2010[3]->Add(var_mc_2010[4]);
  // if(int_exp_2011>0) var_mc_2011[3]->Add(var_mc_2011[4]);
  
  // scale control sample
  vardata->Sumw2();
  //  vardatacs->Sumw2();
  double num_data =  vardata->Integral();
  double num_data_cs = vardatacs->Integral();  
  vardatacs->Scale(num_data/num_data_cs); 

  // stack histograms  
  for (int i=1; i<nbin+1; i++){      
    for (int j=0; j<14; j++){            
      int offset(0);
      if(j>0) offset = 6;    // to add higgs contributions up                                                                        
      if(j>1) offset = 8;    // to add dibosons contribuions up                                                                      
      if(j>6) offset = 9;    // to add Wgg contributions up                                                                          
      if(j>7) offset =11;    // to add Zg contributions up                                                                           
      if(j>8) offset =13;    // to add Wg contributions up                                                                           
      if(j>10) offset =14;   // to add QCD contributions up
      for (int k=0 ; k<28-j-offset; k++){ 
	if(int_exp_2010>0) var[j]->SetBinContent(i,var_mc_2010[k]->GetBinContent(i) + var[j]->GetBinContent(i));
	if(int_exp_2011>0) var[j]->SetBinContent(i,var_mc_2011[k]->GetBinContent(i) + var[j]->GetBinContent(i));
      }	
    }    
  }
 
  //final plots
  char ytitle[100];
  sprintf(ytitle,"%s%d%s","N_{ev}/",int(int_exp_2010+int_exp_2011),"pb^{-1}");
  for(int i=0; i<14; i++){
    var[i]->SetTitle("");
    var[i]->SetStats(0);
    var[i]->SetTitleOffset(1.25,"Y");
    var[i]->SetYTitle(ytitle);
    var[i]->SetXTitle(axis.c_str());
    var[i]->SetLineColor(kBlack);
    var[i]->SetLineWidth(2);
  }

  //legenda
  TLegendEntry *legge;
  TLegend *leg;
  leg = new TLegend(0.6,0.6,0.85,0.85);
  leg->SetFillStyle(0); leg->SetBorderSize(0); leg->SetTextSize(0.05);
  leg->SetFillColor(0);
  sprintf(name,"H M(120)");
  if(boosthiggs!=1) sprintf(name,"%s%2.0f", "H M(120) x ", boosthiggs);

  legge = leg->AddEntry(var[0], name, "f");
  legge = leg->AddEntry(var[1], "VV", "f");
  legge = leg->AddEntry(var[2], "Wjets", "f");
  legge = leg->AddEntry(var[3], "ttjets", "f");
  legge = leg->AddEntry(var[4], "ttgg", "f");
  legge = leg->AddEntry(var[5], "Zgg", "f");
  legge = leg->AddEntry(var[6], "Wgg", "f");
  legge = leg->AddEntry(var[7], "Zgamma", "f");
  legge = leg->AddEntry(var[8], "Wgamma", "f");
  legge = leg->AddEntry(var[9], "DY", "f");
  legge = leg->AddEntry(var[10], "QCD", "f");
  legge = leg->AddEntry(var[11], "#gamma + jets", "f");
  legge = leg->AddEntry(var[12], "di-#gamma + jets", "f");
  legge = leg->AddEntry(var[13], "di-#gamma box", "f");

  //mc only plot
  var[0]->SetFillColor(kYellow);
  var[0]->Draw();
  var[1]->SetFillColor(6);
  var[1]->Draw("same");
  var[2]->SetFillColor(2);
  var[2]->Draw("same");
  var[3]->SetFillColor(38);
  var[3]->Draw("same");
  var[4]->SetFillColor(46);
  var[4]->Draw("same");
  var[5]->SetFillColor(16);
  var[5]->Draw("same");
  var[6]->SetFillColor(4);
  var[6]->Draw("same");
  var[7]->SetFillColor(3);
  var[7]->Draw("same");
  var[8]->SetFillColor(64);
  var[8]->Draw("same");
  var[9]->SetFillColor(40);
  var[9]->Draw("same");
  var[10]->SetFillColor(95);
  var[10]->Draw("same");
  var[11]->SetFillColor(48);
  var[11]->Draw("same");
  var[12]->SetFillColor(52);
  var[12]->Draw("same");
  var[13]->SetFillColor(30);
  var[13]->Draw("same");

  leg->Draw();
  sprintf(name,"%s%s%s%s%s","results_gg/mc_",variableMC.c_str(),"_",allcut,".png");
  c0->SaveAs(name);

  //higgs only plot
  var_mc_2011[21]->Add(var_mc_2011[22]);
  var_mc_2011[21]->Add(var_mc_2011[23]);   
  var_mc_2011[21]->Add(var_mc_2011[24]);
  var_mc_2011[21]->Add(var_mc_2011[25]);
  var_mc_2011[21]->Add(var_mc_2011[26]);
  var_mc_2011[21]->Add(var_mc_2011[27]);
  sprintf(ytitle,"%s%d%s","N_{ev}/",int(int_exp_2011),"pb^{-1}");
  var_mc_2011[21]->SetXTitle(axis.c_str());
  var_mc_2011[21]->SetYTitle(ytitle);
  var_mc_2011[21]->SetTitle("");
  var_mc_2011[21]->SetLineColor(kBlack);
  var_mc_2011[21]->SetLineWidth(2);
  var_mc_2011[21]->SetFillColor(kYellow);
  var_mc_2011[21]->Draw();
  sprintf(name,"%s%s%s%s%s","results_gg/higgs_",variableData.c_str(),"_",allcut,".png");
  c0->SaveAs(name);

  // only VH contribution
  TH1D* myClone2 = (TH1D*) var_mc_2011[23]->Clone("myClone");
  myClone2->Add(var_mc_2011[24]);
  myClone2->Add(var_mc_2011[25]);
  myClone2->Add(var_mc_2011[26]);
  myClone2->Add(var_mc_2011[27]);

  double integralhiggs ;
  double entrieshiggs ;
  double integralbkg ;
  // some calculation for s/sqrt(b)
  if (variableMC == "massgg")
    {  
      sprintf(name,"%s%s%s%d%s","results_gg/optimalcut_",variableMC.c_str(),"_",eb,".txt");
      integralhiggs = var_mc_2011[21]->Integral(43,56);
      entrieshiggs =  var_mc_2011[21]->Integral(0,201);
      integralbkg = var[1]->Integral(30,71)/3.;
      cout << "Fraction in signal box " << integralhiggs/entrieshiggs << endl;
      cout << "Number of signal events " << integralhiggs/boosthiggs << endl;
      cout << "Number of bkg events " << integralbkg << endl;
      cout << "S/sqrt(B) " << integralhiggs/boosthiggs/sqrt(integralbkg) << endl;
  }  

  // here only signal HW, W->lnu                                                                                                   
  TH1D* myClone3 = (TH1D*) var_mc_2011[23]->Clone("myClone3");
  integralhiggs = myClone3->Integral(43,56);
  entrieshiggs  = myClone3->Integral(0,201);
  cout << endl;
  cout << "Analysis considering only W->lnu" << endl;
  cout << "Number of signal events " << integralhiggs/boosthiggs << endl;
  cout << "S/sqrt(B) "   << integralhiggs/boosthiggs/sqrt(integralbkg) << endl;
  cout << "S/B "         << integralhiggs/boosthiggs/integralbkg << endl;
  
  // here only signal HZ, Z->nn                                                                                                    
  TH1D* myClone4 = (TH1D*) var_mc_2011[27]->Clone("myClone4");
  integralhiggs = myClone4->Integral(43,56);
  entrieshiggs  = myClone4->Integral(0,201);
  cout << endl;
  cout << "Analysis considering only Z->nn" << endl;
  cout << "Number of signal events " << integralhiggs/boosthiggs << endl;
  cout << "S/sqrt(B) "   << integralhiggs/boosthiggs/sqrt(integralbkg) << endl;
  cout << "S/B "         << integralhiggs/boosthiggs/integralbkg << endl;

  // data only plot
  vardata->SetXTitle(axis.c_str());
  vardata->SetTitle("");
  vardata->SetStats(0);
  vardata->SetMarkerStyle(8);
  vardata->SetMarkerSize(.9);
  vardata->SetTitleOffset(1.25,"Y");
  vardata->Draw("pe");
  sprintf(name,"%s%s%s%s%s","results_gg/data_",variableData.c_str(),"_",allcut,".png");
  c0->SaveAs(name);

  // data overlaid to mc
  legge = leg->AddEntry(vardata, "data", "p");

  double themax =   vardata->GetMaximum();
  if(var[0]->GetMaximum()>themax) themax = var[0]->GetMaximum();
  if (
      variableMC == "etaphot1" || variableMC == "etaphot2" ||
      variableMC == "phiphot1" || variableMC == "phiphot2" ||
      variableMC == "etajet1" || variableMC == "etajet2" ||
      variableMC == "phijet1" || variableMC == "phijet2"
      )
    vardata->SetMaximum(themax*2.0);
  else
    vardata->SetMaximum(themax*1.1);
  vardata->SetMinimum(0);
  var[0]->Draw("same");
  var[1]->Draw("same");
  var[2]->Draw("same");
  var[3]->Draw("same");
  var[4]->Draw("same");
  var[5]->Draw("same");
  var[6]->Draw("same");
  var[7]->Draw("same");
  var[8]->Draw("same");
  var[9]->Draw("same");
  var[10]->Draw("same");
  var[11]->Draw("same");
  var[12]->Draw("same");
  var[13]->Draw("same");
  leg->Draw();
  vardata->Draw("pesame");
  gPad->RedrawAxis();
  sprintf(name,"%s%s%s%s%s","results_gg/data-mc_",variableMC.c_str(),"_",allcut,".png");
  c0->SaveAs(name);

  //data with control sample
  vardata->Draw("pe");
  vardatacs->SetLineColor(46);
  vardatacs->SetFillColor(42);
  vardatacs->SetLineWidth(3);
  vardatacs->Draw("hsame");
  vardata->Draw("pesame");
  sprintf(name,"%s%s%s%s%s","results_gg/datacs_",variableMC.c_str(),"_",allcut,".png");
  gPad->RedrawAxis();
  TLegendEntry *legge2;
  TLegend *leg2;
  leg2 = new TLegend(0.6,0.65,0.9,0.85);
  leg2->SetFillStyle(0); leg2->SetBorderSize(0); leg2->SetTextSize(0.035);
  leg2->SetFillColor(0);
  legge2 = leg2->AddEntry(vardata, "default sel.", "p");
  legge2 = leg2->AddEntry(vardatacs, "control sample", "f");
  leg2->Draw();
  c0->SaveAs(name);

  // additional plot with larger binnin (only for invariant mass)
  if (variableMC == "massgg"){
//     var[0]->Rebin(4);
//     themax = var[0]->GetMaximum();
//     vardata->Rebin(4);
//     cout << vardata->GetMaximum() << "   " << var[0]->GetMaximum() << endl;
//     if(themax < vardata->GetMaximum()) themax = vardata->GetMaximum();
//     var[0]->SetMaximum(themax*1.1);
//     var[0]->SetMinimum(0.);
//     var[0]->Draw();
//     var[1]->Rebin(4);
//     var[1]->Draw("same");
//     var[2]->Rebin(4);
//     var[2]->Draw("same");
//     var[3]->Rebin(4);
//     var[3]->Draw("same");
//     var[4]->Rebin(4);
//     var[4]->Draw("same");
//     var[5]->Rebin(4);
//     var[5]->Draw("same");
//     leg->Draw();
//     gPad->RedrawAxis();
//     sprintf(name,"%s%s%s%s%s","results_gg/mc_rebin_",variableMC.c_str(),"_",allcut,".png");
//     c0->SaveAs(name);
    
//     vardata->Draw("pesame");
//     sprintf(name,"%s%s%s%s%s","results_gg/data-mc_rebin_",variableMC.c_str(),"_",allcut,".png");
//     c0->SaveAs(name);

//     vardatacs->Rebin(4);
//     vardata->SetMaximum(themax*1.1);
//     vardata->Draw("pe");
//     vardatacs->Draw("same");
//     vardata->Draw("pesame");
//     leg2->Draw();
//     gPad->RedrawAxis();
//     sprintf(name,"%s%s%s%s%s","results_gg/datacs_rebin_",variableMC.c_str(),"_",allcut,".png");
//     c0->SaveAs(name);
    
//     double newmax(0);
//     for (int i=5;i<50;i++){
//       double tempmax = vardata->GetBinContent(i);
//       if(tempmax > newmax) newmax = tempmax;
//     }
//     vardata->SetAxisRange(100.,150.);
//     vardata->SetMaximum(newmax*1.3);
//     vardata->Draw("pe");
//     sprintf(name,"%s%s%s%s%s","results_gg/data_rebin_resize",variableMC.c_str(),"_",allcut,".png");
//     c0->SaveAs(name);

  

    sprintf(name,"%s%s%s","results_gg/yields_",allcut,".txt");
    ofstream outfile(name);  
    outfile << "####################################" << endl;
    outfile << "CUTS " << endl;
    outfile << "####################################" << endl;
    outfile << "ptphot1 : " << pt1 << endl;
    outfile << "ptphot2 : " << pt2 << endl;
    outfile << "ptjet1 : " << ptj1 << endl;
    outfile << "ptjet2 : " << ptj2 << endl;
    outfile << "met : " << misset  << endl;
    outfile << "deltaetacut : " << deltae << endl;
    outfile << "deltaphicut : " << deltap << endl;
    outfile << "deltaPhi(met,jet)cut : "   << jetmet << endl;
    outfile << "deltaPhi(met,phot1)cut : " << p1met << endl;
    outfile << "deltaPhi(met,phot2)cut : " << p2met << endl;
    outfile << "deltaPhi(met,higgs)cut : " << hmet << endl;
    outfile << "deltaPhi(phot1,phot2)cut : " << phigg << endl;
    outfile << "zeppencut : " << zep << endl;
    outfile << "invmassjetcut : " << mjj << endl;
    outfile << "CiC level : " << cic << endl;
    outfile << "ebcat : " << eb << endl;
    outfile << "r9cat : " << r9 << endl;
    outfile << "thirdcat : " << thirdcat << endl;
    outfile << "leptontag : " << leptontag << endl;
    outfile << "leptonveto : " << leptonveto << endl;
    outfile << endl;
    outfile << "####################################" << endl;
    outfile << "N of generated events" << endl;
    outfile << "####################################" << endl;
    outfile << "# events hig_glu2010 =    " << n_mc_2010[21] << endl;
    outfile << "# events hig_vbf2010 =    " << n_mc_2010[22] << endl;
    outfile << "# events hig_w2010, wlept =    " << n_mc_2010[23] << endl;
    outfile << "# events hig_w2010, whad  =    " << n_mc_2010[24] << endl;
    outfile << "# events hig_z2010, zlept =    " << n_mc_2010[25] << endl;
    outfile << "# events hig_z2010, zhad  =    " << n_mc_2010[26] << endl;
    outfile << "# events hig_z2010, zneut =    " << n_mc_2010[27] << endl;
    outfile << "# events dy_2010 =        "  << n_mc_2010[5] << endl;
    outfile << "# events box_2010 =       "  << n_mc_2010[0] << endl;
    outfile << "# events diphotjet_2010 = "  << n_mc_2010[1] << endl;
    outfile << "# events gjet_2010 =      "  << n_mc_2010[2] << endl;
    outfile << "# events qcd_2010 =       "  << n_mc_2010[3] << endl;
    outfile << "# events qcd2_2010 =      "  << n_mc_2010[4] << endl;
    outfile << "# events wgamma_ele_2010 = " << n_mc_2010[6] << endl;
    outfile << "# events wgamma_mu_2010  = " << n_mc_2010[7] << endl;
    outfile << "# events wgamma_tau_2010 = " << n_mc_2010[8] << endl;
    outfile << "# events zgamma_ele_2010 = " << n_mc_2010[9]  << endl;
    outfile << "# events zgamma_mu_2010  = " << n_mc_2010[10] << endl;
    outfile << "# events zgamma_tau_2010 = " << n_mc_2010[11] << endl;
    outfile << "# events ggWminus_2010 = "   << n_mc_2010[12] << endl;
    outfile << "# events ggWplus_2010 = "    << n_mc_2010[13] << endl;
    outfile << "# events ggZ_2010 = "        << n_mc_2010[14] << endl;
    outfile << "# events ggtt_2010 = "       << n_mc_2010[15] << endl;
    outfile << "# events ttjets_2010 = "     << n_mc_2010[16] << endl;
    outfile << "# events Wjets_2010 = "      << n_mc_2010[17] << endl;
    outfile << "# events WW_2010 = "         << n_mc_2010[18] << endl;
    outfile << "# events WZ_2010 = "         << n_mc_2010[19] << endl;
    outfile << "# events ZZ_2010 = "         << n_mc_2010[20] << endl;
    outfile << endl;
    outfile << "# events hig_glu2011 =    " << n_mc_2011[21] << endl;
    outfile << "# events hig_vbf2011 =    " << n_mc_2011[22] << endl;
    outfile << "# events hig_w2011, wlept =    " << n_mc_2011[23] << endl;
    outfile << "# events hig_w2011, whad  =    " << n_mc_2011[24] << endl;
    outfile << "# events hig_z2011, zlept =    " << n_mc_2011[25] << endl;
    outfile << "# events hig_z2011, zhad  =    " << n_mc_2011[26] << endl;
    outfile << "# events hig_z2011, zneut =    " << n_mc_2011[27] << endl;
    outfile << "# events dy_2011 =        "  << n_mc_2011[5] << endl;
    outfile << "# events box_2011 =       "  << n_mc_2011[0] << endl;
    outfile << "# events diphotjet_2011 = "  << n_mc_2011[1] << endl;
    outfile << "# events gjet_2011 =      "  << n_mc_2011[2] << endl;
    outfile << "# events qcd_2011 =       "  << n_mc_2011[3] << endl;
    outfile << "# events qcd2_2011 =      "  << n_mc_2011[4] << endl;
    outfile << "# events wgamma_ele_2011 = " << n_mc_2011[6] << endl;
    outfile << "# events wgamma_mu_2011  = " << n_mc_2011[7] << endl;
    outfile << "# events wgamma_tau_2011 = " << n_mc_2011[8] << endl;
    outfile << "# events zgamma_ele_2011 = " << n_mc_2011[9]  << endl;
    outfile << "# events zgamma_mu_2011  = " << n_mc_2011[10] << endl;
    outfile << "# events zgamma_tau_2011 = " << n_mc_2011[11] << endl;
    outfile << "# events ggWminus_2011 = "   << n_mc_2011[12] << endl;
    outfile << "# events ggWplus_2011 = "    << n_mc_2011[13] << endl;
    outfile << "# events ggZ_2011 = "        << n_mc_2011[14] << endl;
    outfile << "# events ggtt_2011 = "       << n_mc_2011[15] << endl;
    outfile << "# events ttjets_2011 = "     << n_mc_2011[16] << endl;
    outfile << "# events Wjets_2011 = "      << n_mc_2011[17] << endl;
    outfile << "# events WW_2011 = "         << n_mc_2011[18] << endl;
    outfile << "# events WZ_2011 = "         << n_mc_2011[19] << endl;
    outfile << "# events ZZ_2011 = "         << n_mc_2011[20] << endl;
    outfile << endl;
    outfile << "####################################" << endl;
    outfile << "N of selected events and eff." << endl;
    outfile << "####################################" << endl; 
    outfile << "ndata      = " << num_data << endl;
    outfile << endl;
    
    double num_bkg(0), err_num_bkg(0);
    double num_mc_total[28],num_uns_mc_total[28], n_mc_total[28];
    double err_num_mc_total[28],err_num_uns_mc_total[28];
    for (int i=0; i<28; i++){
      if(i<21){
	num_bkg += num_mc_2010[i];
	num_bkg += num_mc_2011[i];
      }
      num_mc_total[i] = num_mc_2010[i] + num_mc_2011[i];
      num_uns_mc_total[i] = num_uns_mc_2011[i];
      err_num_uns_mc_total[i] = sqrt(num_uns_mc_total[i]);
      if(num_uns_mc_total[i]) err_num_mc_total[i] = err_num_uns_mc_total[i] * num_mc_total[i]/num_uns_mc_total[i];
      else err_num_mc_total[i] = 0;
      n_mc_total[i] = n_mc_2010[i] * scale_mc_2010[i] + n_mc_2011[i] * scale_mc_2011[i];
    }
    for (int i=0; i<21; i++){
      err_num_bkg = sqrt(err_num_bkg*err_num_bkg + err_num_mc_total[i]*err_num_mc_total[i]);
    }

    outfile << "nallbkg    = " << num_bkg << " +/- " << err_num_bkg << endl;
    outfile << "nhig glu     = " << num_mc_total[21] << " +/- " << err_num_mc_total[21] << endl;
    outfile << "nhig vbf     = " << num_mc_total[22] << " +/- " << err_num_mc_total[22] << endl;
    outfile << "nhig wlept   = " << num_mc_total[23] << " +/- " << err_num_mc_total[23] << endl;
    outfile << "nhig whad    = " << num_mc_total[24] << " +/- " << err_num_mc_total[24] << endl;
    outfile << "nhig zlept   = " << num_mc_total[25] << " +/- " << err_num_mc_total[25] << endl;
    outfile << "nhig zhad    = " << num_mc_total[26] << " +/- " << err_num_mc_total[26] << endl;
    outfile << "nhig zneut   = " << num_mc_total[27] << " +/- " << err_num_mc_total[27] << endl;
    outfile << "ndy          = " << num_mc_total[5]  << " +/- " << err_num_mc_total[5]  << endl;
    outfile << "nbox         = " << num_mc_total[0]  << " +/- " << err_num_mc_total[0]  << endl;
    outfile << "ndiphot      = " << num_mc_total[1]  << " +/- " << err_num_mc_total[1]  << endl;
    outfile << "ngjet        = " << num_mc_total[2]  << " +/- " << err_num_mc_total[2]  << endl;
    outfile << "nqcd40       = " << num_mc_total[3]  << " +/- " << err_num_mc_total[3]  << endl;
    outfile << "nqcd30-40    = " << num_mc_total[4]  << " +/- " << err_num_mc_total[4]  << endl;
    outfile << "nwg-ele      = " << num_mc_total[6]  << " +/- " << err_num_mc_total[6]  << endl;
    outfile << "nwg-mu       = " << num_mc_total[7]  << " +/- " << err_num_mc_total[7]  << endl;
    outfile << "nwg-tau      = " << num_mc_total[8]  << " +/- " << err_num_mc_total[8]  << endl;
    outfile << "nzg-ele      = " << num_mc_total[9]  << " +/- " << err_num_mc_total[9]  << endl;
    outfile << "nzg-mu       = " << num_mc_total[10] << " +/- " << err_num_mc_total[10] << endl;
    outfile << "nzg-tau      = " << num_mc_total[11] << " +/- " << err_num_mc_total[11] << endl;
    outfile << "nggWminus    = " << num_mc_total[12] << " +/- " << err_num_mc_total[12] << endl;
    outfile << "nggWplus     = " << num_mc_total[13] << " +/- " << err_num_mc_total[13] << endl;
    outfile << "nggZ         = " << num_mc_total[14] << " +/- " << err_num_mc_total[14] << endl;
    outfile << "nggtt        = " << num_mc_total[15] << " +/- " << err_num_mc_total[15] << endl;
    outfile << "nttjets      = " << num_mc_total[16] << " +/- " << err_num_mc_total[16] << endl;
    outfile << "nWjets       = " << num_mc_total[17] << " +/- " << err_num_mc_total[17] << endl;
    outfile << "nWW          = " << num_mc_total[18] << " +/- " << err_num_mc_total[18] << endl;
    outfile << "nWZ          = " << num_mc_total[19] << " +/- " << err_num_mc_total[19] << endl;
    outfile << "nZZ          = " << num_mc_total[20] << " +/- " << err_num_mc_total[20] << endl;
    outfile << endl;
    outfile << "eff nhig     = "
            << (num_mc_total[21] + num_mc_total[22] + num_mc_total[23] + num_mc_total[24] + num_mc_total[25] + num_mc_total[26] + num_mc_total[27]) / (n_mc_total[21] + n_mc_total[22] + n_mc_total[23] + n_mc_total[24] + n_mc_total[25] + n_mc_total[26] + n_mc_total[27]) << endl;

    outfile << "eff nhig glu    = " << num_mc_total[21] / n_mc_total[21] << endl;
    outfile << "eff nhig vbf    = " << num_mc_total[22] / n_mc_total[22] << endl;
    outfile << "eff nhig wlept  = " << num_mc_total[23] / n_mc_total[23] << endl;
    outfile << "eff nhig whad   = " << num_mc_total[24] / n_mc_total[24] << endl;
    outfile << "eff nhig zlept  = " << num_mc_total[25] / n_mc_total[25] << endl;
    outfile << "eff nhig zhad   = " << num_mc_total[26] / n_mc_total[26] << endl;
    outfile << "eff nhig zneut  = " << num_mc_total[27] / n_mc_total[27] << endl;
    outfile << "eff ndy       = " << num_mc_total[5] / n_mc_total[5] << endl;
    outfile << "eff nbox      = " << num_mc_total[0] / n_mc_total[0] << endl;
    outfile << "eff ndiphot   = " << num_mc_total[1] / n_mc_total[1] << endl;
    outfile << "eff ngjet     = " << num_mc_total[2] / n_mc_total[2] << endl;
    outfile << "eff nqcd      = " << num_mc_total[3] / n_mc_total[3] << endl;
    outfile << "eff nqcd30-40 = " << num_mc_total[4] / n_mc_total[4] << endl;
    outfile << "eff wgamma-ele = " << num_mc_total[6] / n_mc_total[6] << endl;
    outfile << "eff wgamma-mu = "  << num_mc_total[7] / n_mc_total[7] << endl;
    outfile << "eff wgamma-tau = " << num_mc_total[8] / n_mc_total[8] << endl;
    outfile << "eff zgamma-ele = " << num_mc_total[9] / n_mc_total[9] << endl;
    outfile << "eff zgamma-mu = "  << num_mc_total[10] / n_mc_total[10] << endl;
    outfile << "eff zgamma-tau = " << num_mc_total[11] / n_mc_total[11] << endl;
    outfile << "eff ggWminus = " << num_mc_total[12] / n_mc_total[12] << endl;
    outfile << "eff ggWplus  = " << num_mc_total[13] / n_mc_total[13] << endl;
    outfile << "eff ggZ  = "     << num_mc_total[14] / n_mc_total[14] << endl;
    outfile << "eff ggtt  = "    << num_mc_total[15] / n_mc_total[15] << endl;
    outfile << "eff ttjets  = "  << num_mc_total[16] / n_mc_total[16] << endl;
    outfile << "eff Wjets  = "   << num_mc_total[17] / n_mc_total[17] << endl;
    outfile << "eff WW = "       << num_mc_total[18] / n_mc_total[18] << endl;
    outfile << "eff WZ = "       << num_mc_total[19] / n_mc_total[19] << endl;
    outfile << "eff ZZ = "       << num_mc_total[20] / n_mc_total[20] << endl;

    for(int i=0; i<7; i++){
      outfile << "eff nhig gluglu "<< h_masses[i] << " = " << num_gluglu_2011[i] / n_gluglu_2011[i] << endl;
      outfile << "eff nhig vbf    "<< h_masses[i] << " = " << num_vbf_2011[i] / n_vbf_2011[i] << endl;
      outfile << "eff nhig wzh    "<< h_masses[i] << " = " << num_wzh_2011[i] / n_wzh_2011[i] << endl;
//       outfile << "eff nhig tth    "<< h_masses[i] << " = " << num_tth_2011[i] / n_tth_2011[i] << endl;
    }
    outfile.close();

  }


 
  hOutputFile->Write() ;
  hOutputFile->Close() ;
  hOutputFile->Delete();

//   delete c0;

  delete data;
  delete datacs;

  for(int i=0; i<28; i++){
    delete mc_2011_fill[i];
    delete mc_2011[i];
  }

  for(int i=0; i<7; i++){
    delete mc_gluglu_2011_fill[i];
    delete mc_vbf_2011_fill[i];
    delete mc_wzh_2011_fill[i];
    //    delete mc_tth_2011_fill[i];
    delete mc_gluglu_2011[i];
    delete mc_vbf_2011[i];
    delete mc_wzh_2011[i];
    //    delete mc_tth_2011[i];
  }

//   delete leg;
//   delete leg2;

  vector<double> values;

  if (variableMC == "massgg"){
    values.push_back(integralhiggs);
    values.push_back(integralbkg);
    values.push_back(integralhiggs/sqrt(integralbkg));
  }

  return values;

}








