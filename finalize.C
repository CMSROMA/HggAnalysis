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
// .x finalize(data int lumi 2010, data int lumi 2011, pt1 cut, pt2 cut, ptj1 cut, ptj2 cut, deltae cut, zep cut, mjj cut, eb cat, r9 cat, trk isocut scale, cic selection, no pixel)
// 
// example:
//
// .x finalize.C(33,230,40,25,20,15,2.5,2.5,300,1,1,4,1) 

vector <double> finalize(double int_exp_2010, double int_exp_2011, double pt1=50, double pt2=30, double pthiggsmin = -100, double pthiggsmax = -100, double ptj1=20, double ptj2=15, double deltae=2.5, double zep=2.5, double mjj=300, int eb = 1, int r9 = 1, int cic = 4, bool pixel = 1, string variable = "massgg", int nbin = 200, double min = 90, double max = 190, string axis = "m(#gamma#gamma)[GeV]", int isolscaletrk = 100., int isolscaleecal = 100., int isolscalehcal = 100., int isolscalehove = 100.){


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
  // legenda for mc inputs 
  // 0 = dy
  // 1 = box
  // 2 = diphoton
  // 3 = gjet
  // 4 = qcd pt>40 
  // 5 = qcd 30<pt<40
  // 6 = higgs gluglu
  // 7 = higgs VBF
  // 8 = higgs WZH    

  string mcnames[9];
  mcnames[0] = "box";
  mcnames[1] = "diphoton";
  mcnames[2] = "gjet";
  mcnames[3] = "qcdpt>40";
  mcnames[4] = "qcd30<pt<40";
  mcnames[5] = "dy";
  mcnames[6] = "higgsgluglu";
  mcnames[7] = "higgsVBF";
  mcnames[8] = "higgsWZH";

  TFile* mc_2010[9];
  TFile* mc_2011[9];

  TFile* mc_gluglu_2011[7];
  TFile* mc_vbf_2011[7];
  TFile* mc_wzh_2011[7];
  TFile* mc_tth_2011[7];

  int h_masses[7] = {100,105,110,115,120,130,140};

  TString redntpDir= "/shome/meridian/software/CMSSW423/src/Analysis/Higgs";

  TString preselectionLevel;


//    if (cic>0)
//      preselectionLevel="cicloose";
//   else
    preselectionLevel="preselectionCS";

  TString preselectionLevelCS="preselectionCS";
  // total data sample
  //TFile* data = new TFile(redntpDir+"/redntp.42xv1_data."+preselectionLevel+".v2/merged/redntp_Photon-Run2011A-DiPhoton-May10ReReco-v1.root");    
  TFile* datacs = new TFile(redntpDir+"/redntp.42xv1_data."+preselectionLevelCS+".v2/merged/redntp_Photon-Run2011A-PromptReco-v4-DiPhotonSkimOnFly-v5_DiPhoton-May10ReReco-v1.root");    
  //  TFile* data = new TFile(redntpDir+"/redntp.41xv10_data."+preselectionLevel+".v1/merged/redntp_Run2011A-May7ReReco-v1-DiPhotonSkim.root");    
  TFile* data = new TFile(redntpDir+"/redntp.42xv1_data."+preselectionLevel+".v2/merged/redntp_Photon-Run2011A-PromptReco-v4-DiPhotonSkimOnFly-v5_DiPhoton-May10ReReco-v1.root");    

  if(int_exp_2010>0){
    // box samples
    mc_2010[0] = new TFile("redntp.39xv7.preselection.v2/redntp_DiPhotonBox_Pt25to250_TrackingParticles_7TeV-pythia6_00.root");
    // diphoton jets samples
    mc_2010[1] = new TFile("redntp.39xv7.preselection.v2/redntp_DiPhotonJets_7TeV-madgraph.root");
    // gjet samples
    mc_2010[2] = new TFile("redntp.39xv7.preselection.v2/redntp_GJet_Pt-20_doubleEMEnriched_TuneZ2_7TeV-pythia6.root");
    // qcd pt>40 samples
    mc_2010[3] = new TFile("redntp.39xv7.preselection.v2/redntp_QCD_Pt-40_doubleEMEnriched_TuneZ2_7TeV-pythia6.root");
    // qcd 30<pt<40 samples
    mc_2010[4] = new TFile("redntp.39xv7.preselection.v2/redntp_QCD_Pt-30to40_doubleEMEnriched_TuneZ2_7TeV-pythia6.root");
    // drell yan samples
    mc_2010[5] = new TFile("redntp.39xv7.preselection.v2/redntp_DYJetsToLL_TuneZ2_M-50_7TeV-madgraph-tauola.root"); 
    // gluglu higgs samples 
    mc_2010[6] = new TFile("redntp.39xv7.preselection.v2/redntp_GluGluToHToGG_M-115_7TeV-pythia6_00.root");
    // vbf higgs samples 
    mc_2010[7] = new TFile("redntp.39xv7.preselection.v2/redntp_VBF_HToGG_M-105_00.root");
    // W/Z/TT H higgs samples 
    mc_2010[8] = new TFile("redntp.39xv7.preselection.v2/redntp_WH_ZH_TTH_HToGG_M-115_00.root");
  }

  if(int_exp_2011>0){
    // box samples
    mc_2011[0] = new TFile(redntpDir+"/redntp.41xv10."+preselectionLevel+".v2/merged/redntp_DiPhotonBox_Pt25to250-41x_ntpv4.root");
    // diphoton jets samples
    mc_2011[1] = new TFile(redntpDir+"/redntp.41xv10."+preselectionLevel+".v2/merged/redntp_DiPhotonJets_7TeV-madgraph.root");
    // gjet samples
    mc_2011[2] = new TFile(redntpDir+"/redntp.41xv10."+preselectionLevel+".v2/merged/redntp_GJet_Pt-20_doubleEMEnriched_TuneZ2_7TeV-pythia6-41x_ntpv1.root");
    // qcd pt>40 samples
    mc_2011[3] = new TFile(redntpDir+"/redntp.41xv10."+preselectionLevel+".v2/merged/redntp_QCD_Pt-40_doubleEMEnriched_TuneZ2_7TeV-pythia6-41x_ntpv1.root");
    // qcd 30<pt<40 samples
    mc_2011[4] = new TFile(redntpDir+"/redntp.41xv10."+preselectionLevel+".v2/merged/redntp_QCD_Pt-30to40_doubleEMEnriched_TuneZ2_7TeV-pythia6-41x_ntpv1.root");
    // drell yan samples
    mc_2011[5] = new TFile(redntpDir+"/redntp.41xv10."+preselectionLevel+".v2/merged/redntp_DYJetsToLL_TuneZ2_M-50_7TeV-madgraph-tauola-PU-winter-newtrkiso-41x_ntpv1.root");
    // gluglu higgs samples 
    mc_2011[6] = new TFile(redntpDir+"/redntp.42xv2."+preselectionLevel+".v2/merged/redntp_GluGluToHToGG_M-115_7TeV-powheg-pythia6.root");
    // vbf higgs samples 
    mc_2011[7] = new TFile(redntpDir+"/redntp.42xv2."+preselectionLevel+".v2/merged/redntp_VBF_HToGG_M-115_7TeV-powheg-pythia6.root");
    // W/Z/TT H higgs samples 
    mc_2011[8] = new TFile(redntpDir+"/redntp.42xv2."+preselectionLevel+".v2/merged/redntp_WH_ZH_HToGG_M-115_7TeV-pythia6.root");
  }
 
  for (int i=0;i<7;i++){
    char hmass[100]; sprintf(hmass,"%d",h_masses[i]);
    // gluglu samples
    mc_gluglu_2011[i] = new TFile(redntpDir+"/redntp.42xv2."+preselectionLevel+".v2/merged/redntp_GluGluToHToGG_M-"+ hmass + "_7TeV-powheg-pythia6.root");
    // vbf higgs samples 
    mc_vbf_2011[i] = new TFile(redntpDir+"/redntp.42xv2."+preselectionLevel+".v2/merged/redntp_VBF_HToGG_M-"+ hmass +"_7TeV-powheg-pythia6.root");
    // W/Z H higgs samples 
    mc_wzh_2011[i] = new TFile(redntpDir+"/redntp.42xv2."+preselectionLevel+".v2/merged/redntp_WH_ZH_HToGG_M-"+ hmass +"_7TeV-pythia6.root");
//     // TT H higgs samples 
//     mc_tth_2011[i] = new TFile(redntpDir+"/redntp.42xv2."+preselectionLevel+".v2/merged/redntp_TTH_HToGG_M-"+ hmass +"_7TeV-pythia6.root");
  }

  // cross sections and scaling
  double boosthiggs(15);
  double cross_mc[9];
  cross_mc[0] = 12.37; // box
  cross_mc[1] = 134; // diphoton jets
  cross_mc[2] = 493.44; // gjet
  cross_mc[3] = 40392; // qcd pt>40
  cross_mc[4] = 9610; // qcd 30<pt<40 
  cross_mc[5] = 2321; // drell yan
  cross_mc[6] = 18.23 * 2.13e-03 * boosthiggs; // glu glu higgs
  cross_mc[7] = 1.332 * 2.13e-03 * boosthiggs; // vbf higgs
  cross_mc[8] = (0.7546 + 0.4107) * 2.13e-03 * boosthiggs; // WHtt higgs

  // getting the number of original events in each sample (processed with CMSSW)
  int n_mc_2010[9], n_mc_2011[9],n_gluglu_2011[7],n_vbf_2011[7],n_wzh_2011[7],n_tth_2011[7];
  for(int i=0; i<9; i++){
    n_mc_2010[i] = n_mc_2011[i] = 0;
    if(int_exp_2010>0) n_mc_2010[i] = ((TH1D*)mc_2010[i]->Get("ptphotgen1"))->GetEntries();
    if(int_exp_2011>0) n_mc_2011[i] = ((TH1D*)mc_2011[i]->Get("ptphotgen1"))->GetEntries();
  }
  for(int i=0; i<7; i++){
    n_gluglu_2011[i] = n_vbf_2011[i] = n_wzh_2011[i] = n_tth_2011[i] = 0;
    n_gluglu_2011[i] = ((TH1D*)mc_gluglu_2011[i]->Get("ptphotgen1"))->GetEntries();
    n_vbf_2011[i] = ((TH1D*)mc_vbf_2011[i]->Get("ptphotgen1"))->GetEntries();
    n_wzh_2011[i] = ((TH1D*)mc_wzh_2011[i]->Get("ptphotgen1"))->GetEntries();
//     n_tth_2011[i] = ((TH1D*)mc_tth_2011[i]->Get("ptphotgen1"))->GetEntries();
  }

  // setting the scaling factor to actual lumi
  double scale_mc_2010[9], scale_mc_2011[9];
  for(int i=0; i<9; i++){
    scale_mc_2010[i] = scale_mc_2011[i] = 0; 
    if(int_exp_2010>0) scale_mc_2010[i] = cross_mc[i] * int_exp_2010 / n_mc_2010[i];
    if(int_exp_2011>0) scale_mc_2011[i] = cross_mc[i] * int_exp_2011 / n_mc_2011[i];
  }

  // char for output name
  char name[1000];
  char allcut[3000];
  sprintf(allcut,"%3.1f%s%3.1f%s%3.1f%s%3.1f%s%3.1f%s%3.1f%s%3.1f%s%3.1f%s%3.1f%s%d%s%d%s%d%s%d%s%d%s%d%s%d%s%d",pt1,"_",pt2,"_",pthiggsmin,"_",pthiggsmax,"_",ptj1,"_",ptj2,"_",deltae,"_",zep,"_",mjj,"_",eb,"_",r9,"_",isolscaletrk,"_",isolscaleecal,"_",isolscalehcal,"_",isolscalehove,"_",pixel,"_",cic);

  // output root file
  sprintf(name,"%s%s%s%s%s","results_gg/histo_",variable.c_str(),"_",allcut,".root");
  TFile * hOutputFile   = new TFile(name, "RECREATE" ) ;

  // histograms needed by the machinery
  TH1D* vardata = new TH1D("vardata","vardata",nbin,min,max);
  TH1D* vardatacs = new TH1D("vardatacs","vardatacs",nbin,min,max);
  TH1D* var_mc_2010[9];
  TH1D* var_mc_2011[9];
  TH1D* var_gluglu_2011[7];
  TH1D* var_vbf_2011[7];
  TH1D* var_wzh_2011[7];
  TH1D* var_tth_2011[7];
  TH1D * var[6];
  for (int i=0; i<6; i++) {
    sprintf(name,"%s%d","var",i);
    var[i] = new TH1D(name,name,nbin,min,max);
  }
  for (int i=0; i<9; i++) {
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
  data_fill.Setcuts(pt1,pt2,pthiggsmin,pthiggsmax,ptj1,ptj2,deltae,zep,mjj,eb,r9,isolscaletrk,isolscaleecal,isolscalehcal,isolscalehove,pixel);
  datacs_fill.Setcuts(pt1,pt2,pthiggsmin,pthiggsmax,ptj1,ptj2,deltae,zep,mjj,eb,r9,isolscaletrk,isolscaleecal,isolscalehcal,isolscalehove,pixel);
  if(cic>0)
    {
      data_fill.setCic(cic);
      datacs_fill.setCic(cic);
    }

  fillPlot* mc_2010_fill[9];
  fillPlot* mc_2011_fill[9];
  fillPlot* mc_gluglu_2011_fill[9];  
  fillPlot* mc_vbf_2011_fill[9];  
  fillPlot* mc_wzh_2011_fill[9];  
  fillPlot* mc_tth_2011_fill[9];  

  for (int i=0; i<9; i++){
    if(int_exp_2010>0) mc_2010_fill[i] = new fillPlot((TTree*)mc_2010[i]->Get("AnaTree"), 1);
    if(int_exp_2011>0) mc_2011_fill[i] = new fillPlot((TTree*)mc_2011[i]->Get("AnaTree"), 1);
    if(int_exp_2010>0) mc_2010_fill[i]->Setcuts(pt1,pt2,pthiggsmin,pthiggsmax,ptj1,ptj2,deltae,zep,mjj,eb,r9,isolscaletrk,isolscaleecal,isolscalehcal,isolscalehove,pixel);
    if(int_exp_2011>0) mc_2011_fill[i]->Setcuts(pt1,pt2,pthiggsmin,pthiggsmax,ptj1,ptj2,deltae,zep,mjj,eb,r9,isolscaletrk,isolscaleecal,isolscalehcal,isolscalehove,pixel);

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
    mc_gluglu_2011_fill[i]->Setcuts(pt1,pt2,pthiggsmin,pthiggsmax,ptj1,ptj2,deltae,zep,mjj,eb,r9,isolscaletrk,isolscaleecal,isolscalehcal,isolscalehove,pixel);
    mc_vbf_2011_fill[i]->Setcuts(pt1,pt2,pthiggsmin,pthiggsmax,ptj1,ptj2,deltae,zep,mjj,eb,r9,isolscaletrk,isolscaleecal,isolscalehcal,isolscalehove,pixel);
    mc_wzh_2011_fill[i]->Setcuts(pt1,pt2,pthiggsmin,pthiggsmax,ptj1,ptj2,deltae,zep,mjj,eb,r9,isolscaletrk,isolscaleecal,isolscalehcal,isolscalehove,pixel);
//     mc_tth_2011_fill[i]->Setcuts(pt1,pt2,pthiggsmin,pthiggsmax,ptj1,ptj2,deltae,zep,mjj,eb,r9,isolscaletrk,isolscaleecal,isolscalehcal,isolscalehove,pixel);
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
//   for (int i=0; i<9; i++){
//     if(int_exp_2010>0) mc_2010_fill[i]->DoSmearing(1.,0.0001);
//     if(int_exp_2011>0) mc_2011_fill[i]->DoSmearing(1.,0.0001);
//   }

  // filling histograms
  std::cout << " ++++++++++++++ DATA ++++++++++++++++" << std::endl;
  cout << "running over " << ((TTree*)data->Get("AnaTree"))->GetEntries("") << " data events" <<  endl;
  if (variable == "massgg") {  
    sprintf(name,"%s%s%s","results_gg/events_",allcut,".txt");
    data_fill.Writetxt(name);
    sprintf(name,"%s%s%s","results_gg/events_",allcut,".root");
    data_fill.WriteRoot(name);
  }
  vardata->Add(data_fill.Plot(variable,"data", nbin, min, max)); 
  std::cout << "Selected events on data " << vardata->GetEntries() << std::endl;
  cout << "running over " << ((TTree*)datacs->Get("AnaTree"))->GetEntries("") << " data events (for cs)" <<  endl; 

  if (variable == "massgg") {  
    datacs_fill.Writetxt(name);
    datacs_fill.WriteRoot(name);
  }
  vardatacs->Add(datacs_fill.Plot(variable,"datacs", nbin, min, max, 1)); 
  std::cout << "Selected events on data cs " << vardatacs->GetEntries() << std::endl;

  std::cout << " ++++++++++++++ MC ++++++++++++++++" << std::endl;
  for (int i=0; i<9; i++){ 
    sprintf(name,"%s%s",mcnames[i].c_str()," 2010");
    if(int_exp_2010>0) {
      cout << "running over " << ((TTree*)mc_2010[i]->Get("AnaTree"))->GetEntries("") << " " << name << " events" <<  endl; 
      sprintf(name,"%s%s%s%s%s","results_gg/events_",mcnames[i].c_str(),"_2010_",allcut,".root");
      if (variable == "massgg") mc_2010_fill[i]->WriteRoot(name);
      var_mc_2010[i]->Add(mc_2010_fill[i]->Plot(variable, name, nbin, min, max));
      std::cout << "Selected events on mc2010 " << name << " " << var_mc_2010[i]->GetEntries() << std::endl;
    }
    sprintf(name,"%s%s",mcnames[i].c_str()," 2011");
    if(int_exp_2011>0) {
      cout << "running over " << ((TTree*)mc_2011[i]->Get("AnaTree"))->GetEntries("") << " " << name << " events" <<  endl; 
      sprintf(name,"%s%s%s%s%s","results_gg/events_",mcnames[i].c_str(),"_2011_",allcut,".root");
      if (variable == "massgg") mc_2011_fill[i]->WriteRoot(name);
      var_mc_2011[i]->Add(mc_2011_fill[i]->Plot(variable, name, nbin, min, max));
      std::cout << "Selected events on mc2011 " << name << " " << var_mc_2011[i]->GetEntries() << std::endl;
    }

  }

  std::cout << " ++++++++++++++ signal MC ++++++++++++++++" << std::endl;
  for (int i=0; i<7; i++){ 
    cout << "running over " << ((TTree*)mc_gluglu_2011[i]->Get("AnaTree"))->GetEntries("") << " gluglu M=" << h_masses[i] << " events" <<  endl; 
    sprintf(name,"%s%d%s%s%s","results_gg/events_gluglu",h_masses[i],"_2011_",allcut,".root");
    if (variable == "massgg") mc_gluglu_2011_fill[i]->WriteRoot(name);
    var_gluglu_2011[i]->Add(mc_gluglu_2011_fill[i]->Plot(variable, name, nbin, min, max));
    std::cout << "Selected events on mc2011 gluglu " << h_masses[i] << " " << var_gluglu_2011[i]->GetEntries() << std::endl;
 
    cout << "running over " << ((TTree*)mc_vbf_2011[i]->Get("AnaTree"))->GetEntries("") << " vbf M=" << h_masses[i] << " events" <<  endl; 
    sprintf(name,"%s%d%s%s%s","results_gg/events_vbf",h_masses[i],"_2011_",allcut,".root");
    if (variable == "massgg") mc_vbf_2011_fill[i]->WriteRoot(name);
    var_vbf_2011[i]->Add(mc_vbf_2011_fill[i]->Plot(variable, name, nbin, min, max));
    std::cout << "Selected events on mc2011 vbf " << h_masses[i] << " " << var_vbf_2011[i]->GetEntries() << std::endl;

    cout << "running over " << ((TTree*)mc_wzh_2011[i]->Get("AnaTree"))->GetEntries("") << " wzh M=" << h_masses[i] << " events" <<  endl; 
    sprintf(name,"%s%d%s%s%s","results_gg/events_wzh",h_masses[i],"_2011_",allcut,".root");
    if (variable == "massgg") mc_wzh_2011_fill[i]->WriteRoot(name);
    var_wzh_2011[i]->Add(mc_wzh_2011_fill[i]->Plot(variable, name, nbin, min, max));
    std::cout << "Selected events on mc2011 wzh " << h_masses[i] << " " << var_wzh_2011[i]->GetEntries() << std::endl;

//     cout << "running over " << ((TTree*)mc_tth_2011[i]->Get("AnaTree"))->GetEntries("") << " tth M=" << h_masses[i] << " events" <<  endl; 
//     sprintf(name,"%s%d%s%s%s","results_gg/events_tth",h_masses[i],"_2011_",allcut,".root");
//     if (variable == "massgg") mc_tth_2011_fill[i]->WriteRoot(name);
//     var_tth_2011[i]->Add(mc_tth_2011_fill[i]->Plot(variable, name, nbin, min, max));
//     std::cout << "Selected events on mc2011 tth " << h_masses[i] << " " << var_tth_2011[i]->GetEntries() << std::endl;
 }

  // scale mc to equivalent lumi
  for (int i=0; i<9; i++){ 
    if(int_exp_2010>0) var_mc_2010[i]->Scale(scale_mc_2010[i]);  
    if(int_exp_2011>0) var_mc_2011[i]->Scale(scale_mc_2011[i]);  
  }

  // counting number of events passing selection (scaled)
  double num_mc_2010[9],num_mc_2011[9],num_gluglu_2011[7],num_vbf_2011[7],num_wzh_2011[7],num_tth_2011[7]; 
  for (int i=0; i<9; i++){ 
    num_mc_2010[i] = num_mc_2011[i] = 0;
    if(int_exp_2010>0) num_mc_2010[i] = var_mc_2010[i]->Integral();  
    if(int_exp_2011>0) num_mc_2011[i] = var_mc_2011[i]->Integral();  
  }
  for(int i=0; i<7; i++){
    num_gluglu_2011[i] = num_vbf_2011[i] = num_wzh_2011[i] = num_tth_2011[i] = 0;
    num_gluglu_2011[i] = num_gluglu_2011[i] = var_gluglu_2011[i]->Integral(); 
    num_vbf_2011[i] = num_vbf_2011[i] = var_vbf_2011[i]->Integral(); 
    num_wzh_2011[i] = num_wzh_2011[i] = var_wzh_2011[i]->Integral(); 
//     num_tth_2011[i] = num_tth_2011[i] = var_tth_2011[i]->Integral(); 
  }

  // add two QCD bins
  if(int_exp_2010>0) var_mc_2010[3]->Add(var_mc_2010[4]);
  if(int_exp_2011>0) var_mc_2011[3]->Add(var_mc_2011[4]);
  
  // scale control sample
  vardata->Sumw2();
  //  vardatacs->Sumw2();
  double num_data =  vardata->Integral();
  double num_data_cs = vardatacs->Integral();  
  vardatacs->Scale(num_data/num_data_cs); 

  // stack histograms  
  for (int i=1; i<nbin+1; i++){      
    for (int j=0; j<6; j++){            
      int offset(0);
      if(j>0) offset = 2; // to add higgs contributions up
      if(j>2) offset = 3; // to add qcd contributions up
      for (int k=0 ; k<9-j-offset; k++){ 
	if(int_exp_2010>0) var[j]->SetBinContent(i,var_mc_2010[k]->GetBinContent(i) + var[j]->GetBinContent(i));
	if(int_exp_2011>0) var[j]->SetBinContent(i,var_mc_2011[k]->GetBinContent(i) + var[j]->GetBinContent(i));
      }	
    }    
  }
 
  //final plots
  char ytitle[100];
  sprintf(ytitle,"%s%d%s","N_{ev}/",int(int_exp_2010+int_exp_2011),"pb^{-1}");
  for(int i=0; i<6; i++){
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
  leg = new TLegend(0.5,0.6,0.75,0.85);
  leg->SetFillStyle(0); leg->SetBorderSize(0); leg->SetTextSize(0.05);
  leg->SetFillColor(0);
  sprintf(name,"H M(115)");
  if(boosthiggs!=1) sprintf(name,"%s%2.0f", "H M(115) x ", boosthiggs);
  legge = leg->AddEntry(var[0], name, "f");
  legge = leg->AddEntry(var[1], "DY", "f");
  legge = leg->AddEntry(var[2], "QCD", "f");
  legge = leg->AddEntry(var[3], "#gamma + jets", "f");
  legge = leg->AddEntry(var[4], "di-#gamma + jets", "f");
  legge = leg->AddEntry(var[5], "di-#gamma box", "f");

  //mc only plot
  var[0]->SetFillColor(kYellow);
  var[0]->Draw();
  var[1]->SetFillColor(32);
  var[1]->Draw("same");
  var[2]->SetFillColor(29);
  var[2]->Draw("same");
  var[3]->SetFillColor(38);
  var[3]->Draw("same");
  var[4]->SetFillColor(46);
  var[4]->Draw("same");
  var[5]->SetFillColor(16);
  var[5]->Draw("same");
  leg->Draw();
  sprintf(name,"%s%s%s%s%s","results_gg/mc_",variable.c_str(),"_",allcut,".gif");
  c0->SaveAs(name);

  //higgs only plot
  var_mc_2011[6]->Add(var_mc_2011[7]);
  var_mc_2011[6]->Add(var_mc_2011[8]);
  sprintf(ytitle,"%s%d%s","N_{ev}/",int(int_exp_2010+int_exp_2011),"pb^{-1}");
  var_mc_2011[6]->SetXTitle(axis.c_str());
  var_mc_2011[6]->SetYTitle(ytitle);
  var_mc_2011[6]->SetTitle("");
  var_mc_2011[6]->SetLineColor(kBlack);
  var_mc_2011[6]->SetLineWidth(2);
  var_mc_2011[6]->SetFillColor(kYellow);
  var_mc_2011[6]->Draw();
  sprintf(name,"%s%s%s%s%s","results_gg/higgs_",variable.c_str(),"_",allcut,".gif");
  c0->SaveAs(name);
  
  double integralhiggs ;
  double entrieshiggs ;
  double integralbkg ;
  // some calculation for s/sqrt(b)
  if (variable == "massgg")
    {  
      sprintf(name,"%s%s%s%d%s","results_gg/optimalcut_",variable.c_str(),"_",eb,".txt");
      integralhiggs = var_mc_2011[6]->Integral(43,56);
      entrieshiggs =  var_mc_2011[6]->Integral(0,201);
      integralbkg = var[1]->Integral(30,71)/3.;
      cout << "Fraction in signal box " << integralhiggs/entrieshiggs << endl;
      cout << "Number of signal events " << integralhiggs/boosthiggs << endl;
      cout << "Number of bkg events " << integralbkg << endl;
      cout << "S/sqrt(B) " << integralhiggs/boosthiggs/sqrt(integralbkg) << endl;
  }  

  // data only plot
  vardata->SetXTitle(axis.c_str());
  vardata->SetTitle("");
  vardata->SetStats(0);
  vardata->SetMarkerStyle(8);
  vardata->SetMarkerSize(.9);
  vardata->SetTitleOffset(1.25,"Y");
  vardata->Draw("pe");
  sprintf(name,"%s%s%s%s%s","results_gg/data_",variable.c_str(),"_",allcut,".gif");
  c0->SaveAs(name);

  // data overlaid to mc
  double themax =   vardata->GetMaximum();
  if(var[0]->GetMaximum()>themax) themax = var[0]->GetMaximum();
  if (
      variable == "etaphot1" || variable == "etaphot2" ||
      variable == "phiphot1" || variable == "phiphot2" ||
      variable == "etajet1" || variable == "etajet2" ||
      variable == "phijet1" || variable == "phijet2"
      )
    vardata->SetMaximum(themax*2.0);
  else
    vardata->SetMaximum(themax*1.1);
  var[0]->Draw("same");
  var[1]->Draw("same");
  var[2]->Draw("same");
  var[3]->Draw("same");
  var[4]->Draw("same");
  var[5]->Draw("same");
  leg->Draw();
  vardata->Draw("pesame");
  gPad->RedrawAxis();
  sprintf(name,"%s%s%s%s%s","results_gg/data-mc_",variable.c_str(),"_",allcut,".gif");
  c0->SaveAs(name);

  //data with control sample
  vardata->Draw("pe");
  vardatacs->SetLineColor(46);
  vardatacs->SetFillColor(42);
  vardatacs->SetLineWidth(3);
  vardatacs->Draw("hsame");
  vardata->Draw("pesame");
  sprintf(name,"%s%s%s%s%s","results_gg/datacs_",variable.c_str(),"_",allcut,".gif");
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
  if (variable == "massgg"){
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
//     sprintf(name,"%s%s%s%s%s","results_gg/mc_rebin_",variable.c_str(),"_",allcut,".gif");
//     c0->SaveAs(name);
    
//     vardata->Draw("pesame");
//     sprintf(name,"%s%s%s%s%s","results_gg/data-mc_rebin_",variable.c_str(),"_",allcut,".gif");
//     c0->SaveAs(name);

    vardatacs->Rebin(4);
    vardata->SetMaximum(themax*1.1);
    vardata->Draw("pe");
    vardatacs->Draw("same");
    vardata->Draw("pesame");
    leg2->Draw();
    gPad->RedrawAxis();
    sprintf(name,"%s%s%s%s%s","results_gg/datacs_rebin_",variable.c_str(),"_",allcut,".gif");
    c0->SaveAs(name);
    
    double newmax(0);
    for (int i=5;i<50;i++){
      double tempmax = vardata->GetBinContent(i);
      if(tempmax > newmax) newmax = tempmax;
    }
    vardata->SetAxisRange(100.,150.);
    vardata->SetMaximum(newmax*1.3);
    vardata->Draw("pe");
    sprintf(name,"%s%s%s%s%s","results_gg/data_rebin_resize",variable.c_str(),"_",allcut,".gif");
    c0->SaveAs(name);

  

    sprintf(name,"%s%s%s","results_gg/yields_",allcut,".txt");
    ofstream outfile(name);  
    outfile << "####################################" << endl;
    outfile << "CUTS " << endl;
    outfile << "####################################" << endl;
    outfile << "ptphot1 : " << pt1 << endl;
    outfile << "ptphot2 : " << pt2 << endl;
    outfile << "ptjet1 : " << ptj1 << endl;
    outfile << "ptjet2 : " << ptj2 << endl;
    outfile << "deltaetacut : " << deltae << endl;
    outfile << "zeppencut : " << zep << endl;
    outfile << "invmassjetcut : " << mjj << endl;
    outfile << "pixelseedcut : " << pixel << endl;
    outfile << "CiC level : " << cic << endl;
    outfile << "ebcat : " << eb << endl;
    outfile << "r9cat : " << r9 << endl;
    outfile << "scaletrk : " << isolscaletrk << endl;
    outfile << "scaleecal : " << isolscaleecal << endl;
    outfile << "scalehcal : " << isolscalehcal << endl;
    outfile << "scalehove : " << isolscalehove << endl;
    outfile << endl;
    outfile << "####################################" << endl;
    outfile << "N of generated events" << endl;
    outfile << "####################################" << endl;
    outfile << "# events hig_glu2010 =    " << n_mc_2010[6] << endl;
    outfile << "# events hig_vbf2010 =    " << n_mc_2010[7] << endl;
    outfile << "# events hig_wzt2010 =    " << n_mc_2010[8] << endl;
    outfile << "# events dy_2010 =        " << n_mc_2010[5] << endl;
    outfile << "# events box_2010 =       " << n_mc_2010[0] << endl;
    outfile << "# events diphotjet_2010 = " << n_mc_2010[1] << endl;
    outfile << "# events gjet_2010 =      " << n_mc_2010[2] << endl;
    outfile << "# events qcd_2010 =       " << n_mc_2010[3] << endl;
    outfile << "# events qcd2_2010 =      " << n_mc_2010[4] << endl;
    outfile << "# events hig_2011 =       " << n_mc_2011[6] << endl;
    outfile << "# events hig_vbf2011 =    " << n_mc_2011[7] << endl;
    outfile << "# events hig_wzt2011 =    " << n_mc_2011[8] << endl;
    outfile << "# events dy_2011 =        " << n_mc_2011[5] << endl;
    outfile << "# events box_2011 =       " << n_mc_2011[0] << endl;
    outfile << "# events diphotjet_2011 = " << n_mc_2011[1] << endl;
    outfile << "# events gjet_2011 =      " << n_mc_2011[2] << endl;
    outfile << "# events qcd_2011 =       " << n_mc_2011[3] << endl;
    outfile << "# events qcd2_2011 =      " << n_mc_2011[4] << endl;
    outfile << endl;
    outfile << "####################################" << endl;
    outfile << "N of selected events and eff." << endl;
    outfile << "####################################" << endl; 
    outfile << "ndata      = " << num_data << endl;
    outfile << endl;
    
    double num_bkg(0);
    double num_mc_total[9], n_mc_total[9];
    for (int i=0; i<9; i++){
      if(i<6){
	num_bkg += num_mc_2010[i];
	num_bkg += num_mc_2011[i];
      }
      num_mc_total[i] = num_mc_2010[i] + num_mc_2011[i];
      n_mc_total[i] = n_mc_2010[i] * scale_mc_2010[i] + n_mc_2011[i] * scale_mc_2011[i];
    }
    
    outfile << "nallbkg    = " << num_bkg << endl;
    outfile << "nhig glu   = " << num_mc_total[6] << endl;
    outfile << "nhig vbf   = " << num_mc_total[7] << endl;
    outfile << "nhig wzt   = " << num_mc_total[8] << endl;
    outfile << "ndy        = " << num_mc_total[5] << endl;
    outfile << "nbox       = " << num_mc_total[0] << endl;
    outfile << "ndiphot    = " << num_mc_total[1] << endl;
    outfile << "ngjet      = " << num_mc_total[2] << endl;
    outfile << "nqcd40     = " << num_mc_total[3] << endl;
    outfile << "nqcd30-40  = " << num_mc_total[4] << endl;
    outfile << endl;
    outfile << "eff nhig      = " << (num_mc_total[6] + num_mc_total[7] + num_mc_total[8]) 
      / (n_mc_total[6] + n_mc_total[7] + n_mc_total[8]) << endl;
    outfile << "eff nhig glu  = " << num_mc_total[6] / n_mc_total[6] << endl;
    outfile << "eff nhig vbf  = " << num_mc_total[7] / n_mc_total[7] << endl;
    outfile << "eff nhig wzt  = " << num_mc_total[8] / n_mc_total[8] << endl;
    outfile << "eff ndy       = " << num_mc_total[5] / n_mc_total[5] << endl;
    outfile << "eff nbox      = " << num_mc_total[0] / n_mc_total[0] << endl;
    outfile << "eff ndiphot   = " << num_mc_total[1] / n_mc_total[1] << endl;
    outfile << "eff ngjet     = " << num_mc_total[2] / n_mc_total[2] << endl;
    outfile << "eff nqcd      = " << num_mc_total[3] / n_mc_total[3] << endl;
    outfile << "eff nqcd30-40 = " << num_mc_total[4] / n_mc_total[4] << endl;
    outfile << endl;
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

  vector<double> values;

  if (variable == "massgg"){
    values.push_back(integralhiggs);
    values.push_back(integralbkg);
    values.push_back(integralhiggs/sqrt(integralbkg));
  }

  return values;

}








