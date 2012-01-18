#include <TTree.h>
#include <TFile.h>
#include <TLegend.h>
#include <TF1.h>
#include <TH2.h>
#include <TStyle.h>
#include <TROOT.h>
#include <TCanvas.h>
#include <TMath.h>
#include <iostream>
#include <fstream>
#include <vector>
#include <TLorentzVector.h>
#include <fillHisto.h>


// USE:
//
// .x optimize(data int lumi 2010, data int lumi 2011, pt1 cut, pt2 cut, ptj1 cut, ptj2 cut, deltae cut, zep cut, mjj cut, eb cat, r9 cat, cic selection)
// 
// example:
//
// .x optimize.C(33,230,40,25,20,15,2.5,2.5,300,1,1,4,1) 

double probability_disc(double obs, double exp){
 double prob=0.;
 TF1 *mypoisson = new TF1("mypoisson","TMath::Poisson(x,[0])",0,10000);
 mypoisson->SetParameter(0, exp);
 if (exp<100) 
   prob= mypoisson->Integral(obs-0.5,999); 
 else 
   prob= mypoisson->Integral(obs-0.5,9999);//  for (int i=0;i <  obs; ++i)
 // std::cout << "p-value to observe " << obs  << " with exp " << exp << ": " << 1-prob  << std::endl;
 // std::cout << "loc. significance " << TMath::NormQuantile(1-(1-prob)) << std::endl;
 return TMath::NormQuantile(1-prob);
} 

double UL(double expsig, double expbkg){
  double prob=0.;
  TF1 *mypoisson = new TF1("mypoisson","TMath::Poisson(x,[0])",0,10000);
  double ul(0);
  for(int i=0; i<100000; i++){
    double exp = i*(0.0002*expsig) + expbkg;
    mypoisson->SetParameter(0, exp);
    if (exp<100) 
      prob= mypoisson->Integral(0,expbkg-0.5); 
    else 
      prob= mypoisson->Integral(0,expbkg-0.5);
    //    std::cout << "p-value to observe " << expbkg  << " with exp " << exp << ": " << prob  << std::endl;    
    if(prob<0.05){
      ul = exp-expbkg;
      break;
    }
  }
  return ul/expsig;
} 


vector <double> optimize(double int_exp, double pt1=50, double pt2=30, double pthiggsmin = -100, double pthiggsmax = -100, double ptj1=20, double ptj2=15, double deltae=2.5, double zep=2.5, double mjj=300, double deltap = 2.6, int eb = 1, int r9 = 1, int cic = 4, bool thirdcat = 0, string variable = "massgg", int nbin = 200, double min = 90, double max = 190, string axis = "m(#gamma#gamma)[GeV]"){


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

  TFile* mc[9];

  TFile* mc_gluglu[7];
  TFile* mc_vbf[7];
  TFile* mc_wzh[7];
  TFile* mc_tth[7];

  int h_masses[7] = {100,105,110,115,120,130,140};

  TString redntpDir= "root://pccmsrm23.cern.ch:1094//u2/xrootd/delre/Higgs/reduced/";
  //TString redntpDir= "root://pccmsrm23.cern.ch:1094//u2/xrootd/meridian/Higgs/reduced_bck/";
  //TString redntpDir= "/Users/delre/";
  TString preselectionLevel;

  preselectionLevel="cicloose";
  TString preselectionLevelCS="preselectionCS";

  // total data sample
  TFile* data = TFile::Open(redntpDir+"/redntp.42xv4_data_new."+preselectionLevel+".regrPho_eCorr.v3/merged/redntp_Photon-Run2011A-03Oct2011-05JulReReco-05AugReReco-Prompt-v1-DiPhotonSkimOnFly-b.root");
  TFile* datacs = TFile::Open(redntpDir+"/redntp.42xv4_data_new.preselectionCS.regrPho_eCorr.v3/merged/redntp_Photon-Run2011A-03Oct2011-05JulReReco-05AugReReco-Prompt-v1-DiPhotonSkimOnFly-b.root");

  if(int_exp>0){
    // box samples
    mc[0] = TFile::Open(redntpDir+"/redntp.42xv4."+preselectionLevel+".regrPho_eCorr.v3/merged/redntp_DiPhotonBox_Pt-25To250_7TeV-pythia6-Fall11-PU_S6_START42_V14B-v1.root");
    // diphoton jets samples
    mc[1] = TFile::Open(redntpDir+"/redntp.42xv4."+preselectionLevel+".regrPho_eCorr.v3/merged/redntp_DiPhotonJets_7TeV-madgraph-Fall11-PU_S6_START42_V14B-v1.root");
    // gjet samples
    mc[2] = TFile::Open(redntpDir+"/redntp.42xv4."+preselectionLevel+".regrPho_eCorr.v3/merged/redntp_GJet_Pt-20_doubleEMEnriched_TuneZ2_7TeV-pythia6-Fall11-PU_S6_START42_V14B-v1.root");
    // qcd pt>40 samples
    mc[3] = TFile::Open(redntpDir+"/redntp.42xv4."+preselectionLevel+".regrPho_eCorr.v3/merged/redntp_QCD_Pt-40_doubleEMEnriched_TuneZ2_7TeV-pythia6-Fall11-PU_S6_START42_V14B-v1.root");
   // qcd 30<pt<40 samples
    mc[4] = TFile::Open(redntpDir+"/redntp.42xv4."+preselectionLevel+".regrPho_eCorr.v3/merged/redntp_QCD_Pt-30to40_doubleEMEnriched_TuneZ2_7TeV-pythia6-Fall11-PU_S6_START42_V14B-v1.root");
    // drell yan samples
    mc[5] = TFile::Open(redntpDir+"/redntp.42xv4."+preselectionLevel+".regrPho_eCorr.v3/merged/redntp_DYJetsToLL_TuneZ2_M-50_7TeV-madgraph-tauola-Fall11-PU_S6_START42_V14B-v1.root");
    // gluglu higgs samples 
    mc[6] = TFile::Open(redntpDir+"/redntp.42xv4."+preselectionLevel+".regrPho_eCorr.v3/merged/redntp_GluGluToHToGG_M-120_7TeV-powheg-pythia6-Fall11-PU_S6_START42_V14B-v1.root");
    // vbf higgs samples
    mc[7] = TFile::Open(redntpDir+"/redntp.42xv4."+preselectionLevel+".regrPho_eCorr.v3/merged/redntp_VBF_HToGG_M-120_7TeV-powheg-pythia6-Fall11-PU_S6_START42_V14B-v1.root");
  // W/Z/TT H higgs samples 
    mc[8] = TFile::Open(redntpDir+"/redntp.42xv4."+preselectionLevel+".regrPho_eCorr.v3/merged/redntp_WH_ZH_HToGG_M-120_7TeV-pythia6-Fall11-PU_S6_START42_V14B-v1.root");
  }

  for (int i=0;i<7;i++){
    char hmass[100]; sprintf(hmass,"%d",h_masses[i]);
    // gluglu samples
    mc_gluglu[i] = TFile::Open(redntpDir+"/redntp.42xv4."+preselectionLevel+".regrPho_eCorr.v3/merged/redntp_GluGluToHToGG_M-"+ hmass + "_7TeV-powheg-pythia6-Fall11-PU_S6_START42_V14B-v1.root");
    // vbf higgs samples 
    mc_vbf[i] = TFile::Open(redntpDir+"/redntp.42xv4."+preselectionLevel+".regrPho_eCorr.v3/merged/redntp_VBF_HToGG_M-"+ hmass +"_7TeV-powheg-pythia6-Fall11-PU_S6_START42_V14B-v1.root");
    // W/Z H higgs samples 
    mc_wzh[i] = TFile::Open(redntpDir+"/redntp.42xv4."+preselectionLevel+".regrPho_eCorr.v3/merged/redntp_WH_ZH_HToGG_M-"+ hmass +"_7TeV-pythia6-Fall11-PU_S6_START42_V14B-v1.root");
    // TT H higgs samples 
    //  mc_tth[i] = TFile::Open(redntpDir+"/redntp.42xv5b."+preselectionLevel+"/merge/redntp_TTH_HToGG_M-"+ hmass +"_7TeV-pythia6.root");
  }

  // k factors  
  double kfactordiphot = 1.3;
  double kfactordiphotmadgraph = 1.15;
  double kfactorgamjet = 1.3;
  double kfactorqcd = 1;
  double kfactordy = 1.15;

  // cross sections and scaling
  double boosthiggs(1);
  double cross_mc[9];
  cross_mc[0] = 12.37 * kfactordiphot; // box
  cross_mc[1] = 134 * kfactordiphotmadgraph; // diphoton jets
  cross_mc[2] = 493.44 * kfactorgamjet; // gjet
  cross_mc[3] = 40392 * kfactorqcd; // qcd pt>40
  cross_mc[4] = 9610 * kfactorqcd; // qcd 30<pt<40 
  cross_mc[5] = 2321 * kfactordy; // drell yan
//   cross_mc[6] = 16.63 * 2.13e-03 * boosthiggs; // glu glu higgs SM 120
//   cross_mc[7] = 1.269 * 2.13e-03 * boosthiggs; // vbf higgs SM 120
//   cross_mc[8] = (0.6561 + 0.3598) * 2.13e-03 * boosthiggs; // WHtt higgs SM 120
  cross_mc[6] = 0; // glu glu FP higgs
  cross_mc[7] = 1.269 * 0.0231 * boosthiggs; // vbf FP higgs
  cross_mc[8] = 1.016 * 0.0231 * boosthiggs; // WHtt FP higgs
 
  // getting the number of original events in each sample (processed with CMSSW)
  int n_mc[9],n_gluglu[7],n_vbf[7],n_wzh[7],n_tth[7];
  for(int i=0; i<9; i++){
    n_mc[i] = 0;
    if(int_exp>0) n_mc[i] = ((TH1D*)mc[i]->Get("ptphotgen1"))->GetEntries();
  }
  for(int i=0; i<7; i++){
    n_gluglu[i] = n_vbf[i] = n_wzh[i] = n_tth[i] = 0;
    n_gluglu[i] = ((TH1D*)mc_gluglu[i]->Get("ptphotgen1"))->GetEntries();
    n_vbf[i] = ((TH1D*)mc_vbf[i]->Get("ptphotgen1"))->GetEntries();
    n_wzh[i] = ((TH1D*)mc_wzh[i]->Get("ptphotgen1"))->GetEntries();
//     n_tth[i] = ((TH1D*)mc_tth[i]->Get("ptphotgen1"))->GetEntries();
  }

  // setting the scaling factor to actual lumi
  double scale_mc[9];
  for(int i=0; i<9; i++){
    scale_mc[i] = 0; 
    if(int_exp>0) scale_mc[i] = cross_mc[i] * int_exp / n_mc[i];
  }

  // char for output name
  char name[1000];
  char allcut[3000];
  sprintf(allcut,"%3.1f%s%3.1f%s%3.1f%s%3.1f%s%3.1f%s%3.1f%s%3.1f%s%3.1f%s%3.1f%s%3.1f%s%d%s%d%s%d%s%d",pt1,"_",pt2,"_",pthiggsmin,"_",pthiggsmax,"_",ptj1,"_",ptj2,"_",deltae,"_",zep,"_",mjj,"_",deltap,"_",eb,"_",r9,"_",thirdcat,"_",cic);

  // output root file
  sprintf(name,"%s%s%s%s%s","results_gg/histo_",variable.c_str(),"_",allcut,".root");
  TFile * hOutputFile   = new TFile(name, "RECREATE" ) ;

  // histograms needed by the machinery
  TH1D* vardata = new TH1D("vardata","vardata",nbin,min,max);
  TH1D* vardatacs = new TH1D("vardatacs","vardatacs",nbin,min,max);
  TH1D* var_mc[9];
  TH1D* var_gluglu[7];
  TH1D* var_vbf[7];
  TH1D* var_wzh[7];
  TH1D* var_tth[7];
  TH1D * var[6];
  for (int i=0; i<6; i++) {
    sprintf(name,"%s%d","var",i);
    var[i] = new TH1D(name,name,nbin,min,max);
  }
  for (int i=0; i<9; i++) {
    sprintf(name,"%s%d","var_mc_",i);
    var_mc[i] = new TH1D(name,name,nbin,min,max);
  }
  for (int i=0; i<7; i++) {
    sprintf(name,"%s%d","var_gluglu_",h_masses[i]);
    var_gluglu[i] = new TH1D(name,name,nbin,min,max);
    sprintf(name,"%s%d","var_vbf_",h_masses[i]);
    var_vbf[i] = new TH1D(name,name,nbin,min,max);
    sprintf(name,"%s%d","var_wzh_",h_masses[i]);
    var_wzh[i] = new TH1D(name,name,nbin,min,max);
//     sprintf(name,"%s%d","var_tth_",h_masses[i]);
//     var_tth[i] = new TH1D(name,name,nbin,min,max);
  }

  // creating the fillers and setting cuts
  fillHisto data_fill((TTree*)data->Get("AnaTree"), 1);
  fillHisto datacs_fill((TTree*)datacs->Get("AnaTree"), 1);
  data_fill.Setcuts(pt1,pt2,pthiggsmin,pthiggsmax,ptj1,ptj2,deltae,zep,mjj,deltap,eb,r9,thirdcat);
  datacs_fill.Setcuts(pt1,pt2,pthiggsmin,pthiggsmax,ptj1,ptj2,deltae,zep,mjj,deltap,eb,r9,thirdcat);
  if(cic>0)
    {
      data_fill.setCic(cic);
      datacs_fill.setCic(cic);
    }

  fillHisto* mc_fill[9];
  fillHisto* mc_gluglu_fill[7];  
  fillHisto* mc_vbf_fill[7];  
  fillHisto* mc_wzh_fill[7];  
  fillHisto* mc_tth_fill[7];  

  for (int i=0; i<9; i++){
    if(int_exp>0) mc_fill[i] = new fillHisto((TTree*)mc[i]->Get("AnaTree"), 1);
    if(int_exp>0) mc_fill[i]->Setcuts(pt1,pt2,pthiggsmin,pthiggsmax,ptj1,ptj2,deltae,zep,mjj,deltap,eb,r9,thirdcat);

    mc_fill[i]->DoPuReweight();
    mc_fill[i]->DoPtReweight();
    
    if(cic>0){
      if(int_exp>0) mc_fill[i]->setCic(cic);
    }
  }
 
  for (int i=0; i<7; i++){
    mc_gluglu_fill[i] = new fillHisto((TTree*)mc_gluglu[i]->Get("AnaTree"), 1);
    mc_vbf_fill[i] = new fillHisto((TTree*)mc_vbf[i]->Get("AnaTree"), 1);
    mc_wzh_fill[i] = new fillHisto((TTree*)mc_wzh[i]->Get("AnaTree"), 1);
//     mc_tth_fill[i] = new fillHisto((TTree*)mc_tth[i]->Get("AnaTree"), 1);
    mc_gluglu_fill[i]->Setcuts(pt1,pt2,pthiggsmin,pthiggsmax,ptj1,ptj2,deltae,zep,mjj,deltap,eb,r9,thirdcat);
    mc_vbf_fill[i]->Setcuts(pt1,pt2,pthiggsmin,pthiggsmax,ptj1,ptj2,deltae,zep,mjj,deltap,eb,r9,thirdcat);
    mc_wzh_fill[i]->Setcuts(pt1,pt2,pthiggsmin,pthiggsmax,ptj1,ptj2,deltae,zep,mjj,deltap,eb,r9,thirdcat);
//     mc_tth_fill[i]->Setcuts(pt1,pt2,pthiggsmin,pthiggsmax,ptj1,ptj2,deltae,zep,mjj,deltap,eb,r9,thirdcat);
    mc_gluglu_fill[i]->DoPuReweight();
    mc_vbf_fill[i]->DoPuReweight();
    mc_wzh_fill[i]->DoPuReweight();
//     mc_tth_fill[i]->DoPuReweight();
    mc_gluglu_fill[i]->DoPtReweight();
    mc_vbf_fill[i]->DoPtReweight();
    mc_wzh_fill[i]->DoPtReweight();
//     mc_tth_fill[i]->DoPtReweight();
    mc_gluglu_fill[i]->setCic(cic);
    mc_vbf_fill[i]->setCic(cic);
    mc_wzh_fill[i]->setCic(cic);
//     mc_tth_fill[i]->setCic(cic);
  }

  TFile *datainput;
  TFile *datacsinput;
  TFile* mcinput[9];
  TFile* mc_glugluinput[7];
  TFile* mc_vbfinput[7];
  TFile* mc_wzhinput[7];
  TFile* mc_tthinput[7];

  // filling histograms
  std::cout << " ++++++++++++++ DATA ++++++++++++++++" << std::endl;
  cout << "running over " << ((TTree*)data->Get("AnaTree"))->GetEntries("") << " data events" <<  endl;
  if (variable == "massgg") {  
    sprintf(name,"%s%s%s","results_gg/events_",allcut,".txt");
    data_fill.Writetxt(name);
    sprintf(name,"%s%s%s","results_gg/events_",allcut,".root");
    //    data_fill.WriteRoot(name);
  }
  datainput = data_fill.File("results_gg/data.root");
  vardata->Add((TH1D*)datainput->Get("massgg")); 
  std::cout << "Selected events on data " << vardata->GetEntries() << std::endl;
  cout << "running over " << ((TTree*)datacs->Get("AnaTree"))->GetEntries("") << " data events (for cs)" <<  endl; 

  if (variable == "massgg") {  
    sprintf(name,"%s%s%s","results_gg/events_",allcut,"_cs.txt");
    datacs_fill.Writetxt(name);
    sprintf(name,"%s%s%s","results_gg/events_",allcut,"_cs.root");
    //    datacs_fill.WriteRoot(name);
  }
  datacsinput = datacs_fill.File("results_gg/datacs.root",1);
  vardatacs->Add((TH1D*)datacsinput->Get("massgg")); 
  std::cout << "Selected events on data cs " << vardatacs->GetEntries() << std::endl;

  std::cout << " ++++++++++++++ MC ++++++++++++++++" << std::endl;
  for (int i=0; i<9; i++){ 
    sprintf(name,"%s%s",mcnames[i].c_str()," 2011");
    if(int_exp>0) {
      cout << "running over " << ((TTree*)mc[i]->Get("AnaTree"))->GetEntries("") << " " << name << " events" <<  endl; 
      sprintf(name,"%s%s%s","results_gg/events_",mcnames[i].c_str(),".root");
      if (variable == "massgg") {
	//	mc_fill[i]->WriteRoot(name);
	// 	sprintf(name,"%s%s%s%s%s","results_gg/events_",mcnames[i].c_str(),"_",allcut,".txt");
	// 	if(i==7) mc_fill[i]->Writetxt(name);
      }
      mcinput[i] = mc_fill[i]->File(name);
      var_mc[i]->Add((TH1D*)mcinput[i]->Get("massgg"));
      std::cout << "Selected events on mc2011 " << name << " " << var_mc[i]->GetEntries() << std::endl;
    }

  }

  std::cout << " ++++++++++++++ signal MC ++++++++++++++++" << std::endl;
  for (int i=0; i<7; i++){ 
    cout << "running over " << ((TTree*)mc_gluglu[i]->Get("AnaTree"))->GetEntries("") << " gluglu M=" << h_masses[i] << " events" <<  endl; 
    sprintf(name,"%s%d%s","results_gg/events_gluglu",h_masses[i],".root");
    //    if (variable == "massgg") mc_gluglu_fill[i]->WriteRoot(name);
    mc_glugluinput[i] = mc_gluglu_fill[i]->File(name);
    var_gluglu[i]->Add((TH1D*)mc_glugluinput[i]->Get("massgg"));
    std::cout << "Selected events on mc2011 gluglu " << h_masses[i] << " " << var_gluglu[i]->GetEntries() << std::endl;
 
    cout << "running over " << ((TTree*)mc_vbf[i]->Get("AnaTree"))->GetEntries("") << " vbf M=" << h_masses[i] << " events" <<  endl; 
    sprintf(name,"%s%d%s","results_gg/events_vbf",h_masses[i],".root");
    //    if (variable == "massgg") mc_vbf_fill[i]->WriteRoot(name);
    mc_vbfinput[i] = mc_vbf_fill[i]->File(name);
    var_vbf[i]->Add((TH1D*)mc_vbfinput[i]->Get("massgg"));
    std::cout << "Selected events on mc2011 vbf " << h_masses[i] << " " << var_vbf[i]->GetEntries() << std::endl;

    cout << "running over " << ((TTree*)mc_wzh[i]->Get("AnaTree"))->GetEntries("") << " wzh M=" << h_masses[i] << " events" <<  endl; 
    sprintf(name,"%s%d%s","results_gg/events_wzh",h_masses[i],".root");
    //    if (variable == "massgg") mc_wzh_fill[i]->WriteRoot(name);
    mc_wzhinput[i] = mc_wzh_fill[i]->File(name);
    var_wzh[i]->Add((TH1D*)mc_wzhinput[i]->Get("massgg"));
    std::cout << "Selected events on mc2011 wzh " << h_masses[i] << " " << var_wzh[i]->GetEntries() << std::endl;

//     cout << "running over " << ((TTree*)mc_tth[i]->Get("AnaTree"))->GetEntries("") << " tth M=" << h_masses[i] << " events" <<  endl; 
//     sprintf(name,"%s%d%s%s%s","results_gg/events_tth",h_masses[i],"_",allcut,".root");
//     if (variable == "massgg") mc_tth_fill[i]->WriteRoot(name);
//     var_tth[i]->Add(mc_tth_fill[i]->Plot(variable, name, nbin, min, max));
//     std::cout << "Selected events on mc2011 tth " << h_masses[i] << " " << var_tth[i]->GetEntries() << std::endl;
 }

  // scale mc to equivalent lumi
  for (int i=0; i<9; i++){ 
    if(int_exp>0) var_mc[i]->Scale(scale_mc[i]);  
  }

  // counting number of events passing selection (scaled)
  double num_mc[9],num_gluglu[7],num_vbf[7],num_wzh[7],num_tth[7]; 
  double num_uns_mc[9],num_uns_gluglu[7],num_uns_vbf[7],num_uns_wzh[7],num_uns_tth[7]; 

  for (int i=0; i<9; i++){ 
    num_mc[i] = 0;
    if(int_exp>0) num_mc[i] = var_mc[i]->Integral();  
    num_uns_mc[i] = 0;
    if(int_exp>0) num_uns_mc[i] = var_mc[i]->GetEntries();  
  }
  for(int i=0; i<7; i++){
    num_gluglu[i] = num_vbf[i] = num_wzh[i] = num_tth[i] = 0;
    num_gluglu[i] = num_gluglu[i] = var_gluglu[i]->Integral(); 
    num_vbf[i] = num_vbf[i] = var_vbf[i]->Integral(); 
    num_wzh[i] = num_wzh[i] = var_wzh[i]->Integral(); 
    num_uns_gluglu[i] = num_uns_vbf[i] = num_uns_wzh[i] = num_uns_tth[i] = 0;
    num_uns_gluglu[i] = num_uns_gluglu[i] = var_gluglu[i]->GetEntries(); 
    num_uns_vbf[i] = num_uns_vbf[i] = var_vbf[i]->GetEntries(); 
    num_uns_wzh[i] = num_uns_wzh[i] = var_wzh[i]->GetEntries(); 
  }

  // add two QCD bins
  if(int_exp>0) var_mc[3]->Add(var_mc[4]);
  
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
	if(int_exp>0) var[j]->SetBinContent(i,var_mc[k]->GetBinContent(i) + var[j]->GetBinContent(i));
      }	
    }    
  }
 
  //final plots
  char ytitle[100];
  sprintf(ytitle,"%s%d%s","N_{ev}/",int(int_exp),"pb^{-1}");
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
  leg = new TLegend(0.6,0.6,0.85,0.85);
  leg->SetFillStyle(0); leg->SetBorderSize(0); leg->SetTextSize(0.05);
  leg->SetFillColor(0);
  sprintf(name,"H M(120)");
  if(boosthiggs!=1) sprintf(name,"%s%2.0f", "H M(120) x ", boosthiggs);
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
  sprintf(name,"%s%s%s%s%s","results_gg/mc_",variable.c_str(),"_",allcut,".png");
  c0->SaveAs(name);

  //higgs only plot
  var_mc[6]->Add(var_mc[7]);
  var_mc[6]->Add(var_mc[8]);
  sprintf(ytitle,"%s%d%s","N_{ev}/",int(int_exp),"pb^{-1}");
  var_mc[6]->SetXTitle(axis.c_str());
  var_mc[6]->SetYTitle(ytitle);
  var_mc[6]->SetTitle("");
  var_mc[6]->SetLineColor(kBlack);
  var_mc[6]->SetLineWidth(2);
  var_mc[6]->SetFillColor(kYellow);
  var_mc[6]->Draw();
  sprintf(name,"%s%s%s%s%s","results_gg/higgs_",variable.c_str(),"_",allcut,".png");
  c0->SaveAs(name);
  
  // data only plot
  vardata->SetXTitle(axis.c_str());
  vardata->SetTitle("");
  vardata->SetStats(0);
  vardata->SetMarkerStyle(8);
  vardata->SetMarkerSize(.9);
  vardata->SetTitleOffset(1.25,"Y");
  vardata->Draw("pe");
  sprintf(name,"%s%s%s%s%s","results_gg/data_",variable.c_str(),"_",allcut,".png");
  c0->SaveAs(name);

  // data overlaid to mc

  legge = leg->AddEntry(vardata, "data", "p");

  double themax =   vardata->GetMaximum();
  if(var[0]->GetMaximum()>themax) themax = var[0]->GetMaximum();
  if (
      variable == "etaphot1" || variable == "etaphot2" ||
      variable == "phiphot1" || variable == "phiphot2" ||
      variable == "etajet1" || variable == "etajet2" ||
      variable == "phijet1" || variable == "phijet2"
      )
    vardata->SetMaximum(themax*2.0);
  else if ( variable == "met" ) {
    c0->SetLogy(1);
    vardata->SetMaximum(themax*5.0);   
  }
  else
    vardata->SetMaximum(themax*1.1);
  vardata->SetMinimum(0.001);
  var[0]->Draw("same");
  var[1]->Draw("same");
  var[2]->Draw("same");
  var[3]->Draw("same");
  var[4]->Draw("same");
  var[5]->Draw("same");
  leg->Draw();
  vardata->Draw("pesame");
  gPad->RedrawAxis();
  sprintf(name,"%s%s%s%s%s","results_gg/data-mc_",variable.c_str(),"_",allcut,".png");
  c0->SaveAs(name);

  //data with control sample
  vardata->Draw("pe");
  vardatacs->SetLineColor(46);
  vardatacs->SetFillColor(42);
  vardatacs->SetLineWidth(3);
  vardatacs->Draw("hsame");
  vardata->Draw("pesame");
  sprintf(name,"%s%s%s%s%s","results_gg/datacs_",variable.c_str(),"_",allcut,".png");
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

  double integralhiggs ;
  double entrieshiggs ;
  double integralbkg ;
  // some calculation for s/sqrt(b)
  TF1 *mypowlaw = new TF1("mypowlaw","[0]*pow(x,[1])",100,180);
  TF1 *mydoublegau = new TF1("mydoublegau","gaus(0)+gaus(3)",100,180);
  double binsize = (max-min)/nbin;
  if (variable == "massgg")
    {  
      sprintf(name,"%s%s%s%d%s","results_gg/optimalcut_",variable.c_str(),"_",eb,".txt");
      mydoublegau->SetParameter(0, 1.);
      mydoublegau->SetParameter(2, 1.5);
      mydoublegau->SetParameter(3, 1.);
      mydoublegau->SetParameter(5, 4.);      
      mydoublegau->FixParameter(1,120.);
      mydoublegau->FixParameter(4,120.);
      mydoublegau->SetParLimits(2,0.5,2.);
      mydoublegau->SetParLimits(5,2.,6.);    
      mypowlaw->SetParameter(0,2.88660e+07);
      mypowlaw->SetParameter(1,-3.42928e+00);
      var_mc[6]->Fit("mydoublegau","l");
      var[1]->Fit("mypowlaw","l");
      integralhiggs = mydoublegau->Integral(117.,123.)/binsize;
      integralbkg = mypowlaw->Integral(117.,123.)/binsize;
      cout << "Number of signal events " << integralhiggs/boosthiggs << endl;
      cout << "Number of bkg events " << integralbkg << endl;
      cout << "S/sqrt(B) " << integralhiggs/boosthiggs/sqrt(integralbkg) << endl;
      cout << "probability "<< probability_disc(integralhiggs/boosthiggs+integralbkg,integralbkg) << endl;
  }  

  var_mc[6]->Draw();
  sprintf(name,"%s%s%s%s%s","results_gg/higgs_fit_",variable.c_str(),"_",allcut,".png");
  c0->SaveAs(name);
 
  var[1]->Draw();
  sprintf(name,"%s%s%s%s%s","results_gg/mc_fit_",variable.c_str(),"_",allcut,".png");
  c0->SaveAs(name);
  


  if (variable == "massgg"){

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
    outfile << "deltaphicut : " << deltap << endl;
    outfile << "zeppencut : " << zep << endl;
    outfile << "invmassjetcut : " << mjj << endl;
    outfile << "CiC level : " << cic << endl;
    outfile << "ebcat : " << eb << endl;
    outfile << "r9cat : " << r9 << endl;
    outfile << "thirdcat : " << thirdcat << endl;
    outfile << endl;
    outfile << "####################################" << endl;
    outfile << "N of generated events" << endl;
    outfile << "####################################" << endl;
    outfile << "# events hig =       " << n_mc[6] << endl;
    outfile << "# events hig_vbf2011 =    " << n_mc[7] << endl;
    outfile << "# events hig_wzt2011 =    " << n_mc[8] << endl;
    outfile << "# events dy =        " << n_mc[5] << endl;
    outfile << "# events box =       " << n_mc[0] << endl;
    outfile << "# events diphotjet = " << n_mc[1] << endl;
    outfile << "# events gjet =      " << n_mc[2] << endl;
    outfile << "# events qcd =       " << n_mc[3] << endl;
    outfile << "# events qcd2 =      " << n_mc[4] << endl;
    outfile << endl;
    outfile << "####################################" << endl;
    outfile << "N of selected events and eff." << endl;
    outfile << "####################################" << endl; 
    outfile << "ndata      = " << num_data << endl;
    outfile << endl;
    
    double num_bkg(0), err_num_bkg(0);
    double num_mc_total[9],num_uns_mc_total[9], n_mc_total[9];
    double err_num_mc_total[9],err_num_uns_mc_total[9];
    for (int i=0; i<9; i++){
      if(i<6){
	num_bkg += num_mc[i];
      }
      num_mc_total[i] = num_mc[i];
      num_uns_mc_total[i] = num_uns_mc[i];
      err_num_uns_mc_total[i] = sqrt(num_uns_mc_total[i]);
      if(num_uns_mc_total[i]) err_num_mc_total[i] = err_num_uns_mc_total[i] * num_mc_total[i]/num_uns_mc_total[i];
      else err_num_mc_total[i] = 0;
      n_mc_total[i] = n_mc[i] * scale_mc[i];
    }
    for (int i=0; i<6; i++){
      err_num_bkg = sqrt(err_num_bkg*err_num_bkg + err_num_mc_total[i]*err_num_mc_total[i]);
    }

    outfile << "nallbkg    = " << num_bkg << " +/- " << err_num_bkg << endl;
    outfile << "nhig glu   = " << num_mc_total[6] << " +/- " << err_num_mc_total[6] << endl;
    outfile << "nhig vbf   = " << num_mc_total[7] << " +/- " << err_num_mc_total[7] << endl;
    outfile << "nhig wzt   = " << num_mc_total[8] << " +/- " << err_num_mc_total[8] << endl;
    outfile << "ndy        = " << num_mc_total[5] << " +/- " << err_num_mc_total[5] << endl;
    outfile << "nbox       = " << num_mc_total[0] << " +/- " << err_num_mc_total[0] << endl;
    outfile << "ndiphot    = " << num_mc_total[1] << " +/- " << err_num_mc_total[1] << endl;
    outfile << "ngjet      = " << num_mc_total[2] << " +/- " << err_num_mc_total[2] << endl;
    outfile << "nqcd40     = " << num_mc_total[3] << " +/- " << err_num_mc_total[3] << endl;
    outfile << "nqcd30-40  = " << num_mc_total[4] << " +/- " << err_num_mc_total[4] << endl;
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
      outfile << "eff nhig gluglu "<< h_masses[i] << " = " << num_gluglu[i] / n_gluglu[i] << endl;
      outfile << "eff nhig vbf    "<< h_masses[i] << " = " << num_vbf[i] / n_vbf[i] << endl;
      outfile << "eff nhig wzh    "<< h_masses[i] << " = " << num_wzh[i] / n_wzh[i] << endl;
//       outfile << "eff nhig tth    "<< h_masses[i] << " = " << num_tth[i] / n_tth[i] << endl;
    }
    outfile.close();

  }

  // CODE FOR OPTIMIZATION
  // ----------------------------
  hOutputFile->cd();
  // VBF
  // double ptgg_range(80), phot1_range(60), phot2_range(30), jet1_range(40), jet2_range(20), deltae_range(2), zep_range(2), mjj_range(800), deltap_range(1.0), met_range(140.);
  // double ptgg_central(40), phot1_central(65), phot2_central(35), jet1_central(40), jet2_central(30), deltae_central(3.5), zep_central(2.5), mjj_central(450), deltap_central(2.6), met_central(70);
  // WZH
  double ptgg_range(100), phot1_range(60), phot2_range(30), jet1_range(40), jet2_range(30), deltae_range(2), zep_range(3), mjj_range(60), deltap_range(1.0),  met_range(140);
  double ptgg_central(50), phot1_central(65), phot2_central(35), jet1_central(40), jet2_central(35), deltae_central(-2.5), zep_central(2.), mjj_central(-35), deltap_central(2.6), met_central(70);
  
  TH1D ptgg_prob_disc("ptgg_prob_disc","",20,ptgg_central-ptgg_range/2.,ptgg_central+ptgg_range/2.);                 
  TH1D phot1_prob_disc("phot1_prob_disc","",20,phot1_central-phot1_range/2.,phot1_central+phot1_range/2.);
  TH1D phot2_prob_disc("phot2_prob_disc","",20,phot2_central-phot2_range/2.,phot2_central+phot2_range/2.);
  TH1D jet1_prob_disc("jet1_prob_disc","",20,jet1_central-jet1_range/2.,jet1_central+jet1_range/2.);
  TH1D jet2_prob_disc("jet2_prob_disc","",20,jet2_central-jet2_range/2.,jet2_central+jet2_range/2.);
  TH1D deltae_prob_disc("deltae_prob_disc","",20,deltae_central-deltae_range/2.,deltae_central+deltae_range/2.);
  TH1D zep_prob_disc("zep_prob_disc","",20,zep_central-zep_range/2.,zep_central+zep_range/2.);
  TH1D mjj_prob_disc("mjj_prob_disc","",20,mjj_central-mjj_range/2.,mjj_central+mjj_range/2.);
  TH1D deltap_prob_disc("deltap_prob_disc","",20,deltap_central-deltap_range/2.,deltap_central+deltap_range/2.);  
  TH1D met_prob_disc("met_prob_disc","",20,met_central-met_range/2.,met_central+met_range/2.);  

  TH1D ptgg_prob_UL("ptgg_prob_UL","",20,ptgg_central-ptgg_range/2.,ptgg_central+ptgg_range/2.);                 
  TH1D phot1_prob_UL("phot1_prob_UL","",20,phot1_central-phot1_range/2.,phot1_central+phot1_range/2.);
  TH1D phot2_prob_UL("phot2_prob_UL","",20,phot2_central-phot2_range/2.,phot2_central+phot2_range/2.);
  TH1D jet1_prob_UL("jet1_prob_UL","",20,jet1_central-jet1_range/2.,jet1_central+jet1_range/2.);
  TH1D jet2_prob_UL("jet2_prob_UL","",20,jet2_central-jet2_range/2.,jet2_central+jet2_range/2.);
  TH1D deltae_prob_UL("deltae_prob_UL","",20,deltae_central-deltae_range/2.,deltae_central+deltae_range/2.);
  TH1D zep_prob_UL("zep_prob_UL","",20,zep_central-zep_range/2.,zep_central+zep_range/2.);
  TH1D mjj_prob_UL("mjj_prob_UL","",20,mjj_central-mjj_range/2.,mjj_central+mjj_range/2.);
  TH1D deltap_prob_UL("deltap_prob_UL","",20,deltap_central-deltap_range/2.,deltap_central+deltap_range/2.);  
  TH1D met_prob_UL("met_prob_UL","",20,met_central-met_range/2.,met_central+met_range/2.);  

  TH1D ptgg_S("ptgg_S","",20,ptgg_central-ptgg_range/2.,ptgg_central+ptgg_range/2.);                 
  TH1D phot1_S("phot1_S","",20,phot1_central-phot1_range/2.,phot1_central+phot1_range/2.);
  TH1D phot2_S("phot2_S","",20,phot2_central-phot2_range/2.,phot2_central+phot2_range/2.);
  TH1D jet1_S("jet1_S","",20,jet1_central-jet1_range/2.,jet1_central+jet1_range/2.);
  TH1D jet2_S("jet2_S","",20,jet2_central-jet2_range/2.,jet2_central+jet2_range/2.);
  TH1D deltae_S("deltae_S","",20,deltae_central-deltae_range/2.,deltae_central+deltae_range/2.);
  TH1D zep_S("zep_S","",20,zep_central-zep_range/2.,zep_central+zep_range/2.);
  TH1D mjj_S("mjj_S","",20,mjj_central-mjj_range/2.,mjj_central+mjj_range/2.);
  TH1D deltap_S("deltap_S","",20,deltap_central-deltap_range/2.,deltap_central+deltap_range/2.);  
  TH1D met_S("met_S","",20,met_central-met_range/2.,met_central+met_range/2.);  

  TH1D ptgg_B("ptgg_B","",20,ptgg_central-ptgg_range/2.,ptgg_central+ptgg_range/2.);                 
  TH1D phot1_B("phot1_B","",20,phot1_central-phot1_range/2.,phot1_central+phot1_range/2.);
  TH1D phot2_B("phot2_B","",20,phot2_central-phot2_range/2.,phot2_central+phot2_range/2.);
  TH1D jet1_B("jet1_B","",20,jet1_central-jet1_range/2.,jet1_central+jet1_range/2.);
  TH1D jet2_B("jet2_B","",20,jet2_central-jet2_range/2.,jet2_central+jet2_range/2.);
  TH1D deltae_B("deltae_B","",20,deltae_central-deltae_range/2.,deltae_central+deltae_range/2.);
  TH1D zep_B("zep_B","",20,zep_central-zep_range/2.,zep_central+zep_range/2.);
  TH1D mjj_B("mjj_B","",20,mjj_central-mjj_range/2.,mjj_central+mjj_range/2.);
  TH1D deltap_B("deltap_B","",20,deltap_central-deltap_range/2.,deltap_central+deltap_range/2.);  
  TH1D met_B("met_B","",20,met_central-met_range/2.,met_central+met_range/2.);  

  TH1D ptgg_SovsqrtB("ptgg_SovsqrtB","",20,ptgg_central-ptgg_range/2.,ptgg_central+ptgg_range/2.);                 
  TH1D phot1_SovsqrtB("phot1_SovsqrtB","",20,phot1_central-phot1_range/2.,phot1_central+phot1_range/2.);
  TH1D phot2_SovsqrtB("phot2_SovsqrtB","",20,phot2_central-phot2_range/2.,phot2_central+phot2_range/2.);
  TH1D jet1_SovsqrtB("jet1_SovsqrtB","",20,jet1_central-jet1_range/2.,jet1_central+jet1_range/2.);
  TH1D jet2_SovsqrtB("jet2_SovsqrtB","",20,jet2_central-jet2_range/2.,jet2_central+jet2_range/2.);
  TH1D deltae_SovsqrtB("deltae_SovsqrtB","",20,deltae_central-deltae_range/2.,deltae_central+deltae_range/2.);
  TH1D zep_SovsqrtB("zep_SovsqrtB","",20,zep_central-zep_range/2.,zep_central+zep_range/2.);
  TH1D mjj_SovsqrtB("mjj_SovsqrtB","",20,mjj_central-mjj_range/2.,mjj_central+mjj_range/2.);
  TH1D deltap_SovsqrtB("deltap_SovsqrtB","",20,deltap_central-deltap_range/2.,deltap_central+deltap_range/2.);  
  TH1D met_SovsqrtB("met_SovsqrtB","",20,met_central-met_range/2.,met_central+met_range/2.);  
             
  TH1D ptgg_SovsqrtSB("ptgg_SovsqrtSB","",20,ptgg_central-ptgg_range/2.,ptgg_central+ptgg_range/2.);                 
  TH1D phot1_SovsqrtSB("phot1_SovsqrtSB","",20,phot1_central-phot1_range/2.,phot1_central+phot1_range/2.);
  TH1D phot2_SovsqrtSB("phot2_SovsqrtSB","",20,phot2_central-phot2_range/2.,phot2_central+phot2_range/2.);
  TH1D jet1_SovsqrtSB("jet1_SovsqrtSB","",20,jet1_central-jet1_range/2.,jet1_central+jet1_range/2.);
  TH1D jet2_SovsqrtSB("jet2_SovsqrtSB","",20,jet2_central-jet2_range/2.,jet2_central+jet2_range/2.);
  TH1D deltae_SovsqrtSB("deltae_SovsqrtSB","",20,deltae_central-deltae_range/2.,deltae_central+deltae_range/2.);
  TH1D zep_SovsqrtSB("zep_SovsqrtSB","",20,zep_central-zep_range/2.,zep_central+zep_range/2.);
  TH1D mjj_SovsqrtSB("mjj_SovsqrtSB","",20,mjj_central-mjj_range/2.,mjj_central+mjj_range/2.);
  TH1D deltap_SovsqrtSB("deltap_SovsqrtSB","",20,deltap_central-deltap_range/2.,deltap_central+deltap_range/2.);  
  TH1D met_SovsqrtSB("met_SovsqrtSB","",20,met_central-met_range/2.,met_central+met_range/2.);  
//   TH1D ptgg_SovsqrtSB("ptgg_SovsqrtSB","",20,0,80);                 
//   TH1D phot1_SovsqrtSB("phot1_SovsqrtSB","",20,35,95);
//   TH1D phot2_SovsqrtSB("phot2_SovsqrtSB","",20,20,50);
//   TH1D jet1_SovsqrtSB("jet1_SovsqrtSB","",20,20,60);
//   TH1D jet2_SovsqrtSB("jet2_SovsqrtSB","",20,20,40);
//   TH1D deltae_SovsqrtSB("deltae_SovsqrtSB","",20,2.5,4.5);
//   TH1D zep_SovsqrtSB("zep_SovsqrtSB","",20,1.5,3.5);
//   TH1D mjj_SovsqrtSB("mjj_SovsqrtSB","",20,50,850);
//   TH1D deltap_SovsqrtSB("deltap_SovsqrtSB","",20,2.1,3.1);  

  for (int kk=0; kk<10; kk++){

    char nametest[100];
    if(kk==0) sprintf(nametest,"ptgg_");
    if(kk==1) sprintf(nametest,"phot1_");
    if(kk==2) sprintf(nametest,"phot2_");
    if(kk==3) sprintf(nametest,"jet1_");
    if(kk==4) sprintf(nametest,"jet2_");
    if(kk==5) sprintf(nametest,"deltae_");
    if(kk==6) sprintf(nametest,"zep_");
    if(kk==7) sprintf(nametest,"mjj_");
    if(kk==8) sprintf(nametest,"deltap_");
    if(kk==9) sprintf(nametest,"met_");
    char nameprob_disc[100]; sprintf(nameprob_disc,"%s%s",nametest,"prob_disc");
    char nameprob_UL[100]; sprintf(nameprob_UL,"%s%s",nametest,"prob_UL");
    char nameS[100]; sprintf(nameS,"%s%s",nametest,"S");
    char nameB[100]; sprintf(nameB,"%s%s",nametest,"B");
    char nameSovsqrtB[100]; sprintf(nameSovsqrtB,"%s%s",nametest,"SovsqrtB");
    char nameSovsqrtSB[100]; sprintf(nameSovsqrtSB,"%s%s",nametest,"SovsqrtSB");

    for (int jj=0; jj<20; jj++){
      
      char namehisto[100];
      sprintf(namehisto,"%s%i",nametest,jj);
      
      for (int i=0; i<9; i++){ 
	var_mc[i]->Reset();
	var_mc[i]->Add((TH1D*)mcinput[i]->Get(namehisto));
      }
      
      for (int i=0; i<7; i++){ 
	var_gluglu[i]->Reset();
	var_gluglu[i]->Add((TH1D*)mc_glugluinput[i]->Get(namehisto));
	var_vbf[i]->Reset();
	var_vbf[i]->Add((TH1D*)mc_vbfinput[i]->Get(namehisto));
	var_wzh[i]->Reset();
	var_wzh[i]->Add((TH1D*)mc_wzhinput[i]->Get(namehisto));
      }
      
      // scale mc to equivalent lumi
      for (int i=0; i<9; i++){ 
	if(int_exp>0) var_mc[i]->Scale(scale_mc[i]);  
      }
      
      //   // counting number of events passing selection (scaled)
      //   double num_mc[9],num_gluglu[7],num_vbf[7],num_wzh[7],num_tth[7]; 
      //   double num_uns_mc[9],num_uns_gluglu[7],num_uns_vbf[7],num_uns_wzh[7],num_uns_tth[7]; 
      
      for (int i=0; i<9; i++){ 
	num_mc[i] = 0;
	if(int_exp>0) num_mc[i] = var_mc[i]->Integral();  
	num_uns_mc[i] = 0;
	if(int_exp>0) num_uns_mc[i] = var_mc[i]->GetEntries();  
      }
      for(int i=0; i<7; i++){
	num_gluglu[i] = num_vbf[i] = num_wzh[i] = num_tth[i] = 0;
	num_gluglu[i] = var_gluglu[i]->Integral(); 
	num_vbf[i] = var_vbf[i]->Integral(); 
	num_wzh[i] = var_wzh[i]->Integral(); 
	num_uns_gluglu[i] = num_uns_vbf[i] = num_uns_wzh[i] = num_uns_tth[i] = 0;
	num_uns_gluglu[i] = var_gluglu[i]->GetEntries(); 
	num_uns_vbf[i] = var_vbf[i]->GetEntries(); 
	num_uns_wzh[i] = var_wzh[i]->GetEntries(); 
      }
      
      // add two QCD bins
      if(int_exp>0) var_mc[3]->Add(var_mc[4]);
      
      for(int i=0; i<6; i++) var[i]->Reset();
      
      // stack histograms  
      for (int i=1; i<nbin+1; i++){      
	for (int j=0; j<6; j++){            
	  int offset(0);
	  if(j>0) offset = 2; // to add higgs contributions up
	  if(j>2) offset = 3; // to add qcd contributions up
	  for (int k=0 ; k<9-j-offset; k++){ 
	    if(int_exp>0) var[j]->SetBinContent(i,var_mc[k]->GetBinContent(i) + var[j]->GetBinContent(i));
	  }	
	}    
      }
      
      //higgs only plot
      var_mc[6]->Add(var_mc[7]);
      var_mc[6]->Add(var_mc[8]);
      
      binsize = (max-min)/nbin;
      sprintf(name,"%s%s%s%d%s","results_gg/optimalcut_",variable.c_str(),"_",eb,".txt");
      mydoublegau->SetParameter(0, 1.);
      mydoublegau->SetParameter(2, 1.5);
      mydoublegau->SetParameter(3, 1.);
      mydoublegau->SetParameter(5, 4.);      
      mydoublegau->FixParameter(1,120.);
      mydoublegau->FixParameter(4,120.);
      mydoublegau->SetParLimits(2,0.5,2.);
      mydoublegau->SetParLimits(5,2.,6.);    
      mypowlaw->SetParameter(0,2.88660e+07);
      mypowlaw->SetParLimits(0,0.,1.e+12);
      mypowlaw->SetParameter(1,-3.42928e+00);
      // TEMP
      // if(var_mc[6].Integral()<5) mypowlaw->FixParameter(1,mypowlaw->GetParameter(1));
      // VBF
      //mypowlaw->FixParameter(1,-3.42928e+00);
      // WZH
      //mypowlaw->FixParameter(1,-2.97321e+00);
      //
      var_mc[6]->Fit("mydoublegau","l");
      var[1]->Fit("mypowlaw","l");
      integralhiggs = mydoublegau->Integral(117.,123.)/binsize;
      integralbkg = mypowlaw->Integral(117.,123.)/binsize;
      cout << "Number of signal events " << integralhiggs/boosthiggs << endl;
      cout << "Number of bkg events " << integralbkg << endl;
      cout << "S/sqrt(B) " << integralhiggs/boosthiggs/sqrt(integralbkg) << endl; 
      cout << "probaility "<< probability_disc(integralhiggs/boosthiggs+integralbkg,integralbkg) << endl;
      if(integralbkg>0){
	((TH1D*)hOutputFile->Get(nameprob_disc))->SetBinContent(jj+1,probability_disc(integralhiggs/boosthiggs+integralbkg,integralbkg) );
	((TH1D*)hOutputFile->Get(nameprob_UL))->SetBinContent(jj+1,UL(integralhiggs/boosthiggs,integralbkg));
	((TH1D*)hOutputFile->Get(nameS))->SetBinContent(jj+1,integralhiggs/boosthiggs);
	((TH1D*)hOutputFile->Get(nameB))->SetBinContent(jj+1,integralbkg);
	((TH1D*)hOutputFile->Get(nameSovsqrtB))->SetBinContent(jj+1,integralhiggs/boosthiggs/sqrt(integralbkg));
	((TH1D*)hOutputFile->Get(nameSovsqrtSB))->SetBinContent(jj+1,integralhiggs/boosthiggs/sqrt(integralbkg+integralhiggs/boosthiggs));
      } else { 
	((TH1D*)hOutputFile->Get(nameprob_disc))->SetBinContent(jj+1,0);
	((TH1D*)hOutputFile->Get(nameprob_UL))->SetBinContent(jj+1,0);
	((TH1D*)hOutputFile->Get(nameS))->SetBinContent(jj+1,0);
	((TH1D*)hOutputFile->Get(nameB))->SetBinContent(jj+1,0);
	((TH1D*)hOutputFile->Get(nameSovsqrtB))->SetBinContent(jj+1,0);
	((TH1D*)hOutputFile->Get(nameSovsqrtSB))->SetBinContent(jj+1,0);
       }
      // ----------------------------
      
    }

  }
 
  hOutputFile->Write() ;
  hOutputFile->Close() ;
  hOutputFile->Delete();

//   delete c0;

  delete data;
  delete datacs;

  cout << "here" << endl;

  for(int i=0; i<9; i++){
    cout << " mc " << i << endl;
    delete mc_fill[i];
    cout << " mc2 " << i << endl;    
    delete mc[i];
  }

  for(int i=0; i<7; i++){
    cout << " signal " << i << endl;
    delete mc_gluglu_fill[i];
    delete mc_vbf_fill[i];
    delete mc_wzh_fill[i];
    //    delete mc_tth_fill[i];
    delete mc_gluglu[i];
    delete mc_vbf[i];
    delete mc_wzh[i];
    //    delete mc_tth[i];
  }

//   delete leg;
//   delete leg2;

  vector<double> values;

  if (variable == "massgg"){
    values.push_back(integralhiggs);
    values.push_back(integralbkg);
    values.push_back(integralhiggs/sqrt(integralbkg));
  }

  return values;

}








