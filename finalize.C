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
// .x finalize(data int lumi 2010, data int lumi 2011, pt1 cut, pt2 cut, ptj1 cut, ptj2 cut, deltae cut, zep cut, mjj cut, eb cat, r9 cat, trk isocut scale, ecal isocut scale, hcal isocut scale, h/e isocut scale)
// 
// example:
//
// .x finalize.C(33,230,40,25,20,15,2.5,2.5,300, 100, 100, 100, 100) 

vector <double> finalize(double int_exp_2010, double int_exp_2011, double pt1=50, double pt2=30, double ptj1=20, double ptj2=15, double deltae=2.5, double zep=2.5, double mjj=300, int eb = 1, int r9 = 1, int isolscaletrk = 100., int isolscaleecal = 100., int isolscalehcal = 100., int isolscalehove = 100., bool pixel = 1, string var = "massgg", int nbin = 200, double min = 90, double max = 190, string axis = "m(#gamma#gamma)[GeV]"){


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

  // total data sample
  TFile data("redntp.41xv7.preselection.v3/redntp_run2010-2011.root");    

  // higgs samples (only gluglu so far)
  TFile hig_2010("redntp.39xv7.preselection.v2/redntp_GluGluToHToGG_M-115_7TeV-pythia6_00.root");  //  2010
  TFile hig_2011("redntp.41xv7.preselection.v3/redntp_GluGluToHToGG_M-115_7TeV-pythia6_00.root");  //  2011

  // drell yan samples
  TFile dy_2010("redntp.39xv7.preselection.v2/redntp_DYJetsToLL_TuneZ2_M-50_7TeV-madgraph-tauola.root"); //  2010
  TFile dy_2011("redntp.41xv7.preselection.v3/redntp_DYJetsToLL_TuneZ2_M-50_7TeV-madgraph-tauola.root"); //  2011

  // box samples
  TFile box_2010("redntp.39xv7.preselection.v2/redntp_DiPhotonBox_Pt25to250_TrackingParticles_7TeV-pythia6_00.root"); //  2010
  TFile box_2011("redntp.41xv7.preselection.v3/redntp_DiPhotonBox_Pt25to250_TrackingParticles_7TeV-pythia6_00.root"); //  2011

  // diphoton jets samples
  TFile diphotjet_2010("redntp.39xv7.preselection.v2/redntp_DiPhotonJets_7TeV-madgraph.root"); //  2010
  TFile diphotjet_2011("redntp.41xv7.preselection.v3/redntp_DiPhotonJets_7TeV-madgraph.root"); //  2011

  // gjet samples
  TFile gjet_2010("redntp.39xv7.preselection.v2/redntp_GJet_Pt-20_doubleEMEnriched_TuneZ2_7TeV-pythia6.root"); //  2010
  TFile gjet_2011("redntp.41xv7.preselection.v3/redntp_GJet_Pt-20_doubleEMEnriched_TuneZ2_7TeV-pythia6.root"); //  2011

  // qcd pt>40 samples
  TFile qcd_2010("redntp.39xv7.preselection.v2/redntp_QCD_Pt-40_doubleEMEnriched_TuneZ2_7TeV-pythia6.root"); //  2010
  TFile qcd_2011("redntp.41xv7.preselection.v3/redntp_QCD_Pt-40_doubleEMEnriched_TuneZ2_7TeV-pythia6.root"); //  2011

  // qcd 30<pt<40 samples
  TFile qcd2_2010("redntp.39xv7.preselection.v2/redntp_QCD_Pt-30to40_doubleEMEnriched_TuneZ2_7TeV-pythia6.root"); //  2010
  TFile qcd2_2011("redntp.41xv7.preselection.v3/redntp_QCD_Pt-30to40_doubleEMEnriched_TuneZ2_7TeV-pythia6.root"); //  2011
 
  // cross sections and scaling
  double cross_hig = 18.23 * 2.13e-03 * 15; // higgs
  double cross_dy = 2321; // drell yan
  double cross_box = 12.37; // box
  double cross_diphotjet = 134; // diphoton jets
  double cross_gjet = 493.44; // gjet
  double cross_qcd = 40392; // qcd pt>40
  double cross_qcd2 = 9610; // qcd 30<pt<40 

  // getting the number of original events in each sample (processed with CMSSW)
  // higgs
  int n_hig_2010 = ((TH1D*)hig_2010.Get("ptphotgen1"))->GetEntries();
  int n_hig_2011 = ((TH1D*)hig_2011.Get("ptphotgen1"))->GetEntries();
  // drell yan
  int n_dy_2010 = ((TH1D*)dy_2010.Get("ptphotgen1"))->GetEntries();
  int n_dy_2011 = ((TH1D*)dy_2011.Get("ptphotgen1"))->GetEntries();
  // box 
  int n_box_2010 = ((TH1D*)box_2010.Get("ptphotgen1"))->GetEntries();
  int n_box_2011 = ((TH1D*)box_2011.Get("ptphotgen1"))->GetEntries();
  // diphoton jets
  int n_diphotjet_2010 = ((TH1D*)diphotjet_2010.Get("ptphotgen1"))->GetEntries();
  int n_diphotjet_2011 = ((TH1D*)diphotjet_2011.Get("ptphotgen1"))->GetEntries();
  // gjet
  int n_gjet_2010 = ((TH1D*)gjet_2010.Get("ptphotgen1"))->GetEntries();
  int n_gjet_2011 = ((TH1D*)gjet_2011.Get("ptphotgen1"))->GetEntries();
  // qcd pt>40
  int n_qcd_2010 = ((TH1D*)qcd_2010.Get("ptphotgen1"))->GetEntries();
  int n_qcd_2011 = ((TH1D*)qcd_2011.Get("ptphotgen1"))->GetEntries();
  // qcd 30<pt<40
  int n_qcd2_2010 = ((TH1D*)qcd2_2010.Get("ptphotgen1"))->GetEntries();
  int n_qcd2_2011 = ((TH1D*)qcd2_2011.Get("ptphotgen1"))->GetEntries();

  // setting the scaling factor to actual lumi
  double scale_hig_2010 = cross_hig * int_exp_2010 / n_hig_2010;
  double scale_dy_2010 = cross_dy * int_exp_2010 / n_dy_2010;
  double scale_box_2010 = cross_box * int_exp_2010 / n_box_2010;
  double scale_diphotjet_2010 = cross_diphotjet * int_exp_2010 / n_diphotjet_2010;
  double scale_gjet_2010 = cross_gjet * int_exp_2010 / n_gjet_2010;
  double scale_qcd_2010 = cross_qcd * int_exp_2010 / n_qcd_2010;
  double scale_qcd2_2010 = cross_qcd2 * int_exp_2010 / n_qcd2_2010;
  double scale_hig_2011 = cross_hig * int_exp_2011 / n_hig_2011;
  double scale_dy_2011 = cross_dy * int_exp_2011 / n_dy_2011;
  double scale_box_2011 = cross_box * int_exp_2011 / n_box_2011;
  double scale_diphotjet_2011 = cross_diphotjet * int_exp_2011 / n_diphotjet_2011;
  double scale_gjet_2011 = cross_gjet * int_exp_2011 / n_gjet_2011;
  double scale_qcd_2011 = cross_qcd * int_exp_2011 / n_qcd_2011;
  double scale_qcd2_2011 = cross_qcd2 * int_exp_2011 / n_qcd2_2011;

  // char for output name
  char name[1000];
  char allcut[3000];
  sprintf(allcut,"%f%s%f%s%f%s%f%s%f%s%f%s%f%s%d%s%d%s%d%s%d%s%d%s%d%s%d",pt1,"_",pt2,"_",ptj1,"_",ptj2,"_",deltae,"_",zep,"_",mjj,"_",eb,"_",r9,"_",isolscaletrk,"_",isolscaleecal,"_",isolscalehcal,"_",isolscalehove,"_",pixel);

  // output root file
  sprintf(name,"%s%s%s%s%s","results_gg/histo_",var.c_str(),"_",allcut,".root");

  TFile * hOutputFile   = new TFile(name, "RECREATE" ) ;

  // histograms needed by the machinery
  TH1D* vardata;
  TH1D* vardatacs;
  TH1D* varhig_2010;
  TH1D* vardy_2010;
  TH1D* varbox_2010;
  TH1D* vardiphotjet_2010;
  TH1D* vargjet_2010;
  TH1D* varqcd_2010;
  TH1D* varqcd2_2010;
  TH1D* varhig_2011;
  TH1D* vardy_2011;
  TH1D* varbox_2011;
  TH1D* vardiphotjet_2011;
  TH1D* vargjet_2011;
  TH1D* varqcd_2011;
  TH1D* varqcd2_2011;
  TH1D* var1 = new TH1D("var1","var1",nbin,min,max);
  TH1D* var2 = new TH1D("var2","var2",nbin,min,max);
  TH1D* var3 = new TH1D("var3","var3",nbin,min,max);
  TH1D* var4 = new TH1D("var4","var4",nbin,min,max);
  TH1D* var5 = new TH1D("var5","var5",nbin,min,max);
  TH1D* var6 = new TH1D("var6","var6",nbin,min,max);
  
  // creating the fillers
  fillPlot data_fill((TTree*)data.Get("AnaTree"));
  fillPlot higgs_2010_fill((TTree*)hig_2010.Get("AnaTree"));
  fillPlot higgs_2011_fill((TTree*)hig_2011.Get("AnaTree"));
  fillPlot dy_2010_fill((TTree*)dy_2010.Get("AnaTree"));
  fillPlot dy_2011_fill((TTree*)dy_2011.Get("AnaTree"));
  fillPlot box_2010_fill((TTree*)box_2010.Get("AnaTree"));
  fillPlot box_2011_fill((TTree*)box_2011.Get("AnaTree"));
  fillPlot diphotjet_2010_fill((TTree*)diphotjet_2010.Get("AnaTree"));
  fillPlot diphotjet_2011_fill((TTree*)diphotjet_2011.Get("AnaTree"));
  fillPlot gjet_2010_fill((TTree*)gjet_2010.Get("AnaTree"));
  fillPlot gjet_2011_fill((TTree*)gjet_2011.Get("AnaTree"));
  fillPlot qcd_2010_fill((TTree*)qcd_2010.Get("AnaTree"));
  fillPlot qcd_2011_fill((TTree*)qcd_2011.Get("AnaTree"));
  fillPlot qcd2_2010_fill((TTree*)qcd2_2010.Get("AnaTree"));
  fillPlot qcd2_2011_fill((TTree*)qcd2_2011.Get("AnaTree"));
 
  // setting cuts
  data_fill.Setcuts(pt1,pt2,ptj1,ptj2,deltae,zep,mjj,eb,r9,isolscaletrk,isolscaleecal,isolscalehcal,isolscalehove,pixel);
  data_fill.Setcuts(pt1,pt2,ptj1,ptj2,deltae,zep,mjj,eb,r9,isolscaletrk,isolscaleecal,isolscalehcal,isolscalehove,pixel);
  higgs_2010_fill.Setcuts(pt1,pt2,ptj1,ptj2,deltae,zep,mjj,eb,r9,isolscaletrk,isolscaleecal,isolscalehcal,isolscalehove,pixel);
  higgs_2011_fill.Setcuts(pt1,pt2,ptj1,ptj2,deltae,zep,mjj,eb,r9,isolscaletrk,isolscaleecal,isolscalehcal,isolscalehove,pixel);
  dy_2010_fill.Setcuts(pt1,pt2,ptj1,ptj2,deltae,zep,mjj,eb,r9,isolscaletrk,isolscaleecal,isolscalehcal,isolscalehove,pixel);
  dy_2011_fill.Setcuts(pt1,pt2,ptj1,ptj2,deltae,zep,mjj,eb,r9,isolscaletrk,isolscaleecal,isolscalehcal,isolscalehove,pixel);
  box_2010_fill.Setcuts(pt1,pt2,ptj1,ptj2,deltae,zep,mjj,eb,r9,isolscaletrk,isolscaleecal,isolscalehcal,isolscalehove,pixel);
  box_2011_fill.Setcuts(pt1,pt2,ptj1,ptj2,deltae,zep,mjj,eb,r9,isolscaletrk,isolscaleecal,isolscalehcal,isolscalehove,pixel);
  diphotjet_2010_fill.Setcuts(pt1,pt2,ptj1,ptj2,deltae,zep,mjj,eb,r9,isolscaletrk,isolscaleecal,isolscalehcal,isolscalehove,pixel);
  diphotjet_2011_fill.Setcuts(pt1,pt2,ptj1,ptj2,deltae,zep,mjj,eb,r9,isolscaletrk,isolscaleecal,isolscalehcal,isolscalehove,pixel);
  gjet_2010_fill.Setcuts(pt1,pt2,ptj1,ptj2,deltae,zep,mjj,eb,r9,isolscaletrk,isolscaleecal,isolscalehcal,isolscalehove,pixel);
  gjet_2011_fill.Setcuts(pt1,pt2,ptj1,ptj2,deltae,zep,mjj,eb,r9,isolscaletrk,isolscaleecal,isolscalehcal,isolscalehove,pixel);
  qcd_2010_fill.Setcuts(pt1,pt2,ptj1,ptj2,deltae,zep,mjj,eb,r9,isolscaletrk,isolscaleecal,isolscalehcal,isolscalehove,pixel);
  qcd_2011_fill.Setcuts(pt1,pt2,ptj1,ptj2,deltae,zep,mjj,eb,r9,isolscaletrk,isolscaleecal,isolscalehcal,isolscalehove,pixel);
  qcd2_2010_fill.Setcuts(pt1,pt2,ptj1,ptj2,deltae,zep,mjj,eb,r9,isolscaletrk,isolscaleecal,isolscalehcal,isolscalehove,pixel);
  qcd2_2011_fill.Setcuts(pt1,pt2,ptj1,ptj2,deltae,zep,mjj,eb,r9,isolscaletrk,isolscaleecal,isolscalehcal,isolscalehove,pixel);

  // filling histograms
  cout << "running over " << ((TTree*)data.Get("AnaTree"))->GetEntries("") << " events of data" <<  endl;
  data_fill.Writetxt(1);
  vardata = data_fill.Plot(var,"data", nbin, min, max); 
  cout << "running over " << ((TTree*)data.Get("AnaTree"))->GetEntries("") << " events of data for cs" <<  endl; 
  data_fill.Writetxt(0);
  vardatacs = data_fill.Plot(var,"datacs", nbin, min, max); 

  cout << "running over " << ((TTree*)hig_2010.Get("AnaTree"))->GetEntries("") << " events of higgs 2010" <<  endl; 
  varhig_2010 = higgs_2010_fill.Plot(var,"massgg_higgs_2010", nbin, min, max);
  cout << "running over " << ((TTree*)hig_2011.Get("AnaTree"))->GetEntries("") << " events of higgs 2011" <<  endl; 
  varhig_2011 = higgs_2011_fill.Plot(var,"massgg_higgs_2011", nbin, min, max);

  cout << "running over " << ((TTree*)dy_2010.Get("AnaTree"))->GetEntries("") << " events of dy 2010" <<  endl;
  vardy_2010 = dy_2010_fill.Plot(var,"massgg_dy_2010", nbin, min, max);  
  cout << "running over " << ((TTree*)dy_2011.Get("AnaTree"))->GetEntries("") << " events of dy 2011" <<  endl; 
  vardy_2011 = dy_2011_fill.Plot(var,"massgg_dy_2011", nbin, min, max);

  cout << "running over " << ((TTree*)box_2010.Get("AnaTree"))->GetEntries("") << " events of box 2010" <<  endl;
  varbox_2010 = box_2010_fill.Plot(var,"massgg_box_2010", nbin, min, max);  
  cout << "running over " << ((TTree*)box_2010.Get("AnaTree"))->GetEntries("") << " events of box 2011" <<  endl; 
  varbox_2011 = box_2011_fill.Plot(var,"massgg_box_2011", nbin, min, max);

  cout << "running over " << ((TTree*)diphotjet_2010.Get("AnaTree"))->GetEntries("") << " events of diphotjet 2010"<<  endl; 
  vardiphotjet_2010 = diphotjet_2010_fill.Plot(var,"massgg_diphotjet_2010", nbin, min, max); 
  cout << "running over " << ((TTree*)diphotjet_2011.Get("AnaTree"))->GetEntries("") << " events of diphotjet 2011" <<  endl; 
  vardiphotjet_2011 = diphotjet_2011_fill.Plot(var,"massgg_diphotjet_2011", nbin, min, max); 

  cout << "running over " << ((TTree*)gjet_2010.Get("AnaTree"))->GetEntries("") << " events of gjet 2010" <<  endl; 
  vargjet_2010 = gjet_2010_fill.Plot(var,"massgg_gjet_2010", nbin, min, max); 
  cout << "running over " << ((TTree*)gjet_2011.Get("AnaTree"))->GetEntries("") << " events of gjet 2011" <<  endl; 
  vargjet_2011 = gjet_2011_fill.Plot(var,"massgg_gjet_2011", nbin, min, max); 

  cout << "running over " << ((TTree*)qcd_2010.Get("AnaTree"))->GetEntries("") << " events of qcd pt40 2010" <<  endl; 
  varqcd_2010 = qcd_2010_fill.Plot(var,"massgg_qcd_2010", nbin, min, max); 
  cout << "running over " << ((TTree*)qcd_2011.Get("AnaTree"))->GetEntries("") << " events of qcd pt40 2011" <<  endl; 
  varqcd_2011 = qcd_2011_fill.Plot(var,"massgg_qcd_2011", nbin, min, max);

  cout << "running over " << ((TTree*)qcd2_2010.Get("AnaTree"))->GetEntries("") << " events of qcd pt30-40 2010" <<  endl; 
  varqcd2_2010 = qcd2_2010_fill.Plot(var,"massgg_qcd2_2010", nbin, min, max); 
  cout << "running over " << ((TTree*)qcd2_2010.Get("AnaTree"))->GetEntries("") << " events of qcd pt30-40 2011" <<  endl; 
  varqcd2_2011 = qcd2_2011_fill.Plot(var,"massgg_qcd2_2011", nbin, min, max); 

  // scale mc to equivalent lumi
  varhig_2010->Scale(scale_hig_2010);
  varhig_2011->Scale(scale_hig_2011);
  vardy_2010->Scale(scale_dy_2010);
  vardy_2011->Scale(scale_dy_2011);
  varbox_2010->Scale(scale_box_2010);
  varbox_2011->Scale(scale_box_2011);
  vardiphotjet_2010->Scale(scale_diphotjet_2010);
  vardiphotjet_2011->Scale(scale_diphotjet_2011);
  vargjet_2010->Scale(scale_gjet_2010);
  vargjet_2011->Scale(scale_gjet_2011);
  varqcd_2010->Scale(scale_qcd_2010);
  varqcd_2011->Scale(scale_qcd_2011);
  varqcd2_2010->Scale(scale_qcd2_2010);
  varqcd2_2011->Scale(scale_qcd2_2011);

  // counting number of events passing selection (scaled)
  double num_hig_2010       = varhig_2010->Integral();  
  double num_hig_2011       = varhig_2011->Integral();	
  double num_dy_2010        = vardy_2010->Integral();	
  double num_dy_2011        = vardy_2011->Integral();	
  double num_box_2010       = varbox_2010->Integral();
  double num_box_2011       = varbox_2011->Integral();
  double num_diphotjet_2010 = vardiphotjet_2010->Integral(); 
  double num_diphotjet_2011 = vardiphotjet_2011->Integral(); 
  double num_gjet_2010      = vargjet_2010->Integral();	  
  double num_gjet_2011      = vargjet_2011->Integral();	  
  double num_qcd_2010       = varqcd_2010->Integral();	  
  double num_qcd_2011       = varqcd_2011->Integral();	  
  double num_qcd2_2010      = varqcd2_2010->Integral();	  
  double num_qcd2_2011      = varqcd2_2011->Integral();	  

  // add two QCD bins
  varqcd_2010->Add(varqcd2_2010);
  varqcd_2011->Add(varqcd2_2011);
  
  // scale control sample
  vardata->Sumw2();
  double num_data =  vardata->GetEntries();
  double num_data_cs = vardatacs->GetEntries();  
  vardatacs->Scale(num_data/num_data_cs); 

  // stack histograms

  for (int j=0; j<6; j++){

     for (int i=1; i<nbin+1; i++){      
      
      int k = i;
      
      if(j==0) {
	var1->SetBinContent(i,varbox_2010->GetBinContent(k) + var1->GetBinContent(i));
	var2->SetBinContent(i,varbox_2010->GetBinContent(k) + var2->GetBinContent(i));
	var3->SetBinContent(i,varbox_2010->GetBinContent(k) + var3->GetBinContent(i));
	var4->SetBinContent(i,varbox_2010->GetBinContent(k) + var4->GetBinContent(i));
	var5->SetBinContent(i,varbox_2010->GetBinContent(k) + var5->GetBinContent(i));
	var6->SetBinContent(i,varbox_2010->GetBinContent(k) + var6->GetBinContent(i));
	var1->SetBinContent(i,varbox_2011->GetBinContent(k) + var1->GetBinContent(i));
	var2->SetBinContent(i,varbox_2011->GetBinContent(k) + var2->GetBinContent(i));
	var3->SetBinContent(i,varbox_2011->GetBinContent(k) + var3->GetBinContent(i));
	var4->SetBinContent(i,varbox_2011->GetBinContent(k) + var4->GetBinContent(i));
	var5->SetBinContent(i,varbox_2011->GetBinContent(k) + var5->GetBinContent(i));
	var6->SetBinContent(i,varbox_2011->GetBinContent(k) + var6->GetBinContent(i));
      }
      if(j==1) {
      	var2->SetBinContent(i,vardiphotjet_2010->GetBinContent(k) + var2->GetBinContent(i));
      	var3->SetBinContent(i,vardiphotjet_2010->GetBinContent(k) + var3->GetBinContent(i));
      	var4->SetBinContent(i,vardiphotjet_2010->GetBinContent(k) + var4->GetBinContent(i));
      	var5->SetBinContent(i,vardiphotjet_2010->GetBinContent(k) + var5->GetBinContent(i));
      	var6->SetBinContent(i,vardiphotjet_2010->GetBinContent(k) + var6->GetBinContent(i));
      	var2->SetBinContent(i,vardiphotjet_2011->GetBinContent(k) + var2->GetBinContent(i));
      	var3->SetBinContent(i,vardiphotjet_2011->GetBinContent(k) + var3->GetBinContent(i));
      	var4->SetBinContent(i,vardiphotjet_2011->GetBinContent(k) + var4->GetBinContent(i));
      	var5->SetBinContent(i,vardiphotjet_2011->GetBinContent(k) + var5->GetBinContent(i));
      	var6->SetBinContent(i,vardiphotjet_2011->GetBinContent(k) + var6->GetBinContent(i));
      }
      if(j==2) {
	var3->SetBinContent(i,vargjet_2010->GetBinContent(k) + var3->GetBinContent(i));
	var4->SetBinContent(i,vargjet_2010->GetBinContent(k) + var4->GetBinContent(i));
	var5->SetBinContent(i,vargjet_2010->GetBinContent(k) + var5->GetBinContent(i));
	var6->SetBinContent(i,vargjet_2010->GetBinContent(k) + var6->GetBinContent(i));
	var3->SetBinContent(i,vargjet_2011->GetBinContent(k) + var3->GetBinContent(i));
	var4->SetBinContent(i,vargjet_2011->GetBinContent(k) + var4->GetBinContent(i));
	var5->SetBinContent(i,vargjet_2011->GetBinContent(k) + var5->GetBinContent(i));
	var6->SetBinContent(i,vargjet_2011->GetBinContent(k) + var6->GetBinContent(i));
      }
      if(j==3) {
	var4->SetBinContent(i,varqcd_2010->GetBinContent(k) + var4->GetBinContent(i));
	var5->SetBinContent(i,varqcd_2010->GetBinContent(k) + var5->GetBinContent(i));
	var6->SetBinContent(i,varqcd_2010->GetBinContent(k) + var6->GetBinContent(i));
	var4->SetBinContent(i,varqcd_2011->GetBinContent(k) + var4->GetBinContent(i));
	var5->SetBinContent(i,varqcd_2011->GetBinContent(k) + var5->GetBinContent(i));
	var6->SetBinContent(i,varqcd_2011->GetBinContent(k) + var6->GetBinContent(i));
      }
      if(j==4) {
	var5->SetBinContent(i,vardy_2010->GetBinContent(k) + var5->GetBinContent(i));
 	var6->SetBinContent(i,vardy_2010->GetBinContent(k) + var6->GetBinContent(i));
	var5->SetBinContent(i,vardy_2011->GetBinContent(k) + var5->GetBinContent(i));
 	var6->SetBinContent(i,vardy_2011->GetBinContent(k) + var6->GetBinContent(i));
      }
      if(j==5) {
 	var6->SetBinContent(i,varhig_2010->GetBinContent(k) + var6->GetBinContent(i));
 	var6->SetBinContent(i,varhig_2011->GetBinContent(k) + var6->GetBinContent(i));
      }
	
    }
    
    
  }
 
  //final plots
  var6->SetTitle("");
  var6->SetStats(0);
  var6->SetTitleOffset(1.25,"Y");
  var6->SetXTitle(axis.c_str());
  char ytitle[100];
  sprintf(ytitle,"%s%d%s","N_{ev}/",int(int_exp_2010+int_exp_2011),"pb^{-1}");
  var6->SetYTitle(ytitle);
  var6->SetLineColor(kBlack);
  var6->SetLineWidth(2);
  var6->SetFillColor(kYellow);
  var6->Draw();

  var5->SetLineColor(kBlack);
  var5->SetLineWidth(2);
  var5->SetFillColor(32);
  var5->Draw("same");

  var4->SetLineColor(kBlack);
  var4->SetLineWidth(2);
  var4->SetFillColor(29);
  var4->Draw("same");

  var3->SetLineColor(kBlack);
  var3->SetLineWidth(2);
  var3->SetFillColor(38);
  var3->Draw("same");

  var2->SetLineColor(kBlack);
  var2->SetLineWidth(2);
  var2->SetFillColor(46);
  var2->Draw("same");

  var1->SetLineColor(kBlack);
  var1->SetLineWidth(2);
  var1->SetFillColor(16);
  var1->Draw("same");

  TLegendEntry *legge;
  TLegend *leg;
  leg = new TLegend(0.5,0.6,0.75,0.85);
  leg->SetFillStyle(0); leg->SetBorderSize(0); leg->SetTextSize(0.05);
  leg->SetFillColor(0);
  legge = leg->AddEntry(var6, "h_{f} M(115) x 15", "f");
  legge = leg->AddEntry(var5, "DY", "f");
  legge = leg->AddEntry(var4, "QCD", "f");
  legge = leg->AddEntry(var3, "#gamma + jets", "f");
  legge = leg->AddEntry(var2, "di-#gamma + jets", "f");
  legge = leg->AddEntry(var1, "di-#gamma box", "f");
  leg->Draw();

  sprintf(name,"%s%s%s%s%s","results_gg/mc_",var.c_str(),"_",allcut,".gif");
  c0->SaveAs(name);

  varhig_2010->Add(varhig_2011);
  sprintf(ytitle,"%s%d%s","N_{ev}/",int(int_exp_2010+int_exp_2011),"pb^{-1}");
  varhig_2010->SetXTitle(axis.c_str());
  varhig_2010->SetYTitle(ytitle);
  varhig_2010->SetTitle("");
  varhig_2010->SetLineColor(kBlack);
  varhig_2010->SetLineWidth(2);
  varhig_2010->SetFillColor(kYellow);
  varhig_2010->Draw();
  sprintf(name,"%s%s%s%s%s","results_gg/higgs_",var.c_str(),"_",allcut,".gif");
  c0->SaveAs(name);
  
  sprintf(name,"%s%s%s%d%s","results_gg/optimalcut_",var.c_str(),"_",eb,".txt");
  double integralhiggs = varhig_2010->Integral(43,56);
  double entrieshiggs =  varhig_2010->Integral(0,201);
  double integralbkg = var5->Integral(30,71)/3.;
  cout << "Fraction in signal box " << integralhiggs/entrieshiggs << endl;
  cout << "Number of signal events " << integralhiggs/15 << endl;
  cout << "Number of bkg events " << integralbkg << endl;
  cout << "S/sqrt(B) " << integralhiggs/15/sqrt(integralbkg) << endl;

  data.cd();
  vardata->SetMinimum(0.01);
  vardata->SetXTitle(axis.c_str());
  vardata->SetTitle("");
  vardata->SetStats(0);
  vardata->SetMarkerStyle(8);
  vardata->SetMarkerSize(.9);
  vardata->SetTitleOffset(1.25,"Y");
  vardata->Draw("pe");
  sprintf(name,"%s%s%s%s%s","results_gg/data_",var.c_str(),"_",allcut,".gif");
  c0->SaveAs(name);

  double themax =   vardata->GetMaximum();
  if(var5->GetMaximum()>max) themax = var5->GetMaximum();
  vardata->SetMaximum(themax*1.3);
  var6->Draw("same");
  var5->Draw("same");
  var4->Draw("same");
  var3->Draw("same");
  var2->Draw("same");
  var1->Draw("same");
  leg->Draw();
  vardata->Draw("pesame");

  gPad->RedrawAxis();

  sprintf(name,"%s%s%s%s%s","results_gg/data-mc_",var.c_str(),"_",allcut,".gif");
  c0->SaveAs(name);

  vardata->Draw("pe");
  vardatacs->Scale(num_data/num_data_cs); 
  vardatacs->SetLineColor(46);
  vardatacs->SetFillColor(42);
  vardatacs->SetLineWidth(3);
  vardatacs->Draw("same");
  vardata->Draw("pesame");
  sprintf(name,"%s%s%s%s%s","results_gg/datacs_",var.c_str(),"_",allcut,".gif");
  gPad->RedrawAxis();
  TLegendEntry *legge2;
  TLegend *leg2;
  leg2 = new TLegend(0.6,0.65,0.9,0.85);
  leg2->SetFillStyle(0); leg2->SetBorderSize(0); leg2->SetTextSize(0.035);
  leg2->SetFillColor(0);
  legge2 = leg2->AddEntry(((TH1D*)gDirectory->Get("vardata")), "default sel.", "p");
  legge2 = leg2->AddEntry(((TH1D*)gDirectory->Get("vardatacs")), "control sample", "f");
  leg2->Draw();
  c0->SaveAs(name);


  if (var == "massgg"){
    var6->Rebin(4);
    var6->Draw();
    var5->Rebin(4);
    var5->Draw("same");
    var4->Rebin(4);
    var4->Draw("same");
    var3->Rebin(4);
    var3->Draw("same");
    var2->Rebin(4);
    var2->Draw("same");
    var1->Rebin(4);
    var1->Draw("same");
    leg->Draw();
    gPad->RedrawAxis();
    
    sprintf(allcut,"%f%s%f%s%f%s%f%s%f%s%f%s%f%s%d%s%d%s%d%s%d%s%d%s%d",pt1,"_",pt2,"_",ptj1,"_",ptj2,"_",deltae,"_",zep,"_",mjj,"_",eb,"_",r9,"_",isolscaletrk,"_",isolscaleecal,"_",isolscalehcal,"_",isolscalehove);
    sprintf(name,"%s%s%s%s%s","results_gg/mc_rebin_",var.c_str(),"_",allcut,".gif");
    c0->SaveAs(name);
    
    vardata->Rebin(4);
    vardata->Draw("pesame");
    
    sprintf(name,"%s%s%s%s%s","results_gg/data-mc_rebin_",var.c_str(),"_",allcut,".gif");
    c0->SaveAs(name);
  }

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
  outfile << "# events hig_2010 = " << n_hig_2010 << endl;
  outfile << "# events dy_2010 = " << n_dy_2010 << endl;
  outfile << "# events box_2010 = " << n_box_2010 << endl;
  outfile << "# events diphotjet_2010 = " << n_diphotjet_2010 << endl;
  outfile << "# events gjet_2010 = " << n_gjet_2010 << endl;
  outfile << "# events qcd_2010 = " << n_qcd_2010 << endl;
  outfile << "# events qcd2_2010 = " << n_qcd2_2010 << endl;
  outfile << "# events hig_2011 = " << n_hig_2011 << endl;
  outfile << "# events dy_2011 = " << n_dy_2011 << endl;
  outfile << "# events box_2011 = " << n_box_2011 << endl;
  outfile << "# events diphotjet_2011 = " << n_diphotjet_2011 << endl;
  outfile << "# events gjet_2011 = " << n_gjet_2011 << endl;
  outfile << "# events qcd_2011 = " << n_qcd_2011 << endl;
  outfile << "# events qcd2_2011 = " << n_qcd2_2011 << endl;
  outfile << endl;
  outfile << "####################################" << endl;
  outfile << "N of selected events and eff." << endl;
  outfile << "####################################" << endl; 
  outfile << "ndata      = " << num_data << endl;
  outfile << endl;
  outfile << "nallbkg    = " << num_box_2011 + num_diphotjet_2011 + num_gjet_2011 + num_qcd_2011 + num_qcd2_2011 + num_dy_2011 + num_box_2011 + num_diphotjet_2011 + num_gjet_2011 + num_qcd_2011 + num_qcd2_2011 + num_dy_2011<< endl;
  outfile << "nhig       = " << num_hig_2010       + num_hig_2011       << endl;
  outfile << "ndy        = " << num_dy_2010        + num_dy_2011        << endl;
  outfile << "nbox       = " << num_box_2010       + num_box_2011       << endl;
  outfile << "ndiphot    = " << num_diphotjet_2010 + num_diphotjet_2011 << endl;
  outfile << "ngjet      = " << num_gjet_2010      + num_gjet_2011      << endl;
  outfile << "nqcd40     = " << num_qcd_2010       + num_qcd_2011       << endl;
  outfile << "nqcd30-40  = " << num_qcd2_2010      + num_qcd2_2011      << endl;
  outfile << endl;
  outfile << "eff nhig      = " << (num_hig_2011       + num_hig_2010)      /(n_hig_2011       * scale_hig_2011       + n_hig_2010       * scale_hig_2010       ) << endl;
  outfile << "eff ndy       = " << (num_dy_2011        + num_dy_2010)       /(n_dy_2011        * scale_dy_2011        + n_dy_2010        * scale_dy_2010        ) << endl;
  outfile << "eff nbox      = " << (num_box_2011       + num_box_2010)      /(n_box_2011       * scale_box_2011       + n_box_2010       * scale_box_2010       ) << endl;
  outfile << "eff ndiphot   = " << (num_diphotjet_2011 + num_diphotjet_2010)/(n_diphotjet_2011 * scale_diphotjet_2011 + n_diphotjet_2010 * scale_diphotjet_2010 ) << endl;
  outfile << "eff ngjet     = " << (num_gjet_2011      + num_gjet_2010)     /(n_gjet_2011      * scale_gjet_2011      + n_gjet_2010      * scale_gjet_2010      ) << endl;
  outfile << "eff nqcd      = " << (num_qcd_2011       + num_qcd_2010)      /(n_qcd_2011       * scale_qcd_2011       + n_qcd_2010       * scale_qcd_2010       ) << endl;
  outfile << "eff nqcd30-40 = " << (num_qcd2_2011      + num_qcd2_2010)     /(n_qcd2_2011      * scale_qcd2_2011      + n_qcd2_2010      * scale_qcd2_2010      ) << endl;
  outfile << endl;

  outfile.close();

 
  hOutputFile->Write() ;
  hOutputFile->Close() ;
  hOutputFile->Delete();

  vector<double> values;
  
  values.push_back(integralhiggs);
  values.push_back(integralbkg);
  values.push_back(integralhiggs/sqrt(integralbkg));
  
  return values;
}








