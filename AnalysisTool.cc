// code written on Sun21 Nov: firs PhId, then ordering photons. Reco_gamma1 and 2. Find jet1 and jet2 not matching Reco_gamma1 and 2.
// generator part to be rewritten (too verbose)
#include "AnalysisTool.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <iostream>
#include <fstream>
#include <math.h>
#include<cmath>
#include<vector>
#define pi 3.141592653589793
#include <TLorentzVector.h>
 
using namespace std;

AnalysisTool::AnalysisTool(TTree *tree, const TString& outname) 
   : tree_reader_V6(tree) { 
     outname_=outname; 
     ptJet1_cut = ptJet2_cut = 20.;
     detaJets_cut = 3.;
     zepp_cut = 3.;
     mjj_cut = 250.; 

     // PhotonID Del Re
     mediumid.hcaliso_rel=         0.05;
     mediumid.hcaliso_abs=         2.4;
     mediumid.ecaliso_rel=         0.05;
     mediumid.ecaliso_abs=         2.4;
     mediumid.tracknb=             3.;
     mediumid.trackiso_rel=        0.10;
     mediumid.sminmin=             0.30;
     mediumid.sminmin_min=         0.15;
     mediumid.smajmaj=             0.35;

     looseid.hcaliso_rel=         0.10;
     looseid.hcaliso_abs=         4.;
     looseid.ecaliso_rel=         0.10;
     looseid.ecaliso_abs=         4.5;
     looseid.tracknb=             5.;
     looseid.trackiso_rel=        0.20;
     looseid.sminmin=             0.50;
     looseid.sminmin_min=         0.15;
     looseid.smajmaj=             0.60; 

    filter2GammaJet = false;
} 

void AnalysisTool::ReadCutsFromFile(const char* fname) {

  ifstream is(fname);
  if(! is.good()) {
     cout << "AnalysisTool::ReadCutsFromFile >> ERROR : file " << fname << " not read" << endl;
     is.close();
     exit(-1);
  }

  cout << "--- setting cuts from file:  " << fname << endl; 

  // order is fixed

  //ptJet1_cut
  is >> ptJet1_cut;
  cout << "ptJet1_cut: " << ptJet1_cut << endl;

  //ptJet2_cut
  is >> ptJet2_cut;
  cout << "ptJet2_cut: " << ptJet2_cut << endl;

  //detaJets_cut
  is >> detaJets_cut;
  cout << "detaJets_cut: " << detaJets_cut << endl;

  //zepp_cut
  is >> zepp_cut;
  cout << "zepp_cut: " << zepp_cut << endl;

  //mjj_cut
  is >> mjj_cut;
  cout << "mjj_cut: " << mjj_cut << endl;

  cout << "--- done setting cuts." << endl;
}


void AnalysisTool::Loop() {



  // open output file
  TFile *output_root = new TFile(outname_,"RECREATE");

  using namespace std;

  double SumP[4],Gen_inv_mass=-100;

  gStyle->SetOptStat(1111111);
  // TFile* plotsfile = TFile::Open("vbf_plotsfile.root","RECREATE");  

/****************************** HISTOGRAMS *******************************/

  int color = hfillcolor;
 
  gStyle->SetFillColor(hfillcolor);

  //add SetMinimum(0)

  TH1F* h1 = new TH1F("h1","Higgs E",100,0,1000.);
  h1->GetXaxis()->SetTitle("E(GeV)"); 
  h1->GetYaxis()->SetTitle("counts");   
  h1->SetFillColor(hfillcolor);
  h1->SetMinimum(0);

  TH1F* h2 = new TH1F("h2","Higgs Pt",100,0,140);
  h2->GetXaxis()->SetTitle("Pt(#frac{GeV}{c})"); 
  h2->GetYaxis()->SetTitle("counts");  
  h2->SetFillColor(hfillcolor);
  h2->SetMinimum(0);
  
  TH1F* h3 = new TH1F("h3","Higgs Eta",100,-6.5,6.5);
  h3->GetXaxis()->SetTitle("#eta");                  
  h3->GetYaxis()->SetTitle("counts");
  h3->SetFillColor(hfillcolor);
  h3->SetMinimum(0);

  TH1F* h4 = new TH1F("h4","Higgs Phi",80,-pi,pi);
  h4->GetYaxis()->SetTitle("counts");             
  h4->GetXaxis()->SetTitle("#phi"); 
  h4->SetFillColor(hfillcolor);
  h4->SetMinimum(0);

/////////////////////////////////////////////////////////////////

  TH1F* h5 = new TH1F("h5","genPhoton1 E",100,0,1000.);
  h5->GetXaxis()->SetTitle("E(GeV)"); 
  h5->GetYaxis()->SetTitle("counts");  
  h5->SetFillColor(hfillcolor);

  TH1F* h6 = new TH1F("h6","genPhoton1 Pt",100,0,140);
  h6->GetXaxis()->SetTitle("Pt(#frac{GeV}{c})"); 
  h6->GetYaxis()->SetTitle("counts"); 
  h6->SetFillColor(hfillcolor);

  TH1F* h7 = new TH1F("h7","genPhoton1 Eta",100,-6.5,6.5);
  h7->GetXaxis()->SetTitle("#eta");                  
  h7->GetYaxis()->SetTitle("counts");  
  h7->SetFillColor(hfillcolor);
  h7->SetMinimum(0);

  TH1F* h8 = new TH1F("h8","genPhoton1 Phi",80,-pi,pi);
  h8->GetXaxis()->SetTitle("#phi");                  
  h8->GetYaxis()->SetTitle("counts");   
  h8->SetFillColor(hfillcolor);
  h8->SetMinimum(0);

////////////////////////////////////////////////////////////////////

  TH1F* h9 = new TH1F("h9","genPhoton2 E",100,0,1000.);
  h9->GetXaxis()->SetTitle("E(GeV)"); 
  h9->GetYaxis()->SetTitle("counts"); 
  h9->SetFillColor(hfillcolor);

  TH1F* h10 = new TH1F("h10","genPhoton2 Pt",100,0,140);
  h10->GetXaxis()->SetTitle("Pt(#frac{GeV}{c})"); 
  h10->GetYaxis()->SetTitle("counts"); 
  h10->SetFillColor(hfillcolor);

  TH1F* h11 = new TH1F("h11","genPhoton2 Eta",100,-6.5,6.5);
  h11->GetXaxis()->SetTitle("#eta");                  
  h11->GetYaxis()->SetTitle("counts");
  h11->SetFillColor(hfillcolor);
  h11->SetMinimum(0);

  TH1F* h12 = new TH1F("h12","genPhoton2 Phi",80,-pi,pi);
  h12->GetXaxis()->SetTitle("#phi");                  
  h12->GetYaxis()->SetTitle("counts"); 
  h12->SetFillColor(hfillcolor);
  h12->SetMinimum(0);

/////////////////////////////////////////////////////////////////

  TH1F* h13 = new TH1F("h13","[NHD]genPhoton E",100,0,150);
  h13->GetXaxis()->SetTitle("E(GeV)"); 
  h13->GetYaxis()->SetTitle("counts"); 
  h13->SetFillColor(hfillcolor);

  TH1F* h14 = new TH1F("h14","[NHD]genPhoton Pt",100,0,10);
  h14->GetXaxis()->SetTitle("Pt(#frac{GeV}{c})"); 
  h14->GetYaxis()->SetTitle("counts"); 
  h14->SetFillColor(hfillcolor);

  TH1F* h15 = new TH1F("h15","[NHD]genPhoton Eta",100,-6.5,6.5);
  h15->GetXaxis()->SetTitle("#eta");                  
  h15->GetYaxis()->SetTitle("counts");
  h15->SetFillColor(hfillcolor);
  h15->SetMinimum(0);

  TH1F* h16 = new TH1F("h16","[NHD]genPhoton Phi",80,pi,pi);
  h16->GetXaxis()->SetTitle("#phi");                  
  h16->GetYaxis()->SetTitle("counts");  
  h16->SetFillColor(hfillcolor);
  h16->SetMinimum(0);

///////////////////////////////////////////////////////////////////////////////////////////////////

  TH1F* h17 = new TH1F ("h17","Istogramma Massa invariante 2 fotoni GEN",100,119.9,120.1);
  h17->GetXaxis()->SetTitle("M(#frac{GeV}{c^{2}})"); 
  h17->GetYaxis()->SetTitle("counts"); 
  h17->SetFillColor(hfillcolor);

  TH1F* h18 = new TH1F ("h18","Istogramma Massa invariante 2 fotoni RECO",100,112.,130.);
  h18->GetXaxis()->SetTitle("M(#frac{GeV}{c^{2}})"); 
  h18->GetYaxis()->SetTitle("counts"); 
  h18->SetFillColor(hfillcolor);

///////////////////////////////////////////////////////////////////////////////////////////////////

  TH1F* h19 = new TH1F ("h19","Associazione: delta_R fotone + energetico",100,0,9);
  h19->GetXaxis()->SetTitle("#DeltaR={#Delta#phi^{2}+#Delta#eta^{2}}^{#frac{1}{2} }");
  h19->GetYaxis()->SetTitle("counts");
  h19->SetFillColor(hfillcolor);

  TH1F* h20 = new TH1F ("h20","Associazione: delta_R fotone - energetico",100,0,9);
  h20->GetXaxis()->SetTitle("#DeltaR={#Delta#phi^{2}+#Delta#eta^{2}}^{#frac{1}{2} }");
  h20->GetYaxis()->SetTitle("counts"); 
  h20->SetFillColor(hfillcolor);

  TH2F* h21=new TH2F("h21","Delta_R1 vs Delta_R2",100,0,0.1,100,0,0.1);
  gStyle->SetPalette(1);
  h21->Draw("COLZ");
  h21->GetXaxis()->SetTitle("#DeltaR1");
  h21->GetYaxis()->SetTitle("#DeltaR2"); 

///////////////////////////////////////////////////////////////////////////////////////////////////

  TH1F* h22 = new TH1F ("h22","RECO Massa Inv 2 fotoni matchati",100,100,150);
  h22->GetXaxis()->SetTitleOffset(1.2);
  h22->GetXaxis()->SetTitle("M(#frac{GeV}{c^{2}})");
  h22->GetYaxis()->SetTitle("counts"); 
  h22->SetFillColor(hfillcolor);

  TH1F* h23 = new TH1F ("h23","RECO Massa Inv 1 fotone matchato",100,100,150);
  h23->GetXaxis()->SetTitleOffset(1.2);
  h23->GetXaxis()->SetTitle("M(#frac{GeV}{c^{2}})");
  h23->GetYaxis()->SetTitle("counts"); 
  h23->SetFillColor(hfillcolor);

  TH1F *h24 = new TH1F("h24","RECO Massa Inv 0 fotoni matchati",100,100,150);
  h24->GetXaxis()->SetTitleOffset(1.2);
  h24->GetXaxis()->SetTitle("M(#frac{GeV}{c^{2}})");
  h24->GetYaxis()->SetTitle("counts"); 
  h24->SetFillColor(hfillcolor);

///////////////////////////////////////////////////////////////////////////////////////////////////

  TH1F* h25 = new TH1F("h25","Jets E",100,0,1000.);                           /* JETS */
  h25->GetXaxis()->SetTitle("E(GeV)"); 
  h25->GetYaxis()->SetTitle("counts");
  h25->SetFillColor(hfillcolor);

  TH1F* h26 = new TH1F("h26","Jets Pt",100,0,160);
  h26->GetXaxis()->SetTitle("Pt(#frac{GeV}{c})"); 
  h26->GetYaxis()->SetTitle("counts");  
  h26->SetFillColor(hfillcolor);

  TH1F* h27 = new TH1F("h27","Jets Eta",100,-6.5,6.5);
  h27->GetXaxis()->SetTitle("#eta");                  
  h27->GetYaxis()->SetTitle("counts");
  h27->SetFillColor(hfillcolor);
  h27->SetMinimum(0);

  TH1F* h28 = new TH1F("h28","Jets Phi",80,-pi,pi);
  h28->GetXaxis()->SetTitle("#phi");                  
  h28->GetYaxis()->SetTitle("counts");
  h28->SetFillColor(hfillcolor);
  h28->SetMinimum(0);

///////////////////////////////////////////////////////////////////////////////////////////////////

  TH1F* h29 = new TH1F("h29","# Jets per evento",100,0,10);
  h29->GetXaxis()->SetTitle("# of jets per event");    
  h29->GetYaxis()->SetTitleOffset(1.5);              
  h29->GetYaxis()->SetTitle("counts");
  h29->SetFillColor(hfillcolor);

  TH1F* h30 = new TH1F("h30","# Jets Eta<1.4",100,0,10);
  h30->GetXaxis()->SetTitle("# of jets per event");   
  h30->GetYaxis()->SetTitleOffset(1.5);                
  h30->GetYaxis()->SetTitle("counts");
  h30->SetFillColor(hfillcolor);

  TH1F* h31 = new TH1F("h31","# Jets Eta>1.4",100,0,10);
  h31->GetXaxis()->SetTitle("# of jets per event");                  
  h31->GetYaxis()->SetTitleOffset(1.5); 
  h31->GetYaxis()->SetTitle("counts");
  h31->SetFillColor(hfillcolor);

//////////////////////////////////////////////////////////////////////////////////////////////////

  TH1F* h32 = new TH1F("h32"," Jet-Jet DeltaEta ",100,-10,10);
  h32->GetXaxis()->SetTitle("#Delta#eta(Jet-Jet)"); 
  h32->GetYaxis()->SetTitle("counts");  
  h32->SetFillColor(hfillcolor);

  TH1F* h33 = new TH1F("h33"," Jet-Higgs DeltaEta ",100,-10,10);
  h33->GetXaxis()->SetTitle("#Delta#eta(Jet-Higgs)"); 
  h33->GetYaxis()->SetTitle("counts");  
  h33->SetFillColor(hfillcolor);

//////////////////////////////////////////////////////////////////////////////////////////////////  

  TH1F* h34 = new TH1F("h34"," Pt Jets in Barrel ",100,0,160);
  h34->GetXaxis()->SetTitle("Pt(#frac{GeV}{c})"); 
  h34->GetYaxis()->SetTitle("counts");
  h34->SetFillColor(hfillcolor);

  TH1F* h35 = new TH1F("h35"," Pt Jets in EndCap ",100,0,160);
  h35->GetXaxis()->SetTitle("Pt(#frac{GeV}{c})"); 
  h35->GetYaxis()->SetTitle("counts");
  h35->SetFillColor(hfillcolor);

////////////////////////////////////////////////////////////////////////////////////////////////// 

  TH1F* h71 = new TH1F("h71"," DeltaR gamma1 Jet1 ",100,0,10);
  h71->GetXaxis()->SetTitle("#DeltaR(#gamma1-Jet1)"); 
  h71->GetYaxis()->SetTitle("counts");
  h71->SetFillColor(hfillcolor);

  TH1F* h72 = new TH1F("h72"," DeltaR gamma1 Jet2 ",100,0,10);
  h72->GetXaxis()->SetTitle("#DeltaR(#gamma1-Jet2)"); 
  h72->GetYaxis()->SetTitle("counts");
  h72->SetFillColor(hfillcolor);

  TH1F* h73 = new TH1F("h73"," DeltaR gamma2 Jet1 ",100,0,10);
  h73->GetXaxis()->SetTitle("#DeltaR(#gamma2-Jet1)"); 
  h73->GetYaxis()->SetTitle("counts");
  h73->SetFillColor(hfillcolor);

  TH1F* h74 = new TH1F("h74"," DeltaR gamma2 Jet2 ",100,0,10);
  h74->GetXaxis()->SetTitle("#DeltaR(#gamma2-Jet2)"); 
  h74->GetYaxis()->SetTitle("counts");
  h74->SetFillColor(hfillcolor);

//////////////////////////////////////////////////////////////////////////////////////////////////  

  TH1F* h38 = new TH1F("h38"," Pt spectrum of 1st Jet in Pt ",100,0,200);
  h38->GetXaxis()->SetTitle("Pt Jet1(#frac{GeV}{c})"); 
  h38->GetYaxis()->SetTitle("counts");
  h38->SetFillColor(hfillcolor);

  TH1F* h39 = new TH1F("h39"," Eta 1st Jet in Pt ",100,-6,6);
  h39->GetXaxis()->SetTitle("#eta Jet1"); 
  h39->GetYaxis()->SetTitle("counts");  
  h39->SetFillColor(hfillcolor);
  h39->SetMinimum(0);

  TH1F* h40 = new TH1F("h40"," Pt spectrum of 2nd Jet in Pt",100,0,200);
  h40->GetXaxis()->SetTitle("Pt Jet2(#frac{GeV}{c})"); 
  h40->GetYaxis()->SetTitle("counts");  
  h40->SetFillColor(hfillcolor);

  TH1F* h41 = new TH1F("h41"," Eta 2nd Jet in Pt ",100,-6,6);
  h41->GetXaxis()->SetTitle("#eta Jet2"); 
  h41->GetYaxis()->SetTitle("counts");
  h41->SetFillColor(hfillcolor);
  h41->SetMinimum(0);

  TH1F* h42 = new TH1F("h42"," DeltaEta 1st-2nd Jets in Pt with cuts applied",100,4,10);
  h42->GetXaxis()->SetTitle("#Delta#eta(Jet1-Jet2)"); 
  h42->GetYaxis()->SetTitle("counts");
  h42->SetFillColor(hfillcolor);
  h42->SetMinimum(0);

//////////////////////////////////////////////////////////////////////////////////////////////////  

  TH1F* h43 = new TH1F("h43"," deltaEta Higgs Jet1 ",100,-10,10);
  h43->GetXaxis()->SetTitle("#Delta#eta(Higgs-Jet1)"); 
  h43->GetYaxis()->SetTitle("counts");
  h43->SetFillColor(hfillcolor);
  h43->SetMinimum(0);

  TH1F* h44 = new TH1F("h44"," DeltaEta Higgs Jet2",100,-10,10);
  h44->GetXaxis()->SetTitle("#Delta#eta(Higgs-Jet2)"); 
  h44->GetYaxis()->SetTitle("counts");
  h44->SetFillColor(hfillcolor);
  h44->SetMinimum(0);

//////////////////////////////////////////////////////////////////////////////////////////////////  

  TH2F* h45 = new TH2F("h45"," Pt Sum of Jets 1 2 vs Pt Higgs",100,0,200,100,0,200);
  gStyle->SetPalette(1);
  //h45->Draw("COLZ");
  h45->GetXaxis()->SetTitle("Pt(Jet1+Jet2)(#frac{GeV}{c})"); 
  h45->GetYaxis()->SetTitle("Pt(Higgs)(#frac{GeV}{c})");

  TH2F* h36 = new TH2F("h36"," Eta Sum of Jets 1 2 vs Eta Higgs",100,-10,10,100,-10,10);  //100,-10,10,100,-10,10  //100,0,500,100,0,500
  gStyle->SetPalette(1);
  //h36->Draw("COLZ");
  h36->GetXaxis()->SetTitle("Eta(Jet1+Jet2)"); 
  h36->GetYaxis()->SetTitle("Eta(Higgs)");
  h36->SetMinimum(0);

  TH2F* h37 = new TH2F("h37"," Phi Sum of Jets 1 2 vs Phi Higgs",80,0,2*pi,80,0,2*pi);
  gStyle->SetPalette(1);
  //h37->Draw("COLZ");
  h37->GetXaxis()->SetTitle("Phi(Jet1+Jet2)"); 
  h37->GetYaxis()->SetTitle("Phi(Higgs) ");
  h37->SetMinimum(0);

  TH1F* h94 = new TH1F("h94"," deltaEta((Jet1+Jet2)-(Gamma1+Gamma2)) ",60,-10.,10.);
  h94->GetXaxis()->SetTitle("#delta#eta((Jet1+Jet2)-(Gamma1+Gamma2))"); 
  h94->GetYaxis()->SetTitle("counts");
  h94->SetFillColor(hfillcolor);
  h94->SetMinimum(0);

  TH2F* h95 = new TH2F("h95"," Eta Sum of Asso Jets 1 2 vs Eta Higgs",100,-10,10,100,-10,10);  //100,-10,10,100,-10,10  //100,0,500,100,0,500
  gStyle->SetPalette(1);
  //h36->Draw("COLZ");
  h95->GetXaxis()->SetTitle("Eta(Jet1+Jet2)"); 
  h95->GetYaxis()->SetTitle("Eta(Higgs)");
  h95->SetMinimum(0);

  TH1F* h96 = new TH1F("h96","Px (g1+g2+j1+j2) distribution",70,-100.,100.);
  h96->GetXaxis()->SetTitle("Px(GeV/c)"); 
  h96->GetYaxis()->SetTitle("counts");  
  h96->SetFillColor(hfillcolor);

  TH1F* h99 = new TH1F("h99","|deltaPhi((Jet1+Jet2)-(Gamma1+Gamma2))|  ",75, 0.,2*pi);   //delta peaked 180 degrees
  h99->GetXaxis()->SetTitle("#delta#phi((Jet1+Jet2)-(Gamma1+Gamma2))"); 
  h99->GetYaxis()->SetTitle("counts");
  h99->SetFillColor(hfillcolor);
  h99->SetMinimum(0);

  TH1F* h97 = new TH1F("h97","Py (g1+g2+j1+j2) distribution",70,-100.,100.);
  h97->GetXaxis()->SetTitle("Py(GeV/c)"); 
  h97->GetYaxis()->SetTitle("counts"); 
  h97->SetFillColor(hfillcolor);

  TH1F* h98 = new TH1F("h98","Pz (g1+g2+j1+j2) distribution",70,-1000.,1000.);
  h98->GetXaxis()->SetTitle("Pz(GeV/c)"); 
  h98->GetYaxis()->SetTitle("counts");
  h98->SetFillColor(hfillcolor);

//////////////////////////////////////////////////////////////////////////////////////////////////  

  TH1F* h46 = new TH1F("h46"," Eta quark1  ",100,-6,6);
  h46->GetXaxis()->SetTitle("#eta(parton1)"); 
  h46->GetYaxis()->SetTitle("counts");
  h46->SetFillColor(hfillcolor);
  h46->SetMinimum(0);

  TH1F* h47 = new TH1F("h47"," Eta quark2  ",100,-6,6);
  h47->GetXaxis()->SetTitle("#eta(parton2)"); 
  h47->GetYaxis()->SetTitle("counts");
  h47->SetFillColor(hfillcolor);
  h47->SetMinimum(0);

  TH1F* h48 = new TH1F("h48"," Pt quark1  ",100,0,200);
  h48->GetXaxis()->SetTitle("Pt(parton1)"); 
  h48->GetYaxis()->SetTitle("counts");
  h48->SetFillColor(hfillcolor);

  TH1F* h49 = new TH1F("h49"," Pt quark2  ",100,0,200);
  h49->GetXaxis()->SetTitle("Pt(parton2)"); 
  h49->GetYaxis()->SetTitle("counts");
  h49->SetFillColor(hfillcolor);

////////////////////////////////////////////////////////////////////////////////////////////////// 

  TH1F* h50 = new TH1F("h50"," DeltaR parton1 Jet ",100,0,10);
  h50->GetXaxis()->SetTitle("#DeltaR(parton1-Jet)"); 
  h50->GetYaxis()->SetTitle("counts");
  h50->SetFillColor(hfillcolor);

  TH1F* h51 = new TH1F("h51"," DeltaR parton2 Jet ",100,0,10);
  h51->GetXaxis()->SetTitle("#DeltaR(parton2-Jet)"); 
  h51->GetYaxis()->SetTitle("counts");
  h51->SetFillColor(hfillcolor);

////////////////////////////////////////////////////////////////////////////////////////////////// 

  TH2F* h52 = new TH2F("h52"," Eta Parton1 vs Eta Jet1 ",100,-6,6,100,-6,6);
  gStyle->SetPalette(1);
  //h52->Draw("COLZ");
  h52->GetXaxis()->SetTitle("#eta(Jet1)"); 
  h52->GetYaxis()->SetTitle("#eta(parton1)"); 
  h52->SetMinimum(0);

  TH2F* h53 = new TH2F("h53"," Eta Parton2 vs Eta Jet2 ",100,-6,6,100,-6,6);
  gStyle->SetPalette(1);
  //h53->Draw("COLZ");
  h53->GetXaxis()->SetTitle("#eta(Jet2)"); 
  h53->GetYaxis()->SetTitle("#eta(parton2)"); 
  h53->SetMinimum(0);

//////////////////////////////////////////////////////////////////////////////////////////////////

  TH2F* h54 = new TH2F("h54"," Phi Parton1 vs Phi Jet1 ",80,-pi,pi,80,-pi,pi);
  gStyle->SetPalette(1);
  //h54->Draw("COLZ");
  h54->GetXaxis()->SetTitle("#phi(Jet1)"); 
  h54->GetYaxis()->SetTitle("#phi(parton1)"); 
  h54->SetMinimum(0);

  TH2F* h55 = new TH2F("h55"," Phi Parton2 vs Phi Jet2 ",80,-pi,pi,80,-pi,pi);
  gStyle->SetPalette(1);
  //h55->Draw("COLZ");
  h55->GetXaxis()->SetTitle("#phi(Jet2)"); 
  h55->GetYaxis()->SetTitle("#phi(parton2)"); 
  h55->SetMinimum(0);

//////////////////////////////////////////////////////////////////////////////////////////////////

  TH1F* h56 = new TH1F("h56"," DeltaR parton1 Jet1 ",100,0,0.5);
  h56->GetXaxis()->SetTitle("#DeltaR(parton1-Jet1)"); 
  h56->GetYaxis()->SetTitle("counts"); 
  h56->SetFillColor(hfillcolor);

  TH1F* h57 = new TH1F("h57"," DeltaR parton2 Jet2 ",100,0,0.5);
  h57->GetXaxis()->SetTitle("#DeltaR(parton2-Jet2)"); 
  h57->GetYaxis()->SetTitle("counts"); 
  h57->SetFillColor(hfillcolor);

  TH2F* h58 = new TH2F("h58"," deltaR2 (parton2jet2) vs delta R1(parton1jet1)",100,0.,1.,100,0.,1.);
  gStyle->SetPalette(1);
  //h58->Draw("COLZ");
  h58->GetXaxis()->SetTitle("#DeltaR(parton1-Jet1)"); 
  h58->GetYaxis()->SetTitle("#DeltaR(parton2-Jet2)");                 

//////////////////////////////////////////////////////////////////////////////////////////////////  
/* NUOVI ISTOGRAMMI : 10 NOVEMBRE */

///////////dopo 1 CUT///////////////
  //1 fotone
  TH1F* h500 = new TH1F("h500","recoPhoton1 E",100,0,1000.);
  h500->GetXaxis()->SetTitle("E(GeV)");
  h500->GetYaxis()->SetTitle("counts");
  h500->SetFillColor(hfillcolor);

  TH1F* h501 = new TH1F("h501","recoPhoton1 Pt",100,0,140);
  h501->GetXaxis()->SetTitle("Pt(#frac{GeV}{c})");
  h501->GetYaxis()->SetTitle("counts");
  h501->SetFillColor(hfillcolor);

  TH1F* h502 = new TH1F("h502","recoPhoton1 Eta",100,-6.5,6.5);
  h502->GetXaxis()->SetTitle("#eta");
  h502->GetYaxis()->SetTitle("counts");
  h502->SetFillColor(hfillcolor);
  h502->SetMinimum(0);

  TH1F* h503 = new TH1F("h503","recoPhoton1 Phi",80,-pi,pi);
  h503->GetXaxis()->SetTitle("#phi");
  h503->GetYaxis()->SetTitle("counts");
  h503->SetFillColor(hfillcolor);
  h503->SetMinimum(0);

  //2 fotone
  TH1F* h504 = new TH1F("h504","recoPhoton2 E",100,0,1000.);
  h504->GetXaxis()->SetTitle("E(GeV)");
  h504->GetYaxis()->SetTitle("counts");
  h504->SetFillColor(hfillcolor);

  TH1F* h505 = new TH1F("h505","recoPhoton2 Pt",100,0,140);
  h505->GetXaxis()->SetTitle("Pt(#frac{GeV}{c})");
  h505->GetYaxis()->SetTitle("counts");
  h505->SetFillColor(hfillcolor);

  TH1F* h506 = new TH1F("h506","recoPhoton2 Eta",100,-6.5,6.5);
  h506->GetXaxis()->SetTitle("#eta");
  h506->GetYaxis()->SetTitle("counts");
  h506->SetFillColor(hfillcolor);
  h506->SetMinimum(0);

  TH1F* h507 = new TH1F("h507","recoPhoton2 Phi",80,-pi,pi);
  h507->GetXaxis()->SetTitle("#phi");
  h507->GetYaxis()->SetTitle("counts");
  h507->SetFillColor(hfillcolor);
  h507->SetMinimum(0);

  //1jet
  TH1F* h508 = new TH1F("h508 ","Jet1 E",100,0,1000.);                           /* RECO JETS */
  h508->GetXaxis()->SetTitle("E(GeV)");
  h508->GetYaxis()->SetTitle("counts");
  h508->SetFillColor(hfillcolor);

  TH1F* h509 = new TH1F("h509","Jet1 Pt",100,0,200);
  h509->GetXaxis()->SetTitle("Pt(#frac{GeV}{c})");
  h509->GetYaxis()->SetTitle("counts");
  h509->SetFillColor(hfillcolor);

  TH1F* h510 = new TH1F("h510","Jet1 Eta",100,-6.5,6.5);
  h510->GetXaxis()->SetTitle("#eta");
  h510->GetYaxis()->SetTitle("counts");
  h510->SetFillColor(hfillcolor);
  h510->SetMinimum(0);

  TH1F* h511 = new TH1F("h511","Jet1 Phi",80,-pi,pi);
  h511->GetXaxis()->SetTitle("#phi");
  h511->GetYaxis()->SetTitle("counts");
  h511->SetFillColor(hfillcolor);
  h511->SetMinimum(0);

  //2 jet
  TH1F* h512 = new TH1F("h512","Jet2 E",100,0,1000.);                           /* RECO JETS */
  h512->GetXaxis()->SetTitle("E(GeV)");
  h512->GetYaxis()->SetTitle("counts");
  h512->SetFillColor(hfillcolor);

  TH1F* h513 = new TH1F("h513","Jet2 Pt",100,0,200);
  h513->GetXaxis()->SetTitle("Pt(#frac{GeV}{c})");
  h513->GetYaxis()->SetTitle("counts");
  h513->SetFillColor(hfillcolor);

  TH1F* h514 = new TH1F("h514","Jet2 Eta",100,-6.5,6.5);
  h514->GetXaxis()->SetTitle("#eta");
  h514->GetYaxis()->SetTitle("counts");
  h514->SetFillColor(hfillcolor);
  h514->SetMinimum(0);

  TH1F* h515 = new TH1F("h515","Jet2 Phi",80,-pi,pi);
  h515->GetXaxis()->SetTitle("#phi");
  h515->GetYaxis()->SetTitle("counts");
  h515->SetFillColor(hfillcolor);
  h515->SetMinimum(0);

  //deltaJets (in modulo)
  TH1F* h516 = new TH1F("h516"," DeltaEta 1st-2nd Jets in Pt with cuts applied",100,0.,8.);
  h516->GetXaxis()->SetTitle("#Delta#eta(Jet1-Jet2)");
  h516->GetYaxis()->SetTitle("counts");
  h516->SetFillColor(hfillcolor);
  h516->SetMinimum(0);

  //eta(jet1)*eta(jet2)<0                                                                                                                                   
  TH1F* h517 = new TH1F("h517"," eta(jet1)*eta(jet2)<0",100,-1.5,1.5);
  h517->GetXaxis()->SetTitle("#Delta#eta(Jet1-Jet2)");
  h517->GetYaxis()->SetTitle("counts");
  h517->SetFillColor(hfillcolor);
  h517->SetMinimum(0);

  //Phi(g1+g2)-Phi(j1+j2)
  TH1F* h518 = new TH1F("h518"," RECO |deltaPhi((Jet1+Jet2)-(Gamma1+Gamma2))|  ", 75, 0.,2*pi);            
  h518->GetXaxis()->SetTitle("#delta#phi((Jet1+Jet2)-(Gamma1+Gamma2))");
  h518->GetYaxis()->SetTitle("counts");
  h518->SetFillColor(hfillcolor);
  h518->SetMinimum(0);

  //Zeppenfeld(g1+g2)
  TH1F* h519 = new TH1F("h519","Zeppenfeld reco ga1+ga2",70,-10.,10.);
  h519->GetXaxis()->SetTitle("Zeppenfeld ga1+ga2");
  h519->GetYaxis()->SetTitle("counts");
  h519->SetFillColor(hfillcolor);

  TH1F* h520 = new TH1F("h520","reco: Zeppenfeld(ga1+ga2)/|#eta(jet1)-#eta(jet2)|",70,-5.,5.);
  h520->GetXaxis()->SetTitle("Zeppenfeld(ga1+ga2)/|#eta(jet1)-#eta(jet2)|");
  h520->GetYaxis()->SetTitle("counts");
  h520->SetFillColor(hfillcolor);

  //Zeppenfeld(j1)                                                                                                                                        
  TH1F* h521 = new TH1F("h521","Zeppenfeld reco jet1",70,-6.,6.);
  h521->GetXaxis()->SetTitle("Zeppenfeld jet1");
  h521->GetYaxis()->SetTitle("counts");
  h521->SetFillColor(hfillcolor);

  TH1F* h522 = new TH1F("h522","reco: Zeppenfeld(jet1)/|#eta(jet1)-#eta(jet2)|",70,-2.,2.);
  h522->GetXaxis()->SetTitle("Zeppenfeld(jet1)/|#eta(jet1)-#eta(jet2)|");
  h522->GetYaxis()->SetTitle("counts");
  h522->SetFillColor(hfillcolor);

  //Zeppenfeld(j2)                                                                                                                                          
  TH1F* h523 = new TH1F("h523","Zeppenfeld reco jet2",70,-6.,6.);
  h523->GetXaxis()->SetTitle("Zeppenfeld jet2");
  h523->GetYaxis()->SetTitle("counts");
  h523->SetFillColor(hfillcolor);

  TH1F* h524 = new TH1F("h524","reco: Zeppenfeld(jet2)/|#eta(jet1)-#eta(jet2)|",70,-5.,5.);
  h524->GetXaxis()->SetTitle("Zeppenfeld(jet2)/|#eta(jet1)-#eta(jet2)|");
  h524->GetYaxis()->SetTitle("counts");
  h524->SetFillColor(hfillcolor);

  TH1F* h525 = new TH1F ("h525","RECO Massa Inv 2 fotoni",100,100,150);
  h525->GetXaxis()->SetTitleOffset(1.2);
  h525->GetXaxis()->SetTitle("M(#frac{GeV}{c^{2}})");
  h525->GetYaxis()->SetTitle("counts");
  h525->SetFillColor(hfillcolor);

  TH1F* h526 = new TH1F ("h526","RECO Massa Inv jets",100,0,1500);
  h526->GetXaxis()->SetTitleOffset(1.2);
  h526->GetXaxis()->SetTitle("M(j1 j2)(#frac{GeV}{c^{2}})");
  h526->GetYaxis()->SetTitle("counts");
  h526->SetFillColor(hfillcolor);

///////////dopo 2 CUT/////////////
  //1 fotone
  TH1F* h600 = new TH1F("h600","recoPhoton1 E",100,0,1000.);
  h600->GetXaxis()->SetTitle("E(GeV)");
  h600->GetYaxis()->SetTitle("counts");
  h600->SetFillColor(hfillcolor);

  TH1F* h601 = new TH1F("h601","recoPhoton1 Pt",100,0,140);
  h601->GetXaxis()->SetTitle("Pt(#frac{GeV}{c})");
  h601->GetYaxis()->SetTitle("counts");
  h601->SetFillColor(hfillcolor);

  TH1F* h602 = new TH1F("h602","recoPhoton1 Eta",100,-6.5,6.5);
  h602->GetXaxis()->SetTitle("#eta");
  h602->GetYaxis()->SetTitle("counts");
  h602->SetFillColor(hfillcolor);
  h602->SetMinimum(0);

  TH1F* h603 = new TH1F("h603","recoPhoton1 Phi",80,-pi,pi);
  h603->GetXaxis()->SetTitle("#phi");
  h603->GetYaxis()->SetTitle("counts");
  h603->SetFillColor(hfillcolor);
  h603->SetMinimum(0);

  //2 fotone
  TH1F* h604 = new TH1F("h604","recoPhoton2 E",100,0,1000.);
  h604->GetXaxis()->SetTitle("E(GeV)");
  h604->GetYaxis()->SetTitle("counts");
  h604->SetFillColor(hfillcolor);

  TH1F* h605 = new TH1F("h605","recoPhoton2 Pt",100,0,140);
  h605->GetXaxis()->SetTitle("Pt(#frac{GeV}{c})");
  h605->GetYaxis()->SetTitle("counts");
  h605->SetFillColor(hfillcolor);

  TH1F* h606 = new TH1F("h606","recoPhoton2 Eta",100,-6.5,6.5);
  h606->GetXaxis()->SetTitle("#eta");
  h606->GetYaxis()->SetTitle("counts");
  h606->SetFillColor(hfillcolor);
  h606->SetMinimum(0);

  TH1F* h607 = new TH1F("h607","recoPhoton2 Phi",80,-pi,pi);
  h607->GetXaxis()->SetTitle("#phi");
  h607->GetYaxis()->SetTitle("counts");
  h607->SetFillColor(hfillcolor);
  h607->SetMinimum(0);

  //1jet
  TH1F* h608 = new TH1F("h608 ","Jet1 E",100,0,1000.);                           /* RECO JETS */
  h608->GetXaxis()->SetTitle("E(GeV)");
  h608->GetYaxis()->SetTitle("counts");
  h608->SetFillColor(hfillcolor);

  TH1F* h609 = new TH1F("h609","Jet1 Pt",100,0,200);
  h609->GetXaxis()->SetTitle("Pt(#frac{GeV}{c})");
  h609->GetYaxis()->SetTitle("counts");
  h609->SetFillColor(hfillcolor);

  TH1F* h610 = new TH1F("h610","Jet1 Eta",100,-6.5,6.5);
  h610->GetXaxis()->SetTitle("#eta");
  h610->GetYaxis()->SetTitle("counts");
  h610->SetFillColor(hfillcolor);
  h610->SetMinimum(0);

  TH1F* h611 = new TH1F("h611","Jet1 Phi",80,-pi,pi);
  h611->GetXaxis()->SetTitle("#phi");
  h611->GetYaxis()->SetTitle("counts");
  h611->SetFillColor(hfillcolor);
  h611->SetMinimum(0);

  //2 jet
  TH1F* h612 = new TH1F("h612","Jet2 E",100,0,1000.);                           /* RECO JETS */
  h612->GetXaxis()->SetTitle("E(GeV)");
  h612->GetYaxis()->SetTitle("counts");
  h612->SetFillColor(hfillcolor);

  TH1F* h613 = new TH1F("h613","Jet2 Pt",100,0,200);
  h613->GetXaxis()->SetTitle("Pt(#frac{GeV}{c})");
  h613->GetYaxis()->SetTitle("counts");
  h613->SetFillColor(hfillcolor);

  TH1F* h614 = new TH1F("h614","Jet2 Eta",100,-6.5,6.5);
  h614->GetXaxis()->SetTitle("#eta");
  h614->GetYaxis()->SetTitle("counts");
  h614->SetFillColor(hfillcolor);
  h614->SetMinimum(0);

  TH1F* h615 = new TH1F("h615","Jet2 Phi",80,-pi,pi);
  h615->GetXaxis()->SetTitle("#phi");
  h615->GetYaxis()->SetTitle("counts");
  h615->SetFillColor(hfillcolor);
  h615->SetMinimum(0);

  //deltaJets (in modulo)
  TH1F* h616 = new TH1F("h616"," DeltaEta 1st-2nd Jets in Pt with cuts applied",100,0.,8.);
  h616->GetXaxis()->SetTitle("#Delta#eta(Jet1-Jet2)");
  h616->GetYaxis()->SetTitle("counts");
  h616->SetFillColor(hfillcolor);
  h616->SetMinimum(0);

  //eta(jet1)*eta(jet2)<0                                                                                                                                   
  TH1F* h617 = new TH1F("h617"," eta(jet1)*eta(jet2)<0",100,-1.5,1.5);
  h617->GetXaxis()->SetTitle("#Delta#eta(Jet1-Jet2)");
  h617->GetYaxis()->SetTitle("counts");
  h617->SetFillColor(hfillcolor);
  h617->SetMinimum(0);

  //Phi(g1+g2)-Phi(j1+j2)
  TH1F* h618 = new TH1F("h618"," RECO |deltaPhi((Jet1+Jet2)-(Gamma1+Gamma2))|  ", 75, 0.,2*pi);            
  h618->GetXaxis()->SetTitle("#delta#phi((Jet1+Jet2)-(Gamma1+Gamma2))");
  h618->GetYaxis()->SetTitle("counts");
  h618->SetFillColor(hfillcolor);
  h618->SetMinimum(0);

  //Zeppenfeld(g1+g2)
  TH1F* h619 = new TH1F("h619","Zeppenfeld reco ga1+ga2",70,-10.,10.);
  h619->GetXaxis()->SetTitle("Zeppenfeld ga1+ga2");
  h619->GetYaxis()->SetTitle("counts");
  h619->SetFillColor(hfillcolor);

  TH1F* h620 = new TH1F("h620","reco: Zeppenfeld(ga1+ga2)/|#eta(jet1)-#eta(jet2)|",70,-5.,5.);
  h620->GetXaxis()->SetTitle("Zeppenfeld(ga1+ga2)/|#eta(jet1)-#eta(jet2)|");
  h620->GetYaxis()->SetTitle("counts");
  h620->SetFillColor(hfillcolor);

  //Zeppenfeld(j1)                                                                                                                                        
  TH1F* h621 = new TH1F("h621","Zeppenfeld reco jet1",70,-6.,6.);
  h621->GetXaxis()->SetTitle("Zeppenfeld jet1");
  h621->GetYaxis()->SetTitle("counts");
  h621->SetFillColor(hfillcolor);

  TH1F* h622 = new TH1F("h622","reco: Zeppenfeld(jet1)/|#eta(jet1)-#eta(jet2)|",70,-2.,2.);
  h622->GetXaxis()->SetTitle("Zeppenfeld(jet1)/|#eta(jet1)-#eta(jet2)|");
  h622->GetYaxis()->SetTitle("counts");
  h622->SetFillColor(hfillcolor);

  //Zeppenfeld(j2)                                                                                                                                          
  TH1F* h623 = new TH1F("h623","Zeppenfeld reco jet2",70,-6.,6.);
  h623->GetXaxis()->SetTitle("Zeppenfeld jet2");
  h623->GetYaxis()->SetTitle("counts");
  h623->SetFillColor(hfillcolor);

  TH1F* h624 = new TH1F("h624","reco: Zeppenfeld(jet2)/|#eta(jet1)-#eta(jet2)|",70,-2.,2.);
  h624->GetXaxis()->SetTitle("Zeppenfeld(jet2)/|#eta(jet1)-#eta(jet2)|");
  h624->GetYaxis()->SetTitle("counts");
  h624->SetFillColor(hfillcolor);

  TH1F* h625 = new TH1F ("h625","RECO Massa Inv 2 fotoni",100,100,150);
  h625->GetXaxis()->SetTitleOffset(1.2);
  h625->GetXaxis()->SetTitle("M(#frac{GeV}{c^{2}})");
  h625->GetYaxis()->SetTitle("counts");
  h625->SetFillColor(hfillcolor);

  TH1F* h626 = new TH1F ("h626","RECO Massa Inv jets",100,0,1500);
  h626->GetXaxis()->SetTitleOffset(1.2);
  h626->GetXaxis()->SetTitle("M(j1 j2)(#frac{GeV}{c^{2}})");
  h626->GetYaxis()->SetTitle("counts");
  h626->SetFillColor(hfillcolor);


///////////dopo 3 CUT/////////////
  //1 fotone
  TH1F* h700 = new TH1F("h700","recoPhoton1 E",100,0,1000.);
  h700->GetXaxis()->SetTitle("E(GeV)");
  h700->GetYaxis()->SetTitle("counts");
  h700->SetFillColor(hfillcolor);

  TH1F* h701 = new TH1F("h701","recoPhoton1 Pt",100,0,140);
  h701->GetXaxis()->SetTitle("Pt(#frac{GeV}{c})");
  h701->GetYaxis()->SetTitle("counts");
  h701->SetFillColor(hfillcolor);

  TH1F* h702 = new TH1F("h702","recoPhoton1 Eta",100,-6.5,6.5);
  h702->GetXaxis()->SetTitle("#eta");
  h702->GetYaxis()->SetTitle("counts");
  h702->SetFillColor(hfillcolor);
  h702->SetMinimum(0);

  TH1F* h703 = new TH1F("h703","recoPhoton1 Phi",80,-pi,pi);
  h703->GetXaxis()->SetTitle("#phi");
  h703->GetYaxis()->SetTitle("counts");
  h703->SetFillColor(hfillcolor);
  h703->SetMinimum(0);

  //2 fotone
  TH1F* h704 = new TH1F("h704","recoPhoton2 E",100,0,1000.);
  h704->GetXaxis()->SetTitle("E(GeV)");
  h704->GetYaxis()->SetTitle("counts");
  h704->SetFillColor(hfillcolor);

  TH1F* h705 = new TH1F("h705","recoPhoton2 Pt",100,0,140);
  h705->GetXaxis()->SetTitle("Pt(#frac{GeV}{c})");
  h705->GetYaxis()->SetTitle("counts");
  h705->SetFillColor(hfillcolor);

  TH1F* h706 = new TH1F("h706","recoPhoton2 Eta",100,-6.5,6.5);
  h706->GetXaxis()->SetTitle("#eta");
  h706->GetYaxis()->SetTitle("counts");
  h706->SetFillColor(hfillcolor);
  h706->SetMinimum(0);

  TH1F* h707 = new TH1F("h707","recoPhoton2 Phi",80,-pi,pi);
  h707->GetXaxis()->SetTitle("#phi");
  h707->GetYaxis()->SetTitle("counts");
  h707->SetFillColor(hfillcolor);
  h707->SetMinimum(0);

  //1jet
  TH1F* h708 = new TH1F("h708 ","Jet1 E",100,0,1000.);                           /* RECO JETS */
  h708->GetXaxis()->SetTitle("E(GeV)");
  h708->GetYaxis()->SetTitle("counts");
  h708->SetFillColor(hfillcolor);

  TH1F* h709 = new TH1F("h709","Jet1 Pt",100,0,200);
  h709->GetXaxis()->SetTitle("Pt(#frac{GeV}{c})");
  h709->GetYaxis()->SetTitle("counts");
  h709->SetFillColor(hfillcolor);

  TH1F* h710 = new TH1F("h710","Jet1 Eta",100,-6.5,6.5);
  h710->GetXaxis()->SetTitle("#eta");
  h710->GetYaxis()->SetTitle("counts");
  h710->SetFillColor(hfillcolor);
  h710->SetMinimum(0);

  TH1F* h711 = new TH1F("h711","Jet1 Phi",80,-pi,pi);
  h711->GetXaxis()->SetTitle("#phi");
  h711->GetYaxis()->SetTitle("counts");
  h711->SetFillColor(hfillcolor);
  h711->SetMinimum(0);

  //2 jet
  TH1F* h712 = new TH1F("h712","Jet2 E",100,0,1000.);                           /* RECO JETS */
  h712->GetXaxis()->SetTitle("E(GeV)");
  h712->GetYaxis()->SetTitle("counts");
  h712->SetFillColor(hfillcolor);

  TH1F* h713 = new TH1F("h713","Jet2 Pt",100,0,200);
  h713->GetXaxis()->SetTitle("Pt(#frac{GeV}{c})");
  h713->GetYaxis()->SetTitle("counts");
  h713->SetFillColor(hfillcolor);

  TH1F* h714 = new TH1F("h714","Jet2 Eta",100,-6.5,6.5);
  h714->GetXaxis()->SetTitle("#eta");
  h714->GetYaxis()->SetTitle("counts");
  h714->SetFillColor(hfillcolor);
  h714->SetMinimum(0);

  TH1F* h715 = new TH1F("h715","Jet2 Phi",80,-pi,pi);
  h715->GetXaxis()->SetTitle("#phi");
  h715->GetYaxis()->SetTitle("counts");
  h715->SetFillColor(hfillcolor);
  h715->SetMinimum(0);

  //deltaJets (in modulo)
  TH1F* h716 = new TH1F("h716"," DeltaEta 1st-2nd Jets in Pt with cuts applied",100,0.,8.);
  h716->GetXaxis()->SetTitle("#Delta#eta(Jet1-Jet2)");
  h716->GetYaxis()->SetTitle("counts");
  h716->SetFillColor(hfillcolor);
  h716->SetMinimum(0);

  //eta(jet1)*eta(jet2)<0                                                                                                                                   
  TH1F* h717 = new TH1F("h717"," eta(jet1)*eta(jet2)<0",100,-1.5,1.5);
  h717->GetXaxis()->SetTitle("#Delta#eta(Jet1-Jet2)");
  h717->GetYaxis()->SetTitle("counts");
  h717->SetFillColor(hfillcolor);
  h717->SetMinimum(0);

  //Phi(g1+g2)-Phi(j1+j2)
  TH1F* h718 = new TH1F("h718"," RECO |deltaPhi((Jet1+Jet2)-(Gamma1+Gamma2))|  ", 75, 0.,2*pi);            
  h718->GetXaxis()->SetTitle("#delta#phi((Jet1+Jet2)-(Gamma1+Gamma2))");
  h718->GetYaxis()->SetTitle("counts");
  h718->SetFillColor(hfillcolor);
  h718->SetMinimum(0);

  //Zeppenfeld(g1+g2)
  TH1F* h719 = new TH1F("h719","Zeppenfeld reco ga1+ga2",70,-10.,10.);
  h719->GetXaxis()->SetTitle("Zeppenfeld ga1+ga2");
  h719->GetYaxis()->SetTitle("counts");
  h719->SetFillColor(hfillcolor);

  TH1F* h720 = new TH1F("h720","reco: Zeppenfeld(ga1+ga2)/|#eta(jet1)-#eta(jet2)|",70,-5.,5.);
  h720->GetXaxis()->SetTitle("Zeppenfeld(ga1+ga2)/|#eta(jet1)-#eta(jet2)|");
  h720->GetYaxis()->SetTitle("counts");
  h720->SetFillColor(hfillcolor);

  //Zeppenfeld(j1)                                                                                                                                        
  TH1F* h721 = new TH1F("h721","Zeppenfeld reco jet1",70,-6.,6.);
  h721->GetXaxis()->SetTitle("Zeppenfeld jet1");
  h721->GetYaxis()->SetTitle("counts");
  h721->SetFillColor(hfillcolor);

  TH1F* h722 = new TH1F("h722","reco: Zeppenfeld(jet1)/|#eta(jet1)-#eta(jet2)|",70,-2.,2.);
  h722->GetXaxis()->SetTitle("Zeppenfeld(jet1)/|#eta(jet1)-#eta(jet2)|");
  h722->GetYaxis()->SetTitle("counts");
  h722->SetFillColor(hfillcolor);

  //Zeppenfeld(j2)                                                                                                                                          
  TH1F* h723 = new TH1F("h723","Zeppenfeld reco jet2",70,-6.,6.);
  h723->GetXaxis()->SetTitle("Zeppenfeld jet2");
  h723->GetYaxis()->SetTitle("counts");
  h723->SetFillColor(hfillcolor);

  TH1F* h724 = new TH1F("h724","reco: Zeppenfeld(jet2)/|#eta(jet1)-#eta(jet2)|",70,-2.,2.);
  h724->GetXaxis()->SetTitle("Zeppenfeld(jet2)/|#eta(jet1)-#eta(jet2)|");
  h724->GetYaxis()->SetTitle("counts");
  h724->SetFillColor(hfillcolor);

  TH1F* h725 = new TH1F ("h725","RECO Massa Inv 2 fotoni",100,100,150);
  h725->GetXaxis()->SetTitleOffset(1.2);
  h725->GetXaxis()->SetTitle("M(#frac{GeV}{c^{2}})");
  h725->GetYaxis()->SetTitle("counts");
  h725->SetFillColor(hfillcolor);

  TH1F* h726 = new TH1F ("h726","RECO Massa Inv jets",100,0,1500);
  h726->GetXaxis()->SetTitleOffset(1.2);
  h726->GetXaxis()->SetTitle("M(j1 j2)(#frac{GeV}{c^{2}})");
  h726->GetYaxis()->SetTitle("counts");
  h726->SetFillColor(hfillcolor);
////////////dopo 4 CUT///////////////////
  //1 fotone
  TH1F* h800 = new TH1F("h800","recoPhoton1 E",100,0,1000.);
  h800->GetXaxis()->SetTitle("E(GeV)");
  h800->GetYaxis()->SetTitle("counts");
  h800->SetFillColor(hfillcolor);

  TH1F* h801 = new TH1F("h801","recoPhoton1 Pt",100,0,140);
  h801->GetXaxis()->SetTitle("Pt(#frac{GeV}{c})");
  h801->GetYaxis()->SetTitle("counts");
  h801->SetFillColor(hfillcolor);

  TH1F* h802 = new TH1F("h802","recoPhoton1 Eta",100,-6.5,6.5);
  h802->GetXaxis()->SetTitle("#eta");
  h802->GetYaxis()->SetTitle("counts");
  h802->SetFillColor(hfillcolor);
  h802->SetMinimum(0);

  TH1F* h803 = new TH1F("h803","recoPhoton1 Phi",80,-pi,pi);
  h803->GetXaxis()->SetTitle("#phi");
  h803->GetYaxis()->SetTitle("counts");
  h803->SetFillColor(hfillcolor);
  h803->SetMinimum(0);

  //2 fotone
  TH1F* h804 = new TH1F("h804","recoPhoton2 E",100,0,1000.);
  h804->GetXaxis()->SetTitle("E(GeV)");
  h804->GetYaxis()->SetTitle("counts");
  h804->SetFillColor(hfillcolor);

  TH1F* h805 = new TH1F("h805","recoPhoton2 Pt",100,0,140);
  h805->GetXaxis()->SetTitle("Pt(#frac{GeV}{c})");
  h805->GetYaxis()->SetTitle("counts");
  h805->SetFillColor(hfillcolor);

  TH1F* h806 = new TH1F("h806","recoPhoton2 Eta",100,-6.5,6.5);
  h806->GetXaxis()->SetTitle("#eta");
  h806->GetYaxis()->SetTitle("counts");
  h806->SetFillColor(hfillcolor);
  h806->SetMinimum(0);

  TH1F* h807 = new TH1F("h807","recoPhoton2 Phi",80,-pi,pi);
  h807->GetXaxis()->SetTitle("#phi");
  h807->GetYaxis()->SetTitle("counts");
  h807->SetFillColor(hfillcolor);
  h807->SetMinimum(0);

  //1jet
  TH1F* h808 = new TH1F("h808 ","Jet1 E",100,0,1000.);                           /* RECO JETS */
  h808->GetXaxis()->SetTitle("E(GeV)");
  h808->GetYaxis()->SetTitle("counts");
  h808->SetFillColor(hfillcolor);

  TH1F* h809 = new TH1F("h809","Jet1 Pt",100,0,200);
  h809->GetXaxis()->SetTitle("Pt(#frac{GeV}{c})");
  h809->GetYaxis()->SetTitle("counts");
  h809->SetFillColor(hfillcolor);

  TH1F* h810 = new TH1F("h810","Jet1 Eta",100,-6.5,6.5);
  h810->GetXaxis()->SetTitle("#eta");
  h810->GetYaxis()->SetTitle("counts");
  h810->SetFillColor(hfillcolor);
  h810->SetMinimum(0);

  TH1F* h811 = new TH1F("h811","Jet1 Phi",80,-pi,pi);
  h811->GetXaxis()->SetTitle("#phi");
  h811->GetYaxis()->SetTitle("counts");
  h811->SetFillColor(hfillcolor);
  h811->SetMinimum(0);

  //2 jet
  TH1F* h812 = new TH1F("h812","Jet2 E",100,0,1000.);                           /* RECO JETS */
  h812->GetXaxis()->SetTitle("E(GeV)");
  h812->GetYaxis()->SetTitle("counts");
  h812->SetFillColor(hfillcolor);

  TH1F* h813 = new TH1F("h813","Jet2 Pt",100,0,200);
  h813->GetXaxis()->SetTitle("Pt(#frac{GeV}{c})");
  h813->GetYaxis()->SetTitle("counts");
  h813->SetFillColor(hfillcolor);

  TH1F* h814 = new TH1F("h814","Jet2 Eta",100,-6.5,6.5);
  h814->GetXaxis()->SetTitle("#eta");
  h814->GetYaxis()->SetTitle("counts");
  h814->SetFillColor(hfillcolor);
  h814->SetMinimum(0);

  TH1F* h815 = new TH1F("h815","Jet2 Phi",80,-pi,pi);
  h815->GetXaxis()->SetTitle("#phi");
  h815->GetYaxis()->SetTitle("counts");
  h815->SetFillColor(hfillcolor);
  h815->SetMinimum(0);

  //deltaJets (in modulo)
  TH1F* h816 = new TH1F("h816"," DeltaEta 1st-2nd Jets in Pt with cuts applied",100,0.,8.);
  h816->GetXaxis()->SetTitle("#Delta#eta(Jet1-Jet2)");
  h816->GetYaxis()->SetTitle("counts");
  h816->SetFillColor(hfillcolor);
  h816->SetMinimum(0);

  //eta(jet1)*eta(jet2)<0                                                                                                                                   
  TH1F* h817 = new TH1F("h817"," eta(jet1)*eta(jet2)<0",100,-1.5,1.5);
  h817->GetXaxis()->SetTitle("#Delta#eta(Jet1-Jet2)");
  h817->GetYaxis()->SetTitle("counts");
  h817->SetFillColor(hfillcolor);
  h817->SetMinimum(0);

  //Phi(g1+g2)-Phi(j1+j2)
  TH1F* h818 = new TH1F("h818"," RECO |deltaPhi((Jet1+Jet2)-(Gamma1+Gamma2))|  ", 75, 0.,2*pi);            
  h818->GetXaxis()->SetTitle("#delta#phi((Jet1+Jet2)-(Gamma1+Gamma2))");
  h818->GetYaxis()->SetTitle("counts");
  h818->SetFillColor(hfillcolor);
  h818->SetMinimum(0);

  //Zeppenfeld(g1+g2)
  TH1F* h819 = new TH1F("h819","Zeppenfeld reco ga1+ga2",70,-10.,10.);
  h819->GetXaxis()->SetTitle("Zeppenfeld ga1+ga2");
  h819->GetYaxis()->SetTitle("counts");
  h819->SetFillColor(hfillcolor);

  TH1F* h820 = new TH1F("h820","reco: Zeppenfeld(ga1+ga2)/|#eta(jet1)-#eta(jet2)|",70,-5.,5.);
  h820->GetXaxis()->SetTitle("Zeppenfeld(ga1+ga2)/|#eta(jet1)-#eta(jet2)|");
  h820->GetYaxis()->SetTitle("counts");
  h820->SetFillColor(hfillcolor);

  //Zeppenfeld(j1)                                                                                                                                        
  TH1F* h821 = new TH1F("h821","Zeppenfeld reco jet1",70,-6.,6.);
  h821->GetXaxis()->SetTitle("Zeppenfeld jet1");
  h821->GetYaxis()->SetTitle("counts");
  h821->SetFillColor(hfillcolor);

  TH1F* h822 = new TH1F("h822","reco: Zeppenfeld(jet1)/|#eta(jet1)-#eta(jet2)|",70,-2.,2.);
  h822->GetXaxis()->SetTitle("Zeppenfeld(jet1)/|#eta(jet1)-#eta(jet2)|");
  h822->GetYaxis()->SetTitle("counts");
  h822->SetFillColor(hfillcolor);

  //Zeppenfeld(j2)                                                                                                                                          
  TH1F* h823 = new TH1F("h823","Zeppenfeld reco jet2",70,-6.,6.);
  h823->GetXaxis()->SetTitle("Zeppenfeld jet2");
  h823->GetYaxis()->SetTitle("counts");
  h823->SetFillColor(hfillcolor);

  TH1F* h824 = new TH1F("h824","reco: Zeppenfeld(jet2)/|#eta(jet1)-#eta(jet2)|",70,-2.,2.);
  h824->GetXaxis()->SetTitle("Zeppenfeld(jet2)/|#eta(jet1)-#eta(jet2)|");
  h824->GetYaxis()->SetTitle("counts");
  h824->SetFillColor(hfillcolor);

  TH1F* h825 = new TH1F ("h825","RECO Massa Inv 2 fotoni",100,100,150);
  h825->GetXaxis()->SetTitleOffset(1.2);
  h825->GetXaxis()->SetTitle("M(#frac{GeV}{c^{2}})");
  h825->GetYaxis()->SetTitle("counts");
  h825->SetFillColor(hfillcolor);

  TH1F* h826 = new TH1F ("h826","RECO Massa Inv jets",100,0,1500);
  h826->GetXaxis()->SetTitleOffset(1.2);
  h826->GetXaxis()->SetTitle("M(j1 j2)(#frac{GeV}{c^{2}})");
  h826->GetYaxis()->SetTitle("counts");
  h826->SetFillColor(hfillcolor);

//extra histogram

  TH2F* h827=new TH2F("h827","M higgs vs DeltaPhi",75,0.,2*pi,100,100,150);
  gStyle->SetPalette(1);
  h827->Draw("COLZ");
  h827->GetXaxis()->SetTitle("#Delta#phi");
  h827->GetYaxis()->SetTitle("M#gamma#gamma"); 
  

  ////////////dopo 5 CUT///////////////////                                                                                                                    
  //1 fotone                                                                                                                                                 
  TH1F* h900 = new TH1F("h900","recoPhoton1 E",100,0,1000.);
  h900->GetXaxis()->SetTitle("E(GeV)");
  h900->GetYaxis()->SetTitle("counts");
  h900->SetFillColor(hfillcolor);

  TH1F* h901 = new TH1F("h901","recoPhoton1 Pt",100,0,140);
  h901->GetXaxis()->SetTitle("Pt(#frac{GeV}{c})");
  h901->GetYaxis()->SetTitle("counts");
  h901->SetFillColor(hfillcolor);

  TH1F* h902 = new TH1F("h902","recoPhoton1 Eta",100,-6.5,6.5);
  h902->GetXaxis()->SetTitle("#eta");
  h902->GetYaxis()->SetTitle("counts");
  h902->SetFillColor(hfillcolor);
  h902->SetMinimum(0);

  TH1F* h903 = new TH1F("h903","recoPhoton1 Phi",80,-pi,pi);
  h903->GetXaxis()->SetTitle("#phi");
  h903->GetYaxis()->SetTitle("counts");
  h903->SetFillColor(hfillcolor);
  h903->SetMinimum(0);

  //2 fotone                                                                                                                                                 
  TH1F* h904 = new TH1F("h904","recoPhoton2 E",100,0,1000.);
  h904->GetXaxis()->SetTitle("E(GeV)");
  h904->GetYaxis()->SetTitle("counts");
  h904->SetFillColor(hfillcolor);

  TH1F* h905 = new TH1F("h905","recoPhoton2 Pt",100,0,140);
  h905->GetXaxis()->SetTitle("Pt(#frac{GeV}{c})");
  h905->GetYaxis()->SetTitle("counts");
  h905->SetFillColor(hfillcolor);

  TH1F* h906 = new TH1F("h906","recoPhoton2 Eta",100,-6.5,6.5);
  h906->GetXaxis()->SetTitle("#eta");
  h906->GetYaxis()->SetTitle("counts");
  h906->SetFillColor(hfillcolor);
  h906->SetMinimum(0);

  TH1F* h907 = new TH1F("h907","recoPhoton2 Phi",80,-pi,pi);
  h907->GetXaxis()->SetTitle("#phi");
  h907->GetYaxis()->SetTitle("counts");
  h907->SetFillColor(hfillcolor);
  h907->SetMinimum(0);

  //1jet                                                                                                                                                     
  TH1F* h908 = new TH1F("h908 ","Jet1 E",100,0,1000.);                           /* RECO JETS */
  h908->GetXaxis()->SetTitle("E(GeV)");
  h908->GetYaxis()->SetTitle("counts");
  h908->SetFillColor(hfillcolor);

  TH1F* h909 = new TH1F("h909","Jet1 Pt",100,0,200);
  h909->GetXaxis()->SetTitle("Pt(#frac{GeV}{c})");
  h909->GetYaxis()->SetTitle("counts");
  h909->SetFillColor(hfillcolor);

  TH1F* h910 = new TH1F("h910","Jet1 Eta",100,-6.5,6.5);
  h910->GetXaxis()->SetTitle("#eta");
  h910->GetYaxis()->SetTitle("counts");
  h910->SetFillColor(hfillcolor);
  h910->SetMinimum(0);

  TH1F* h911 = new TH1F("h911","Jet1 Phi",80,-pi,pi);
  h911->GetXaxis()->SetTitle("#phi");
  h911->GetYaxis()->SetTitle("counts");
  h911->SetFillColor(hfillcolor);
  h911->SetMinimum(0);

  //2jet                                                                                                                                                    
  TH1F* h912 = new TH1F("h912","Jet2 E",100,0,1000.);                           /* RECO JETS */
  h912->GetXaxis()->SetTitle("E(GeV)");
  h912->GetYaxis()->SetTitle("counts");
  h912->SetFillColor(hfillcolor);

  TH1F* h913 = new TH1F("h913","Jet2 Pt",100,0,200);
  h913->GetXaxis()->SetTitle("Pt(#frac{GeV}{c})");
  h913->GetYaxis()->SetTitle("counts");
  h913->SetFillColor(hfillcolor);

  TH1F* h914 = new TH1F("h914","Jet2 Eta",100,-6.5,6.5);
  h914->GetXaxis()->SetTitle("#eta");
  h914->GetYaxis()->SetTitle("counts");
  h914->SetFillColor(hfillcolor);
  h914->SetMinimum(0);

  TH1F* h915 = new TH1F("h915","Jet2 Phi",80,-pi,pi);
  h915->GetXaxis()->SetTitle("#phi");
  h915->GetYaxis()->SetTitle("counts");
  h915->SetFillColor(hfillcolor);
  h915->SetMinimum(0);

  //deltaJets (in modulo)                                                                                                                                    
  TH1F* h916 = new TH1F("h916"," DeltaEta 1st-2nd Jets in Pt with cuts applied",100,0.,8.);
  h916->GetXaxis()->SetTitle("#Delta#eta(Jet1-Jet2)");
  h916->GetYaxis()->SetTitle("counts");
  h916->SetFillColor(hfillcolor);
  h916->SetMinimum(0);

  //eta(jet1)*eta(jet2)<0                                                                                                                                   
  TH1F* h917 = new TH1F("h917"," eta(jet1)*eta(jet2)<0",100,-1.5,1.5);
  h917->GetXaxis()->SetTitle("#Delta#eta(Jet1-Jet2)");
  h917->GetYaxis()->SetTitle("counts");
  h917->SetFillColor(hfillcolor);
  h917->SetMinimum(0);

  //Phi(g1+g2)-Phi(j1+j2)                                                                                                                                    
  TH1F* h918 = new TH1F("h918"," RECO |deltaPhi((Jet1+Jet2)-(Gamma1+Gamma2))|  ", 75, 0.,2*pi);                                                           
  h918->GetXaxis()->SetTitle("#delta#phi((Jet1+Jet2)-(Gamma1+Gamma2))");
  h918->GetYaxis()->SetTitle("counts");
  h918->SetFillColor(hfillcolor);
  h918->SetMinimum(0);

  //Zeppenfeld(g1+g2)                                                                                                                                        
  TH1F* h919 = new TH1F("h919","Zeppenfeld reco ga1+ga2",70,-10.,10.);
  h919->GetXaxis()->SetTitle("Zeppenfeld ga1+ga2");
  h919->GetYaxis()->SetTitle("counts");
  h919->SetFillColor(hfillcolor);
  
  TH1F* h920 = new TH1F("h920","reco: Zeppenfeld(ga1+ga2)/|#eta(jet1)-#eta(jet2)|",70,-5.,5.);
  h920->GetXaxis()->SetTitle("Zeppenfeld(ga1+ga2)/|#eta(jet1)-#eta(jet2)|");
  h920->GetYaxis()->SetTitle("counts");
  h920->SetFillColor(hfillcolor);
  
  //Zeppenfeld(j1)                                                                                                                                           
  TH1F* h921 = new TH1F("h921","Zeppenfeld reco jet1",70,-6.,6.);
  h921->GetXaxis()->SetTitle("Zeppenfeld jet1");
  h921->GetYaxis()->SetTitle("counts");
  h921->SetFillColor(hfillcolor);
  
  TH1F* h922 = new TH1F("h922","reco: Zeppenfeld(jet1)/|#eta(jet1)-#eta(jet2)|",70,-2.,2.);
  h922->GetXaxis()->SetTitle("Zeppenfeld(jet1)/|#eta(jet1)-#eta(jet2)|");
  h922->GetYaxis()->SetTitle("counts");
  h922->SetFillColor(hfillcolor);
  
  //Zeppenfeld(j2)                                                                                                                                          
  TH1F* h923 = new TH1F("h923","Zeppenfeld reco jet2",70,-6.,6.);
  h923->GetXaxis()->SetTitle("Zeppenfeld jet2");
  h923->GetYaxis()->SetTitle("counts");
  h923->SetFillColor(hfillcolor);
  
  TH1F* h924 = new TH1F("h924","reco: Zeppenfeld(jet2)/|#eta(jet1)-#eta(jet2)|",70,-2.,2.);
  h924->GetXaxis()->SetTitle("Zeppenfeld(jet2)/|#eta(jet1)-#eta(jet2)|");
  h924->GetYaxis()->SetTitle("counts");
  h924->SetFillColor(hfillcolor);
  
  TH1F* h925 = new TH1F ("h925","RECO Massa Inv 2 fotoni",100,100,150);
  h925->GetXaxis()->SetTitleOffset(1.2);
  h925->GetXaxis()->SetTitle("M(#frac{GeV}{c^{2}})");
  h925->GetYaxis()->SetTitle("counts");
  h925->SetFillColor(hfillcolor);

  TH1F* h927 = new TH1F ("h927","RECO Massa Inv 2 fotoni",1000,0.,1000.);
  h927->GetXaxis()->SetTitleOffset(1.2);
  h927->GetXaxis()->SetTitle("M(#frac{GeV}{c^{2}})");
  h927->GetYaxis()->SetTitle("counts");
  h927->SetFillColor(hfillcolor);



  TH1F* h926 = new TH1F ("h926","RECO Massa Inv jets",100,0,1500);
  h926->GetXaxis()->SetTitleOffset(1.2);
  h926->GetXaxis()->SetTitle("M(j1 j2)(#frac{GeV}{c^{2}})");
  h926->GetYaxis()->SetTitle("counts");
  h926->SetFillColor(hfillcolor);

  ///////////dopo6CUT///////////////////                                                                                                                    
  //1fotone                                                                                                                                                 
  TH1F* h1000 = new TH1F("h1000","recoPhoton1 E",100,0,1000.);
  h1000->GetXaxis()->SetTitle("E(GeV)");
  h1000->GetYaxis()->SetTitle("counts");
  h1000->SetFillColor(hfillcolor);

  TH1F* h1001 = new TH1F("h1001","recoPhoton1 Pt",100,0,140);
  h1001->GetXaxis()->SetTitle("Pt(#frac{GeV}{c})");
  h1001->GetYaxis()->SetTitle("counts");
  h1001->SetFillColor(hfillcolor);

  TH1F* h1002 = new TH1F("h1002","recoPhoton1 Eta",100,-6.5,6.5);
  h1002->GetXaxis()->SetTitle("#eta");
  h1002->GetYaxis()->SetTitle("counts");
  h1002->SetFillColor(hfillcolor);
  h1002->SetMinimum(0);

  TH1F* h1003 = new TH1F("h1003","recoPhoton1 Phi",80,-pi,pi);
  h1003->GetXaxis()->SetTitle("#phi");
  h1003->GetYaxis()->SetTitle("counts");
  h1003->SetFillColor(hfillcolor);
  h1003->SetMinimum(0);

  //2 fotone                                                                                                                                                 
  TH1F* h1004 = new TH1F("h1004","recoPhoton2 E",100,0,1000.);
  h1004->GetXaxis()->SetTitle("E(GeV)");
  h1004->GetYaxis()->SetTitle("counts");
  h1004->SetFillColor(hfillcolor);

  TH1F* h1005 = new TH1F("h1005","recoPhoton2 Pt",100,0,140);
  h1005->GetXaxis()->SetTitle("Pt(#frac{GeV}{c})");
  h1005->GetYaxis()->SetTitle("counts");
  h1005->SetFillColor(hfillcolor);

  TH1F* h1006 = new TH1F("h1006","recoPhoton2 Eta",100,-6.5,6.5);
  h1006->GetXaxis()->SetTitle("#eta");
  h1006->GetYaxis()->SetTitle("counts");
  h1006->SetFillColor(hfillcolor);
  h1006->SetMinimum(0);

  TH1F* h1007 = new TH1F("h1007","recoPhoton2 Phi",80,-pi,pi);
  h1007->GetXaxis()->SetTitle("#phi");
  h1007->GetYaxis()->SetTitle("counts");
  h1007->SetFillColor(hfillcolor);
  h1007->SetMinimum(0);

  //1jet                                                                                                                                                     
  TH1F* h1008 = new TH1F("h1008 ","Jet1 E",100,0,1000.);                           /* RECO JETS */
  h1008->GetXaxis()->SetTitle("E(GeV)");
  h1008->GetYaxis()->SetTitle("counts");
  h1008->SetFillColor(hfillcolor);

  TH1F* h1009 = new TH1F("h1009","Jet1 Pt",100,0,200);
  h1009->GetXaxis()->SetTitle("Pt(#frac{GeV}{c})");
  h1009->GetYaxis()->SetTitle("counts");
  h1009->SetFillColor(hfillcolor);

  TH1F* h1010 = new TH1F("h1010","Jet1 Eta",100,-6.5,6.5);
  h1010->GetXaxis()->SetTitle("#eta");
  h1010->GetYaxis()->SetTitle("counts");
  h1010->SetFillColor(hfillcolor);
  h1010->SetMinimum(0);

  TH1F* h1011 = new TH1F("h1011","Jet1 Phi",80,-pi,pi);
  h1011->GetXaxis()->SetTitle("#phi");
  h1011->GetYaxis()->SetTitle("counts");
  h1011->SetFillColor(hfillcolor);
  h1011->SetMinimum(0);

  //2jet                                                                                                                                                    
  TH1F* h1012 = new TH1F("h1012","Jet2 E",100,0,1000.);                           /* RECO JETS */
  h1012->GetXaxis()->SetTitle("E(GeV)");
  h1012->GetYaxis()->SetTitle("counts");
  h1012->SetFillColor(hfillcolor);

  TH1F* h1013 = new TH1F("h1013","Jet2 Pt",100,0,200);
  h1013->GetXaxis()->SetTitle("Pt(#frac{GeV}{c})");
  h1013->GetYaxis()->SetTitle("counts");
  h1013->SetFillColor(hfillcolor);

  TH1F* h1014 = new TH1F("h1014","Jet2 Eta",100,-6.5,6.5);
  h1014->GetXaxis()->SetTitle("#eta");
  h1014->GetYaxis()->SetTitle("counts");
  h1014->SetFillColor(hfillcolor);
  h1014->SetMinimum(0);

  TH1F* h1015 = new TH1F("h1015","Jet2 Phi",80,-pi,pi);
  h1015->GetXaxis()->SetTitle("#phi");
  h1015->GetYaxis()->SetTitle("counts");
  h1015->SetFillColor(hfillcolor);
  h1015->SetMinimum(0);

  //deltaJets (in modulo)                                                                                                                                    
  TH1F* h1016 = new TH1F("h1016"," DeltaEta 1st-2nd Jets in Pt with cuts applied",100,0.,8.);
  h1016->GetXaxis()->SetTitle("#Delta#eta(Jet1-Jet2)");
  h1016->GetYaxis()->SetTitle("counts");
  h1016->SetFillColor(hfillcolor);
  h1016->SetMinimum(0);

  //eta(jet1)*eta(jet2)<0                                                                                                                                   
  TH1F* h1017 = new TH1F("h1017"," eta(jet1)*eta(jet2)<0",100,-1.5,1.5);
  h1017->GetXaxis()->SetTitle("#Delta#eta(Jet1-Jet2)");
  h1017->GetYaxis()->SetTitle("counts");
  h1017->SetFillColor(hfillcolor);
  h1017->SetMinimum(0);

  //Phi(g1+g2)-Phi(j1+j2)                                                                                                                                    
  TH1F* h1018 = new TH1F("h1018"," RECO |deltaPhi((Jet1+Jet2)-(Gamma1+Gamma2))|  ", 75, 0.,2*pi);                                                           
  h1018->GetXaxis()->SetTitle("#delta#phi((Jet1+Jet2)-(Gamma1+Gamma2))");
  h1018->GetYaxis()->SetTitle("counts");
  h1018->SetFillColor(hfillcolor);
  h1018->SetMinimum(0);

  //Zeppenfeld(g1+g2)                                                                                                                                        
  TH1F* h1019 = new TH1F("h1019","Zeppenfeld reco ga1+ga2",70,-10.,10.);
  h1019->GetXaxis()->SetTitle("Zeppenfeld ga1+ga2");
  h1019->GetYaxis()->SetTitle("counts");
  h1019->SetFillColor(hfillcolor);
  
  TH1F* h1020 = new TH1F("h1020","reco: Zeppenfeld(ga1+ga2)/|#eta(jet1)-#eta(jet2)|",70,-5.,5.);
  h1020->GetXaxis()->SetTitle("Zeppenfeld(ga1+ga2)/|#eta(jet1)-#eta(jet2)|");
  h1020->GetYaxis()->SetTitle("counts");
  h1020->SetFillColor(hfillcolor);
  
  //Zeppenfeld(j1)                                                                                                                                           
  TH1F* h1021 = new TH1F("h1021","Zeppenfeld reco jet1",70,-6.,6.);
  h1021->GetXaxis()->SetTitle("Zeppenfeld jet1");
  h1021->GetYaxis()->SetTitle("counts");
  h1021->SetFillColor(hfillcolor);
  
  TH1F* h1022 = new TH1F("h1022","reco: Zeppenfeld(jet1)/|#eta(jet1)-#eta(jet2)|",70,-2.,2.);
  h1022->GetXaxis()->SetTitle("Zeppenfeld(jet1)/|#eta(jet1)-#eta(jet2)|");
  h1022->GetYaxis()->SetTitle("counts");
  h1022->SetFillColor(hfillcolor);
  
  //Zeppenfeld(j2)                                                                                                                                          
  TH1F* h1023 = new TH1F("h1023","Zeppenfeld reco jet2",70,-6.,6.);
  h1023->GetXaxis()->SetTitle("Zeppenfeld jet2");
  h1023->GetYaxis()->SetTitle("counts");
  h1023->SetFillColor(hfillcolor);
  
  TH1F* h1024 = new TH1F("h1024","reco: Zeppenfeld(jet2)/|#eta(jet1)-#eta(jet2)|",70,-2.,2.);
  h1024->GetXaxis()->SetTitle("Zeppenfeld(jet2)/|#eta(jet1)-#eta(jet2)|");
  h1024->GetYaxis()->SetTitle("counts");
  h1024->SetFillColor(hfillcolor);
  
  TH1F* h1025 = new TH1F ("h1025","RECO Massa Inv 2 fotoni",100,100,150);
  h1025->GetXaxis()->SetTitleOffset(1.2);
  h1025->GetXaxis()->SetTitle("M(#frac{GeV}{c^{2}})");
  h1025->GetYaxis()->SetTitle("counts");
  h1025->SetFillColor(hfillcolor);  

  TH1F* h1026 = new TH1F ("h1026","RECO Massa Inv jets",100,0,1500);
  h1026->GetXaxis()->SetTitleOffset(1.2);
  h1026->GetXaxis()->SetTitle("M(j1 j2)(#frac{GeV}{c^{2}})");
  h1026->GetYaxis()->SetTitle("counts");
  h1026->SetFillColor(hfillcolor);

/************************************************/
  if (fChain == 0) return;

   Long64_t nentries = fChain->GetEntriesFast();
   Long64_t nbytes = 0, nb = 0;

/************************************************/

   float count_ga1_match=0,count_ga2_match=0,count_tot=0, count_ga1andga2_match=0;
   float count_ga1_match40=0,count_ga2_match20=0; // Pt reco gamma1 >40 GeV

   int recocount_evtga1andga2=0;  //numero di eventi con due fotoni reco energetici
   int count_and_match=0,count_xor_match=0, count_nor_match=0;
   int howmany_subcount_ga1=0,howmany_subcount_ga2=0;

   // dati in input

   float Pt_Cut = 0.;
   float checkcount_jetsforward=0; 
   float count_jet1jet2_Cut2=0;    // eventi con gen jet1 e jet2 che soddisfano al CUT2 (G1)
   float recojcount_jet1jet2_Cut2=0;    // eventi con recon jet1 e jet2 che soddisfano al CUT2 (R1) 
   float count_jet1jet2_Cut2_NMgengammas=0; // eventi con gen jet1 e jet2 che soddisfano al CUT2 e che non matchano gamma gen (G2)
   float recojcount_jet1jet2_Cut2_NMgengammas=0; // eventi con reco jet1 e jet2 che soddisfano al CUT2 e che non matchano gamma gen (R2)
   float recojcount_jet1jet2_Cut2_NMrecogammas=0; // eventi con reco jet1 e jet2 che soddisfano al CUT2 e che non matchano reco gen (R2)

   float count_match_q1Jet1 = 0;
   float count_match_q2Jet2 = 0;
   float count_match_q1Jet1q2Jet2 = 0;
   float recojcount_match_q1Jet1 = 0;
   float recojcount_match_q2Jet2 = 0;
   float recojcount_match_q1Jet1q2Jet2 = 0;
   int checkclusters=0, checkgammas=0;

   float RECOphotmismatch=0;
   int num_evts=0;
   int num_evts_beforePHID=0; //aggiunto 17 nov
   int num_evts_after_C1=0, num_evts_after_C2=0, num_evts_after_C3=0, num_evts_after_C4=0, num_evts_after_C5=0, num_evts_after_C6=0  ; 
   //# evts dopo i tagli
   int num_evts_2recophotons_ptmorethan30=0; //num di evts con reco photons pt>30 GeV prima associazione
   int morethanoneHiggs=0;
   int totnumbofHiggs=0;
   int ntotHgg=0;
   int n1totHgg_cuteta=0;

   int spike_gamma1=0,spike_gamma2=0, spike_jet1=0,spike_jet2=0;

   struct momento{
     float x;
     float y;
     float z;
     float pt;
     float p;
     float E;
     float eta;
     float phi;
     float dR, dEta, dPhi;
     float match_g2, match_clu2;
   } ;

   momento v2_GA, v2_R;  
   // sono i due momenti che calcolo nel caso 'teorico', con grandezze 
   //generate associate e nel caso 'reale', con grandezze reco 
   //[vai dopo RECO Jets]
   // v2_GA.match_g2=0; 
   v2_R.match_clu2=0;    // contatori che incrementano nell'evento  

   float PtRecoPho1Pho2Cut=0,PtEtaRecoPho1Pho2Cut=0;  

   /**************************************************************/
   int nev(0);     // equivale a nev=0 inizializzazione  
   // temp varables to ckeep track of the file being processed
   TString foldname("");
   TString currfilename("");
   int ifile(0);
   int nfiles = ((TChain*)fChain)->GetListOfFiles()->GetEntries();


   // loop over events
   for (Long64_t jentry=0; jentry<nentries;jentry++) { /* 2 */
     Long64_t ientry = LoadTree(jentry);

     if (ientry < 0) break;
     nb = fChain->GetEntry(jentry);   nbytes += nb;
     // if (Cut(ientry) < 0) continue;
     /**************************************************************/

     int count_photonsMC=0;
     num_evts++;

     if(int(num_evts)%1000 == 0 ) cout << "Evento: " << num_evts << endl;
      // print name of crrent file
      currfilename = TString(fChain->GetCurrentFile()->GetName());
      if(currfilename != foldname) {
         ifile++;
         cout << "Opening file " << ifile << " of "  << nfiles << "\n"
              << currfilename  << "\n"
              << "------------------------------"
              << endl;
         foldname = currfilename;
      }


     //GENERATOR PART
     int subcount_ga1_match=0,subcount_ga2_match=0,subcount_ga1andga2_match=0;
     int numbofHiggs=0; //per evento
     //   int subcount_and_match=0; // count_xor_match=0, count_nor_match=0;
     int N1_geng1g2=0;  //segnala attivazione cut1 per gen gamma Higgs (fabs(eta)|1,2 < 2.5  
     /***********************************************************************/
     //            ANALISI PARTICELLE MC: GEN PHOTONS, GEN PARTONS
     /***********************************************************************/

     int ig1=-11, ig2=-11, iH=-11; /* indici selezione fotoni figli e Higgs */  
     float gamma[2][4]={ {-100, -100, -100, -100 } , {-100, -100, -100, -100} };
     //float BBgamma[2][4]={ {-3, -3,-3, -3 } , {-3, -3,-3, -3 }}; /* fotoni Born Box */  
     int  BBgammaidx1=-9, BBgammaidx2=-9;
     int BBphotons=0;
     float temp_energy1=0,temp_energy2=0;
     /* array 2 fotoni figli Higgs g1 g2 per ogni evento */      
     //float reco[2][4]={ {-3, -3,-3, -3 } , {-3, -3,-3, -3 }};   
     /* array 2 fotoni reco associati a Higgs g1 g2 per ogni evento */    
     float delta_eta1;
     float delta_eta2;
     float delta_phi1;
     float delta_phi2;
     float delta_R1=100;
     float delta_R2=100;  
     float delta_R3=100;
     float delta_R4=100; 
     int quark=0;

     float tempparton[20][4]; 
     // di questi vorro' solo due righe (i due quark adronizzanti in jet), in status 3, a Pt maggiore
     //4: E, Pt, Eta, Phi partoni
     for(int i=0;i<20;i++){      
       for(int j=0;j<4;j++){ tempparton[i][j]=-1000;      
       }
     }

     for(int i=0; i<nMC;++i) {     /* 3MC */


       if(pdgIdMC[i]==22 && statusMC[i]==3) count_photonsMC++; //activate for g+jet only!!
       /***********variables from GammaJetAnalyzer.cc***********/
       //   Int_t nMC
       //   Int_t pdgIdMC 
       //   Int_t statusMC
       //   Int_t motherIDMC
       // Most MC particles have mass, but why do photons (status=1 and 3) have mass?
       //m_tree->Branch("massMC ",&massMC ,"massMC[nMC]/F");
       //to a good approximation, m = sqrt(e^2 - (pt*cosh(eta))^2), when m>1e-6 GeV
       //   Float_t ptMC
       //   Float_t eMC
       //   Float_t etaMC
       //   Float_t phiMC

       if(pdgIdMC[i]==25) {    /* E Pt eta phi Higgs*/
  
	 numbofHiggs++;
	 totnumbofHiggs++;
	 iH= i;

	 h1->Fill(eMC[iH]);
	 h2->Fill(ptMC[iH]);
	 h3->Fill(etaMC[iH]);
	 h4->Fill(phiMC[iH]);

	 //          std::cout<< "Daughters indices:"<< Dau1[0]<<" "<<Dau2[0]<<std::endl;
	 //          std::cout<< "Higgs array index:"<< i <<std::endl;
	 
	 int Dau[1][2];
	 
	 Dau[0][0]=-1;
	 Dau[0][1]=-1;
	 
	 int k=0;
	 
	 for(int j=0;j<nMC;j++){

           if( pdgIdMC[j] == 2212 ) continue;
           if( motherIDMC[j]>149) continue;
	   if(pdgIdMC[motherIDMC[j]]==25 && statusMC[j]==3 && pdgIdMC[j]==22){
	     
	     /** motherIDMC è l'indice nell'array MC della madre della particella **/
        
	     Dau[0][k]=j;
	     k++;
	     //cout<<"Dau[0]["<<k<<"]="<<j<<endl;

	   }
	 }
               
	 if(k>2){cout<<"more than 2 photons in Status 3 for this sample. "<<endl;}
               
	 ig1=Dau[0][0];                   
	 ig2=Dau[0][1];
	 
	 if(ptMC[ig1]<ptMC[ig2]){   //correzione (26 ottobre): anche i due fotoni Higgs ordinati in Pt 

	   ig1=Dau[0][1];                   
	   ig2=Dau[0][0];
	 }

	 if(ig1!=-11 && ig2!=-11){

	   gamma[0][0]= eMC[ig1];
	   gamma[0][1]= ptMC[ig1];
	   gamma[0][2]= etaMC[ig1]; 
	   gamma[0][3]= phiMC[ig1];
	   h5->Fill(gamma[0][0]);   /* 1st highest Pt photon  */
	   h6->Fill(gamma[0][1]);
	   h7->Fill(gamma[0][2]);
	   h8->Fill(gamma[0][3]);      
	
	   gamma[1][0]= eMC[ig2];
	   gamma[1][1]= ptMC[ig2];
	   gamma[1][2]= etaMC[ig2]; 
	   gamma[1][3]= phiMC[ig2];
	   h9->Fill(gamma[1][0]);    /* 2nd highest Pt photon */
	   h10->Fill(gamma[1][1]);
	   h11->Fill(gamma[1][2]);
	   h12->Fill(gamma[1][3]); 
	 }          

	 ntotHgg++; 
  
	 /*if(fabs(gamma[0][2])<2.5 && fabs(gamma[1][2])<2.5){
	   n1totHgg_cuteta++;
	   N1_geng1g2=1; 
	   }*/ // N1 cut1 sui fotoni attivato

       } /* id= Higgs */

 
       if(pdgIdMC[i]==22 && pdgIdMC[motherIDMC[i]]!=25 && statusMC[i]==3){         /* BORN BOX */

         BBphotons++; //check 2 fotoni in status 3 Born Box
	 // cout<<" i:"<<i<<" Energy:"<<eMC[i]<<endl;
	 if(eMC[i]>temp_energy1)     {
	   
	   BBgammaidx2=BBgammaidx1;
	   BBgammaidx1=i;
	   temp_energy2=temp_energy1;
	   temp_energy1=eMC[i];
	 }
	 
	 else if(eMC[i]>temp_energy2 && eMC[i]<temp_energy1){
	   
	   BBgammaidx2=i;
	   temp_energy2=eMC[i]; 
	 }                     
 
         // cout<<"BBgammaidx1: "<<BBgammaidx1<<", BBgammaidx2: "<<BBgammaidx2<<endl;
	 
	 /* additional photons in Higgs events OR background events born,box,gamma+jet */
	 h13->Fill(eMC[i]);
	 h14->Fill(ptMC[i]);
	 h15->Fill(etaMC[i]);
	 h16->Fill(phiMC[i]);
       }                              /* BORN BOX */
                   
       for(int j=i+1; j<nMC;++j) {/* 4MC */   //MASSA INVARIANTE FOTONI GEN

	 if(pdgIdMC[i]==22 && pdgIdMC[j]==22 && pdgIdMC[motherIDMC[i]]==25){ /* 5MC */ //tutti i fotoni
	   
	   SumP[1]= ptMC[i]*cos(phiMC[i])+ ptMC[j]*cos(phiMC[j]);
	   SumP[2]= ptMC[i]*sin(phiMC[i])+ ptMC[j]*sin(phiMC[j]);
	   SumP[3]= ptMC[i]*sinh(etaMC[i])+ ptMC[j]*sinh(etaMC[j]);
	   SumP[0]= eMC[i]+ eMC[j];

	   Gen_inv_mass= (SumP[0]*SumP[0]-(SumP[1]*SumP[1]+SumP[2]*SumP[2]+SumP[3]*SumP[3]));

	   if(Gen_inv_mass>0){
	     Gen_inv_mass=sqrt(Gen_inv_mass);
	     h17->Fill(Gen_inv_mass);     
	   } 

	 } /* 5MC */ //tutti i fotoni

	 if(pdgIdMC[i]==22 && statusMC[i]==3 && pdgIdMC[j]==22 && statusMC[j]==3 && pdgIdMC[motherIDMC[i]]!=25){ /* 5MC */ 
	   //tutti fotoni Born e Box

	   SumP[1]= ptMC[i]*cos(phiMC[i])+ ptMC[j]*cos(phiMC[j]);
	   SumP[2]= ptMC[i]*sin(phiMC[i])+ ptMC[j]*sin(phiMC[j]);
	   SumP[3]= ptMC[i]*sinh(etaMC[i])+ ptMC[j]*sinh(etaMC[j]);
	   SumP[0]= eMC[i]+ eMC[j];

	   Gen_inv_mass= (SumP[0]*SumP[0]-(SumP[1]*SumP[1]+SumP[2]*SumP[2]+SumP[3]*SumP[3]));
	   
	   if(Gen_inv_mass>0){
	     Gen_inv_mass=sqrt(Gen_inv_mass);
	     h17->Fill(Gen_inv_mass);     
	   } 

	 } /* 5MC */ //tutti i fotoni
	 
       }   /* 4MC */ 

       //cout << "i : " << i << " pdgIdMC[i]: " << pdgIdMC[i] << " statusMC[i]: " << statusMC[i] << " mother: " << motherIDMC[i] << endl;
       if( motherIDMC[i] > 149 ) continue;
       if(pdgIdMC[i]!=25 && pdgIdMC[i]!=2212 && pdgIdMC[i]!=22 && statusMC[i]==3 
          && pdgIdMC[motherIDMC[i]]!=25 && pdgIdMC[motherIDMC[i]]!=22 && statusMC[motherIDMC[i]]==3){ /* partoni generati, check listdraw.txt */
	 
	 tempparton[quark][0]= eMC[i];
	 tempparton[quark][1]= ptMC[i];
	 tempparton[quark][2]= etaMC[i]; 
	 tempparton[quark][3]= phiMC[i];
        
	 quark++;
       }
                                                                                                                 
     }     /* 3MC */  // for (array MC dell'ntupla di un "evento")   

     if(BBgammaidx1!=-9 && pdgIdMC[motherIDMC[BBgammaidx1]]!=25) {   /* BORN || BOX || gamma+jet */
       
       gamma[0][0]= eMC[BBgammaidx1];   /* fotone piu' energetico */
       gamma[0][1]= ptMC[BBgammaidx1];
       gamma[0][2]= etaMC[BBgammaidx1];
       gamma[0][3]= phiMC[BBgammaidx1];
       
       h5->Fill(gamma[0][0]);    
       h6->Fill(gamma[0][1]);
       h7->Fill(gamma[0][2]);
       h8->Fill(gamma[0][3]);      
     }
     
     if(BBgammaidx2!=-9 && pdgIdMC[motherIDMC[BBgammaidx2]]!=25) {  /* BORN || BOX  */
       
       gamma[1][0]= eMC[BBgammaidx2];   /* fotone meno energetico */
       gamma[1][1]= ptMC[BBgammaidx2];
       gamma[1][2]= etaMC[BBgammaidx2];
       gamma[1][3]= phiMC[BBgammaidx2];  
       
       h9->Fill(gamma[1][0]);    
       h10->Fill(gamma[1][1]);
       h11->Fill(gamma[1][2]);
       h12->Fill(gamma[1][3]);     
     }

     /*******************************************/

     const int nMCMax = 150;
     float vect[nMCMax];
     int indices_partons[nMCMax]; //4NOV
     
     for(int i=0;i<nMCMax;i++){ vect[i]=0;}

     int tot_quark=0;
     tot_quark=quark;
     
     float parton[2][4];  //variabili in cui salvero' E Pt Eta Phi dei due partoni a Pt maggiore
     
     for(int i=0;i<2;i++){//inizializzione
       
       for(int j=0;j<4;j++){parton[i][j]=-1000;}      
     }    
     
     for(int i=0; i<tot_quark;i++) { vect[i]=tempparton[i][1]; }

     int index_parton1=-9, index_parton2=-9;

     if(tot_quark>0){
       ordinamento(vect,indices_partons,nMCMax);  //Pt quark ordinati in ordine decrescente
                 
       if(vect[0]>0){index_parton1=indices_partons[0];}
       if(vect[1]>0){index_parton2=indices_partons[1];}
     }
     
     //if(index_parton1==-9 || index_parton2==-9){cout<<"meno di 2 partoni gen adronizzanti nell'evento"<<endl;}
     
     if(index_parton1!=-9){

       parton[0][0]=tempparton[index_parton1][0];
       parton[0][1]=tempparton[index_parton1][1];
       parton[0][2]=tempparton[index_parton1][2];
       parton[0][3]=tempparton[index_parton1][3];
     }
     
     if(index_parton2!=-9){
       
       parton[1][0]=tempparton[index_parton2][0];
       parton[1][1]=tempparton[index_parton2][1];
       parton[1][2]=tempparton[index_parton2][2];
       parton[1][3]=tempparton[index_parton2][3];
     }            

     if(index_parton1!=-9){                           //Pt Eta due partoni a Pt maggiore

       h48->Fill(parton[0][1]);
       h46->Fill(parton[0][2]);
     }
     
     if(index_parton2!=-9){

       h49->Fill(parton[1][1]);  
       h47->Fill(parton[1][2]); 
     }

     if(filter2GammaJet && count_photonsMC>1){continue;} //active only for g+jet sample

     /***********************************************************************/
     //                           RECO PHOTONS AND JETS
     /***********************************************************************/            //from here: 21 NOVEMBER

     float threshold_deltaR=0.05;

      vector<bool> photassocMC, photIdentif, jetassocphot;

      for(int i=0; i<nPhot; i++){      

          //reco photons geometrical associated to MC photons
  
          delta_eta1= etaPhot[i]-gamma[0][2];            
          delta_eta2= etaPhot[i]-gamma[1][2];           
          delta_phi1= fabs(phiPhot[i]-gamma[0][3]);     
          delta_phi2= fabs(phiPhot[i]-gamma[1][3]);      
     
          if( (2*pi-delta_phi1)< delta_phi1 ){delta_phi1=2*pi-delta_phi1; }
          if( (2*pi-delta_phi2)< delta_phi2 ){delta_phi2=2*pi-delta_phi2; }
     
          delta_R1= sqrt( delta_eta1*delta_eta1+delta_phi1*delta_phi1 );      /* + energetico*/
          delta_R2= sqrt( delta_eta2*delta_eta2+delta_phi2*delta_phi2 );      /* - energetico*/
     
          h19->Fill(delta_R1);
          h20->Fill(delta_R2);

         if(delta_R1<threshold_deltaR) count_ga1_match++;
         if(delta_R2<threshold_deltaR) count_ga2_match++;
     
          if(delta_R1<threshold_deltaR || delta_R2<threshold_deltaR) photassocMC.push_back(1);
          else photassocMC.push_back(0);
	  
          //photonID [|eta|<2.5 included]

      //RECO PART
        vector<bool> idpass(7);
        //previous PHID
	//if(cutID(i, mediumid, &idpass)==1) photIdentif.push_back(1);  //cutID select gammas with |eta|<2.5, inside detector
	if(cutID(i, looseid, &idpass)==1) photIdentif.push_back(1);  //cutID select gammas with |eta|<2.5, inside detector
	else photIdentif.push_back(0);
             
	//if(pid_isTight[i]) photIdentif.push_back(1);
	//else photIdentif.push_back(0);
      
      }
 
      vector<int> firsttwophotassocMC = firsttwo(ptPhot,&photassocMC); //do I need it???
      vector<int> firsttwophotIdentif = firsttwo(ptPhot,&photIdentif); //find reco_gamma1, reco_gamma2

      int reco_gamma1=-9,reco_gamma2=-9;  

      if( firsttwophotIdentif.at(0)>-1 && firsttwophotIdentif.at(1)>-1 ){
           reco_gamma1= firsttwophotIdentif.at(0);
           reco_gamma2= firsttwophotIdentif.at(1);
      }

     // associazione reco jets reco_gamma1 reco_gamma2

     float recoga1Jet_deltaR,recoga2Jet_deltaR;
       
     for(int i=0;i<nJet_pfakt5;i++){/* 3 */  
	 
     recoga1Jet_deltaR=500;
     recoga2Jet_deltaR=500; //each time by default

       if(reco_gamma1!=-9){
	 
	  delta_eta1=etaJet_pfakt5[i]-etaPhot[reco_gamma1];
	  delta_phi1=fabs(phiJet_pfakt5[i]-phiPhot[reco_gamma1]);              
	 
	  if( (2*pi-delta_phi1)< delta_phi1 )delta_phi1=2*pi-delta_phi1;               
	  recoga1Jet_deltaR= sqrt(delta_eta1*delta_eta1 + delta_phi1*delta_phi1);
       }

       if(reco_gamma2!=-9){	 
	  delta_eta2=etaJet_pfakt5[i]-etaPhot[reco_gamma2];
	  delta_phi2=fabs(phiJet_pfakt5[i]-phiPhot[reco_gamma2]);
	 
	  if( (2*pi-delta_phi2)< delta_phi2 )delta_phi2=2*pi-delta_phi2; 	           
	  recoga2Jet_deltaR= sqrt(delta_eta2*delta_eta2 + delta_phi2*delta_phi2);  
       }
       
       if(recoga1Jet_deltaR<0.25 || recoga2Jet_deltaR<0.25 || fabs(etaJet_pfakt5[i])>5.) jetassocphot.push_back(0);  //rejected cases     
       else jetassocphot.push_back(1); //accepted cases: jets NOT associated either to reco_gamma1 or to reco_gamma2 and with |eta|<5. 
       
     }/* 3 */
 
      vector<int> firsttwojetassocphot = firsttwo(ptJet_pfakt5,&jetassocphot);    //trovo idxjet1, idxjet2

     int idxjet1=-9,idxjet2=-9;  

      if( firsttwojetassocphot.at(0)>-1 && firsttwojetassocphot.at(1)>-1 ){  
      idxjet1= firsttwojetassocphot.at(0);
      idxjet2= firsttwojetassocphot.at(1);      
      }
  
     /************************************************ observables to plot *************************************************/

     ///////////////////////// plot eta1*eta2<0
     float cond_a=5;
     if(idxjet1!=-9 && idxjet2!=-9) cond_a= etaJet_pfakt5[idxjet1]*etaJet_pfakt5[idxjet2];    //non viene sfruttata?

     ///////////////////////// plot deltaPhi
     //jets
     double Jets_inv_mass=-100;
     float PhiTotJets=1000, EtaTotJets=1000;

     TLorentzVector jet1, jet2;        
     jet1.SetPtEtaPhiE(ptJet_pfakt5[firsttwojetassocphot.at(0)],etaJet_pfakt5[firsttwojetassocphot.at(0)],
                       phiJet_pfakt5[firsttwojetassocphot.at(0)], eJet_pfakt5[firsttwojetassocphot.at(0)]);
     jet2.SetPtEtaPhiE(ptJet_pfakt5[firsttwojetassocphot.at(1)],etaJet_pfakt5[firsttwojetassocphot.at(1)],
                       phiJet_pfakt5[firsttwojetassocphot.at(1)], eJet_pfakt5[firsttwojetassocphot.at(1)]);

     TLorentzVector jetsum = jet1 + jet2;
 
     if(idxjet1!=-9 && idxjet2!=-9){
     PhiTotJets=jetsum.Phi();
         if(PhiTotJets<0){PhiTotJets= pi*2.+ PhiTotJets;}
     EtaTotJets=jetsum.Eta();
     Jets_inv_mass=jetsum.M();
     }

     //photons
     double Reco_inv_mass=-100; 
     float PhiTotgammas=1000, EtaTotgammas=1000;

     TLorentzVector phot1, phot2;        
     phot1.SetPtEtaPhiE(ptPhot[firsttwophotIdentif.at(0)],etaPhot[firsttwophotIdentif.at(0)],phiPhot[firsttwophotIdentif.at(0)],ePhot[firsttwophotIdentif.at(0)]);
     phot2.SetPtEtaPhiE(ptPhot[firsttwophotIdentif.at(1)],etaPhot[firsttwophotIdentif.at(1)],phiPhot[firsttwophotIdentif.at(1)],ePhot[firsttwophotIdentif.at(1)]);

     TLorentzVector higgs = phot1 + phot2;
 
     if(reco_gamma1!=-9 && reco_gamma2!=-9){
     PhiTotgammas=higgs.Phi();
             if(PhiTotgammas<0){PhiTotgammas= pi*2.+ PhiTotgammas;}   
     EtaTotgammas=higgs.Eta();
     Reco_inv_mass=higgs.M();     ////////////////////plot Reco_inv_mass
     }

     /////
     float deltaPhi=1000;
     if(PhiTotJets!=1000 && PhiTotgammas!=1000){ 
       deltaPhi=fabs(PhiTotJets-PhiTotgammas);
     }
     
     ////////////////////plot Zeppenfeld
     float Zepp_ga1ga2=1000, Zepp_jet1=1000, Zepp_jet2=1000; 
     if(idxjet1!=-9 && idxjet2!=-9 && reco_gamma1!=-9 && reco_gamma2!=-9){ 
     Zepp_ga1ga2= EtaTotgammas-(etaJet_pfakt5[idxjet1]+etaJet_pfakt5[idxjet2])/2.; }

     if(idxjet1!=-9 && idxjet2!=-9){
     Zepp_jet1= etaJet_pfakt5[idxjet1] -(etaJet_pfakt5[idxjet1]+etaJet_pfakt5[idxjet2])/2.;
     Zepp_jet2= etaJet_pfakt5[idxjet2] -(etaJet_pfakt5[idxjet1]+etaJet_pfakt5[idxjet2])/2.;
     }

     /*************************************************************** CUTS ********************************************************/

     /*******/
     // 1CUT     //INSERIRE QUI 1 CUT 
     /*******/ 
    
     if(reco_gamma1==-9 || reco_gamma2==-9 || ptPhot[reco_gamma1]<50.0 || ptPhot[reco_gamma2]<30.0 ){continue;}

     num_evts_after_C1++;
       
     //plots after 1CUT
     /*********************************************/
     if(reco_gamma1==-9){spike_gamma1++;}
     if(reco_gamma2==-9){spike_gamma2++;}
     if(idxjet1==-9){spike_jet1++;}
     if(idxjet2==-9){spike_jet2++;}

     //1 photon    
     if(reco_gamma1!=-9){
     h500->Fill(ePhot[reco_gamma1]); 
     h501->Fill(ptPhot[reco_gamma1]);
     h502->Fill(etaPhot[reco_gamma1]);
     h503->Fill(phiPhot[reco_gamma1]);
     }
     //2 photon
     if(reco_gamma2!=-9){
     h504->Fill(ePhot[reco_gamma2]);
     h505->Fill(ptPhot[reco_gamma2]);
     h506->Fill(etaPhot[reco_gamma2]);
     h507->Fill(phiPhot[reco_gamma2]);
     }
     //1 jet 
     if(idxjet1!=-9){ 
     h508->Fill(eJet_pfakt5[idxjet1]);
     h509->Fill(ptJet_pfakt5[idxjet1]);
     h510->Fill(etaJet_pfakt5[idxjet1]);
     h511->Fill(phiJet_pfakt5[idxjet1]);
     }
     //2 jet              
     if(idxjet2!=-9){
     h512->Fill(eJet_pfakt5[idxjet2]);
     h513->Fill(ptJet_pfakt5[idxjet2]);
     h514->Fill(etaJet_pfakt5[idxjet2]);
     h515->Fill(phiJet_pfakt5[idxjet2]);
     }
     //deltaEta jets1,2
     if(idxjet1!=-9 && idxjet2!=-9){
       h516->Fill(fabs(etaJet_pfakt5[idxjet1]-etaJet_pfakt5[idxjet2]));
       //eta_j1*eta_j2<0
       h517->Fill(etaJet_pfakt5[idxjet1]*etaJet_pfakt5[idxjet2] /fabs(etaJet_pfakt5[idxjet1]*etaJet_pfakt5[idxjet2]));
     }
     if(idxjet1!=-9 && idxjet2!=-9 && reco_gamma1!=-9 && reco_gamma2!=-9){
       //Phi(g1+g2)-Phi(j1+j2)
       if(deltaPhi!=1000) h518->Fill(deltaPhi); 

       //Zepp(g1+g2)
       h519->Fill(Zepp_ga1ga2);
       h520->Fill(Zepp_ga1ga2/fabs(etaJet_pfakt5[idxjet1]-etaJet_pfakt5[idxjet2]));
     }
     if(idxjet1!=-9 && idxjet2!=-9){
       //Zepp(jet1)
       h521->Fill(Zepp_jet1);
       h522->Fill(Zepp_jet1/fabs(etaJet_pfakt5[idxjet1]-etaJet_pfakt5[idxjet2]));
       //Zepp(jet2)
       h523->Fill(Zepp_jet2);
       h524->Fill(Zepp_jet2/fabs(etaJet_pfakt5[idxjet1]-etaJet_pfakt5[idxjet2]));
       //jets inv mass
       h526->Fill(Jets_inv_mass);
     }
     if( reco_gamma1!=-9 && reco_gamma2!=-9){
       //Reco_inv_mass
       if(Reco_inv_mass>0){ h525->Fill(Reco_inv_mass); }
       // cout<<"inv_mass: "<<Reco_inv_mass<<endl;
     }
     /**********************************************/
     /*******/
     //2CUT     //INSERIRE QUI 2 CUT 
     /*******/


     if(idxjet1==-9 || idxjet2==-9 || ptJet_pfakt5[idxjet1]<ptJet1_cut || ptJet_pfakt5[idxjet2]<ptJet2_cut){continue;}     
     
     num_evts_after_C2++;
     //inserire plots dopo 2CUT
     /*********************************************/
     //1 fotone
     h600->Fill(ePhot[reco_gamma1]); 
     h601->Fill(ptPhot[reco_gamma1]);
     h602->Fill(etaPhot[reco_gamma1]);
     h603->Fill(phiPhot[reco_gamma1]);
     //2 fotone
     h604->Fill(ePhot[reco_gamma2]);
     h605->Fill(ptPhot[reco_gamma2]);
     h606->Fill(etaPhot[reco_gamma2]);
     h607->Fill(phiPhot[reco_gamma2]);
     //1 jet              
     h608->Fill(eJet_pfakt5[idxjet1]);
     h609->Fill(ptJet_pfakt5[idxjet1]);
     h610->Fill(etaJet_pfakt5[idxjet1]);
     h611->Fill(phiJet_pfakt5[idxjet1]);
     //2 jet              
     h612->Fill(eJet_pfakt5[idxjet2]);
     h613->Fill(ptJet_pfakt5[idxjet2]);
     h614->Fill(etaJet_pfakt5[idxjet2]);
     h615->Fill(phiJet_pfakt5[idxjet2]);
     //deltaEta jets1,2
     h616->Fill(fabs(etaJet_pfakt5[idxjet1]-etaJet_pfakt5[idxjet2]));
     //eta_j1*eta_j2<0
     h617->Fill(etaJet_pfakt5[idxjet1]*etaJet_pfakt5[idxjet2] /fabs(etaJet_pfakt5[idxjet1]*etaJet_pfakt5[idxjet2]));
     //Phi(g1+g2)-Phi(j1+j2)
     if(deltaPhi!=1000) h618->Fill(deltaPhi);
     //Zepp(g1+g2)
     h619->Fill(Zepp_ga1ga2);
     h620->Fill(Zepp_ga1ga2/fabs(etaJet_pfakt5[idxjet1]-etaJet_pfakt5[idxjet2]));
     //Zepp(jet1)
     h621->Fill(Zepp_jet1);
     h622->Fill(Zepp_jet1/fabs(etaJet_pfakt5[idxjet1]-etaJet_pfakt5[idxjet2]));
     //Zepp(jet2)
     h623->Fill(Zepp_jet2);
     h624->Fill(Zepp_jet2/fabs(etaJet_pfakt5[idxjet1]-etaJet_pfakt5[idxjet2]));
     //Reco_inv_mass
     if(Reco_inv_mass>0){ h625->Fill(Reco_inv_mass); }
     //jets inv mass
     h626->Fill(Jets_inv_mass);
     /**********************************************/
     /*******/
     //3CUT     //INSERIRE QUI 3 CUT 
     /*******/     
     if(fabs(etaJet_pfakt5[idxjet1]-etaJet_pfakt5[idxjet2])<detaJets_cut){continue;} 

     num_evts_after_C3++;
     
     //inserire plots dopo 3CUT
     /*********************************************/
     //1 fotone
     h700->Fill(ePhot[reco_gamma1]); 
     h701->Fill(ptPhot[reco_gamma1]);
     h702->Fill(etaPhot[reco_gamma1]);
     h703->Fill(phiPhot[reco_gamma1]);
     //2 fotone
     h704->Fill(ePhot[reco_gamma2]);
     h705->Fill(ptPhot[reco_gamma2]);
     h706->Fill(etaPhot[reco_gamma2]);
     h707->Fill(phiPhot[reco_gamma2]);
     //1 jet              
     h708->Fill(eJet_pfakt5[idxjet1]);
     h709->Fill(ptJet_pfakt5[idxjet1]);
     h710->Fill(etaJet_pfakt5[idxjet1]);
     h711->Fill(phiJet_pfakt5[idxjet1]);
     //2 jet              
     h712->Fill(eJet_pfakt5[idxjet2]);
     h713->Fill(ptJet_pfakt5[idxjet2]);
     h714->Fill(etaJet_pfakt5[idxjet2]);
     h715->Fill(phiJet_pfakt5[idxjet2]);
     //deltaEta jets1,2
     h716->Fill(fabs(etaJet_pfakt5[idxjet1]-etaJet_pfakt5[idxjet2]));
     //eta_j1*eta_j2<0
     h717->Fill(etaJet_pfakt5[idxjet1]*etaJet_pfakt5[idxjet2] /fabs(etaJet_pfakt5[idxjet1]*etaJet_pfakt5[idxjet2]));
     //Phi(g1+g2)-Phi(j1+j2)
     if(deltaPhi!=1000) h718->Fill(deltaPhi);
     //Zepp(g1+g2)
     h719->Fill(Zepp_ga1ga2);
     h720->Fill(Zepp_ga1ga2/fabs(etaJet_pfakt5[idxjet1]-etaJet_pfakt5[idxjet2]));
     //Zepp(jet1)
     h721->Fill(Zepp_jet1);
     h722->Fill(Zepp_jet1/fabs(etaJet_pfakt5[idxjet1]-etaJet_pfakt5[idxjet2]));
     //Zepp(jet2)
     h723->Fill(Zepp_jet2);
     h724->Fill(Zepp_jet2/fabs(etaJet_pfakt5[idxjet1]-etaJet_pfakt5[idxjet2]));
     //Reco_inv_mass
     if(Reco_inv_mass>0){ h725->Fill(Reco_inv_mass); }
     /**********************************************/
     //jets inv mass
     h726->Fill(Jets_inv_mass);

/*******/
//4CUT     //INSERIRE QUI 4 CUT 
/*******/

    if(fabs(Zepp_ga1ga2)>zepp_cut){continue;}

    num_evts_after_C4++;

     //inserire plots dopo 4CUT
     /*********************************************/
     //1 fotone
     h800->Fill(ePhot[reco_gamma1]); 
     h801->Fill(ptPhot[reco_gamma1]);
     h802->Fill(etaPhot[reco_gamma1]);
     h803->Fill(phiPhot[reco_gamma1]);
     //2 fotone
     h804->Fill(ePhot[reco_gamma2]);
     h805->Fill(ptPhot[reco_gamma2]);
     h806->Fill(etaPhot[reco_gamma2]);
     h807->Fill(phiPhot[reco_gamma2]);
     //1 jet              
     h808->Fill(eJet_pfakt5[idxjet1]);
     h809->Fill(ptJet_pfakt5[idxjet1]);
     h810->Fill(etaJet_pfakt5[idxjet1]);
     h811->Fill(phiJet_pfakt5[idxjet1]);
     //2 jet              
     h812->Fill(eJet_pfakt5[idxjet2]);
     h813->Fill(ptJet_pfakt5[idxjet2]);
     h814->Fill(etaJet_pfakt5[idxjet2]);
     h815->Fill(phiJet_pfakt5[idxjet2]);
     //deltaEta jets1,2
     h816->Fill(fabs(etaJet_pfakt5[idxjet1]-etaJet_pfakt5[idxjet2]));
     //eta_j1*eta_j2<0
     h817->Fill(etaJet_pfakt5[idxjet1]*etaJet_pfakt5[idxjet2] /fabs(etaJet_pfakt5[idxjet1]*etaJet_pfakt5[idxjet2]));
     //Phi(g1+g2)-Phi(j1+j2)
     if(deltaPhi!=1000) h818->Fill(deltaPhi);
     //Zepp(g1+g2)
     h819->Fill(Zepp_ga1ga2);
     h820->Fill(Zepp_ga1ga2/fabs(etaJet_pfakt5[idxjet1]-etaJet_pfakt5[idxjet2]));
     //Zepp(jet1)
     h821->Fill(Zepp_jet1);
     h822->Fill(Zepp_jet1/fabs(etaJet_pfakt5[idxjet1]-etaJet_pfakt5[idxjet2]));
     //Zepp(jet2)
     h823->Fill(Zepp_jet2);
     h824->Fill(Zepp_jet2/fabs(etaJet_pfakt5[idxjet1]-etaJet_pfakt5[idxjet2]));
     //Reco_inv_mass
     if(Reco_inv_mass>0){ h825->Fill(Reco_inv_mass); }
     //jets inv mass
     h826->Fill(Jets_inv_mass);
     /**********************************************/
     if(deltaPhi!=1000 && Reco_inv_mass>0) h827->Fill(deltaPhi,Reco_inv_mass);

/*******/
//5CUT     //INSERIRE QUI 5 CUT 
/*******/

     if(fabs(Jets_inv_mass)<mjj_cut){continue;} //right cut??

     num_evts_after_C5++;

     //inserire plots dopo 5CUT
     /*********************************************/
     //1 fotone
     h900->Fill(ePhot[reco_gamma1]); 
     h901->Fill(ptPhot[reco_gamma1]);
     h902->Fill(etaPhot[reco_gamma1]);
     h903->Fill(phiPhot[reco_gamma1]);
     //2 fotone
     h904->Fill(ePhot[reco_gamma2]);
     h905->Fill(ptPhot[reco_gamma2]);
     h906->Fill(etaPhot[reco_gamma2]);
     h907->Fill(phiPhot[reco_gamma2]);
     //1 jet              
     h908->Fill(eJet_pfakt5[idxjet1]);
     h909->Fill(ptJet_pfakt5[idxjet1]);
     h910->Fill(etaJet_pfakt5[idxjet1]);
     h911->Fill(phiJet_pfakt5[idxjet1]);
     //2 jet              
     h912->Fill(eJet_pfakt5[idxjet2]);
     h913->Fill(ptJet_pfakt5[idxjet2]);
     h914->Fill(etaJet_pfakt5[idxjet2]);
     h915->Fill(phiJet_pfakt5[idxjet2]);
     //deltaEta jets1,2
     h916->Fill(fabs(etaJet_pfakt5[idxjet1]-etaJet_pfakt5[idxjet2]));
     //eta_j1*eta_j2<0
     h917->Fill(etaJet_pfakt5[idxjet1]*etaJet_pfakt5[idxjet2] /fabs(etaJet_pfakt5[idxjet1]*etaJet_pfakt5[idxjet2]));
     //Phi(g1+g2)-Phi(j1+j2)
     if(deltaPhi!=1000) h918->Fill(deltaPhi);
     //Zepp(g1+g2)
     h919->Fill(Zepp_ga1ga2);
     h920->Fill(Zepp_ga1ga2/fabs(etaJet_pfakt5[idxjet1]-etaJet_pfakt5[idxjet2]));
     //Zepp(jet1)
     h921->Fill(Zepp_jet1);
     h922->Fill(Zepp_jet1/fabs(etaJet_pfakt5[idxjet1]-etaJet_pfakt5[idxjet2]));
     //Zepp(jet2)
     h923->Fill(Zepp_jet2);
     h924->Fill(Zepp_jet2/fabs(etaJet_pfakt5[idxjet1]-etaJet_pfakt5[idxjet2]));
     //Reco_inv_mass
     if(Reco_inv_mass>0){ 
         h925->Fill(Reco_inv_mass); 
         h927->Fill(Reco_inv_mass); 
     }

     //jets inv mass
     h926->Fill(Jets_inv_mass);

     /**********************************************/

     /*******/
     //6CUT     //INSERIRE QUI 6 CUT
     /*******/

    if(fabs(Reco_inv_mass-120)>5){continue;}

     num_evts_after_C6++;

     //inserire plots dopo 5CUT
     /*********************************************/
     //1 fotone
     h1000->Fill(ePhot[reco_gamma1]); 
     h1001->Fill(ptPhot[reco_gamma1]);
     h1002->Fill(etaPhot[reco_gamma1]);
     h1003->Fill(phiPhot[reco_gamma1]);
     //2 fotone
     h1004->Fill(ePhot[reco_gamma2]);
     h1005->Fill(ptPhot[reco_gamma2]);
     h1006->Fill(etaPhot[reco_gamma2]);
     h1007->Fill(phiPhot[reco_gamma2]);
     //1 jet              
     h1008->Fill(eJet_pfakt5[idxjet1]);
     h1009->Fill(ptJet_pfakt5[idxjet1]);
     h1010->Fill(etaJet_pfakt5[idxjet1]);
     h1011->Fill(phiJet_pfakt5[idxjet1]);
     //2 jet              
     h1012->Fill(eJet_pfakt5[idxjet2]);
     h1013->Fill(ptJet_pfakt5[idxjet2]);
     h1014->Fill(etaJet_pfakt5[idxjet2]);
     h1015->Fill(phiJet_pfakt5[idxjet2]);
     //deltaEta jets1,2
     h1016->Fill(fabs(etaJet_pfakt5[idxjet1]-etaJet_pfakt5[idxjet2]));
     //eta_j1*eta_j2<0
     h1017->Fill(etaJet_pfakt5[idxjet1]*etaJet_pfakt5[idxjet2] /fabs(etaJet_pfakt5[idxjet1]*etaJet_pfakt5[idxjet2]));
     //Phi(g1+g2)-Phi(j1+j2)
     if(deltaPhi!=1000) h1018->Fill(deltaPhi);
     //Zepp(g1+g2)
     h1019->Fill(Zepp_ga1ga2);
     h1020->Fill(Zepp_ga1ga2/fabs(etaJet_pfakt5[idxjet1]-etaJet_pfakt5[idxjet2]));
     //Zepp(jet1)
     h1021->Fill(Zepp_jet1);
     h1022->Fill(Zepp_jet1/fabs(etaJet_pfakt5[idxjet1]-etaJet_pfakt5[idxjet2]));
     //Zepp(jet2)
     h1023->Fill(Zepp_jet2);
     h1024->Fill(Zepp_jet2/fabs(etaJet_pfakt5[idxjet1]-etaJet_pfakt5[idxjet2]));
     //Reco_inv_mass
     if(Reco_inv_mass>0){ h1025->Fill(Reco_inv_mass); }
     //jets inv mass
     h1026->Fill(Jets_inv_mass);

     /**********************************************/


     if(numbofHiggs>1){
       cout<<"num. of Higgs nell'evento: " <<numbofHiggs<<endl;
       morethanoneHiggs++;
     }
     
     /*     std::cout<< " Fine evento. " << std::endl;
	    std::cout<< " " << std::endl;*/
     
     nev++;
     using namespace std;
   }  /* 2 */ // fine for dimensione albero ("eventi") 
   
   // float fraction;
   
   cout<<endl;
   cout<<"N0 (n.tot eventi): "<<num_evts<<endl;
   cout<<"Num.evts con 2 fotoni reco pt>30 e eta<2.5 (CUT1): "<<num_evts_beforePHID<<endl;
   cout<<"N1 (n.eventi dopo CUT1 && PHID): "<<num_evts_after_C1<<endl;
   cout<<"N2 (n.eventi dopo CUT2): "<<num_evts_after_C2<<endl;
   cout<<"N3 (n.eventi dopo CUT3): "<<num_evts_after_C3<<endl;
   cout<<"N4 (n.eventi dopo CUT4): "<<num_evts_after_C4<<endl;
   cout<<"N5 (n.eventi dopo CUT5): "<<num_evts_after_C5<<endl;
   cout<<"N6 (n.eventi dopo CUT6): "<<num_evts_after_C6<<endl;
   cout<<"num. reco gamma1 matchati "<<count_ga1_match<<endl;
   cout<<"num. reco gamma2 matchati "<<count_ga2_match<<endl;
   cout<<"num. di reco_gamma1=-9  "<<spike_gamma1<<endl;
   cout<<"num. di reco_gamma2=-9  "<<spike_gamma2<<endl;
   cout<<"num. di reco_jet1=-9  "<<spike_jet1<<endl;
   cout<<"num. di reco_jet2=-9  "<<spike_jet2<<endl;
   cout<<"num. eventi con 2 reco photons pt>30GeV: "<<num_evts_2recophotons_ptmorethan30<<endl;

   output_root->Write();
   output_root->Close();
   
}
   
   
/*****************************************************************************/
   
using namespace std;

void AnalysisTool::ordinamento(float *n, int*array, int dim) { 
  
  float W,tmp;
  
  for(int i=0;i!=dim;i++) {            
    array[i]=i;
  } 
  
  for (int j=0; j!=dim; j++)  
    {
      for (int k=0; k!=dim; k++)
	{
	  if(n[j]>n[k])
	    {
	      W=n[j];     
	      tmp=array[j];       
	      n[j]=n[k];
	      array[j]=array[k];
	      n[k]=W;
	      array[k]=tmp;
	    }
	}
    }
}                 //implementa funzione bubblesort da pc Desk/codes


//


//first tow indices Del Re  //implementare nel.h così come per PhotonID

vector<int>  AnalysisTool::firsttwo(Float_t *vec, vector<bool> *asso){

  double max(-999); int idmax(-999);
  double secondmax(-999); int idsecondmax(-999);

  for (int i=0; i<int(asso->size()); i++) {

    if ( vec[i] > max && asso->at(i)) {                             //save boolean = 1 case
      max = vec[i];
      idmax = i;
    }

  }
  for (int i=0; i<int(asso->size()); i++) {

    if ( vec[i] > secondmax && asso->at(i) && i!= idmax) { 
      secondmax = vec[i];
      idsecondmax = i;
    }

  }
  
  //  int themaxtemp[] = {idmax,idsecondmax};

  vector<int> themax;
  
  themax.push_back(idmax);
  themax.push_back(idsecondmax);

  return themax;

}



//PhotonID Del Re

bool AnalysisTool::cutID(int i, photonidcuts const& pid, vector<bool> *vpass) {

  // Use photon supercluster energy (would be e5x5 if r9>0.93 otherwise)
  bool ntrkiso = ntrkiso035Phot[i] < pid.tracknb;
  bool ptiso = (ptiso035Phot[i] / ptPhot[i] < pid.trackiso_rel);
  bool ecaliso = (ecaliso04Phot[i] / ePhot[i] < pid.ecaliso_rel ||
                  ecaliso04Phot[i] < pid.ecaliso_abs);
  double fhcal = hcalovecal04Phot[i];
  bool hcaliso = (fhcal < pid.hcaliso_rel ||
                  fhcal*ptPhot[i] < pid.hcaliso_abs);
  bool smaj = sMajMajPhot[i] < pid.smajmaj;
  bool smin = sMinMinPhot[i] < pid.sminmin;
  bool smin_min = sMinMinPhot[i] > pid.sminmin_min;
  bool eta = TMath::Abs(etaPhot[i]) < 2.5; 


  if(TMath::Abs(etaPhot[i]) > 1.44) {
    smaj = 1; smin = 1; smin_min = 1;
  }
  
  if (vpass) {
    // assert((*vpass).size()==7);
    (*vpass)[0] = ntrkiso;
    (*vpass)[1] = ptiso;
    (*vpass)[2] = hcaliso;
    (*vpass)[3] = ecaliso;
    (*vpass)[4] = smaj;
    (*vpass)[5] = smin;
    (*vpass)[6] = smin_min; 
  }

  return (ntrkiso && ptiso && hcaliso && ecaliso && smaj && smin && smin_min && eta);
}

