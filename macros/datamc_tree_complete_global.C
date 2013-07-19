inline double delta_phi(double phi1, double phi2) {

  double dphi = TMath::Abs(phi1 - phi2);
  return (dphi <= TMath::Pi())? dphi : TMath::TwoPi() - dphi;

}

// lumi runA 3.1pb-1 runB 30.9pb-1

Double_t gausconst(Double_t *x, Double_t *par)
{

  double exponent = 0.0;
  double gau1 = 0.0;

  exponent = (x[0] - par[0])/par[1];
  gau1 = exp(-exponent*exponent/2.);
  gau1 = gau1/(sqrt(2.*3.14)*par[1]);

  //  Double_t total = par[3] + par[4] * x[0] + par[2] * gau1;

  Double_t total = par[3] / exp(par[4] * x[0]) + par[2] * gau1;  

  return total;
}

void datamc_tree_complete_global(int rebin, double int_exp, double int_exp_2011, double pt1=50, double pt2=30, double ptj1=20, double ptj2=15, double deltae=2.5, double zep=2.5, double mjj=300,double isol1=1, double isol2=1, double isom1=0, double isom2=0){

  gROOT->SetStyle("Plain");
   
  gStyle->SetPalette(1);
  gStyle->SetOptStat(1111111);  // Show overflow, underflow + SumOfWeights 
  gStyle->SetOptFit(111110); 
  gStyle->SetOptFile(1); 
  
  gStyle->SetMarkerStyle(20);
  gStyle->SetMarkerSize(.3);
  gStyle->SetMarkerColor(1);

  c0 = new TCanvas("c1"," ",200,10,500,500);
  c0->Clear();

  //  double int_exp = 150;

  char name[100];
  sprintf(name,"higgsmassjustisocutrecofull");

  // LOOSE EG PRESEL
  // old
  //TFile data("redntp.V27.41x.eg/redntp_run2011.root");  
  //TFile data("redntp.V27.41x.eg/redntp_run2010-2011.root");  
  //TFile data("redntp.V27.41x.eg/redntp_run2011-upto163334.root");
  //TFile data("redntp.V27.41x.eg/redntp_run2010-2011-upto163334.root");
  // prompt
  //TFile data("redntp.V27.41x.eg/redntp_run2011_prompt0425.root");
  //TFile data("redntp.V27.41x.eg/redntp_run2010-2011_prompt0425.root");
  // 
  //TFile data("redntp.V27.39x.eg/redntp_run2010_time.root");
  TFile data("redntp.V27.41x.eg/redntp_run2011_upto-163340-ESfix_time.root");
  //TFile data("redntp.V27.41x.eg/redntp_run2010-2011_upto-163340-ESfix_time.root");
  // HGG presel
  //TFile data("redntp.V27.39x.hgg/redntp_run2010.root");
  //TFile data("redntp.V27.41x.eg/redntp_run2011_upto-163289-ESfix.root");
  //TFile data("redntp.V27.41x.hgg/redntp_run2010-2011_upto-163340-ESfix.root");

  TH1D* higgsmassdata = new TH1D("higgsmassdata","higgsmassdata", ((TH1D*)gDirectory->Get(name))->GetNbinsX(), ((TH1D*)gDirectory->Get(name))->GetBinLowEdge(1),((TH1D*)gDirectory->Get(name))->GetBinLowEdge(((TH1D*)gDirectory->Get(name))->GetNbinsX()+1));
  TH1D* higgsmassdatacs = new TH1D("higgsmassdatacs","higgsmassdatacs", ((TH1D*)gDirectory->Get(name))->GetNbinsX(), ((TH1D*)gDirectory->Get(name))->GetBinLowEdge(1),((TH1D*)gDirectory->Get(name))->GetBinLowEdge(((TH1D*)gDirectory->Get(name))->GetNbinsX()+1));
  TH1D* higgsmassdata2010 = new TH1D("higgsmassdata2010","higgsmassdata2010", ((TH1D*)gDirectory->Get(name))->GetNbinsX(), ((TH1D*)gDirectory->Get(name))->GetBinLowEdge(1),((TH1D*)gDirectory->Get(name))->GetBinLowEdge(((TH1D*)gDirectory->Get(name))->GetNbinsX()+1));


  TFile vbf("redntp.preegloose.V12/redntp_GluGluToHToGG_M-115_7TeV-powheg-pythia6_00.root");
  int n_vbf = ptphotgen1->GetEntries();
  cout << "# events vbf = " << n_vbf << endl;
  //  double cross_vbf = 1.3062 * 0.028;
  //  double cross_vbf = 18.23 * 2.13e-03 * 20;
  double cross_vbf = 18.23 * 2.13e-03 * 30;
  TH1D* higgsmassvbf = new TH1D("higgsmassvbf","higgsmassvbf", ((TH1D*)gDirectory->Get(name))->GetNbinsX(), ((TH1D*)gDirectory->Get(name))->GetBinLowEdge(1),((TH1D*)gDirectory->Get(name))->GetBinLowEdge(((TH1D*)gDirectory->Get(name))->GetNbinsX()+1));

  //  TFile vbf_2011("redntp.preegloose.V12/redntp_GluGluToHToGG_M-115_7TeV-powheg-pythia6_00.root");
  TFile vbf_2011("redntp.V27.41x.eg/redntp_GluGluToHToGG_M-115_7TeV-pythia6_00.root");
  int n_vbf_2011 = ptphotgen1->GetEntries();
  cout << "# events vbf_2011 = " << n_vbf_2011 << endl;
  //  double cross_vbf_2011 = 1.3062 * 0.028;
  //  double cross_vbf_2011 = 18.23 * 2.13e-03 * 20;
  double cross_vbf_2011 = 18.23 * 2.13e-03 * 30;
  TH1D* higgsmassvbf_2011 = new TH1D("higgsmassvbf_2011","higgsmassvbf_2011", ((TH1D*)gDirectory->Get(name))->GetNbinsX(), ((TH1D*)gDirectory->Get(name))->GetBinLowEdge(1),((TH1D*)gDirectory->Get(name))->GetBinLowEdge(((TH1D*)gDirectory->Get(name))->GetNbinsX()+1));

  TFile dy("redntp.V27.39x.eg/redntp_DYJetsToLL_TuneZ2_M-50_7TeV-madgraph-tauola.root");
  int n_dy = ptphotgen1->GetEntries();
  cout << "# events dy = " << n_dy << endl;
  double cross_dy = 2321.;
  TH1D* higgsmassdy = new TH1D("higgsmassdy","higgsmassdy", ((TH1D*)gDirectory->Get(name))->GetNbinsX(), ((TH1D*)gDirectory->Get(name))->GetBinLowEdge(1),((TH1D*)gDirectory->Get(name))->GetBinLowEdge(((TH1D*)gDirectory->Get(name))->GetNbinsX()+1));

  TFile dy_2011("redntp.V27.41x.eg/redntp_DYJetsToLL_TuneZ2_M-50_7TeV-madgraph-tauola.root");
  int n_dy_2011 = ptphotgen1->GetEntries();
  cout << "# events dy_2011 = " << n_dy_2011 << endl;
  double cross_dy_2011 = 2321.;
  TH1D* higgsmassdy_2011 = new TH1D("higgsmassdy_2011","higgsmassdy_2011", ((TH1D*)gDirectory->Get(name))->GetNbinsX(), ((TH1D*)gDirectory->Get(name))->GetBinLowEdge(1),((TH1D*)gDirectory->Get(name))->GetBinLowEdge(((TH1D*)gDirectory->Get(name))->GetNbinsX()+1));

  TFile box("redntp.V27.39x.eg/redntp_DiPhotonBox_Pt25to250_TrackingParticles_7TeV-pythia6_00.root");
  int n_box = ptphotgen1->GetEntries();
  cout << "# events box = " << n_box << endl;
  double cross_box = 12.37;
  TH1D* higgsmassbox = new TH1D("higgsmassbox","higgsmassbox", ((TH1D*)gDirectory->Get(name))->GetNbinsX(), ((TH1D*)gDirectory->Get(name))->GetBinLowEdge(1),((TH1D*)gDirectory->Get(name))->GetBinLowEdge(((TH1D*)gDirectory->Get(name))->GetNbinsX()+1));

  TFile box_2011("redntp.V27.41x.eg/redntp_DiPhotonBox_Pt25to250_TrackingParticles_7TeV-pythia6_00.root");
  int n_box_2011 = ptphotgen1->GetEntries();
  cout << "# events box_2011 = " << n_box_2011 << endl;
  double cross_box_2011 = 12.37;
  TH1D* higgsmassbox_2011 = new TH1D("higgsmassbox_2011","higgsmassbox_2011", ((TH1D*)gDirectory->Get(name))->GetNbinsX(), ((TH1D*)gDirectory->Get(name))->GetBinLowEdge(1),((TH1D*)gDirectory->Get(name))->GetBinLowEdge(((TH1D*)gDirectory->Get(name))->GetNbinsX()+1));

  TFile diphotjet("redntp.V27.39x.eg/redntp_DiPhotonJets_7TeV-madgraph.root");
  int n_diphotjet = ptphotgen1->GetEntries();
  cout << "# events diphotjet = " << n_diphotjet << endl;
  double cross_diphotjet = 134;
  TH1D* higgsmassdiphotjet = new TH1D("higgsmassdiphotjet","higgsmassdiphotjet", ((TH1D*)gDirectory->Get(name))->GetNbinsX(), ((TH1D*)gDirectory->Get(name))->GetBinLowEdge(1),((TH1D*)gDirectory->Get(name))->GetBinLowEdge(((TH1D*)gDirectory->Get(name))->GetNbinsX()+1));

  TFile diphotjet_2011("redntp.V27.39x.eg/redntp_DiPhotonJets_7TeV-madgraph.root");
  int n_diphotjet_2011 = ptphotgen1->GetEntries();
  cout << "# events diphotjet_2011 = " << n_diphotjet_2011 << endl;
  double cross_diphotjet_2011 = 134  * 2./3.;
  TH1D* higgsmassdiphotjet_2011 = new TH1D("higgsmassdiphotjet_2011","higgsmassdiphotjet_2011", ((TH1D*)gDirectory->Get(name))->GetNbinsX(), ((TH1D*)gDirectory->Get(name))->GetBinLowEdge(1),((TH1D*)gDirectory->Get(name))->GetBinLowEdge(((TH1D*)gDirectory->Get(name))->GetNbinsX()+1));

  TFile gjet("redntp.V27.39x.eg/redntp_GJet_Pt-20_doubleEMEnriched_TuneZ2_7TeV-pythia6.root");
  int n_gjet = ptphotgen1->GetEntries();
  cout << "# events gjet = " << n_gjet << endl;
  double cross_gjet = 493.44;
  TH1D* higgsmassgjet = new TH1D("higgsmassgjet","higgsmassgjet", ((TH1D*)gDirectory->Get(name))->GetNbinsX(), ((TH1D*)gDirectory->Get(name))->GetBinLowEdge(1),((TH1D*)gDirectory->Get(name))->GetBinLowEdge(((TH1D*)gDirectory->Get(name))->GetNbinsX()+1));

  TFile gjet_2011("redntp.V27.41x.eg/redntp_GJet_Pt-20_doubleEMEnriched_TuneZ2_7TeV-pythia6.root");
  int n_gjet_2011 = ptphotgen1->GetEntries();
  cout << "# events gjet_2011 = " << n_gjet_2011 << endl;
  double cross_gjet_2011 = 493.44;
  TH1D* higgsmassgjet_2011 = new TH1D("higgsmassgjet_2011","higgsmassgjet_2011", ((TH1D*)gDirectory->Get(name))->GetNbinsX(), ((TH1D*)gDirectory->Get(name))->GetBinLowEdge(1),((TH1D*)gDirectory->Get(name))->GetBinLowEdge(((TH1D*)gDirectory->Get(name))->GetNbinsX()+1));

  TFile qcd("redntp.V27.39x.eg/redntp_QCD_Pt-40_doubleEMEnriched_TuneZ2_7TeV-pythia6.root");
  int n_qcd = ptphotgen1->GetEntries();
  cout << "# events qcd = " << n_qcd << endl;
  double cross_qcd = 40392;
  TH1D* higgsmassqcd = new TH1D("higgsmassqcd","higgsmassqcd", ((TH1D*)gDirectory->Get(name))->GetNbinsX(), ((TH1D*)gDirectory->Get(name))->GetBinLowEdge(1),((TH1D*)gDirectory->Get(name))->GetBinLowEdge(((TH1D*)gDirectory->Get(name))->GetNbinsX()+1));

  TFile qcd_2011("redntp.V27.41x.eg/redntp_QCD_Pt-40_doubleEMEnriched_TuneZ2_7TeV-pythia6.root");
  int n_qcd_2011 = ptphotgen1->GetEntries();
  cout << "# events qcd_2011 = " << n_qcd_2011 << endl;
  double cross_qcd_2011 = 40392;
  TH1D* higgsmassqcd_2011 = new TH1D("higgsmassqcd_2011","higgsmassqcd_2011", ((TH1D*)gDirectory->Get(name))->GetNbinsX(), ((TH1D*)gDirectory->Get(name))->GetBinLowEdge(1),((TH1D*)gDirectory->Get(name))->GetBinLowEdge(((TH1D*)gDirectory->Get(name))->GetNbinsX()+1));

  TFile qcd2("redntp.V27.39x.eg/redntp_QCD_Pt-30to40_doubleEMEnriched_TuneZ2_7TeV-pythia6.root");
  int n_qcd2 = ptphotgen1->GetEntries();
  cout << "# events qcd2 = " << n_qcd2 << endl;
  double cross_qcd2 = 9610;
  TH1D* higgsmassqcd2 = new TH1D("higgsmassqcd2","higgsmassqcd2", ((TH1D*)gDirectory->Get(name))->GetNbinsX(), ((TH1D*)gDirectory->Get(name))->GetBinLowEdge(1),((TH1D*)gDirectory->Get(name))->GetBinLowEdge(((TH1D*)gDirectory->Get(name))->GetNbinsX()+1));

  TFile qcd2_2011("redntp.V27.41x.eg/redntp_QCD_Pt-30to40_doubleEMEnriched_TuneZ2_7TeV-pythia6.root");
  int n_qcd2_2011 = ptphotgen1->GetEntries();
  cout << "# events qcd2_2011 = " << n_qcd2_2011 << endl;
  double cross_qcd2_2011 = 9610;
  TH1D* higgsmassqcd2_2011 = new TH1D("higgsmassqcd2_2011","higgsmassqcd2_2011", ((TH1D*)gDirectory->Get(name))->GetNbinsX(), ((TH1D*)gDirectory->Get(name))->GetBinLowEdge(1),((TH1D*)gDirectory->Get(name))->GetBinLowEdge(((TH1D*)gDirectory->Get(name))->GetNbinsX()+1));

  char allcut[3000], allcut2[3000];
  char cut1[150],cut2[150],cut3[150],cut4[150],cut5[150],cut6[150],cut7[150],cut8[150],cut9[150],cut9[150],cut10[150],cut11[150], cutextra[200], cut100[200];
  sprintf(cut1,"%s%f","ptphot1>",pt1);
  sprintf(cut2,"%s%f","ptphot2>",pt2);
  sprintf(cut3,"%s%f","ptjet1>",ptj1);
  sprintf(cut4,"%s%f","ptjet2>",ptj2);
  sprintf(cut5,"%s%f","abs(deltaeta)>",deltae);
  sprintf(cut6,"%s%f","abs(zeppenjet)<",zep);
  sprintf(cut7,"%s%f","invmassjet>",mjj);
  // TEMP
  sprintf(cut100,"ptgg<1000000000");
  //printf(cut100,"ptjet1<15");
  //sprintf(cut100,"nvtx<5");
  //sprintf(cut100,"delta_phi(phiphot1,phiphot2)>3.");
  //sprintf(cut100,"!(abs(etaphot1)<1.4442&&abs(etaphot2)<1.4442)&&delta_phi(phiphot1,phiphot2)>3.");
  //sprintf(cut100,"!(abs(etaphot1)<1.4442&&abs(etaphot2)<1.4442)&&ptjet1<20");
  //sprintf(cut100,"!(abs(etaphot1)<1.4442&&abs(etaphot2)<1.4442)");
  // END TEMP
  if(isol1 && isol2) sprintf(cutextra,"((idloosephot1==0&&idloosephot2)||(idloosephot1&&idloosephot2==0))");
  if(isol1 && isom2) sprintf(cutextra,"((idloosephot1==0&&idmediumphot2)||(idloosephot1&&idmediumphot2==0))");
  if(isom1 && isol2) sprintf(cutextra,"((idmediumphot1==0&&idloosephot2)||(idmediumphot1&&idloosephot2==0))");
  if(isom1 && isom2) sprintf(cutextra,"((idmediumphot1==0&&idmediumphot2)||(idmediumphot1&&idmediumphot2==0))");
  if(isol1)  sprintf(cut8,"idtightnewEGphot1&&!pid_haspixelseedphot1");
  else   sprintf(cut8,"massgg>0");
  if(isol2)  sprintf(cut9,"idtightnewEGphot2&&!pid_haspixelseedphot2");
  else   sprintf(cut9,"massgg>0");
  if(isom1)  sprintf(cut10,"idhggtightnewEGphot1&&!pid_haspixelseedphot1");
  else   sprintf(cut10,"massgg>0");
  if(isom2)  sprintf(cut11,"idhggtightnewEGphot2&&!pid_haspixelseedphot2");
  else   sprintf(cut11,"massgg>0");
  sprintf(allcut,"%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s",cut1,"&&!(abs(etaphot1)>1.4442&&abs(etaphot1)<1.566)&&!(abs(etaphot2)>1.4442&&abs(etaphot2)<1.566)&&abs(etaphot1)<2.5&&abs(etaphot2)<2.5&&",cut2,"&&",cut3,"&&",cut4,"&&",cut5,"&&",cut6,"&&",cut7,"&&",cut8,"&&",cut9,"&&",cut10,"&&",cut11,"&&",cut100);
  sprintf(allcut2,"%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s",cut1,"&&!(abs(etaphot1)>1.4442&&abs(etaphot1)<1.566)&&!(abs(etaphot2)>1.4442&&abs(etaphot2)<1.566)&&abs(etaphot1)<2.5&&abs(etaphot2)<2.5&&",cut2,"&&",cut3,"&&",cut4,"&&",cut5,"&&",cut6,"&&",cut7,"&&",cut100,"&&",cutextra);

  cout << allcut << endl;
  cout << allcut2 << endl;
  
  data.cd();
  char tempcut[10000]; 
  double scale = 1./rebin;
  sprintf(name,"higgsmassdata");
  AnaTree->Project(name,"massgg",allcut);
  ((TH1D*)gDirectory->Get(name))->Sumw2();
  ((TH1D*)gDirectory->Get(name))->Rebin(rebin);
  double num_data = ((TH1D*)gDirectory->Get(name))->Integral(101,150); 
  sprintf(name,"higgsmassdata2010");
  sprintf(tempcut,"%s%s",allcut,"&&run<160000");
  AnaTree->Project(name,"massgg",tempcut);
  ((TH1D*)gDirectory->Get(name))->Rebin(rebin);
  sprintf(tempcut,"%s%s",allcut,"&&massgg>114&&massgg<119");
  AnaTree->Scan("ptphot1:etaphot1:timephot1:pid_hlwTrackphot1:pid_etawidphot1:pid_HoverEphot1:pid_jurECALphot1:pid_twrHCALphot1:ptphot2:etaphot2:timephot2:pid_hlwTrackphot2:pid_etawidphot2:pid_jurECALphot2:pid_twrHCALphot2:pid_HoverEphot2",tempcut);
  TH1D timesigphot1("timesigphot1","",20,-2.,2.);
  TH1D timesigphot2("timesigphot2","",20,-2.,2.);
  AnaTree->Draw("timephot1>>timesigphot1",tempcut); 
  AnaTree->Draw("timephot2>>timesigphot2",tempcut); 
  timesigphot1.Add(&timesigphot2);
  timesigphot1.Draw("pe");
  c0->SaveAs("timesigphot.gif");
  TH1D timebkgphot1("timebkgphot1","",20,-2.,2.);
  TH1D timebkgphot2("timebkgphot2","",20,-2.,2.);
   sprintf(tempcut,"%s%s",allcut,"&&!(massgg>114&&massgg<119)");
  AnaTree->Draw("timephot1>>timebkgphot1",tempcut); 
  AnaTree->Draw("timephot2>>timebkgphot2",tempcut); 
  timebkgphot1.Add(&timebkgphot2);
  timebkgphot1.Draw("pe");
  c0->SaveAs("timebkgphot.gif");
  sprintf(tempcut,"%s%s",allcut,"&&massgg>114&&massgg<119");
  double entries_sb = AnaTree->GetEntries(tempcut);  
  AnaTree->Scan("ptphot1:ptphot2:ptjet1:ptjet2:met:nvtx:massgg:run:lumi:event",tempcut);
  cout << endl; cout << endl;
  sprintf(tempcut,"%s%s",allcut,"&&massgg>100&&massgg<150");
  //  AnaTree->Scan("ptphot1:etaphot1:timephot1:pid_hlwTrackphot1:pid_etawidphot1:pid_HoverEphot1:pid_jurECALphot1:pid_twrHCALphot1:ptphot2:etaphot2:timephot2:pid_hlwTrackphot2:pid_etawidphot2:pid_jurECALphot2:pid_twrHCALphot2:pid_HoverEphot2",tempcut);
  AnaTree->Scan("ptphot1:ptphot2:massgg:run:lumi:event",tempcut);
  sprintf(name,"higgsmassdatacs");
  AnaTree->Project(name,"massgg",allcut2);
  sprintf(tempcut,"%s%s",allcut2,"&&massgg>123&&massgg<131");
  double entries_sb_cs = AnaTree->GetEntries(tempcut);  
  ((TH1D*)gDirectory->Get(name))->Rebin(rebin);

  vbf.cd();
  scale = cross_vbf * int_exp / n_vbf;
  sprintf(name,"higgsmassvbf");
  AnaTree->Project(name,"massgg",allcut);
  double num_vbf = ((TH1D*)gDirectory->Get(name))->Integral(101,150)*scale;
  double num_vbf_unscaled = ((TH1D*)gDirectory->Get(name))->Integral(101,150);
  cout << "unscaled " << num_vbf_unscaled <<endl;
  ((TH1D*)gDirectory->Get(name))->Rebin(rebin);
  ((TH1D*)gDirectory->Get(name))->Scale(scale);

  vbf_2011.cd();
  scale = cross_vbf_2011 * int_exp_2011 / n_vbf_2011;
  sprintf(name,"higgsmassvbf_2011");
  AnaTree->Project(name,"massgg",allcut);
  double num_vbf_2011 = ((TH1D*)gDirectory->Get(name))->Integral(101,150)*scale;
  double num_vbf_2011_unscaled = ((TH1D*)gDirectory->Get(name))->Integral(101,150);
  cout << "unscaled " << num_vbf_2011_unscaled <<endl;
  ((TH1D*)gDirectory->Get(name))->Rebin(rebin);
  ((TH1D*)gDirectory->Get(name))->Scale(scale);

  dy.cd();
  scale = cross_dy * int_exp / n_dy;
  sprintf(name,"higgsmassdy");
  AnaTree->Project(name,"massgg",allcut);
  double num_dy = ((TH1D*)gDirectory->Get(name))->Integral(101,150)*scale;
  double num_dy_unscaled = ((TH1D*)gDirectory->Get(name))->Integral(101,150);
  cout << "unscaled " << num_dy_unscaled <<endl;
  ((TH1D*)gDirectory->Get(name))->Rebin(rebin);
  ((TH1D*)gDirectory->Get(name))->Scale(scale);

  dy_2011.cd();
  scale = cross_dy_2011 * int_exp_2011 / n_dy_2011;
  sprintf(name,"higgsmassdy_2011");
  AnaTree->Project(name,"massgg",allcut);
  double num_dy_2011 = ((TH1D*)gDirectory->Get(name))->Integral(101,150)*scale;
  double num_dy_2011_unscaled = ((TH1D*)gDirectory->Get(name))->Integral(101,150);
  cout << "unscaled " << num_dy_2011_unscaled <<endl;
  ((TH1D*)gDirectory->Get(name))->Rebin(rebin);
  ((TH1D*)gDirectory->Get(name))->Scale(scale);

  box.cd();
  scale = cross_box * int_exp / n_box;
  sprintf(name,"higgsmassbox");
  AnaTree->Project(name,"massgg",allcut);
  double num_box = ((TH1D*)gDirectory->Get(name))->Integral(101,150)*scale;
  double num_box_unscaled = ((TH1D*)gDirectory->Get(name))->Integral(101,150);
  ((TH1D*)gDirectory->Get(name))->Rebin(rebin);
  ((TH1D*)gDirectory->Get(name))->Scale(scale);
  ((TH1D*)gDirectory->Get(name))->SetLineColor(kRed);

  box_2011.cd();
  scale = cross_box_2011 * int_exp_2011 / n_box_2011;
  sprintf(name,"higgsmassbox_2011");
  AnaTree->Project(name,"massgg",allcut);
  double num_box_2011 = ((TH1D*)gDirectory->Get(name))->Integral(101,150)*scale;
  double num_box_2011_unscaled = ((TH1D*)gDirectory->Get(name))->Integral(101,150);
  ((TH1D*)gDirectory->Get(name))->Rebin(rebin);
  ((TH1D*)gDirectory->Get(name))->Scale(scale);
  ((TH1D*)gDirectory->Get(name))->SetLineColor(kRed);

  diphotjet.cd();
  scale = cross_diphotjet * int_exp / n_diphotjet;
  sprintf(name,"higgsmassdiphotjet");
  AnaTree->Project(name,"massgg",allcut);
  double num_diphotjet = ((TH1D*)gDirectory->Get(name))->Integral(101,150)*scale;
  double num_diphotjet_unscaled = ((TH1D*)gDirectory->Get(name))->Integral(101,150);
  ((TH1D*)gDirectory->Get(name))->Rebin(rebin);
  ((TH1D*)gDirectory->Get(name))->Scale(scale);
  ((TH1D*)gDirectory->Get(name))->SetLineColor(kBlue);

  diphotjet_2011.cd();
  scale = cross_diphotjet_2011 * int_exp_2011 / n_diphotjet_2011;
  sprintf(name,"higgsmassdiphotjet_2011");
  AnaTree->Project(name,"massgg",allcut);
  double num_diphotjet_2011 = ((TH1D*)gDirectory->Get(name))->Integral(101,150)*scale;
  double num_diphotjet_2011_unscaled = ((TH1D*)gDirectory->Get(name))->Integral(101,150);
  ((TH1D*)gDirectory->Get(name))->Rebin(rebin);
  ((TH1D*)gDirectory->Get(name))->Scale(scale);
  ((TH1D*)gDirectory->Get(name))->SetLineColor(kBlue);

  gjet.cd();
  scale = cross_gjet * int_exp / n_gjet;
  sprintf(name,"higgsmassgjet");
  AnaTree->Project(name,"massgg",allcut);
  double num_gjet = ((TH1D*)gDirectory->Get(name))->Integral(101,150)*scale;
  double num_gjet_unscaled = ((TH1D*)gDirectory->Get(name))->Integral(101,150);
  ((TH1D*)gDirectory->Get(name))->Rebin(rebin);
  ((TH1D*)gDirectory->Get(name))->Scale(scale);
  ((TH1D*)gDirectory->Get(name))->SetLineColor(kMagenta);

  gjet_2011.cd();
  scale = cross_gjet_2011 * int_exp_2011 / n_gjet_2011;
  sprintf(name,"higgsmassgjet_2011");
  AnaTree->Project(name,"massgg",allcut);
  double num_gjet_2011 = ((TH1D*)gDirectory->Get(name))->Integral(101,150)*scale;
  double num_gjet_2011_unscaled = ((TH1D*)gDirectory->Get(name))->Integral(101,150);
  ((TH1D*)gDirectory->Get(name))->Rebin(rebin);
  ((TH1D*)gDirectory->Get(name))->Scale(scale);
  ((TH1D*)gDirectory->Get(name))->SetLineColor(kMagenta);

  qcd.cd();
  scale = cross_qcd * int_exp / n_qcd;
  sprintf(name,"higgsmassqcd");
  AnaTree->Project(name,"massgg",allcut);
  double num_qcd = ((TH1D*)gDirectory->Get(name))->Integral(101,150)*scale;
  double num_qcd_unscaled = ((TH1D*)gDirectory->Get(name))->Integral(101,150);
  ((TH1D*)gDirectory->Get(name))->Rebin(rebin);
  ((TH1D*)gDirectory->Get(name))->Scale(scale);
  ((TH1D*)gDirectory->Get(name))->SetLineColor(kGreen);

  qcd_2011.cd();
  scale = cross_qcd_2011 * int_exp_2011 / n_qcd_2011;
  sprintf(name,"higgsmassqcd_2011");
  AnaTree->Project(name,"massgg",allcut);
  double num_qcd_2011 = ((TH1D*)gDirectory->Get(name))->Integral(101,150)*scale;
  double num_qcd_2011_unscaled = ((TH1D*)gDirectory->Get(name))->Integral(101,150);
  ((TH1D*)gDirectory->Get(name))->Rebin(rebin);
  ((TH1D*)gDirectory->Get(name))->Scale(scale);
  ((TH1D*)gDirectory->Get(name))->SetLineColor(kGreen);

  qcd2.cd();
  scale = cross_qcd2 * int_exp / n_qcd2;
  sprintf(name,"higgsmassqcd2");
  AnaTree->Project(name,"massgg",allcut);
  double num_qcd2 = ((TH1D*)gDirectory->Get(name))->Integral(101,150)*scale;
  double num_qcd2_unscaled = ((TH1D*)gDirectory->Get(name))->Integral(101,150);
  ((TH1D*)gDirectory->Get(name))->Rebin(rebin);
  ((TH1D*)gDirectory->Get(name))->Scale(scale);
  ((TH1D*)gDirectory->Get(name))->SetLineColor(kGreen);
  higgsmassqcd->Add(((TH1D*)gDirectory->Get(name)));

  qcd2_2011.cd();
  scale = cross_qcd2_2011 * int_exp_2011 / n_qcd2_2011;
  sprintf(name,"higgsmassqcd2_2011");
  AnaTree->Project(name,"massgg",allcut);
  double num_qcd2_2011 = ((TH1D*)gDirectory->Get(name))->Integral(101,150)*scale;
  double num_qcd2_2011_unscaled = ((TH1D*)gDirectory->Get(name))->Integral(101,150);
  ((TH1D*)gDirectory->Get(name))->Rebin(rebin);
  ((TH1D*)gDirectory->Get(name))->Scale(scale);
  ((TH1D*)gDirectory->Get(name))->SetLineColor(kGreen);
  higgsmassqcd->Add(((TH1D*)gDirectory->Get(name)));


  vbf.cd();
  sprintf(name,"higgsmassvbf");
  TH1D* higgsmass1 = new TH1D("higgsmass1","higgsmass1", ((TH1D*)gDirectory->Get(name))->GetNbinsX(), ((TH1D*)gDirectory->Get(name))->GetBinLowEdge(1),((TH1D*)gDirectory->Get(name))->GetBinLowEdge(((TH1D*)gDirectory->Get(name))->GetNbinsX()+1));
  TH1D* higgsmass2 = new TH1D("higgsmass2","higgsmass2", ((TH1D*)gDirectory->Get(name))->GetNbinsX(), ((TH1D*)gDirectory->Get(name))->GetBinLowEdge(1),((TH1D*)gDirectory->Get(name))->GetBinLowEdge(((TH1D*)gDirectory->Get(name))->GetNbinsX()+1));
  TH1D* higgsmass3 = new TH1D("higgsmass3","higgsmass3", ((TH1D*)gDirectory->Get(name))->GetNbinsX(), ((TH1D*)gDirectory->Get(name))->GetBinLowEdge(1),((TH1D*)gDirectory->Get(name))->GetBinLowEdge(((TH1D*)gDirectory->Get(name))->GetNbinsX()+1));
  TH1D* higgsmass4 = new TH1D("higgsmass4","higgsmass4", ((TH1D*)gDirectory->Get(name))->GetNbinsX(), ((TH1D*)gDirectory->Get(name))->GetBinLowEdge(1),((TH1D*)gDirectory->Get(name))->GetBinLowEdge(((TH1D*)gDirectory->Get(name))->GetNbinsX()+1));
  TH1D* higgsmass5 = new TH1D("higgsmass5","higgsmass5", ((TH1D*)gDirectory->Get(name))->GetNbinsX(), ((TH1D*)gDirectory->Get(name))->GetBinLowEdge(1),((TH1D*)gDirectory->Get(name))->GetBinLowEdge(((TH1D*)gDirectory->Get(name))->GetNbinsX()+1));
  TH1D* higgsmass6 = new TH1D("higgsmass6","higgsmass6", ((TH1D*)gDirectory->Get(name))->GetNbinsX(), ((TH1D*)gDirectory->Get(name))->GetBinLowEdge(1),((TH1D*)gDirectory->Get(name))->GetBinLowEdge(((TH1D*)gDirectory->Get(name))->GetNbinsX()+1));

  for (int j=0; j<6; j++){

     for (int i=1; i<((TH1D*)gDirectory->Get(name))->GetNbinsX()+1; i++){      
      
      int k = i;
      // if(j!=3) k=int((i-1)/10.+1);      
      
      if(j==0) {
	box.cd();
	sprintf(name,"higgsmassbox");
	higgsmass1->SetBinContent(i,((TH1D*)gDirectory->Get(name))->GetBinContent(k) + higgsmass1->GetBinContent(i));
	higgsmass2->SetBinContent(i,((TH1D*)gDirectory->Get(name))->GetBinContent(k) + higgsmass2->GetBinContent(i));
	higgsmass3->SetBinContent(i,((TH1D*)gDirectory->Get(name))->GetBinContent(k) + higgsmass3->GetBinContent(i));
	higgsmass4->SetBinContent(i,((TH1D*)gDirectory->Get(name))->GetBinContent(k) + higgsmass4->GetBinContent(i));
	higgsmass5->SetBinContent(i,((TH1D*)gDirectory->Get(name))->GetBinContent(k) + higgsmass5->GetBinContent(i));
	higgsmass6->SetBinContent(i,((TH1D*)gDirectory->Get(name))->GetBinContent(k) + higgsmass6->GetBinContent(i));
 	box_2011.cd();
	sprintf(name,"higgsmassbox_2011");
	higgsmass1->SetBinContent(i,((TH1D*)gDirectory->Get(name))->GetBinContent(k) + higgsmass1->GetBinContent(i));
	higgsmass2->SetBinContent(i,((TH1D*)gDirectory->Get(name))->GetBinContent(k) + higgsmass2->GetBinContent(i));
	higgsmass3->SetBinContent(i,((TH1D*)gDirectory->Get(name))->GetBinContent(k) + higgsmass3->GetBinContent(i));
	higgsmass4->SetBinContent(i,((TH1D*)gDirectory->Get(name))->GetBinContent(k) + higgsmass4->GetBinContent(i));
	higgsmass5->SetBinContent(i,((TH1D*)gDirectory->Get(name))->GetBinContent(k) + higgsmass5->GetBinContent(i));
	higgsmass6->SetBinContent(i,((TH1D*)gDirectory->Get(name))->GetBinContent(k) + higgsmass6->GetBinContent(i));
      }
      if(j==1) {
	diphotjet.cd();
	sprintf(name,"higgsmassdiphotjet");
	higgsmass2->SetBinContent(i,((TH1D*)gDirectory->Get(name))->GetBinContent(k) + higgsmass2->GetBinContent(i));
	higgsmass3->SetBinContent(i,((TH1D*)gDirectory->Get(name))->GetBinContent(k) + higgsmass3->GetBinContent(i));
	higgsmass4->SetBinContent(i,((TH1D*)gDirectory->Get(name))->GetBinContent(k) + higgsmass4->GetBinContent(i));
	higgsmass5->SetBinContent(i,((TH1D*)gDirectory->Get(name))->GetBinContent(k) + higgsmass5->GetBinContent(i));
	higgsmass6->SetBinContent(i,((TH1D*)gDirectory->Get(name))->GetBinContent(k) + higgsmass6->GetBinContent(i));
 	diphotjet_2011.cd();
	sprintf(name,"higgsmassdiphotjet_2011");
	higgsmass2->SetBinContent(i,((TH1D*)gDirectory->Get(name))->GetBinContent(k) + higgsmass2->GetBinContent(i));
	higgsmass3->SetBinContent(i,((TH1D*)gDirectory->Get(name))->GetBinContent(k) + higgsmass3->GetBinContent(i));
	higgsmass4->SetBinContent(i,((TH1D*)gDirectory->Get(name))->GetBinContent(k) + higgsmass4->GetBinContent(i));
	higgsmass5->SetBinContent(i,((TH1D*)gDirectory->Get(name))->GetBinContent(k) + higgsmass5->GetBinContent(i));
	higgsmass6->SetBinContent(i,((TH1D*)gDirectory->Get(name))->GetBinContent(k) + higgsmass6->GetBinContent(i));
      }
      if(j==2) {
	gjet.cd();
	sprintf(name,"higgsmassgjet");
	higgsmass3->SetBinContent(i,((TH1D*)gDirectory->Get(name))->GetBinContent(k) + higgsmass3->GetBinContent(i));
	higgsmass4->SetBinContent(i,((TH1D*)gDirectory->Get(name))->GetBinContent(k) + higgsmass4->GetBinContent(i));
	higgsmass5->SetBinContent(i,((TH1D*)gDirectory->Get(name))->GetBinContent(k) + higgsmass5->GetBinContent(i));
	higgsmass6->SetBinContent(i,((TH1D*)gDirectory->Get(name))->GetBinContent(k) + higgsmass6->GetBinContent(i));
	gjet_2011.cd();
	sprintf(name,"higgsmassgjet_2011");
	higgsmass3->SetBinContent(i,((TH1D*)gDirectory->Get(name))->GetBinContent(k) + higgsmass3->GetBinContent(i));
	higgsmass4->SetBinContent(i,((TH1D*)gDirectory->Get(name))->GetBinContent(k) + higgsmass4->GetBinContent(i));
	higgsmass5->SetBinContent(i,((TH1D*)gDirectory->Get(name))->GetBinContent(k) + higgsmass5->GetBinContent(i));
	higgsmass6->SetBinContent(i,((TH1D*)gDirectory->Get(name))->GetBinContent(k) + higgsmass6->GetBinContent(i));
      }
      if(j==3) {
	qcd.cd();
	sprintf(name,"higgsmassqcd");
	higgsmass4->SetBinContent(i,((TH1D*)gDirectory->Get(name))->GetBinContent(k) + higgsmass4->GetBinContent(i));
	higgsmass5->SetBinContent(i,((TH1D*)gDirectory->Get(name))->GetBinContent(k) + higgsmass5->GetBinContent(i));
	higgsmass6->SetBinContent(i,((TH1D*)gDirectory->Get(name))->GetBinContent(k) + higgsmass6->GetBinContent(i));
	qcd_2011.cd();
	sprintf(name,"higgsmassqcd_2011");
	higgsmass4->SetBinContent(i,((TH1D*)gDirectory->Get(name))->GetBinContent(k) + higgsmass4->GetBinContent(i));
	higgsmass5->SetBinContent(i,((TH1D*)gDirectory->Get(name))->GetBinContent(k) + higgsmass5->GetBinContent(i));
	higgsmass6->SetBinContent(i,((TH1D*)gDirectory->Get(name))->GetBinContent(k) + higgsmass6->GetBinContent(i));
      }
      if(j==4) {
	dy.cd();
	sprintf(name,"higgsmassdy");
	higgsmass5->SetBinContent(i,((TH1D*)gDirectory->Get(name))->GetBinContent(k) + higgsmass5->GetBinContent(i));
	higgsmass6->SetBinContent(i,((TH1D*)gDirectory->Get(name))->GetBinContent(k) + higgsmass6->GetBinContent(i));
 	dy_2011.cd();
	sprintf(name,"higgsmassdy_2011");
	higgsmass5->SetBinContent(i,((TH1D*)gDirectory->Get(name))->GetBinContent(k) + higgsmass5->GetBinContent(i));
	higgsmass6->SetBinContent(i,((TH1D*)gDirectory->Get(name))->GetBinContent(k) + higgsmass6->GetBinContent(i));
      }
      if(j==5) {
	vbf.cd();
	sprintf(name,"higgsmassvbf");
	higgsmass6->SetBinContent(i,((TH1D*)gDirectory->Get(name))->GetBinContent(k) + higgsmass6->GetBinContent(i));
	vbf_2011.cd();
	sprintf(name,"higgsmassvbf_2011");
	higgsmass6->SetBinContent(i,((TH1D*)gDirectory->Get(name))->GetBinContent(k) + higgsmass6->GetBinContent(i));
      }
	
    }
    
    
  }
 
  higgsmass6->SetAxisRange(100.,150.);
  higgsmass6->SetTitle("");
  higgsmass6->SetStats(0);
  higgsmass6->SetTitleOffset(1.25,"Y");
  higgsmass6->SetXTitle("m(#gamma#gamma) [GeV]");
  char ytitle[100];
  sprintf(ytitle,"%s%d%s","N_{ev}/",int_exp,"pb^{-1}");
  higgsmass6->SetYTitle(ytitle);
  higgsmass6->SetLineColor(kBlack);
  higgsmass6->SetLineWidth(2);
  higgsmass6->SetFillColor(kYellow);
  higgsmass6->SetMinimum(0);
  higgsmass6->Draw();

  higgsmass5->SetLineColor(kBlack);
  higgsmass5->SetLineWidth(2);
  higgsmass5->SetFillColor(32);
  higgsmass5->Draw("same");

  higgsmass4->SetLineColor(kBlack);
  higgsmass4->SetLineWidth(2);
  higgsmass4->SetFillColor(29);
  higgsmass4->Draw("same");

  higgsmass3->SetLineColor(kBlack);
  higgsmass3->SetLineWidth(2);
  higgsmass3->SetFillColor(38);
  higgsmass3->Draw("same");

  higgsmass2->SetLineColor(kBlack);
  higgsmass2->SetLineWidth(2);
  higgsmass2->SetFillColor(46);
  higgsmass2->Draw("same");

  higgsmass1->SetLineColor(kBlack);
  higgsmass1->SetLineWidth(2);
  higgsmass1->SetFillColor(16);
  higgsmass1->Draw("same");

  TLegendEntry *legge;
  TLegend *leg;
  leg = new TLegend(0.6,0.6,0.85,0.85);
  leg->SetFillStyle(0); leg->SetBorderSize(0); leg->SetTextSize(0.05);
  leg->SetFillColor(0);
  legge = leg->AddEntry(higgsmass6, "h_{gg}(115)", "f");
  legge = leg->AddEntry(higgsmass5, "DY", "f");
  legge = leg->AddEntry(higgsmass4, "QCD", "f");
  legge = leg->AddEntry(higgsmass3, "#gamma + jets", "f");
  legge = leg->AddEntry(higgsmass2, "di-#gamma + jets", "f");
  legge = leg->AddEntry(higgsmass1, "di-#gamma box", "f");
  leg->Draw();

  sprintf(allcut,"%f%s%f%s%f%s%f%s%f%s%f%s%f%s%d%s%d%s%d%s%d",pt1,"_",pt2,"_",ptj1,"_",ptj2,"_",deltae,"_",zep,"_",mjj,"_",isol1,"_",isol2,"_",isom1,"_",isom2);
  char newname[1000];
  sprintf(newname,"%s%s%s","results_gg/mc_",allcut,".gif");
  c0->SaveAs(newname);

  data.cd();
  sprintf(name,"higgsmassdata");
  ((TH1D*)gDirectory->Get(name))->SetXTitle("m(#gamma#gamma) [GeV]");
  //  ((TH1D*)gDirectory->Get(name))->SetAxisRange(90.,299.);
  c0->SetLogy();
  ((TH1D*)gDirectory->Get(name))->SetAxisRange(89.,999.);
  ((TH1D*)gDirectory->Get(name))->SetTitle("");
  ((TH1D*)gDirectory->Get(name))->SetStats(0);
  ((TH1D*)gDirectory->Get(name))->SetMarkerStyle(8);
  ((TH1D*)gDirectory->Get(name))->SetMarkerSize(.9);
  ((TH1D*)gDirectory->Get(name))->SetTitleOffset(1.25,"Y");
  ((TH1D*)gDirectory->Get(name))->SetMinimum(0.002);
  ((TH1D*)gDirectory->Get(name))->Draw("pe");
  sprintf(newname,"%s%s%s","results_gg/data_",allcut,".gif");
  c0->SaveAs(newname);

  c0->SetLogy(0);

  ((TH1D*)gDirectory->Get(name))->SetAxisRange(90.,180.);
  ((TH1D*)gDirectory->Get(name))->Draw("pe");
  sprintf(newname,"%s%s%s","results_gg/data_",allcut,".gif");
  c0->SaveAs(newname);

  ((TH1D*)gDirectory->Get(name))->SetAxisRange(80.,150.);
  ((TH1D*)gDirectory->Get(name))->Draw("pe");
  sprintf(newname,"%s%s%s","results_gg/data_",allcut,"_2.gif");
  c0->SaveAs(newname);

  double norm = ((TH1D*)gDirectory->Get(name))->Integral()-entries_sb;
  sprintf(name,"higgsmassdatacs");
  ((TH1D*)gDirectory->Get(name))->SetAxisRange(90.,299.);
  double normcs = ((TH1D*)gDirectory->Get(name))->Integral()-entries_sb_cs;
  cout << norm << "  " << normcs << endl;
  ((TH1D*)gDirectory->Get(name))->Scale(norm/normcs); 
  ((TH1D*)gDirectory->Get(name))->SetLineColor(46);
  ((TH1D*)gDirectory->Get(name))->SetFillColor(42);
  ((TH1D*)gDirectory->Get(name))->SetLineWidth(3);
  ((TH1D*)gDirectory->Get(name))->Draw("same");
  sprintf(name,"higgsmassdata");
  ((TH1D*)gDirectory->Get(name))->Draw("pesame");
  sprintf(newname,"%s%s%s","results_gg/datacs_",allcut,".gif");
  gPad->RedrawAxis();
  TLegendEntry *legge2;
  TLegend *leg2;
  leg2 = new TLegend(0.6,0.65,0.9,0.85);
  leg2->SetFillStyle(0); leg2->SetBorderSize(0); leg2->SetTextSize(0.035);
  leg2->SetFillColor(0);
  legge2 = leg2->AddEntry(((TH1D*)gDirectory->Get("higgsmassdata")), "default sel.", "p");
  legge2 = leg2->AddEntry(((TH1D*)gDirectory->Get("higgsmassdatacs")), "control sample", "f");
  leg2->Draw();
  c0->SaveAs(newname);



  ((TH1D*)gDirectory->Get(name))->Draw("pe");
  higgsmass6->Draw("same");
  higgsmass5->Draw("same");
  higgsmass4->Draw("same");
  higgsmass3->Draw("same");
  higgsmass2->Draw("same");
  higgsmass1->Draw("same");
  leg->Draw();
  ((TH1D*)gDirectory->Get(name))->Draw("pesame");

  gPad->RedrawAxis();

  sprintf(newname,"%s%s%s","results_gg/data-mc_",allcut,".gif");
  c0->SaveAs(newname);

  c0->SetLogy(0);
  ((TH1D*)gDirectory->Get(name))->SetAxisRange(90.,180.);
  ((TH1D*)gDirectory->Get(name))->Draw("pe");
  higgsmass6->Draw("same");
  higgsmass5->Draw("same");
  higgsmass4->Draw("same");
  higgsmass3->Draw("same");
  higgsmass2->Draw("same");
  higgsmass1->Draw("same");
  leg->Draw();
  ((TH1D*)gDirectory->Get(name))->Draw("pesame");

  gPad->RedrawAxis();

  sprintf(newname,"%s%s%s","results_gg/data-mc_",allcut,"_2.gif");
  c0->SaveAs(newname);

  ((TH1D*)gDirectory->Get(name))->SetAxisRange(110.,149.);
  ((TH1D*)gDirectory->Get(name))->Draw("pe");

  sprintf(newname,"%s%s%s","results_gg/data_",allcut,"_red.gif");
  c0->SaveAs(newname);

  higgsmass6->Draw("same");
  higgsmass5->Draw("same");
  higgsmass4->Draw("same");
  higgsmass3->Draw("same");
  higgsmass2->Draw("same");
  higgsmass1->Draw("same");
  leg->Draw();
  ((TH1D*)gDirectory->Get(name))->Draw("pesame");

  gPad->RedrawAxis();

  sprintf(newname,"%s%s%s","results_gg/data-mc_",allcut,"_red.gif");
  c0->SaveAs(newname);

  ((TH1D*)gDirectory->Get(name))->Draw("pe");
  sprintf(name,"higgsmassdatacs");
  ((TH1D*)gDirectory->Get(name))->Draw("same");
  sprintf(name,"higgsmassdata");
  ((TH1D*)gDirectory->Get(name))->Draw("pesame");
  gPad->RedrawAxis();
  sprintf(newname,"%s%s%s","results_gg/datacs_",allcut,"_red.gif");
  leg2->Draw();
  c0->SaveAs(newname);

  ((TH1D*)gDirectory->Get(name))->SetAxisRange(150.,299.);
  ((TH1D*)gDirectory->Get(name))->Draw("pe");
  higgsmass6->Draw("same");
  higgsmass5->Draw("same");
  higgsmass4->Draw("same");
  higgsmass3->Draw("same");
  higgsmass2->Draw("same");
  higgsmass1->Draw("same");
  leg->Draw();
  ((TH1D*)gDirectory->Get(name))->Draw("pesame");

  gPad->RedrawAxis();


  sprintf(newname,"%s%s%s","results_gg/data-mc_",allcut,"_red-150-300.gif");
  c0->SaveAs(newname);


  ((TH1D*)gDirectory->Get(name))->Draw("pe");

  sprintf(newname,"%s%s%s","results_gg/data_",allcut,"_red-150-300.gif");
  c0->SaveAs(newname);


  ((TH1D*)gDirectory->Get(name))->SetAxisRange(299.,1000.);
  ((TH1D*)gDirectory->Get(name))->Draw("pe");
  higgsmass6->Draw("same");
  higgsmass5->Draw("same");
  higgsmass4->Draw("same");
  higgsmass3->Draw("same");
  higgsmass2->Draw("same");
  higgsmass1->Draw("same");
  leg->Draw();
  ((TH1D*)gDirectory->Get(name))->Draw("pesame");

  gPad->RedrawAxis();


  sprintf(newname,"%s%s%s","results_gg/data-mc_",allcut,"_red-3000-1000.gif");
  c0->SaveAs(newname);


  ((TH1D*)gDirectory->Get(name))->Draw("pe");

  sprintf(newname,"%s%s%s","results_gg/data_",allcut,"_red-3000-1000.gif");
  c0->SaveAs(newname);


  sprintf(newname,"%s%s%s","results_gg/yields_",allcut,".txt");
  ofstream outfile(newname);  
  outfile << "ndata   (100-150GeV) = " << num_data << endl;
  outfile << "nallbkg (100-150GeV) = " << num_box + num_diphotjet + num_gjet + num_qcd + num_qcd2 + num_dy << endl;
  outfile << "nvbf    (100-150GeV) = " << num_vbf << endl;
  outfile << "ndy     (100-150GeV) = " << num_dy << endl;
  outfile << "nbox    (100-150GeV) = " << num_box << endl;
  outfile << "ndiphot (100-150GeV) = " << num_diphotjet << endl;
  outfile << "ngjet   (100-150GeV) = " << num_gjet << endl;
  outfile << "nqcd    (100-150GeV) = " << num_qcd << endl;
  outfile << "nqcd30-40(100-150GeV)= " << num_qcd2 << endl;
  outfile << endl;
  outfile << "eff nvbf    (100-150GeV) = " << num_vbf_unscaled/n_vbf << endl;
  outfile << "eff ndy     (100-150GeV) = " << num_dy_unscaled/n_dy << endl;
  outfile << "eff nbox    (100-150GeV) = " << num_box_unscaled/n_box << endl;
  outfile << "eff ndiphot (100-150GeV) = " << num_diphotjet_unscaled/n_diphotjet << endl;
  outfile << "eff ngjet   (100-150GeV) = " << num_gjet_unscaled/n_gjet << endl;
  outfile << "eff nqcd    (100-150GeV) = " << num_qcd_unscaled/n_qcd << endl;
  outfile << "eff nqcd30-40(100-150GeV)= " << num_qcd2_unscaled/n_qcd2 << endl;
  outfile << endl;
  outfile << "nallbkg_2011 (100-150GeV) = " << num_box_2011 + num_diphotjet_2011 + num_gjet_2011 + num_qcd_2011 + num_qcd2_2011 + num_dy_2011 << endl;
  outfile << "nvbf_2011    (100-150GeV) = " << num_vbf_2011 << endl;
  outfile << "ndy_2011     (100-150GeV) = " << num_dy_2011 << endl;
  outfile << "nbox_2011    (100-150GeV) = " << num_box_2011 << endl;
  outfile << "ndiphot_2011 (100-150GeV) = " << num_diphotjet_2011 << endl;
  outfile << "ngjet_2011   (100-150GeV) = " << num_gjet_2011 << endl;
  outfile << "nqcd_2011    (100-150GeV) = " << num_qcd_2011 << endl;
  outfile << "nqcd30-40_2011(100-150GeV)= " << num_qcd2_2011 << endl;
  outfile << endl;
  outfile << "eff nvbf_2011    (100-150GeV) = " << num_vbf_2011_unscaled/n_vbf_2011 << endl;
  outfile << "eff ndy_2011     (100-150GeV) = " << num_dy_2011_unscaled/n_dy_2011 << endl;
  outfile << "eff nbox_2011    (100-150GeV) = " << num_box_2011_unscaled/n_box_2011 << endl;
  outfile << "eff ndiphot_2011 (100-150GeV) = " << num_diphotjet_2011_unscaled/n_diphotjet_2011 << endl;
  outfile << "eff ngjet_2011   (100-150GeV) = " << num_gjet_2011_unscaled/n_gjet_2011 << endl;
  outfile << "eff nqcd_2011    (100-150GeV) = " << num_qcd_2011_unscaled/n_qcd_2011 << endl;
  outfile << "eff nqcd30-40_2011(100-150GeV)= " << num_qcd2_2011_unscaled/n_qcd2_2011 << endl;
  outfile << endl;

  outfile.close();

 

  data.cd();
  sprintf(name,"higgsmassdata");
  ((TH1D*)gDirectory->Get(name))->SetAxisRange(90.,180.);
  ((TH1D*)gDirectory->Get(name))->Draw("pe");
  sprintf(name,"higgsmassdata2010");
  ((TH1D*)gDirectory->Get(name))->SetFillColor(kRed);
  ((TH1D*)gDirectory->Get(name))->Draw("same");
  gPad->RedrawAxis();

  sprintf(newname,"%s%s%s","results_gg/data_",allcut,"_2010-2011.gif");
  c0->SaveAs(newname);


  TF1 *gaussian_p;
  gaussian_p = new TF1("gausconst",gausconst,90.,180., 5) ;
  gaussian_p->SetLineColor(kBlue);
  gaussian_p->SetParNames ("Mean1","Sigma1","Norm1","p0","p1");
  gaussian_p->SetParameter(0, 118);
  gaussian_p->SetParLimits(0, 110, 125);
  gaussian_p->SetParameter(1, 1.5);
  gaussian_p->SetParLimits(1, .8,3.5);
  gaussian_p->SetParameter(2, 7.);
  gaussian_p->SetParLimits(2, .1,1000.);
  gaussian_p->SetParameter(3, 1.);
  gaussian_p->SetParLimits(3, 0.,1000.);  
  gaussian_p->SetParameter(4, 1.);
  gaussian_p->SetParLimits(4, -10.,10.);  

  higgsmass6->Fit ("gausconst","L","",90,180) ;
  gaussian_p->FixParameter(1, gaussian_p->GetParameter(1));

  sprintf(name,"higgsmassdata");
  ((TH1D*)gDirectory->Get(name))->SetAxisRange(90.,180.);
 
  ((TH1D*)gDirectory->Get(name))->Fit ("gausconst","L","",90,180) ;
  // ((TH1D*)gDirectory->Get(name))->Fit ("gausconst","L","",110,300) ;
  //  ((TH1D*)gDirectory->Get(name))->SetAxisRange(110.,150.);

  // higgsmass6->Draw("same");
  // higgsmass5->Draw("same");
  // higgsmass4->Draw("same");
  // higgsmass3->Draw("same");
  // higgsmass2->Draw("same");
  // higgsmass1->Draw("same");
  // ((TH1D*)gDirectory->Get(name))->Draw("pesame");
  sprintf(newname,"%s%s%s","results_gg/data-mc_",allcut,"_fit.gif");
  double mpi0= gaussian_p->GetParameter(0);
  double errmpi0= gaussian_p->GetParError(0);
  
  double smpi0= gaussian_p->GetParameter(1);
  double errsmpi0= gaussian_p->GetParError(1);
  
  double par2= gaussian_p->GetParameter(2);
  double events= gaussian_p->GetParameter(2)/((TH1D*)gDirectory->Get(name))->GetBinWidth(1);
  double errevents = gaussian_p->GetParError(2)/((TH1D*)gDirectory->Get(name))->GetBinWidth(1);

  TLatex tl;
  char line[100];
  double max = ((TH1D*)gDirectory->Get(name))->GetMaximum();

  tl.SetTextSize(.04);
  sprintf(line, "N=%3.1f#pm%3.1f", events, errevents); 
  tl.DrawLatex(145, max*0.7, line);
  sprintf(line, "#mu=(%3.1f#pm%3.1f)MeV", mpi0, errmpi0); 
  tl.DrawLatex(145, max*0.6, line);
  sprintf(line, "#sigma=(%3.1f#pm%3.1f)MeV", smpi0, errsmpi0); 
  tl.DrawLatex(145, max*0.5, line);
  
  leg->Draw("same");

  c0->SaveAs(newname);
  
}
