
// lumi runA 3.1pb-1 runB 30.9pb-1
//
// USE:
//
// .x var_tree_off("var name", "axis name", n bin in the histo, axis min, axis max, data int lumi, pt1 cut, pt2 cut, ptj1 cut, ptj2 cut, deltae cut, zep cut, mjj cut, 1st phot isol loose bool, 2nd phot isol loose bool, 1st phot isol medium bool, 2nd phot isol medium bool)
// 
// example:
//
// .x var_tree_off.C("ptjet1","p_{T}(jet1)[GeV]",9,20.,150.,34,20,20,20,15,2.5,2.5,300,1,1) 

void var_tree_off(char *var, char *axis, int nbin, double min, double max, double int_exp, double pt1=50, double pt2=30, double ptj1=20, double ptj2=15, double deltae=2.5, double zep=2.5, double mjj=300,double isol1=1, double isol2=1, double isom1=0, double isom2=0){

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

  char name[100];

  TFile data("redntp.V10/redntp_run2010.root");  
  TH1D* vardata = new TH1D("vardata","vardata",nbin,min,max);
  TH1D* vardatacs = new TH1D("vardatacs","vardatacs",nbin,min,max);

  TFile vbf("redntp.V10/redntp_VBF_HToGG_M-115_7TeV-powheg-pythia6_00.root");
  
  int n_vbf = ptphotgen1->GetEntries();
  cout << "# events vbf = " << n_vbf << endl;
  double cross_vbf = 1.3062 * 0.028;
  TH1D* varvbf = new TH1D("varvbf","varvbf",nbin,min,max);

  TFile box("redntp.V10/redntp_DiPhotonBox_Pt25to250_TrackingParticles_7TeV-pythia6_00.root");
  int n_box = ptphotgen1->GetEntries();
  cout << "# events box = " << n_box << endl;
  double cross_box = 12.37;
  TH1D* varbox = new TH1D("varbox","varbox",nbin,min,max);

  TFile diphotjet("redntp.V10/redntp_DiPhotonJets_7TeV-madgraph.root");
  int n_diphotjet = ptphotgen1->GetEntries();
  cout << "# events diphotjet = " << n_diphotjet << endl;
  double cross_diphotjet = 134;
  TH1D* vardiphotjet = new TH1D("vardiphojet","vardiphotjet",nbin,min,max);

  TFile gjet("redntp.V10/redntp_GJet_Pt-20_doubleEMEnriched_TuneZ2_7TeV-pythia6.root");
  int n_gjet = ptphotgen1->GetEntries();
  cout << "# events gjet = " << n_gjet << endl;
  double cross_gjet = 493.44;
  TH1D* vargjet = new TH1D("vargjet","vargjet",nbin,min,max);

  TFile qcd("redntp.V10/redntp_QCD_Pt-40_doubleEMEnriched_TuneZ2_7TeV-pythia6.root");
  int n_qcd = ptphotgen1->GetEntries();
  cout << "# events qcd = " << n_qcd << endl;
  double cross_qcd = 40392;
  TH1D* varqcd = new TH1D("varqcd","varqcd",nbin,min,max);

  TFile qcd2("redntp.V10/redntp_QCD_Pt-30to40_doubleEMEnriched_TuneZ2_7TeV-pythia6.root");
  int n_qcd2 = ptphotgen1->GetEntries();
  cout << "# events qcd2 = " << n_qcd2 << endl;
  double cross_qcd2 = 9610;
  TH1D* varqcd2 = new TH1D("varqcd2","varqcd2",nbin,min,max);


  char allcut[3000], allcut2[3000];
  char cut1[50],cut2[50],cut3[50],cut4[50],cut5[50],cut6[50],cut7[50],cut8[50],cut9[50],cut9[50],cut10[50],cut11[50], cutextra[200];
  sprintf(cut1,"%s%f","ptphot1>",pt1);
  sprintf(cut2,"%s%f","ptphot2>",pt2);
  sprintf(cut3,"%s%f","ptjet1>",ptj1);
  sprintf(cut4,"%s%f","ptjet2>",ptj2);
  sprintf(cut5,"%s%f","abs(deltaeta)>",deltae);
  sprintf(cut6,"%s%f","abs(zeppenjet)<",zep);
  // TEMP
  // to select events with low jet inv mass
  // sprintf(cut7,"%s%f","invmassjet<",mjj);
  // 
  sprintf(cut7,"%s%f","invmassjet>",mjj);
  if(isol1)  sprintf(cut8,"idloosephot1");
  else   sprintf(cut8,"massgg>0");
  if(isol2)  sprintf(cut9,"idloosephot2");
  else   sprintf(cut9,"massgg>0");
  if(isom1)  sprintf(cut10,"idmediumphot1");
  else   sprintf(cut10,"massgg>0");
  if(isom2)  sprintf(cut11,"idmediumphot2");
  else   sprintf(cut11,"massgg>0");
  if(isol1 && isol2) sprintf(cutextra,"((idloosephot1==0&&idloosephot2)||(idloosephot1&&idloosephot2==0))");
  if(isol1 && isom2) sprintf(cutextra,"((idloosephot1==0&&idmediumphot2)||(idloosephot1&&idmediumphot2==0))");
  if(isom1 && isol2) sprintf(cutextra,"((idmediumphot1==0&&idloosephot2)||(idmediumphot1&&idloosephot2==0))");
  if(isom1 && isom2) sprintf(cutextra,"((idmediumphot1==0&&idmediumphot2)||(idmediumphot1&&idmediumphot2==0))");
  sprintf(allcut,"%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s",cut1,"&&massgg>100&&massgg<150000&&abs(etaphot1)<2.5&&abs(etaphot2)<2.5&&",cut2,"&&",cut3,"&&",cut4,"&&",cut5,"&&",cut6,"&&",cut7,"&&",cut8,"&&",cut9,"&&",cut10,"&&",cut11);
  sprintf(allcut2,"%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s",cut1,"&&massgg>100&&massgg<150000&&abs(etaphot1)<2.5&&abs(etaphot2)<2.5&&",cut2,"&&",cut3,"&&",cut4,"&&",cut5,"&&",cut6,"&&",cut7,"&&",cutextra);

  // TEMP
  // EGAMMA isolation
  //
  // if(isol1)  sprintf(cut8,"idtightnewEGphot1&&!pid_haspixelseedphot1");
  // else   sprintf(cut8,"massgg>0");
  // if(isol2)  sprintf(cut9,"idtightnewEGphot2&&!pid_haspixelseedphot2");
  // else   sprintf(cut9,"massgg>0");
  // if(isom1)  sprintf(cut10,"idhggtightnewEGphot1&&!pid_haspixelseedphot1");
  // else   sprintf(cut10,"massgg>0");
  // if(isom2)  sprintf(cut11,"idhggtightnewEGphot2&&!pid_haspixelseedphot2");
  // else   sprintf(cut11,"massgg>0");
  // sprintf(allcut,"%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s",cut1,"&&massgg>123&&massgg<131&&!(abs(etaphot1)>1.4442&&abs(etaphot1)<1.566)&&!(abs(etaphot2)>1.4442&&abs(etaphot2)<1.566)&&abs(etaphot1)<2.5&&abs(etaphot2)<2.5&&",cut2,"&&",cut3,"&&",cut4,"&&",cut5,"&&",cut6,"&&",cut7,"&&",cut8,"&&",cut9,"&&",cut10,"&&",cut11);
  // sprintf(allcut2,"%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s",cut1,"&&massgg>123&&massgg<131&&!(abs(etaphot1)>1.4442&&abs(etaphot1)<1.566)&&!(abs(etaphot2)>1.4442&&abs(etaphot2)<1.566)&&abs(etaphot1)<2.5&&abs(etaphot2)<2.5&&",cut2,"&&",cut3,"&&",cut4,"&&",cut5,"&&",cut6,"&&",cut7,"&&",cutextra);
  //
  cout << allcut << endl;
  cout << allcut2 << endl;

  data.cd();
  double scale;
  sprintf(name,"vardata");
  AnaTree->Project(name,var,allcut);
  char tempcut[10000]; sprintf(tempcut,"%s%s",allcut,"&&massgg>123&&massgg<131");
  sprintf(tempcut,"%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s",cut1,"&&!(massgg>123&&massgg<131)&&massgg>90&&massgg<300&&abs(etaphot1)<2.5&&abs(etaphot2)<2.5&&",cut2,"&&",cut3,"&&",cut4,"&&",cut5,"&&",cut6,"&&",cut7,"&&",cut8,"&&",cut9,"&&",cut10,"&&",cut11);
  double entries_nosb = AnaTree->GetEntries(tempcut);  
  double num_data = ((TH1D*)gDirectory->Get(name))->Integral(101,150);
  ((TH1D*)gDirectory->Get(name))->Sumw2();
  sprintf(name,"vardatacs");
  AnaTree->Project(name,var,allcut2);
  sprintf(tempcut,"%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s",cut1,"&&!(massgg>123&&massgg<131)&&massgg>90&&massgg<300&&abs(etaphot1)<2.5&&abs(etaphot2)<2.5&&",cut2,"&&",cut3,"&&",cut4,"&&",cut5,"&&",cut6,"&&",cut7,"&&",cutextra);
  double entries_nosb_cs = AnaTree->GetEntries(tempcut);  

  vbf.cd();
  scale = cross_vbf * int_exp / n_vbf;
  sprintf(name,"varvbf");
  char allcut2[1000];
  AnaTree->Project(name,var,allcut);
  double num_vbf = ((TH1D*)gDirectory->Get(name))->Integral(101,150)*scale;
  double num_vbf_unscaled = ((TH1D*)gDirectory->Get(name))->Integral(101,150);
  ((TH1D*)gDirectory->Get(name))->Scale(scale);

  box.cd();
  scale = cross_box * int_exp / n_box;
  sprintf(name,"varbox");
  AnaTree->Project(name,var,allcut);
  double num_box = ((TH1D*)gDirectory->Get(name))->Integral(101,150)*scale;
  double num_box_unscaled = ((TH1D*)gDirectory->Get(name))->Integral(101,150);
  ((TH1D*)gDirectory->Get(name))->Scale(scale);
  ((TH1D*)gDirectory->Get(name))->SetLineColor(kRed);

  diphotjet.cd();
  scale = cross_diphotjet * int_exp / n_diphotjet;
  sprintf(name,"vardiphojet");
  AnaTree->Project(name,var,allcut);
  double num_diphotjet = ((TH1D*)gDirectory->Get(name))->Integral(101,150)*scale;
  double num_diphotjet_unscaled = ((TH1D*)gDirectory->Get(name))->Integral(101,150);
  ((TH1D*)gDirectory->Get(name))->Scale(scale);
  ((TH1D*)gDirectory->Get(name))->SetLineColor(kBlue);

  gjet.cd();
  scale = cross_gjet * int_exp / n_gjet;
  sprintf(name,"vargjet");
  AnaTree->Project(name,var,allcut);
  double num_gjet = ((TH1D*)gDirectory->Get(name))->Integral(101,150)*scale;
  double num_gjet_unscaled = ((TH1D*)gDirectory->Get(name))->Integral(101,150);
  ((TH1D*)gDirectory->Get(name))->Scale(scale);
  ((TH1D*)gDirectory->Get(name))->SetLineColor(kMagenta);

  qcd.cd();
  scale = cross_qcd * int_exp / n_qcd;
  sprintf(name,"varqcd");
  AnaTree->Project(name,var,allcut);
  double num_qcd = ((TH1D*)gDirectory->Get(name))->Integral(101,150)*scale;
  double num_qcd_unscaled = ((TH1D*)gDirectory->Get(name))->Integral(101,150);
  ((TH1D*)gDirectory->Get(name))->Scale(scale);
  ((TH1D*)gDirectory->Get(name))->SetLineColor(kGreen);

  qcd2.cd();
  scale = cross_qcd2 * int_exp / n_qcd2;
  sprintf(name,"varqcd2");
  AnaTree->Project(name,"massgg",allcut);
  double num_qcd2 = ((TH1D*)gDirectory->Get(name))->Integral(101,150)*scale;
  double num_qcd2_unscaled = ((TH1D*)gDirectory->Get(name))->Integral(101,150);
  ((TH1D*)gDirectory->Get(name))->Scale(scale);
  ((TH1D*)gDirectory->Get(name))->SetLineColor(kGreen);
  varqcd->Add(((TH1D*)gDirectory->Get(name)));


  vbf.cd();
  sprintf(name,"varvbf");
  TH1D* var1 = new TH1D("var1","var1",nbin,min,max);
  TH1D* var2 = new TH1D("var2","var2",nbin,min,max);
  TH1D* var3 = new TH1D("var3","var3",nbin,min,max);
  TH1D* var4 = new TH1D("var4","var4",nbin,min,max);
  TH1D* var5 = new TH1D("var5","var5",nbin,min,max);

  for (int j=0; j<5; j++){

     for (int i=1; i<((TH1D*)gDirectory->Get(name))->GetNbinsX()+1; i++){      
      
      int k = i;
      // if(j!=3) k=int((i-1)/10.+1);      
      
      if(j==0) {
	box.cd();
	sprintf(name,"varbox");
	var1->SetBinContent(i,((TH1D*)gDirectory->Get(name))->GetBinContent(k) + var1->GetBinContent(i));
	var2->SetBinContent(i,((TH1D*)gDirectory->Get(name))->GetBinContent(k) + var2->GetBinContent(i));
	var3->SetBinContent(i,((TH1D*)gDirectory->Get(name))->GetBinContent(k) + var3->GetBinContent(i));
	var4->SetBinContent(i,((TH1D*)gDirectory->Get(name))->GetBinContent(k) + var4->GetBinContent(i));
	var5->SetBinContent(i,((TH1D*)gDirectory->Get(name))->GetBinContent(k) + var4->GetBinContent(i));
      }
      if(j==1) {
	diphotjet.cd();
	sprintf(name,"vardiphojet");
	var2->SetBinContent(i,((TH1D*)gDirectory->Get(name))->GetBinContent(k) + var2->GetBinContent(i));
	var3->SetBinContent(i,((TH1D*)gDirectory->Get(name))->GetBinContent(k) + var3->GetBinContent(i));
	var4->SetBinContent(i,((TH1D*)gDirectory->Get(name))->GetBinContent(k) + var4->GetBinContent(i));
	var5->SetBinContent(i,((TH1D*)gDirectory->Get(name))->GetBinContent(k) + var4->GetBinContent(i));
      }
      if(j==2) {
	gjet.cd();
	sprintf(name,"vargjet");
	var3->SetBinContent(i,((TH1D*)gDirectory->Get(name))->GetBinContent(k) + var3->GetBinContent(i));
	var4->SetBinContent(i,((TH1D*)gDirectory->Get(name))->GetBinContent(k) + var4->GetBinContent(i));
	var5->SetBinContent(i,((TH1D*)gDirectory->Get(name))->GetBinContent(k) + var4->GetBinContent(i));
      }
      if(j==3) {
	qcd.cd();
	sprintf(name,"varqcd");
	var4->SetBinContent(i,((TH1D*)gDirectory->Get(name))->GetBinContent(k) + var4->GetBinContent(i));
	var5->SetBinContent(i,((TH1D*)gDirectory->Get(name))->GetBinContent(k) + var4->GetBinContent(i));
      }
      if(j==4) {
	vbf.cd();
	sprintf(name,"varvbf");
	var5->SetBinContent(i,((TH1D*)gDirectory->Get(name))->GetBinContent(k) + var4->GetBinContent(i));
      }
	
    }
    
    
  }
 
  var5->SetTitle("");
  var5->SetStats(0);
  var5->SetTitleOffset(1.25,"Y");
  var5->SetXTitle(axis);
  char ytitle[100];
  sprintf(ytitle,"%s%d%s","N_{ev}/",int_exp,"pb^{-1}");
  var5->SetYTitle(ytitle);
  var5->SetLineColor(kBlack);
  var5->SetLineWidth(2);
  var5->SetFillColor(kYellow);
  var5->SetMinimum(0);
  var5->Draw();

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
  leg = new TLegend(0.6,0.65,0.85,0.85);
  leg->SetFillStyle(0); leg->SetBorderSize(0); leg->SetTextSize(0.05);
  leg->SetFillColor(0);
  legge = leg->AddEntry(var5, "h_{f} M(130)", "f");
  legge = leg->AddEntry(var4, "QCD", "f");
  legge = leg->AddEntry(var3, "#gamma + jets", "f");
  legge = leg->AddEntry(var2, "di-#gamma + jets", "f");
  legge = leg->AddEntry(var1, "di-#gamma box", "f");
  leg->Draw();

  sprintf(allcut,"%f%s%f%s%f%s%f%s%f%s%f%s%f%s%d%s%d%s%d%s%d",pt1,"_",pt2,"_",ptj1,"_",ptj2,"_",deltae,"_",zep,"_",mjj,"_",isol1,"_",isol2,"_",isom1,"_",isom2);
  char newname[1000];
  sprintf(newname,"%s%s%s%s%s","results_gg/mc_",var,"_",allcut,".gif");
  cout << newname << endl;
  c0->SaveAs(newname);

  data.cd();
  sprintf(name,"vardata");
  ((TH1D*)gDirectory->Get(name))->SetXTitle(axis);
  ((TH1D*)gDirectory->Get(name))->SetTitle("");
  ((TH1D*)gDirectory->Get(name))->SetStats(0);
  ((TH1D*)gDirectory->Get(name))->SetMarkerStyle(8);
  ((TH1D*)gDirectory->Get(name))->SetMarkerSize(.9);
  ((TH1D*)gDirectory->Get(name))->SetTitleOffset(1.25,"Y");
  ((TH1D*)gDirectory->Get(name))->SetMinimum(0);
  ((TH1D*)gDirectory->Get(name))->Draw("pe");
  sprintf(newname,"%s%s%s%s%s","results_gg/data_",var,"_",allcut,".gif");
  c0->SaveAs(newname);

  double max =   ((TH1D*)gDirectory->Get(name))->GetMaximum();
  if(var5->GetMaximum()>max) max = var5->GetMaximum();
  ((TH1D*)gDirectory->Get(name))->SetMaximum(max*1.3);
  var5->Draw("same");
  var4->Draw("same");
  var3->Draw("same");
  var2->Draw("same");
  var1->Draw("same");
  leg->Draw();
  ((TH1D*)gDirectory->Get(name))->Draw("pesame");

  gPad->RedrawAxis();

  sprintf(newname,"%s%s%s%s%s","results_gg/data-mc_",var,"_",allcut,".gif");
  c0->SaveAs(newname);

  sprintf(name,"vardata");
  ((TH1D*)gDirectory->Get(name))->Draw("pe");
  double norm = entries_nosb;
  sprintf(name,"vardatacs");
  double normcs = entries_nosb_cs;
  cout << norm << "  " << normcs << endl;
  ((TH1D*)gDirectory->Get(name))->Scale(norm/normcs); 
  ((TH1D*)gDirectory->Get(name))->SetLineColor(46);
  ((TH1D*)gDirectory->Get(name))->SetFillColor(42);
  ((TH1D*)gDirectory->Get(name))->SetLineWidth(3);
  ((TH1D*)gDirectory->Get(name))->Draw("same");
  sprintf(name,"vardata");
  ((TH1D*)gDirectory->Get(name))->Draw("pesame");
  sprintf(newname,"%s%s%s%s%s","results_gg/datacs_",var,"_",allcut,".gif");
  gPad->RedrawAxis();
  TLegendEntry *legge2;
  TLegend *leg2;
  leg2 = new TLegend(0.6,0.65,0.9,0.85);
  leg2->SetFillStyle(0); leg2->SetBorderSize(0); leg2->SetTextSize(0.035);
  leg2->SetFillColor(0);
  legge2 = leg2->AddEntry(((TH1D*)gDirectory->Get("vardata")), "default sel.", "p");
  legge2 = leg2->AddEntry(((TH1D*)gDirectory->Get("vardatacs")), "control sample", "f");
  leg2->Draw();
  c0->SaveAs(newname);

}
