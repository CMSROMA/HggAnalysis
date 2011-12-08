void jecplots_glu(){

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

  TFile g("/tmp/delre/glu_official.root");
  TFile g1("/tmp/delre/glu_jec_1.root");
  TFile g2("/tmp/delre/glu_jec_2.root");
  TFile g3("/tmp/delre/glu_jec_3.root");
  TFile g4("/tmp/delre/glu_jec_4.root");
 
  g.cd();
  JECresovbf->Rebin(4);
  JECresovbf->SetTitle("");
  JECresovbf->SetStats(0);
  JECresovbf->SetXTitle("(p_{T}^{meas}-p_{T}^{gen}))/p_{T}^{gen}");
  JECresovbf->Draw();
  g1.cd();
  JECresovbf->Rebin(4);
  JECresovbf->SetLineColor(kRed);
  JECresovbf->Draw("same");
  g2.cd();
  JECresovbf->Rebin(4);
  JECresovbf->SetLineColor(kBlue);
  JECresovbf->Draw("same");
  
  c0->SaveAs("JECsyst_vbf.png");
  
 //  h.cd();
//   JECresovh->SetTitle("");
//   JECresovh->SetStats(0);
//   JECresovh->SetXTitle("(p_{T}^{meas}-p_{T}^{gen}))/p_{T}^{gen}");
//   JECresovh->Draw();
//   h1.cd();
//   JECresovh->SetLineColor(kRed);
//   JECresovh->Draw("same");
//   h2.cd();
//   JECresovh->SetLineColor(kBlue);
//   JECresovh->Draw("same");
  
//   c0->SaveAs("JECsyst_vh.png");
   
  g4.cd();
  JECresovbf->Rebin(4);
  JECresovbf->SetLineColor(kRed);
  JECresovbf->SetTitle("");
  JECresovbf->SetStats(0);
  JECresovbf->SetXTitle("(p_{T}^{meas}-p_{T}^{gen}))/p_{T}^{gen}");
  JECresovbf->Draw();
  g.cd();
  JECresovbf->Draw("same");
  g3.cd();
  JECresovbf->Rebin(4);
  JECresovbf->SetLineColor(kBlue);
  JECresovbf->Draw("same");
  
  c0->SaveAs("JERsyst_vbf.png");
  
//   h4.cd();
//   JECresovh->SetTitle("");
//   JECresovh->SetLineColor(kRed);
//   JECresovh->SetStats(0);
//   JECresovh->SetXTitle("(p_{T}^{meas}-p_{T}^{gen}))/p_{T}^{gen}");
//   JECresovh->Draw();
//   h.cd();
//   JECresovh->Draw("same");
//   h3.cd();
//   JECresovh->SetLineColor(kBlue);
//   JECresovh->Draw("same");
  
//   c0->SaveAs("JERsyst_vh.png");
   

}
