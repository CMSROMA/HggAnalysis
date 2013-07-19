void jecplots(){

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

  TFile g("/tmp/delre/vbf.root");
  TFile g1("/tmp/delre/vbf1.root");
  TFile g2("/tmp/delre/vbf2.root");
  TFile g3("/tmp/delre/vbf3.root");
  TFile g4("/tmp/delre/vbf4.root");
  TFile h("/tmp/delre/vh.root");
  TFile h1("/tmp/delre/vh1.root");
  TFile h2("/tmp/delre/vh2.root");
  TFile h3("/tmp/delre/vh3.root");
  TFile h4("/tmp/delre/vh4.root");

  g.cd();
  JECresovbf->SetTitle("");
  JECresovbf->SetStats(0);
  JECresovbf->SetXTitle("(p_{T}^{meas}-p_{T}^{gen}))/p_{T}^{gen}");
  JECresovbf->Draw();
  g1.cd();
  JECresovbf->SetLineColor(kRed);
  JECresovbf->Draw("same");
  g2.cd();
  JECresovbf->SetLineColor(kBlue);
  JECresovbf->Draw("same");
  
  c0->SaveAs("JECsyst_vbf.png");
  
  h.cd();
  JECresovh->SetTitle("");
  JECresovh->SetStats(0);
  JECresovh->SetXTitle("(p_{T}^{meas}-p_{T}^{gen}))/p_{T}^{gen}");
  JECresovh->Draw();
  h1.cd();
  JECresovh->SetLineColor(kRed);
  JECresovh->Draw("same");
  h2.cd();
  JECresovh->SetLineColor(kBlue);
  JECresovh->Draw("same");
  
  c0->SaveAs("JECsyst_vh.png");
   
  g4.cd();
  JECresovbf->SetLineColor(kRed);
  JECresovbf->SetTitle("");
  JECresovbf->SetStats(0);
  JECresovbf->SetXTitle("(p_{T}^{meas}-p_{T}^{gen}))/p_{T}^{gen}");
  JECresovbf->Draw();
  g.cd();
  JECresovbf->Draw("same");
  g3.cd();
  JECresovbf->SetLineColor(kBlue);
  JECresovbf->Draw("same");
  
  c0->SaveAs("JERsyst_vbf.png");
  
  h4.cd();
  JECresovh->SetTitle("");
  JECresovh->SetLineColor(kRed);
  JECresovh->SetStats(0);
  JECresovh->SetXTitle("(p_{T}^{meas}-p_{T}^{gen}))/p_{T}^{gen}");
  JECresovh->Draw();
  h.cd();
  JECresovh->Draw("same");
  h3.cd();
  JECresovh->SetLineColor(kBlue);
  JECresovh->Draw("same");
  
  c0->SaveAs("JERsyst_vh.png");
   

}
