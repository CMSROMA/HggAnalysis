
Double_t fun(Double_t *x, Double_t *par)
{

  Double_t total = par[0] * ( 1 + par[1] * x[0] + par[2] * x[0] * x[0] );

  return total;
}

void fit(int p2){

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
  sprintf(name,"%i",p2);

  TString thename = name;

  TFile vbf("results_gg/histo_massgg_55.0_25.0_-10000.0_-10000.0_30.0_20.0_3.5_2.0_400.0_-1_-1_0_4.root");

  double scale = vardatacs->GetEntries()/vardatacs->Integral();

  vardatacs->Scale(scale);

  TF1 *myfun;
  myfun = new TF1("myfun",fun, 90.,190., 3) ;
  myfun->SetParameter(0,100);
  myfun->SetParameter(1,-0.001);
  myfun->SetParLimits(0,0.,10000000);
  myfun->SetParLimits(1,-10,10);
  myfun->SetParLimits(2,-100.,100);
  if (!p2) myfun->FixParameter(2,0.);

  myfun->SetLineColor(kBlue);

  vardatacs->SetMinimum(0);
  vardatacs->SetMarkerStyle(8);
  vardatacs->SetMarkerSize(.9);
  vardatacs->SetLineColor(1);
  //  vardatacs->SetLineSize(1);
  vardatacs->Fit ("myfun","L","",90.,190.) ;
  vardatacs->Draw("pe");

  c0->SaveAs("fitvbfcs"+thename+".png");

  myfun->FixParameter(0,myfun->GetParameter(0) / scale);
  
  double maxbin = vardata->GetMaximumBin();
  double max = vardata->GetBinContent(maxbin);
  vardata->SetMaximum( max * 1.8 );

  vardata->SetAxisRange(90,190);
  vardata->Draw("pe");
  myfun->Draw("same");
  
  c0->SaveAs("fitvbfsig"+thename+".png");

  myfun->SetParLimits(0,-10.,1000);
  myfun->SetParLimits(1,-10,10);

  vardata->SetAxisRange(90,190);
  vardata->Fit ("myfun","L","",90.,190.) ;
  
  c0->SaveAs("fitvbfsigfloat"+thename+".png");

  TFile hstra("results_gg/histo_massgg_65.0_25.0_-10000.0_-10000.0_35.0_20.0_-2.5_1.5_-30.0_-1_-1_0_4.root");

  scale = vardatacs->GetEntries()/vardatacs->Integral();
  vardatacs->Scale(scale);

  if(p2) myfun->SetParLimits(2,-100.,100);
  myfun->SetParLimits(0,0,100000);
 
  vardatacs->SetMinimum(0);
  vardatacs->SetMarkerStyle(8);
  vardatacs->SetMarkerSize(.9);
  vardatacs->SetLineColor(1);
  //  vardatacs->SetLineSize(1);
  vardatacs->Fit ("myfun","L","",90.,190.) ;
  vardatacs->Draw("pe");

  c0->SaveAs("fithstracs"+thename+".png");

  myfun->FixParameter(0,myfun->GetParameter(0) / scale);
  
  vardata->SetAxisRange(90,190);
  vardata->Draw("pe");
  myfun->Draw("same");
  
  c0->SaveAs("fithstrasig"+thename+".png");

 
  myfun->SetParLimits(0,-100.,1000);
  myfun->SetParLimits(1,-10,10);
  if(p2) myfun->SetParLimits(2,-100.,100);
  
  vardata->SetAxisRange(90,190);
  vardata->Fit ("myfun","L","",90.,190.) ;
  
  c0->SaveAs("fithstrasigfloat"+thename+".png");


}
