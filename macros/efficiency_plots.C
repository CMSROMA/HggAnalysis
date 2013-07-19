{
  gROOT->SetStyle("Plain");
   
  gStyle->SetPalette(1);
  gStyle->SetOptStat(1111111);  // Show overflow, underflow + SumOfWeights 
  gStyle->SetOptFit(111110); 
  gStyle->SetOptFile(1); 
  
  gStyle->SetMarkerStyle(20);
  gStyle->SetMarkerSize(.3);
  gStyle->SetMarkerColor(1);
  gStyle->SetMarkerColor(1);
  TCanvas* c0 = new TCanvas("c0"," ",200,10,500,500);
 
  TString redntpDir= "root://pccmsrm23.cern.ch:1094//u2/xrootd/delre/Higgs/reduced/";
  TString preselectionLevel;
  preselectionLevel="cicloose";
  TFile *glu = TFile::Open(redntpDir+"/redntp.42xv4."+preselectionLevel+".regrPho_eCorr.v3/merged/redntp_GluGluToHToGG_M-120_7TeV-powheg-pythia6-Fall11-PU_S6_START42_V14B-v1.root");
  TFile *vbf = TFile::Open(redntpDir+"/redntp.42xv4."+preselectionLevel+".regrPho_eCorr.v3/merged/redntp_VBF_HToGG_M-120_7TeV-powheg-pythia6-Fall11-PU_S6_START42_V14B-v1.root");  
  TFile *data = TFile::Open(redntpDir+"/redntp.42xv4_data_new."+preselectionLevel+".regrPho_eCorr.v3/merged/redntp_Photon-Run2011A-03Oct2011-05JulReReco-05AugReReco-Prompt-v1-DiPhotonSkimOnFly-b.root");

  glu->cd();
  TH1D npu_afterall("npu_afterall","",12,0,24);
  TH1D npu_afterphot("npu_afterphot","",12,0,24);
  AnaTree.Draw("npu>>npu_afterall","(ptphot1>55&&ptphot2>25&&ptcorrjet1>30&&ptcorrjet2>20&&abs(deltaeta)>3.5&&abs(zeppenjet)<2.5&&invmassjet>350&&abs(deltaphi)>2.6&&abs(etascphot1)<2.5&&abs(etascphot2)<2.5&&idcicphot1>=4&&idcicphot2>=4&&!((TMath::Abs(etascphot1)>1.4442&&TMath::Abs(etascphot1)<1.566)||(TMath::Abs(etascphot2)>1.4442&&TMath::Abs(etascphot2)<1.566))&&massgg>90&&massgg<190)")/ptphotgen1.GetEntries();
  AnaTree.Draw("npu>>npu_afterphot","(ptphot1>55&&ptphot2>25&&abs(etascphot1)<2.5&&abs(etascphot2)<2.5&&idcicphot1>=4&&idcicphot2>=4&&!((TMath::Abs(etascphot1)>1.4442&&TMath::Abs(etascphot1)<1.566)||(TMath::Abs(etascphot2)>1.4442&&TMath::Abs(etascphot2)<1.566))&&massgg>90&&massgg<190)")/ptphotgen1.GetEntries();
  npu_afterall.Sumw2();
  npu_afterall->SetStats(0);
  npu_afterall->SetMarkerStyle(8);
  npu_afterall->SetMarkerSize(.9);
  npu_afterall.SetMinimum(0);
  npu_afterall.SetMaximum(.04);
  npu_afterall.Divide(&npu_afterphot);
  npu_afterall.SetXTitle("npu");
  npu_afterall.SetYTitle("N(all cuts)/N(phot cuts)");
  npu_afterall.Fit("pol0");
  c0->SaveAs("efficiency_MC_glu.png");

  vbf->cd();
  TH1D npu_afterall("npu_afterall","",12,0,24);
  TH1D npu_afterphot("npu_afterphot","",12,0,24);
  AnaTree.Draw("npu>>npu_afterall","(ptphot1>55&&ptphot2>25&&ptcorrjet1>30&&ptcorrjet2>20&&abs(deltaeta)>3.5&&abs(zeppenjet)<2.5&&invmassjet>350&&abs(deltaphi)>2.6&&abs(etascphot1)<2.5&&abs(etascphot2)<2.5&&idcicphot1>=4&&idcicphot2>=4&&!((TMath::Abs(etascphot1)>1.4442&&TMath::Abs(etascphot1)<1.566)||(TMath::Abs(etascphot2)>1.4442&&TMath::Abs(etascphot2)<1.566))&&massgg>90&&massgg<190)")/ptphotgen1.GetEntries();
  AnaTree.Draw("npu>>npu_afterphot","(ptphot1>55&&ptphot2>25&&abs(etascphot1)<2.5&&abs(etascphot2)<2.5&&idcicphot1>=4&&idcicphot2>=4&&!((TMath::Abs(etascphot1)>1.4442&&TMath::Abs(etascphot1)<1.566)||(TMath::Abs(etascphot2)>1.4442&&TMath::Abs(etascphot2)<1.566))&&massgg>90&&massgg<190)")/ptphotgen1.GetEntries();
  npu_afterall.Sumw2();
  npu_afterall->SetStats(0);
  npu_afterall->SetMarkerStyle(8);
  npu_afterall->SetMarkerSize(.9);
  npu_afterall.SetMinimum(0);
  npu_afterall.SetMaximum(.5);
  npu_afterall.Divide(&npu_afterphot);
  npu_afterall.SetXTitle("npu");
  npu_afterall.SetYTitle("N(all cuts)/N(phot cuts)");
  npu_afterall.Fit("pol0");
  c0->SaveAs("efficiency_MC_vbf.png");

  data->cd();
  TH1D npu_afterall_data("npu_afterall_data","",12,0,24);
  TH1D npu_afterphot_data("npu_afterphot_data","",12,0,24);
  AnaTree.Draw("nvtx>>npu_afterall_data","(ptphot1>55&&ptphot2>25&&ptcorrjet1>30&&ptcorrjet2>20&&abs(deltaeta)>3.5&&abs(zeppenjet)<2.5&&invmassjet>350&&abs(deltaphi)>2.6&&abs(etascphot1)<2.5&&abs(etascphot2)<2.5&&idcicphot1>=4&&idcicphot2>=4&&!((TMath::Abs(etascphot1)>1.4442&&TMath::Abs(etascphot1)<1.566)||(TMath::Abs(etascphot2)>1.4442&&TMath::Abs(etascphot2)<1.566))&&massgg>90&&massgg<190)")/ptphotgen1.GetEntries();
  AnaTree.Draw("nvtx>>npu_afterphot_data","(ptphot1>55&&ptphot2>25&&abs(etascphot1)<2.5&&abs(etascphot2)<2.5&&idcicphot1>=4&&idcicphot2>=4&&!((TMath::Abs(etascphot1)>1.4442&&TMath::Abs(etascphot1)<1.566)||(TMath::Abs(etascphot2)>1.4442&&TMath::Abs(etascphot2)<1.566))&&massgg>90&&massgg<190)")/ptphotgen1.GetEntries();
  npu_afterall_data.Sumw2();
  npu_afterall_data.SetStats(0);
  npu_afterall_data.SetMarkerStyle(8);
  npu_afterall_data.SetMarkerSize(.2);
  npu_afterall_data.SetMinimum(0);
  npu_afterall_data.SetMaximum(.035);
  npu_afterall_data.Divide(&npu_afterphot_data);
  npu_afterall_data.SetXTitle("nvtx");
  npu_afterall_data.SetYTitle("N(all cuts)/N(phot cuts)");
  npu_afterall_data.Fit("pol0");
  c0->SaveAs("efficiency_data.png");



}
