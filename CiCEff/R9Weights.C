{
  bool doVariableR9binning=true;

  gROOT->SetStyle("Plain");
  gROOT->ProcessLine(".L FindVariableBins.C");

  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);

  
  TFile *_file0 = TFile::Open("/cmsrm/pc24_2/meridian/data/redntp_DYJetsToLL_TuneZ2_M-50_7TeV-madgraph-tauola.root");
  TFile* f= TFile::Open("root://pccmsrm23.cern.ch:1094//u2/xrootd/meridian/Higgs/reduced/redntp.42xv3.cicloose.-1.v2/redntp_GluGluToHToGG_M-90_7TeV-powheg-pythia6_00.root");

  TString etacatName[7]={"All","EB","EE","EBlowEta","EBhighEta","EElowEta","EEhighEta"};
  float etaMin[7]={0.,0.,1.5,0.,1.,1.5,2.};
  float etaMax[7]={2.5,1.5,2.5,1.,1.5,,2.,2.5};
  
  TString r9Cat[3]={"All","Gold","Bad"};
  float r9CutMin[3]={-999.,0.94,-999.};
  float r9CutMax[3]={999.,999.,0.94};

  Int_t R9bins[3]={50,50,50};
  Float_t R9BinValues[3][500];

  for (int ieta=0;ieta<7;++ieta)
    {
      int etaRange=0;
      if (etacatName[ieta].Contains("EB"))
	etaRange=1;
      else if (etacatName[ieta].Contains("EE"))
	etaRange=2;
      for (int ir9=0;ir9<3;++ir9)
	{
	  TTree* a=(TTree*)_file0->Get("AnaTree");
	  TTree* b=(TTree*)f->Get("AnaTree");
	  TH1F* r9Ele;
	  if ((ieta==0 || ieta==1 || ieta==2) && ir9==0)
	      r9Ele=new TH1F("r9Ele"+etacatName[ieta]+r9Cat[ir9],"r9Ele"+etacatName[ieta]+r9Cat[ir9],2100,0.,1.05);
	  else 
	    {
	      if (doVariableR9binning)
		{
		  r9Ele=new TH1F("r9Ele"+etacatName[ieta]+r9Cat[ir9],"r9Ele"+etacatName[ieta]+r9Cat[ir9],R9bins[etaRange],&R9BinValues[etaRange][0]);
		}
	      else
	      r9Ele=new TH1F("r9Ele"+etacatName[ieta]+r9Cat[ir9],"r9Ele"+etacatName[ieta]+r9Cat[ir9],210,0.,1.05);
	    }
	  
	  TH1F* etaEle=new TH1F("etaEle"+etacatName[ieta]+r9Cat[ir9],"etaEle"+etacatName[ieta]+r9Cat[ir9],60,0.,3.);
	  TH1F* etaEleR9=new TH1F("etaEleR9"+etacatName[ieta]+r9Cat[ir9],"etaEleR9"+etacatName[ieta]+r9Cat[ir9],60,0.,3.);
	
	  TString accCut=" && !(abs(etascphot1)>1.4442 && abs(etascphot1)<1.566) &&";
	  accCut+=" abs(etascphot1)<";
	  accCut+=etaMax[ieta];
	  accCut+=" && abs(etascphot1)>";
	  accCut+=etaMin[ieta];
	  
	  TString r9Cut=" && r9phot1>";
	  r9Cut+=r9CutMin[ir9];
	  r9Cut+=" && r9phot1<";
	  r9Cut+=r9CutMax[ir9];

	  a->Draw("r9phot1>>r9Ele"+etacatName[ieta]+r9Cat[ir9],"ptphot1>30. && idcicnoelvetophot1>=4 && pid_hasMatchedPromptElephot1==1"+accCut+r9Cut); 
	  a->Draw("abs(etascphot1)>>etaEle"+etacatName[ieta]+r9Cat[ir9],"ptphot1>30 && idcicnoelvetophot1>=4 && pid_hasMatchedPromptElephot1==1"+accCut); 
	  a->Draw("abs(etascphot1)>>etaEleR9"+etacatName[ieta]+r9Cat[ir9],"ptphot1>30 && idcicnoelvetophot1>=4 && pid_hasMatchedPromptElephot1==1"+accCut+r9Cut); 
	  if ((ieta==0 || ieta==1 || ieta==2) && ir9==0 && doVariableR9binning)
	    FindVariableBins(r9Ele,R9bins[etaRange],&R9BinValues[etaRange][0]);


	  TH1F* r9Pho;
	  if ((ieta==0 || ieta==1 || ieta==2) && ir9==0)
	    r9Pho= new TH1F("r9Pho"+etacatName[ieta]+r9Cat[ir9],"r9Pho"+etacatName[ieta]+r9Cat[ir9],2100,0.,1.05);
	  else
	    {
	      if (doVariableR9binning)
		r9Pho= new TH1F("r9Pho"+etacatName[ieta]+r9Cat[ir9],"r9Pho"+etacatName[ieta]+r9Cat[ir9],R9bins[etaRange],&R9BinValues[etaRange][0]);
	      else
		r9Pho= new TH1F("r9Pho"+etacatName[ieta]+r9Cat[ir9],"r9Pho"+etacatName[ieta]+r9Cat[ir9],210,0.,1.05);
	    }

	  TH1F* etaPho= new TH1F("etaPho"+etacatName[ieta]+r9Cat[ir9],"etaPho"+etacatName[ieta]+r9Cat[ir9],60,0.,3.);
	  TH1F* etaPhoR9= new TH1F("etaPhoR9"+etacatName[ieta]+r9Cat[ir9],"etaPhoR9"+etacatName[ieta]+r9Cat[ir9],60,0.,3.);
		
	  b->Draw("r9phot1>>r9Pho"+etacatName[ieta]+r9Cat[ir9],"ptphot1>30. && idcicphot1>=4"+accCut+r9Cut);  
	  b->Draw("abs(etascphot1)>>etaPho"+etacatName[ieta]+r9Cat[ir9],"ptphot1>30 && idcicphot1>=4"+accCut);
	  b->Draw("abs(etascphot1)>>etaPhoR9"+etacatName[ieta]+r9Cat[ir9],"ptphot1>30 && idcicphot1>=4"+accCut+r9Cut);
	
	  TLegend l(0.65,0.72,0.89,0.89);
	  l.SetFillColor(0);
	  l.SetBorderSize(0);
	  l.SetHeader("#eta after #gammaId, R9, p_{T}>30");
	  //	  l.SetHeader(accCut+r9Cut);
	

	  TH1F* etaEleEff= (TH1F*) etaEleR9->Clone("etaEleEff");
	  etaEleEff->Divide(etaEle);
	  etaEleR9->Scale(1./etaEleR9->Integral());
	  etaEleR9->Draw();
	  etaEleR9->GetXaxis()->SetTitle("#eta");
	  etaEleR9->SetLineWidth(2);
	  l.AddEntry(etaEle,"Z#rightarrow ee","l");

	
	  TH1F* etaPhoEff= (TH1F*) etaPhoR9->Clone("etaPhoEff");
	  etaPhoEff->Divide(etaPho);
	  etaPhoR9->Scale(1./etaPhoR9->Integral());

	  etaPhoR9->SetLineColor(kBlue);
	  etaPhoR9->SetLineStyle(2);
	  etaPhoR9->SetLineWidth(2);
	  l.AddEntry(etaPhoR9,"H(90)#rightarrow #gamma#gamma","l");
	  etaPhoR9->Draw("SAME");

	  float maximum=TMath::Max(etaEleR9->GetMaximum(),etaPhoR9->GetMaximum())*1.7;
	  etaEle->SetMaximum(maximum);
	  gPad->Update();
	  l.Draw();
	  gPad->SaveAs("etaAfterR9"+etacatName[ieta]+r9Cat[ir9]+".png");
	
	
	  TLegend l(0.65,0.72,0.89,0.89);
	  l.SetFillColor(0);
	  l.SetBorderSize(0);

	
	  etaEleEff->Draw();
	  etaEleEff->GetXaxis()->SetTitle("#eta");
	  etaEleEff->GetYaxis()->SetTitle("fraction of e/#gamma");
	  etaEleEff->SetLineWidth(2);
	  l.AddEntry(etaEleEff,"Z#rightarrow ee","l");

	  etaPhoEff->SetLineColor(kBlue);
	  etaPhoEff->SetLineStyle(2);
	  etaPhoEff->SetLineWidth(2);
	  l.AddEntry(etaPhoEff,"H(90)#rightarrow #gamma#gamma","l");
	  etaPhoEff->Draw("SAME");
	
	
	  float maximum=TMath::Max(etaEleEff->GetMaximum(),etaPhoEff->GetMaximum())*1.3;
	  etaEleEff.SetMaximum(maximum);
  
	  gPad->Update();
	  l.Draw();
	  gPad->SaveAs("etaEff"+etacatName[ieta]+r9Cat[ir9]+".png");
	
	  TLegend l(0.65,0.72,0.89,0.89);
	  l.SetFillColor(0);
	  l.SetBorderSize(0);
	  //	  l.SetHeader(accCut+r9Cut);
  
	  etaPhoEff->Divide(etaEleEff);
	  etaPhoEff->GetXaxis()->SetTitle("#eta");
	  etaPhoEff->GetYaxis()->SetTitle("ratio of #gamma/e eff");
	  etaPhoEff->SetLineWidth(2);
	  etaPhoEff->SetLineColor(1);
	  etaPhoEff->SetLineStyle(1);
	  TH1F* etaEleReweighted=(TH1F*)etaEleR9->Clone(TString(etaEleR9->GetName())+"_reweighted");
	  etaEleReweighted->Multiply(etaPhoEff);
	  
	  std::cout <<  etaEleR9->Integral() << " " << etaEleReweighted->Integral() << std::endl;

	  etaPhoEff->Scale( etaEleR9->Integral()/etaEleReweighted->Integral() );
	  etaPhoEff->Draw();
	  l.AddEntry(etaPhoEff,"H(90)#rightarrow #gamma#gamma/Z#rightarrow ee","l");
	
	  float maximum=etaPhoEff->GetMaximum()*1.3;
	  etaPhoEff->SetMaximum(maximum);
  
	  gPad->Update();
	  l.Draw();
	  gPad->SaveAs("etaEffGammaOverEle"+etacatName[ieta]+r9Cat[ir9]+".png");
	  etaPhoEff->SetName("R9EtaWeight"+etacatName[ieta]+r9Cat[ir9]);
	  etaPhoEff->SetTitle("R9EtaWeight"+etacatName[ieta]+r9Cat[ir9]);
	  etaPhoEff->SaveAs("R9EtaWeight"+etacatName[ieta]+r9Cat[ir9]+".root");

	  TLegend l(0.65,0.72,0.89,0.89);
	  l.SetFillColor(0);
	  l.SetBorderSize(0);
	  l.SetHeader("#eta after #gammaId  p_{T}>30 ");
	
	  etaEle->Scale(1./etaEle->Integral());
	  cout << etaEle->Integral() << endl;
	  etaEle->Draw();
	  etaEle->GetXaxis()->SetTitle("#eta");
	  etaEle->SetLineWidth(2);
	  l.AddEntry(etaEle,"Z#rightarrow ee","l");

	
	  etaPho->Scale(1./etaPho->Integral());
	  cout << etaPho->Integral() << endl;
	  etaPho->SetLineColor(kBlue);
	  etaPho->SetLineStyle(2);
	  etaPho->SetLineWidth(2);
	  l.AddEntry(etaPho,"H(90)#rightarrow #gamma#gamma","l");
	  etaPho->Draw("SAME");
	
	  float maximum=TMath::Max(etaEle->GetMaximum(),etaPho->GetMaximum())*1.3;
	  etaEle->SetMaximum(maximum);
	  etaEle->SetMinimum(0.);
	  gPad->Update();
	  l.Draw();
	  gPad->SaveAs("eta"+etacatName[ieta]+r9Cat[ir9]+".png");
	
	
	  TLegend l(0.15,0.7,0.45,0.87);
	  l.SetFillColor(0);
	  l.SetBorderSize(0);
	  l.SetHeader("R_{9} after #gammaId  p_{T}>30 ");
	
	  r9Ele->Scale(1./r9Ele->Integral());
	  cout << r9Ele->Integral() << endl;
	  TH1F* r9EleDraw=(TH1F*)r9Ele->Clone("r9EleDraw");
	  if (doVariableR9binning)
	    for (int ibin=1;ibin<r9EleDraw->GetNbinsX()+2;++ibin)
	      r9EleDraw->SetBinContent(ibin,r9EleDraw->GetBinContent(ibin)/r9EleDraw->GetBinWidth(ibin));
	  r9EleDraw->Draw();
	  r9EleDraw->GetXaxis()->SetTitle("R_{9}");
	  r9EleDraw->GetXaxis()->SetRangeUser(0.55,1.03);
	  r9EleDraw->SetLineWidth(2);
	  l.AddEntry(r9EleDraw,"Z#rightarrow ee","l");
	
	  r9Pho->Scale(1./r9Pho->Integral());
	  cout << r9Pho->Integral() << endl;
	  TH1F* r9PhoDraw=(TH1F*)r9Pho->Clone("r9PhoDraw");
	  if (doVariableR9binning)
	    for (int ibin=1;ibin<r9PhoDraw->GetNbinsX()+2;++ibin)
	      r9PhoDraw->SetBinContent(ibin,r9PhoDraw->GetBinContent(ibin)/r9PhoDraw->GetBinWidth(ibin));
	  r9PhoDraw->SetLineColor(kBlue+2);
	  r9PhoDraw->SetLineStyle(2);
	  r9PhoDraw->SetLineWidth(2);
	  r9PhoDraw->SetFillColor(kBlue);
	  r9PhoDraw->SetFillStyle(3002);
	  l.AddEntry(r9PhoDraw,"H(90)#rightarrow #gamma#gamma","F");
	  r9PhoDraw->Draw("SAME");
	
	  float maximum=TMath::Max(r9EleDraw->GetMaximum(),r9PhoDraw->GetMaximum())*1.1;
	  r9EleDraw->SetMaximum(maximum);
	  r9EleDraw->SetMinimum(0.);
	  r9EleDraw->Draw("SAME");
	  gPad->Update();
	  l.Draw();
	  gPad->SaveAs("r9"+etacatName[ieta]+r9Cat[ir9]+".png");

	  r9Pho->Divide(r9Ele);
	  r9Ele->Multiply(r9Pho);
	  r9Pho->Scale(1./r9Ele->Integral());

	  TLegend l(0.15,0.72,0.45,0.89);
	  l.SetFillColor(0);
	  l.SetBorderSize(0);
	  r9Pho->SetLineColor(1);
	  r9Pho->SetLineStyle(1);
	  r9Pho->SetLineWidth(2);
	  r9Pho->SetFillStyle(0);
	  if (r9Cat[ir9] == "Gold")
	    r9Pho->GetXaxis()->SetRangeUser(0.85,1.03);
	  else
	    r9Pho->GetXaxis()->SetRangeUser(0.55,1.03);

	  r9Pho->Draw();
	  l.AddEntry(r9Pho,"H(90)#rightarrow #gamma#gamma/Z#rightarrow ee","l");
	
	  float maximum=r9Pho->GetMaximum()*1.3;
	  r9Pho->SetMaximum(maximum);
  
	  gPad->Update();
	  l.Draw();
	  gPad->SaveAs("r9EffGammaOverEle"+etacatName[ieta]+r9Cat[ir9]+".png");
	  r9Pho->SetName("R9Weight"+etacatName[ieta]+r9Cat[ir9]);
	  r9Pho->SetTitle("R9Weight"+etacatName[ieta]+r9Cat[ir9]);
	  r9Pho->SaveAs("R9Weight"+etacatName[ieta]+r9Cat[ir9]+".root");

	}
    }
}
