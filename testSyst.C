{
  gROOT->ProcessLine(".L fillPlot.C++");
  gROOT->SetStyle("Plain");
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);

  TFile* file0=TFile::Open("redntp_GluGluToHToGG_M-120_7TeV-powheg-pythia6-Fall11-PU_S6_START42_V14B-v1_00_systcent.root");
  TFile* file1=TFile::Open("redntp_GluGluToHToGG_M-120_7TeV-powheg-pythia6-Fall11-PU_S6_START42_V14B-v1_00_systdown.root");
  TFile* file2=TFile::Open("redntp_GluGluToHToGG_M-120_7TeV-powheg-pythia6-Fall11-PU_S6_START42_V14B-v1_00_systup.root");

  fillPlot cent_fill((TTree*)file0->Get("AnaTree"), 1);
  fillPlot down_fill((TTree*)file1->Get("AnaTree"), 1);
  fillPlot up_fill((TTree*)file2->Get("AnaTree"), 1);
  float pt1=40;
  float pt2=30;
  float eb[3]={-1,0,1};
  float r9[3]={-1,0,1};

  cent_fill.setCic(4);
  down_fill.setCic(4);
  up_fill.setCic(4);
  cent_fill.DoPuReweight();
  down_fill.DoPuReweight();
  up_fill.DoPuReweight();
  
  for (int i=0;i<3;++i)
    {
      for (int j=0;j<3;++j)
	{
	  cent_fill.Setcuts(pt1,pt2,-100,-100,-10000,-10000,0,0,0,eb[i],r9[j],1);
	  down_fill.Setcuts(pt1,pt2,-100,-100,-10000,-10000,0,0,0,eb[i],r9[j],1);
	  up_fill.Setcuts(pt1,pt2,-100,-100,-10000,-10000,0,0,0,eb[i],r9[j],1);
	  
	  TLegend* leg=new TLegend(0.6,0.7,0.9,0.87);
	  leg->SetBorderSize(0);
	  leg->SetFillColor(0);
	  TH1D* cent = new TH1D("cent","cent",80,100,140);
	  cent->Sumw2();
	  cent->Add(cent_fill.Plot("massgg","masscent", 80,100,140));
	  cent->SetLineColor(1);
	  cent->SetFillColor(1);
	  cent->SetFillStyle(3001);
	  cent->SetLineWidth(2);
	  float mass_cent=cent->GetMean();
	  leg->AddEntry(cent,TString("default scale "),"FL");
	  cent->GetXaxis()->SetTitle("#gamma#gamma mass (GeV)");
	  cent->DrawNormalized("HIST");
	
	  TH1D* down = new TH1D("down","down",80,100,140);
	  down->Sumw2();
	  down->Add(down_fill.Plot("massgg","massdown", 80,100,140));
	  down->SetLineColor(2);
	  down->SetLineWidth(2);
	  float mass_down=down->GetMean();
	  leg->AddEntry(down,TString("-1 #sigma scale all categories"),"L");
	  down->DrawNormalized("HISTSAME");
	
	  TH1D* up = new TH1D("up","up",80,100,140);
	  up->Sumw2();
	  up->SetLineColor(3);
	  up->SetLineWidth(2);
	  up->GetMean();
	  up->Add(up_fill.Plot("massgg","massup", 80,100,140));
	  leg->AddEntry(up,TString("+1 #sigma scale all categories"),"L");
	  float mass_up=up->GetMean();
	  up->DrawNormalized("HISTSAME");
	
	
	  TLatex a;
	  char text[300];
	  sprintf(text,"Total syst: %.2f", (mass_up-mass_down)/2/mass_cent*100);
	  a.SetTextAlign(23);
	  a.SetTextSize(0.08);
	  a.DrawLatex(0.5,0.5,text);
	  leg->Draw();
	  c1->Update();
	  std::cout << "SYST ERROR ON MASS EB " << eb[i] << " R9 " << r9[j] <<  ":"  << (mass_up-mass_down)/2/mass_cent*100 << " %" << std::endl;
	  TString fileName("results_gg/energyScaleSyst_");
	  fileName+=eb[i];
	  fileName+="_";
	  fileName+=r9[j];
	  fileName+=".png";
	  c1->SaveAs(fileName);
	}
    }
}

  
