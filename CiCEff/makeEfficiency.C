{
  gROOT->SetStyle("Plain");
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);
  gStyle->SetPadGridY(1);

  TString ciclevels[4]={"Loose","Medium","Tight","SuperTight"};
  //  TFile f("redntp.42xv1.preselection.v1/redntp_GluGluToHToGG_M-115_7TeV-powheg-pythia6_00.root");
  TFile f("../redntp_GluGluToHToGG_M-115_7TeV-powheg-pythia6_00.root");

  TH1F den_lead(  "den_lead", "den_lead",70,30.,100.);
  TH1F num1_lead("num1_lead","num1_lead",70,30.,100.);
  TH1F num2_lead("num2_lead","num2_lead",70,30.,100.);
  TH1F num3_lead("num3_lead","num3_lead",70,30.,100.);
  TH1F num4_lead("num4_lead","num4_lead",70,30.,100.);

  TH1F den_sublead(  "den_sublead", "den_sublead",70,30.,100.);
  TH1F num1_sublead("num1_sublead","num1_sublead",70,30.,100.);
  TH1F num2_sublead("num2_sublead","num2_sublead",70,30.,100.);
  TH1F num3_sublead("num3_sublead","num3_sublead",70,30.,100.);
  TH1F num4_sublead("num4_sublead","num4_sublead",70,30.,100.);

  AnaTree->Draw("ptphot1>>den_lead","abs(etascphot1)>1.566 && ptphot1>30. && abs(etascphot1)<2.5");
  AnaTree->Draw("ptphot1>>num1_lead","abs(etascphot1)>1.566 && ptphot1>30. && abs(etascphot1)<2.5 && idcicphot1>0");
  AnaTree->Draw("ptphot1>>num2_lead","abs(etascphot1)>1.566 && ptphot1>30. && abs(etascphot1)<2.5 && idcicphot1>1");
  AnaTree->Draw("ptphot1>>num3_lead","abs(etascphot1)>1.566 && ptphot1>30. && abs(etascphot1)<2.5 && idcicphot1>2");
  AnaTree->Draw("ptphot1>>num4_lead","abs(etascphot1)>1.566 && ptphot1>30. && abs(etascphot1)<2.5 && idcicphot1>3");

  AnaTree->Draw("ptphot2>>den_sublead","abs(etascphot2)>1.566 && ptphot2>30. && abs(etascphot2)<2.5");
  AnaTree->Draw("ptphot2>>num1_sublead","abs(etascphot2)>1.566 && ptphot2>30. && abs(etascphot2)<2.5 && idcicphot2>0");
  AnaTree->Draw("ptphot2>>num2_sublead","abs(etascphot2)>1.566 && ptphot2>30. && abs(etascphot2)<2.5 && idcicphot2>1");
  AnaTree->Draw("ptphot2>>num3_sublead","abs(etascphot2)>1.566 && ptphot2>30. && abs(etascphot2)<2.5 && idcicphot2>2");
  AnaTree->Draw("ptphot2>>num4_sublead","abs(etascphot2)>1.566 && ptphot2>30. && abs(etascphot2)<2.5 && idcicphot2>3");

  TH1F* den= (TH1F*) den_lead.Clone("den");
  TH1F* num1=(TH1F*) num1_lead.Clone("num1");
  TH1F* num2=(TH1F*) num2_lead.Clone("num2");
  TH1F* num3=(TH1F*) num3_lead.Clone("num3");
  TH1F* num4=(TH1F*) num4_lead.Clone("num4");

  num1_lead.Divide(&den_lead);
  num2_lead.Divide(&den_lead);
  num3_lead.Divide(&den_lead);
  num4_lead.Divide(&den_lead);
 
  TH2F a("a","a",10,30.,100.,10,0.4,1.);
  a.Draw();
  a.GetXaxis()->SetTitle("p_{T} [GeV]");

  num1_lead.SetLineWidth(1.3);
  num1_lead.Draw("SAME");
  num2_lead.SetLineWidth(1.3);
  num2_lead.SetLineColor(2);
  num2_lead.Draw("SAME");
  num3_lead.SetLineWidth(1.3);
  num3_lead.SetLineColor(3);
  num3_lead.Draw("SAME");
  num4_lead.SetLineWidth(1.3);
  num4_lead.SetLineColor(4);
  num4_lead.Draw("SAME");

  leg=TLegend(0.8,0.1,0.99,0.3);
  leg.SetFillColor(0);
  leg.SetBorderSize(0);
  leg.AddEntry(&num1_lead,ciclevels[0],"l");
  leg.AddEntry(&num2_lead,ciclevels[1],"l");
  leg.AddEntry(&num3_lead,ciclevels[2],"l");
  leg.AddEntry(&num4_lead,ciclevels[3],"l");
  leg.Draw();
  c1->SaveAs("lead_EE_CiC_eff.png");

  den->Add(&den_sublead);
  num1->Add(&num1_sublead);
  num2->Add(&num2_sublead);
  num3->Add(&num3_sublead);
  num4->Add(&num4_sublead);

  num1_sublead.Divide(&den_sublead);
  num2_sublead.Divide(&den_sublead);
  num3_sublead.Divide(&den_sublead);
  num4_sublead.Divide(&den_sublead);
 
  a.Draw();
  num1_sublead.SetLineWidth(1.3);
  num1_sublead.Draw("SAME");
  num2_sublead.SetLineWidth(1.3);
  num2_sublead.SetLineColor(2);
  num2_sublead.Draw("SAME");
  num3_sublead.SetLineWidth(1.3);
  num3_sublead.SetLineColor(3);
  num3_sublead.Draw("SAME");
  num4_sublead.SetLineWidth(1.3);
  num4_sublead.SetLineColor(4);
  num4_sublead.Draw("SAME");

  leg.Draw();
  c1->SaveAs("sublead_EE_CiC_eff.png");

  num1->Divide(den);
  num2->Divide(den);
  num3->Divide(den);
  num4->Divide(den);

  a.Draw();
  num1->SetLineWidth(1.3);
  num1->Draw("SAME");
  num2->SetLineWidth(1.3);
  num2->SetLineColor(2);
  num2->Draw("SAME");
  num3->SetLineWidth(1.3);
  num3->SetLineColor(3);
  num3->Draw("SAME");
  num4->SetLineWidth(1.3);
  num4->SetLineColor(4);
  num4->Draw("SAME");

  leg.Draw();
  c1->SaveAs("EE_CiC_eff.png");

  TH1F den_lead(  "den_lead", "den_lead",70,30.,100.);
  TH1F num1_lead("num1_lead","num1_lead",70,30.,100.);
  TH1F num2_lead("num2_lead","num2_lead",70,30.,100.);
  TH1F num3_lead("num3_lead","num3_lead",70,30.,100.);
  TH1F num4_lead("num4_lead","num4_lead",70,30.,100.);

  TH1F den_sublead(  "den_sublead", "den_sublead",70,30.,100.);
  TH1F num1_sublead("num1_sublead","num1_sublead",70,30.,100.);
  TH1F num2_sublead("num2_sublead","num2_sublead",70,30.,100.);
  TH1F num3_sublead("num3_sublead","num3_sublead",70,30.,100.);
  TH1F num4_sublead("num4_sublead","num4_sublead",70,30.,100.);


  AnaTree->Draw("ptphot1>>den_lead", "abs(etascphot1)<1.4442 && ptphot1>30. ");
  AnaTree->Draw("ptphot1>>num1_lead","abs(etascphot1)<1.4442 && ptphot1>30.  && idcicphot1>0");
  AnaTree->Draw("ptphot1>>num2_lead","abs(etascphot1)<1.4442 && ptphot1>30.  && idcicphot1>1");
  AnaTree->Draw("ptphot1>>num3_lead","abs(etascphot1)<1.4442 && ptphot1>30.  && idcicphot1>2");
  AnaTree->Draw("ptphot1>>num4_lead","abs(etascphot1)<1.4442 && ptphot1>30.  && idcicphot1>3");

  AnaTree->Draw("ptphot2>>den_sublead", "abs(etascphot2)<1.4442 && ptphot2>30. ");
  AnaTree->Draw("ptphot2>>num1_sublead","abs(etascphot2)<1.4442 && ptphot2>30.  && idcicphot2>0");
  AnaTree->Draw("ptphot2>>num2_sublead","abs(etascphot2)<1.4442 && ptphot2>30.  && idcicphot2>1");
  AnaTree->Draw("ptphot2>>num3_sublead","abs(etascphot2)<1.4442 && ptphot2>30.  && idcicphot2>2");
  AnaTree->Draw("ptphot2>>num4_sublead","abs(etascphot2)<1.4442 && ptphot2>30.  && idcicphot2>3");

  TH1F* den= (TH1F*) den_lead.Clone("den");
  TH1F* num1=(TH1F*) num1_lead.Clone("num1");
  TH1F* num2=(TH1F*) num2_lead.Clone("num2");
  TH1F* num3=(TH1F*) num3_lead.Clone("num3");
  TH1F* num4=(TH1F*) num4_lead.Clone("num4");

  num1_lead.Divide(&den_lead);
  num2_lead.Divide(&den_lead);
  num3_lead.Divide(&den_lead);
  num4_lead.Divide(&den_lead);
 
  //  TH2F a("a","a",10,30.,100.,10,0.,1.);
  a.Draw();
  num1_lead.SetLineWidth(1.3);
  num1_lead.Draw("SAME");
  num2_lead.SetLineWidth(1.3);
  num2_lead.SetLineColor(2);
  num2_lead.Draw("SAME");
  num3_lead.SetLineWidth(1.3);
  num3_lead.SetLineColor(3);
  num3_lead.Draw("SAME");
  num4_lead.SetLineWidth(1.3);
  num4_lead.SetLineColor(4);
  num4_lead.Draw("SAME");

  leg.Draw();
  c1->SaveAs("lead_EB_CiC_eff.png");

  den->Add(&den_sublead);
  num1->Add(&num1_sublead);
  num2->Add(&num2_sublead);
  num3->Add(&num3_sublead);
  num4->Add(&num4_sublead);

  num1_sublead.Divide(&den_sublead);
  num2_sublead.Divide(&den_sublead);
  num3_sublead.Divide(&den_sublead);
  num4_sublead.Divide(&den_sublead);
 
  //  TH2F a("a","a",10,30.,100.,10,0.,1.);
  a.Draw();

  num1_sublead.SetLineWidth(1.3);
  num1_sublead.Draw("SAME");
  num2_sublead.SetLineWidth(1.3);
  num2_sublead.SetLineColor(2);
  num2_sublead.Draw("SAME");
  num3_sublead.SetLineWidth(1.3);
  num3_sublead.SetLineColor(3);
  num3_sublead.Draw("SAME");
  num4_sublead.SetLineWidth(1.3);
  num4_sublead.SetLineColor(4);
  num4_sublead.Draw("SAME");

  leg.Draw();
  c1->SaveAs("sublead_EB_CiC_eff.png");

  num1->Divide(den);
  num2->Divide(den);
  num3->Divide(den);
  num4->Divide(den);

  a.Draw();
  num1->SetLineWidth(1.3);
  num1->Draw("SAME");
  num2->SetLineWidth(1.3);
  num2->SetLineColor(2);
  num2->Draw("SAME");
  num3->SetLineWidth(1.3);
  num3->SetLineColor(3);
  num3->Draw("SAME");
  num4->SetLineWidth(1.3);
  num4->SetLineColor(4);
  num4->Draw("SAME");

  leg.Draw();
  c1->SaveAs("EB_CiC_eff.png");

  TH1F  den_leadRhoPF( "den_leadRhoPF", "den_leadRhoPF",30,0.,10.);
  TH1F num1_leadRhoPF("num1_leadRhoPF","num1_leadRhoPF",30,0.,10.);
  TH1F num2_leadRhoPF("num2_leadRhoPF","num2_leadRhoPF",30,0.,10.);
  TH1F num3_leadRhoPF("num3_leadRhoPF","num3_leadRhoPF",30,0.,10.);
  TH1F num4_leadRhoPF("num4_leadRhoPF","num4_leadRhoPF",30,0.,10.);

  AnaTree->Draw("rhoPF>>den_leadRhoPF","abs(etascphot1)>1.566 && ptphot1>30. && abs(etascphot1)<2.5");
  AnaTree->Draw("rhoPF>>num1_leadRhoPF","abs(etascphot1)>1.566 && ptphot1>30. && abs(etascphot1)<2.5 && idcicphot1>0");
  AnaTree->Draw("rhoPF>>num2_leadRhoPF","abs(etascphot1)>1.566 && ptphot1>30. && abs(etascphot1)<2.5 && idcicphot1>1");
  AnaTree->Draw("rhoPF>>num3_leadRhoPF","abs(etascphot1)>1.566 && ptphot1>30. && abs(etascphot1)<2.5 && idcicphot1>2");
  AnaTree->Draw("rhoPF>>num4_leadRhoPF","abs(etascphot1)>1.566 && ptphot1>30. && abs(etascphot1)<2.5 && idcicphot1>3");

  TH1F* den= (TH1F*) den_leadRhoPF.Clone("den");
  TH1F* num1=(TH1F*) num1_leadRhoPF.Clone("num1");
  TH1F* num2=(TH1F*) num2_leadRhoPF.Clone("num2");
  TH1F* num3=(TH1F*) num3_leadRhoPF.Clone("num3");
  TH1F* num4=(TH1F*) num4_leadRhoPF.Clone("num4");

  num1_leadRhoPF.Divide(&den_leadRhoPF);
  num2_leadRhoPF.Divide(&den_leadRhoPF);
  num3_leadRhoPF.Divide(&den_leadRhoPF);
  num4_leadRhoPF.Divide(&den_leadRhoPF);

  TH2F aRho("aRho","aRho",10,0.,10.,10,0.4,1.);
  aRho.Draw();
  aRho.GetXaxis()->SetTitle("#rho [GeV]"); 

  num1_leadRhoPF.SetLineWidth(1.3);
  num1_leadRhoPF.Draw("SAME");
  num2_leadRhoPF.SetLineWidth(1.3);
  num2_leadRhoPF.SetLineColor(2);
  num2_leadRhoPF.Draw("SAME");
  num3_leadRhoPF.SetLineWidth(1.3);
  num3_leadRhoPF.SetLineColor(3);
  num3_leadRhoPF.Draw("SAME");
  num4_leadRhoPF.SetLineWidth(1.3);
  num4_leadRhoPF.SetLineColor(4);
  num4_leadRhoPF.Draw("SAME");

  TLegend leg1=TLegend(0.8,0.1,0.99,0.3);
  leg1.SetFillColor(0);
  leg1.SetBorderSize(0);
  leg1.AddEntry(&num1_leadRhoPF,ciclevels[0],"l");
  leg1.AddEntry(&num2_leadRhoPF,ciclevels[1],"l");
  leg1.AddEntry(&num3_leadRhoPF,ciclevels[2],"l");
  leg1.AddEntry(&num4_leadRhoPF,ciclevels[3],"l");
  leg1.Draw();

  c1->SaveAs("leadRhoPF_EE_CiC_eff.png");


  TH1F  den_leadRhoPF( "den_leadRhoPF", "den_leadRhoPF",30,0.,10.);
  TH1F num1_leadRhoPF("num1_leadRhoPF","num1_leadRhoPF",30,0.,10.);
  TH1F num2_leadRhoPF("num2_leadRhoPF","num2_leadRhoPF",30,0.,10.);
  TH1F num3_leadRhoPF("num3_leadRhoPF","num3_leadRhoPF",30,0.,10.);
  TH1F num4_leadRhoPF("num4_leadRhoPF","num4_leadRhoPF",30,0.,10.);

  AnaTree->Draw( "rhoPF>>den_leadRhoPF","ptphot1>30. && abs(etascphot1)<1.4442");
  AnaTree->Draw("rhoPF>>num1_leadRhoPF","ptphot1>30. && abs(etascphot1)<1.4442 && idcicphot1>0");
  AnaTree->Draw("rhoPF>>num2_leadRhoPF","ptphot1>30. && abs(etascphot1)<1.4442 && idcicphot1>1");
  AnaTree->Draw("rhoPF>>num3_leadRhoPF","ptphot1>30. && abs(etascphot1)<1.4442 && idcicphot1>2");
  AnaTree->Draw("rhoPF>>num4_leadRhoPF","ptphot1>30. && abs(etascphot1)<1.4442 && idcicphot1>3");

  TH1F* den= (TH1F*) den_leadRhoPF.Clone("den");
  TH1F* num1=(TH1F*) num1_leadRhoPF.Clone("num1");
  TH1F* num2=(TH1F*) num2_leadRhoPF.Clone("num2");
  TH1F* num3=(TH1F*) num3_leadRhoPF.Clone("num3");
  TH1F* num4=(TH1F*) num4_leadRhoPF.Clone("num4");

  num1_leadRhoPF.Divide(&den_leadRhoPF);
  num2_leadRhoPF.Divide(&den_leadRhoPF);
  num3_leadRhoPF.Divide(&den_leadRhoPF);
  num4_leadRhoPF.Divide(&den_leadRhoPF);
 
  aRho.Draw();
  aRho.GetXaxis()->SetTitle("#rho [GeV]"); 

  num1_leadRhoPF.SetLineWidth(1.3);
  num1_leadRhoPF.Draw("SAME");
  num2_leadRhoPF.SetLineWidth(1.3);
  num2_leadRhoPF.SetLineColor(2);
  num2_leadRhoPF.Draw("SAME");
  num3_leadRhoPF.SetLineWidth(1.3);
  num3_leadRhoPF.SetLineColor(3);
  num3_leadRhoPF.Draw("SAME");
  num4_leadRhoPF.SetLineWidth(1.3);
  num4_leadRhoPF.SetLineColor(4);
  num4_leadRhoPF.Draw("SAME");

  leg1.Draw();
  c1->SaveAs("leadRhoPF_EB_CiC_eff.png");
}
