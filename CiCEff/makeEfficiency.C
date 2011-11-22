{
  gROOT->SetStyle("Plain");
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);
  gStyle->SetPadGridY(1);

  TString ciclevels[4]={"Loose","Medium","Tight","SuperTight"};
  //  TFile f("redntp.42xv1.preselection.v1/redntp_GluGluToHToGG_M-115_7TeV-powheg-pythia6_00.root");
  //TFile f("/xrootdfs/u2/xrootd/meridian/Higgs/reduced/redntp.42xv4_b.mcass.regrPho_eCorr.v1/redntp_GluGluToHToGG_M-120_7TeV-powheg-pythia6-Fall11-PU_S6_START42_V14B-v1_00.root");
  TFile f("/xrootdfs/u2/xrootd/meridian/Higgs/reduced/redntp.42xv4.mcass.regrPho_eCorr.v1/redntp_GluGluToHToGG_M-120_7TeV-powheg-pythia6-Fall11-PU_S6_START42_V14B-v1_00.root");
  ///xrootdfs/u2/xrootd/meridian/Higgs/reduced/redntp.42xv4.nosel.regrPho_eCorr.v1/merged/redntp_GluGluToHToGG_M-120_7TeV-powheg-pythia6-Fall11-PU_S6_START42_V14B-v1.root");
  //  TFile f("/xrootdfs/u2/xrootd/meridian/Higgs/reduced/redntp.42xv3.cicloose.-1.v2/merged/redntp_GluGluToHToGG_M-120_7TeV-powheg-pythia6.root");
  //  TFile f("/xrootdfs/u2/xrootd/meridian/Higgs/reduced/redntp.42xv3.preselectionCS.-1.v2/merged/redntp_GluGluToHToGG_M-120_7TeV-powheg-pythia6.root");


  TString category[9]= {
    "All",
    "All_highR9",
    "All_lowR9",
    "EB",
    "EE",
    "EB_highR9",
    "EB_lowR9",
    "EE_highR9",
    "EE_lowR9"
  }
  
  TString cut_category_1[9]= {
    "ptphot1>30.",
    "ptphot1>30. && r9phot1>=0.94",
    "ptphot1>30. && r9phot1<0.94",
    "abs(etascphot1)<1.4442 && ptphot1>30.",
    "abs(etascphot1)>1.566 && ptphot1>30. && abs(etascphot1)<2.5",
    "abs(etascphot1)<1.4442 && ptphot1>30. && r9phot1>=0.94",
    "abs(etascphot1)<1.4442 && ptphot1>30. && r9phot1<0.94",
    "abs(etascphot1)>1.566 && ptphot1>30. && abs(etascphot1)<2.5 && r9phot1>=0.94",
    "abs(etascphot1)>1.566 && ptphot1>30. && abs(etascphot1)<2.5 && r9phot1<0.94"
  }

  TString cut_category_2[9]= {
    "ptphot2>30.",
    "ptphot2>30. && r9phot2>=0.94",
    "ptphot2>30. && r9phot2<0.94",
    "abs(etascphot2)<1.4442 && ptphot2>30.",
    "abs(etascphot2)>1.566 && ptphot2>30. && abs(etascphot2)<2.5",
    "abs(etascphot2)<1.4442 && ptphot2>30. && r9phot2>=0.94",
    "abs(etascphot2)<1.4442 && ptphot2>30. && r9phot2<0.94",
    "abs(etascphot2)>1.566 && ptphot2>30. && abs(etascphot2)<2.5 && r9phot2>=0.94",
    "abs(etascphot2)>1.566 && ptphot2>30. && abs(etascphot2)<2.5 && r9phot2<0.94"
  }

  TH2F a("a","a",10,30.,100.,10,0.,1.);
  a.GetXaxis()->SetTitle("p_{T} [GeV]");  

  TH2F b("b","b",10,-3.,3.,10,0.,1.);
  b.GetXaxis()->SetTitle("#eta");

  TH2F aRho("aRho","aRho",10,0.,30.,10,0.,1.);
  aRho.GetXaxis()->SetTitle("#rho [GeV]");
 
  std::cout << "INIT DONE" << std::endl;
  for (int i=0;i<9;++i)
    {
      TH1F den_lead(  "den_lead_"+category[i], "den_lead_"+category[i],35,30.,100.);
      TH1F num1_lead("num1_lead_"+category[i],"num1_lead_"+category[i],35,30.,100.);
      TH1F num2_lead("num2_lead_"+category[i],"num2_lead_"+category[i],35,30.,100.);
      TH1F num3_lead("num3_lead_"+category[i],"num3_lead_"+category[i],35,30.,100.);
      TH1F num4_lead("num4_lead_"+category[i],"num4_lead_"+category[i],35,30.,100.);

      TH1F den_sublead(  "den_sublead_"+category[i], "den_sublead_"+category[i],35,30.,100.);
      TH1F num1_sublead("num1_sublead_"+category[i],"num1_sublead_"+category[i],35,30.,100.);
      TH1F num2_sublead("num2_sublead_"+category[i],"num2_sublead_"+category[i],35,30.,100.);
      TH1F num3_sublead("num3_sublead_"+category[i],"num3_sublead_"+category[i],35,30.,100.);
      TH1F num4_sublead("num4_sublead_"+category[i],"num4_sublead_"+category[i],35,30.,100.);


       den_lead.Sumw2();
      num1_lead.Sumw2();
      num2_lead.Sumw2();
      num3_lead.Sumw2();
      num4_lead.Sumw2();
       den_sublead.Sumw2();
      num1_sublead.Sumw2();
      num2_sublead.Sumw2();
      num3_sublead.Sumw2();
      num4_sublead.Sumw2();

      AnaTree->Draw("ptphot1>>den_lead_"+category[i], "("+cut_category_1[i]+")*pu_weight");
      AnaTree->Draw("ptphot1>>num1_lead_"+category[i],"("+cut_category_1[i]+" && idcicphot1>0)*pu_weight");
      AnaTree->Draw("ptphot1>>num2_lead_"+category[i],"("+cut_category_1[i]+" && idcicphot1>1)*pu_weight");
      AnaTree->Draw("ptphot1>>num3_lead_"+category[i],"("+cut_category_1[i]+" && idcicphot1>2)*pu_weight");
      AnaTree->Draw("ptphot1>>num4_lead_"+category[i],"("+cut_category_1[i]+" && idcicphot1>3)*pu_weight");
      
      AnaTree->Draw("ptphot2>>den_sublead_"+category[i], "("+cut_category_2[i]+")*pu_weight");
      AnaTree->Draw("ptphot2>>num1_sublead_"+category[i],"("+cut_category_2[i]+" && idcicphot2>0)*pu_weight");
      AnaTree->Draw("ptphot2>>num2_sublead_"+category[i],"("+cut_category_2[i]+" && idcicphot2>1)*pu_weight");
      AnaTree->Draw("ptphot2>>num3_sublead_"+category[i],"("+cut_category_2[i]+" && idcicphot2>2)*pu_weight");
      AnaTree->Draw("ptphot2>>num4_sublead_"+category[i],"("+cut_category_2[i]+" && idcicphot2>3)*pu_weight");

      TH1F* den= (TH1F*) den_lead.Clone("den");
      TH1F* num1=(TH1F*) num1_lead.Clone("num1");
      TH1F* num2=(TH1F*) num2_lead.Clone("num2");
      TH1F* num3=(TH1F*) num3_lead.Clone("num3");
      TH1F* num4=(TH1F*) num4_lead.Clone("num4");
      
      num1_lead.Divide(&den_lead);
      num2_lead.Divide(&den_lead);
      num3_lead.Divide(&den_lead);
      num4_lead.Divide(&den_lead);
      

      a.Draw();      

      num1_lead.SetLineWidth(2);
      num1_lead.Draw("SAME");
      num2_lead.SetLineWidth(2);
      num2_lead.SetLineColor(2);
      num2_lead.Draw("SAME");
      num3_lead.SetLineWidth(2);
      num3_lead.SetLineColor(3);
      num3_lead.Draw("SAME");
      num4_lead.SetLineWidth(2);
      num4_lead.SetLineColor(4);
      num4_lead.Draw("SAME");
      
      TLegend leg(0.8,0.1,0.99,0.3);
      leg.SetFillColor(0);
      leg.SetBorderSize(0);
      leg.AddEntry(&num1_lead,ciclevels[0],"l");
      leg.AddEntry(&num2_lead,ciclevels[1],"l");
      leg.AddEntry(&num3_lead,ciclevels[2],"l");
      leg.AddEntry(&num4_lead,ciclevels[3],"l");
      leg.Draw();
      c1->SaveAs("lead_"+category[i]+"_CiC_eff.png");
      
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
      num1_sublead.SetLineWidth(2);
      num1_sublead.Draw("SAME");
      num2_sublead.SetLineWidth(2);
      num2_sublead.SetLineColor(2);
      num2_sublead.Draw("SAME");
      num3_sublead.SetLineWidth(2);
      num3_sublead.SetLineColor(3);
      num3_sublead.Draw("SAME");
      num4_sublead.SetLineWidth(2);
      num4_sublead.SetLineColor(4);
      num4_sublead.Draw("SAME");
      
      leg.Draw();
      c1->SaveAs("sublead_"+category[i]+"_CiC_eff.png");

      num1->Divide(den);
      num2->Divide(den);
      num3->Divide(den);
      num4->Divide(den);
      
      a.Draw();
      num1->SetLineWidth(2);
      num1->Draw("SAME");
      num2->SetLineWidth(2);
      num2->SetLineColor(2);
      num2->Draw("SAME");
      num3->SetLineWidth(2);
      num3->SetLineColor(3);
      num3->Draw("SAME");
      num4->SetLineWidth(2);
      num4->SetLineColor(4);
      num4->Draw("SAME");
      
      leg.Draw();
      c1->SaveAs(category[i]+"_CiC_eff.png");


      TH1F den_eta_lead(  "den_eta_lead_"+category[i], "den_eta_lead_"+category[i],50,-2.5,2.5);
      TH1F num_eta1_lead("num_eta1_lead_"+category[i],"num_eta1_lead_"+category[i],50,-2.5,2.5);
      TH1F num_eta2_lead("num_eta2_lead_"+category[i],"num_eta2_lead_"+category[i],50,-2.5,2.5);
      TH1F num_eta3_lead("num_eta3_lead_"+category[i],"num_eta3_lead_"+category[i],50,-2.5,2.5);
      TH1F num_eta4_lead("num_eta4_lead_"+category[i],"num_eta4_lead_"+category[i],50,-2.5,2.5);

      TH1F den_eta_sublead(  "den_eta_sublead_"+category[i], "den_eta_sublead_"+category[i],50,-2.5,2.5);
      TH1F num_eta1_sublead("num_eta1_sublead_"+category[i],"num_eta1_sublead_"+category[i],50,-2.5,2.5);
      TH1F num_eta2_sublead("num_eta2_sublead_"+category[i],"num_eta2_sublead_"+category[i],50,-2.5,2.5);
      TH1F num_eta3_sublead("num_eta3_sublead_"+category[i],"num_eta3_sublead_"+category[i],50,-2.5,2.5);
      TH1F num_eta4_sublead("num_eta4_sublead_"+category[i],"num_eta4_sublead_"+category[i],50,-2.5,2.5);

      den_eta_lead.Sumw2();
      num_eta1_lead.Sumw2();
      num_eta2_lead.Sumw2();
      num_eta3_lead.Sumw2();
      num_eta4_lead.Sumw2();
      den_eta_sublead.Sumw2();
      num_eta1_sublead.Sumw2();
      num_eta2_sublead.Sumw2();
      num_eta3_sublead.Sumw2();
      num_eta4_sublead.Sumw2();

      AnaTree->Draw("etaphot1>>den_eta_lead_"+category[i], "("+cut_category_1[i]+")*pu_weight");
      AnaTree->Draw("etaphot1>>num_eta1_lead_"+category[i],"("+cut_category_1[i]+" && idcicphot1>0)*pu_weight");
      AnaTree->Draw("etaphot1>>num_eta2_lead_"+category[i],"("+cut_category_1[i]+" && idcicphot1>1)*pu_weight");
      AnaTree->Draw("etaphot1>>num_eta3_lead_"+category[i],"("+cut_category_1[i]+" && idcicphot1>2)*pu_weight");
      AnaTree->Draw("etaphot1>>num_eta4_lead_"+category[i],"("+cut_category_1[i]+" && idcicphot1>3)*pu_weight");
      
      AnaTree->Draw("etaphot2>>den_eta_sublead_"+category[i], "("+cut_category_2[i]+")*pu_weight");
      AnaTree->Draw("etaphot2>>num_eta1_sublead_"+category[i],"("+cut_category_2[i]+" && idcicphot2>0)*pu_weight");
      AnaTree->Draw("etaphot2>>num_eta2_sublead_"+category[i],"("+cut_category_2[i]+" && idcicphot2>1)*pu_weight");
      AnaTree->Draw("etaphot2>>num_eta3_sublead_"+category[i],"("+cut_category_2[i]+" && idcicphot2>2)*pu_weight");
      AnaTree->Draw("etaphot2>>num_eta4_sublead_"+category[i],"("+cut_category_2[i]+" && idcicphot2>3)*pu_weight");

      TH1F* den_eta= (TH1F*) den_eta_lead.Clone("den_eta");
      TH1F* num_eta1=(TH1F*) num_eta1_lead.Clone("num_eta1");
      TH1F* num_eta2=(TH1F*) num_eta2_lead.Clone("num_eta2");
      TH1F* num_eta3=(TH1F*) num_eta3_lead.Clone("num_eta3");
      TH1F* num_eta4=(TH1F*) num_eta4_lead.Clone("num_eta4");
      
      num_eta1_lead.Divide(&den_eta_lead);
      num_eta2_lead.Divide(&den_eta_lead);
      num_eta3_lead.Divide(&den_eta_lead);
      num_eta4_lead.Divide(&den_eta_lead);
      

      b.Draw();      
      num_eta1_lead.SetLineWidth(2);
      num_eta1_lead.Draw("SAME");
      num_eta2_lead.SetLineWidth(2);
      num_eta2_lead.SetLineColor(2);
      num_eta2_lead.Draw("SAME");
      num_eta3_lead.SetLineWidth(2);
      num_eta3_lead.SetLineColor(3);
      num_eta3_lead.Draw("SAME");
      num_eta4_lead.SetLineWidth(2);
      num_eta4_lead.SetLineColor(4);
      num_eta4_lead.Draw("SAME");
      
      TLegend leg(0.8,0.1,0.99,0.3);
      leg.SetFillColor(0);
      leg.SetBorderSize(0);
      leg.AddEntry(&num_eta1_lead,ciclevels[0],"l");
      leg.AddEntry(&num_eta2_lead,ciclevels[1],"l");
      leg.AddEntry(&num_eta3_lead,ciclevels[2],"l");
      leg.AddEntry(&num_eta4_lead,ciclevels[3],"l");
      leg.Draw();
      c1->SaveAs("lead_eta_"+category[i]+"_CiC_eff.png");
      
      den_eta->Add(&den_eta_sublead);
      num_eta1->Add(&num_eta1_sublead);
      num_eta2->Add(&num_eta2_sublead);
      num_eta3->Add(&num_eta3_sublead);
      num_eta4->Add(&num_eta4_sublead);
      
      num_eta1_sublead.Divide(&den_eta_sublead);
      num_eta2_sublead.Divide(&den_eta_sublead);
      num_eta3_sublead.Divide(&den_eta_sublead);
      num_eta4_sublead.Divide(&den_eta_sublead);
      
      b.Draw();
      num_eta1_sublead.SetLineWidth(2);
      num_eta1_sublead.Draw("SAME");
      num_eta2_sublead.SetLineWidth(2);
      num_eta2_sublead.SetLineColor(2);
      num_eta2_sublead.Draw("SAME");
      num_eta3_sublead.SetLineWidth(2);
      num_eta3_sublead.SetLineColor(3);
      num_eta3_sublead.Draw("SAME");
      num_eta4_sublead.SetLineWidth(2);
      num_eta4_sublead.SetLineColor(4);
      num_eta4_sublead.Draw("SAME");
      
      leg.Draw();
      c1->SaveAs("sublead_eta_"+category[i]+"_CiC_eff.png");

      num_eta1->Divide(den_eta);
      num_eta2->Divide(den_eta);
      num_eta3->Divide(den_eta);
      num_eta4->Divide(den_eta);
      
      b.Draw();
      num_eta1->SetLineWidth(2);
      num_eta1->Draw("SAME");
      num_eta2->SetLineWidth(2);
      num_eta2->SetLineColor(2);
      num_eta2->Draw("SAME");
      num_eta3->SetLineWidth(2);
      num_eta3->SetLineColor(3);
      num_eta3->Draw("SAME");
      num_eta4->SetLineWidth(2);
      num_eta4->SetLineColor(4);
      num_eta4->Draw("SAME");
      
      leg.Draw();
      c1->SaveAs(category[i]+"_eta_CiC_eff.png");


      TH1F  den_leadRhoPF( "den_leadRhoPF_"+category[i], "den_leadRhoPF_"+category[i],60,0.,30.);
      TH1F num1_leadRhoPF("num1_leadRhoPF_"+category[i],"num1_leadRhoPF_"+category[i],60,0.,30.);
      TH1F num2_leadRhoPF("num2_leadRhoPF_"+category[i],"num2_leadRhoPF_"+category[i],60,0.,30.);
      TH1F num3_leadRhoPF("num3_leadRhoPF_"+category[i],"num3_leadRhoPF_"+category[i],60,0.,30.);
      TH1F num4_leadRhoPF("num4_leadRhoPF_"+category[i],"num4_leadRhoPF_"+category[i],60,0.,30.);

       den_leadRhoPF.Sumw2();
      num1_leadRhoPF.Sumw2();
      num2_leadRhoPF.Sumw2();
      num3_leadRhoPF.Sumw2();
      num4_leadRhoPF.Sumw2();
      
      AnaTree->Draw("rhoPF>>den_leadRhoPF_"+category[i], "("+cut_category_1[i]+")*pu_weight");
      AnaTree->Draw("rhoPF>>num1_leadRhoPF_"+category[i],"("+cut_category_1[i]+" && idcicphot1>0)*pu_weight");
      AnaTree->Draw("rhoPF>>num2_leadRhoPF_"+category[i],"("+cut_category_1[i]+" && idcicphot1>1)*pu_weight");
      AnaTree->Draw("rhoPF>>num3_leadRhoPF_"+category[i],"("+cut_category_1[i]+" && idcicphot1>2)*pu_weight");
      AnaTree->Draw("rhoPF>>num4_leadRhoPF_"+category[i],"("+cut_category_1[i]+" && idcicphot1>3)*pu_weight");
      
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
      num1_leadRhoPF.SetLineWidth(2);
      num1_leadRhoPF.Draw("SAME");
      num2_leadRhoPF.SetLineWidth(2);
      num2_leadRhoPF.SetLineColor(2);
      num2_leadRhoPF.Draw("SAME");
      num3_leadRhoPF.SetLineWidth(2);
      num3_leadRhoPF.SetLineColor(3);
      num3_leadRhoPF.Draw("SAME");
      num4_leadRhoPF.SetLineWidth(2);
      num4_leadRhoPF.SetLineColor(4);
      num4_leadRhoPF.Draw("SAME");
      
      TLegend leg1(0.8,0.1,0.99,0.3);
      leg1.SetFillColor(0);
      leg1.SetBorderSize(0);
      leg1.AddEntry(&num1_leadRhoPF,ciclevels[0],"l");
      leg1.AddEntry(&num2_leadRhoPF,ciclevels[1],"l");
      leg1.AddEntry(&num3_leadRhoPF,ciclevels[2],"l");
      leg1.AddEntry(&num4_leadRhoPF,ciclevels[3],"l");
      leg1.Draw();
      
      c1->SaveAs("leadRhoPF_"+category[i]+"_CiC_eff.png");
    }
}
