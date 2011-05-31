{
  //  TFile f("redntp.42xv1.preselection.v1/redntp_GluGluToHToGG_M-115_7TeV-powheg-pythia6_00.root");
  TFile f(" redntp_GluGluToHToGG_M-115_7TeV-powheg-pythia6_00.root");

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

  AnaTree->Draw("ptphot1>>den_lead","abs(etascphot1)>1.5 && ptphot1>30. && abs(etascphot1)<2.5");
  AnaTree->Draw("ptphot1>>num1_lead","abs(etascphot1)>1.5 && ptphot1>30. && abs(etascphot1)<2.5 && idcicphot1>0");
  AnaTree->Draw("ptphot1>>num2_lead","abs(etascphot1)>1.5 && ptphot1>30. && abs(etascphot1)<2.5 && idcicphot1>1");
  AnaTree->Draw("ptphot1>>num3_lead","abs(etascphot1)>1.5 && ptphot1>30. && abs(etascphot1)<2.5 && idcicphot1>2");
  AnaTree->Draw("ptphot1>>num4_lead","abs(etascphot1)>1.5 && ptphot1>30. && abs(etascphot1)<2.5 && idcicphot1>3");

  AnaTree->Draw("ptphot2>>den_sublead","abs(etascphot2)>1.5 && ptphot2>30. && abs(etascphot2)<2.5");
  AnaTree->Draw("ptphot2>>num1_sublead","abs(etascphot2)>1.5 && ptphot2>30. && abs(etascphot2)<2.5 && idcicphot2>0");
  AnaTree->Draw("ptphot2>>num2_sublead","abs(etascphot2)>1.5 && ptphot2>30. && abs(etascphot2)<2.5 && idcicphot2>1");
  AnaTree->Draw("ptphot2>>num3_sublead","abs(etascphot2)>1.5 && ptphot2>30. && abs(etascphot2)<2.5 && idcicphot2>2");
  AnaTree->Draw("ptphot2>>num4_sublead","abs(etascphot2)>1.5 && ptphot2>30. && abs(etascphot2)<2.5 && idcicphot2>3");

  num1_lead.Divide(&den_lead);
  num2_lead.Divide(&den_lead);
  num3_lead.Divide(&den_lead);
  num4_lead.Divide(&den_lead);
 
  TH2F a("a","a",10,30.,100.,10,0.,1.);
  a.Draw();
  num1_lead.Draw("SAME");
  num2_lead.SetLineColor(2);
  num2_lead.Draw("SAME");
  num3_lead.SetLineColor(3);
  num3_lead.Draw("SAME");
  num4_lead.SetLineColor(4);
  num4_lead.Draw("SAME");

  c1->SaveAs("lead_EE_CiC_eff.png");

  num1_sublead.Divide(&den_sublead);
  num2_sublead.Divide(&den_sublead);
  num3_sublead.Divide(&den_sublead);
  num4_sublead.Divide(&den_sublead);
 
  TH2F a("a","a",10,30.,100.,10,0.,1.);
  a.Draw();
  num1_sublead.Draw("SAME");
  num2_sublead.SetLineColor(2);
  num2_sublead.Draw("SAME");
  num3_sublead.SetLineColor(3);
  num3_sublead.Draw("SAME");
  num4_sublead.SetLineColor(4);
  num4_sublead.Draw("SAME");

  c1->SaveAs("sublead_EE_CiC_eff.png");

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


  AnaTree->Draw("ptphot1>>den_lead", "abs(etascphot1)<1.479 && ptphot1>30. ");
  AnaTree->Draw("ptphot1>>num1_lead","abs(etascphot1)<1.479 && ptphot1>30.  && idcicphot1>0");
  AnaTree->Draw("ptphot1>>num2_lead","abs(etascphot1)<1.479 && ptphot1>30.  && idcicphot1>1");
  AnaTree->Draw("ptphot1>>num3_lead","abs(etascphot1)<1.479 && ptphot1>30.  && idcicphot1>2");
  AnaTree->Draw("ptphot1>>num4_lead","abs(etascphot1)<1.479 && ptphot1>30.  && idcicphot1>3");

  AnaTree->Draw("ptphot2>>den_sublead", "abs(etascphot2)<1.479 && ptphot2>30. ");
  AnaTree->Draw("ptphot2>>num1_sublead","abs(etascphot2)<1.479 && ptphot2>30.  && idcicphot2>0");
  AnaTree->Draw("ptphot2>>num2_sublead","abs(etascphot2)<1.479 && ptphot2>30.  && idcicphot2>1");
  AnaTree->Draw("ptphot2>>num3_sublead","abs(etascphot2)<1.479 && ptphot2>30.  && idcicphot2>2");
  AnaTree->Draw("ptphot2>>num4_sublead","abs(etascphot2)<1.479 && ptphot2>30.  && idcicphot2>3");

  num1_lead.Divide(&den_lead);
  num2_lead.Divide(&den_lead);
  num3_lead.Divide(&den_lead);
  num4_lead.Divide(&den_lead);
 
  TH2F a("a","a",10,30.,100.,10,0.,1.);
  a.Draw();
  num1_lead.Draw("SAME");
  num2_lead.SetLineColor(2);
  num2_lead.Draw("SAME");
  num3_lead.SetLineColor(3);
  num3_lead.Draw("SAME");
  num4_lead.SetLineColor(4);
  num4_lead.Draw("SAME");

  c1->SaveAs("lead_EB_CiC_eff.png");

  num1_sublead.Divide(&den_sublead);
  num2_sublead.Divide(&den_sublead);
  num3_sublead.Divide(&den_sublead);
  num4_sublead.Divide(&den_sublead);
 
  TH2F a("a","a",10,30.,100.,10,0.,1.);
  a.Draw();
  num1_sublead.Draw("SAME");
  num2_sublead.SetLineColor(2);
  num2_sublead.Draw("SAME");
  num3_sublead.SetLineColor(3);
  num3_sublead.Draw("SAME");
  num4_sublead.SetLineColor(4);
  num4_sublead.Draw("SAME");

  c1->SaveAs("sublead_EB_CiC_eff.png");

}
