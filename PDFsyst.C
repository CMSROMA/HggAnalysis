void PDFsyst(int type){

  // results
  // 120
  // CL10     +0.0077105 -0.0066569
  // con as   +0.0093669 -0.007862
  // MSTW2008 +0.0052567 -0.0037918
  // con  as  +0.0136158 -0.0116928
  // NNPDF2.0 +0.0056882 -0.0054511
  // con as   +0.0063479 -0.0049924


  TFile g("/tmp/delre/glu_fullas.root");
  char allcut[3000];
  char cut1[500],cut2[500],cut2ap[500],cut2am[500];
  sprintf(cut1,"%s","(ptphot1>55&&ptphot2>25&&ptcorrjet1>30&&ptcorrjet2>20&&abs(deltaeta)>3.5&&abs(zeppenjet)<2.5&&invmassjet>350&&abs(deltaphi)>2.6&&abs(etascphot1)<2.5&&abs(etascphot2)<2.5&&idcicphot1>=4&&idcicphot2>=4&&!((TMath::Abs(etascphot1)>1.4442&&TMath::Abs(etascphot1)<1.566)||(TMath::Abs(etascphot2)>1.4442&&TMath::Abs(etascphot2)<1.566))&&massgg>90&&massgg<190)");
  
  TH1D mass("mass","mass",1000,90.,190.);

  double PDF[450];
  double norm[450];

  for(int i=0;i<150; i++){

    norm[i] = 0;
    if(type == 1) {
      norm[i] = nPDFweight1->GetBinContent(i+1);
      sprintf(cut2,"%s%i%s"," * PDFweight1[",i,"]");
      norm[i+150] = nPDFweight2->GetBinContent(i+1);
      sprintf(cut2ap,"%s%i%s"," * PDFweight2[",i,"]");
      norm[i+300] =  0;
   }
    if(type == 2) {
      norm[i] = nPDFweight3->GetBinContent(i+1);
      sprintf(cut2,"%s%i%s"," * PDFweight3[",i,"]");
      norm[i+150] = nPDFweight4->GetBinContent(i+1);
      sprintf(cut2ap,"%s%i%s"," * PDFweight4[",i,"]");
      norm[i+300] = nPDFweight5->GetBinContent(i+1);
      sprintf(cut2am,"%s%i%s"," * PDFweight5[",i,"]");
    }
    if(type == 3) {
      norm[i] = nPDFweight6->GetBinContent(i+1);
      sprintf(cut2,"%s%i%s"," * PDFweight6[",i,"]");
      norm[i+150] = nPDFweight7->GetBinContent(i+1);
      sprintf(cut2ap,"%s%i%s"," * PDFweight7[",i,"]");
      norm[i+300] = nPDFweight8->GetBinContent(i+1);
      sprintf(cut2am,"%s%i%s"," * PDFweight8[",i,"]");
    }

    sprintf(allcut,"%s%s",cut1,cut2);    
    cout << allcut << endl;
    AnaTree->Draw("massgg>>mass",allcut);
    if(norm[i]) PDF[i] = mass.Integral();
    if(norm[i]) cout << i << "  " << PDF[i] << "  " << norm[i] << " " << AnaTree->GetEntries(allcut) << "  " << ptphotgen1->GetEntries() << "  " << PDF[i]/norm[i] << "  " << "  " << (PDF[i]/norm[i]-double(AnaTree->GetEntries(allcut))/norm[0])/ (PDF[i]/norm[i]) << endl;

    sprintf(allcut,"%s%s",cut1,cut2ap);    
    cout << allcut << endl;
    AnaTree->Draw("massgg>>mass",allcut);
    if(norm[i+150]) PDF[i+150] = mass.Integral();
    if(norm[i+150]) cout << i+150 << "  " << PDF[i+150] << "  " << norm[i+150] << " " << AnaTree->GetEntries(allcut) << "  " << ptphotgen1->GetEntries() << "  " << PDF[i+150]/norm[i+150] << "  " << "  " << (PDF[i+150]/norm[i+150]-double(AnaTree->GetEntries(allcut))/norm[0])/ (PDF[i+150]/norm[i+150]) << endl;

    if(type == 1) continue;

    sprintf(allcut,"%s%s",cut1,cut2am);    
    cout << allcut << endl;
    AnaTree->Draw("massgg>>mass",allcut);
    if(norm[i+300]) PDF[i+300] = mass.Integral();
    if(norm[i+300]) cout << i << "  " << PDF[i+300] << "  " << norm[i+300] << " " << AnaTree->GetEntries(allcut) << "  " << ptphotgen1->GetEntries() << "  " << PDF[i+300]/norm[i+300] << "  " << "  " << (PDF[i+300]/norm[i+300]-double(AnaTree->GetEntries(allcut))/norm[0])/ (PDF[i+300]/norm[i+300]) << endl;
    
  }
  
  double totalsystp(0);
  double totalsystm(0);
  int counterp(0);
  int counterm(0);

  for (int i=1; i<450; i++){

    if(norm[i] && (PDF[i]/norm[i] - PDF[0]/norm[0]) > 0) {
      totalsystp = sqrt(totalsystp * totalsystp + ((PDF[i]/norm[i] - PDF[0]/norm[0]) * (PDF[i]/norm[i] - PDF[0]/norm[0]))/((PDF[0]/norm[0])*(PDF[0]/norm[0])));
      counterp++;
      cout << " p " << i << "  " << totalsystp << "  " << counterp << endl;
    }
    if(norm[i] && (PDF[i]/norm[i] - PDF[0]/norm[0]) < 0) {
      totalsystm = sqrt(totalsystm * totalsystm + ((PDF[i]/norm[i] - PDF[0]/norm[0]) * (PDF[i]/norm[i] - PDF[0]/norm[0]))/((PDF[0]/norm[0])*(PDF[0]/norm[0])));
      counterm++;
      cout << " m " << i << "  " << totalsystm << "  " << counterm << endl;      
    }
  }

  if(type == 1) {
    totalsystp = totalsystp / 1.645;
    totalsystm = totalsystm / 1.645;
  }
  else if(type == 3) {
    totalsystp = totalsystp / sqrt(counterp-1);
    totalsystm = totalsystm / sqrt(counterm-1);
  }

  cout << "+" << totalsystp << " -" << totalsystm << endl;

}
