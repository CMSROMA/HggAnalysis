{
  gROOT->ProcessLine(".L fillPlot.C++");
  gROOT->ProcessLine(".L finalize.C++");

  //finalize
  //(data int lumi 2010, data int lumi 2011, pt1 cut, pt2 cut, ptj1 cut, ptj2 cut, deltae cut, zep cut, mjj cut, eb cat, r9 cat, trk isocut scale, cic selection, no pixel)

  float lumi2011=214;
  float lumi2010=0;
  float ptlead=40.;
  float ptsublead=30.;
  float pthiggscut=40.;
  //CiC Supertight is 4
  int ciclevels[2]={4,0};
  //  int ciclevels[1]={4};

  /*

      _____ _  _____ 
     / ____(_)/ ____|
    | |     _| |     
    | |    | | |     
    | |____| | |____ 
     \_____|_|\_____|

  */                 

  for (int icic=0;icic<1;++icic)
    {
      //
      // Plots in 2PT 2R9 2ETA categories
      //

      std::cout << "******* Doing CiC " << ciclevels[icic] << " 2PT 2R9 2ETA categories plots ******* " << std::endl;

      // Cat0 EBEB R9R9 low pt
      finalize(lumi2010,lumi2011,ptlead,ptsublead,-100,pthiggscut,-10000,-10000,0,0,0,1,1,ciclevels[icic],1);
      // Cat0 EBEB R9R9 low pt
      finalize(lumi2010,lumi2011,ptlead,ptsublead,pthiggscut,-100,-10000,-10000,0,0,0,1,1,ciclevels[icic],1);
      // Cat1 EBEB !R9R9 low pt
      finalize(lumi2010,lumi2011,ptlead,ptsublead,-100,pthiggscut,-10000,-10000,0,0,0,1,0,ciclevels[icic],1);
      // Cat1 EBEB !R9R9 high pt
      finalize(lumi2010,lumi2011,ptlead,ptsublead,pthiggscut,-100,-10000,-10000,0,0,0,1,0,ciclevels[icic],1);
      // Cat2 !EBEB R9R9 low pt
      finalize(lumi2010,lumi2011,ptlead,ptsublead,-100,pthiggscut,-10000,-10000,0,0,0,0,1,ciclevels[icic],1);
      // Cat2 !EBEB R9R9 higgs pt
      finalize(lumi2010,lumi2011,ptlead,ptsublead,pthiggscut,-100,-10000,-10000,0,0,0,0,1,ciclevels[icic],1);
      // Cat3 !EBEB !R9R9 low pt
      finalize(lumi2010,lumi2011,ptlead,ptsublead,-100,pthiggscut,-10000,-10000,0,0,0,0,0,ciclevels[icic],1);
      // Cat3 !EBEB !R9R9 high pt
      finalize(lumi2010,lumi2011,ptlead,ptsublead,pthiggscut,-100,-10000,-10000,0,0,0,0,0,ciclevels[icic],1);

      //
      // Plots in 2R9 2ETA categories
      //
      // Cat0 EBEB R9R9
      std::cout << "******* Doing CiC " << ciclevels[icic] << " 2R9 2ETA categories plots ******* " << std::endl;
      finalize(lumi2010,lumi2011,ptlead,ptsublead,-100,-100,-10000,-10000,0,0,0,1,1,ciclevels[icic],1);
      // Cat1 EBEB !R9R9
      finalize(lumi2010,lumi2011,ptlead,ptsublead,-100,-100,-10000,-10000,0,0,0,1,0,ciclevels[icic],1);
      // Cat2 !EBEB R9R9
      finalize(lumi2010,lumi2011,ptlead,ptsublead,-100,-100,-10000,-10000,0,0,0,0,1,ciclevels[icic],1);
      // Cat3 !EBEB !R9R9
      finalize(lumi2010,lumi2011,ptlead,ptsublead,-100,-100,-10000,-10000,0,0,0,0,0,ciclevels[icic],1);
      
      //
      // Plots in 2ETA categories
      //
      std::cout << "******* Doing CiC " << ciclevels[icic] << " 2ETA categories plots ******* " << std::endl;
      //EBEB NoR9Cat
      finalize(lumi2010,lumi2011,ptlead,ptsublead,-100,-100,-10000,-10000,0,0,0,1,-1,ciclevels[icic],1);
      //!(EBEB) NoR9Cat
      finalize(lumi2010,lumi2011,ptlead,ptsublead,-100,-100,-10000,-10000,0,0,0,0,-1,ciclevels[icic],1);

      //
      // All together
      //
      std::cout << "******* Doing CiC " << ciclevels[icic] << " all together plots ******* " << std::endl;
      finalize(lumi2010,lumi2011,ptlead,ptsublead,-100,-100,-10000,-10000,0,0,0,-1,-1,ciclevels[icic],1);

//      finalize(lumi2010,lumi2011,ptlead,ptsublead,-10000,-10000,0,0,0,-1,-1,ciclevels[icic],1,"nvtx",30,-0.5,29.5,"nRecoVtx");
//      finalize(lumi2010,lumi2011,ptlead,ptsublead,-10000,-10000,0,0,0,-1,-1,ciclevels[icic],1,"npu",30,-0.5,29.5,"nSimVtx");
//       finalize(lumi2010,lumi2011,ptlead,ptsublead,-10000,-10000,0,0,0,-1,-1,ciclevels[icic],1,"ptphot1",75,0.,150.,"pt #Gamma_{1} [GeV]");
//       finalize(lumi2010,lumi2011,ptlead,ptsublead,-10000,-10000,0,0,0,-1,-1,ciclevels[icic],1,"ptphot2",75,0.,150.,"pt #Gamma_{2} [GeV]");
//      finalize(lumi2010,lumi2011,ptlead,ptsublead,-10000,-10000,0,0,0,-1,-1,ciclevels[icic],1,"etaphot1",60,-3.,3.,"#eta #Gamma_{1}");
//       finalize(lumi2010,lumi2011,ptlead,ptsublead,-10000,-10000,0,0,0,-1,-1,ciclevels[icic],1,"etaphot2",60,-3.,3.,"#eta #Gamma_{2}");
    }
}
