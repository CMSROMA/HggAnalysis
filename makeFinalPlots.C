{
  gROOT->ProcessLine(".L fillPlot.C++");
  gROOT->ProcessLine(".L finalize.C++");

  //finalize
  //(data int lumi 2010, data int lumi 2011, pt1 cut, pt2 cut, ptj1 cut, ptj2 cut, deltae cut, zep cut, mjj cut, eb cat, r9 cat, trk isocut scale, cic selection, no pixel)

  float lumi2011=214;
  float lumi2010=0;
  float ptlead=40.;
  float ptsublead=30.;
  //CiC Supertight is 4
  int ciclevel=4;

  /*

      _____ _  _____ 
     / ____(_)/ ____|
    | |     _| |     
    | |    | | |     
    | |____| | |____ 
     \_____|_|\_____|

  */                 

  //
  // Plots in 2R9 2ETA categories
  //

  // Cat0 EBEB R9R9

  std::cout << "******* Doing CiC " << ciclevel << " 2R9 2ETA categories plots ******* " << std::endl;
  finalize(lumi2010,lumi2011,ptlead,ptsublead,-10000,-10000,0,0,0,1,1,ciclevel,1);
  // Cat1 EBEB !R9R9
  finalize(lumi2010,lumi2011,ptlead,ptsublead,-10000,-10000,0,0,0,1,0,ciclevel,1);
  // Cat2 !EBEB R9R9
  finalize(lumi2010,lumi2011,ptlead,ptsublead,-10000,-10000,0,0,0,0,1,ciclevel,1);
  // Cat3 !EBEB !R9R9
  finalize(lumi2010,lumi2011,ptlead,ptsublead,-10000,-10000,0,0,0,0,0,ciclevel,1);

  //
  // Plots in 2ETA categories
  //
  std::cout << "******* Doing CiC " << ciclevel << " 2ETA categories plots ******* " << std::endl;
  //EBEB NoR9Cat
  finalize(lumi2010,lumi2011,ptlead,ptsublead,-10000,-10000,0,0,0,1,-1,ciclevel,1);
  //!(EBEB) NoR9Cat
  finalize(lumi2010,lumi2011,ptlead,ptsublead,-10000,-10000,0,0,0,0,-1,ciclevel,1);

  //
  // All together
  //
  std::cout << "******* Doing CiC " << ciclevel << " all together plots ******* " << std::endl;
  finalize(lumi2010,lumi2011,ptlead,ptsublead,-10000,-10000,0,0,0,-1,-1,ciclevel,1);
  
}
