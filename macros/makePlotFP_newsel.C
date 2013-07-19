
{
  gROOT->ProcessLine(".L fillPlot.C++");
  gROOT->ProcessLine(".L finalize.C++");

  //finalize
  //(data int lumi 2010, data int lumi 2011, pt1 cut, pt2 cut, ptj1 cut, ptj2 cut, deltae cut, zep cut, mjj cut, eb cat, r9 cat, trk isocut scale, cic selection, no pixel)

  // float lumi2011=1090;
//   float lumi2011=1660;
  float lumi2011=4660;
  float lumi2010=0;

  float ptleadvbf=55.;
  float ptsubleadvbf=25.;
  float ptjetleadvbf=30.;
  float ptjetsubleadvbf=20.;
  float deltaetavbf=3.5;
  float zeppenvbf=2.5;
  float mjjvbf=300.;

  float ptleadwzh=65.;
  float ptsubleadwzh=25.;
  float ptjetleadwzh=35.;
  float ptjetsubleadwzh=20.;
  float deltaetawzh=-2.5;
  float zeppenwzh=1.5;
  float mjjwzh=-30.;

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

  //
  // All together
  //
                                                                                                 
  // All together                                                                                                                                                                                                 
  //         

//   finalize(lumi2010,lumi2011,40,30,-10000, -10000, 0,0,0,0,0,-1,-1,ciclevel,1, "massgg",50,80,180);
                                                                                                                                                                                                     
//   std::cout << "******* Doing CiC " << ciclevel << " all together plots for VBF ******* " << std::endl;
  finalize(lumi2010,lumi2011,ptleadvbf,ptsubleadvbf,-10000, -10000, ptjetleadvbf,ptjetsubleadvbf,deltaetavbf,zeppenvbf,mjjvbf,-1,-1,ciclevel,0, "massgg",50,90.,190.);
//   //  finalize(lumi2010,lumi2011,ptleadvbf,ptsubleadvbf,-10000, -10000, ptjetleadvbf,ptjetsubleadvbf,deltaetavbf,zeppenvbf,mjjvbf,-1,-1,ciclevel,0, "massgg",19,92.5,187.5);

//   finalize(lumi2010,lumi2011,-10000,ptsubleadvbf,-10000, -10000, -100000,-10000,0,0,-10000,-1,-1,ciclevel,0, "ptphot1",10,45.,145.,"p_{T}(#gamma1)[GeV]");
//   finalize(lumi2010,lumi2011,ptleadvbf,-10000,-10000, -10000, -100000,-10000,0,0,-10000,-1,-1,ciclevel,0, "ptphot2",10,25.,125.,"p_{T}(#gamma2)[GeV]");
   finalize(lumi2010,lumi2011,ptleadvbf,ptsubleadvbf,-10000, -10000, -100000,-10000,0,0,-10000,-1,-1,ciclevel,0, "ptjet1",10,25.,200.,"p_{T}(jet1)[GeV]");
   finalize(lumi2010,lumi2011,ptleadvbf,ptsubleadvbf,-10000, -10000, -100000,-10000,0,0,-10000,-1,-1,ciclevel,0, "ptjet2",10,25.,200.,"p_{T}(jet2)[GeV]");
   finalize(lumi2010,lumi2011,ptleadvbf,ptsubleadvbf,-10000, -10000, ptjetleadvbf,ptjetsubleadvbf,0,0,mjjvbf,-1,-1,ciclevel,0, "deltaeta",15,0.,9.,"#Delta#eta");
   finalize(lumi2010,lumi2011,ptleadvbf,ptsubleadvbf,-10000, -10000, ptjetleadvbf,ptjetsubleadvbf,deltaetavbf,0,mjjvbf,-1,-1,ciclevel,0, "zeppenjet",15,0.,9.,"zeppen");
   finalize(lumi2010,lumi2011,ptleadvbf,ptsubleadvbf,-10000, -10000, ptjetleadvbf,ptjetsubleadvbf,0,0,-10000,-1,-1,ciclevel,0, "invmassjet",8,0.,1600.,"m(jj)[GeV]");
//   finalize(lumi2010,lumi2011,ptleadvbf,ptsubleadvbf,-10000, -10000, ptjetleadvbf,ptjetsubleadvbf,0,0,-10000,-1,-1,ciclevel,0, "etajet1",35,0.,7.,"|#eta|(jet1)");
//   finalize(lumi2010,lumi2011,ptleadvbf,ptsubleadvbf,-10000, -10000, ptjetleadvbf,ptjetsubleadvbf,0,0,-10000,-1,-1,ciclevel,0, "etajet2",35,0.,7.,"|#eta|(jet2)");

// //   //
// //   // Plots in 2ETA categories
// //   //
// //   std::cout << "******* Doing CiC " << ciclevel << " 2ETA categories plots for VBF ******* " << std::endl;
//   //EBEB NoR9Cat
//   finalize(lumi2010,lumi2011,ptleadvbf,ptsubleadvbf,-10000, -10000,ptjetleadvbf,ptjetsubleadvbf,deltaetavbf,zeppenvbf,mjjvbf,1,-1,ciclevel,0, "massgg",50,90.,190.);
//   //!(EBEB) NoR9Cat
//   finalize(lumi2010,lumi2011,ptleadvbf,ptsubleadvbf,-10000, -10000,ptjetleadvbf,ptjetsubleadvbf,deltaetavbf,zeppenvbf,mjjvbf,0,-1,ciclevel,0, "massgg",50,90.,190.);

//   std::cout << "******* Doing CiC " << ciclevel << " all together plots for WZH ******* " << std::endl;
//  finalize(lumi2010,lumi2011,ptleadwzh,ptsubleadwzh,-10000, -10000, ptjetleadwzh,ptjetsubleadwzh,deltaetawzh,zeppenwzh,mjjwzh,-1,-1,ciclevel,0, "massgg",50,90,190);
//   finalize(lumi2010,lumi2011,ptleadwzh,ptsubleadwzh,-10000, -10000, ptjetleadwzh,ptjetsubleadwzh,deltaetawzh,zeppenwzh,mjjwzh,-1,-1,ciclevel,0, "massgg",19,92.5,187.5);

//  finalize(lumi2010,lumi2011,-10000,ptsubleadwzh,-10000, -10000, -100000,-10000,0,0,-10000,-1,-1,ciclevel,0, "ptphot1",10,45.,145.,"p_{T}(#gamma1)[GeV]");
//   finalize(lumi2010,lumi2011,ptleadwzh,-10000,-10000, -10000, -100000,-10000,0,0,-10000,-1,-1,ciclevel,0, "ptphot2",10,25.,125.,"p_{T}(#gamma2)[GeV]");
//   finalize(lumi2010,lumi2011,ptleadwzh,ptsubleadwzh,-10000, -10000, -100000,-10000,0,0,-10000,-1,-1,ciclevel,0, "ptjet1",10,25.,200.,"p_{T}(jet1)[GeV]");
//   finalize(lumi2010,lumi2011,ptleadwzh,ptsubleadwzh,-10000, -10000, -100000,-10000,0,0,-10000,-1,-1,ciclevel,0, "ptjet2",10,25.,200.,"p_{T}(jet2)[GeV]");
//   finalize(lumi2010,lumi2011,ptleadwzh,ptsubleadwzh,-10000, -10000, ptjetleadwzh,ptjetsubleadwzh,0,0,-10000,-1,-1,ciclevel,0, "deltaeta",15,-9.,9.,"#Delta#eta");
//   finalize(lumi2010,lumi2011,ptleadwzh,ptsubleadwzh,-10000, -10000, ptjetleadwzh,ptjetsubleadwzh,0,0,-10000,-1,-1,ciclevel,0, "zeppenjet",15,-9.,9.,"zeppen");
//   finalize(lumi2010,lumi2011,ptleadwzh,ptsubleadwzh,-10000, -10000, ptjetleadwzh,ptjetsubleadwzh,deltaetawzh,zeppenwzh,-10000,-1,-1,ciclevel,0, "invmassjet",10,0.,300.,"m(jj)[GeV]");

  

  
//   std::cout << "******* Doing CiC " << ciclevel << " 2ETA categories plots for WZH ******* " << std::endl;
//   //EBEB NoR9Cat
//   finalize(lumi2010,lumi2011,ptleadwzh,ptsubleadwzh,-10000, -10000,ptjetleadwzh,ptjetsubleadwzh,deltaetawzh,zeppenwzh,mjjwzh,1,-1,ciclevel,0, "massgg",19,92.5,187.5);
//   //!(EBEB) NoR9Cat
//   finalize(lumi2010,lumi2011,ptleadwzh,ptsubleadwzh,-10000, -10000,ptjetleadwzh,ptjetsubleadwzh,deltaetawzh,zeppenwzh,mjjwzh,0,-1,ciclevel,0, "massgg",19,92.5,187.5);
  
}
