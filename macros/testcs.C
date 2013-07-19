{
  gROOT->ProcessLine(".L fillPlot.C++");
  gROOT->ProcessLine(".L finalize.C++");

  //finalize
  //(data int lumi 2010, data int lumi 2011, pt1 cut, pt2 cut, ptj1 cut, ptj2 cut, deltae cut, zep cut, mjj cut, eb cat, r9 cat, trk isocut scale, cic selection, no pixel)

  // float lumi2011=1090;
  float lumi2011=4700;
  float lumi2010=0;

  float ptleadvbf=55.;
  float ptsubleadvbf=25.;
  float ptjetleadvbf=30.;
  float ptjetsubleadvbf=20.;
  float deltaetavbf=3.5;
  float zeppenvbf=2.5;
  float mjjvbf=350.;
//   float zeppenvbf=2.;
//   float mjjvbf=400.;
  float deltaphivbf=2.6;

  float ptleadwzh=60.;
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

//   finalize(lumi2010,lumi2011,ptleadvbf,ptsubleadvbf,-10000, -10000, ptjetleadvbf,ptjetsubleadvbf,deltaetavbf,zeppenvbf,mjjvbf,deltaphivbf,-1,-1,ciclevel,0, "massgg",40,100,180   );

  //   finalize(lumi2010,lumi2011,ptleadvbf,ptsubleadvbf,-10000, -10000, -100000,-10000,0,0,-10000,0,-1,-1,ciclevel,0, "massgg",40,100,180);
//    finalize(lumi2010,lumi2011,ptleadvbf,ptsubleadvbf,-10000, -10000, ptjetleadvbf,-100000,0,0,-10000,0,-1,-1,ciclevel,0, "massgg",40,100,180);
//    finalize(lumi2010,lumi2011,ptleadvbf,ptsubleadvbf,-10000, -10000, ptjetleadvbf,ptjetsubleadvbf,0,0,-10000,0,-1,-1,ciclevel,0, "massgg",40,100,180);
//    finalize(lumi2010,lumi2011,ptleadvbf,ptsubleadvbf,-10000, -10000, ptjetleadvbf,ptjetsubleadvbf,deltaetavbf,0,-10000,0,-1,-1,ciclevel,0, "massgg",40,100,180);
//    finalize(lumi2010,lumi2011,ptleadvbf,ptsubleadvbf,-10000, -10000, ptjetleadvbf,ptjetsubleadvbf,deltaetavbf,zeppenvbf,-10000.,0,-1,-1,ciclevel,0, "massgg",40,100,180);
//    finalize(lumi2010,lumi2011,ptleadvbf,ptsubleadvbf,-10000, -10000, ptjetleadvbf,ptjetsubleadvbf,deltaetavbf,zeppenvbf,mjjvbf ,0,-1,-1,ciclevel,0, "massgg",40,100,180);


//    finalize(lumi2010,lumi2011,ptleadvbf,ptsubleadvbf,-10000, -10000, -100000,-10000,0,0,-10000,0,-1,-1,ciclevel,0, "ptphot1",40,55.,200.,"p_{T}(#gamma1)[GeV]");
//    finalize(lumi2010,lumi2011,ptleadvbf,ptsubleadvbf,-10000, -10000, -100000,-10000,0,0,-10000,0,-1,-1,ciclevel,0, "ptphot2",40,25.,150.,"p_{T}(#gamma2)[GeV]");
//    finalize(lumi2010,lumi2011,ptleadvbf,ptsubleadvbf,-10000, -10000, -100000,-10000,0,0,-10000,0,-1,-1,ciclevel,0, "ptgg",40,0.,300.,"p_{T}(#gamma#gamma)[GeV]");
   finalize(lumi2010,lumi2011,ptleadvbf,ptsubleadvbf,-10000, -10000, ptjetleadvbf,ptjetsubleadvbf,0,0,-10000,0,-1,-1,ciclevel,0, "ptphot1",40,55.,200.,"p_{T}(#gamma1)[GeV]");
   finalize(lumi2010,lumi2011,ptleadvbf,ptsubleadvbf,-10000, -10000, ptjetleadvbf,ptjetsubleadvbf,0,0,-10000,0,-1,-1,ciclevel,0, "ptphot2",40,25.,150.,"p_{T}(#gamma2)[GeV]");
//    finalize(lumi2010,lumi2011,ptleadvbf,ptsubleadvbf,-10000, -10000, ptjetleadvbf,ptjetsubleadvbf,0,0,-10000,0,-1,-1,ciclevel,0, "ptgg",40,0.,300.,"p_{T}(#gamma#gamma)[GeV]");
//    finalize(lumi2010,lumi2011,ptleadvbf,ptsubleadvbf,-10000, -10000, -100000,-10000,0,0,-10000,0,-1,-1,ciclevel,0, "ptjet1",40,10.,200.,"p_{T}(jet1)[GeV]");
//    finalize(lumi2010,lumi2011,ptleadvbf,ptsubleadvbf,-10000, -10000, ptjetleadvbf,-10000,0,0,-10000,0,-1,-1,ciclevel,0, "ptjet2",40,10.,200.,"p_{T}(jet2)[GeV]");
//    finalize(lumi2010,lumi2011,ptleadvbf,ptsubleadvbf,-10000, -10000, ptjetleadvbf,ptjetsubleadvbf,0,0,mjjvbf,0,-1,-1,ciclevel,0, "deltaeta",28,0.,7.,"#Delta#eta");
//    finalize(lumi2010,lumi2011,ptleadvbf,ptsubleadvbf,-10000, -10000, ptjetleadvbf,ptjetsubleadvbf,deltaetavbf,0,mjjvbf,0,-1,-1,ciclevel,0, "zeppenjet",32,0.,8.,"zeppen");
//    finalize(lumi2010,lumi2011,ptleadvbf,ptsubleadvbf,-10000, -10000, ptjetleadvbf,ptjetsubleadvbf,0,0,-10000,0,-1,-1,ciclevel,0, "invmassjet",40,0.,1000.,"m(jj)[GeV]");
//    finalize(lumi2010,lumi2011,ptleadvbf,ptsubleadvbf,-10000, -10000, ptjetleadvbf,ptjetsubleadvbf,deltaetavbf,zeppenvbf,mjjvbf,0,-1,-1,ciclevel,0, "deltaphi",16,0,3.143,"#Delta#phi");
//    finalize(lumi2010,lumi2011,ptleadvbf,ptsubleadvbf,-10000, -10000, ptjetleadvbf,ptjetsubleadvbf,0,0,-10000,0,-1,-1,ciclevel,0, "etajet1",35,0.,7.,"|#eta|(jet1)");
//    finalize(lumi2010,lumi2011,ptleadvbf,ptsubleadvbf,-10000, -10000, ptjetleadvbf,ptjetsubleadvbf,0,0,-10000,0,-1,-1,ciclevel,0, "etajet2",35,0.,7.,"|#eta|(jet2)");


   finalize(lumi2010,lumi2011,ptleadwzh,ptsubleadwzh,-10000, -10000, ptjetleadwzh,ptjetsubleadwzh,deltaetawzh,zeppenwzh,mjjwzh,0,-1,-1,ciclevel,0, "massgg",40,100,180);

   finalize(lumi2010,lumi2011,ptleadwzh,ptsubleadwzh,-10000, -10000, -100000,-10000,0,0,-10000,0,-1,-1,ciclevel,0, "massgg",40,100,180);
   finalize(lumi2010,lumi2011,ptleadwzh,ptsubleadwzh,-10000, -10000,  ptjetleadwzh,-10000,0,0,-10000,0,-1,-1,ciclevel,0, "massgg",40,100,180);
   finalize(lumi2010,lumi2011,ptleadwzh,ptsubleadwzh,-10000, -10000,  ptjetleadwzh,ptjetsubleadwzh,0,0,-10000,0,-1,-1,ciclevel,0, "massgg",40,100,180);
   finalize(lumi2010,lumi2011,ptleadwzh,ptsubleadwzh,-10000, -10000,  ptjetleadwzh,ptjetsubleadwzh,deltaetawzh,0,-10000,0,-1,-1,ciclevel,0, "massgg",40,100,180);
   finalize(lumi2010,lumi2011,ptleadwzh,ptsubleadwzh,-10000, -10000,  ptjetleadwzh,ptjetsubleadwzh,deltaetawzh,zeppenwzh,-10000,0,-1,-1,ciclevel,0, "massgg",40,100,180);

//    finalize(lumi2010,lumi2011,ptleadwzh,ptsubleadwzh,-10000, -10000, -100000,-10000,0,0,-10000,0,-1,-1,ciclevel,0, "ptjet1",40,10.,200.,"p_{T}(jet1)[GeV]");
//    finalize(lumi2010,lumi2011,ptleadwzh,ptsubleadwzh,-10000, -10000, ptjetleadwzh,-10000,0,0,-10000,0,-1,-1,ciclevel,0, "ptjet2",40,10.,200.,"p_{T}(jet2)[GeV]");
//    finalize(lumi2010,lumi2011,ptleadwzh,ptsubleadwzh,-10000, -10000, ptjetleadwzh,ptjetsubleadwzh,0,0,mjjwzh,0,-1,-1,ciclevel,0, "deltaeta",28,0.,7.,"#Delta#eta");
//    finalize(lumi2010,lumi2011,ptleadwzh,ptsubleadwzh,-10000, -10000, ptjetleadwzh,ptjetsubleadwzh,deltaetawzh,0,mjjwzh,0,-1,-1,ciclevel,0, "zeppenjet",32,0.,8.,"zeppen");
//    finalize(lumi2010,lumi2011,ptleadwzh,ptsubleadwzh,-10000, -10000, ptjetleadwzh,ptjetsubleadwzh,0,0,-10000,0,-1,-1,ciclevel,0, "invmassjet",20,40.,140.,"m(jj)[GeV]");

//    finalize(lumi2010,lumi2011,40,30,-10000, -10000, -100000,-10000,0,0,-10000,0,-1,-1,ciclevel,0, "massgg",40,100,180);

}
