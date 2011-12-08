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





   finalize(lumi2010,lumi2011,33,25,-10000, -10000, -100000,-10000,0,0,-10000,0,-1,-1,ciclevel,0, "massgg",40,100,180);
   finalize(lumi2010,lumi2011,33,25,-10000, -10000, -100000,-10000,0,0,-10000,0,1,1,ciclevel,0, "massgg",40,100,180);
   finalize(lumi2010,lumi2011,33,25,-10000, -10000, -100000,-10000,0,0,-10000,0,1,0,ciclevel,0, "massgg",40,100,180);
   finalize(lumi2010,lumi2011,33,25,-10000, -10000, -100000,-10000,0,0,-10000,0,0,1,ciclevel,0, "massgg",40,100,180);
   finalize(lumi2010,lumi2011,33,25,-10000, -10000, -100000,-10000,0,0,-10000,0,0,0,ciclevel,0, "massgg",40,100,180);

}
