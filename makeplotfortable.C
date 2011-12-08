



{
  gROOT->ProcessLine(".L fillPlot.C++");
  gROOT->ProcessLine(".L finalize.C++");

  //finalize
  //(data int lumi 2010, data int lumi 2011, pt1 cut, pt2 cut, ptj1 cut, ptj2 cut, deltae cut, zep cut, mjj cut, eb cat, r9 cat, trk isocut scale, cic selection, no pixel)

  // float lumi2011=1090;
  float lumi2011=1660;
  float lumi2010=0;

  float ptleadvbf=55.;
  float ptsubleadvbf=25.;
  float ptjetleadvbf=30.;
  float ptjetsubleadvbf=20.;
  float deltaetavbf=3.5;
  float zeppenvbf=2.;
  float mjjvbf=400.;

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

  finalize(lumi2010,lumi2011,ptleadvbf,ptsubleadvbf,-10000, -10000, ptjetleadvbf,ptjetsubleadvbf,deltaetavbf,zeppenvbf,mjjvbf,-1,-1,ciclevel,0, "massgg",19,92.5,187.5);

  finalize(lumi2010,lumi2011,ptleadvbf,ptsubleadvbf,-10000, -10000, -100000,-10000,0,0,-10000,-1,-1,ciclevel,0, "massgg",19,92.5,187.5);
  finalize(lumi2010,lumi2011,ptleadvbf,ptsubleadvbf,-10000, -10000, ptjetleadvbf,ptjetsubleadvbf,0,0,-10000,-1,-1,ciclevel,0, "massgg",19,92.5,187.5);
  finalize(lumi2010,lumi2011,ptleadvbf,ptsubleadvbf,-10000, -10000, ptjetleadvbf,ptjetsubleadvbf,deltaetavbf,0,-10000,-1,-1,ciclevel,0, "massgg",19,92.5,187.5);
  finalize(lumi2010,lumi2011,ptleadvbf,ptsubleadvbf,-10000, -10000, ptjetleadvbf,ptjetsubleadvbf,deltaetavbf,zeppenvbf,-10000,-1,-1,ciclevel,0, "massgg",19,92.5,187.5);


  finalize(lumi2010,lumi2011,ptleadwzh,ptsubleadwzh,-10000, -10000, ptjetleadwzh,ptjetsubleadwzh,deltaetawzh,zeppenwzh,mjjwzh,-1,-1,ciclevel,0, "massgg",19,92.5,187.5);

  finalize(lumi2010,lumi2011,ptleadwzh,ptsubleadwzh,-10000, -10000, -100000,-10000,0,0,-10000,-1,-1,ciclevel,0, "massgg",19,92.5,187.5);
  finalize(lumi2010,lumi2011,ptleadwzh,ptsubleadwzh,-10000, -10000,  ptjetleadwzh,ptjetsubleadwzh,0,0,-10000,-1,-1,ciclevel,0, "massgg",19,92.5,187.5);
  finalize(lumi2010,lumi2011,ptleadwzh,ptsubleadwzh,-10000, -10000,  ptjetleadwzh,ptjetsubleadwzh,deltaetawzh,0,-10000,-1,-1,ciclevel,0, "massgg",19,92.5,187.5);
  finalize(lumi2010,lumi2011,ptleadwzh,ptsubleadwzh,-10000, -10000,  ptjetleadwzh,ptjetsubleadwzh,deltaetawzh,zeppenwzh,-10000,-1,-1,ciclevel,0, "massgg",19,92.5,187.5);

}
