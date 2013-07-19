{
  gROOT->ProcessLine(".L fillHisto.C++");
  gROOT->ProcessLine(".L optimize.C++");

  //optimize
  //(data int lumi 2010, data int lumi 2011, pt1 cut, pt2 cut, ptj1 cut, ptj2 cut, deltae cut, zep cut, mjj cut, eb cat, r9 cat, trk isocut scale, cic selection, no pixel)

  // float lumi2011=1090;
  float lumi2011=4700;

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

//   float ptleadvbf=60.;
//   float ptsubleadvbf=25.;
//   float ptjetleadvbf=30.;
//   float ptjetsubleadvbf=20.;
//   float deltaetavbf=3.;
//   float zeppenvbf=2.5;
//   float mjjvbf=250.;
// //   float zeppenvbf=2.;
// //   float mjjvbf=400.;
//   float deltaphivbf=2.6;

  float ptleadwzh=60.;
  float ptsubleadwzh=40.;
  float ptjetleadwzh=40.;
  float ptjetsubleadwzh=25.;
  float deltaetawzh=0.0;
  float zeppenwzh=1.5;
  float mjjwzh=-30.;
  float deltaphiwzh=0.0;

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
 
  //  optimize(lumi2011,ptleadvbf,ptsubleadvbf,-10000, -10000, ptjetleadvbf,ptjetsubleadvbf,deltaetavbf,zeppenvbf,mjjvbf,deltaphivbf,-1,-1,ciclevel,0, "massgg",80,100,180);
  //  optimize(lumi2011,ptleadvbf,ptsubleadvbf,-10000, -10000, ptjetleadvbf,ptjetsubleadvbf,deltaetavbf,zeppenvbf,mjjvbf,deltaphivbf,1,-1,ciclevel,0, "massgg",80,100,180);
  //  optimize(lumi2011,ptleadvbf,ptsubleadvbf,-10000, -10000, ptjetleadvbf,ptjetsubleadvbf,deltaetavbf,zeppenvbf,mjjvbf,deltaphivbf,1,1,ciclevel,0, "massgg",80,100,180);
  //  optimize(lumi2011,ptleadvbf,ptsubleadvbf,-10000, -10000, ptjetleadvbf,ptjetsubleadvbf,deltaetavbf,zeppenvbf,mjjvbf,deltaphivbf,0,-1,ciclevel,0, "massgg",80,100,180);
//   optimize(lumi2011,ptleadvbf,ptsubleadvbf,-10000, -10000, ptjetleadvbf,ptjetsubleadvbf,deltaetavbf,zeppenvbf,mjjvbf,deltaphivbf,0,1,ciclevel,0, "massgg",80,100,180);

  optimize(lumi2011,ptleadwzh,ptsubleadwzh,-10000, -10000, ptjetleadwzh,ptjetsubleadwzh,deltaetawzh,zeppenwzh,mjjwzh,deltaphiwzh,-1,-1,ciclevel,0, "massgg",80,100,180);


}
