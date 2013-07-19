{
  gROOT->ProcessLine(".L fillPlot.C++");
  gROOT->ProcessLine(".L finalize.C++");
  
  // float lumi2011=4700;
  float lumi2011=5084;
  float lumi2010=0;

  // CiC Supertight is 4
  int ciclevel=4;

  // per met tag:
  finalize(lumi2010,lumi2011,40.,25.,-10000,-10000,-10000,-10000,70.,0,0,-10000,0,0.,0.,0.,30.,170.,-1,-1,ciclevel,0,0,1,"massgg","massgg",200,95.,195.,"m(#gamma#gamma)[GeV]"); 

  // per lepton tag
  // finalize(lumi2010,lumi2011,45.,25.,-10000,-10000,-10000,-10000,0.,0,0,-10000,0,0.,0.,0.,0.,300.,-1,-1,ciclevel,0,1,0,"massgg","massgg",12,100.,160.,"m(#gamma#gamma)[GeV]"); 
}
