#ifndef LeptonIdCuts_H
#define LeptonIdCuts_H

struct electronidcuts {

  float eta;
  float crack1;
  float crack2;
  float pt;
  float setaetaEB;
  float setaetaEE;
  float dphiEB;
  float dphiEE;
  float detaEB;
  float detaEE;
  float minhitsEB;
  float minhitsEE;
  float dcotEB;
  float dcotEE;
  float distEB;
  float distEE;
  float d0EB;
  float d0EE;
  float dzEB;
  float dzEE;
  float iso_relEB;
  float iso_relEE;
};

struct muonidcuts {

  float eta;
  float pt;
  float pixhits;
  float tkhits;
  float hits;
  float chi2;
  float match;
  float d0;
  float dz;
  float iso_rel;
};

#endif
