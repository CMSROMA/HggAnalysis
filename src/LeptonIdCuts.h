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

struct electronidcuts2012 {

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
  float hoeEB;
  float hoeEE;
  float oemopEB;
  float oemopEE;
  float minhitsEB;
  float minhitsEE;
  float d0EB;
  float d0EE;
  float dzEB;
  float dzEE;
  float iso_relEB;
  float iso_relEE;
};

struct electronidcutsMva2012 {

  float eta;
  float crack1;
  float crack2;
  float pt;
  float mvaCentEB;
  float mvaOutEB;
  float mvaEE;
  float minhitsEB;
  float minhitsEE;
  float d0EB;
  float d0EE;
  float dzEB;
  float dzEE;
  float iso_relCentEB;
  float iso_relOutEB;
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

struct muonidcuts2012 {

  float eta;
  float pt;
  float chi2;
  float hits;
  float match;
  float pixhits;
  float withm;
  float d0;
  float dz;
  float iso_rel;
};

#endif
