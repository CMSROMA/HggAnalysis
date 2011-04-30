#ifndef PhotonIdCuts_H
#define PhotonIdCuts_H
//Photon ID DelRe                                                            
struct photonidcuts {
  int tracknb;
  float trackiso_rel;
  float ecaliso_rel;
  float ecaliso_abs;
  float hcaliso_rel;
  float hcaliso_abs;
  float sminmin;
  float sminmin_min;
  float smajmaj;
};

struct photonidegcuts {
  float trackiso_rel;
  float trackiso_abs;
  float ecaliso_rel;
  float ecaliso_abs;
  float hovereiso;
  float hcaliso_rel;
  float hcaliso_abs;
  float setaetaEB;
  float setaetaEE;
};

struct photonidelecuts {
  float trackiso_relEB;
  float ecaliso_relEB;
  float hovereisoEB;
  float hcaliso_relEB;
  float setaetaEB;
  float detaEB;
  float dphiEB;
  int minhitsEB;
  float dcotEB;
  float distEB;
  float trackiso_relEE;
  float ecaliso_relEE;
  float hovereisoEE;
  float hcaliso_relEE;
  float setaetaEE;
  float detaEE;
  float dphiEE;
  int minhitsEE;
  float dcotEE;
  float distEE;
};

#endif
