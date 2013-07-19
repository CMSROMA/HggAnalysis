#ifndef EnergyScaleCorrection_hh
#define EnergyScaleCorrection_hh
#include <TString.h>
#include <iostream>
#include <map>
#include <math.h>
#include<fstream>

typedef std::map < TString, std::pair<double, double> > correction_map_t;
typedef std::map < TString, std::pair<double, double> > smearing_map_t;
typedef struct {
  int runMax;
  correction_map_t correction_map;
  smearing_map_t smearing_map;
} correction_t;

class EnergyScaleCorrection{

 public:
  std::map <int, correction_t > runMin_map;
  bool isHggCat;
  bool noCorrections;
  TString correctionType;
 public:
  //  EnergyScaleCorrection(TString correctionFileName, bool isHggCat_=false);
  EnergyScaleCorrection(TString correctionFileName, TString correctionType_);
  ~EnergyScaleCorrection(void);
  void Add(TString category_, int runMin_, int runMax_, double deltaP_, double err_deltaP_, double s, double es);
  void ReadFromFile(TString filename);
  float getScaleOffset(int runNumber, bool isEBEle, double R9Ele, double etaSCEle);
  float getSmearing(int runNumber, bool isEBEle, double R9Ele, double etaSCEle);

 private:
  TString GetElectronCategory(bool isEBEle, double R9Ele, double etaSCEle);
  std::map< int, correction_t>::const_iterator runCorrection_itr;
  std::map< int, correction_t>::const_iterator FindRunCorrection_itr(int runNumber);
  
};

#endif
