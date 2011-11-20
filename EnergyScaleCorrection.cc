#include "EnergyScaleCorrection.h"
//#define DEBUG

//EnergyScaleCorrection::EnergyScaleCorrection(TString correctionFileName, bool isHggCat_):
EnergyScaleCorrection::EnergyScaleCorrection(TString correctionFileName, TString correctionType_):
  runMin_map(), noCorrections(false), correctionType(correctionType_)
{
  if(correctionFileName.Sizeof()>1) noCorrections=false;
  else noCorrections=true;

#ifdef DEBUG
  std::cout << "[DEBUG] Using correctionType: " << correctionType << std::endl;
#endif

  if(correctionType.CompareTo("noCalib")==0){
    noCorrections=true;
  }

  if(correctionType.Contains("Hgg")){
    isHggCat=true;
    noCorrections=false;
  }

  if(!noCorrections) ReadFromFile(correctionFileName);
  runCorrection_itr=runMin_map.begin();

  return;
}

EnergyScaleCorrection::~EnergyScaleCorrection(void){
  return;
}

TString EnergyScaleCorrection::GetElectronCategory(bool isEBEle, double R9Ele, double etaSCEle){

  TString category;

  if(isEBEle) category="EB";
  else category="EE";
  
  if(isHggCat){
    if(correctionType.Contains("Hgg")){
      if(correctionType.Contains("eta")){
	
	//    if(correctionType.CompareTo("Hgg_eta_runNumber")==0){
	if(isEBEle){
	  if (fabs(etaSCEle)>1.) category+="highEta";
	  else category+="lowEta";
	} else {
	  if (fabs(etaSCEle) > 2) category+="highEta";
	  else category+="lowEta";
	}
      }
      if(correctionType.Contains("residual")||correctionType.Contains("R9")){
	// 	   if(R9Ele>0.94) category+="_gold";
	// 	   else category+="_bad";
	if(R9Ele>0.94) category+="Gold";
	else category+="Bad";
      }
    }
  }
#ifdef DEBUG
  //  std::cout << "[DEBUG] category: " << category << std::endl;
#endif  
  return category;
  
}

std::map< int, correction_t>::const_iterator EnergyScaleCorrection::FindRunCorrection_itr(int runNumber){
  std::map< int, correction_t>::const_iterator itr=runMin_map.begin();

  while (runNumber >= itr->first && itr!=runMin_map.end()){
#ifdef DEBUG
    std::cout << "runNumber - itr->first: " << runNumber << " - " << itr->first << std::endl;
#endif
    itr++;
  }
  itr--;
#ifdef DEBUG
  std::cout << "[DEBUG] runNumber - [itr->first: itr->second.runMax] " << runNumber << " - [" << itr->first << ":" << itr->second.runMax << "]" << std::endl;
#endif

  if(itr->second.runMax < runNumber){
    //    std::cerr << "[WARNING] runNumber > runs in recalib file: " << runNumber << std::endl;
    return runMin_map.end();
  }
  //  if(runNumber <= itr->second.runMax) // c'e' un altro errore ma non e' possibile
  return itr;
}

  
float EnergyScaleCorrection::getScaleOffset(int runNumber, bool isEBEle, double R9Ele, double etaSCEle){
  if(noCorrections) return 1;
  if(runMin_map.empty()){
    std::cerr << "[ERROR] runMin_map empty and calibration required" << std::endl;
    std::cerr << "[ERROR] turning into noCorrection mode" << std::endl;
    noCorrections=true;
    return 1;
  }


  if(!(runNumber >= runCorrection_itr->first && runNumber <= runCorrection_itr->second.runMax)){
#ifdef DEBUG
    std::cout << "[DEBUG] starting FindRunCorrection: " << runNumber << "\t" << runCorrection_itr->first << "\t" << runCorrection_itr->second.runMax << std::endl;
#endif    
    runCorrection_itr=FindRunCorrection_itr(runNumber);
    if (runCorrection_itr==runMin_map.end()){
      runCorrection_itr=runMin_map.begin(); // lo reinizializzo per i prossimi run
      return 0; // se per esempio il runNumber non e' definito nel recalib file
    }
  }

  //  if (runMin_map.count(runNumber)==1){
#ifdef DDEBUG
  std::cout << "[DEBUG] " << runNumber << "\t" << runCorrection_itr->first << "\t" << runCorrection_itr->second.runMax << std::endl;
#endif
    
    const correction_map_t *correction_map_p = &(runCorrection_itr->second.correction_map);
    
#ifdef HggCatTry
  // questo se voglio provare a cercare prima le HggCat nel file 
  // e poi se non c'e' usare quelle EB/EE standard
  if(correction_map.count(GetElectronCategory(isEBEle, R9Ele, true))==1)
    return correction_map[GetElectronCategory(isEBEle, R9Ele, true)].first;
  else if(correction_map.count(GetElectronCategory(isEBEle, R9Ele, false))==1)
    return correction_map[GetElectronCategory(isEBEle, R9Ele, false)].first;
  else {
    std::cerr << "[ERROR] Electron category not found!!!!  Scale offset not applied" << std::endl;
    return 1;
  }
#endif
  TString category=GetElectronCategory(isEBEle, R9Ele, etaSCEle);
#ifdef DDEBUG
  std::cout << "[DEBUG] category(isEBEle, R9Ele) = " << category << "(" << isEBEle << ", " << R9Ele << ")" <<   correction_map_p->find(category)->second.first << std::endl;

#endif

  if(correction_map_p->count(category)==1)
    return correction_map_p->find(category)->second.first;
  else{
    std::cerr << "[ERROR] Electron category not found!!!!  Scale offset not applied" << std::endl;
    return 1;
  }

  
}

float EnergyScaleCorrection::getSmearing(int runNumber, bool isEBEle, double R9Ele, double etaSCEle){
  if(noCorrections) return 1;
  if(runMin_map.empty()){
    std::cerr << "[ERROR] runMin_map empty and calibration required" << std::endl;
    std::cerr << "[ERROR] turning into noCorrection mode" << std::endl;
    noCorrections=true;
    return 1;
  }


  if(!(runNumber >= runCorrection_itr->first && runNumber <= runCorrection_itr->second.runMax)){
#ifdef DEBUG
    std::cout << "[DEBUG] starting FindRunCorrection: " << runNumber << "\t" << runCorrection_itr->first << "\t" << runCorrection_itr->second.runMax << std::endl;
#endif    
    runCorrection_itr=FindRunCorrection_itr(runNumber);
    if (runCorrection_itr==runMin_map.end()){
      runCorrection_itr=runMin_map.begin(); // lo reinizializzo per i prossimi run
      return 0; // se per esempio il runNumber non e' definito nel recalib file
    }
  }

  //  if (runMin_map.count(runNumber)==1){
#ifdef DDEBUG
  std::cout << "[DEBUG] " << runNumber << "\t" << runCorrection_itr->first << "\t" << runCorrection_itr->second.runMax << std::endl;
#endif
    
    const smearing_map_t *smearing_map_p = &(runCorrection_itr->second.smearing_map);
    
#ifdef HggCatTry
  // questo se voglio provare a cercare prima le HggCat nel file 
  // e poi se non c'e' usare quelle EB/EE standard
  if(smearing_map.count(GetElectronCategory(isEBEle, R9Ele, true))==1)
    return smearing_map[GetElectronCategory(isEBEle, R9Ele, true)].first;
  else if(smearing_map.count(GetElectronCategory(isEBEle, R9Ele, false))==1)
    return smearing_map[GetElectronCategory(isEBEle, R9Ele, false)].first;
  else {
    std::cerr << "[ERROR] Electron category not found!!!!  Scale offset not applied" << std::endl;
    return 1;
  }
#endif
  TString category=GetElectronCategory(isEBEle, R9Ele, etaSCEle);
#ifdef DDEBUG
  std::cout << "[DEBUG] category(isEBEle, R9Ele) = " << category << "(" << isEBEle << ", " << R9Ele << ")" <<   smearing_map_p->find(category)->second.first << std::endl;

#endif

  if(smearing_map_p->count(category)==1)
    return smearing_map_p->find(category)->second.first;
  else{
    std::cerr << "[ERROR] Electron category not found!!!!  Scale offset not applied" << std::endl;
    return 1;
  }

  
}

void EnergyScaleCorrection::Add(TString category_, int runMin_, int runMax_, double deltaP_, double err_deltaP_, double sigmaE_, double err_sigmaE_){

  // se le corrections non sono definite per Hgg (gold e bad) allora ignoro le righe che si riferiscono a queste categorie
  if( (!isHggCat) && (category_.Contains("gold") || category_.Contains("bad"))) {
#ifdef DEBUG
    std::cout << "[DEBUG] line of HggCat recalib: " << category_ << " ignore" << std::endl;
#endif
    return;
  }
  //  if( isHggCat && (!(category_.Contains("gold") || category_.Contains("bad")))){
#ifdef DEBUG
  //    std::cout << "[DEBUG] line of not HggCat recalib: " << category_ << " ignore" << std::endl;
#endif
  //   return;
  // }

  std::map <int, correction_t >::iterator itr=runMin_map.find(runMin_);
  if (itr!=runMin_map.end()){ // if exists
    if(itr->second.runMax != runMax_){
      std::cerr << "[ERROR] two run ranges with same runMin and different runMax" << std::endl;
      std::cerr << "        value not added" << std::endl;
      std::cerr << "        " << category_ << "\t" << runMin_ << "\t" << runMax_ << "\t" << deltaP_ << std::endl;
      std::cerr << "        " << itr->first << "\t" << itr->second.runMax << std::endl;
      return;
    }

    // probabilmente la categoria non e' ancora stata definita e la aggiungo
    //1-deltaP_; // deltaP_ e' lo shift, devo correggerlo (quindi -deltaP_)
    itr->second.correction_map[category_].first=deltaP_; // allo stato attuale passo come input il file con le correzioni (1-deltaP/100)
    itr->second.correction_map[category_].second=err_deltaP_; // inefficiente
    itr->second.smearing_map[category_].first=sigmaE_; // allo stato attuale passo come input il file con le correzioni (1-deltaP/100)
    itr->second.smearing_map[category_].second=err_sigmaE_; // inefficiente
    //#ifdef DEBUG
    std::cout << "[DEBUG] " << category_ << "\t" 
	      << itr->first << "-" << itr->second.runMax 
	      << "\t" << itr->second.correction_map.size() 
	      << "\t" << itr->second.correction_map[category_].first
	      << "\t" << itr->second.smearing_map.size() 
	      << "\t" << itr->second.smearing_map[category_].first << std::endl;
    //"\t" << deltaP_ << "\t" << itr->second.correction_map[category_].first << std::endl;
    //std::cout << "[DEBUG] " << err_deltaP_ << "\t" << itr->second.correction_map[category_].second << std::endl;
    //#endif
    
    
  } else {
    correction_t corr;
    corr.runMax=runMax_;
    corr.correction_map[category_]=std::make_pair<double,double>(deltaP_, err_deltaP_); //(1-deltaP_, err_deltaP_);
    corr.smearing_map[category_]=std::make_pair<double,double>(sigmaE_,err_sigmaE_); //(1-deltaP_, err_deltaP_);
    //    corr.correction_map[category_].second=err_deltaP_;
    runMin_map[runMin_]=corr; // aggiungo il nuovo
    itr=runMin_map.find(runMin_);
    
    std::cout << "[DEBUG] " << category_ << "\t" 
	      << itr->first << "-" << itr->second.runMax 
	      << "\t" << itr->second.correction_map.size() 
	      << "\t" << itr->second.correction_map[category_].first
	      << "\t" << itr->second.smearing_map.size() 
	      << "\t" << itr->second.smearing_map[category_].first << std::endl;

  }
  return;
}

void EnergyScaleCorrection::ReadFromFile(TString filename){
  std::cout << "[STATUS] Reading recalibration values from file: " << filename << std::endl;
  ifstream f_in(filename);
  if(!f_in.good()){
    std::cerr << "[ERROR] file " << filename << " not readable" << std::endl;
    return;
  }
  
  int runMin, runMax;
#ifdef SHERVIN
  double deltaM_data, err_deltaM_data;
  double deltaM_MC, err_deltaM_MC;
#endif
  TString category, region2;
  double deltaP, err_deltaP;
  double sigmaE, err_sigmaE;


  for(f_in >> category; f_in.good(); f_in >> category){
    f_in >> region2 
	 >> runMin >> runMax 
#ifdef SHERVIN
      	 >> deltaM_data >> err_deltaM_data 
      	 >> deltaM_MC >> err_deltaM_MC;
    deltaP=(deltaM_data - deltaM_MC)/91.188;
    err_deltaP=sqrt(err_deltaM_data*err_deltaM_data + err_deltaM_MC * err_deltaM_MC)/91.188;

#else
    >> deltaP >> err_deltaP >> sigmaE >> err_sigmaE;
#endif
    
#ifdef DEBUG
    std::cout << "[DEBUG]" << "\t" << "category" << "\t" << "runMin" << "\t" <<"deltaM_data" << "\t" << "deltaM_MC" << "\t" << "deltaP" << std::endl;
#ifdef SHERVIN
    std::cout << "[DEBUG]" << "\t" << category << "\t" << runMin << "\t" << runMax << "\t" << deltaM_data << "\t" << deltaM_MC << "\t" << deltaP << std::endl;
#else
    std::cout << "[DEBUG]" << "\t" << category << "\t" << runMin << "\t" << runMax << "\t" << "\t" << deltaP << "\t" << err_deltaP<< "\t" << sigmaE << "\t" << err_sigmaE<< std::endl;
#endif
#endif
    Add(category, runMin, runMax, deltaP, err_deltaP,sigmaE,err_sigmaE);
  }
  
  f_in.close();
  runCorrection_itr=runMin_map.begin();
#ifdef DEBUG
  std::cout << "[STATUS DEBUG] runCorrection_itr->first: " << runCorrection_itr->first << "\t runMin_map.begin(): " << runMin_map.begin()->first << std::endl;
#endif

  return;
}

