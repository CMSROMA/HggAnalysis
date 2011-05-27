#include "JSON.h"

#include <iostream>
#include <sstream>
#include <stdexcept>
using std::cout;
using std::endl;
using std::vector;
using std::map;
using std::pair;
using std::stringstream;
using std::runtime_error;


JSON::JSON(const char* json) {

  goodLS_ = LSRange();
  goodLSCache_ = goodLS_.end();

  cout << "Reading JSON file of good runs " << json << endl;
  FILE* iff = fopen(json,"r");

  if(iff == 0) {
    cout << "cannot open JSON file " << json << " ... now exiting." << endl;
    throw std::runtime_error("JSON file does not exist");
  }

  char c1, c2, c3;
  int run1, run2, LS1, LS2;

  cout << "Following LS will be used" << endl;
  cout << "-------------------------" << endl;
  while( fscanf(iff,"%c%d:%d-%d:%d%c%c,",&c1,&run1,&LS1,&run2,&LS2,&c2,&c3) != EOF ) {
      cout << "run: " << run1 << "  LS range: " << LS1
	   << " --> " << LS2 << "  " << c1 << "   " << c2<< "  " << c3<<endl;
      goodLS_[run1].push_back(  pair<int,int>(LS1,LS2) );
  }
  fclose(iff);
}

//========================
bool JSON::isGoodLS(int run, int lumi) {
//========================
     //if(!filterGoodRuns_) return true; // if filtered not requested all events are good

      // 
      if( oldRun != run ) {
        oldRun = run;
        goodLSCache_ = goodLS_.find( run );
      }

     // check whether this run is part of the good runs. else retrun false
     if( goodLSCache_ != goodLS_.end() ) {

        // get list of LS intervals
        const GoodLSVector& lsvector =   goodLSCache_->second; 
        // loop over good LS intervals and return as soon as one interval contains this event
        for(GoodLSVector::const_iterator iLS = lsvector.begin(); iLS != lsvector.end(); iLS++) {
           if(lumi >= iLS->first && lumi <= iLS->second ) {
             //cout << "Accepting run: " << Run << " LS: " << LumiSection << endl;
             return true;
           } // check current LS being in the interval
        } // loop over good LS for this run
     }
     return false;
}



