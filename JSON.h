#ifndef JSON_H
#define JSON_H

#include <vector>
#include <map>

class JSON {
 public:
   JSON(const char* file);

   void ReadJSONFile(const char* json);
   bool isGoodLS(int run, int lumi);

 private:
   int oldRun;
   typedef std::vector< std::pair<int,int> > GoodLSVector;
   typedef std::map< int, GoodLSVector  >    LSRange ;
   LSRange goodLS_;
   LSRange::const_iterator goodLSCache_; // ptr to list of good LS for last run

};
#endif
