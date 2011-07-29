#ifndef __ENERGYSCALECORRECTION__
#define __ENERGYSCALECORRECTION__

#include <string>
#include <vector>
#include <iostream>
#include <map>


// ------------------------------------------------------------------------------------
class EnergyScaleOffset {
 public:
  EnergyScaleOffset(int first, int  last) : firstrun(first), lastrun(last) {};
    bool operator == (int run) const { return run>=firstrun && ( lastrun<0 || run<=lastrun); };

    int firstrun,  lastrun;

    std::map<std::string,float> scale_offset;
    std::map<std::string,float> scale_offset_error;
    std::map<std::string,float> smearing;
    std::map<std::string,float> smearing_error;

};

// ------------------------------------------------------------------------------------
class EnergyScaleCorrection 
{
 public:

     struct energyScaleParameters
     {
       int n_categories;
       std::string categoryType;
       std::string parameterSetName;

       typedef std::vector<EnergyScaleOffset> eScaleVector;
       typedef std::vector<EnergyScaleOffset>::iterator eScaleVectorIt;
       typedef std::vector<EnergyScaleOffset>::const_iterator eScaleVectorConstIt;

       typedef std::map<std::string,float> parameterMap;
       typedef std::map<std::string,float>::iterator parameterMapIt;
       typedef std::map<std::string,float>::const_iterator parameterMapConstIt;

       eScaleVector scale_offset_byrun;

     };

     EnergyScaleCorrection( const energyScaleParameters& par );
     ~EnergyScaleCorrection();

     float getScaleOffset(int run, const std::string & category) const;
     float getSmearing(int run, const std::string & category) const;

     std::string category(bool isEB, float r9) const;

     energyScaleParameters  myParameters_;

};

#endif
