#include "EnergyScaleCorrection.h"
#include <algorithm>
#include <assert.h>


EnergyScaleCorrection::EnergyScaleCorrection(const energyScaleParameters& par) : myParameters_(par)
{
  //Checking consistency of input parameters
  std::cerr << myParameters_.categoryType << " " <<  myParameters_.n_categories << std::endl;
  for(energyScaleParameters::eScaleVectorIt it=myParameters_.scale_offset_byrun.begin(); it!=myParameters_.scale_offset_byrun.end();
      ++it ) {
    assert( myParameters_.n_categories == (int) it->scale_offset.size() );
    assert( myParameters_.n_categories == (int) it->scale_offset_error.size() );
    assert( myParameters_.n_categories == (int) it->smearing.size() );
    assert( myParameters_.n_categories == (int) it->smearing_error.size() );
  }
}
  
std::string EnergyScaleCorrection::category(bool isEB, float r9) const
{
  std::string myCategory="";
  if (myParameters_.categoryType=="2CatR9_EBEE" )
    {
      if (isEB)
	myCategory+="EB";
      else
	myCategory+="EE";
      
      if (r9>=0.94)
	myCategory+="HighR9";
      else
	myCategory+="LowR9";
    }
  else if (myParameters_.categoryType=="EBEE")
    {
      if (isEB)
	myCategory+="EB";
	else
	  myCategory+="EE";
    }
  else
      {
	std::cout << "Unknown categorization. No category name is returned" << std::endl;
      }
  return myCategory;
}

float EnergyScaleCorrection::getScaleOffset(int run, const std::string & category) const
{
  const std::map<std::string, float> * scale_offset =  0;

  scale_offset = &(find(myParameters_.scale_offset_byrun.begin(),myParameters_.scale_offset_byrun.end(),run)->scale_offset) ;

  energyScaleParameters::parameterMapConstIt it=scale_offset->find(category);

  if ( it == scale_offset->end())
    {
      std::cout << "Category was not found in the configuration. Giving Up" << std::endl;
      return false;
    }

  return it->second;

}

float EnergyScaleCorrection::getSmearing(int run, const std::string & category) const
{
  const std::map<std::string, float> * smearing =  0;

  smearing = &(find(myParameters_.scale_offset_byrun.begin(),myParameters_.scale_offset_byrun.end(),run)->smearing) ;

  energyScaleParameters::parameterMapConstIt it=smearing->find(category);

  if ( it == smearing->end())
    {
      std::cout << "Category was not found in the configuration. Giving Up" << std::endl;
      return false;
    }

  return it->second;

}
