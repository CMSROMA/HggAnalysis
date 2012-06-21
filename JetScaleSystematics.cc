#include "JetScaleSystematics.h"
#include <algorithm>
#include <assert.h>
#include <fstream>
#include <iostream>

//----------------------------------------------------------------------
float getFloat(const std::string& token) 
{
  char* endptr;
  float result = strtod (token.c_str(), &endptr);
  if (endptr == token.c_str()) 
    {
      std::cout <<"can't convert token "<<token<<" to float value"<<std::endl;
    }
  return result;
} 

//----------------------------------------------------------------------
unsigned getUnsigned(const std::string& token) 
{
  char* endptr;
  unsigned result = strtoul (token.c_str(), &endptr, 0);
  if (endptr == token.c_str()) 
    {
      std::cout<<"can't convert token "<<token<<" to unsigned value"<<std::endl;
    }
  return result;
}

std::vector<std::string> getTokens(const std::string& fLine)
{
  std::vector<std::string> tokens;
  std::string currentToken;
  for (unsigned ipos = 0; ipos < fLine.length (); ++ipos) 

    {
      char c = fLine[ipos];
      //      cout << c <<endl;
      if (c == '#') break; // ignore comments
      else if (c == ' ') 
	{ // flush current token if any
	  if (!currentToken.empty()) 
	    {
	      tokens.push_back(currentToken);
	      currentToken.clear();
	    }
	}
      else
	currentToken += c;
      }
  if (!currentToken.empty()) tokens.push_back(currentToken); // flush end 
  return tokens;
}


JetScaleSystematics::JetScaleSystematics(TString inputfile)
{
  std::ifstream input(inputfile);
  std::string line;
  int counteta(-1);
  while (std::getline(input,line)) 
    {
      std::vector<std::string> tokens = getTokens(line); 
      if (!tokens.empty())
	{ 
	  counteta++;
	  if (tokens.size() < 6) 
	    {
	      std::cout<<"(line "<<inputfile<<"): less than 6 expected tokens:"<<tokens.size()<<std::endl;
	    }
	  
	  //	  cout << tokens.size() << "   " << tokens[0] << endl;
	  
	  for(int i=0; i<39; i++){
	    double value = getFloat(tokens[i*3+4]);
	    //	    cout << value << " " ;
	    unc[counteta][i] = value;
	    value = getFloat(tokens[i*3+3]);
	    minpt[counteta][i] = value;
	  }
	  //	  cout << endl;
	  mineta[counteta] = getFloat(tokens[0]);	 	  

	}
    }
  for (int j=0;j<38; j++){
    std::cout << mineta[j] << "  " ;
    for (int i=0;i<39; i++){
      std::cout << minpt[j][i] << "   " << unc[j][i] << "   " ;
    }
    std::cout << std::endl;
  }

  return;

}


  

float JetScaleSystematics::getJESUncertainty(double eta, double pt) 
{
  
  float uncertainty(0);

  for(int j=0; j<38; j++)

    if(eta > mineta[j] && eta < mineta[j+1]){

      for(int i=0; i<39; i++){

	if(pt > minpt[j][i] && pt<minpt[j][i+1]){
	  uncertainty = unc[j][i];
	  continue;
	}

      }

      continue;

    }

  //  std::cout << eta << "   " << pt << "   " << uncertainty << std::endl;

  return uncertainty;

}

