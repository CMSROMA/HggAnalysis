void handleError(const std::string& fClass, const std::string& fMessage)
  {
#ifdef STANDALONE 
    std::stringstream sserr;
    sserr<<fClass<<" ERROR: "<<fMessage;
    throw std::runtime_error(sserr.str());
#else
    //throw cms::Exception(fClass)<<fMessage;
#endif
  }

//----------------------------------------------------------------------
float getFloat(const std::string& token) 
{
  char* endptr;
  float result = strtod (token.c_str(), &endptr);
  if (endptr == token.c_str()) 
    {
      std::stringstream sserr; 
      sserr<<"can't convert token "<<token<<" to float value";
      handleError("getFloat",sserr.str());
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
      std::stringstream sserr; 
      sserr<<"can't convert token "<<token<<" to unsigned value";
      handleError("getUnsigned",sserr.str());
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

//---------------------------------------------------------------------- 

void readfile(const std::string& fLine)
{

  cout << fLine<< endl;
  std::ifstream input(fLine.c_str());
  std::string line;
  double mineta[28];  
  double minpt[28][39];
  double unc[28][39];
  int counteta(-1);
  while (std::getline(input,line)) 
    {
      std::vector<std::string> tokens = getTokens(line); 
      if (!tokens.empty())
	{ 
	  counteta++;
	  if (tokens.size() < 6) 
	    {
	      std::stringstream sserr;
	      sserr<<"(line "<<fLine<<"): less than 6 expected tokens:"<<tokens.size();
	      //handleError("JetCorrectorParameters::Definitions",sserr.str());
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
  for (int j=0;j<28; j++){
    cout << mineta[j] << "  " ;
    for (int i=0;i<39; i++){
      cout << minpt[j][i] << "   " << unc[j][i] << "   " ;
    }
    cout << endl;
  }
}
