void FindVariableBins(TH1F* h, int nbins, float* binValues)
{
  binValues[0]=h->GetBinLowEdge(1);
  binValues[nbins]=h->GetBinLowEdge(h->GetNbinsX())+h->GetBinWidth(h->GetNbinsX());
  std::cout << binValues[0] << " " << binValues[nbins] << " " <<   (float)h->Integral()
	    << std::endl;
  for(int ibin=0;ibin<nbins;++ibin)
    {
      
      std::cout << "+++++++++" << std::endl;
      std::cout << ((float)(ibin+1))/((float)nbins) << std::endl;
      float cumulativeThreshold=(((float)(ibin+1))/((float)nbins))*(float)h->Integral();
      std::cout << cumulativeThreshold << std::endl;
      for(int ib=1;ib<h->GetNbinsX();++ib)
	//	std::cout << h->Integral(1,ib+1) << std::endl;
	if (h->Integral(1,ib+1)>cumulativeThreshold)
	  {
	    std::cout << h->GetBinLowEdge(ib+1) << std::endl; 
	    binValues[ibin+1]=h->GetBinLowEdge(ib+1);
	    std::cout << binValues[ibin+1] << std::endl;
	    break;
	  }
    }
}
