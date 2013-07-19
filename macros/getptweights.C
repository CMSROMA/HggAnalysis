{

  TFile *f  = new TFile("ptreweight.root","READ");

  double weights_[10][10];
  double minptsublead(20), maxptsublead(100);
  double minptlead(30), maxptlead(160);

  TH1D *puweights = 0;
  
  puweights= (TH1D*) f->Get("pt2d");
    
  for (int i = 0; i<10; i++) {
    for (int j = 0; j<10; j++) {
      float weight=1.;
      weight=puweights->GetBinContent(i+1,j+1);
      weights_[i][j] =  weight;
      cout << i << "  " << "  " << j << "   " << weight << endl; 
    }
  }
  
  //std::cout << "weights sum is " << sumPuWeights << std::endl;

}
