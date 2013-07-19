{
  gROOT->ProcessLine(".L makereweighthisto.C++");
  
  TString redntpDir= "root://pccmsrm23.cern.ch:1094//u2/xrootd/delre/Higgs/reduced/";
  TFile* data = TFile::Open(redntpDir+"redntp.42xv4_data_new.cicloose.regrPho_eCorr.v3/merged/redntp_Photon-Run2011A-03Oct2011-05JulReReco-05AugReReco-Prompt-v1-DiPhotonSkimOnFly-b.root");
  TFile* datacs = TFile::Open(redntpDir+"redntp.42xv4_data_new.preselectionCS.regrPho_eCorr.v3/merged/redntp_Photon-Run2011A-03Oct2011-05JulReReco-05AugReReco-Prompt-v1-DiPhotonSkimOnFly-b.root");
  //  TFile* datacs = TFile::Open(redntpDir+"redntp.42xv4_data_new.cicloose.regrPho_eCorr.v3/merged/redntp_Photon-Run2011A-03Oct2011-05JulReReco-05AugReReco-Prompt-v1-DiPhotonSkimOnFly-b.root");
  TFile * hOutputFile   = new TFile("ptreweight.root", "RECREATE" ) ;

  makereweighthisto data_fill((TTree*)data->Get("AnaTree"), 1);
  makereweighthisto datacs_fill((TTree*)datacs->Get("AnaTree"), 1);
  data_fill.Setcuts(55,25,-10000,-10000,-10000,-10000,0,0,-100000,0,-1,-1,0);
  data_fill.setCic(4);
  datacs_fill.Setcuts(55,25,-10000,-10000,-10000,-10000,0,0,-100000,0,-1,-1,0);
  datacs_fill.setCic(4);


  TH2D* pt2d_nocs =  data_fill.Plot("ptphot2","ptphot1","pt2d_nocs",15,25,100,20,40,160);
  TH2D* pt2d_cs =  datacs_fill.Plot("ptphot2","ptphot1","pt2d_cs",15,25,100,20,40,160,1);
  TH2D* pt2d = new TH2D("pt2d","pt2d",15,25,100,20,40,160);

  int norm1 = pt2d_nocs.Integral();
  int norm2 = pt2d_cs.Integral();
  pt2d_nocs.Scale(1./norm1);
  pt2d_cs.Scale(1./norm2);
  pt2d_nocs.Draw("box");
  pt2d_cs.Draw("samebox");
  c1.SaveAs("plot.png");

  for (int i=0; i<15; i++){
    for (int j=0; j<20; j++){
      double bin_nocs = pt2d_nocs.GetBinContent(i+1,j+1);
      double bin_cs = pt2d_cs.GetBinContent(i+1,j+1);
      double ratio(1);
      if(bin_nocs!=0 && bin_cs!=0) ratio = bin_nocs/bin_cs;
      if(ratio>100) ratio = 1;

      pt2d.SetBinContent(i+1,j+1,ratio);
    }
  }


  hOutputFile->Write() ;
  hOutputFile->Close() ;
  hOutputFile->Delete();

}
