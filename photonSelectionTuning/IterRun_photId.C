#include <vector>
#include <map>
#include <TROOT.h>
#include "SobIter_photId.cc"

int IterRun()
{ 
  // you have to include vector because the interpreter 
  // does not recognize vector by default 
  // if you have enabled FWLite or you have compiled STL in root
  // then the following line is not necessary 
  //  gROOT->LoadMacro("++"); 
  // define your filters 
  Filter signalfilter, bkgfilter; 
  signalfilter.NAME =  "signalFilter"; 
  signalfilter.ET1Min = 40; 
  signalfilter.ET2Min = 25; 
  signalfilter.mGGMin = 110;
  signalfilter.mGGMax = 120;

  bkgfilter.NAME =  "bkgFilter"; 
  bkgfilter.ET1Min = 40; 
  bkgfilter.ET2Min = 25; 
  bkgfilter.mGGMin = 100;
  bkgfilter.mGGMax = 150;

  float step = 0.01; 
  bool relative = true; 
  float min_step = 0.05; 
  float lowest_eff = 0.7; 

  // Set the weights 
  const float Nw=22; 
  float weights[22]; 
  weights[0]=0.1667; // type 0
  weights[1]=0.5352; // type 1
  weights[2]=0.8271; // type 2
  weights[3]=0.2422; // type 3
  weights[4]=0.3935; // type 4
  weights[5]=0.5670; // type 5
  weights[6]=0.07881; // type 6
  weights[7]=10.566; // type 7
  weights[8]=6.0976; // type 8
  weights[9]=3.7908; // type 9
  weights[10]=0.2482; // type 10
  weights[11]=0.04026; // type 11
  weights[12]=0.01367; // type 12
  weights[13]=0.003460; // type 13
  weights[14]=0.1667; // type 14
  weights[15]=0.1667; // type 15
  weights[16]=0.1667; // type 16
  weights[17]=0.1667; // type 17
  weights[18]=0.065330; // type 18
  weights[19]=0.1; // type 19
  weights[20]=0.1; // type 20
  weights[21]=0.1; // type 21

  // Set the files 
  vector<TString> files; 
  files.push_back("results/redntp_GluGluToHToGG_M-115-41x_ntpv1_00.root");
  files.push_back("results/redntp_GJet_Pt-20_doubleEMEnriched_TuneZ2_7TeV-pythia6-41x_ntpv1_00.root");
  files.push_back("results/redntp_GJet_Pt-20_doubleEMEnriched_TuneZ2_7TeV-pythia6-41x_ntpv1_01.root");
  files.push_back("results/redntp_GJet_Pt-20_doubleEMEnriched_TuneZ2_7TeV-pythia6-41x_ntpv1_02.root");

  // Set the variable related info 
  std::map<int, vector<TString> > branches;
  std::map<int, vector<TString> > names;
  std::map<int, vector<float> > initCut;
  std::map<int, vector<float> > minCut;
  std::map<int, vector<float> > minAllowed;
  std::map<int, vector<int> > points;

  for (int icat=0;icat<4;++icat)
    {
      vector<TString> branches_cat;
      vector<TString> names_cat;
      vector<float> initCut_cat;
      vector<float> minCut_cat;
      vector<float> minAllowed_cat;
      vector<int> points_cat;

      branches_cat.push_back("pid_hlwTrackNoDzphot2");
      names_cat.push_back("pid_hlwTrackNoDzphot2");
      initCut_cat.push_back(7.);
      minCut_cat.push_back(0.);
      minAllowed_cat.push_back(1.);
      points_cat.push_back(70);

      branches_cat.push_back("pid_jurECALphot2");
      names_cat.push_back("pid_jurECALphot2");
      initCut_cat.push_back(7.);
      minCut_cat.push_back(0.);
      minAllowed_cat.push_back(1.);
      points_cat.push_back(70);

      branches_cat.push_back("pid_twrHCALphot2");
      names_cat.push_back("pid_twrHCALphot2");
      initCut_cat.push_back(7.);
      minCut_cat.push_back(0.);
      minAllowed_cat.push_back(1.);
      points_cat.push_back(70);

      branches_cat.push_back("pid_HoverEphot2");
      names_cat.push_back("pid_HoverEphot2");
      initCut_cat.push_back(0.15);
      minCut_cat.push_back(0.);
      minAllowed_cat.push_back(0.01);
      points_cat.push_back(45);

      branches_cat.push_back("pid_etawidphot2");
      names_cat.push_back("pid_etawidphot2");
      if (icat<2) //Barrel categories
	{
	  initCut_cat.push_back(0.02);
	  minCut_cat.push_back(0.005);
	  minAllowed_cat.push_back(0.008);
	  points_cat.push_back(45);
	}
      else
	{
	  initCut_cat.push_back(0.04);
	  minCut_cat.push_back(0.01);
	  minAllowed_cat.push_back(0.025);
	  points_cat.push_back(50);
	}
      branches.insert(std::pair<int,vector<TString> >(icat,branches_cat));
      names.insert(std::pair<int,vector<TString> >(icat,names_cat));
      initCut.insert(std::pair<int,vector<float> >(icat,initCut_cat));
      minCut.insert(std::pair<int,vector<float> >(icat,minCut_cat));
      minAllowed.insert(std::pair<int,vector<float> >(icat,minAllowed_cat));
      points.insert(std::pair<int,vector<int> >(icat,points_cat));
    }

  SobIter iter(step, relative, lowest_eff, min_step,
               branches, names, initCut, minCut, points, minAllowed, files, signalfilter, bkgfilter); 
  iter.Analysis();
  cout << "Selection tuning has finished " << endl; 
  return 0;
} 
