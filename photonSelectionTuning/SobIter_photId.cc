//
//
/*
  class SobIter:
  root implementation of the S/B Iterative method
  for selection tuning

  Version 2.0 Released 30 April 09

  Further Information/Inquiries:
  Nikolaos.Rompotis@Cern.ch

  Nikolaos Rompotis - 23 March 2009 - initial release
  26 March 09: Detector specific tuning implemented
  03 April 09: bug corrected in TChain implementation that was stopping the
               code when more than one files for bkg were specified
  30 April 09: version 2.0 This implementation gives identical results with
               the previous version, but has some more flexibility in the
               definition of signal and bkg. Instead of giving a signal and
               bkg files a signal filter and a bkg filter is defined in the
               input configuration and the electrons are read from a single
               dataset
  11 May 09:   version 2.1 enables to use an artificially lower cut in any
               variable that we want
  22 June 09:  version 3.1  bug  is corrected in sig/bkg mc values
               code extended to receive many different types as input
  22 Sept 09:  modification for the october exersize
  15 Marc 10:  met type chosen in cfg, preselection option in cfg
               some refinements in order to be able to calculate the S/B init
               old cfg still work, but the numbers may not have the same
               interpretation
               TO BE DONE: proper integration of Zee tuning, inverted cuts
*/

#include <iostream>
#include <vector>
#include <algorithm>
#include "TString.h"
#include "TTree.h"
#include "TBranch.h"
#include "TFile.h"
#include "TGraph.h"
#include "TMath.h"
#include "TH1F.h"
#include "TChain.h"
#include "TDatime.h"
#include "TCanvas.h"
#include "TROOT.h"


using namespace std;

struct Filter {
  // this structure allows the application of a filter to an event
  // the filter contains ET cut, MET cut and particular type of event cut
  // see function PassFilter in SobIter class
  float ET1Min;
  float ET2Min;

  float mGGMin;
  float mGGMax;

  TString NAME;
};

class SobIter
{
public:
  SobIter( float step, bool relative, float lowest_eff,
	   float min_step, 
	   std::map<int,vector<TString> > vars,
	   std::map<int,vector<TString> > names,
	   std::map<int, vector<float> > initcuts,
	   std::map<int, vector<float> > mins,
	   std::map<int, vector<int> > points,
	   std::map<int, vector<float> > lowerAllowed, 
	   vector<TString> files, 
	   Filter signalfilter, Filter bkgfilter
	  );

  typedef std::map<int,vector<TString> >  stringVectorMap;
  typedef std::map<int,vector<float> >  floatVectorMap;
  typedef std::map<int,vector<int> >  intVectorMap;

  typedef std::pair<int,vector<TString> >  stringVectorMapElement;
  typedef std::pair<int,vector<float> >  floatVectorMapElement;
  typedef std::pair<int,vector<int> >  intVectorMapElement;

  typedef std::map<int,vector<TString> >::const_iterator  stringVectorMapIter;
  typedef std::map<int,vector<float> >::const_iterator  floatVectorMapIter;
  typedef std::map<int,vector<int> >::const_iterator  intVectorMapIter;

  void Analysis();
private:
  bool RELATIVE_STEP_;
  float LOWEST_EFF_;
  std::map<int, vector<float> > InitCuts_;
  std::map<int, vector<float> > Mins_;
  std::map<int, vector<int> > Points_; 
  std::map<int, vector<TString> > Branches_;
  std::map<int, vector<TString> > Names_;
  std::map<int, vector<float> > LowerAllowed_;
  std::map<int, vector<float> > CurrentCuts_;
  std::map<int, int > nvars_;
  std::map<int, vector<int> > AllBranchesToVariable_;
  vector<float> LowerAllowedEndcaps_;
  vector<TString> AllBranches_;


  vector<TString> DATAFILES_;
  float minStep_;
//   float *w_;
//   int Nw_;
  float Step_;
  bool isAlive_;
  Filter SignalFilter_;
  Filter BkgFilter_;
  //

//   int n_barrel_;
//   int n_endcaps_;
  int n_variables_;
  int n_categories_;

  //
  bool CheckCuts( float *values, int category);
  int  GetCategory(float etaValue, float r9);
  bool PassFilter(Filter myfilter, float massgg, float et1, float et2);
  void  PrintFilter(Filter filter);
  int FindName(TString name, vector<TString> vec);
  vector<float>  ReturnSeffForSoB(float target, TH1F* h_signal, TH1F* h_bkg, float rest_signal, 
   float rest_bkg, float SIGNAL_INIT, float tolerance);
  float Bisection(TGraph *g, float a_, float b_, float target, float tolerance);
};


SobIter::SobIter( float step, bool relative, float lowest_eff,
		  float min_step,
		  std::map<int, vector<TString> > vars, 
		  std::map<int, vector<TString> > names, 
		  std::map<int, vector<float> > initcuts,  
		  std::map<int, vector<float> > mins, 
		  std::map<int, vector<int> > points,
		  std::map<int, vector<float> > lower,
		  vector<TString> files, 
		  Filter signalfilter, Filter bkgfilter)
{
  isAlive_ =true;
  // check for consistency
  if (vars.size() != names.size() ||
      vars.size() != initcuts.size() ||
      vars.size() != mins.size() ||
      vars.size() != points.size() ||
      vars.size() != lower.size()
      )
    {
      cout << "Inconsistent Number of Categories " << endl;
      isAlive_ = false;
      return;
    }
  
  
  for (stringVectorMapIter it=vars.begin();it!=vars.end();++it)
    if (it->second.size() != names[it->first].size())
	{
	  cout << "Inconsistent Number of Variables for category" << it->first << endl;
	  isAlive_ = false;
	  return;
	}

  for (stringVectorMapIter it=vars.begin();it!=vars.end();++it)
    if (it->second.size() != initcuts[it->first].size())
	{
	  cout << "Inconsistent Number of Variables for category" << it->first << endl;
	  isAlive_ = false;
	  return;
	}

  for (stringVectorMapIter it=vars.begin();it!=vars.end();++it)
    if (it->second.size() != mins[it->first].size())
	{
	  cout << "Inconsistent Number of Variables for category" << it->first << endl;
	  isAlive_ = false;
	  return;
	}

  for (stringVectorMapIter it=vars.begin();it!=vars.end();++it)
    if (it->second.size() != lower[it->first].size())
	{
	  cout << "Inconsistent Number of Variables for category" << it->first << endl;
	  isAlive_ = false;
	  return;
	}

  for (stringVectorMapIter it=vars.begin();it!=vars.end();++it)
    if (it->second.size() != points[it->first].size())
      {
	cout << "Inconsistent Number of Variables for category" << it->first << endl;
	isAlive_ = false;
	return;
      }
  
  //
  if ( (int) files.size() == 0){
    cout << "You must have at least one sample" << endl;
    isAlive_ = false;
    return;
  }
//   if (nw == 0) {
//     cout << "You should insert at least one weight" << endl;
//     isAlive_ = false;
//     return;
//   }
  for (stringVectorMapIter it=vars.begin();it!=vars.end();++it)
      if (it->second.size() == 0)
	{
	  cout << "You should have at least one variable" << endl;
	  isAlive_ = false;
	  return;
	}

  InitCuts_ = initcuts;
  Mins_ = mins;
  Points_ = points;
  Branches_ = vars;
  Names_ = names;
  LowerAllowed_ = lower;

  //  w_ = weights;
//   Nw_ = nw;
  DATAFILES_ = files;
  Step_ = step;
  RELATIVE_STEP_ = relative;
  LOWEST_EFF_ = lowest_eff;
  minStep_ = min_step;

  for (stringVectorMapIter it=vars.begin();it!=vars.end();++it)
    {
      nvars_.insert(std::make_pair(it->first,it->second.size()));
    }

  for (stringVectorMapIter it=vars.begin();it!=vars.end();++it)
    {
      std::vector<int> branchesToVariable;
      for (unsigned int i=0;i!=it->second.size();++i)
	{
	  if (std::find(AllBranches_.begin(),AllBranches_.end(),it->second[i])==AllBranches_.end())
	    {
	      AllBranches_.push_back(it->second[i]);
	      branchesToVariable.push_back(AllBranches_.size()-1);
	      cout << "didn't find " <<  it->second[i] << " assign ind: "
		   <<  AllBranches_.size()-1 << endl;
	    }
	  else
	    {
	      int ind=std::find(AllBranches_.begin(),AllBranches_.end(),it->second[i])-AllBranches_.begin();
	      branchesToVariable.push_back(ind);
	      cout << "Found " <<  it->second[i] << " and assign id: "
		   << ind << endl;
	    }
	}
      AllBranchesToVariable_.insert(std::make_pair(it->first,branchesToVariable));
    }

  n_variables_ = (int)  AllBranches_.size();
  n_categories_ = (int)  vars.size();
  SignalFilter_ = signalfilter;
  BkgFilter_ = bkgfilter;
  PrintFilter(SignalFilter_);
  PrintFilter(BkgFilter_);
  cout << "SobIter object was created..." << endl;
}

void SobIter::Analysis()
{
  if (not isAlive_) return;
  TString relative = "yes";
  if (not RELATIVE_STEP_) relative = "no";
  cout << "==============================" << endl;
  cout << "Iterative Method starts ...." << endl;
  cout << "Step = " << Step_ << ", relative? "<< relative
       << " Min step = " << minStep_       << endl;
  cout << "Variables and initial values: " << endl;

  //Printing values and initizialing the cuts values
  for (int icat=0;icat<n_categories_;++icat)
    {
      cout << "CATEGORY: " << icat << endl;
      std::vector<float> initCuts;
//       initCuts.reserve(InitCuts_[icat].size());
//       std::copy(InitCuts_[icat].begin(),InitCuts_[icat].end(),initCuts.begin());
      for (unsigned int i=0; i<InitCuts_[icat].size(); ++i) {
	cout << i << ". "<< Names_[icat][i] << ": "<<  InitCuts_[icat][i] << endl;
	initCuts.push_back(InitCuts_[icat][i]);
      }
      CurrentCuts_.insert(std::make_pair(icat,initCuts));
    }

//   cout << "WEIGHTS IN USE: " << endl;
//   for (int i=0; i<Nw_; ++i) cout << "w_[" << i << "]=" << w_[i] << endl;
  cout << "==============================" << endl;
  //
  //
  //
  //
  // for book keeping
  vector<float> SOB;   
  //vector<float> SOB_MC;
  vector<float> SEFF;  vector<float> SEFF_MC;
  vector<float> SOB_ERROR;  
  //vector<float> SOB_MC_ERROR;
  vector<TString> SELECTION;
  //
  // prepare the input files
  TChain  tree("AnaTree");
//   TChain  bkg("T1");
//   signal.Add(DATAFILES_[0]);
//   sA.Add(DATAFILES_[1]);
  vector<TString>::iterator it = DATAFILES_.begin();
  for ( ; it !=  DATAFILES_.end(); ++it) {
    tree.Add(*it); 
  }
  
  //
  // 
  const int N_signal = tree.GetEntries();
  cout << "Total contains " << N_signal << " events" << endl;
  //return;
  // / / / / / / / / / / / / / / / / / / / / / / / / / / / / / / / 
  /////////////////////////////////////////////////////////////////
  //                                         //////////////////////
  //   Sig/Bkg Branches   ...........................//////////////
  // ..................................................////////////
  // Common branch:
  float massgg;  // eta of the supercluster
  float et1;
  float et2;
  float r9phot1;
  float r9phot2;
  float etascphot1;
  float etascphot2;
  float rhoPF;
  float weight;
//   short int type;  // type of the event
//   short int exclusion; // ==1 when electron belongs to event with 2 elecs
  //
  tree.SetBranchAddress("massgg", &massgg);  
  tree.SetBranchAddress("ptphot1", &et1);  
  tree.SetBranchAddress("ptphot2", &et2);    
  tree.SetBranchAddress("r9phot1", &r9phot1);    
  tree.SetBranchAddress("r9phot2", &r9phot2);
  tree.SetBranchAddress("etascphot1", &etascphot1);    
  tree.SetBranchAddress("etascphot2", &etascphot2);    
  tree.SetBranchAddress("rhoPF", &rhoPF);    
  tree.SetBranchAddress("weight", &weight);  
  
  float *vars = new float[n_variables_];
  for (int i =0; i< n_variables_; ++i) { 
    tree.SetBranchAddress(AllBranches_[i], &vars[i]);
  }

  //////////////////////////////////////////////////////////////////////
  ///////////////////////////////////////////////
  ///////////////////////////////////////////
  ////////////////////////////////////
  ////////////////////////
  //
  //   ITERATION LOOP
  //   ^^^^^^^^^^^^^^
  //
  //

  float signal_init = 0;    // needed for the efficiency definition
  float signal_init_cat[n_categories_];
//   float signal_init_mc = 0;
//   float bkg_init_mc = 0;
  float bkg_init    = 0;
  float signal_init2 = 0, bkg_init2 = 0;
//   float signal_init_mc2 = 0, bkg_init_mc2 = 0;
  float sob_previous=-1;
  float sobInit = 0; 
    //    sobInitMC = 0;
  float sobInitError = 0;
  //sobInitMCError = 0;
  // .....................
  for (int Iter = 0; Iter<100000; ++Iter) {  // S/B Iteration loop
  //  for (int Iter = 0; Iter<1; ++Iter) {  // S/B Iteration loop
    // define your histograms .........................................

    vector <TH1F> h[n_categories_]; // signal
    vector <TH1F> b[n_categories_]; // bkg
    for (int icat=0; icat<n_categories_; ++icat) {
	char catName[50];
	sprintf(catName,"cat%d",icat);
	for (int i=0; i<nvars_[icat]; ++i) {
	  TString name = "h_";
	  name += Names_[icat][i];
	  name += "_";
	  name += TString(catName);
	  TH1F sig(name, name, Points_[icat][i], Mins_[icat][i], 
		 CurrentCuts_[icat][i]);
	  h[icat].push_back(sig);
	  TString b_name = "b_";
	  b_name += Names_[icat][i];
	  b_name += "_";
	  b_name += TString(catName);
	  TH1F bkg(b_name, b_name,Points_[icat][i], Mins_[icat][i], 
		   CurrentCuts_[icat][i]);
	  b[icat].push_back(bkg);
	}
    }

    float current_signal = 0;
    float current_bkg = 0;
    float current_signal_cat[n_categories_];    
    float current_bkg_cat[n_categories_];
    //    float current_actual_signal = 0; // assuming type 0 is the signal
    //    float current_actual_bkg = 0; // assuming type!= 0 is the bkg
    float current_signal_in_sigSample_cat[n_categories_];
    float current_signal_in_bkgSample_cat[n_categories_];
    
    for (int i=0;i<n_categories_;++i)
      {
	signal_init_cat[i]=0;
	current_signal_cat[i]=0.;
	current_bkg_cat[i]=0.;
	current_signal_in_sigSample_cat[i]=0.;
	current_signal_in_bkgSample_cat[i]=0.;
      }
  
    // statistical errors on s/b
    float dS2 = 0. , dB2 = 0.;
    //    float dS2_mc = 0. , dB2_mc = 0.;

    //
    // LOOP OVER THE TREE ENTRIES ........................................
    //
    // Loop over the signal 
    for (int iT = 0; iT < N_signal; ++iT) { // loop over the signal
      tree.GetEntry(iT);
      //Hardcoded absvalue
       for (int iv =0; iv< n_variables_; ++iv) 
	 {
	   if (TMath::Abs(etascphot2) < 1.4442)
	     {
	       if(AllBranches_[iv].Contains("pid_hlwTrackNoDz"))
		 vars[iv]=vars[iv]-(et2 * 0.001 + 5.48136e-01*rhoPF);
	       if(AllBranches_[iv].Contains("pid_jurECAL"))
		 vars[iv]=vars[iv]-(et2 * 0.006 + 2.98677e-01*rhoPF);
	       if(AllBranches_[iv].Contains("pid_twrHCAL"))
		 vars[iv]=vars[iv]-(et2 * 0.0025 + 2.44899e-01*rhoPF);
	       if(AllBranches_[iv].Contains("pid_HoverE"))
		 vars[iv]=vars[iv]-(1.00859e-03*rhoPF);
	     } 
	   else 
	     {
	       if(AllBranches_[iv].Contains("pid_hlwTrackNoDz"))
		 vars[iv]=vars[iv]-(et2 * 0.001 +  5.25491e-01*rhoPF);
	       if(AllBranches_[iv].Contains("pid_jurECAL"))
		 vars[iv]=vars[iv]-(et2 * 0.006 + 1.91840e-01*rhoPF);
	       if(AllBranches_[iv].Contains("pid_twrHCAL"))
		 vars[iv]=vars[iv]-(et2 * 0.0025 + 2.74598e-01*rhoPF);
	       if(AllBranches_[iv].Contains("pid_HoverE"))
		 vars[iv]=vars[iv]-(1.14826e-03*rhoPF);
	     }
	 }

      ///////// this is signal related
      bool isSignal=(TString(tree.GetCurrentFile()->GetName()).Contains("HToGG")); //nasty fix to avoid type as a branch  
      if (isSignal)
	weight= 300. * 18.23 * 2.13e-03 * 15 / 84058.; //Hardcoded weights for 300 pb-1
      else
	weight= 300. * 493.44 / 626717. / 5.; //Hardcoded weights for 300 pb-1. Dividing by 5 to normalize in the signal region (using 10 GeV)

      int icat=GetCategory(etascphot2,r9phot2);
      if (icat == -1)
	continue;
      bool pass = CheckCuts(vars, icat);

      //

      if (PassFilter(SignalFilter_, massgg, et1, et2) && isSignal) 
	{
	/////////////////////////////////////////////////////////////
	  //	  cout << "met=" << event_MET << ", et=" << probe_sc_et 
	  if (Iter == 0) { 
	    signal_init2 += (weight*weight);
	    signal_init_cat[icat] += weight;
	  }
	  //
	  if (pass) 
	    {
	      current_signal_cat[icat] += weight;
	      dS2 +=  weight*weight;
	      current_signal_in_sigSample_cat[icat] += weight;
	      for (int i=0; i<nvars_[icat]; ++i) {
		h[icat][i].Fill( vars[ AllBranchesToVariable_[icat][i] ], weight);
	      }
	    }
	  /////////////////////////////////////////////////////////////////////
	}
      //// bkg related ......................................................
      if (PassFilter(BkgFilter_, massgg, et1, et2) && !isSignal) {
	if (Iter == 0) {
	  bkg_init += weight; bkg_init2 += (weight*weight);
	}
	//
	if (pass) {
	  current_bkg_cat[icat] += weight;
	  dB2 +=  weight*weight;
	  current_signal_in_bkgSample_cat[icat] += weight;
	  for (int i=0; i<nvars_[icat]; ++i)
	    b[icat][i].Fill( vars[ AllBranchesToVariable_[icat][i] ] , weight);
	}
	
      }
      // MC values for signal and bkg -- to be changed -- this is for wenu
      // only
//       bool Ex =  true;
// //       if (SignalFilter_.EXCLUDE) Ex = exclusion != 1;
//       if (probe_sc_et> SignalFilter_.ET && Ex && somewhere) {
// 	if (isSignal) {
// 	  if (Iter==0) {
// 	    signal_init_mc += weight; signal_init_mc2+=(weight*weight);
// 	  }
// 	  if (pass) {
// 	    current_actual_signal += weight; dS2_mc += weight*weight;
// 	  }
// 	}
// 	//	if (type >=1  && type <=13) {
// 	else { // all the rest are bkg
// 	  if (Iter==0) {
// 	    bkg_init_mc += weight;  bkg_init_mc2 += (weight*weight);
// 	  }
// 	  if (pass) {
// 	    current_actual_bkg += weight; dB2_mc += weight*weight;
// 	  }
// 	}

    }
    //////////////////////////////////////////////////////////////////////
   // end loop over the samples
    // ...................................................................
    // -------------------------------------------------------------------
    // List of the counters that you have available now ------------------
    // ...................................................................
    // signal_init_barrel           signal sample contents at Iter == 0 
    // signal_init_endcaps, signal_init          
    // ...................................................................
    // signal_init_mc             number of events with type == 0 that
    //                            pass the kinematic restrictions of the
    //                            signal sample (ET, fiducial, exclusion)
    //                            at Iter == 0
    // bkg_init_mc                number of events with type != 0 + 
    //                            signal sample kinematic restrictions
    //                         >> may need to be modified for Zee-driven
    // ...................................................................
    // current_signal_barrel      number of events in signal sample that
    // current_signal_endcaps     pass kinematics+selection
    // current_bkg_barrel         number of events in bkg sample that pass
    // current_bkg_endcaps        kinematics+selection
    // ...................................................................
    // current_actual_signal      number of events that pass the kinematics
    // current_actual_bkg         + selection of the signal sample
    // ...................................................................
    // current_signal_in_sigSample  number of events with type == 0 in the
    // current_signal_in_bkgSample  signal and bkg sample
    // ...................................................................
    // to be done - slight redefinitions needed for the Zee-driven tuning
    // -------------------------------------------------------------------
    if (Iter==0) {
      
      TFile f("selectionVariables.root","RECREATE");
      f.cd();
      for (int icat=0;icat<n_categories_;icat++)
	for (int i=0;i<nvars_[icat];i++)
	  {
	    h[icat][i].Write();
	    b[icat][i].Write();
	  }
      f.Write();
      f.Close();

      for (int i=0;i<n_categories_;i++)
	signal_init += signal_init_cat[i];
      sobInit = signal_init/bkg_init;
      sobInitError = sobInit * sqrt(signal_init2/(signal_init*signal_init)+
				    bkg_init2/(bkg_init*bkg_init));
//       sobInitMC = signal_init_mc/bkg_init_mc;
//       sobInitMCError = sobInitMC * 
// 	sqrt(signal_init_mc2/(signal_init_mc*signal_init_mc)+
// 	     bkg_init_mc2/(bkg_init_mc*bkg_init_mc));
    }
    float sob_cat[n_categories_];
    float seff_cat[n_categories_]; 
    float current_signal_in_sigSample=0;
    float current_signal_in_bkgSample=0; 

    for (int i=0;i<n_categories_;i++)
      {
	current_signal += current_signal_cat[i];
	current_bkg += current_bkg_cat[i];
	current_signal_in_sigSample += current_signal_in_sigSample_cat[i];
	current_signal_in_bkgSample += current_signal_in_bkgSample_cat[i];
	sob_cat[i]= (current_bkg_cat[i]>0)?(current_signal_cat[i]/current_bkg_cat[i]):0;
	seff_cat[i] = (signal_init_cat[i]>0)?(current_signal_cat[i]/signal_init_cat[i]):0;
      }
    if (current_bkg == 0 || current_signal == 0) break;

    float sob = current_signal / current_bkg;
    float sob_error = sob * sqrt(dS2/(current_signal*current_signal)+
				  dB2/(current_bkg*current_bkg));

    if (Iter==0) {
      // this is in order to be able to have available the S/B after presel
      cout << "initial S/B= " << sobInit << " +- " << sobInitError
	//	   << " fromMC: " << sobInitMC << " +- " << sobInitMCError
	   << endl;
      cout << "Current signal = " << current_signal
	   << " Current bkg = " << current_bkg
	   << ", init=" << signal_init 
	//<< " sigMC=" << current_actual_signal
	// 	   << ", bkgMC=" << bkg_init_mc << ", sigInitMC="<<signal_init_mc
	   << endl;
    }
    // debugging...

    cout << "Current signal = " << current_signal 
	 << " Current bkg = " << current_bkg
	 << ", init=" << signal_init 
      //<< " sigMC=" << current_actual_signal
      //          << ", bkgMC=" << current_actual_bkg << ", sigInitMC="<<signal_init_mc
	 << ", CurrentSigInSig=" << current_signal_in_sigSample
	 << ", CurrentSigInBkg=" << current_signal_in_bkgSample
	 << endl;
    
    float seff = current_signal/signal_init; 

//     float sob_mc = 
//       (current_actual_bkg>0)?(current_actual_signal/current_actual_bkg):0;
//     float sob_mc_error = 
//       sob_mc * sqrt(dS2_mc/(current_actual_signal*current_actual_signal)+
// 		    dB2_mc/(current_actual_bkg*current_actual_bkg));
//     float seff_mc = 
//       (signal_init_mc>0)?(current_actual_signal/signal_init_mc):0;
    //
    TString selection = "(";
    for (int icat=0;icat<n_categories_;++icat)
      {
	for (int i=0; i<nvars_[icat]-1; ++i) {
	  selection +=  CurrentCuts_[icat][i];
	  selection += ", ";
	}
	if (nvars_[icat]>0)
	  selection +=  CurrentCuts_[icat][nvars_[icat]-1];
	else
	  {
	    selection +=  "(no class ";
	    selection +=  icat;
	    selection += " variables)";
	  }
	selection += " // ";
      }
    if (Iter == 0) {
      cout << "+++++++ Initial Iteration: S/B=" << sob << "  Initial Signal=" << signal_init;
      for (int icat=0;icat<n_categories_;++icat)
	cout  << "\nin Cat" << icat << " : " << signal_init_cat[icat] ;
      std::cout << "\nInitial Conditions: " << selection << endl;
    }
    else {
      cout << "+++++++ Iteration #" << Iter << " S/B=" << sob << "+/-" << sob_error
	   << "; seff=" << seff  << "; S="<< current_signal << "; B="
	   << current_bkg;
      for (int icat=0;icat<n_categories_;++icat)
	std::cout << "\nCat" << icat << " : S/B=" << sob_cat[icat] << "; seff="
		  << seff_cat[icat];
      std::cout << "\nCuts: " << selection << endl;
    }
    SOB.push_back(sob);    SEFF.push_back(seff);
    SOB_ERROR.push_back(sob_error);
    SELECTION.push_back(selection);
//     SOB_MC.push_back(sob_mc);    SEFF_MC.push_back(seff_mc);
//     SOB_MC_ERROR.push_back(sob_mc_error);
    // -----------------------------------------------------------------------
    // variable  tuning  to be performed here --------------------------------
    //
    // set the S/B target ....................................................
    float abs_step;
    if (RELATIVE_STEP_) {
      abs_step = Step_ * sob;
      ///if (Iter == 0) abs_step = 0.01;
      if (abs_step > minStep_) abs_step = minStep_;
    }
    else abs_step = Step_;
    float target = sob + abs_step;
    float tolerance = abs_step/100.;
    cout << "Targeting S/B " << target << ", S/B step is " << abs_step << endl;
    //
    // for each of the variables calculate the seff for the target S/B
    std::map<int,std::vector<float> > tmp_cuts;
    std::map<int,std::vector<float> >  tmp_seffs;
    std::map<int,std::vector<float> > tmp_sobs;
    float error_flag = 0;

    for (int icat=0; icat<n_categories_; ++icat) 
      {
	std::vector<float> cuts,seffs,sobs;
	for (int i=0; i<nvars_[icat]; ++i) {
	  cout << "Working on Cat" << icat << " " << Names_[icat][i] << endl;
	  vector<float> a = ReturnSeffForSoB
	    (target, &h[icat][i], &b[icat][i], 
	     current_signal - current_signal_cat[icat], current_bkg-current_bkg_cat[icat], signal_init, tolerance);
	  cuts.push_back(a[0]);    seffs.push_back(a[1]);    sobs.push_back(a[2]);
	  error_flag += a[3];
	  if (a[0]< LowerAllowed_[icat][i]) {
	    cout << "Cut Value Lower than Limit: cut will not be moved..."<< endl;
	    seffs[i] = 0.;
	  }
	}
	tmp_cuts.insert(std::make_pair(icat,cuts));
	tmp_seffs.insert(std::make_pair(icat,seffs));
	tmp_sobs.insert(std::make_pair(icat,sobs));
      }
    
    std::pair<int,int> bestCut;
    float maxeff=-2.;
    for (int icat=0; icat<n_categories_; ++icat) 
      for (int i=0; i<nvars_[icat]; ++i) 
	if (tmp_seffs[icat][i]>=maxeff)
	  {
	    bestCut=std::make_pair(icat,i);
	    maxeff=tmp_seffs[icat][i];
	  }
    
    if (maxeff < -1) {
      cout << "S/B cannot increased more ..." << endl;
      break;
    }
    // update the variable cuts now
    float movement=CurrentCuts_[bestCut.first][bestCut.second]-tmp_cuts[bestCut.first][bestCut.second];
    CurrentCuts_[bestCut.first][bestCut.second]=tmp_cuts[bestCut.first][bestCut.second];

    cout << "======= Iter #" << Iter << ": Variable Chosen: Cat" <<  bestCut.first
	 << " " <<  Names_[bestCut.first][bestCut.second] << " to " << tmp_cuts[bestCut.first][bestCut.second]
	 << " moved by " << movement << endl;
  
    
  
    //for (int i=0; i< n_endcaps_; ++i)
    //  cout << CurrentCutsEndcaps_[i] << endl;
    // release memory:
    // -----------------------------------------------------------------------
    //cout << "Current Signal is: " << current_signal << endl
    //	 << "Current Bkg is: " << current_bkg << endl;
    if (seff < LOWEST_EFF_)  break;
    if (movement == 0) break;
    if (sob_previous == sob || sob_previous > sob || error_flag != 0) break;
    sob_previous = sob;
  } // end S/B Iteration loop
  //
  // print out the summary of the method:
  const int Niter = (int) SOB.size();
  cout << "************** Tuning Summary ******************" << endl;
  cout << "Cuts appear in Category order:" << endl;
  for (int icat=0; icat<n_categories_; ++icat)
    {
      std::cout <<"\nCat" << icat;
      for (int i=0; i<nvars_[icat]-1; ++i)
	{
	  cout << Names_[icat][i] << ",";
	}
      if (nvars_[icat]>0)
	cout << Names_[icat][ nvars_[icat]-1 ] << " / ";
      else
	cout << "(no Cat " << icat << " variables)" << " / ";
    }
  cout << "\nfloat sob[" << Niter << "], seff[" << Niter <<  "]," << "sob_error["
       << Niter << "];"   << endl;
  for (int i=0;i<Niter; ++i)
    cout << "sob[" << i << "]=" << SOB[i] << "; seff[" << i << "]=" << SEFF[i]
	 << "; sob_error[" << i << "]=" << SOB_ERROR[i] 
	 << ";     //  "  << SELECTION[i] << endl;
  delete [] vars;
}


bool SobIter::CheckCuts( float *values, int category)
{
  bool res = true;
  for (int i=0; i< nvars_[category]; ++i) {
    bool cut = (values[ AllBranchesToVariable_[category][i] ]< 
		CurrentCuts_[category][i]);
    res = res && cut;
    if ( not cut) { 
      // 	cout << "failed #" << i << Names_endcaps_[i] << " val=" 
      // 	     << fabs(values[ BranchToVariable_[i+n_barrel_] ])
      // 	     << " current cut =" << CurrentCutsEndcaps_[i]
// 	     << endl;
      break;
    }
  }
  return res;
}

//
// try to find the string name in a vector of strings vec
// return the index of vec where you have found it, or -1
// if it is not there
int SobIter::FindName(TString name, vector<TString> vec)
{
  int res = -1;
  for (int i=0; i< (int) vec.size(); ++i)
    if (name == vec[i]) { res = i; break;}

  return res;
}


vector<float>  SobIter::ReturnSeffForSoB
(float target, TH1F* h_signal, TH1F* h_bkg,  
 float rest_signal, float rest_bkg, float SIGNAL_INIT, float tolerance)
{
   int points = h_signal->GetNbinsX();
  //cout << "DEBUG: S/B =" 
  //     << h_signal->Integral(0,nbins)/h_bkg->Integral(0,nbins) << endl;
  //
  float *sob = new float[points];
  float *seff= new float[points];
  float *cut = new float[points];
  float max_sob=0; float max_cut=0; float min_sob = 10000;
  int RealPoints = 0;
  for (int i=1; i < points+1; ++i) {
    float signal = h_signal->Integral(0,i) + rest_signal;
    float bkg = h_bkg->Integral(0,i)+rest_bkg;
    //cout << "signal " << signal << " signal_init " << signal_init << endl;
    seff[i-1] = signal/SIGNAL_INIT;
    if (bkg == 0) {
      cout << "ERROR!!!! Bkg is zero!!!" << endl;
      continue;
    }
    RealPoints++;
    sob[RealPoints-1] = signal/bkg;  //cout << signal/bkg << endl;
    cut[RealPoints-1] = h_signal->GetBinLowEdge(i) + h_signal->GetBinWidth(i);
    if (sob[RealPoints-1] > max_sob && sob[RealPoints-1] < 1000. ) { 
      // for more stable results
      max_sob = sob[RealPoints-1]; max_cut = cut[RealPoints-1];}
    if (sob[RealPoints-1] < min_sob) { min_sob = sob[RealPoints-1];}
  }
  TGraph g_sob(RealPoints, cut, sob); // this is to store the sob
  TGraph g_seff(RealPoints, cut, seff); // this is to store the sob

  //
  // bisection now to find the needed cut value
  int error_flag = 0;
  if (max_sob > 1000){
    cout << "ERROR: Your lower limit in this variable"
			  << " is too low: max sob > 1000!!!!!" << endl;
    error_flag = 1;
  }
  //
  if (sob[0]<sob[RealPoints-1]) {
    min_sob = sob[RealPoints-1];
  }
  //
  //
  // target too small
  if (min_sob > target) {
    cout << "ERROR:  min S/B is "<< min_sob <<" > target "<< target  << endl;
    cout << "Last Bin: Cut["<< RealPoints-1 <<"] = " << cut[RealPoints-1] 
	 << " Seff["<< RealPoints-1 <<"] = " << seff[RealPoints-1]
	 << " S/B["<< RealPoints-1  <<"] = " << sob[RealPoints-1]  << endl;
    error_flag = 2;
  }
  //
  //
  //
  // what happens if the target is too big
  if (max_sob < target) {
    if ( sob[0] < max_sob ) {
      cout << "NOTICE: S/B maximum "<< max_sob <<" is located at cut " 
	   << max_cut
	   << " and it is smaller than the target " << target << endl;
    }
    else {
      cout << "ERROR: S/B max is "<< max_sob  <<" < target " << target << endl;
      cout << "First Bin: Cut[0] = " << cut[0] << " Seff[0] = " 
	   << seff[0] << " S/B[0] = "  << sob[0]  << endl;
    }
    error_flag = 3;
  }

  //  float tolerance = target / 1000.; ==> inserted as an argument
  float working_point = 
    (error_flag == 2)?cut[RealPoints-1]:Bisection(&g_sob, max_cut, 
						  cut[RealPoints-1], 
						  target, tolerance);
  float working_seff = g_seff.Eval(working_point,0,"S");
  if (working_seff>1. && error_flag == 0) {
    cout << "SERIOUS ERROR: Fitting/Bisection gives efficiency " 
	 << working_seff	 << " manual replacement to 1" << endl;
    working_seff = 1.;
  }
  float working_sob = g_sob.Eval(working_point,0,"S");
  //
  // release memory:
  delete [] sob;   delete [] seff;   delete [] cut; 
  vector<float> res;
  res.push_back(working_point); // first entry
  if (error_flag == 2) {
    cout  << "Serious Error detected: program will exit..." << endl;
    res.push_back(-1.);
  }
  else if (error_flag == 3) {
    res.push_back(-3.);
    cout << "No active cut here " << endl;
  }
  else if (working_point > 999.) {
    cout << "Your actual values is seff = " << working_seff 
	 << " but this working point will be excluded because" 
	 <<  " the output is 1000 - Bisection Failure!!!!"
	 << endl;
    res.push_back(-1.);
  }
  else {
    res.push_back(working_seff); // sendond entry
    cout << "Working cut "<< working_point
	 << " Working seff is " << working_seff 
	 << " Working Sob is " << working_sob 
	 << endl;
  }
  res.push_back(working_sob); // 3rd entry
  if (error_flag == 2) {
    res.push_back(1);  // 4th entry
  } 
  else
    res.push_back(0);
  return res;


}
//
// this function is Bisection_Method2 from SelectionOptimizer class
// written by NR - December 08
//
// the meaning of tolerance  is wrt the target value
// i.e. if |f(c) - Target| < tolerance => stop
float SobIter::Bisection
(TGraph *g, float a_, float b_, float target, float tolerance)
{
  float a = a_;
  float b = b_;
  float c = (a+b)/2.;
  float fa = g->Eval(a,0,"S") - target;
  float fb = g->Eval(b,0,"S") - target;
  if (fa*fb < 0) {
    for (int i =0; i < 100000; ++i) {
      c = (a + b)/2.;
      float fc = g->Eval(c,0,"S") - target;
      if (fabs(fc) < tolerance) break;
      if (fc*fa < 0) {	fb = fc; b =c;      }
      else if (fc*fb < 0) {fa = fc; a=c;}
      else break;
    }

  }
  else if (fa*fb == 0) return (fa == 0.0 ? a:b);
  else return 10000;

  return c;
}

//
//
//
bool SobIter::PassFilter
(Filter myfilter, float massgg, float et1, float et2)
{
  if (myfilter.ET1Min>=0) 
    if (et1<myfilter.ET1Min) return false;

  if (myfilter.ET2Min>=0) 
    if (et2<myfilter.ET2Min) return false;

  if (myfilter.mGGMin>=0) 
    if (massgg<myfilter.mGGMin) return false;

  if (myfilter.mGGMax>=0) 
    if (massgg>myfilter.mGGMax) return false;

  return true;
}

int  SobIter::GetCategory(float etaValue, float r9)
{
  int category=-1;
  if (fabs(etaValue)<1.4442)
    {
      if (r9>0.93)
	category=0;
      else
	category=1;
    }
  else if (fabs(etaValue)>1.56 && fabs(etaValue)<2.5) 
    {
      if (r9>0.9)
	category=2;
      else
	category=3;
    }
  return category;
}

void SobIter::PrintFilter(Filter filter)
{
  
  cout << "Filter Summary of filter: " << filter.NAME << endl;
  
  if (filter.ET1Min >=0) 
    cout << "Cuts: et1>" << filter.ET1Min << endl;
  
  if (filter.ET2Min >=0) 
    cout << "Cuts: et2>" << filter.ET2Min << endl;
  
  if (filter.mGGMin >=0) 
    cout << "Cuts: mGG>" << filter.mGGMin << endl;
    
  if (filter.mGGMax >=0) 
    cout << "Cuts: mGG<" << filter.mGGMax << endl;
    
  cout << "End Summary for filter " << filter.NAME << endl;
}
