#define fillPlot_cxx
#include "fillPlot.h"
#include <TH2.h>
#include <TStyle.h>
#include <TRandom3.h>
#include <TCanvas.h>
#include <iostream>
#include <iomanip>
#include <fstream>

#define MAX_PU_REWEIGHT 22

TH1D * fillPlot::Plot(string var, string name, int nbin, double min, double max, bool cs)
{
//   In a ROOT session, you can do:
//      Root > .L fillPlot.C
//      Root > fillPlot t
//      Root > t.GetEntry(12); // Fill t data members with entry number 12
//      Root > t.Show();       // Show values of entry 12
//      Root > t.Show(16);     // Read and show values of entry 16
//      Root > t.Loop();       // Loop on all entries
//

//     This is the loop skeleton where:
//    jentry is the global entry number in the chain
//    ientry is the entry number in the current Tree
//  Note that the argument to GetEntry must be:
//    jentry for TChain::GetEntry
//    ientry for TTree::GetEntry and TBranch::GetEntry
//
//       To read only selected branches, Insert statements like:
// METHOD1:
//    fChain->SetBranchStatus("*",0);  // disable all branches
//    fChain->SetBranchStatus("branchname",1);  // activate branchname
// METHOD2: replace line
//    fChain->GetEntry(jentry);       //read all branches
//by  b_branchname->GetEntry(ientry); //read only this branch

   Long64_t nentries = fChain->GetEntriesFast();

   Long64_t nbytes = 0, nb = 0;
   TH1D * tempplot = new TH1D(name.c_str(),name.c_str(),nbin,min,max);
   
   ofstream outfile;

   TFile* fOut=0;
   TTree* myTree=0;

   if (var == "massgg" && writeRoot != "")
     {
       string filename(writeRoot);
       if (cs)
	   filename+=".cs";

       fOut=TFile::Open(filename.c_str(),"RECREATE"); 
       fOut->cd();
       myTree = new TTree("diPhotonEvents","");
       TString treeVariables = "run/I:lumi/I:event/I:massgg/F";    
       myTree->Branch("diPhotonEvents",&(tree_.run),treeVariables);
     }

   if (writetxt != "") 
     {
       string filename(writetxt);
       if (cs)
	   filename+=".cs";
       outfile.open(filename.c_str()); 
     }

   for (Long64_t jentry=0; jentry<nentries;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;

      // analysis cuts

      if((TMath::Abs(etascphot1)>1.4442&&TMath::Abs(etascphot1)<1.566)||(TMath::Abs(etascphot2)>1.4442&&TMath::Abs(etascphot2)<1.566)
	 || TMath::Abs(etascphot1)>2.5 || TMath::Abs(etascphot2)>2.5) continue;  // acceptance

      if(ptphot1<ptphot1cut) continue; //pt first photon
      if(ptphot2<ptphot2cut) continue; //pt second photon

      if(pthiggsmincut>0 && ptgg<pthiggsmincut) continue; //pt higgs min
      if(pthiggsmaxcut>0 && ptgg>=pthiggsmaxcut) continue; //pt higgs max

      if(ptjet1cut>0 && ptcorrjet1<ptjet1cut) continue; //pt first jet
      if(ptjet2cut>0 && ptcorrjet2<ptjet2cut) continue; //pt second jet

      //delteta
      if(deltaetacut!=0){
	if(deltaetacut>0){
	  if(TMath::Abs(deltaeta)<deltaetacut) continue;  // vbf selection 
	}else{
	  if(TMath::Abs(deltaeta)>-deltaetacut) continue;  // WZH selection
	}
      }

      //zeppenfeld
      if(zeppencut!=0) 
	if(TMath::Abs(zeppenjet)>zeppencut) continue; 

      //inv mass of jets
      if(invmassjetcut!=0){
	if(invmassjetcut>0){
	  if(invmassjet<invmassjetcut) continue; // vbf selection 
	}else{
	  if(TMath::Abs(invmassjet-85)>-invmassjetcut) continue; // WZH selection
	}
      }

      if(ebcat == 1) { // EB EE categories
	if((TMath::Abs(etascphot1)>1.4442||TMath::Abs(etascphot2)>1.4442)) continue; 
      } else if(ebcat == 0){
	if((TMath::Abs(etascphot1)<1.4442&&TMath::Abs(etascphot2)<1.4442)) continue; 
      }

      // r9 categories
      bool isr9phot1(0), isr9phot2(0);

      if(TMath::Abs(etascphot1)<1.4442 && r9phot1>.94) isr9phot1 = 1;
      if(TMath::Abs(etascphot2)<1.4442 && r9phot2>.94) isr9phot2 = 1;
      if(TMath::Abs(etascphot1)>1.4442 && r9phot1>.94) isr9phot1 = 1;
      if(TMath::Abs(etascphot2)>1.4442 && r9phot2>.94) isr9phot2 = 1;

      if(r9cat == 1) {
	if(!isr9phot1 || !isr9phot2) continue;
      } else if (r9cat == 0){
	if(isr9phot1 && isr9phot2) continue;
      } 

      // photon id
      bool idphot1(0), idphot2(0), looseidphot1(0), looseidphot2(0), pxlphot1(1), pxlphot2(1);

      if(pixelseedcut) { 
	pxlphot1 = !pid_haspixelseedphot1;
	pxlphot2 = !pid_haspixelseedphot2;
      }


      if(cicselection>0) {
	idphot1 = (idcicphot1 >= cicselection);
	idphot2 = (idcicphot2 >= cicselection);
      }else{	
	idphot1 = cutIDEG(ptphot1, etascphot1, pid_hlwTrackNoDzphot1, pid_jurECALphot1, pid_twrHCALphot1, pid_HoverEphot1, pid_etawidphot1, scaletrk, scaleecal, scalehcal, scalehove);
	idphot2 = cutIDEG(ptphot2, etascphot2, pid_hlwTrackNoDzphot2, pid_jurECALphot2, pid_twrHCALphot2, pid_HoverEphot2, pid_etawidphot2, scaletrk, scaleecal, scalehcal, scalehove);
      }


      if(!cs){ // photon id no control sample

	if(cicselection>0) {
	  if(!(idphot1)) continue;
	  if(!(idphot2)) continue;
	}else{
	  if(!(idphot1 && pxlphot1)) continue;
	  if(!(idphot2 && pxlphot2)) continue;
	}

      }else{ // photon id for control sample

	if(cicselection>0) {
	  looseidphot1 = (idcicphot1 > 0 );
	  looseidphot2 = (idcicphot2 > 0 );
	  //	  if( !( (idphot1 && looseidphot2 && !idphot2) || (idphot2 && looseidphot1 && !idphot1) ) ) continue; 
	  // Not perfect should be using the same electronVeto wrt CiC selection (now using matchedPromptEle veto)
	  if( !( (idphot1 && !idphot2 && !pid_hasMatchedPromptElephot2) || (idphot2 && !idphot1 && !pid_hasMatchedPromptElephot1) ) ) continue; 
	}else{
	  looseidphot1 = cutIDEG(ptphot1, etascphot1, pid_hlwTrackNoDzphot1, pid_jurECALphot1, pid_twrHCALphot1, pid_HoverEphot1, pid_etawidphot1, scaletrk*10, scaleecal*10, scalehcal*10, scalehove*10);
	  looseidphot2 = cutIDEG(ptphot2, etascphot2, pid_hlwTrackNoDzphot2, pid_jurECALphot2, pid_twrHCALphot2, pid_HoverEphot2, pid_etawidphot2, scaletrk*10, scaleecal*10, scalehcal*10, scalehove*10);
	  
	  if( !( (idphot1 && pxlphot1 && looseidphot2 && !idphot2) || (idphot2 && pxlphot2 && looseidphot1 && !idphot1) ) ) continue;
	}

      }

      // finding variable to be plotted
      double variable(0);
      if (var == "massgg")  variable = massgg;
      else if (var == "ptphot1")  variable = ptphot1;
      else if (var == "ptphot2")  variable = ptphot2;
      else if (var == "ptjet1")  variable = ptcorrjet1;
      else if (var == "ptjet2")  variable = ptcorrjet2;
      else if (var == "etajet1")  variable = etajet1;
      else if (var == "etajet2")  variable = etajet2;
      else if (var == "phijet1")  variable = phijet1;
      else if (var == "phijet2")  variable = phijet2;
      else if (var == "etaphot1")  variable = etaphot1;
      else if (var == "etaphot2")  variable = etaphot2;
      else if (var == "phiphot1")  variable = phiphot1;
      else if (var == "phiphot2")  variable = phiphot2;
      else if (var == "deltaeta")  variable = deltaeta;
      else if (var == "zeppenjet")  variable = zeppenjet;
      else if (var == "invmassjet")  variable = invmassjet;
      else if (var == "nvtx")  variable = nvtx;
      else if (var == "npu")  variable = npu;
      else if (var == "met")  variable = met;
      else if (var == "pid_haspixelseedphot1")  variable = pid_haspixelseedphot1;
      else if (var == "pid_jurECALphot1")  variable = pid_jurECALphot1;
      else if (var == "pid_twrHCALphot1")  variable = pid_twrHCALphot1;
      else if (var == "pid_HoverEphot1")  variable = pid_HoverEphot1;
      else if (var == "pid_hlwTrackNoDzphot1")  variable = pid_hlwTrackNoDzphot1;
      else if (var == "pid_etawidphot1")  variable = pid_etawidphot1;
      else if (var == "pid_haspixelseedphot2")  variable = pid_haspixelseedphot2;
      else if (var == "pid_jurECALphot2")  variable = pid_jurECALphot2;
      else if (var == "pid_twrHCALphot2")  variable = pid_twrHCALphot2;
      else if (var == "pid_HoverEphot2")  variable = pid_HoverEphot2;
      else if (var == "pid_hlwTrackNoDzphot2")  variable = pid_hlwTrackNoDzphot2;
      else if (var == "pid_etawidphot2")  variable = pid_etawidphot2;
      else{
	cout << "NO SUCH VARIABLE IMPLEMENTED!" << endl;
	break;
      }

      // energy smearing
      if(dosmear){
	TRandom3 smearing(565656);
	variable *= smearing.Gaus(meansmear,spreadsmear);
      }

      // pu reweighting
      if(npu<MAX_PU_REWEIGHT && puweights_.size()>0) 
	tempplot->Fill(variable, puweights_[npu]);
      else{
	//Using weight 1. for DATA and for dataOutside reweighting window
	tempplot->Fill(variable,1.);
	//	cout << "Event outside the reweighting range N<MAX_PU_REWEIGHT" << endl;
      } 

      if (var == "massgg" && writeRoot != "")
	{
	  tree_.run=run;
	  tree_.lumi=lumi;
	  tree_.event=event;
	  tree_.massgg=massgg;
	  fOut->cd();
	  myTree->Fill();
	}

      if(writetxt != "") 
	outfile << "run " << run << "\t lumi "  << std::setw(4) << lumi << "\t event " << std::setw(12) << event  << "\t massgg " << std::setprecision(6) << massgg << endl;      

   }
   
   if (var == "massgg" && writeRoot != "")
     fOut->Write();

   if(writetxt != "") 
     outfile.close();

   return tempplot;

}

void  fillPlot::Setcuts(double pt1, double pt2, double higgsptcutmin, double higgsptcutmax, double ptj1, double ptj2, double deltae, double zep, double mjj, int eb, int r9, int isolscaletrk, int isolscaleecal, int isolscalehcal, int isolscalehove, bool pixelseedveto)
{
  ptphot1cut = pt1;
  ptphot2cut = pt2;
  pthiggsmincut = higgsptcutmin;
  pthiggsmaxcut = higgsptcutmax;
  ptjet1cut = ptj1;
  ptjet2cut = ptj2;
  deltaetacut = deltae;
  zeppencut = zep;
  invmassjetcut = mjj;
  pixelseedcut = pixelseedveto;
  ebcat = eb;
  r9cat = r9;
  scaletrk = isolscaletrk;
  scaleecal = isolscaleecal;
  scalehcal = isolscalehcal;
  scalehove = isolscalehove;
  
}


void fillPlot::setCic(int cic) {
  cicselection = cic;
}


bool fillPlot::cutIDEG(double ptPhot, double etascphot, double pid_hlwTrackNoDz, double pid_jurECAL, double pid_twrHCAL, double pid_HoverE, double pid_etawid, int scatrk, int scaecal, int scahcal, int scahove) {

  double hovereiso_eb=           0.01725*double(scahove)/100.;
  double hovereiso_ee=           0.022  *double(scahove)/100.;
  double hcaliso_rel=            0.0025;
  double hcaliso_abs_eb=         1.475*double(scahcal)/100.;
  double hcaliso_abs_ee=         1.7  *double(scahcal)/100.;
  double ecaliso_rel=            0.006;
  double ecaliso_abs_eb=         2.84*double(scaecal)/100.;
  double ecaliso_abs_ee=         1.64*double(scaecal)/100.;
  double trackiso_rel=           0.001;
  double trackiso_abs_eb=        3.85*double(scatrk)/100.;
  double trackiso_abs_ee=        3.55*double(scatrk)/100.;
  double setaetaEB=              0.010;
  double setaetaEE=              0.028;

  bool ptiso; 
  bool ecaliso;
  bool hcaliso; 
  bool hoveiso; 
  bool setaeta; 

  if(TMath::Abs(etascphot) < 1.479) {
    ptiso = (pid_hlwTrackNoDz < ptPhot * trackiso_rel + 8.34071e-01 + 5.48136e-01*rhoPF - 1.5 + trackiso_abs_eb);
    ecaliso = (pid_jurECAL < ptPhot * ecaliso_rel + 1.58995 + 2.98677e-01*rhoPF - 2.0 + ecaliso_abs_eb );
    hcaliso = (pid_twrHCAL < ptPhot * hcaliso_rel + 1.49628 + 2.44899e-01*rhoPF - 2.0 + hcaliso_abs_eb );
    hoveiso = (pid_HoverE < 1.96440e-02 + 1.00859e-03*rhoPF - 0.02 + hovereiso_eb);
  }else{
    ptiso = (pid_hlwTrackNoDz < ptPhot * trackiso_rel + 8.86732e-01 + 5.25491e-01*rhoPF  - 1.5 + trackiso_abs_ee);
    ecaliso = (pid_jurECAL < ptPhot * ecaliso_rel + 8.32333e-01 + 1.91840e-01*rhoPF - 2.0 + ecaliso_abs_ee );
    hcaliso = (pid_twrHCAL < ptPhot * hcaliso_rel + 1.24901 + 2.74598e-01*rhoPF - 2.0 + hcaliso_abs_ee );
    hoveiso = (pid_HoverE < 1.95369e-02 + 1.14826e-03*rhoPF - 0.02 + hovereiso_ee);
  }

  setaeta = pid_etawid < setaetaEB;

  if(TMath::Abs(etascphot) > 1.479) {
    setaeta = pid_etawid < setaetaEE;
  }  


  return (ptiso && hcaliso && ecaliso && hoveiso && setaeta);
}

void fillPlot::Writetxt(char * filename)
{
  writetxt=std::string(filename);
}

void fillPlot::WriteRoot(char * filename)
{
  writeRoot=std::string(filename);
}

void fillPlot::SetPuWeights(bool isData,std::string puWeightFile)
{
  if (puWeightFile == "")
    {
      std::cout << "you need a weights file to use this function" << std::endl;
      return;
    }

  if (!isData)
    std::cout << "PU REWEIGHTING:: Using file " << puWeightFile << std::endl;

  TFile *f_pu  = new TFile(puWeightFile.c_str(),"READ");

  TH1D *puweights = 0;
  TH1D *gen_pu = 0;

  gen_pu= (TH1D*) f_pu->Get("generated_pu");
  puweights= (TH1D*) f_pu->Get("weights");

  if (!puweights || !gen_pu)
    {
      std::cout << "weights histograms  not found in file " << puWeightFile << std::endl;
      return;
    }
  
  
  TH1D* weightedPU= (TH1D*)gen_pu->Clone("weightedPU");
  weightedPU->Multiply(puweights);
  //Rescaling weights in order to preserve same integral of events
  TH1D* weights= (TH1D*)puweights->Clone("rescaledWeights");
  weights->Scale( gen_pu->Integral(1,MAX_PU_REWEIGHT) / weightedPU->Integral(1,MAX_PU_REWEIGHT) );
		  
  float sumPuWeights=0.;

  for (int i = 0; i<MAX_PU_REWEIGHT; i++) {
    float weight=1.;
    if( !isData ) 
	weight=weights->GetBinContent(i+1);
    
    sumPuWeights+=weight;
    puweights_.push_back(weight);
  }
  
  //  std::cout << "weights sum is " << sumPuWeights << std::endl;
}

void fillPlot::DoSmearing(double mean, double spread)
{
  dosmear = 1;
  meansmear = mean;
  spreadsmear = spread;
}
