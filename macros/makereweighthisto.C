#define makereweighthisto_cxx
#include "makereweighthisto.h"
#include <TH2.h>
#include <TStyle.h>
#include <TRandom3.h>
#include <TCanvas.h>
#include <iostream>
#include <iomanip>
#include <fstream>

#define MAX_PU_REWEIGHT 22

TH2D * makereweighthisto::Plot(string var1, string var2, string name, int nbin1, double min1, double max1, int nbin2, double min2, double max2, bool cs)
{
//   In a ROOT session, you can do:
//      Root > .L makereweighthisto.C
//      Root > makereweighthisto t
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
   TH2D * tempplot = new TH2D(name.c_str(),name.c_str(),nbin1,min1,max1,nbin2,min2,max2);
   
   ofstream outfile;

   TFile* fOut=0;
   TTree* myTree=0;

   if (var1 == "massgg" && writeRoot != "")
     {
       string filename(writeRoot);
       fOut=TFile::Open(filename.c_str(),"RECREATE"); 
       fOut->cd();
       myTree = new TTree("diPhotonEvents","");
       TString treeVariables = "run/I:lumi/I:event/I:ptgg/F:ebeb/I:massggnewvtx/F:weight/F";    
       myTree->Branch("diPhotonEvents",&(tree_.run),treeVariables);
     }

   if (writetxt != "") 
     {
       string filename(writetxt);
       outfile.open(filename.c_str()); 
     }

   for (Long64_t jentry=0; jentry<nentries;jentry++) {
      Long64_t ientry = LoadTree(jentry);

      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;

      // analysis cuts

      if(massggnewvtx<90 || massggnewvtx>190) continue;
      //if(massggnewvtx<100 || massggnewvtx>180) continue;

      if((TMath::Abs(etascphot1)>1.4442&&TMath::Abs(etascphot1)<1.566)||(TMath::Abs(etascphot2)>1.4442&&TMath::Abs(etascphot2)<1.566)
	 || TMath::Abs(etascphot1)>2.5 || TMath::Abs(etascphot2)>2.5) continue;  // acceptance

      if(ptphot1<ptphot1cut* massggnewvtx/120.) continue; //pt first photon
      //      if(ptphot1<ptphot1cut) continue; //pt first photon
      if(ptphot2<ptphot2cut) continue; //pt second photon

// TEMPPPP

//       if(ptphot1<ptphot1cut* massggnewvtx/100.) continue; //pt first photon
//       if(ptphot2<ptphot2cut* massggnewvtx/100.) continue; //pt second photon

      if(pthiggsmincut>0 && ptgg<pthiggsmincut) continue; //pt higgs min
      if(pthiggsmaxcut>0 && ptgg>=pthiggsmaxcut) continue; //pt higgs max


      if(ptjet1cut>0 && (ptcorrjet1<ptjet1cut || TMath::Abs(etajet1)>4.7)) continue; //pt first jet
      if(ptjet2cut>0 && (ptcorrjet2<ptjet2cut || TMath::Abs(etajet2)>4.7)) continue; //pt second jet

 
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

     //deltaphi
      if(deltaphicut!=0)
	if(TMath::Abs(deltaphi)<deltaphicut) continue;  // vbf selection 

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

//       if(pixelseedcut) { 
// 	pxlphot1 = !pid_haspixelseedphot1;
// 	pxlphot2 = !pid_haspixelseedphot2;
//       }
      
      idphot1 = (idcicphot1 >= cicselection);
      idphot2 = (idcicphot2 >= cicselection);

      if(!cs){ // photon id no control sample

	if(cicselection>0) {
	  if(!(idphot1)) continue;
	  if(!(idphot2)) continue;
	}else{
	  if(!(idphot1 && pxlphot1)) continue;
	  if(!(idphot2 && pxlphot2)) continue;
	}

      }else{ // photon id for control sample
	
	looseidphot1 = (idcicphot1 > 0 );
	looseidphot2 = (idcicphot2 > 0 );
	//	  if( !( (idphot1 && looseidphot2 && !idphot2) || (idphot2 && looseidphot1 && !idphot1) ) ) continue; 
	// Not perfect should be using the same electronVeto wrt CiC selection (now using matchedPromptEle veto)
	if( !( (idphot1 && !idphot2 && !pid_hasMatchedPromptElephot2) || (idphot2 && !idphot1 && !pid_hasMatchedPromptElephot1) ) ) continue; 

      }

      if(thirdcat && exclSel()) continue; 

      // finding variable1 to be plotted
      double variable1(0);
      if (var1== "massgg")  variable1 = massggnewvtx;
      else if (var1== "ptphot1")  variable1 = ptphot1;
      else if (var1== "ptphot2")  variable1 = ptphot2;
      else if (var1== "ptjet1")  variable1 = ptcorrjet1;
      else if (var1== "ptjet2")  variable1 = ptcorrjet2;
      else if (var1== "etajet1")  variable1 = TMath::Abs(etajet1);
      else if (var1== "etajet2")  variable1 = TMath::Abs(etajet2);
      else if (var1== "phijet1")  variable1 = phijet1;
      else if (var1== "phijet2")  variable1 = phijet2;
      else if (var1== "etaphot1")  variable1 = TMath::Abs(etaphot1);
      else if (var1== "etaphot2")  variable1 = TMath::Abs(etaphot2);
      else if (var1== "phiphot1")  variable1 = phiphot1;
      else if (var1== "phiphot2")  variable1 = phiphot2;
      else if (var1== "deltaeta")  variable1 = TMath::Abs(deltaeta);
      else if (var1== "deltaphi")  variable1 = TMath::Abs(deltaphi);
      else if (var1== "zeppenjet")  variable1 = TMath::Abs(zeppenjet);
      else if (var1== "invmassjet")  variable1 = invmassjet;
//       else if (var1== "deltaeta")  variable1 = deltaeta;
//       else if (var1== "zeppenjet")  variable1 = zeppenjet;
      else if (var1== "nvtx")  variable1 = nvtx;
      else if (var1== "npu")  variable1 = npu;
      else if (var1== "met")  variable1 = met;
      else if (var1== "ptgg")  variable1 = ptgg;
      else if (var1== "pid_haspixelseedphot1")  variable1 = pid_haspixelseedphot1;
      else if (var1== "pid_jurECALphot1")  variable1 = pid_jurECALphot1;
      else if (var1== "pid_twrHCALphot1")  variable1 = pid_twrHCALphot1;
      else if (var1== "pid_HoverEphot1")  variable1 = pid_HoverEphot1;
      else if (var1== "pid_hlwTrackNoDzphot1")  variable1 = pid_hlwTrackNoDzphot1;
      else if (var1== "pid_etawidphot1")  variable1 = pid_etawidphot1;
      else if (var1== "pid_haspixelseedphot2")  variable1 = pid_haspixelseedphot2;
      else if (var1== "pid_jurECALphot2")  variable1 = pid_jurECALphot2;
      else if (var1== "pid_twrHCALphot2")  variable1 = pid_twrHCALphot2;
      else if (var1== "pid_HoverEphot2")  variable1 = pid_HoverEphot2;
      else if (var1== "pid_hlwTrackNoDzphot2")  variable1 = pid_hlwTrackNoDzphot2;
      else if (var1== "pid_etawidphot2")  variable1 = pid_etawidphot2;
      else if (var1== "idcicphot1")  variable1 = idcicphot1;
      else if (var1== "idcicphot2")  variable1 = idcicphot2;      
      else{
	cout << "NO SUCH VARIABLE IMPLEMENTED!" << endl;
	break;
      }

      // finding variable2 to be plotted
      double variable2(0);
      if (var2== "massgg")  variable2 = massggnewvtx;
      else if (var2== "ptphot1")  variable2 = ptphot1;
      else if (var2== "ptphot2")  variable2 = ptphot2;
      else if (var2== "ptjet1")  variable2 = ptcorrjet1;
      else if (var2== "ptjet2")  variable2 = ptcorrjet2;
      else if (var2== "etajet1")  variable2 = TMath::Abs(etajet1);
      else if (var2== "etajet2")  variable2 = TMath::Abs(etajet2);
      else if (var2== "phijet1")  variable2 = phijet1;
      else if (var2== "phijet2")  variable2 = phijet2;
      else if (var2== "etaphot1")  variable2 = TMath::Abs(etaphot1);
      else if (var2== "etaphot2")  variable2 = TMath::Abs(etaphot2);
      else if (var2== "phiphot1")  variable2 = phiphot1;
      else if (var2== "phiphot2")  variable2 = phiphot2;
      else if (var2== "deltaeta")  variable2 = TMath::Abs(deltaeta);
      else if (var2== "deltaphi")  variable2 = TMath::Abs(deltaphi);
      else if (var2== "zeppenjet")  variable2 = TMath::Abs(zeppenjet);
      else if (var2== "invmassjet")  variable2 = invmassjet;
//       else if (var2== "deltaeta")  variable2 = deltaeta;
//       else if (var2== "zeppenjet")  variable2 = zeppenjet;
      else if (var2== "nvtx")  variable2 = nvtx;
      else if (var2== "npu")  variable2 = npu;
      else if (var2== "met")  variable2 = met;
      else if (var2== "ptgg")  variable2 = ptgg;
      else if (var2== "pid_haspixelseedphot1")  variable2 = pid_haspixelseedphot1;
      else if (var2== "pid_jurECALphot1")  variable2 = pid_jurECALphot1;
      else if (var2== "pid_twrHCALphot1")  variable2 = pid_twrHCALphot1;
      else if (var2== "pid_HoverEphot1")  variable2 = pid_HoverEphot1;
      else if (var2== "pid_hlwTrackNoDzphot1")  variable2 = pid_hlwTrackNoDzphot1;
      else if (var2== "pid_etawidphot1")  variable2 = pid_etawidphot1;
      else if (var2== "pid_haspixelseedphot2")  variable2 = pid_haspixelseedphot2;
      else if (var2== "pid_jurECALphot2")  variable2 = pid_jurECALphot2;
      else if (var2== "pid_twrHCALphot2")  variable2 = pid_twrHCALphot2;
      else if (var2== "pid_HoverEphot2")  variable2 = pid_HoverEphot2;
      else if (var2== "pid_hlwTrackNoDzphot2")  variable2 = pid_hlwTrackNoDzphot2;
      else if (var2== "pid_etawidphot2")  variable2 = pid_etawidphot2;
      else if (var2== "idcicphot1")  variable2 = idcicphot1;
      else if (var2== "idcicphot2")  variable2 = idcicphot2;      
      else{
	cout << "NO SUCH VARIABLE2 IMPLEMENTED!" << endl;
	break;
      }
//       // energy smearing
//       if(dosmear){
// 	TRandom3 smearing(565656);
// 	variable *= smearing.Gaus(meansmear,spreadsmear);
//       }

      // pu/pt reweighting
      float weight(1);
      if(dopureweight) weight *= pu_weight;
      if(doptreweight) weight *= pt_weight;

      int ebeb(0);
      if((TMath::Abs(etascphot1)<1.4442&&TMath::Abs(etascphot2)<1.4442)) ebeb = 1;	

      tempplot->Fill(variable1, variable2, weight);
      
      if (var1== "massgg" && writeRoot != "")
	{
	  tree_.run=run;
 	  tree_.lumi=lumi;
	  tree_.event=event;
 	  tree_.ptgg=ptgg;
	  tree_.ebeb=ebeb;
 	  tree_.massggnewvtx=massggnewvtx;
 	  tree_.weight=weight;
	  fOut->cd();
	  myTree->Fill();
	}

      if(writetxt != "") 
	outfile << "run " << run << "\t lumi "  << std::setw(4) << lumi << "\t event " << std::setw(12) << event  << "\t massgg " << std::setprecision(6) << massggnewvtx << endl;      

   }
   
   if (var1== "massgg" && writeRoot != ""){
     fOut->Write();
     delete fOut;
   }

   if(writetxt != "") 
     outfile.close();

   return tempplot;

}

void  makereweighthisto::Setcuts(double pt1, double pt2, double higgsptcutmin, double higgsptcutmax, double ptj1, double ptj2, double deltae, double zep, double mjj, double deltap, int eb, int r9, bool third)
{
  ptphot1cut = pt1;
  ptphot2cut = pt2;
  pthiggsmincut = higgsptcutmin;
  pthiggsmaxcut = higgsptcutmax;
  ptjet1cut = ptj1;
  ptjet2cut = ptj2;
  deltaetacut = deltae;
  deltaphicut = deltap;
  zeppencut = zep;
  invmassjetcut = mjj;
  ebcat = eb;
  r9cat = r9;
  thirdcat = third;
  
}

bool makereweighthisto::exclSel(){

  bool selvbf(1), selwzh(1);
  float ptleadvbf=55.;
  float ptsubleadvbf=25.;
  float ptjetleadvbf=30.;
  float ptjetsubleadvbf=20.;
  float deltaetavbf=3.5;
  float zeppenvbf=2.;
  float mjjvbf=400.;

  float ptleadwzh=65.;
  float ptsubleadwzh=25.;
  float ptjetleadwzh=35.;
  float ptjetsubleadwzh=20.;
  float deltaetawzh=-2.5;
  float zeppenwzh=1.5;
  float mjjwzh=-30.;

  if(ptphot1<ptleadvbf) selvbf =0; //pt first photon
  if(ptphot2<ptsubleadvbf) selvbf =0; //pt second photon
  if(ptcorrjet1<ptjetleadvbf) selvbf =0; //pt first jet
  if(ptcorrjet2<ptjetsubleadvbf) selvbf =0; //pt second jet
  if(TMath::Abs(deltaeta)<deltaetavbf) selvbf =0;  // vbf selection 
  if(TMath::Abs(zeppenjet)>zeppenvbf) selvbf =0; 
  if(invmassjet<mjjvbf) selvbf =0; // vbf selection 

  if(ptphot1<ptleadwzh) selwzh =0; //pt first photon
  if(ptphot2<ptsubleadwzh) selwzh =0; //pt second photon
  if(ptcorrjet1<ptjetleadwzh) selwzh =0; //pt first jet
  if(ptcorrjet2<ptjetsubleadwzh) selwzh =0; //pt second jet
  if(TMath::Abs(deltaeta)>deltaetawzh) selwzh =0;  // WZH selection
  if(TMath::Abs(zeppenjet)>zeppenwzh) selwzh =0; 
  if(TMath::Abs(invmassjet-85)>mjjwzh) selwzh =0; // WZH selection

  if(selvbf || selwzh) return 1;
  else return 0;

}

void makereweighthisto::setCic(int cic) {
  cicselection = cic;
}


bool makereweighthisto::cutIDEG(double ptPhot, double etascphot, double pid_hlwTrackNoDz, double pid_jurECAL, double pid_twrHCAL, double pid_HoverE, double pid_etawid, int scatrk, int scaecal, int scahcal, int scahove) {

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

void makereweighthisto::Writetxt(char * filename)
{
  writetxt=std::string(filename);
}

void makereweighthisto::WriteRoot(char * filename)
{
  writeRoot=std::string(filename);
}

void makereweighthisto::DoPuReweight(){
  dopureweight = 1;
}

void makereweighthisto::DoPtReweight(){
  doptreweight = 1;
}

void makereweighthisto::SetPuWeights(bool isData,std::string puWeightFile)
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

void makereweighthisto::DoSmearing(double mean, double spread)
{
  dosmear = 1;
  meansmear = mean;
  spreadsmear = spread;
}
