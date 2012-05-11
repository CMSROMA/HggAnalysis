#define fillPlot_cxx
#include "fillPlot.h"
#include <TH2.h>
#include <TVector3.h>
#include <TLorentzVector.h>
#include <TStyle.h>
#include <TRandom3.h>
#include <TCanvas.h>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <TMath.h>

#define MAX_PU_REWEIGHT 22

inline double delta_phi(double phi1, double phi2) {

  double dphi = TMath::Abs(phi1 - phi2);
  return (dphi <= TMath::Pi())? dphi : TMath::TwoPi() - dphi;
}

TH1D * fillPlot::Plot(string var, string name, int nbin, double min, double max, bool cs, int signal)
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

   if(cs) getweights();

   if (var == "massgg" && writeRoot != "")
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
      
      // for VH, select only the production mechanism I want to study                                                                          
      if (signal==23 && gen_custom_processId!=3101 && gen_custom_processId!=3201 && gen_custom_processId!=3301) continue;
      if (signal==24 && gen_custom_processId!=3401) continue;
      if (signal==25 && gen_custom_processId!=4101 && gen_custom_processId!=4201 && gen_custom_processId!=4301) continue;
      if (signal==26 && gen_custom_processId!=4401) continue;
      if (signal==27 && gen_custom_processId!=4501) continue;

      // analysis cuts

      if(npu>=30) continue;

      // if(massggnewvtx<90 || massggnewvtx>190) continue;
      if(massggnewvtx<100 || massggnewvtx>180) continue;

      if((TMath::Abs(etascphot1)>1.4442&&TMath::Abs(etascphot1)<1.566)||(TMath::Abs(etascphot2)>1.4442&&TMath::Abs(etascphot2)<1.566)
	 || TMath::Abs(etascphot1)>2.5 || TMath::Abs(etascphot2)>2.5) continue;  // acceptance

      //     if(ptphot1<ptphot1cut) continue; //pt first photon
      if(ptphot2<ptphot2cut) continue; //pt second photon

// TEMPPPP

      if(ptphot1<ptphot1cut* massggnewvtx/120.) continue; //pt first photon
//       if(ptphot2<ptphot2cut* massggnewvtx/120.) continue; //pt second photon

      if(pthiggsmincut>0 && ptgg<pthiggsmincut) continue; //pt higgs min
      if(pthiggsmaxcut>0 && ptgg>=pthiggsmaxcut) continue; //pt higgs max


      if(ptjet1cut>0 && (ptcorrjet1<ptjet1cut || TMath::Abs(etajet1)>4.7)) continue; //pt first jet
      if(ptjet2cut>0 && (ptcorrjet2<ptjet2cut || TMath::Abs(etajet2)>4.7)) continue; //pt second jet

      // smeared / shifted met
      if(signal==100){ // data                                                                                                         
	if(metcut>0 && eShiftedScaledMet<metcut) continue;
      }
      if(signal!=100){ // MC                                                                                                           
	if(metcut>0 && eSmearedShiftedMet<metcut) continue;
      }

      // angular variable
      TVector3 t3met, t3jet1, t3phot1, t3phot2;
      if(signal==100) t3met.SetPtEtaPhi(eShiftedScaledMet,0,phiShiftedScaledMet);      // data                                         
      if(signal!=100) t3met.SetPtEtaPhi(eSmearedShiftedMet,0,phiSmearedShiftedMet);    // MC                                           
      // t3met.SetPtEtaPhi(epfMet,0,phipfMet);           // unsmeared                                                                  
      t3jet1.SetPtEtaPhi(ptcorrjet1,etajet1,phijet1);
      t3phot1.SetPtEtaPhi(ptphot1,etaphot1,phiphot1);
      t3phot2.SetPtEtaPhi(ptphot2,etaphot2,phiphot2);
      TVector3 t3higgs(t3phot1.Px()+t3phot2.Px(),t3phot1.Py()+t3phot2.Py(),t3phot1.Pz()+t3phot2.Pz());
      float deltaPhiMetJet1 = 10000.;
      if(ptcorrjet1>15.0) deltaPhiMetJet1 = fabs(180./3.1415927 * t3jet1.DeltaPhi(t3met));
      float deltaPhiMetHiggs = fabs(180./3.1415927 * t3higgs.DeltaPhi(t3met));
      float deltaPhiMetPhot1 = fabs(180./3.1415927 * t3phot1.DeltaPhi(t3met));
      float deltaPhiMetPhot2 = fabs(180./3.1415927 * t3phot2.DeltaPhi(t3met));
      float deltaPhiPho1Pho2 = fabs(180./3.1415927 * t3phot1.DeltaPhi(t3phot2));
      if (deltaPhiMetJet1>0  && deltaPhiMetJet1<phijetmetcut)    continue;
      if (deltaPhiMetHiggs>0 && deltaPhiMetHiggs<phihiggsmetcut) continue;
      if (deltaPhiMetPhot1>0 && deltaPhiMetPhot1<phiphot1metcut) continue;
      if (deltaPhiMetPhot2>0 && deltaPhiMetPhot2<phiphot2metcut) continue;
      if (deltaPhiPho1Pho2>0 && deltaPhiPho1Pho2>phipho1pho2cut) continue;
 
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

      // for the lepton tag analysis
      if (leptontag && !muonTagSelection() && !electronTagSelection()) continue;                                                    
      // if (leptontag && !electronTagSelection()) continue;
      // if (leptontag && !muonTagSelection()) continue;                                                                               

      // to apply the lepton tag veto                                                                                          
      if (leptonveto && (muonTagSelection() || electronTagSelection())) continue;

      // if (signal==100) cout << "data: event = " << event << ", run = " << run << ", mgg = " << massggnewvtx << endl;

      // finding variable to be plotted
      double variable(0);
      if (var == "massgg")  variable = massggnewvtx;
      else if (var == "ptphot1")  variable = ptphot1;
      else if (var == "ptphot2")  variable = ptphot2;
      else if (var == "ptjet1")  variable = ptcorrjet1;
      else if (var == "ptjet2")  variable = ptcorrjet2;
      else if (var == "etajet1")  variable = TMath::Abs(etajet1);
      else if (var == "etajet2")  variable = TMath::Abs(etajet2);
      else if (var == "phijet1")  variable = phijet1;
      else if (var == "phijet2")  variable = phijet2;
      else if (var == "etaphot1")  variable = TMath::Abs(etaphot1);
      else if (var == "etaphot2")  variable = TMath::Abs(etaphot2);
      else if (var == "phiphot1")  variable = phiphot1;
      else if (var == "phiphot2")  variable = phiphot2;
      else if (var == "deltaeta")  variable = TMath::Abs(deltaeta);
      else if (var == "deltaphi")  variable = TMath::Abs(deltaphi);
      else if (var == "zeppenjet")  variable = TMath::Abs(zeppenjet);
      else if (var == "invmassjet")  variable = invmassjet;
//       else if (var == "deltaeta")  variable = deltaeta;
//       else if (var == "zeppenjet")  variable = zeppenjet;
      else if (var == "nvtx")  variable = nvtx;
      else if (var == "npu")  variable = npu;
      else if (var == "rhopfnvtx")  variable = (rhoPF-0.9*nvtx)/nvtx;
      else if (var == "met")  variable = epfMet;
      else if (var == "smearedMet" && signal==100) variable = eShiftedScaledMet;   // data                                             
      else if (var == "smearedMet" && signal!=100) variable = eSmearedShiftedMet;  // MC                                               
      else if (var == "smearedMetPhi" && signal==100) variable = phiShiftedScaledMet;   // data                                        
      else if (var == "smearedMetPhi" && signal!=100) variable = phiSmearedShiftedMet;  // MC       
      else if (var == "rhopf")  variable = rhoPF;
      else if (var == "ptgg")  variable = ptgg;
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
      else if (var == "idcicphot1")  variable = idcicphot1;
      else if (var == "idcicphot2")  variable = idcicphot2;      
      else if (var == "methiggs") variable = deltaPhiMetHiggs;
      else if (var == "dphigg")   variable = deltaPhiPho1Pho2;
      else{
	cout << "NO SUCH VARIABLE IMPLEMENTED!" << endl;
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

      double minptsublead(25), maxptsublead(100);
      double minptlead(40), maxptlead(160);
      double sizex = (maxptsublead - minptsublead)/15.;
      double sizey = (maxptlead - minptlead)/15.;
      int i = int((ptphot2-minptsublead)/sizex);
      int j = int((ptphot1-minptlead)/sizey);
      if(i<0) i=0;
      if(j>14) j=14;
      if(i<0) i=0;
      if(j>19) j=19;
 
      if(cs) 
	if(i>-1 && i<16 && j>-1 && j<20) weight *= weights_[i][j];
    
      int ebeb(0);
      if((TMath::Abs(etascphot1)<1.4442&&TMath::Abs(etascphot2)<1.4442)) ebeb = 1;	

      tempplot->Fill(variable, weight);
      
      if (var == "massgg" && writeRoot != "")
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
   
   if (var == "massgg" && writeRoot != ""){
     fOut->Write();
     delete fOut;
   }

   if(writetxt != "") 
     outfile.close();

   return tempplot;
}

void  fillPlot::Setcuts(double pt1, double pt2, double higgsptcutmin, double higgsptcutmax, double ptj1, double ptj2, double misset,double deltae, double zep, double mjj, double deltap, double jetmet, double p1met, double p2met, double hmet, double delphigg, int eb, int r9, bool third, bool leptag, bool lepveto) {

  ptphot1cut = pt1;
  ptphot2cut = pt2;
  pthiggsmincut = higgsptcutmin;
  pthiggsmaxcut = higgsptcutmax;
  ptjet1cut = ptj1;
  ptjet2cut = ptj2;
  metcut = misset;
  deltaetacut = deltae;
  deltaphicut = deltap;
  phijetmetcut   = jetmet;
  phiphot1metcut = p1met;
  phiphot2metcut = p2met;
  phihiggsmetcut = hmet;
  phipho1pho2cut = delphigg;
  zeppencut = zep;
  invmassjetcut = mjj;
  ebcat = eb;
  r9cat = r9;
  thirdcat = third;
  leptontag  = leptag;
  leptonveto = lepveto;
}


bool fillPlot::exclSel(){

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

void fillPlot::DoPuReweight(){
  dopureweight = 1;
}

void fillPlot::DoPtReweight(){
  doptreweight = 1;
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



void fillPlot::getweights()
{

  TFile *f  = new TFile("ptreweight.root","READ");

  TH1D *puweights = 0;
  
  puweights= (TH1D*) f->Get("pt2d");

 for (int i = 0; i<15; i++) {
    for (int j = 0; j<20; j++) {
      float weight=1.;
      weight=puweights->GetBinContent(i+1,j+1);
      weights_[i][j] =  weight;
      // cout << i << "  " << "  " << j << "   " << weight << endl; 
    }
  }
  
  //std::cout << "weights sum is " << sumPuWeights << std::endl;

}

bool fillPlot::muonTagSelection(){

  bool passMuTag = true;

  // at least one full reco/id/isol muon                                                                                             
  if (ptmu1<0) passMuTag = false;

  // non implementazione migliore perche' andrebbe applicata a tutti... ma vabbe'                                                    
  if (deltaRToTrackphot1<1) passMuTag = false;
  if (deltaRToTrackphot2<1) passMuTag = false;

  if (passMuTag) {

    // lepton/photon distance                                                                                                        
    TVector3 t3muon, t3phot1, t3phot2;
    t3muon.SetPtEtaPhi(ptmu1,etamu1,phimu1);
    t3phot1.SetPtEtaPhi(ptphot1,etaphot1,phiphot1);
    t3phot2.SetPtEtaPhi(ptphot2,etaphot2,phiphot2);
    float deltaR1 = t3muon.DeltaR(t3phot1);
    float deltaR2 = t3muon.DeltaR(t3phot2);
    if (deltaR1<1) passMuTag = false;
    if (deltaR2<1) passMuTag = false;
  }

  return passMuTag;
}

bool fillPlot::electronTagSelection(){

  bool passEleTag = true;

  // not considered if passed the muon tag                                                                                           
  if (muonTagSelection()) passEleTag = false;
  
  // at least one full reco/id/isol electron                                                                                         
  if (ptele1<0) passEleTag = false;

  // non implementazione migliore perche' andrebbe applicata a tutti... ma vabbe'                                                    
  if (deltaRToTrackphot1<1) passEleTag = false;
  if (deltaRToTrackphot2<1) passEleTag = false;

  if (passEleTag) {

    // lepton/photon distance                                                                                                        
    TVector3 t3ele1, t3phot1, t3phot2;
    t3ele1.SetPtEtaPhi(ptele1,etaele1,phiele1);
    t3phot1.SetPtEtaPhi(ptphot1,etaphot1,phiphot1);
    t3phot2.SetPtEtaPhi(ptphot2,etaphot2,phiphot2);
    float deltaR1 = t3ele1.DeltaR(t3phot1);
    float deltaR2 = t3ele1.DeltaR(t3phot2);
    if (deltaR1<1) passEleTag = false;
    if (deltaR2<1) passEleTag = false;

    // invariant mass                                                                                                                
    float enePhot1 = ptphot1/(fabs(sin(t3phot1.Theta())));
    float enePhot2 = ptphot2/(fabs(sin(t3phot2.Theta())));
    float eneEle   = ptele1/(fabs(sin(t3ele1.Theta())));
    TLorentzVector t4ele1, t4phot1, t4phot2;

    // t4ele1.SetPtEtaPhiE(ptele1,etaele1,phiele1,eneEle);
    // t4phot1.SetPtEtaPhiE(ptphot1,etaphot1,phiphot1,enePhot1);
    // t4phot2.SetPtEtaPhiE(ptphot2,etaphot2,phiphot2,enePhot2);

    t4ele1.SetPtEtaPhiM(ptele1,etaele1,phiele1,0.000511); 
    t4phot1.SetPtEtaPhiM(ptphot1,etaphot1,phiphot1,0.);
    t4phot2.SetPtEtaPhiM(ptphot2,etaphot2,phiphot2,0.);

    float invMass1 = (t4phot1 + t4ele1).M();
    float invMass2 = (t4phot2 + t4ele1).M();
    float diff1 = fabs(invMass1-91.19);
    float diff2 = fabs(invMass2-91.19);

    if ( diff1<5. ) passEleTag = false;
    if ( diff2<5. ) passEleTag = false;
  }

  return passEleTag;
}
