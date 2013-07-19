#define testweight_cxx
#include "testweight.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <iostream>

using std::cout;
using std::endl;


void testweight::Loop(bool cs)
{
//   In a ROOT session, you can do:
//      Root > .L testweight.C
//      Root > testweight t
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
   if (fChain == 0) return;

   TFile * hOutputFile   = new TFile("spectrumrew.root", "RECREATE" ) ;
   TH1D* pt_rew =  new TH1D("pt_rew","pt_rew",50,20,100);
   TH1D* pt_norew =   new TH1D("pt_norew","pt_norew",50,20,100);
   TH1D* mass_rew =  new TH1D("mass_rew","mass_rew",40,100,180);
   TH1D* mass_norew =  new TH1D("mass_norew","mass_norew",40,100,180);
   TH2D* pt2d_norew = new TH2D("pt2d_norew","pt2d_norew",15,25,100,15,55,160);
   TH2D* pt2d_rew = new TH2D("pt2d_rew","pt2d_rew",15,25,100,15,55,160);

   getweights();

   Long64_t nentries = fChain->GetEntriesFast();

   Long64_t nbytes = 0, nb = 0;
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      // if (Cut(ientry) < 0) continue;
      
      if(massggnewvtx<90 || massggnewvtx>190) continue;
      //if(massggnewvtx<100 || massggnewvtx>180) continue;
      
      if((TMath::Abs(etascphot1)>1.4442&&TMath::Abs(etascphot1)<1.566)||(TMath::Abs(etascphot2)>1.4442&&TMath::Abs(etascphot2)<1.566)
	 || TMath::Abs(etascphot1)>2.5 || TMath::Abs(etascphot2)>2.5) continue;  // acceptance
      
      if(ptphot1<55) continue; //pt first photon
      if(ptphot2<25) continue; //pt second photon
      
      //      if(ptcorrjet1<30 || TMath::Abs(etajet1)>4.7) continue; //pt first jet

      bool idphot1 = (idcicphot1 >= 4);
      bool idphot2 = (idcicphot2 >= 4);

      if(!cs){ // photon id no control sample

	if(!(idphot1)) continue;
	if(!(idphot2)) continue;
  
      }else{ // photon id for control sample
	
	if( !( (idphot1 && !idphot2 && !pid_hasMatchedPromptElephot2) || (idphot2 && !idphot1 && !pid_hasMatchedPromptElephot1) ) ) continue; 

      }

      pt_norew->Fill(ptphot2);
      mass_norew->Fill(massggnewvtx);

      double minptsublead(25), maxptsublead(100);
      double minptlead(55), maxptlead(160);
      double sizex = (maxptsublead - minptsublead)/15.;
      double sizey = (maxptlead - minptlead)/15.;
      int i = int((ptphot2-25)/sizex);
      int j = int((ptphot1-55)/sizey);
      if(i<0) i=0;
      if(j>14) j=14;
      if(i<0) i=0;
      if(j>14) j=14;
      if(i>-1 && i<15 && j>-1 && j<15){
	pt_rew->Fill(ptphot2,weights_[i][j]);
	mass_rew->Fill(massggnewvtx,weights_[i][j]);
	pt2d_norew->Fill(ptphot2,ptphot1,1);
	pt2d_rew->Fill(ptphot2,ptphot1,weights_[i][j]);
      }
      //      else  pt_rew->Fill(ptphot2,1.);
   }
   
   hOutputFile->Write() ;
   hOutputFile->Close() ;
   hOutputFile->Delete();
     
}




void testweight::getweights()
{

  TFile *f  = new TFile("ptreweight.root","READ");

  TH1D *puweights = 0;
  
  puweights= (TH1D*) f->Get("pt2d");

  for (int i = 0; i<15; i++) {
    for (int j = 0; j<15; j++) {
      float weight=1.;
      weight=puweights->GetBinContent(i+1,j+1);
      weights_[i][j] =  weight;
      cout << i << "  " << "  " << j << "   " << weight << endl; 
    }
  }
  
  //std::cout << "weights sum is " << sumPuWeights << std::endl;

}
