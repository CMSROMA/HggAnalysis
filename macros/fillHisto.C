#define fillHisto_cxx
#include "fillHisto.h"
#include <TH2.h>
#include <TStyle.h>
#include <TRandom3.h>
#include <TCanvas.h>
#include <iostream>
#include <iomanip>
#include <fstream>

#define MAX_PU_REWEIGHT 22

TFile * fillHisto::File(char* writeRoot, bool cs)
{
//   In a ROOT session, you can do:
//      Root > .L fillHisto.C
//      Root > fillHisto t
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
   
   ofstream outfile;

   if(cs) getweights();

   TFile* fOut=0;
   TTree* myTree=0;
   string filename(writeRoot);
   fOut=TFile::Open(filename.c_str(),"RECREATE"); 
   fOut->cd();
   TH1D * mass_ptgg[20];
   TH1D * mass_phot1[20];
   TH1D * mass_phot2[20];
   TH1D * mass_jet1[20];
   TH1D * mass_jet2[20];
   TH1D * mass_deltae[20];
   TH1D * mass_zep[20];
   TH1D * mass_mjj[20];
   TH1D * mass_deltap[20];
   TH1D * mass_met[20];

   char namehisto[100], titlehisto[100]; 

   for(int i=0; i<20; i++){

     sprintf(namehisto,"%s%i","ptgg_",i);sprintf(titlehisto,"%s%i","variation of ptgg ",i);
     mass_ptgg[i] = new TH1D(namehisto,titlehisto,80,100,180);
     sprintf(namehisto,"%s%i","phot1_",i);sprintf(titlehisto,"%s%i","variation of phot1 ",i);
     mass_phot1[i] = new TH1D(namehisto,titlehisto,80,100,180);
     sprintf(namehisto,"%s%i","phot2_",i);sprintf(titlehisto,"%s%i","variation of phot2 ",i);
     mass_phot2[i] = new TH1D(namehisto,titlehisto,80,100,180);
     sprintf(namehisto,"%s%i","jet1_",i);sprintf(titlehisto,"%s%i","variation of jet1 ",i);
     mass_jet1[i] = new TH1D(namehisto,titlehisto,80,100,180);
     sprintf(namehisto,"%s%i","jet2_",i);sprintf(titlehisto,"%s%i","variation of jet2 ",i);
     mass_jet2[i] = new TH1D(namehisto,titlehisto,80,100,180);
     sprintf(namehisto,"%s%i","deltae_",i);sprintf(titlehisto,"%s%i","variation of deltaeta ",i);
     mass_deltae[i] = new TH1D(namehisto,titlehisto,80,100,180);
     sprintf(namehisto,"%s%i","zep_",i);sprintf(titlehisto,"%s%i","variation of zeppenfeld ",i);
     mass_zep[i] = new TH1D(namehisto,titlehisto,80,100,180);
     sprintf(namehisto,"%s%i","mjj_",i);sprintf(titlehisto,"%s%i","variation of mjj ",i);
     mass_mjj[i] = new TH1D(namehisto,titlehisto,80,100,180);
     sprintf(namehisto,"%s%i","deltap_",i);sprintf(titlehisto,"%s%i","variation of deltap ",i);
     mass_deltap[i] = new TH1D(namehisto,titlehisto,80,100,180);
     sprintf(namehisto,"%s%i","met_",i);sprintf(titlehisto,"%s%i","variation of met ",i);
     mass_met[i] = new TH1D(namehisto,titlehisto,80,100,180);
     
   }
  
   TH1D * tempplot = new TH1D("massgg","massgg",80,100,180);
   myTree = new TTree("diPhotonEvents","");
   TString treeVariables = "run/I:lumi/I:event/I:ptgg/F:ebeb/I:massggnewvtx/F:weight/F";    
   myTree->Branch("diPhotonEvents",&(tree_.run),treeVariables);

   if (writetxt != "") 
     {
       string filename2(writetxt);
       outfile.open(filename2.c_str()); 
     }

   for (Long64_t jentry=0; jentry<nentries;jentry++) {
      Long64_t ientry = LoadTree(jentry);

      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;

      // analysis cuts

      if(npu>=30) continue;

      if(massggnewvtx<100 || massggnewvtx>180) continue;

      if((TMath::Abs(etascphot1)>1.4442&&TMath::Abs(etascphot1)<1.566)||(TMath::Abs(etascphot2)>1.4442&&TMath::Abs(etascphot2)<1.566)
	 || TMath::Abs(etascphot1)>2.5 || TMath::Abs(etascphot2)>2.5) continue;  // acceptance

      bool ptgg_cut_default(1), phot1_cut_default(1), phot2_cut_default(1), jet1_cut_default(1), jet2_cut_default(1), deltae_cut_default(1), zep_cut_default(1), mjj_cut_default(1), deltap_cut_default(1), met_cut_default(1);
      bool ptgg_cut[20], phot1_cut[20], phot2_cut[20], jet1_cut[20], jet2_cut[20], deltae_cut[20], zep_cut[20], mjj_cut[20], deltap_cut[20], met_cut[20];
      for (int i=0; i<20; i++){
	ptgg_cut[i]=1; phot1_cut[i]=1; phot2_cut[i]=1; jet1_cut[i]=1; jet2_cut[i]=1; deltae_cut[i]=1; zep_cut[i]=1; mjj_cut[i]=1; deltap_cut[i]=1, met_cut[i]=1;
      }

      ptgg_cut_default = ptgg>pthiggsmincut;//pt higgs

      phot1_cut_default = ptphot1>ptphot1cut* massggnewvtx/120.;//pt first photon

      phot2_cut_default = ptphot2>ptphot2cut;//pt second photon

      jet1_cut_default = (ptcorrjet1>ptjet1cut && TMath::Abs(etajet1)<4.7);//pt first jet

      jet2_cut_default = (ptcorrjet2>ptjet2cut && TMath::Abs(etajet1)<4.7);//pt second jet

      if(deltaetacut!=0){//delteta
	if(deltaetacut>0){
	  deltae_cut_default = TMath::Abs(deltaeta)>deltaetacut;  // vbf selection 
	}else{
	  deltae_cut_default = TMath::Abs(deltaeta)<-deltaetacut;  // WZH selection
	}
      }

      if(zeppencut!=0) //zeppenfeld
	zep_cut_default = TMath::Abs(zeppenjet)<zeppencut;       

      if(invmassjetcut!=0){//inv mass of jets
	if(invmassjetcut>0){
	  mjj_cut_default = invmassjet>invmassjetcut; // vbf selection 
	}else{
	  mjj_cut_default = TMath::Abs(invmassjet-85)<-invmassjetcut; // WZH selection
	}
      }

      deltap_cut_default = TMath::Abs(deltaphi)>deltaphicut; 

      //      met_cut_default = met>70;
      met_cut_default = 1;
      
      // CUT VARIATION
      // VBF
      // double ptgg_range(80), phot1_range(60), phot2_range(30), jet1_range(40), jet2_range(20), deltae_range(2), zep_range(2), mjj_range(800), deltap_range(1.0),  met_range(140);
      // double ptgg_central(40), phot1_central(65), phot2_central(35), jet1_central(40), jet2_central(30), deltae_central(3.5), zep_central(2.5), mjj_central(450), deltap_central(2.6), met_central(70);
      // WZH
     double ptgg_range(100), phot1_range(60), phot2_range(30), jet1_range(40), jet2_range(30), deltae_range(2), zep_range(3), mjj_range(60), deltap_range(1.0),  met_range(140);
     double ptgg_central(50), phot1_central(65), phot2_central(35), jet1_central(40), jet2_central(35), deltae_central(-2.5), zep_central(2.), mjj_central(-35), deltap_central(2.6), met_central(70);
      double pthiggsmincut_variation(1), ptphot1cut_variation(1), ptphot2cut_variation(1), ptjet1cut_variation(1), ptjet2cut_variation(1), deltaetacut_variation(1), zeppencut_variation(1), invmassjetcut_variation(1), deltaphicut_variation(1), metcut_variation(1);

      bool extramjjcut;
      //      extramjjcut = invmassjet<500;

      for(int i=0; i<20; i++){
	pthiggsmincut_variation = ptgg_central + ptgg_range * ( i - 10. ) / 20.;
	ptgg_cut[i] = ptgg>pthiggsmincut_variation;//pt higgs
	
	ptphot1cut_variation = phot1_central + phot1_range * ( i - 10. ) / 20.;
	phot1_cut[i] = ptphot1>ptphot1cut_variation* massggnewvtx/120.;//pt first photon
	
	ptphot2cut_variation = phot2_central + phot2_range * ( i - 10. ) / 20.;
	phot2_cut[i] = ptphot2>ptphot2cut_variation;//pt second photon
	
	ptjet1cut_variation = jet1_central + jet1_range * ( i - 10. ) / 20.;
	jet1_cut[i] = (ptcorrjet1>ptjet1cut_variation && TMath::Abs(etajet1)<4.7);//pt first jet
	
	ptjet2cut_variation = jet2_central + jet2_range * ( i - 10. ) / 20.;
	jet2_cut[i] = (ptcorrjet2>ptjet2cut_variation && TMath::Abs(etajet1)<4.7);//pt second jet
	
	deltaetacut_variation = deltae_central + deltae_range * ( i - 10. ) / 20.;
	if(deltaetacut_variation!=0){//delteta
	  if(deltaetacut>0){
	    deltae_cut[i] = TMath::Abs(deltaeta)>deltaetacut_variation;  // vbf selection 
	  }else{
	    deltae_cut[i] = TMath::Abs(deltaeta)<-deltaetacut_variation;  // WZH selection
	  }
	}
	
	zeppencut_variation = zep_central + zep_range * ( i - 10. ) / 20.;
	if(zeppencut_variation!=0){ //zeppenfeld
	  zep_cut[i] = TMath::Abs(zeppenjet)<zeppencut_variation;       
	}

	invmassjetcut_variation = mjj_central + mjj_range * ( i - 10. ) / 20.;
	if(invmassjetcut_variation!=0){//inv mass of jets
	  if(invmassjetcut>0){
	    mjj_cut[i] = invmassjet>invmassjetcut_variation; // vbf selection 
	  }else{
	    mjj_cut[i] = TMath::Abs(invmassjet-85)<-invmassjetcut_variation; // WZH selection
	  }
	}

	deltaphicut_variation = deltap_central + deltap_range * ( i - 10. ) / 20.;
	deltap_cut[i] = TMath::Abs(deltaphi)>deltaphicut_variation; 

	metcut_variation = met_central + met_range * ( i - 10. ) / 20.;
	met_cut[i] = met>metcut_variation; 

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

      double variable(0);
      variable = massggnewvtx;

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

       // applying all analysis cuts
      if( ptgg_cut_default && phot1_cut_default && phot2_cut_default && jet1_cut_default && jet2_cut_default && deltae_cut_default && zep_cut_default && mjj_cut_default && deltap_cut_default && met_cut_default) {

	tempplot->Fill(variable, weight);
	
	tree_.run=run;
	tree_.lumi=lumi;
	tree_.event=event;
	tree_.ptgg=ptgg;
	tree_.ebeb=ebeb;
	tree_.massggnewvtx=massggnewvtx;
	tree_.weight=weight;
	fOut->cd();
	myTree->Fill();
	
	if(writetxt != "") 
	  outfile << "run " << run << "\t lumi "  << std::setw(4) << lumi << "\t event " << std::setw(12) << event  << "\t massgg " << std::setprecision(6) << massggnewvtx << endl;      

      }
      for(int ii=0; ii<20; ii++){

	//       	if(!extramjjcut) continue;
	if( ptgg_cut[ii] && phot1_cut_default && phot2_cut_default && jet1_cut_default && jet2_cut_default && deltae_cut_default && zep_cut_default && mjj_cut_default && deltap_cut_default && met_cut_default) mass_ptgg[ii]->Fill(variable, weight);

	if( ptgg_cut_default && phot1_cut[ii] && phot2_cut_default && jet1_cut_default && jet2_cut_default && deltae_cut_default && zep_cut_default && mjj_cut_default && deltap_cut_default && met_cut_default) mass_phot1[ii]->Fill(variable, weight);
	
	if( ptgg_cut_default && phot1_cut_default && phot2_cut[ii] && jet1_cut_default && jet2_cut_default && deltae_cut_default && zep_cut_default && mjj_cut_default && deltap_cut_default && met_cut_default) mass_phot2[ii]->Fill(variable, weight);

	if( ptgg_cut_default && phot1_cut_default && phot2_cut_default && jet1_cut[ii] && jet2_cut_default && deltae_cut_default && zep_cut_default && mjj_cut_default && deltap_cut_default && met_cut_default) mass_jet1[ii]->Fill(variable, weight);

	if( ptgg_cut_default && phot1_cut_default && phot2_cut_default && jet1_cut_default && jet2_cut[ii] && deltae_cut_default && zep_cut_default && mjj_cut_default && deltap_cut_default && met_cut_default) mass_jet2[ii]->Fill(variable, weight);

	if( ptgg_cut_default && phot1_cut_default && phot2_cut_default && jet1_cut_default && jet2_cut_default && deltae_cut[ii] && zep_cut_default && mjj_cut_default && deltap_cut_default && met_cut_default) mass_deltae[ii]->Fill(variable, weight);

	if( ptgg_cut_default && phot1_cut_default && phot2_cut_default && jet1_cut_default && jet2_cut_default && deltae_cut_default && zep_cut[ii] && mjj_cut_default && deltap_cut_default && met_cut_default) mass_zep[ii]->Fill(variable, weight);

	if( ptgg_cut_default && phot1_cut_default && phot2_cut_default && jet1_cut_default && jet2_cut_default && deltae_cut_default && zep_cut_default && mjj_cut[ii] && deltap_cut_default && met_cut_default) mass_mjj[ii]->Fill(variable, weight);

	if( ptgg_cut_default && phot1_cut_default && phot2_cut_default && jet1_cut_default && jet2_cut_default && deltae_cut_default && zep_cut_default && mjj_cut_default && deltap_cut[ii] && met_cut_default) mass_deltap[ii]->Fill(variable, weight);

	if( ptgg_cut_default && phot1_cut_default && phot2_cut_default && jet1_cut_default && jet2_cut_default && deltae_cut_default && zep_cut_default && mjj_cut_default && deltap_cut_default && met_cut[ii]) mass_met[ii]->Fill(variable, weight);

      }

   }
   
   fOut->Write();

   if(writetxt != "") 
     outfile.close();

   return fOut;

}

void  fillHisto::Setcuts(double pt1, double pt2, double higgsptcutmin, double higgsptcutmax, double ptj1, double ptj2, double deltae, double zep, double mjj, double deltap, int eb, int r9, bool third)
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

bool fillHisto::exclSel(){

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

void fillHisto::setCic(int cic) {
  cicselection = cic;
}


bool fillHisto::cutIDEG(double ptPhot, double etascphot, double pid_hlwTrackNoDz, double pid_jurECAL, double pid_twrHCAL, double pid_HoverE, double pid_etawid, int scatrk, int scaecal, int scahcal, int scahove) {

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

void fillHisto::Writetxt(char * filename)
{
  writetxt=std::string(filename);
}

void fillHisto::DoPuReweight(){
  dopureweight = 1;
}

void fillHisto::DoPtReweight(){
  doptreweight = 1;
}

void fillHisto::SetPuWeights(bool isData,std::string puWeightFile)
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

void fillHisto::DoSmearing(double mean, double spread)
{
  dosmear = 1;
  meansmear = mean;
  spreadsmear = spread;
}



void fillHisto::getweights()
{

  TFile *f  = new TFile("ptreweight.root","READ");

  TH1D *puweights = 0;
  
  puweights= (TH1D*) f->Get("pt2d");

 for (int i = 0; i<15; i++) {
    for (int j = 0; j<20; j++) {
      float weight=1.;
      weight=puweights->GetBinContent(i+1,j+1);
      weights_[i][j] =  weight;
      cout << i << "  " << "  " << j << "   " << weight << endl; 
    }
  }
  
  //std::cout << "weights sum is " << sumPuWeights << std::endl;

}
