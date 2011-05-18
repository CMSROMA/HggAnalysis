#define fillPlot_cxx
#include "fillPlot.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <iostream>
#include <fstream>

TH1D * fillPlot::Plot(string var, string name, int nbin, double min, double max)
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
   if (writetxt) outfile.open("events.txt"); 


   for (Long64_t jentry=0; jentry<nentries;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;

      // analysis cuts

      if((TMath::Abs(etaphot1)>1.4442&&TMath::Abs(etaphot1)<1.566)||(TMath::Abs(etaphot2)>1.4442&&TMath::Abs(etaphot2)<1.566)) continue;  // acceptance

      if(ptphot1<ptphot1cut) continue; //pt first photon
      if(ptphot2<ptphot2cut) continue; //pt second photon

      if(ptjet1<ptjet1cut) continue; //pt first photon
      if(ptjet2<ptjet2cut) continue; //pt second photon

      if(TMath::Abs(deltaeta)<deltaetacut) continue; //delteta
      if(TMath::Abs(zeppenjet)>zeppencut) continue; //zeppenfeld

      if(invmassjet<invmassjetcut) continue; //inv mass of jets

      if((pid_haspixelseedphot1||pid_haspixelseedphot2) && pixelseedcut) continue; // pixel seed
      if(ebcat) { // EB EE categories
	if((TMath::Abs(etaphot1)>1.4442||TMath::Abs(etaphot2)>1.4442)) continue; 
      } else {
	if((TMath::Abs(etaphot1)<1.4442&&TMath::Abs(etaphot2)<1.4442)) continue; 
      }

      bool isr9phot1(0), isr9phot2(0);
      if(TMath::Abs(etaphot1)<1.4442 && r9phot1>.93) isr9phot1 = 1;
      if(TMath::Abs(etaphot2)<1.4442 && r9phot2>.93) isr9phot2 = 1;
      if(TMath::Abs(etaphot1)>1.4442 && r9phot1>.9) isr9phot1 = 1;
      if(TMath::Abs(etaphot2)>1.4442 && r9phot2>.9) isr9phot2 = 1;
      if(r9cat) { // r9 categories
	if(!isr9phot1 || !isr9phot2) continue;
      } else {
	if(isr9phot1 && isr9phot2) continue;
      }

      // photon id
      if(!cutIDEG(ptphot1, etaphot1, pid_hlwTrackNoDzphot1, pid_jurECALphot1, pid_twrHCALphot1, pid_HoverEphot1, pid_etawidphot1, scaletrk, scaleecal, scalehcal, scalehove)) 
	continue;
      if(!cutIDEG(ptphot2, etaphot2, pid_hlwTrackNoDzphot2, pid_jurECALphot2, pid_twrHCALphot2, pid_HoverEphot2, pid_etawidphot2, scaletrk, scaleecal, scalehcal, scalehove)) 
	continue;

      double variable(0);
      if (var == "massgg")  variable = massgg;
      else if (var == "ptphot1")  variable = ptphot1;
      else if (var == "ptphot2")  variable = ptphot2;
      else if (var == "ptjet1")  variable = ptjet1;
      else if (var == "ptjet2")  variable = ptjet2;
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
      tempplot->Fill(variable);
      if(writetxt) 
	outfile << run << " : " << lumi << " : " << event << " : "<< ptphot1 << " : " << ptphot2<< " : " << massgg<< endl;      

   }
   
   if(writetxt) outfile.close();

   return tempplot;

}

void  fillPlot::Setcuts(double pt1, double pt2, double ptj1, double ptj2, double deltae, double zep, double mjj, bool eb, bool r9, int isolscaletrk, int isolscaleecal, int isolscalehcal, int isolscalehove, bool pixelseedveto)
{
  ptphot1cut = pt1;
  ptphot2cut = pt2;
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

bool fillPlot::cutIDEG(double ptPhot, double etaPhot, double pid_hlwTrackNoDz, double pid_jurECAL, double pid_twrHCAL, double pid_HoverE, double pid_etawid, int scatrk, int scaecal, int scahcal, int scahove) {

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

  if(TMath::Abs(etaPhot) < 1.4442) {
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

  if(TMath::Abs(etaPhot) > 1.4442) {
    setaeta = pid_etawid < setaetaEE;
  }  


  return (ptiso && hcaliso && ecaliso && hoveiso && setaeta);
}

void fillPlot::Writetxt(bool value)
{
  writetxt = value;
}
