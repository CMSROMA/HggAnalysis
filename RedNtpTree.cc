#include "RedNtpTree.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <iostream>
#include <vector>
#include <TLorentzVector.h>
using std::cout;
using std::endl;


inline double delta_phi(double phi1, double phi2) {

  double dphi = fabs(phi1 - phi2);
  return (dphi <= TMath::Pi())? dphi : TMath::TwoPi() - dphi;

}

inline double delta_eta(double eta1, double eta2) {

  return (eta2 >= 0 ? eta1 - eta2 : eta2 - eta1);
}

RedNtpTree::RedNtpTree(TTree *tree, const TString& outname) : tree_reader_V2(tree) 
{
   hOutputFile   = new TFile(outname, "RECREATE" ) ;
}

RedNtpTree::~RedNtpTree()
{
   hOutputFile->Write() ;
   hOutputFile->Close() ;
   hOutputFile->Delete();
}



vector<int>  RedNtpTree::firsttwo(Float_t *vec, vector<bool> *asso){

  double max(-999); int idmax(-999);
  double secondmax(-999); int idsecondmax(-999);

  for (int i=0; i<int(asso->size()); i++) {

    if ( vec[i] > max && asso->at(i)) {
      max = vec[i];
      idmax = i;
    }

  }
  for (int i=0; i<int(asso->size()); i++) {

    if ( vec[i] > secondmax && asso->at(i) && i!= idmax) {
      secondmax = vec[i];
      idsecondmax = i;
    }

  }
  
  //  int themaxtemp[] = {idmax,idsecondmax};

  vector<int> themax;
  
  themax.push_back(idmax);
  themax.push_back(idsecondmax);

  return themax;

}

bool RedNtpTree::cutID(int i, photonidcuts const& pid, vector<bool> *vpass) {

  // Use photon supercluster energy (would be e5x5 if r9>0.93 otherwise)
  bool ntrkiso = ntrkiso035Phot[i] < pid.tracknb;
  bool ptiso = (ptiso035Phot[i] / ptPhot[i] < pid.trackiso_rel);
  bool ecaliso = (ecaliso04Phot[i] / ePhot[i] < pid.ecaliso_rel ||
                  ecaliso04Phot[i] < pid.ecaliso_abs);
  double fhcal = hcalovecal04Phot[i];
  bool hcaliso = (fhcal < pid.hcaliso_rel ||
                  fhcal*ptPhot[i] < pid.hcaliso_abs);
  bool smaj = sMajMajPhot[i] < pid.smajmaj;
  bool smin = sMinMinPhot[i] < pid.sminmin;
  bool smin_min = sMinMinPhot[i] > pid.sminmin_min;
  bool eta = TMath::Abs(etaPhot[i]) < 2.5; 


  if(TMath::Abs(etaPhot[i]) > 1.44) {
    smaj = 1; smin = 1; smin_min = 1;
  }
  
  if (vpass) {
    //assert((*vpass).size()==7);
    if((*vpass).size()!=7) { cout << "major failure! (*vpass).size()!=7.. die!" << endl; exit(0) ; }
    (*vpass)[0] = ntrkiso;
    (*vpass)[1] = ptiso;
    (*vpass)[2] = hcaliso;
    (*vpass)[3] = ecaliso;
    (*vpass)[4] = smaj;
    (*vpass)[5] = smin;
    (*vpass)[6] = smin_min; 
  }

  return (ntrkiso && ptiso && hcaliso && ecaliso && smaj && smin && smin_min && eta);
}


void RedNtpTree::Loop(int isgjet)
{
//   In a ROOT session, you can do:
//      Root > .L RedNtpTree.C
//      Root > RedNtpTree t
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

   Long64_t nentries = fChain->GetEntriesFast();

   Long64_t nbytes = 0, nb = 0;

   // hOutputFile = new TFile("output.root" , "RECREATE" ) ;
   TH1D higgsmasshiggsassreco("higgsmasshiggsassreco","higgsmasshiggsassreco", 100, 100.,150.);
   TH1D higgsmassassreco("higgsmassassreco","higgsmassassreco", 100, 100.,150.);
   TH1D higgsmasscutreco("higgsmasscutreco","higgsmasscutreco", 100, 100.,150.);
   TH1D higgsmasscutzeppreco("higgsmasscutzeppreco","higgsmasscutzeppreco", 100, 100.,150.);
   TH1D higgsmasscutzeppdijetreco("higgsmasscutzeppdijetreco","higgsmasscutzeppdijetreco", 100, 100.,150.);
   TH1D higgsmassisocutreco("higgsmassisocutreco","higgsmassisocutreco", 100, 100.,150.);
   TH1D higgsmassjustisocutreco("higgsmassjustisocutreco","higgsmassjustisocutreco", 100, 100.,150.);
   TH1D higgsmassisojetptcutreco("higgsmassisojetptcutreco","higgsmassisojetptcutreco", 100, 100.,150.);
   TH1D higgsmassisocutzeppreco("higgsmassisocutzeppreco","higgsmassisocutzeppreco", 100, 100.,150.);
   TH1D higgsmassisocutzeppdijetreco("higgsmassisocutzeppdijetreco","higgsmassisocutzeppdijetreco", 100, 100.,150.);
   TH1D higgsmassreco("higgsmassreco","higgsmassreco", 100, 100.,150.);
   TH1D higgsmassisoreco("higgsmassisoreco","higgsmassisoreco", 100, 100.,150.);
   TH1D higgsmassisorecocheck("higgsmassisorecocheck","higgsmassisorecocheck", 1000, 0.,1000.);
   TH1D higgsmassisocutrecofull("higgsmassisocutrecofull","higgsmassisocutrecofull",1000, 0.,1000.);
   TH1D higgsmassjustisocutrecofull("higgsmassjustisocutrecofull","higgsmassjustisocutrecofull",1000, 0.,1000.);
   TH1D higgsmassisojetptcutrecofull("higgsmassisojetptcutrecofull","higgsmassisojetptcutrecofull",1000, 0.,1000.);
   TH1D higgsmassisocutzepprecofull("higgsmassisocutzepprecofull","higgsmassisocutzepprecofull",1000, 0.,1000.);
   TH1D higgsmassisocutzeppdijetrecofull("higgsmassisocutzeppdijetrecofull","higgsmassisocutzeppdijetrecofull",1000, 0.,1000.);
   TH1D higgsmassrecofull("higgsmassrecofull","higgsmassrecofull",1000, 0.,1000.);
   TH1D higgsmassisorecofull("higgsmassisorecofull","higgsmassisorecofull",1000, 0.,1000.);
   TH1D higgsmassisorecocheckfull("higgsmassisorecocheckfull","higgsmassisorecocheckfull",1000, 0.,1000.);
   TH1D pthiggshiggsassreco("pthiggshiggsassreco","pthiggshiggsassreco", 100, 0.,250.);
   TH1D pthiggsassreco("pthiggsassreco","pthiggsassreco", 100, 0.,250.);
   TH1D pthiggsisoreco("pthiggsisoreco","pthiggsisoreco", 100, 0.,250.);
   TH1D ptphotgen1("ptphotgen1","ptphotgen1", 100, 0.,300.);
   TH1D ptphotgen2("ptphotgen2","ptphotgen2", 100, 0.,300.);
   TH1D ptphothiggsgen1("ptphothiggsgen1","ptphothiggsgen1", 100, 0.,300.);
   TH1D ptphothiggsgen2("ptphothiggsgen2","ptphothiggsgen2", 100, 0.,300.);
   TH1D ptphothiggsassreco1("ptphothiggsassreco1","ptphothiggsassreco1", 100, 0.,300.);
   TH1D ptphothiggsassreco2("ptphothiggsassreco2","ptphothiggsassreco2", 100, 0.,300.);
   TH1D ptphotassreco("ptphotassreco","ptphotassreco", 100, 0.,300.);
   TH1D ptphotassreco1("ptphotassreco1","ptphotassreco1", 100, 0.,300.);
   TH1D ptphotassreco2("ptphotassreco2","ptphotassreco2", 100, 0.,300.);
   TH1D ptphotisoreco("ptphotisoreco","ptphotisoreco", 100, 0.,300.);
   TH1D ptphotisoreco1("ptphotisoreco1","ptphotisoreco1", 100, 0.,300.);
   TH1D ptphotisoreco2("ptphotisoreco2","ptphotisoreco2", 100, 0.,300.);
   TH1D ptphotisoassreco("ptphotisoassreco","ptphotisoassreco", 100, 0.,300.);
   TH1D ptphotjetreco("ptphotjetreco","ptphotjetreco", 100, 0.,300.);
   TH1D ptphotisonotassreco("ptphotisonotassreco","ptphotisonotassreco", 100, 0.,300.);
   TH1D ptphotisojetreco("ptphotisojetreco","ptphotisojetreco", 100, 0.,300.);
   TH1D ptphotnotassreco("ptphotnotassreco","ptphotnotassreco", 100, 0.,300.);
   TH1D etaphotgen1("etaphotgen1","etaphotgen1", 100, -5.5,5.5);
   TH1D etaphotgen2("etaphotgen2","etaphotgen2", 100, -5.5,5.5);
   TH1D etaphothiggsgen1("etaphothiggsgen1","etaphothiggsgen1", 100, -5.5,5.5);
   TH1D etaphothiggsgen2("etaphothiggsgen2","etaphothiggsgen2", 100, -5.5,5.5);
   TH1D etaphothiggsassreco1("etaphothiggsassreco1","etaphothiggsassreco1", 100, -5.5,5.5);
   TH1D etaphothiggsassreco2("etaphothiggsassreco2","etaphothiggsassreco2", 100, -5.5,5.5);
   TH1D etaphotassreco1("etaphotassreco1","etaphotassreco1", 100, -5.5,5.5);
   TH1D etaphotassreco2("etaphotassreco2","etaphotassreco2", 100, -5.5,5.5);
   TH1D etaphotisoreco1("etaphotisoreco1","etaphotisoreco1", 100, -5.5,5.5);
   TH1D etaphotisoreco2("etaphotisoreco2","etaphotisoreco2", 100, -5.5,5.5);
   TH1D etaphotassreco("etaphotassreco","etaphotassreco", 100, -5.5,5.5);
   TH1D etaphotisoreco("etaphotisoreco","etaphotisoreco", 100, -5.5,5.5);
   TH1D etaphotisoassreco("etaphotisoassreco","etaphotisoassreco", 100, -5.5,5.5);
   TH1D etaphotisonotassreco("etaphotisonotassreco","etaphotisonotassreco", 100, -5.5,5.5);
   TH1D etaphotnotassreco("etaphotnotassreco","etaphotnotassreco", 100, -5.5,5.5);
   TH1D etaphotjetreco("etaphotjetreco","etaphotjetreco", 100, -5.5,5.5);
   TH1D etaphotisojetreco("etaphotisojetreco","etaphotisojetreco", 100, -5.5,5.5);
   TH1D ptjetgen1("ptjetgen1","ptjetgen1", 100, 0.,300.);
   TH1D ptjetgen2("ptjetgen2","ptjetgen2", 100, 0.,300.);
   TH1D ptjethiggsassreco1("ptjethiggsassreco1","ptjethiggsassreco1", 100, 0.,300.);
   TH1D ptjethiggsassreco2("ptjethiggsassreco2","ptjethiggsassreco2", 100, 0.,300.);
   TH1D ptjetassreco1("ptjetassreco1","ptjetassreco1", 100, 0.,300.);
   TH1D ptjetassreco2("ptjetassreco2","ptjetassreco2", 100, 0.,300.);
   TH1D ptjetreco1("ptjetreco1","ptjetreco1", 100, 0.,300.);
   TH1D ptjetreco2("ptjetreco2","ptjetreco2", 100, 0.,300.);
   TH1D ptjetisoreco1("ptjetisoreco1","ptjetisoreco1", 100, 0.,300.);
   TH1D ptjetisoreco2("ptjetisoreco2","ptjetisoreco2", 100, 0.,300.);
   TH1D deltaetajetgen("deltaetajetgen","deltaetajetgen", 100, -7.5,7.5);
   TH1D deltaetajetgencut("deltaetajetgencut","deltaetajetgencut", 100, -7.5,7.5);
   TH1D etajetgen1("etajetgen1","etajetgen1", 100, -5.5,5.5);
   TH1D etajetgen2("etajetgen2","etajetgen2", 100, -5.5,5.5);
   TH1D deltaetajethiggsassreco("deltaetajethiggsassreco","deltaetajethiggsassreco", 100, -7.5,7.5);
   TH1D etajethiggsassreco1("etajethiggsassreco1","etajethiggsassreco1", 100, -5.5,5.5);
   TH1D etajethiggsassreco2("etajethiggsassreco2","etajethiggsassreco2", 100, -5.5,5.5);
   TH1D zeppenjethiggsassreco1("zeppenjethiggsassreco1","zeppenjethiggsassreco1", 100, -5.5,5.5);
   TH1D zeppenjethiggsassreco2("zeppenjethiggsassreco2","zeppenjethiggsassreco2", 100, -5.5,5.5);
   TH1D deltaetajetassreco("deltaetajetassreco","deltaetajetassreco", 100, -7.5,7.5);
   TH1D etajetassreco1("etajetassreco1","etajetassreco1", 100, -5.5,5.5);
   TH1D etajetassreco2("etajetassreco2","etajetassreco2", 100, -5.5,5.5);
   TH1D zeppenjetassreco1("zeppenjetassreco1","zeppenjetassreco1", 100, -5.5,5.5);
   TH1D zeppenjetassreco2("zeppenjetassreco2","zeppenjetassreco2", 100, -5.5,5.5);
   TH1D zeppenhiggsassreco("zeppenhiggsassreco","zeppenhiggsassreco", 100, -5.5,5.5);
   TH1D dijetmassassreco("dijetmassassreco","dijetmassassreco", 100, 50.,1500.);
   TH1D deltaetajetreco("deltaetajetreco","deltaetajetreco", 100, -7.5,7.5);
   TH1D etajetreco1("etajetreco1","etajetreco1", 100, -5.5,5.5);
   TH1D etajetreco2("etajetreco2","etajetreco2", 100, -5.5,5.5);
   TH1D zeppenjetreco1("zeppenjetreco1","zeppenjetreco1", 100, -5.5,5.5);
   TH1D zeppenjetreco2("zeppenjetreco2","zeppenjetreco2", 100, -5.5,5.5);
   TH1D deltaetajetisoreco("deltaetajetisoreco","deltaetajetisoreco", 100, -7.5,7.5);
   TH1D etajetisoreco1("etajetisoreco1","etajetisoreco1", 100, -5.5,5.5);
   TH1D etajetisoreco2("etajetisoreco2","etajetisoreco2", 100, -5.5,5.5);
   TH1D zeppenjetisoreco1("zeppenjetisoreco1","zeppenjetisoreco1", 100, -5.5,5.5);
   TH1D zeppenjetisoreco2("zeppenjetisoreco2","zeppenjetisoreco2", 100, -5.5,5.5);
   TH1D zeppenhiggsisoreco("zeppenhiggsisoreco","zeppenhiggsisoreco", 100, -5.5,5.5);
   TH1D dijetmassisoreco("dijetmassisoreco","dijetmassisoreco", 100, 50.,1500.);

   // isolation variables
   TH1D hcalisoassphot_EB("hcalisoassphot_EB","hcalisoassphot_EB",100,0.,1.);
   TH1D ecalisoassphot_EB("ecalisoassphot_EB","ecalisoassphot_EB",100,0.,1.);
   TH1D ptisoassphot_EB("ptisoassphot_EB","ptisoassphot_EB",100,0.,1.);
   TH1D ntrkisoassphot_EB("ntrkisoassphot_EB","ntrkisoassphot_EB",10,0.,10);
   TH1D sminminclusassphot_EB("sminminclusassphot_EB","sminminclusassphot_EB",100,0.,1.);
   TH1D smaxmaxclusassphot_EB("smaxmaxclusassphot_EB","smaxmaxclusassphot_EB",100,0.,1.);
   TH1D hcalisoassjet_EB("hcalisoassjet_EB","hcalisoassjet_EB",100,0.,1.);
   TH1D ecalisoassjet_EB("ecalisoassjet_EB","ecalisoassjet_EB",100,0.,1.);
   TH1D ptisoassjet_EB("ptisoassjet_EB","ptisoassjet_EB",100,0.,1.);
   TH1D ntrkisoassjet_EB("ntrkisoassjet_EB","ntrkisoassjet_EB",10,0.,10);
   TH1D sminminclusassjet_EB("sminminclusassjet_EB","sminminclusassjet_EB",100,0.,1.);
   TH1D smaxmaxclusassjet_EB("smaxmaxclusassjet_EB","smaxmaxclusassjet_EB",100,0.,1.);
   TH1D hcalisoassphot_EE("hcalisoassphot_EE","hcalisoassphot_EE",100,0.,1.);
   TH1D ecalisoassphot_EE("ecalisoassphot_EE","ecalisoassphot_EE",100,0.,1.);
   TH1D ptisoassphot_EE("ptisoassphot_EE","ptisoassphot_EE",100,0.,1.);
   TH1D ntrkisoassphot_EE("ntrkisoassphot_EE","ntrkisoassphot_EE",10,0.,10);
   TH1D sminminclusassphot_EE("sminminclusassphot_EE","sminminclusassphot_EE",100,0.,1.);
   TH1D smaxmaxclusassphot_EE("smaxmaxclusassphot_EE","smaxmaxclusassphot_EE",100,0.,1.);
   TH1D hcalisoassjet_EE("hcalisoassjet_EE","hcalisoassjet_EE",100,0.,1.);
   TH1D ecalisoassjet_EE("ecalisoassjet_EE","ecalisoassjet_EE",100,0.,1.);
   TH1D ptisoassjet_EE("ptisoassjet_EE","ptisoassjet_EE",100,0.,1.);
   TH1D ntrkisoassjet_EE("ntrkisoassjet_EE","ntrkisoassjet_EE",10,0.,10);
   TH1D sminminclusassjet_EE("sminminclusassjet_EE","sminminclusassjet_EE",100,0.,1.);
   TH1D smaxmaxclusassjet_EE("smaxmaxclusassjet_EE","smaxmaxclusassjet_EE",100,0.,1.);
      
   ana_tree = new TTree ("AnaTree","Reduced tree for final analysis") ;
   ana_tree->Branch("run",&runRN,"run/I");
   ana_tree->Branch("event",&eventRN,"event/I");
   ana_tree->Branch("massgg",&massgg,"massgg/F");
   ana_tree->Branch("ptphot1",&ptphot1,"ptphot1/F");
   ana_tree->Branch("ptphot2",&ptphot2,"ptphot2/F");
   ana_tree->Branch("etaphot1",&etaphot1,"etaphot1/F");
   ana_tree->Branch("etaphot2",&etaphot2,"etaphot2/F");
   ana_tree->Branch("E1phot1",&E1phot1,"E1phot1/F");
   ana_tree->Branch("E1phot2",&E1phot2,"E1phot2/F");
   ana_tree->Branch("E9phot1",&E9phot1,"E9phot1/F");
   ana_tree->Branch("E9phot2",&E9phot2,"E9phot2/F");
   ana_tree->Branch("idlooseEGphot1",&idlooseEGphot1,"idlooseEGphot1/I");
   ana_tree->Branch("idlooseEGphot2",&idlooseEGphot2,"idlooseEGphot2/I");
   ana_tree->Branch("idtightEGphot1",&idtightEGphot1,"idtightEGphot1/I");
   ana_tree->Branch("idtightEGphot2",&idtightEGphot2,"idtightEGphot2/I");
   ana_tree->Branch("idloosephot1",&idloosephot1,"idloosephot1/I");
   ana_tree->Branch("idloosephot2",&idloosephot2,"idloosephot2/I");
   ana_tree->Branch("idmediumphot1",&idmediumphot1,"idmediumphot1/I");
   ana_tree->Branch("idmediumphot2",&idmediumphot2,"idmediumphot2/I");
   ana_tree->Branch("ptjet1",&ptjet1,"ptjet1/F");
   ana_tree->Branch("ptjet2",&ptjet2,"ptjet2/F");
   ana_tree->Branch("ptcorrjet1",&ptcorrjet1,"ptcorrjet1/F");
   ana_tree->Branch("ptcorrjet2",&ptcorrjet2,"ptcorrjet2/F");
   ana_tree->Branch("etajet1",&etajet1,"etajet1/F");
   ana_tree->Branch("etajet2",&etajet2,"etajet2/F");
   ana_tree->Branch("deltaeta",&deltaeta,"deltaeta/F");
   ana_tree->Branch("zeppenjet",&zeppenjet,"zeppenjet/F");
   ana_tree->Branch("invmassjet",&invmassjet,"invmassjet/F");
   ana_tree->Branch("nvtx",&nvtx,"nvtx/F");
   ana_tree->Branch("met",&met,"met/F");


   photonidcuts mediumid;
   mediumid.hcaliso_rel=         0.05;
   mediumid.hcaliso_abs=         2.4;
   mediumid.ecaliso_rel=         0.05;
   mediumid.ecaliso_abs=         2.4;
   mediumid.tracknb=             3.;
   mediumid.trackiso_rel=        0.10;
   mediumid.sminmin=             0.30;
   mediumid.sminmin_min=         0.15;
   mediumid.smajmaj=             0.35;
  
   photonidcuts looseid;
   looseid.hcaliso_rel=         0.10;
   looseid.hcaliso_abs=         4.;
   looseid.ecaliso_rel=         0.10;
   looseid.ecaliso_abs=         4.5;
   looseid.tracknb=             5.;
   looseid.trackiso_rel=        0.20;
   looseid.sminmin=             0.50;
   looseid.sminmin_min=         0.15;
   looseid.smajmaj=             0.60; 

   photonidcuts superlooseid;
   superlooseid.hcaliso_rel=         0.15;
   superlooseid.hcaliso_abs=         6.;
   superlooseid.ecaliso_rel=         0.15;
   superlooseid.ecaliso_abs=         7;
   superlooseid.tracknb=             7.;
   superlooseid.trackiso_rel=        0.30;
   superlooseid.sminmin=             0.60;
   superlooseid.sminmin_min=         0.15;
   superlooseid.smajmaj=             0.70; 


   //event based cuts
   double ptphot1cut = 50;
   double ptphot2cut = 30;
   double ptjet1cut = 20;
   double ptjet2cut = 15;
   double deltaetacut = 2.5;
   double zeppencut = 2.5;
   double dijetmasscut = 300;
   double deltaphicut = 2.;

   for (Long64_t jentry=0; jentry<nentries&& jentry<100000;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      // if (Cut(ientry) < 0) continue;
      if (jentry%1000 == 0) cout << "Event " << jentry << endl;

      vector<bool> photassocMC, photassocMChiggs;
      
      int counter(0);

      for(int i=0; i<nMC; i++){      

	if(pdgIdMC[i] == 22 && statusMC[i] == 3){
	  photassocMC.push_back(1);	
	  counter++;
	}
	else
	  photassocMC.push_back(0);

	if(pdgIdMC[i] == 22 && statusMC[i] == 3 && pdgIdMC[motherIDMC[i]] == 25)
	  photassocMChiggs.push_back(1);
	else
	  photassocMChiggs.push_back(0);
	
      }
      
      if(isgjet && counter > 1) continue; 

      vector<int> firsttwogenphot = firsttwo(ptMC,&photassocMC);
      vector<int> firsttwohiggsgenphot = firsttwo(ptMC,&photassocMChiggs);

      ptphotgen1.Fill(ptMC[firsttwogenphot.at(0)]);
      ptphotgen2.Fill(ptMC[firsttwogenphot.at(1)]);
      etaphotgen1.Fill(etaMC[firsttwogenphot.at(0)]);
      etaphotgen2.Fill(etaMC[firsttwogenphot.at(1)]);

      ptphothiggsgen1.Fill(ptMC[firsttwohiggsgenphot.at(0)]);
      ptphothiggsgen2.Fill(ptMC[firsttwohiggsgenphot.at(1)]);
      etaphothiggsgen1.Fill(etaMC[firsttwohiggsgenphot.at(0)]);
      etaphothiggsgen2.Fill(etaMC[firsttwohiggsgenphot.at(1)]);
     
      ptjetgen1.Fill(ptJetGen_akt5[0]);
      ptjetgen2.Fill(ptJetGen_akt5[1]);
      etajetgen1.Fill(etaJetGen_akt5[0]);
      etajetgen2.Fill(etaJetGen_akt5[1]);
   
      deltaetajetgen.Fill(etaJetGen_akt5[0]-etaJetGen_akt5[1]);
      if(etaJetGen_akt5[0]*etaJetGen_akt5[1]<0) deltaetajetgencut.Fill(etaJetGen_akt5[0]-etaJetGen_akt5[1]);
   
      vector<bool> assophothiggs;
      vector<bool> assophot;
      vector<bool> jetphot;
      vector<bool> isophot;
      vector<bool> isophotloose;
      vector<bool> isophotmedium;

      for(int i=0; i<nPhot; i++){
     
	bool assh(0);
	bool assp(0);
	bool assj(0);

	for(int j=0; j<nMC; j++){
	  
	  double DR;
	  if(photassocMChiggs.at(j)){
	    DR = sqrt(delta_eta(etaPhot[i],etaMC[j])*delta_eta(etaPhot[i],etaMC[j]) + 
		      delta_phi(phiPhot[i],phiMC[j])*delta_phi(phiPhot[i],phiMC[j]) ) ;
	    if( DR < .01 )  assh = 1; 
	  }

	  if(photassocMC.at(j)){
	    DR = sqrt(delta_eta(etaPhot[i],etaMC[j])*delta_eta(etaPhot[i],etaMC[j]) + 
		      delta_phi(phiPhot[i],phiMC[j])*delta_phi(phiPhot[i],phiMC[j]) ) ;
	    if(DR < .01 ) assp = 1; 
	  }
	  
	}

	for(int j=0; j<nJet_pfakt5; j++){

	  double DR = sqrt(delta_eta(etaPhot[i],etaJetGen_akt5[j])*delta_eta(etaPhot[i],etaJetGen_akt5[j]) + 
			   delta_phi(phiPhot[i],phiJetGen_akt5[j])*delta_phi(phiPhot[i],phiJetGen_akt5[j]) ) ;
	  if(DR < .1) assj = 1; 
       
	}

	if(assh) assophothiggs.push_back(1);
	else assophothiggs.push_back(0); 

	if(assp) assophot.push_back(1); 
        else assophot.push_back(0);  

	if(assj) jetphot.push_back(1); 
        else jetphot.push_back(0);  

	vector<bool> idpass(7);
	//	if(cutID(i, mediumid, &idpass)) isophot.push_back(1); 
	//if(cutID(i, looseid, &idpass)) isophot.push_back(1); 
	// TEMP
  	if(cutID(i, superlooseid, &idpass)) isophot.push_back(1); 
	// END TEMP
        else isophot.push_back(0);  

	if(cutID(i, looseid, &idpass)) isophotloose.push_back(1); 
        else isophotloose.push_back(0);  

	if(cutID(i, mediumid, &idpass)) isophotmedium.push_back(1); 
        else isophotmedium.push_back(0);  
	
	if( assp )	{
	  ptphotassreco.Fill(ptPhot[i]);
	  etaphotassreco.Fill(etaPhot[i]);
	  if(ptPhot[i]>30) {
	    if(etaPhot[i]<1.47){
	      hcalisoassphot_EB.Fill(hcalovecal04Phot[i]);
	      ecalisoassphot_EB.Fill(ecaliso04Phot[i] / ePhot[i]);
	      ptisoassphot_EB.Fill(ptiso035Phot[i] / ptPhot[i]);
	      ntrkisoassphot_EB.Fill(ntrkiso035Phot[i]);
	      sminminclusassphot_EB.Fill(sMinMinPhot[i]);
	      smaxmaxclusassphot_EB.Fill(sMajMajPhot[i]);
	    }else if(etaPhot[i]<2.5){
	      hcalisoassphot_EE.Fill(hcalovecal04Phot[i]);
	      ecalisoassphot_EE.Fill(ecaliso04Phot[i] / ePhot[i]);
	      ptisoassphot_EE.Fill(ptiso035Phot[i] / ptPhot[i]);
	      ntrkisoassphot_EE.Fill(ntrkiso035Phot[i]);
	      sminminclusassphot_EE.Fill(sMinMinPhot[i]);
	      smaxmaxclusassphot_EE.Fill(sMajMajPhot[i]);	    
	    }
          }
	}
	if( isophot.at(i) )	{
	  ptphotisoreco.Fill(ptPhot[i]);
	  etaphotisoreco.Fill(etaPhot[i]);
	}
	if( assp && isophot.at(i) )	{
	  ptphotisoassreco.Fill(ptPhot[i]);
	  etaphotisoassreco.Fill(etaPhot[i]);
	}
	if( !assp )	{
	  ptphotnotassreco.Fill(ptPhot[i]);
	  etaphotnotassreco.Fill(etaPhot[i]);
	}
	if( !assp && isophot.at(i) )	{
	  ptphotisonotassreco.Fill(ptPhot[i]);
	  etaphotisonotassreco.Fill(etaPhot[i]);
	}
	if( assj )	{
	  ptphotjetreco.Fill(ptPhot[i]);
	  etaphotjetreco.Fill(etaPhot[i]);
	  if(ptPhot[i]>30) {
	    if(etaPhot[i]<1.47){
	      hcalisoassjet_EB.Fill(hcalovecal04Phot[i]);
	      ecalisoassjet_EB.Fill(ecaliso04Phot[i] / ePhot[i]);
	      ptisoassjet_EB.Fill(ptiso035Phot[i] / ptPhot[i]);
	      ntrkisoassjet_EB.Fill(ntrkiso035Phot[i]);
	      sminminclusassjet_EB.Fill(sMinMinPhot[i]);
	      smaxmaxclusassjet_EB.Fill(sMajMajPhot[i]);
	    }else if(etaPhot[i]<2.5){
	      hcalisoassjet_EE.Fill(hcalovecal04Phot[i]);
	      ecalisoassjet_EE.Fill(ecaliso04Phot[i] / ePhot[i]);
	      ptisoassjet_EE.Fill(ptiso035Phot[i] / ptPhot[i]);
	      ntrkisoassjet_EE.Fill(ntrkiso035Phot[i]);
	      sminminclusassjet_EE.Fill(sMinMinPhot[i]);
	      smaxmaxclusassjet_EE.Fill(sMajMajPhot[i]);	    
	    }
          }
	}
	if( assj && isophot.at(i) )	{
	  ptphotisojetreco.Fill(ptPhot[i]);
	  etaphotisojetreco.Fill(etaPhot[i]);
	}
      }  

      vector<int> firsttwohiggsassphot = firsttwo(ptPhot,&assophothiggs);
      vector<int> firsttwoassphot = firsttwo(ptPhot,&assophot);
      vector<int> firsttwoisophot = firsttwo(ptPhot,&isophot);      
      
      vector<bool> jetnohiggsphot;
      vector<bool> jetnoassphot;
      vector<bool> jetnoisophot;

      for(int i=0; i<nJet_pfakt5; i++){

	bool assh(0);
	bool assp(0);
	bool assi(0);

	double DR;

	for(int k=0; k<2; k++){

	  DR = sqrt(delta_eta(etaJet_pfakt5[i],etaPhot[firsttwohiggsassphot.at(k)])*delta_eta(etaJet_pfakt5[i],etaPhot[firsttwohiggsassphot.at(k)]) + 
		    delta_phi(phiJet_pfakt5[i],phiPhot[firsttwohiggsassphot.at(k)])*delta_phi(phiJet_pfakt5[i],phiPhot[firsttwohiggsassphot.at(k)]) ) ;
	  if( DR < .25 ) assh = 1; 
	  
	  DR = sqrt(delta_eta(etaJet_pfakt5[i],etaPhot[firsttwoassphot.at(k)])*delta_eta(etaJet_pfakt5[i],etaPhot[firsttwoassphot.at(k)]) + 
		    delta_phi(phiJet_pfakt5[i],phiPhot[firsttwoassphot.at(k)])*delta_phi(phiJet_pfakt5[i],phiPhot[firsttwoassphot.at(k)]) ) ;
	  if( DR < .25 ) assp = 1; 
	  
	  DR = sqrt(delta_eta(etaJet_pfakt5[i],etaPhot[firsttwoisophot.at(k)])*delta_eta(etaJet_pfakt5[i],etaPhot[firsttwoisophot.at(k)]) + 
		    delta_phi(phiJet_pfakt5[i],phiPhot[firsttwoisophot.at(k)])*delta_phi(phiJet_pfakt5[i],phiPhot[firsttwoisophot.at(k)]) ) ;
	  if( DR < .25 ) assi = 1; 

	}

	if(!assh) jetnohiggsphot.push_back(1);
	else jetnohiggsphot.push_back(0); 

	if(!assp) jetnoassphot.push_back(1); 
        else jetnoassphot.push_back(0);  

	if(!assi) jetnoisophot.push_back(1); 
        else jetnoisophot.push_back(0);  

      }

      vector<int> firsttwonohiggsjet = firsttwo(ptJet_pfakt5,&jetnohiggsphot);
      vector<int> firsttwonoassjet = firsttwo(ptJet_pfakt5,&jetnoassphot);
      vector<int> firsttwonoisojet = firsttwo(ptJet_pfakt5,&jetnoisophot);      

      if( firsttwohiggsassphot.at(0)>-1 && firsttwohiggsassphot.at(1)>-1 ) { 

	TLorentzVector phot1, phot2;	
	phot1.SetPtEtaPhiE(ptPhot[firsttwohiggsassphot.at(0)],etaPhot[firsttwohiggsassphot.at(0)],phiPhot[firsttwohiggsassphot.at(0)],ePhot[firsttwohiggsassphot.at(0)]);
	phot2.SetPtEtaPhiE(ptPhot[firsttwohiggsassphot.at(1)],etaPhot[firsttwohiggsassphot.at(1)],phiPhot[firsttwohiggsassphot.at(1)],ePhot[firsttwohiggsassphot.at(1)]);

	TLorentzVector higgs = phot1 + phot2;
	
	higgsmasshiggsassreco.Fill(higgs.M());
	pthiggshiggsassreco.Fill(higgs.Pt());
  
	ptphothiggsassreco1.Fill(phot1.Pt());
	ptphothiggsassreco2.Fill(phot2.Pt());
	etaphothiggsassreco1.Fill(etaPhot[firsttwohiggsassphot.at(0)]);
	etaphothiggsassreco2.Fill(etaPhot[firsttwohiggsassphot.at(1)]);

      }

      double higgsmass(0), etahiggs(-999);

      if( firsttwoassphot.at(0)>-1 && firsttwoassphot.at(1)>-1 ) { 

	TLorentzVector phot1, phot2;	
	phot1.SetPtEtaPhiE(ptPhot[firsttwoassphot.at(0)],etaPhot[firsttwoassphot.at(0)],phiPhot[firsttwoassphot.at(0)],ePhot[firsttwoassphot.at(0)]);
	phot2.SetPtEtaPhiE(ptPhot[firsttwoassphot.at(1)],etaPhot[firsttwoassphot.at(1)],phiPhot[firsttwoassphot.at(1)],ePhot[firsttwoassphot.at(1)]);

	TLorentzVector higgs = phot1 + phot2;

	higgsmass = higgs.M();
	etahiggs = higgs.Eta();
	
	higgsmassassreco.Fill(higgs.M());
	pthiggsassreco.Fill(higgs.Pt());
  
	ptphotassreco1.Fill(phot1.Pt());
	ptphotassreco2.Fill(phot2.Pt());
	etaphotassreco1.Fill(etaPhot[firsttwoassphot.at(0)]);
	etaphotassreco2.Fill(etaPhot[firsttwoassphot.at(1)]);

      }	

      double higgsisomass(0), etahiggsiso(-999);

      if( firsttwoisophot.at(0)>-1 && firsttwoisophot.at(1)>-1 ) { 

	TLorentzVector phot1, phot2;	
	phot1.SetPtEtaPhiE(ptPhot[firsttwoisophot.at(0)],etaPhot[firsttwoisophot.at(0)],phiPhot[firsttwoisophot.at(0)],ePhot[firsttwoisophot.at(0)]);
	phot2.SetPtEtaPhiE(ptPhot[firsttwoisophot.at(1)],etaPhot[firsttwoisophot.at(1)],phiPhot[firsttwoisophot.at(1)],ePhot[firsttwoisophot.at(1)]);

	TLorentzVector higgs = phot1 + phot2;

	higgsisomass = higgs.M();
	etahiggsiso = higgs.Eta();
	
	higgsmassisoreco.Fill(higgs.M());
	higgsmassisorecofull.Fill(higgs.M());
	pthiggsisoreco.Fill(higgs.Pt());
  
	ptphotisoreco1.Fill(phot1.Pt());
	ptphotisoreco2.Fill(phot2.Pt());
	etaphotisoreco1.Fill(etaPhot[firsttwoisophot.at(0)]);
	etaphotisoreco2.Fill(etaPhot[firsttwoisophot.at(1)]);

      }	

      if( firsttwohiggsassphot.at(0)>-1 && firsttwohiggsassphot.at(1)>-1 ) { 

	if( firsttwonohiggsjet.at(0) > -1) {
	  ptjethiggsassreco1.Fill(ptJet_pfakt5[firsttwonohiggsjet.at(0)]);
	  etajethiggsassreco1.Fill(etaJet_pfakt5[firsttwonohiggsjet.at(0)]);
	}
	if( firsttwonohiggsjet.at(1) > -1) {
	  ptjethiggsassreco2.Fill(ptJet_pfakt5[firsttwonohiggsjet.at(1)]);
	  etajethiggsassreco2.Fill(etaJet_pfakt5[firsttwonohiggsjet.at(1)]);
	}
	if( firsttwonohiggsjet.at(0) > -1 && firsttwonohiggsjet.at(1) > -1) {
	  deltaetajethiggsassreco.Fill(etaJet_pfakt5[firsttwonohiggsjet.at(0)]-etaJet_pfakt5[firsttwonohiggsjet.at(1)]);
	  double aveeta = (etaJet_pfakt5[firsttwonohiggsjet.at(0)]+etaJet_pfakt5[firsttwonohiggsjet.at(1)])/2;
	  double zeppen1 = etaJet_pfakt5[firsttwonohiggsjet.at(0)] - aveeta;
	  double zeppen2 = etaJet_pfakt5[firsttwonohiggsjet.at(1)] - aveeta;
	  zeppenjethiggsassreco1.Fill(zeppen1);
	  zeppenjethiggsassreco2.Fill(zeppen2);	
	}

      }
      
      double twojetsmass(0), etatwojets(-999);

      if( firsttwoassphot.at(0)>-1 && firsttwoassphot.at(1)>-1 ) { 

	if( firsttwonoassjet.at(0) > -1) {
	  ptjetassreco1.Fill(ptJet_pfakt5[firsttwonoassjet.at(0)]);
	  etajetassreco1.Fill(etaJet_pfakt5[firsttwonoassjet.at(0)]);
	}
	if( firsttwonoassjet.at(1) > -1) {
	  ptjetassreco2.Fill(ptJet_pfakt5[firsttwonoassjet.at(1)]);
	  etajetassreco2.Fill(etaJet_pfakt5[firsttwonoassjet.at(1)]);
	}
	if( firsttwonoassjet.at(0) > -1 && firsttwonoassjet.at(1) > -1) {
	  TLorentzVector jet1, jet2;	
	  jet1.SetPtEtaPhiE(ptJet_pfakt5[firsttwonoassjet.at(0)],etaJet_pfakt5[firsttwonoassjet.at(0)],phiJet_pfakt5[firsttwonoassjet.at(0)],eJet_pfakt5[firsttwonoassjet.at(0)]);
	  jet2.SetPtEtaPhiE(ptJet_pfakt5[firsttwonoassjet.at(1)],etaJet_pfakt5[firsttwonoassjet.at(1)],phiJet_pfakt5[firsttwonoassjet.at(1)],eJet_pfakt5[firsttwonoassjet.at(1)]);
	  
	  TLorentzVector sum = jet1 + jet2;
	  
	  twojetsmass = sum.M();
	  etatwojets = sum.Eta();
	 
	  deltaetajetassreco.Fill(etaJet_pfakt5[firsttwonoassjet.at(0)]-etaJet_pfakt5[firsttwonoassjet.at(1)]);
	  double aveeta = (etaJet_pfakt5[firsttwonoassjet.at(0)]+etaJet_pfakt5[firsttwonoassjet.at(1)])/2;
	  double zeppen1 = etaJet_pfakt5[firsttwonoassjet.at(0)] - aveeta;
	  double zeppen2 = etaJet_pfakt5[firsttwonoassjet.at(1)] - aveeta;
	  zeppenjetassreco1.Fill(zeppen1);
	  zeppenjetassreco2.Fill(zeppen2);	
	}

	if(  ptJet_pfakt5[firsttwonoassjet.at(0)] > ptjet1cut && ptJet_pfakt5[firsttwonoassjet.at(1)] > ptjet2cut 
	    && ptPhot[firsttwoassphot.at(0)] > ptphot1cut && ptPhot[firsttwoassphot.at(1)] > ptphot2cut ){

	  //	  deltaetajetreco.Fill(etaJet_pfakt5[firsttwonoassjet.at(0)]-etaJet_pfakt5[firsttwonoassjet.at(1)]);
	  //	  double zeppen = higgsreco_pt - aveeta;
	  if(TMath::Abs(etaJet_pfakt5[firsttwonoassjet.at(0)]-etaJet_pfakt5[firsttwonoassjet.at(1)])>deltaetacut){
	    higgsmasscutreco.Fill(higgsmass);
	    if(isophot.at(firsttwoassphot.at(0)) && isophot.at(firsttwoassphot.at(1)))  higgsmassisorecocheck.Fill(higgsisomass);	    
	    double aveeta = (etaJet_pfakt5[firsttwonoassjet.at(0)]+etaJet_pfakt5[firsttwonoassjet.at(1)])/2;
	    double zeppen = etahiggs - aveeta;
	    zeppenhiggsassreco.Fill(zeppen);
	    if(TMath::Abs(zeppen)<zeppencut) {
	      higgsmasscutzeppreco.Fill(higgsmass);
	      dijetmassassreco.Fill(twojetsmass);
	      if(twojetsmass>250){
		higgsmasscutzeppdijetreco.Fill(higgsmass);
	      }
	    }
	  }
	}

      }

      double twojetsmassiso(0), etatwojetsiso(-999);
      
      if( firsttwoisophot.at(0)>-1 && firsttwoisophot.at(1)>-1 ) { 
	
	if( firsttwonoisojet.at(0) > -1) {
	  ptjetisoreco1.Fill(ptJet_pfakt5[firsttwonoisojet.at(0)]);
	  etajetisoreco1.Fill(etaJet_pfakt5[firsttwonoisojet.at(0)]);
	}
	if( firsttwonoisojet.at(1) > -1) {
	  ptjetisoreco2.Fill(ptJet_pfakt5[firsttwonoisojet.at(1)]);
	  etajetisoreco2.Fill(etaJet_pfakt5[firsttwonoisojet.at(1)]);
	}
	if( firsttwonoisojet.at(0) > -1 && firsttwonoisojet.at(1) > -1) {
	  
	  TLorentzVector jet1, jet2;	
	  jet1.SetPtEtaPhiE(ptJet_pfakt5[firsttwonoisojet.at(0)],etaJet_pfakt5[firsttwonoisojet.at(0)],phiJet_pfakt5[firsttwonoisojet.at(0)],eJet_pfakt5[firsttwonoisojet.at(0)]);
	  jet2.SetPtEtaPhiE(ptJet_pfakt5[firsttwonoisojet.at(1)],etaJet_pfakt5[firsttwonoisojet.at(1)],phiJet_pfakt5[firsttwonoisojet.at(1)],eJet_pfakt5[firsttwonoisojet.at(1)]);
	  
	  TLorentzVector sum = jet1 + jet2;
	  
	  twojetsmassiso = sum.M();
	  etatwojetsiso = sum.Eta();

	  deltaetajetisoreco.Fill(etaJet_pfakt5[firsttwonoisojet.at(0)]-etaJet_pfakt5[firsttwonoisojet.at(1)]);
	  double aveeta = (etaJet_pfakt5[firsttwonoisojet.at(0)]+etaJet_pfakt5[firsttwonoisojet.at(1)])/2;
	  double zeppen1 = etaJet_pfakt5[firsttwonoisojet.at(0)] - aveeta;
	  double zeppen2 = etaJet_pfakt5[firsttwonoisojet.at(1)] - aveeta;
	  zeppenjetisoreco1.Fill(zeppen1);
	  zeppenjetisoreco2.Fill(zeppen2);	
	}

	massgg = higgsisomass;
	ptphot1 = ptPhot[firsttwoisophot.at(0)];
	ptphot2 = ptPhot[firsttwoisophot.at(1)];  
	etaphot1 = etaPhot[firsttwoisophot.at(0)];
	etaphot2 = etaPhot[firsttwoisophot.at(1)];  
	E1phot1 = E1Phot[firsttwoisophot.at(0)];
	E1phot2 = E1Phot[firsttwoisophot.at(1)];  
	E9phot1 = E9Phot[firsttwoisophot.at(0)];
	E9phot2 = E9Phot[firsttwoisophot.at(1)];  
	idlooseEGphot1 = pid_isLoose[firsttwoisophot.at(0)];
	idlooseEGphot2 = pid_isLoose[firsttwoisophot.at(1)];
	idtightEGphot1 = pid_isTight[firsttwoisophot.at(0)];
	idtightEGphot2 = pid_isTight[firsttwoisophot.at(1)];
	idloosephot1 = isophotloose.at(firsttwoisophot.at(0));
	idloosephot2 = isophotloose.at(firsttwoisophot.at(1));
	idmediumphot1 = isophotmedium.at(firsttwoisophot.at(0));
	idmediumphot2 = isophotmedium.at(firsttwoisophot.at(1));
	if( firsttwonoisojet.at(0) > -1) {
	  ptjet1 = ptJet_pfakt5[firsttwonoisojet.at(0)];
	  ptcorrjet1 = ptCorrJet_pfakt5[firsttwonoisojet.at(0)];	  
	  etajet1 = etaJet_pfakt5[firsttwonoisojet.at(0)];
	}else{
	  ptjet1 = -999;
	  ptcorrjet1 = -999;
	  etajet1 = -999;	 
	}
	if( firsttwonoisojet.at(1) > -1) {
	  ptjet2 = ptJet_pfakt5[firsttwonoisojet.at(1)];
	  ptcorrjet2 = ptCorrJet_pfakt5[firsttwonoisojet.at(1)];	  
	  etajet2 = etaJet_pfakt5[firsttwonoisojet.at(1)];
	}else{
	  ptjet2 = -999;
	  ptcorrjet2 = -999;
	  etajet2 = -999;	 
	}
	if( firsttwonoisojet.at(0) > -1 && firsttwonoisojet.at(1) > -1) {
	  deltaeta = etaJet_pfakt5[firsttwonoisojet.at(0)]-etaJet_pfakt5[firsttwonoisojet.at(1)];
	  double aveeta = (etaJet_pfakt5[firsttwonoisojet.at(0)]+etaJet_pfakt5[firsttwonoisojet.at(1)])/2;
	  zeppenjet = etahiggsiso - aveeta;
	  invmassjet = twojetsmassiso;
	}
	met = epfMet;
	nvtx = nvertex;

        runRN = run;
        eventRN = event;
	
	ana_tree->Fill();

	if(ptPhot[firsttwoisophot.at(0)] > ptphot1cut && ptPhot[firsttwoisophot.at(1)] > ptphot2cut){
	  higgsmassjustisocutreco.Fill(higgsisomass);
	  higgsmassjustisocutrecofull.Fill(higgsisomass);
	}

	if( ptJet_pfakt5[firsttwonoisojet.at(0)] > ptjet1cut && ptJet_pfakt5[firsttwonoisojet.at(1)] > ptjet2cut 
	    && ptPhot[firsttwoisophot.at(0)] > ptphot1cut && ptPhot[firsttwoisophot.at(1)] > ptphot2cut){
	  higgsmassisojetptcutreco.Fill(higgsisomass);
	  higgsmassisojetptcutrecofull.Fill(higgsisomass);
	  deltaetajetreco.Fill(etaJet_pfakt5[firsttwonoisojet.at(0)]-etaJet_pfakt5[firsttwonoisojet.at(1)]);
	  //	  double zeppen = higgsreco_pt - aveeta;
	  if(TMath::Abs(etaJet_pfakt5[firsttwonoisojet.at(0)]-etaJet_pfakt5[firsttwonoisojet.at(1)])>deltaetacut){
	    higgsmassisocutreco.Fill(higgsisomass);
	    higgsmassisocutrecofull.Fill(higgsisomass);
	    double aveeta = (etaJet_pfakt5[firsttwonoisojet.at(0)]+etaJet_pfakt5[firsttwonoisojet.at(1)])/2;
	    double zeppen = etahiggsiso - aveeta;
	    zeppenhiggsisoreco.Fill(zeppen);
	    if(TMath::Abs(zeppen)<zeppencut) {
	      higgsmassisocutzeppreco.Fill(higgsisomass);
	      higgsmassisocutzepprecofull.Fill(higgsisomass);
	      dijetmassisoreco.Fill(twojetsmassiso);
	      if(twojetsmassiso>dijetmasscut){
		higgsmassisocutzeppdijetreco.Fill(higgsisomass);
		higgsmassisocutzeppdijetrecofull.Fill(higgsisomass);
	      }
	    }
	  }
	}

      }


   }

   hOutputFile->Write() ;

}
