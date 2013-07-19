//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Wed May 18 00:04:04 2011 by ROOT version 5.27/06
// from TTree AnaTree/Reduced tree for final analysis
// found on file: redntp.41xv7.preselection.v3/redntp_run2010-2011.root
//////////////////////////////////////////////////////////

#ifndef fillHisto_h
#define fillHisto_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TH1D.h>

#include <string>

class fillHisto {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

   struct diPhotonTree_structure_ {
     int run;
     int lumi;
     int event;
     float ptgg;
     int ebeb;
     float massggnewvtx;
     float weight;
   };
   
   diPhotonTree_structure_ tree_;

   //Cuts values
   double ptphot1cut;
   double ptphot2cut;
   double pthiggsmincut;
   double pthiggsmaxcut;
   double ptjet1cut;
   double ptjet2cut;
   double deltaetacut;
   double deltaphicut;
   double zeppencut;
   double invmassjetcut;
   int ebcat;
   int r9cat;
   int cicselection;
   bool thirdcat;
   double ptphot1cut_2;
   double ptphot2cut_2;
   double pthiggsmincut_2;
   double pthiggsmaxcut_2;

   // bool to decide if we want to write output txt file
   std::string writetxt;
   
   // bool for switching on smearing and smearing parameters
   bool dosmear,dopureweight,doptreweight;
   double meansmear, spreadsmear;

   // vector of pu weights
   std::vector<Double_t> puweights_;

   double weights_[15][15];

   // Declaration of leaf types
   Int_t           run;
   Int_t           event;
   Int_t           lumi;
   Float_t         rhoPF;
   Float_t         massgg;
   Float_t         massggnewvtx;
   Float_t         ptgg;
   Float_t         ptphot1;
   Float_t         ptphot2;
   Float_t         timephot1;
   Float_t         timephot2;
   Float_t         etaphot1;
   Float_t         etaphot2;
   Float_t         phiphot1;
   Float_t         phiphot2;
   Float_t         etascphot1;
   Float_t         etascphot2;
   Float_t         phiscphot1;
   Float_t         phiscphot2;
   Float_t         E1phot1;
   Float_t         E1phot2;
   Float_t         E9phot1;
   Float_t         E9phot2;
   Float_t         r9phot1;
   Float_t         r9phot2;
   Int_t           isemEGphot1;
   Int_t           isemEGphot2;
   Int_t           idloosenewEGphot1;
   Int_t           idloosenewEGphot2;
   Int_t           idloose006newEGphot1;
   Int_t           idloose006newEGphot2;
   Int_t           idtightnewEGphot1;
   Int_t           idtightnewEGphot2;
   Int_t           idhggtightnewEGphot1;
   Int_t           idhggtightnewEGphot2;
   Int_t           idloosenewpuEGphot1;
   Int_t           idloosenewpuEGphot2;
   Int_t           idtightnewpuEGphot1;
   Int_t           idtightnewpuEGphot2;
   Int_t           idhggtightnewpuEGphot1;
   Int_t           idhggtightnewpuEGphot2;
   Int_t           idcicphot1;
   Int_t           idcicphot2;
   Int_t           idlooseEGphot1;
   Int_t           idlooseEGphot2;
   Int_t           idtightEGphot1;
   Int_t           idtightEGphot2;
   Int_t           idloosephot1;
   Int_t           idloosephot2;
   Int_t           idmediumphot1;
   Int_t           idmediumphot2;
   Int_t           idloosecsphot1;
   Int_t           idloosecsphot2;
   Int_t           idmediumcsphot1;
   Int_t           idmediumcsphot2;
   Int_t           idelephot1;
   Int_t           idelephot2;
   Int_t           pid_isEMphot1;
   Int_t           pid_isEMphot2;
   Int_t           pid_haspixelseedphot1;
   Int_t           pid_haspixelseedphot2;
   Float_t         pid_jurECALphot1;
   Float_t         pid_jurECALphot2;
   Float_t         pid_twrHCALphot1;
   Float_t         pid_twrHCALphot2;
   Float_t         pid_HoverEphot1;
   Float_t         pid_HoverEphot2;
   Float_t         pid_hlwTrackphot1;
   Float_t         pid_hlwTrackphot2;
   Float_t         pid_etawidphot1;
   Float_t         pid_etawidphot2;
   Float_t         pid_hlwTrackNoDzphot1;
   Float_t         pid_hlwTrackNoDzphot2;
   Float_t         pid_hasMatchedConvphot1;
   Float_t         pid_hasMatchedConvphot2;
   Float_t         pid_hasMatchedPromptElephot1;
   Float_t         pid_hasMatchedPromptElephot2;
   Float_t         pid_sminphot1;
   Float_t         pid_sminphot2;
   Float_t         pid_smajphot1;
   Float_t         pid_smajphot2;
   Int_t           pid_ntrkphot1;
   Int_t           pid_ntrkphot2;
   Float_t         pid_ptisophot1;
   Float_t         pid_ptisophot2;
   Int_t           pid_ntrkcsphot1;
   Int_t           pid_ntrkcsphot2;
   Float_t         pid_ptisocsphot1;
   Float_t         pid_ptisocsphot2;
   Float_t         pid_ecalisophot1;
   Float_t         pid_ecalisophot2;
   Float_t         pid_hcalisophot1;
   Float_t         pid_hcalisophot2;
   Float_t         ptjet1;
   Float_t         ptjet2;
   Float_t         ptcorrjet1;
   Float_t         ptcorrjet2;
   Float_t         etajet1;
   Float_t         etajet2;
   Float_t         phijet1;
   Float_t         phijet2;
   Float_t         deltaeta;
   Float_t         deltaphi;
   Float_t         deltaphinewvtx;
   Float_t         zeppenjet;
   Float_t         invmassjet;
   Float_t         invmass2g1j;
   Float_t         invmass2g2j;
   Float_t         nvtx;
   Int_t           npu;
   Float_t         met;
   Int_t           NtotEvents;
   Float_t         xsection;
   Float_t         EquivLumi;
   Int_t           SampleID;
   Float_t         pu_weight;
   Float_t         pt_weight;

   // List of branches
   TBranch        *b_run;   //!
   TBranch        *b_event;   //!
   TBranch        *b_lumi;   //!
   TBranch        *b_rhoPF;   //!
   TBranch        *b_massgg;   //!                                                                                                                         
   TBranch        *b_massggnewvtx;   //!
   TBranch        *b_ptgg;   //!
   TBranch        *b_ptphot1;   //!
   TBranch        *b_ptphot2;   //!
   TBranch        *b_timephot1;   //!
   TBranch        *b_timephot2;   //!
   TBranch        *b_etaphot1;   //!
   TBranch        *b_etaphot2;   //!
   TBranch        *b_phiphot1;   //!
   TBranch        *b_phiphot2;   //!
   TBranch        *b_etascphot1;   //!
   TBranch        *b_etascphot2;   //!
   TBranch        *b_phiscphot1;   //!
   TBranch        *b_phiscphot2;   //!
   TBranch        *b_E1phot1;   //!
   TBranch        *b_E1phot2;   //!
   TBranch        *b_E9phot1;   //!
   TBranch        *b_E9phot2;   //!
   TBranch        *b_r9phot1;   //!
   TBranch        *b_r9phot2;   //!
   TBranch        *b_isemEGphot1;   //!
   TBranch        *b_isemEGphot2;   //!
   TBranch        *b_idloosenewEGphot1;   //!
   TBranch        *b_idloosenewEGphot2;   //!
   TBranch        *b_idloose006newEGphot1;   //!
   TBranch        *b_idloose006newEGphot2;   //!
   TBranch        *b_idtightnewEGphot1;   //!
   TBranch        *b_idtightnewEGphot2;   //!
   TBranch        *b_idhggtightnewEGphot1;   //!
   TBranch        *b_idhggtightnewEGphot2;   //!
   TBranch        *b_idloosenewpuEGphot1;   //!
   TBranch        *b_idloosenewpuEGphot2;   //!
   TBranch        *b_idtightnewpuEGphot1;   //!
   TBranch        *b_idtightnewpuEGphot2;   //!
   TBranch        *b_idhggtightnewpuEGphot1;   //!
   TBranch        *b_idhggtightnewpuEGphot2;   //!
   TBranch        *b_idcicphot1;   //!
   TBranch        *b_idcicphot2;   //!
   TBranch        *b_idlooseEGphot1;   //!
   TBranch        *b_idlooseEGphot2;   //!
   TBranch        *b_idtightEGphot1;   //!
   TBranch        *b_idtightEGphot2;   //!
   TBranch        *b_idloosephot1;   //!
   TBranch        *b_idloosephot2;   //!
   TBranch        *b_idmediumphot1;   //!
   TBranch        *b_idmediumphot2;   //!
   TBranch        *b_idloosecsphot1;   //!
   TBranch        *b_idloosecsphot2;   //!
   TBranch        *b_idmediumcsphot1;   //!
   TBranch        *b_idmediumcsphot2;   //!
   TBranch        *b_idelephot1;   //!
   TBranch        *b_idelephot2;   //!
   TBranch        *b_pid_isEMphot1;   //!
   TBranch        *b_pid_isEMphot2;   //!
   TBranch        *b_pid_haspixelseedphot1;   //!
   TBranch        *b_pid_haspixelseedphot2;   //!
   TBranch        *b_pid_jurECALphot1;   //!
   TBranch        *b_pid_jurECALphot2;   //!
   TBranch        *b_pid_twrHCALphot1;   //!
   TBranch        *b_pid_twrHCALphot2;   //!
   TBranch        *b_pid_HoverEphot1;   //!
   TBranch        *b_pid_HoverEphot2;   //!
   TBranch        *b_pid_hlwTrackphot1;   //!
   TBranch        *b_pid_hlwTrackphot2;   //!
   TBranch        *b_pid_etawidphot1;   //!
   TBranch        *b_pid_etawidphot2;   //!
   TBranch        *b_pid_hlwTrackNoDzphot1;   //!
   TBranch        *b_pid_hlwTrackNoDzphot2;   //!
   TBranch        *b_pid_hasMatchedConvphot1;   //!
   TBranch        *b_pid_hasMatchedConvphot2;   //!
   TBranch        *b_pid_hasMatchedPromptElephot1;   //!
   TBranch        *b_pid_hasMatchedPromptElephot2;   //!
   TBranch        *b_pid_sminphot1;   //!
   TBranch        *b_pid_sminphot2;   //!
   TBranch        *b_pid_smajphot1;   //!
   TBranch        *b_pid_smajphot2;   //!
   TBranch        *b_pid_ntrkphot1;   //!
   TBranch        *b_pid_ntrkphot2;   //!
   TBranch        *b_pid_ptisophot1;   //!
   TBranch        *b_pid_ptisophot2;   //!
   TBranch        *b_pid_ntrkcsphot1;   //!
   TBranch        *b_pid_ntrkcsphot2;   //!
   TBranch        *b_pid_ptisocsphot1;   //!
   TBranch        *b_pid_ptisocsphot2;   //!
   TBranch        *b_pid_ecalisophot1;   //!
   TBranch        *b_pid_ecalisophot2;   //!
   TBranch        *b_pid_hcalisophot1;   //!
   TBranch        *b_pid_hcalisophot2;   //!
   TBranch        *b_ptjet1;   //!
   TBranch        *b_ptjet2;   //!
   TBranch        *b_ptcorrjet1;   //!
   TBranch        *b_ptcorrjet2;   //!
   TBranch        *b_etajet1;   //!
   TBranch        *b_etajet2;   //!
   TBranch        *b_phijet1;   //!
   TBranch        *b_phijet2;   //!
   TBranch        *b_deltaeta;   //!
   TBranch        *b_deltaphi;   //!
   TBranch        *b_deltaphinewvtx;   //!
   TBranch        *b_zeppenjet;   //!
   TBranch        *b_invmassjet;   //!
   TBranch        *b_invmass2g1j;   //!
   TBranch        *b_invmass2g2j;   //!
   TBranch        *b_nvtx;   //!
   TBranch        *b_npu;   //!
   TBranch        *b_met;   //!
   TBranch        *b_NtotEvents;   //!
   TBranch        *b_xsection;   //!
   TBranch        *b_EquivLumi;   //!
   TBranch        *b_SampleID;   //!
   TBranch        *b_pu_weight;   //!
   TBranch        *b_pt_weight;   //!

   fillHisto(TTree *tree=0, bool isData=0);
   virtual ~fillHisto();
   virtual void     getweights();
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Setcuts(double pt1=50, double pt2=30, double pthiggsmin=-100, double pthiggsmax=-100, double ptj1=20, double ptj2=15, double deltae=2.5, double zep=2.5, double mjj=300, double deltap=2.6, int eb = 1, int r9 = 1, bool thirdcat = 0);
   virtual TFile*   File(char* writeRoot, bool cs=0);
   virtual bool     cutIDEG(double ptPhot, double etaPhot, double pid_hlwTrackNoDz, double pid_jurECAL, double pid_twrHCAL, double pid_HoverE, double pid_etawid, int scaletrk=100, int scaleecal=100, int scalehcal=100, int scalehove=100);
   virtual bool     exclSel();   
   virtual void     setCic(int cic=5);
   virtual void     Writetxt(char * filename);
   virtual void     SetPuWeights(bool isData = 0,std::string file = "");
   virtual void     DoSmearing(double mean, double spread);   
   virtual void     DoPuReweight();   
   virtual void     DoPtReweight();   
   virtual Bool_t   Notify();
};

#endif

#ifdef fillHisto_cxx
fillHisto::fillHisto(TTree *tree, bool isData)
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("redntp.41xv7.preselection.v3/redntp_run2010-2011.root");
      if (!f) {
         f = new TFile("redntp.41xv7.preselection.v3/redntp_run2010-2011.root");
      }
      tree = (TTree*)gDirectory->Get("AnaTree");

   }
   Init(tree);
   writetxt = "";
   dosmear = 0;
   dopureweight = 0;
   doptreweight = 0;
   cicselection = -1;
//   SetPuWeights(isData);
}

fillHisto::~fillHisto()
{

   if (!fChain) return;
   //   delete fChain->GetCurrentFile();   
}

Int_t fillHisto::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t fillHisto::LoadTree(Long64_t entry)
{
// Set the environment to read one entry
   if (!fChain) return -5;
   Long64_t centry = fChain->LoadTree(entry);
   if (centry < 0) return centry;
   if (!fChain->InheritsFrom(TChain::Class()))  return centry;
   TChain *chain = (TChain*)fChain;
   if (chain->GetTreeNumber() != fCurrent) {
      fCurrent = chain->GetTreeNumber();
      Notify();
   }
   return centry;
}

void fillHisto::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("run", &run, &b_run);
   fChain->SetBranchAddress("event", &event, &b_event);
   fChain->SetBranchAddress("lumi", &lumi, &b_lumi);
   fChain->SetBranchAddress("rhoPF", &rhoPF, &b_rhoPF);
   fChain->SetBranchAddress("massgg", &massgg, &b_massgg);
   fChain->SetBranchAddress("massggnewvtx", &massggnewvtx, &b_massggnewvtx);
   fChain->SetBranchAddress("ptgg", &ptgg, &b_ptgg);
   fChain->SetBranchAddress("ptphot1", &ptphot1, &b_ptphot1);
   fChain->SetBranchAddress("ptphot2", &ptphot2, &b_ptphot2);
   fChain->SetBranchAddress("timephot1", &timephot1, &b_timephot1);
   fChain->SetBranchAddress("timephot2", &timephot2, &b_timephot2);
   fChain->SetBranchAddress("etaphot1", &etaphot1, &b_etaphot1);
   fChain->SetBranchAddress("etaphot2", &etaphot2, &b_etaphot2);
   fChain->SetBranchAddress("phiphot1", &phiphot1, &b_phiphot1);
   fChain->SetBranchAddress("phiphot2", &phiphot2, &b_phiphot2);
   fChain->SetBranchAddress("etascphot1", &etascphot1, &b_etascphot1);
   fChain->SetBranchAddress("etascphot2", &etascphot2, &b_etascphot2);
   fChain->SetBranchAddress("phiscphot1", &phiscphot1, &b_phiscphot1);
   fChain->SetBranchAddress("phiscphot2", &phiscphot2, &b_phiscphot2);
   fChain->SetBranchAddress("E1phot1", &E1phot1, &b_E1phot1);
   fChain->SetBranchAddress("E1phot2", &E1phot2, &b_E1phot2);
   fChain->SetBranchAddress("E9phot1", &E9phot1, &b_E9phot1);
   fChain->SetBranchAddress("E9phot2", &E9phot2, &b_E9phot2);
   fChain->SetBranchAddress("r9phot1", &r9phot1, &b_r9phot1);
   fChain->SetBranchAddress("r9phot2", &r9phot2, &b_r9phot2);
   fChain->SetBranchAddress("isemEGphot1", &isemEGphot1, &b_isemEGphot1);
   fChain->SetBranchAddress("isemEGphot2", &isemEGphot2, &b_isemEGphot2);
   fChain->SetBranchAddress("idloosenewEGphot1", &idloosenewEGphot1, &b_idloosenewEGphot1);
   fChain->SetBranchAddress("idloosenewEGphot2", &idloosenewEGphot2, &b_idloosenewEGphot2);
   fChain->SetBranchAddress("idloose006newEGphot1", &idloose006newEGphot1, &b_idloose006newEGphot1);
   fChain->SetBranchAddress("idloose006newEGphot2", &idloose006newEGphot2, &b_idloose006newEGphot2);
   fChain->SetBranchAddress("idtightnewEGphot1", &idtightnewEGphot1, &b_idtightnewEGphot1);
   fChain->SetBranchAddress("idtightnewEGphot2", &idtightnewEGphot2, &b_idtightnewEGphot2);
   fChain->SetBranchAddress("idhggtightnewEGphot1", &idhggtightnewEGphot1, &b_idhggtightnewEGphot1);
   fChain->SetBranchAddress("idhggtightnewEGphot2", &idhggtightnewEGphot2, &b_idhggtightnewEGphot2);
   fChain->SetBranchAddress("idloosenewpuEGphot1", &idloosenewpuEGphot1, &b_idloosenewpuEGphot1);
   fChain->SetBranchAddress("idloosenewpuEGphot2", &idloosenewpuEGphot2, &b_idloosenewpuEGphot2);
   fChain->SetBranchAddress("idtightnewpuEGphot1", &idtightnewpuEGphot1, &b_idtightnewpuEGphot1);
   fChain->SetBranchAddress("idtightnewpuEGphot2", &idtightnewpuEGphot2, &b_idtightnewpuEGphot2);
   fChain->SetBranchAddress("idhggtightnewpuEGphot1", &idhggtightnewpuEGphot1, &b_idhggtightnewpuEGphot1);
   fChain->SetBranchAddress("idhggtightnewpuEGphot2", &idhggtightnewpuEGphot2, &b_idhggtightnewpuEGphot2);
   fChain->SetBranchAddress("idcicphot1", &idcicphot1, &b_idcicphot1);
   fChain->SetBranchAddress("idcicphot2", &idcicphot2, &b_idcicphot2);
   fChain->SetBranchAddress("idlooseEGphot1", &idlooseEGphot1, &b_idlooseEGphot1);
   fChain->SetBranchAddress("idlooseEGphot2", &idlooseEGphot2, &b_idlooseEGphot2);
   fChain->SetBranchAddress("idtightEGphot1", &idtightEGphot1, &b_idtightEGphot1);
   fChain->SetBranchAddress("idtightEGphot2", &idtightEGphot2, &b_idtightEGphot2);
   fChain->SetBranchAddress("idloosephot1", &idloosephot1, &b_idloosephot1);
   fChain->SetBranchAddress("idloosephot2", &idloosephot2, &b_idloosephot2);
   fChain->SetBranchAddress("idmediumphot1", &idmediumphot1, &b_idmediumphot1);
   fChain->SetBranchAddress("idmediumphot2", &idmediumphot2, &b_idmediumphot2);
   fChain->SetBranchAddress("idloosecsphot1", &idloosecsphot1, &b_idloosecsphot1);
   fChain->SetBranchAddress("idloosecsphot2", &idloosecsphot2, &b_idloosecsphot2);
   fChain->SetBranchAddress("idmediumcsphot1", &idmediumcsphot1, &b_idmediumcsphot1);
   fChain->SetBranchAddress("idmediumcsphot2", &idmediumcsphot2, &b_idmediumcsphot2);
   fChain->SetBranchAddress("idelephot1", &idelephot1, &b_idelephot1);
   fChain->SetBranchAddress("idelephot2", &idelephot2, &b_idelephot2);
   fChain->SetBranchAddress("pid_isEMphot1", &pid_isEMphot1, &b_pid_isEMphot1);
   fChain->SetBranchAddress("pid_isEMphot2", &pid_isEMphot2, &b_pid_isEMphot2);
   fChain->SetBranchAddress("pid_haspixelseedphot1", &pid_haspixelseedphot1, &b_pid_haspixelseedphot1);
   fChain->SetBranchAddress("pid_haspixelseedphot2", &pid_haspixelseedphot2, &b_pid_haspixelseedphot2);
   fChain->SetBranchAddress("pid_jurECALphot1", &pid_jurECALphot1, &b_pid_jurECALphot1);
   fChain->SetBranchAddress("pid_jurECALphot2", &pid_jurECALphot2, &b_pid_jurECALphot2);
   fChain->SetBranchAddress("pid_twrHCALphot1", &pid_twrHCALphot1, &b_pid_twrHCALphot1);
   fChain->SetBranchAddress("pid_twrHCALphot2", &pid_twrHCALphot2, &b_pid_twrHCALphot2);
   fChain->SetBranchAddress("pid_HoverEphot1", &pid_HoverEphot1, &b_pid_HoverEphot1);
   fChain->SetBranchAddress("pid_HoverEphot2", &pid_HoverEphot2, &b_pid_HoverEphot2);
   fChain->SetBranchAddress("pid_hlwTrackphot1", &pid_hlwTrackphot1, &b_pid_hlwTrackphot1);
   fChain->SetBranchAddress("pid_hlwTrackphot2", &pid_hlwTrackphot2, &b_pid_hlwTrackphot2);
   fChain->SetBranchAddress("pid_etawidphot1", &pid_etawidphot1, &b_pid_etawidphot1);
   fChain->SetBranchAddress("pid_etawidphot2", &pid_etawidphot2, &b_pid_etawidphot2);
   fChain->SetBranchAddress("pid_hlwTrackNoDzphot1", &pid_hlwTrackNoDzphot1, &b_pid_hlwTrackNoDzphot1);
   fChain->SetBranchAddress("pid_hlwTrackNoDzphot2", &pid_hlwTrackNoDzphot2, &b_pid_hlwTrackNoDzphot2);
   fChain->SetBranchAddress("pid_hasMatchedConvphot1", &pid_hasMatchedConvphot1, &b_pid_hasMatchedConvphot1);
   fChain->SetBranchAddress("pid_hasMatchedConvphot2", &pid_hasMatchedConvphot2, &b_pid_hasMatchedConvphot2);
   fChain->SetBranchAddress("pid_hasMatchedPromptElephot1", &pid_hasMatchedPromptElephot1, &b_pid_hasMatchedPromptElephot1);
   fChain->SetBranchAddress("pid_hasMatchedPromptElephot2", &pid_hasMatchedPromptElephot2, &b_pid_hasMatchedPromptElephot2);
   fChain->SetBranchAddress("pid_sminphot1", &pid_sminphot1, &b_pid_sminphot1);
   fChain->SetBranchAddress("pid_sminphot2", &pid_sminphot2, &b_pid_sminphot2);
   fChain->SetBranchAddress("pid_smajphot1", &pid_smajphot1, &b_pid_smajphot1);
   fChain->SetBranchAddress("pid_smajphot2", &pid_smajphot2, &b_pid_smajphot2);
   fChain->SetBranchAddress("pid_ntrkphot1", &pid_ntrkphot1, &b_pid_ntrkphot1);
   fChain->SetBranchAddress("pid_ntrkphot2", &pid_ntrkphot2, &b_pid_ntrkphot2);
   fChain->SetBranchAddress("pid_ptisophot1", &pid_ptisophot1, &b_pid_ptisophot1);
   fChain->SetBranchAddress("pid_ptisophot2", &pid_ptisophot2, &b_pid_ptisophot2);
   fChain->SetBranchAddress("pid_ntrkcsphot1", &pid_ntrkcsphot1, &b_pid_ntrkcsphot1);
   fChain->SetBranchAddress("pid_ntrkcsphot2", &pid_ntrkcsphot2, &b_pid_ntrkcsphot2);
   fChain->SetBranchAddress("pid_ptisocsphot1", &pid_ptisocsphot1, &b_pid_ptisocsphot1);
   fChain->SetBranchAddress("pid_ptisocsphot2", &pid_ptisocsphot2, &b_pid_ptisocsphot2);
   fChain->SetBranchAddress("pid_ecalisophot1", &pid_ecalisophot1, &b_pid_ecalisophot1);
   fChain->SetBranchAddress("pid_ecalisophot2", &pid_ecalisophot2, &b_pid_ecalisophot2);
   fChain->SetBranchAddress("pid_hcalisophot1", &pid_hcalisophot1, &b_pid_hcalisophot1);
   fChain->SetBranchAddress("pid_hcalisophot2", &pid_hcalisophot2, &b_pid_hcalisophot2);
   fChain->SetBranchAddress("ptjet1", &ptjet1, &b_ptjet1);
   fChain->SetBranchAddress("ptjet2", &ptjet2, &b_ptjet2);
   fChain->SetBranchAddress("ptcorrjet1", &ptcorrjet1, &b_ptcorrjet1);
   fChain->SetBranchAddress("ptcorrjet2", &ptcorrjet2, &b_ptcorrjet2);
   fChain->SetBranchAddress("etajet1", &etajet1, &b_etajet1);
   fChain->SetBranchAddress("etajet2", &etajet2, &b_etajet2);
   fChain->SetBranchAddress("phijet1", &phijet1, &b_phijet1);
   fChain->SetBranchAddress("phijet2", &phijet2, &b_phijet2);
   fChain->SetBranchAddress("deltaeta", &deltaeta, &b_deltaeta);
   fChain->SetBranchAddress("deltaphi", &deltaphi, &b_deltaphi);
   fChain->SetBranchAddress("deltaphinewvtx", &deltaphinewvtx, &b_deltaphinewvtx);
   fChain->SetBranchAddress("zeppenjet", &zeppenjet, &b_zeppenjet);
   fChain->SetBranchAddress("invmassjet", &invmassjet, &b_invmassjet);
   fChain->SetBranchAddress("invmass2g1j", &invmass2g1j, &b_invmass2g1j);
   fChain->SetBranchAddress("invmass2g2j", &invmass2g2j, &b_invmass2g2j);
   fChain->SetBranchAddress("nvtx", &nvtx, &b_nvtx);
   fChain->SetBranchAddress("npu", &npu, &b_npu);
   fChain->SetBranchAddress("met", &met, &b_met);
   fChain->SetBranchAddress("NtotEvents", &NtotEvents, &b_NtotEvents);
   fChain->SetBranchAddress("xsection", &xsection, &b_xsection);
   fChain->SetBranchAddress("EquivLumi", &EquivLumi, &b_EquivLumi);
   fChain->SetBranchAddress("SampleID", &SampleID, &b_SampleID);
   fChain->SetBranchAddress("pu_weight", &pu_weight, &b_pu_weight);
   fChain->SetBranchAddress("pt_weight", &pt_weight, &b_pt_weight);
   Notify();
}

Bool_t fillHisto::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

#endif // #ifdef fillHisto_cxx
