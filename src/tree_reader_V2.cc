#define tree_reader_V2_cxx
#include "tree_reader_V2.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>

void tree_reader_V2::Loop()
{
//   In a ROOT session, you can do:
//      Root > .L tree_reader_V2.C
//      Root > tree_reader_V2 t
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
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      // if (Cut(ientry) < 0) continue;
   }
}

tree_reader_V2::tree_reader_V2(TTree *tree)
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("dcap://cmsrm-se01.roma1.infn.it//pnfs/roma1.infn.it/data/cms/store/user/rahatlou/data/Cert_132440-149442_7TeV_StreamExpress_Collisions10_JSON_v3_NOVRERECO_2010B_smaller_new2/output_100_1_4fg.root");
      if (!f) {
         f = new TFile("dcap://cmsrm-se01.roma1.infn.it//pnfs/roma1.infn.it/data/cms/store/user/rahatlou/data/Cert_132440-149442_7TeV_StreamExpress_Collisions10_JSON_v3_NOVRERECO_2010B_smaller_new2/output_100_1_4fg.root");
         f->cd("dcap://cmsrm-se01.roma1.infn.it//pnfs/roma1.infn.it/data/cms/store/user/rahatlou/data/Cert_132440-149442_7TeV_StreamExpress_Collisions10_JSON_v3_NOVRERECO_2010B_smaller_new2/output_100_1_4fg.root:/myanalysis");
      }
      tree = (TTree*)gDirectory->Get("pippo");

   }
   Init(tree);
}

tree_reader_V2::~tree_reader_V2()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t tree_reader_V2::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t tree_reader_V2::LoadTree(Long64_t entry)
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

void tree_reader_V2::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set object pointer
   HLTNames = 0;
   HLTResults = 0;
   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("genpt", &genpt, &b_genpt);
   fChain->SetBranchAddress("isMC", &isMC, &b_isMC);
   fChain->SetBranchAddress("store", &store, &b_store);
   fChain->SetBranchAddress("lbn", &lbn, &b_lbn);
   fChain->SetBranchAddress("bx", &bx, &b_bx);
   fChain->SetBranchAddress("orbit", &orbit, &b_orbit);
   fChain->SetBranchAddress("run", &run, &b_run);
   fChain->SetBranchAddress("event", &event, &b_event);
   fChain->SetBranchAddress("nMC", &nMC, &b_nMC);
   fChain->SetBranchAddress("pdgIdMC", &pdgIdMC, &b_pdgIdMC);
   fChain->SetBranchAddress("statusMC", &statusMC, &b_statusMC);
   fChain->SetBranchAddress("motherIDMC", &motherIDMC, &b_motherIDMC);
   fChain->SetBranchAddress("ptMC ", &ptMC , &b_ptMC );
   fChain->SetBranchAddress("eMC  ", &eMC  , &b_eMC  );
   fChain->SetBranchAddress("etaMC", &etaMC, &b_etaMC);
   fChain->SetBranchAddress("phiMC", &phiMC, &b_phiMC);
   fChain->SetBranchAddress("nPhot", &nPhot, &b_nPhot);
   fChain->SetBranchAddress("ptPhot ", ptPhot , &b_ptPhot );
   fChain->SetBranchAddress("ePhot  ", ePhot  , &b_ePhot  );
   fChain->SetBranchAddress("escPhot  ", escPhot  , &b_escPhot  );
   fChain->SetBranchAddress("eseedPhot  ", eseedPhot  , &b_eseedPhot  );
   fChain->SetBranchAddress("etaPhot", etaPhot, &b_etaPhot);
   fChain->SetBranchAddress("phiPhot", phiPhot, &b_phiPhot);
   fChain->SetBranchAddress("timePhot", timePhot, &b_timePhot);
   fChain->SetBranchAddress("e4SwissCrossPhot", e4SwissCrossPhot, &b_e4SwissCrossPhot);
   fChain->SetBranchAddress("hasPixelSeedPhot", hasPixelSeedPhot, &b_hasPixelSeedPhot);
   fChain->SetBranchAddress("pid_isEM", pid_isEM, &b_pid_isEM);
   fChain->SetBranchAddress("pid_isLoose", pid_isLoose, &b_pid_isLoose);
   fChain->SetBranchAddress("pid_isTight", pid_isTight, &b_pid_isTight);
   fChain->SetBranchAddress("pid_jurECAL", pid_jurECAL, &b_pid_jurECAL);
   fChain->SetBranchAddress("pid_twrHCAL", pid_twrHCAL, &b_pid_twrHCAL);
   fChain->SetBranchAddress("pid_HoverE", pid_HoverE, &b_pid_HoverE);
   fChain->SetBranchAddress("pid_hlwTrack", pid_hlwTrack, &b_pid_hlwTrack);
   fChain->SetBranchAddress("pid_etawid", pid_etawid, &b_pid_etawid);
   fChain->SetBranchAddress("ptiso0015Phot", ptiso0015Phot, &b_ptiso0015Phot);
   fChain->SetBranchAddress("ntrkiso0015Phot", ntrkiso0015Phot, &b_ntrkiso0015Phot);
   fChain->SetBranchAddress("ptiso035Phot", ptiso035Phot, &b_ptiso035Phot);
   fChain->SetBranchAddress("ntrkiso035Phot", ntrkiso035Phot, &b_ntrkiso035Phot);
   fChain->SetBranchAddress("ptiso04Phot", ptiso04Phot, &b_ptiso04Phot);
   fChain->SetBranchAddress("ntrkiso04Phot", ntrkiso04Phot, &b_ntrkiso04Phot);
   fChain->SetBranchAddress("hcalovecal04Phot", hcalovecal04Phot, &b_hcalovecal04Phot);
   fChain->SetBranchAddress("ecaliso04Phot", ecaliso04Phot, &b_ecaliso04Phot);
   fChain->SetBranchAddress("sMajMajPhot", sMajMajPhot, &b_sMajMajPhot);
   fChain->SetBranchAddress("sMinMinPhot", sMinMinPhot, &b_sMinMinPhot);
   fChain->SetBranchAddress("alphaPhot", alphaPhot, &b_alphaPhot);
   fChain->SetBranchAddress("sEtaEtaPhot", sEtaEtaPhot, &b_sEtaEtaPhot);
   fChain->SetBranchAddress("sEtaPhiPhot", sEtaPhiPhot, &b_sEtaPhiPhot);
   fChain->SetBranchAddress("sPhiPhiPhot", sPhiPhiPhot, &b_sPhiPhiPhot);
   fChain->SetBranchAddress("E1Phot", E1Phot, &b_E1Phot);
   fChain->SetBranchAddress("E9Phot", E9Phot, &b_E9Phot);
   fChain->SetBranchAddress("E25Phot", E25Phot, &b_E25Phot);
   fChain->SetBranchAddress("ieleassocPhot", ieleassocPhot, &b_ieleassocPhot);
   fChain->SetBranchAddress("nElePhot", &nElePhot, &b_nElePhot);
   fChain->SetBranchAddress("pid_jurECALElePhot ", pid_jurECALElePhot , &b_pid_jurECALElePhot );
   fChain->SetBranchAddress("pid_twrHCALElePhot ", pid_twrHCALElePhot , &b_pid_twrHCALElePhot );
   fChain->SetBranchAddress("pid_HoverEElePhot ", pid_HoverEElePhot , &b_pid_HoverEElePhot );
   fChain->SetBranchAddress("pid_hlwTrackElePhot ", pid_hlwTrackElePhot , &b_pid_hlwTrackElePhot );
   fChain->SetBranchAddress("pid_etawidElePhot ", pid_etawidElePhot , &b_pid_etawidElePhot );
   fChain->SetBranchAddress("pid_dphivtxElePhot ", pid_dphivtxElePhot , &b_pid_dphivtxElePhot );
   fChain->SetBranchAddress("pid_detavtxElePhot ", pid_detavtxElePhot , &b_pid_detavtxElePhot );
   fChain->SetBranchAddress("pid_dcotElePhot ", pid_dcotElePhot , &b_pid_dcotElePhot );
   fChain->SetBranchAddress("pid_distElePhot ", pid_distElePhot , &b_pid_distElePhot );
   fChain->SetBranchAddress("pid_mishitsElePhot ", pid_mishitsElePhot , &b_pid_mishitsElePhot );
   fChain->SetBranchAddress("pid_ptElePhot ", pid_ptElePhot , &b_pid_ptElePhot );
//   fChain->SetBranchAddress("nJet_akt5", &nJet_akt5, &b_nJet_akt5);
//   fChain->SetBranchAddress("ptJet_akt5 ", ptJet_akt5 , &b_ptJet_akt5 );
//   fChain->SetBranchAddress("ptCorrJet_akt5 ", ptCorrJet_akt5 , &b_ptCorrJet_akt5 );
//   fChain->SetBranchAddress("eJet_akt5  ", eJet_akt5  , &b_eJet_akt5  );
//   fChain->SetBranchAddress("etaJet_akt5", etaJet_akt5, &b_etaJet_akt5);
//   fChain->SetBranchAddress("phiJet_akt5", phiJet_akt5, &b_phiJet_akt5);
//   fChain->SetBranchAddress("emfJet_akt5", emfJet_akt5, &b_emfJet_akt5);
//   fChain->SetBranchAddress("n90Jet_akt5", n90Jet_akt5, &b_n90Jet_akt5);
//   fChain->SetBranchAddress("n90HitsJet_akt5", n90HitsJet_akt5, &b_n90HitsJet_akt5);
//   fChain->SetBranchAddress("fHPDJet_akt5", fHPDJet_akt5, &b_fHPDJet_akt5);
//   fChain->SetBranchAddress("fRBXJet_akt5", fRBXJet_akt5, &b_fRBXJet_akt5);
//   fChain->SetBranchAddress("nJet_akt7", &nJet_akt7, &b_nJet_akt7);
//   fChain->SetBranchAddress("ptJet_akt7 ", ptJet_akt7 , &b_ptJet_akt7 );
//   fChain->SetBranchAddress("ptCorrJet_akt7 ", ptCorrJet_akt7 , &b_ptCorrJet_akt7 );
//   fChain->SetBranchAddress("eJet_akt7  ", eJet_akt7  , &b_eJet_akt7  );
//   fChain->SetBranchAddress("etaJet_akt7", etaJet_akt7, &b_etaJet_akt7);
//   fChain->SetBranchAddress("phiJet_akt7", phiJet_akt7, &b_phiJet_akt7);
//   fChain->SetBranchAddress("emfJet_akt7", emfJet_akt7, &b_emfJet_akt7);
//   fChain->SetBranchAddress("n90Jet_akt7", n90Jet_akt7, &b_n90Jet_akt7);
//   fChain->SetBranchAddress("n90HitsJet_akt7", n90HitsJet_akt7, &b_n90HitsJet_akt7);
//   fChain->SetBranchAddress("fHPDJet_akt7", fHPDJet_akt7, &b_fHPDJet_akt7);
//   fChain->SetBranchAddress("fRBXJet_akt7", fRBXJet_akt7, &b_fRBXJet_akt7);
   fChain->SetBranchAddress("nJet_jptak5", &nJet_jptak5, &b_nJet_jptak5);
   fChain->SetBranchAddress("ptJet_jptak5 ", ptJet_jptak5 , &b_ptJet_jptak5 );
   fChain->SetBranchAddress("eJet_jptak5  ", eJet_jptak5  , &b_eJet_jptak5  );
   fChain->SetBranchAddress("etaJet_jptak5", etaJet_jptak5, &b_etaJet_jptak5);
   fChain->SetBranchAddress("phiJet_jptak5", phiJet_jptak5, &b_phiJet_jptak5);
   fChain->SetBranchAddress("emfJet_jptak5", emfJet_jptak5, &b_emfJet_jptak5);
   fChain->SetBranchAddress("nJet_pfkt4", &nJet_pfkt4, &b_nJet_pfkt4);
   fChain->SetBranchAddress("ptJet_pfkt4 ", ptJet_pfkt4 , &b_ptJet_pfkt4 );
   fChain->SetBranchAddress("eJet_pfkt4  ", eJet_pfkt4  , &b_eJet_pfkt4  );
   fChain->SetBranchAddress("etaJet_pfkt4", etaJet_pfkt4, &b_etaJet_pfkt4);
   fChain->SetBranchAddress("phiJet_pfkt4", phiJet_pfkt4, &b_phiJet_pfkt4);
   fChain->SetBranchAddress("nJet_pfakt5", &nJet_pfakt5, &b_nJet_pfakt5);
   fChain->SetBranchAddress("ptJet_pfakt5 ", ptJet_pfakt5 , &b_ptJet_pfakt5 );
   fChain->SetBranchAddress("ptCorrJet_pfakt5 ", ptCorrJet_pfakt5 , &b_ptCorrJet_pfakt5 );
   fChain->SetBranchAddress("eJet_pfakt5  ", eJet_pfakt5  , &b_eJet_pfakt5  );
   fChain->SetBranchAddress("etaJet_pfakt5", etaJet_pfakt5, &b_etaJet_pfakt5);
   fChain->SetBranchAddress("phiJet_pfakt5", phiJet_pfakt5, &b_phiJet_pfakt5);
   fChain->SetBranchAddress("nChargedHadrons_pfakt5", nChargedHadrons_pfakt5, &b_nChargedHadrons_pfakt5);
   fChain->SetBranchAddress("nPhotons_pfakt5", nPhotons_pfakt5, &b_nPhotons_pfakt5);
   fChain->SetBranchAddress("nMuons_pfakt5", nMuons_pfakt5, &b_nMuons_pfakt5);
   fChain->SetBranchAddress("nElectrons_pfakt5", nElectrons_pfakt5, &b_nElectrons_pfakt5);
   fChain->SetBranchAddress("nNeutralHadrons_pfakt5", nNeutralHadrons_pfakt5, &b_nNeutralHadrons_pfakt5);
   fChain->SetBranchAddress("nHFHadrons_pfakt5", nHFHadrons_pfakt5, &b_nHFHadrons_pfakt5);
   fChain->SetBranchAddress("nHFEM_pfakt5", nHFEM_pfakt5, &b_nHFEM_pfakt5);
   fChain->SetBranchAddress("eChargedHadrons_pfakt5", eChargedHadrons_pfakt5, &b_eChargedHadrons_pfakt5);
   fChain->SetBranchAddress("ePhotons_pfakt5", ePhotons_pfakt5, &b_ePhotons_pfakt5);
   fChain->SetBranchAddress("eMuons_pfakt5", eMuons_pfakt5, &b_eMuons_pfakt5);
   fChain->SetBranchAddress("eElectrons_pfakt5", eElectrons_pfakt5, &b_eElectrons_pfakt5);
   fChain->SetBranchAddress("eNeutralHadrons_pfakt5", eNeutralHadrons_pfakt5, &b_eNeutralHadrons_pfakt5);
   fChain->SetBranchAddress("eHFHadrons_pfakt5", eHFHadrons_pfakt5, &b_eHFHadrons_pfakt5);
   fChain->SetBranchAddress("eHFEM_pfakt5", eHFEM_pfakt5, &b_eHFEM_pfakt5);
   fChain->SetBranchAddress("nJet_pfakt7", &nJet_pfakt7, &b_nJet_pfakt7);
   fChain->SetBranchAddress("ptJet_pfakt7 ", ptJet_pfakt7 , &b_ptJet_pfakt7 );
   fChain->SetBranchAddress("ptCorrJet_pfakt7 ", ptCorrJet_pfakt7 , &b_ptCorrJet_pfakt7 );
   fChain->SetBranchAddress("eJet_pfakt7  ", eJet_pfakt7  , &b_eJet_pfakt7  );
   fChain->SetBranchAddress("etaJet_pfakt7", etaJet_pfakt7, &b_etaJet_pfakt7);
   fChain->SetBranchAddress("phiJet_pfakt7", phiJet_pfakt7, &b_phiJet_pfakt7);
//   fChain->SetBranchAddress("nJet_pfkt6", &nJet_pfkt6, &b_nJet_pfkt6);
//   fChain->SetBranchAddress("ptJet_pfkt6 ", ptJet_pfkt6 , &b_ptJet_pfkt6 );
//   fChain->SetBranchAddress("eJet_pfkt6  ", eJet_pfkt6  , &b_eJet_pfkt6  );
//   fChain->SetBranchAddress("etaJet_pfkt6", etaJet_pfkt6, &b_etaJet_pfkt6);
//   fChain->SetBranchAddress("phiJet_pfkt6", phiJet_pfkt6, &b_phiJet_pfkt6);
//   fChain->SetBranchAddress("nJetGen_akt5", &nJetGen_akt5, &b_nJetGen_akt5);
//   fChain->SetBranchAddress("ptJetGen_akt5 ", &ptJetGen_akt5 , &b_ptJetGen_akt5 );
//   fChain->SetBranchAddress("eJetGen_akt5  ", &eJetGen_akt5  , &b_eJetGen_akt5  );
//   fChain->SetBranchAddress("etaJetGen_akt5", &etaJetGen_akt5, &b_etaJetGen_akt5);
//   fChain->SetBranchAddress("phiJetGen_akt5", &phiJetGen_akt5, &b_phiJetGen_akt5);
//   fChain->SetBranchAddress("nMuonsGen_akt5", &nMuonsGen_akt5, &b_nMuonsGen_akt5);
//   fChain->SetBranchAddress("nElectronsGen_akt5", &nElectronsGen_akt5, &b_nElectronsGen_akt5);
//   fChain->SetBranchAddress("nPhotonsGen_akt5", &nPhotonsGen_akt5, &b_nPhotonsGen_akt5);
//   fChain->SetBranchAddress("nTracksGen_akt5", &nTracksGen_akt5, &b_nTracksGen_akt5);
//   fChain->SetBranchAddress("nNeutralHadronsGen_akt5", &nNeutralHadronsGen_akt5, &b_nNeutralHadronsGen_akt5);
//   fChain->SetBranchAddress("nHFHadronsGen_akt5", &nHFHadronsGen_akt5, &b_nHFHadronsGen_akt5);
//   fChain->SetBranchAddress("nHFEMGen_akt5", &nHFEMGen_akt5, &b_nHFEMGen_akt5);
//   fChain->SetBranchAddress("nNeutronsGen_akt5", &nNeutronsGen_akt5, &b_nNeutronsGen_akt5);
//   fChain->SetBranchAddress("nK0LGen_akt5", &nK0LGen_akt5, &b_nK0LGen_akt5);
//   fChain->SetBranchAddress("nK0SGen_akt5", &nK0SGen_akt5, &b_nK0SGen_akt5);
//   fChain->SetBranchAddress("nLambdasGen_akt5", &nLambdasGen_akt5, &b_nLambdasGen_akt5);
//   fChain->SetBranchAddress("nCsiGen_akt5", &nCsiGen_akt5, &b_nCsiGen_akt5);
//   fChain->SetBranchAddress("nOtherNeutralHadronsGen_akt5", &nOtherNeutralHadronsGen_akt5, &b_nOtherNeutralHadronsGen_akt5);
//   fChain->SetBranchAddress("eMuonsGen_akt5", &eMuonsGen_akt5, &b_eMuonsGen_akt5);
//   fChain->SetBranchAddress("eElectronsGen_akt5", &eElectronsGen_akt5, &b_eElectronsGen_akt5);
//   fChain->SetBranchAddress("ePhotonsGen_akt5", &ePhotonsGen_akt5, &b_ePhotonsGen_akt5);
//   fChain->SetBranchAddress("eTracksGen_akt5", &eTracksGen_akt5, &b_eTracksGen_akt5);
//   fChain->SetBranchAddress("eNeutralHadronsGen_akt5", &eNeutralHadronsGen_akt5, &b_eNeutralHadronsGen_akt5);
//   fChain->SetBranchAddress("eHFHadronsGen_akt5", &eHFHadronsGen_akt5, &b_eHFHadronsGen_akt5);
//   fChain->SetBranchAddress("eHFEMGen_akt5", &eHFEMGen_akt5, &b_eHFEMGen_akt5);
//   fChain->SetBranchAddress("eNeutronsGen_akt5", &eNeutronsGen_akt5, &b_eNeutronsGen_akt5);
//   fChain->SetBranchAddress("eK0LGen_akt5", &eK0LGen_akt5, &b_eK0LGen_akt5);
//   fChain->SetBranchAddress("eK0SGen_akt5", &eK0SGen_akt5, &b_eK0SGen_akt5);
//   fChain->SetBranchAddress("eLambdasGen_akt5", &eLambdasGen_akt5, &b_eLambdasGen_akt5);
//   fChain->SetBranchAddress("eCsiGen_akt5", &eCsiGen_akt5, &b_eCsiGen_akt5);
//   fChain->SetBranchAddress("eOtherNeutralHadronsGen_akt5", &eOtherNeutralHadronsGen_akt5, &b_eOtherNeutralHadronsGen_akt5);
//   fChain->SetBranchAddress("nJetGen_akt7", &nJetGen_akt7, &b_nJetGen_akt7);
//   fChain->SetBranchAddress("ptJetGen_akt7 ", &ptJetGen_akt7 , &b_ptJetGen_akt7 );
//   fChain->SetBranchAddress("eJetGen_akt7  ", &eJetGen_akt7  , &b_eJetGen_akt7  );
//   fChain->SetBranchAddress("etaJetGen_akt7", &etaJetGen_akt7, &b_etaJetGen_akt7);
//   fChain->SetBranchAddress("phiJetGen_akt7", &phiJetGen_akt7, &b_phiJetGen_akt7);
//   fChain->SetBranchAddress("nJetGen_kt4", &nJetGen_kt4, &b_nJetGen_kt4);
//   fChain->SetBranchAddress("ptJetGen_kt4 ", &ptJetGen_kt4 , &b_ptJetGen_kt4 );
//   fChain->SetBranchAddress("eJetGen_kt4  ", &eJetGen_kt4  , &b_eJetGen_kt4  );
//   fChain->SetBranchAddress("etaJetGen_kt4", &etaJetGen_kt4, &b_etaJetGen_kt4);
//   fChain->SetBranchAddress("phiJetGen_kt4", &phiJetGen_kt4, &b_phiJetGen_kt4);
//   fChain->SetBranchAddress("nJetGen_kt6", &nJetGen_kt6, &b_nJetGen_kt6);
//   fChain->SetBranchAddress("ptJetGen_kt6 ", &ptJetGen_kt6 , &b_ptJetGen_kt6 );
//   fChain->SetBranchAddress("eJetGen_kt6  ", &eJetGen_kt6  , &b_eJetGen_kt6  );
//   fChain->SetBranchAddress("etaJetGen_kt6", &etaJetGen_kt6, &b_etaJetGen_kt6);
//   fChain->SetBranchAddress("phiJetGen_kt6", &phiJetGen_kt6, &b_phiJetGen_kt6);
   fChain->SetBranchAddress("sMet  ", &sMet  , &b_sMet);
   fChain->SetBranchAddress("eMet  ", &eMet  , &b_eMet);
   fChain->SetBranchAddress("phiMet", &phiMet, &b_phiMet);
   fChain->SetBranchAddress("signifMet", &signifMet, &b_signifMet);
   fChain->SetBranchAddress("sCorrMet  ", &sCorrMet  , &b_sCorrMet);
   fChain->SetBranchAddress("eCorrMet  ", &eCorrMet  , &b_eCorrMet);
   fChain->SetBranchAddress("phiCorrMet", &phiCorrMet, &b_phiCorrMet);
   fChain->SetBranchAddress("signifCorrMet", &signifCorrMet, &b_signifCorrMet);
   fChain->SetBranchAddress("smuCorrMet  ", &smuCorrMet  , &b_smuCorrMet);
   fChain->SetBranchAddress("emuCorrMet  ", &emuCorrMet  , &b_emuCorrMet);
   fChain->SetBranchAddress("phimuCorrMet", &phimuCorrMet, &b_phimuCorrMet);
   fChain->SetBranchAddress("signifmuCorrMet", &signifmuCorrMet, &b_signifmuCorrMet);
   fChain->SetBranchAddress("sNoHFMet  ", &sNoHFMet  , &b_sNoHFMet);
   fChain->SetBranchAddress("eNoHFMet  ", &eNoHFMet  , &b_eNoHFMet);
   fChain->SetBranchAddress("phiNoHFMet", &phiNoHFMet, &b_phiNoHFMet);
   fChain->SetBranchAddress("signifNoHFMet", &signifNoHFMet, &b_signifNoHFMet);
   fChain->SetBranchAddress("stcMet  ", &stcMet  , &b_stcMet);
   fChain->SetBranchAddress("etcMet  ", &etcMet  , &b_etcMet);
   fChain->SetBranchAddress("phitcMet", &phitcMet, &b_phitcMet);
   fChain->SetBranchAddress("signiftcMet", &signiftcMet, &b_signiftcMet);
   fChain->SetBranchAddress("spfMet  ", &spfMet  , &b_spfMet);
   fChain->SetBranchAddress("epfMet  ", &epfMet  , &b_epfMet);
   fChain->SetBranchAddress("phipfMet", &phipfMet, &b_phipfMet);
   fChain->SetBranchAddress("signifpfMet", &signifpfMet, &b_signifpfMet);
   fChain->SetBranchAddress("sMetGen  ", &sMetGen  , &b_sMetGen);
   fChain->SetBranchAddress("eMetGen  ", &eMetGen  , &b_eMetGen);
   fChain->SetBranchAddress("phiMetGen", &phiMetGen, &b_phiMetGen);
   fChain->SetBranchAddress("signifMetGen", &signifMetGen, &b_signifMetGen);
   fChain->SetBranchAddress("sMetGen2  ", &sMetGen2  , &b_sMetGen2);
   fChain->SetBranchAddress("eMetGen2  ", &eMetGen2  , &b_eMetGen2);
   fChain->SetBranchAddress("phiMetGen2", &phiMetGen2, &b_phiMetGen2);
   fChain->SetBranchAddress("nvertex", &nvertex, &b_nvertex);
   fChain->SetBranchAddress("vxMC", &vxMC, &b_vxMC);
   fChain->SetBranchAddress("vyMC", &vyMC, &b_vyMC);
   fChain->SetBranchAddress("vzMC", &vzMC, &b_vzMC);
   fChain->SetBranchAddress("vx", vx, &b_vx);
   fChain->SetBranchAddress("vy", vy, &b_vy);
   fChain->SetBranchAddress("vz", vz, &b_vz);
   fChain->SetBranchAddress("vntracks", vntracks, &b_vntracks);
   fChain->SetBranchAddress("vchi2", vchi2, &b_vchi2);
   fChain->SetBranchAddress("vndof", vndof, &b_vndof);
   //fChain->SetBranchAddress("nHLT", &nHLT, &b_nHLT);
   //fChain->SetBranchAddress("hltNamesLen", &hltNamesLen, &b_hltNamesLen);
   //fChain->SetBranchAddress("HLTNames", &HLTNames, &b_HLTNames);
   //fChain->SetBranchAddress("HLTResults", &HLTResults, &b_HLTResults);
   fChain->SetBranchAddress("Xsec", &Xsec, &b_Xsec);
   Notify();
}

Bool_t tree_reader_V2::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void tree_reader_V2::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t tree_reader_V2::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
