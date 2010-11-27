#include "tree_reader_V1.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>


tree_reader_V1::tree_reader_V1(TTree *tree)
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("/cmshome/fanellic/CMSSW_3_6_3/src/rootfiles/crab_VBF/output_1_1_NnK.root");
      if (!f) {
         f = new TFile("/cmshome/fanellic/CMSSW_3_6_3/src/rootfiles/crab_VBF/output_1_1_NnK.root");
         f->cd("/cmshome/fanellic/CMSSW_3_6_3/src/rootfiles/crab_VBF/output_1_1_NnK.root:/myanalysis");
      }
      tree = (TTree*)gDirectory->Get("pippo");

   }
   Init(tree);
}

tree_reader_V1::~tree_reader_V1()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t tree_reader_V1::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t tree_reader_V1::LoadTree(Long64_t entry)
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

void tree_reader_V1::Init(TTree *tree)
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

   fChain->SetBranchAddress("genpt", &genpt, &b_genpt);
   fChain->SetBranchAddress("isMC", &isMC, &b_isMC);
   fChain->SetBranchAddress("store", &store, &b_store);
   fChain->SetBranchAddress("lbn", &lbn, &b_lbn);
   fChain->SetBranchAddress("bx", &bx, &b_bx);
   fChain->SetBranchAddress("orbit", &orbit, &b_orbit);
   fChain->SetBranchAddress("run", &run, &b_run);
   fChain->SetBranchAddress("event", &event, &b_event);
   fChain->SetBranchAddress("nMC", &nMC, &b_nMC);
   fChain->SetBranchAddress("pdgIdMC", pdgIdMC, &b_pdgIdMC);
   fChain->SetBranchAddress("statusMC", statusMC, &b_statusMC);
   fChain->SetBranchAddress("motherIDMC", motherIDMC, &b_motherIDMC);
   fChain->SetBranchAddress("ptMC ", ptMC , &b_ptMC );
   fChain->SetBranchAddress("eMC  ", eMC  , &b_eMC  );
   fChain->SetBranchAddress("etaMC", etaMC, &b_etaMC);
   fChain->SetBranchAddress("phiMC", phiMC, &b_phiMC);
   fChain->SetBranchAddress("nSIM", &nSIM, &b_nSIM);
   fChain->SetBranchAddress("pdgIdSIM", pdgIdSIM, &b_pdgIdSIM);
   fChain->SetBranchAddress("statusSIM", statusSIM, &b_statusSIM);
   fChain->SetBranchAddress("ptSIM ", ptSIM , &b_ptSIM );
   fChain->SetBranchAddress("eSIM  ", eSIM  , &b_eSIM  );
   fChain->SetBranchAddress("etaSIM", etaSIM, &b_etaSIM);
   fChain->SetBranchAddress("phiSIM", phiSIM, &b_phiSIM);
   fChain->SetBranchAddress("rSIM", rSIM, &b_rSIM);
   fChain->SetBranchAddress("zSIM", zSIM, &b_zSIM);
   fChain->SetBranchAddress("nPF", &nPF, &b_nPF);
   fChain->SetBranchAddress("pdgIdPF", pdgIdPF, &b_pdgIdPF);
   fChain->SetBranchAddress("ptPF ", ptPF , &b_ptPF );
   fChain->SetBranchAddress("ePF  ", ePF  , &b_ePF  );
   fChain->SetBranchAddress("etaPF", etaPF, &b_etaPF);
   fChain->SetBranchAddress("phiPF", phiPF, &b_phiPF);
   fChain->SetBranchAddress("nPhot", &nPhot, &b_nPhot);
   fChain->SetBranchAddress("ptPhot ", ptPhot , &b_ptPhot );
   fChain->SetBranchAddress("ePhot  ", ePhot  , &b_ePhot  );
   fChain->SetBranchAddress("escPhot  ", escPhot  , &b_escPhot  );
   fChain->SetBranchAddress("eseedPhot  ", eseedPhot  , &b_eseedPhot  );
   fChain->SetBranchAddress("etaPhot", etaPhot, &b_etaPhot);
   fChain->SetBranchAddress("phiPhot", phiPhot, &b_phiPhot);
   fChain->SetBranchAddress("timePhot", timePhot, &b_timePhot);
   fChain->SetBranchAddress("e4SwissCrossPhot", e4SwissCrossPhot, &b_e4SwissCrossPhot);
   fChain->SetBranchAddress("nconvPhot", &nconvPhot, &b_nconvPhot);
   fChain->SetBranchAddress("chi2convPhot", chi2convPhot, &b_chi2convPhot);
   fChain->SetBranchAddress("ndofconvPhot", ndofconvPhot, &b_ndofconvPhot);
   fChain->SetBranchAddress("rconvPhot", rconvPhot, &b_rconvPhot);
   fChain->SetBranchAddress("phiconvPhot", phiconvPhot, &b_phiconvPhot);
   fChain->SetBranchAddress("zconvPhot", zconvPhot, &b_zconvPhot);
   fChain->SetBranchAddress("ntrkconvPhot", ntrkconvPhot, &b_ntrkconvPhot);
   fChain->SetBranchAddress("eovpconvPhot", eovpconvPhot, &b_eovpconvPhot);
   fChain->SetBranchAddress("etaecalconvPhot", etaecalconvPhot, &b_etaecalconvPhot);
   fChain->SetBranchAddress("phiecalconvPhot", phiecalconvPhot, &b_phiecalconvPhot);
   fChain->SetBranchAddress("eecalconvPhot", eecalconvPhot, &b_eecalconvPhot);
   fChain->SetBranchAddress("algoconvPhot", algoconvPhot, &b_algoconvPhot);
   fChain->SetBranchAddress("d0convPhot", d0convPhot, &b_d0convPhot);
   fChain->SetBranchAddress("detaecalconvPhot", detaecalconvPhot, &b_detaecalconvPhot);
   fChain->SetBranchAddress("dphiecalconvPhot", dphiecalconvPhot, &b_dphiecalconvPhot);
   fChain->SetBranchAddress("dphivtxconvPhot", dphivtxconvPhot, &b_dphivtxconvPhot);
   fChain->SetBranchAddress("pairsepconvPhot", pairsepconvPhot, &b_pairsepconvPhot);
   fChain->SetBranchAddress("pairmassconvPhot", pairmassconvPhot, &b_pairmassconvPhot);
   fChain->SetBranchAddress("trchi21convPhot", trchi21convPhot, &b_trchi21convPhot);
   fChain->SetBranchAddress("trndof1convPhot", trndof1convPhot, &b_trndof1convPhot);
   fChain->SetBranchAddress("trqual1convPhot", trqual1convPhot, &b_trqual1convPhot);
   fChain->SetBranchAddress("trpt1convPhot", trpt1convPhot, &b_trpt1convPhot);
   fChain->SetBranchAddress("trerr1convPhot", trerr1convPhot, &b_trerr1convPhot);
   fChain->SetBranchAddress("trchi22convPhot", trchi22convPhot, &b_trchi22convPhot);
   fChain->SetBranchAddress("trndof2convPhot", trndof2convPhot, &b_trndof2convPhot);
   fChain->SetBranchAddress("trqual2convPhot", trqual2convPhot, &b_trqual2convPhot);
   fChain->SetBranchAddress("trpt2convPhot", trpt2convPhot, &b_trpt2convPhot);
   fChain->SetBranchAddress("trerr2convPhot", trerr2convPhot, &b_trerr2convPhot);
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
   fChain->SetBranchAddress("ptiso05Phot", ptiso05Phot, &b_ptiso05Phot);
   fChain->SetBranchAddress("ntrkiso05Phot", ntrkiso05Phot, &b_ntrkiso05Phot);
   fChain->SetBranchAddress("ptiso07Phot", ptiso07Phot, &b_ptiso07Phot);
   fChain->SetBranchAddress("ntrkiso07Phot", ntrkiso07Phot, &b_ntrkiso07Phot);
   fChain->SetBranchAddress("ptiso1Phot", ptiso1Phot, &b_ptiso1Phot);
   fChain->SetBranchAddress("ntrkiso1Phot", ntrkiso1Phot, &b_ntrkiso1Phot);
   fChain->SetBranchAddress("hcalovecal01Phot", hcalovecal01Phot, &b_hcalovecal01Phot);
   fChain->SetBranchAddress("hcalovecal015Phot", hcalovecal015Phot, &b_hcalovecal015Phot);
   fChain->SetBranchAddress("hcalovecal04Phot", hcalovecal04Phot, &b_hcalovecal04Phot);
   fChain->SetBranchAddress("hcalovecal05Phot", hcalovecal05Phot, &b_hcalovecal05Phot);
   fChain->SetBranchAddress("hcalovecal07Phot", hcalovecal07Phot, &b_hcalovecal07Phot);
   fChain->SetBranchAddress("hcalovecal1Phot", hcalovecal1Phot, &b_hcalovecal1Phot);
   fChain->SetBranchAddress("ecaliso01Phot", ecaliso01Phot, &b_ecaliso01Phot);
   fChain->SetBranchAddress("ecaliso015Phot", ecaliso015Phot, &b_ecaliso015Phot);
   fChain->SetBranchAddress("ecaliso04Phot", ecaliso04Phot, &b_ecaliso04Phot);
   fChain->SetBranchAddress("ecaliso05Phot", ecaliso05Phot, &b_ecaliso05Phot);
   fChain->SetBranchAddress("ecaliso07Phot", ecaliso07Phot, &b_ecaliso07Phot);
   fChain->SetBranchAddress("ecaliso1Phot", ecaliso1Phot, &b_ecaliso1Phot);
   fChain->SetBranchAddress("LATPhot", LATPhot, &b_LATPhot);
   fChain->SetBranchAddress("sMajMajPhot", sMajMajPhot, &b_sMajMajPhot);
   fChain->SetBranchAddress("sMinMinPhot", sMinMinPhot, &b_sMinMinPhot);
   fChain->SetBranchAddress("alphaPhot", alphaPhot, &b_alphaPhot);
   fChain->SetBranchAddress("sEtaEtaPhot", sEtaEtaPhot, &b_sEtaEtaPhot);
   fChain->SetBranchAddress("sEtaPhiPhot", sEtaPhiPhot, &b_sEtaPhiPhot);
   fChain->SetBranchAddress("sPhiPhiPhot", sPhiPhiPhot, &b_sPhiPhiPhot);
   fChain->SetBranchAddress("E1Phot", E1Phot, &b_E1Phot);
   fChain->SetBranchAddress("E9Phot", E9Phot, &b_E9Phot);
   fChain->SetBranchAddress("E25Phot", E25Phot, &b_E25Phot);
   fChain->SetBranchAddress("FisherPhot", FisherPhot, &b_FisherPhot);
   fChain->SetBranchAddress("nJet_ite", &nJet_ite, &b_nJet_ite);
   fChain->SetBranchAddress("ptJet_ite ", ptJet_ite , &b_ptJet_ite );
   fChain->SetBranchAddress("eJet_ite  ", eJet_ite  , &b_eJet_ite  );
   fChain->SetBranchAddress("etaJet_ite", etaJet_ite, &b_etaJet_ite);
   fChain->SetBranchAddress("phiJet_ite", phiJet_ite, &b_phiJet_ite);
   fChain->SetBranchAddress("emfJet_ite", emfJet_ite, &b_emfJet_ite);
   fChain->SetBranchAddress("nJet_kt4", &nJet_kt4, &b_nJet_kt4);
   fChain->SetBranchAddress("ptJet_kt4 ", ptJet_kt4 , &b_ptJet_kt4 );
   fChain->SetBranchAddress("eJet_kt4  ", eJet_kt4  , &b_eJet_kt4  );
   fChain->SetBranchAddress("etaJet_kt4", etaJet_kt4, &b_etaJet_kt4);
   fChain->SetBranchAddress("phiJet_kt4", phiJet_kt4, &b_phiJet_kt4);
   fChain->SetBranchAddress("emfJet_kt4", emfJet_kt4, &b_emfJet_kt4);
   fChain->SetBranchAddress("nJet_kt6", &nJet_kt6, &b_nJet_kt6);
   fChain->SetBranchAddress("ptJet_kt6 ", ptJet_kt6 , &b_ptJet_kt6 );
   fChain->SetBranchAddress("eJet_kt6  ", eJet_kt6  , &b_eJet_kt6  );
   fChain->SetBranchAddress("etaJet_kt6", etaJet_kt6, &b_etaJet_kt6);
   fChain->SetBranchAddress("phiJet_kt6", phiJet_kt6, &b_phiJet_kt6);
   fChain->SetBranchAddress("emfJet_kt6", emfJet_kt6, &b_emfJet_kt6);
   fChain->SetBranchAddress("nJet_akt5", &nJet_akt5, &b_nJet_akt5);
   fChain->SetBranchAddress("ptJet_akt5 ", ptJet_akt5 , &b_ptJet_akt5 );
   fChain->SetBranchAddress("ptCorrJet_akt5 ", ptCorrJet_akt5 , &b_ptCorrJet_akt5 );
   fChain->SetBranchAddress("eJet_akt5  ", eJet_akt5  , &b_eJet_akt5  );
   fChain->SetBranchAddress("etaJet_akt5", etaJet_akt5, &b_etaJet_akt5);
   fChain->SetBranchAddress("phiJet_akt5", phiJet_akt5, &b_phiJet_akt5);
   fChain->SetBranchAddress("emfJet_akt5", emfJet_akt5, &b_emfJet_akt5);
   fChain->SetBranchAddress("nJet_sis5", &nJet_sis5, &b_nJet_sis5);
   fChain->SetBranchAddress("ptJet_sis5 ", &ptJet_sis5 , &b_ptJet_sis5 );
   fChain->SetBranchAddress("eJet_sis5  ", &eJet_sis5  , &b_eJet_sis5  );
   fChain->SetBranchAddress("etaJet_sis5", &etaJet_sis5, &b_etaJet_sis5);
   fChain->SetBranchAddress("phiJet_sis5", &phiJet_sis5, &b_phiJet_sis5);
   fChain->SetBranchAddress("emfJet_sis5", &emfJet_sis5, &b_emfJet_sis5);
   fChain->SetBranchAddress("nJet_sis7", &nJet_sis7, &b_nJet_sis7);
   fChain->SetBranchAddress("ptJet_sis7 ", &ptJet_sis7 , &b_ptJet_sis7 );
   fChain->SetBranchAddress("eJet_sis7  ", &eJet_sis7  , &b_eJet_sis7  );
   fChain->SetBranchAddress("etaJet_sis7", &etaJet_sis7, &b_etaJet_sis7);
   fChain->SetBranchAddress("phiJet_sis7", &phiJet_sis7, &b_phiJet_sis7);
   fChain->SetBranchAddress("emfJet_sis7", &emfJet_sis7, &b_emfJet_sis7);
   fChain->SetBranchAddress("nJet_jptak5", &nJet_jptak5, &b_nJet_jptak5);
   fChain->SetBranchAddress("ptJet_jptak5 ", ptJet_jptak5 , &b_ptJet_jptak5 );
   fChain->SetBranchAddress("eJet_jptak5  ", eJet_jptak5  , &b_eJet_jptak5  );
   fChain->SetBranchAddress("etaJet_jptak5", etaJet_jptak5, &b_etaJet_jptak5);
   fChain->SetBranchAddress("phiJet_jptak5", phiJet_jptak5, &b_phiJet_jptak5);
   fChain->SetBranchAddress("emfJet_jptak5", emfJet_jptak5, &b_emfJet_jptak5);
   fChain->SetBranchAddress("nJet_pfite", &nJet_pfite, &b_nJet_pfite);
   fChain->SetBranchAddress("ptJet_pfite ", ptJet_pfite , &b_ptJet_pfite );
   fChain->SetBranchAddress("eJet_pfite  ", eJet_pfite  , &b_eJet_pfite  );
   fChain->SetBranchAddress("etaJet_pfite", etaJet_pfite, &b_etaJet_pfite);
   fChain->SetBranchAddress("phiJet_pfite", phiJet_pfite, &b_phiJet_pfite);
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
   fChain->SetBranchAddress("ptChargedHadrons_pfakt5", ptChargedHadrons_pfakt5, &b_ptChargedHadrons_pfakt5);
   fChain->SetBranchAddress("ptPhotons_pfakt5", ptPhotons_pfakt5, &b_ptPhotons_pfakt5);
   fChain->SetBranchAddress("ptMuons_pfakt5", ptMuons_pfakt5, &b_ptMuons_pfakt5);
   fChain->SetBranchAddress("ptElectrons_pfakt5", ptElectrons_pfakt5, &b_ptElectrons_pfakt5);
   fChain->SetBranchAddress("ptNeutralHadrons_pfakt5", ptNeutralHadrons_pfakt5, &b_ptNeutralHadrons_pfakt5);
   fChain->SetBranchAddress("ptHFHadrons_pfakt5", ptHFHadrons_pfakt5, &b_ptHFHadrons_pfakt5);
   fChain->SetBranchAddress("ptHFEM_pfakt5", ptHFEM_pfakt5, &b_ptHFEM_pfakt5);
   fChain->SetBranchAddress("phiChargedHadrons_pfakt5", phiChargedHadrons_pfakt5, &b_phiChargedHadrons_pfakt5);
   fChain->SetBranchAddress("phiPhotons_pfakt5", phiPhotons_pfakt5, &b_phiPhotons_pfakt5);
   fChain->SetBranchAddress("phiMuons_pfakt5", phiMuons_pfakt5, &b_phiMuons_pfakt5);
   fChain->SetBranchAddress("phiElectrons_pfakt5", phiElectrons_pfakt5, &b_phiElectrons_pfakt5);
   fChain->SetBranchAddress("phiNeutralHadrons_pfakt5", phiNeutralHadrons_pfakt5, &b_phiNeutralHadrons_pfakt5);
   fChain->SetBranchAddress("phiHFHadrons_pfakt5", phiHFHadrons_pfakt5, &b_phiHFHadrons_pfakt5);
   fChain->SetBranchAddress("phiHFEM_pfakt5", phiHFEM_pfakt5, &b_phiHFEM_pfakt5);
   fChain->SetBranchAddress("etaChargedHadrons_pfakt5", etaChargedHadrons_pfakt5, &b_etaChargedHadrons_pfakt5);
   fChain->SetBranchAddress("etaPhotons_pfakt5", etaPhotons_pfakt5, &b_etaPhotons_pfakt5);
   fChain->SetBranchAddress("etaMuons_pfakt5", etaMuons_pfakt5, &b_etaMuons_pfakt5);
   fChain->SetBranchAddress("etaElectrons_pfakt5", etaElectrons_pfakt5, &b_etaElectrons_pfakt5);
   fChain->SetBranchAddress("etaNeutralHadrons_pfakt5", etaNeutralHadrons_pfakt5, &b_etaNeutralHadrons_pfakt5);
   fChain->SetBranchAddress("etaHFHadrons_pfakt5", etaHFHadrons_pfakt5, &b_etaHFHadrons_pfakt5);
   fChain->SetBranchAddress("etaHFEM_pfakt5", etaHFEM_pfakt5, &b_etaHFEM_pfakt5);
   fChain->SetBranchAddress("nJet_pfakt7", &nJet_pfakt7, &b_nJet_pfakt7);
   fChain->SetBranchAddress("ptJet_pfakt7 ", ptJet_pfakt7 , &b_ptJet_pfakt7 );
   fChain->SetBranchAddress("ptCorrJet_pfakt7 ", ptCorrJet_pfakt7 , &b_ptCorrJet_pfakt7 );
   fChain->SetBranchAddress("eJet_pfakt7  ", eJet_pfakt7  , &b_eJet_pfakt7  );
   fChain->SetBranchAddress("etaJet_pfakt7", etaJet_pfakt7, &b_etaJet_pfakt7);
   fChain->SetBranchAddress("phiJet_pfakt7", phiJet_pfakt7, &b_phiJet_pfakt7);
   fChain->SetBranchAddress("nChargedHadrons_pfakt7", nChargedHadrons_pfakt7, &b_nChargedHadrons_pfakt7);
   fChain->SetBranchAddress("nPhotons_pfakt7", nPhotons_pfakt7, &b_nPhotons_pfakt7);
   fChain->SetBranchAddress("nMuons_pfakt7", nMuons_pfakt7, &b_nMuons_pfakt7);
   fChain->SetBranchAddress("nElectrons_pfakt7", nElectrons_pfakt7, &b_nElectrons_pfakt7);
   fChain->SetBranchAddress("nNeutralHadrons_pfakt7", nNeutralHadrons_pfakt7, &b_nNeutralHadrons_pfakt7);
   fChain->SetBranchAddress("nHFHadrons_pfakt7", nHFHadrons_pfakt7, &b_nHFHadrons_pfakt7);
   fChain->SetBranchAddress("nHFEM_pfakt7", nHFEM_pfakt7, &b_nHFEM_pfakt7);
   fChain->SetBranchAddress("eChargedHadrons_pfakt7", eChargedHadrons_pfakt7, &b_eChargedHadrons_pfakt7);
   fChain->SetBranchAddress("ePhotons_pfakt7", ePhotons_pfakt7, &b_ePhotons_pfakt7);
   fChain->SetBranchAddress("eMuons_pfakt7", eMuons_pfakt7, &b_eMuons_pfakt7);
   fChain->SetBranchAddress("eElectrons_pfakt7", eElectrons_pfakt7, &b_eElectrons_pfakt7);
   fChain->SetBranchAddress("eNeutralHadrons_pfakt7", eNeutralHadrons_pfakt7, &b_eNeutralHadrons_pfakt7);
   fChain->SetBranchAddress("eHFHadrons_pfakt7", eHFHadrons_pfakt7, &b_eHFHadrons_pfakt7);
   fChain->SetBranchAddress("eHFEM_pfakt7", eHFEM_pfakt7, &b_eHFEM_pfakt7);
   fChain->SetBranchAddress("ptChargedHadrons_pfakt7", ptChargedHadrons_pfakt7, &b_ptChargedHadrons_pfakt7);
   fChain->SetBranchAddress("ptPhotons_pfakt7", ptPhotons_pfakt7, &b_ptPhotons_pfakt7);
   fChain->SetBranchAddress("ptMuons_pfakt7", ptMuons_pfakt7, &b_ptMuons_pfakt7);
   fChain->SetBranchAddress("ptElectrons_pfakt7", ptElectrons_pfakt7, &b_ptElectrons_pfakt7);
   fChain->SetBranchAddress("ptNeutralHadrons_pfakt7", ptNeutralHadrons_pfakt7, &b_ptNeutralHadrons_pfakt7);
   fChain->SetBranchAddress("ptHFHadrons_pfakt7", ptHFHadrons_pfakt7, &b_ptHFHadrons_pfakt7);
   fChain->SetBranchAddress("ptHFEM_pfakt7", ptHFEM_pfakt7, &b_ptHFEM_pfakt7);
   fChain->SetBranchAddress("phiChargedHadrons_pfakt7", phiChargedHadrons_pfakt7, &b_phiChargedHadrons_pfakt7);
   fChain->SetBranchAddress("phiPhotons_pfakt7", phiPhotons_pfakt7, &b_phiPhotons_pfakt7);
   fChain->SetBranchAddress("phiMuons_pfakt7", phiMuons_pfakt7, &b_phiMuons_pfakt7);
   fChain->SetBranchAddress("phiElectrons_pfakt7", phiElectrons_pfakt7, &b_phiElectrons_pfakt7);
   fChain->SetBranchAddress("phiNeutralHadrons_pfakt7", phiNeutralHadrons_pfakt7, &b_phiNeutralHadrons_pfakt7);
   fChain->SetBranchAddress("phiHFHadrons_pfakt7", phiHFHadrons_pfakt7, &b_phiHFHadrons_pfakt7);
   fChain->SetBranchAddress("phiHFEM_pfakt7", phiHFEM_pfakt7, &b_phiHFEM_pfakt7);
   fChain->SetBranchAddress("etaChargedHadrons_pfakt7", etaChargedHadrons_pfakt7, &b_etaChargedHadrons_pfakt7);
   fChain->SetBranchAddress("etaPhotons_pfakt7", etaPhotons_pfakt7, &b_etaPhotons_pfakt7);
   fChain->SetBranchAddress("etaMuons_pfakt7", etaMuons_pfakt7, &b_etaMuons_pfakt7);
   fChain->SetBranchAddress("etaElectrons_pfakt7", etaElectrons_pfakt7, &b_etaElectrons_pfakt7);
   fChain->SetBranchAddress("etaNeutralHadrons_pfakt7", etaNeutralHadrons_pfakt7, &b_etaNeutralHadrons_pfakt7);
   fChain->SetBranchAddress("etaHFHadrons_pfakt7", etaHFHadrons_pfakt7, &b_etaHFHadrons_pfakt7);
   fChain->SetBranchAddress("etaHFEM_pfakt7", etaHFEM_pfakt7, &b_etaHFEM_pfakt7);
   fChain->SetBranchAddress("nJet_pfsis5", &nJet_pfsis5, &b_nJet_pfsis5);
   fChain->SetBranchAddress("ptJet_pfsis5 ", &ptJet_pfsis5 , &b_ptJet_pfsis5 );
   fChain->SetBranchAddress("eJet_pfsis5  ", &eJet_pfsis5  , &b_eJet_pfsis5  );
   fChain->SetBranchAddress("etaJet_pfsis5", &etaJet_pfsis5, &b_etaJet_pfsis5);
   fChain->SetBranchAddress("phiJet_pfsis5", &phiJet_pfsis5, &b_phiJet_pfsis5);
   fChain->SetBranchAddress("nJet_pfkt6", &nJet_pfkt6, &b_nJet_pfkt6);
   fChain->SetBranchAddress("ptJet_pfkt6 ", ptJet_pfkt6 , &b_ptJet_pfkt6 );
   fChain->SetBranchAddress("eJet_pfkt6  ", eJet_pfkt6  , &b_eJet_pfkt6  );
   fChain->SetBranchAddress("etaJet_pfkt6", etaJet_pfkt6, &b_etaJet_pfkt6);
   fChain->SetBranchAddress("phiJet_pfkt6", phiJet_pfkt6, &b_phiJet_pfkt6);
   fChain->SetBranchAddress("nJet_pfsis7", &nJet_pfsis7, &b_nJet_pfsis7);
   fChain->SetBranchAddress("ptJet_pfsis7 ", &ptJet_pfsis7 , &b_ptJet_pfsis7 );
   fChain->SetBranchAddress("eJet_pfsis7  ", &eJet_pfsis7  , &b_eJet_pfsis7  );
   fChain->SetBranchAddress("etaJet_pfsis7", &etaJet_pfsis7, &b_etaJet_pfsis7);
   fChain->SetBranchAddress("phiJet_pfsis7", &phiJet_pfsis7, &b_phiJet_pfsis7);
   fChain->SetBranchAddress("nJetGen_ite", &nJetGen_ite, &b_nJetGen_ite);
   fChain->SetBranchAddress("ptJetGen_ite ", ptJetGen_ite , &b_ptJetGen_ite );
   fChain->SetBranchAddress("eJetGen_ite  ", eJetGen_ite  , &b_eJetGen_ite  );
   fChain->SetBranchAddress("etaJetGen_ite", etaJetGen_ite, &b_etaJetGen_ite);
   fChain->SetBranchAddress("phiJetGen_ite", phiJetGen_ite, &b_phiJetGen_ite);
   fChain->SetBranchAddress("nJetGen_akt5", &nJetGen_akt5, &b_nJetGen_akt5);
   fChain->SetBranchAddress("ptJetGen_akt5 ", ptJetGen_akt5 , &b_ptJetGen_akt5 );
   fChain->SetBranchAddress("eJetGen_akt5  ", eJetGen_akt5  , &b_eJetGen_akt5  );
   fChain->SetBranchAddress("etaJetGen_akt5", etaJetGen_akt5, &b_etaJetGen_akt5);
   fChain->SetBranchAddress("phiJetGen_akt5", phiJetGen_akt5, &b_phiJetGen_akt5);
   fChain->SetBranchAddress("nMuonsGen_akt5", nMuonsGen_akt5, &b_nMuonsGen_akt5);
   fChain->SetBranchAddress("nElectronsGen_akt5", nElectronsGen_akt5, &b_nElectronsGen_akt5);
   fChain->SetBranchAddress("nPhotonsGen_akt5", nPhotonsGen_akt5, &b_nPhotonsGen_akt5);
   fChain->SetBranchAddress("nTracksGen_akt5", nTracksGen_akt5, &b_nTracksGen_akt5);
   fChain->SetBranchAddress("nNeutralHadronsGen_akt5", nNeutralHadronsGen_akt5, &b_nNeutralHadronsGen_akt5);
   fChain->SetBranchAddress("nHFHadronsGen_akt5", nHFHadronsGen_akt5, &b_nHFHadronsGen_akt5);
   fChain->SetBranchAddress("nHFEMGen_akt5", nHFEMGen_akt5, &b_nHFEMGen_akt5);
   fChain->SetBranchAddress("nNeutronsGen_akt5", nNeutronsGen_akt5, &b_nNeutronsGen_akt5);
   fChain->SetBranchAddress("nK0LGen_akt5", nK0LGen_akt5, &b_nK0LGen_akt5);
   fChain->SetBranchAddress("nK0SGen_akt5", nK0SGen_akt5, &b_nK0SGen_akt5);
   fChain->SetBranchAddress("nLambdasGen_akt5", nLambdasGen_akt5, &b_nLambdasGen_akt5);
   fChain->SetBranchAddress("nCsiGen_akt5", nCsiGen_akt5, &b_nCsiGen_akt5);
   fChain->SetBranchAddress("nOtherNeutralHadronsGen_akt5", nOtherNeutralHadronsGen_akt5, &b_nOtherNeutralHadronsGen_akt5);
   fChain->SetBranchAddress("eMuonsGen_akt5", eMuonsGen_akt5, &b_eMuonsGen_akt5);
   fChain->SetBranchAddress("eElectronsGen_akt5", eElectronsGen_akt5, &b_eElectronsGen_akt5);
   fChain->SetBranchAddress("ePhotonsGen_akt5", ePhotonsGen_akt5, &b_ePhotonsGen_akt5);
   fChain->SetBranchAddress("eTracksGen_akt5", eTracksGen_akt5, &b_eTracksGen_akt5);
   fChain->SetBranchAddress("eNeutralHadronsGen_akt5", eNeutralHadronsGen_akt5, &b_eNeutralHadronsGen_akt5);
   fChain->SetBranchAddress("eHFHadronsGen_akt5", eHFHadronsGen_akt5, &b_eHFHadronsGen_akt5);
   fChain->SetBranchAddress("eHFEMGen_akt5", eHFEMGen_akt5, &b_eHFEMGen_akt5);
   fChain->SetBranchAddress("eNeutronsGen_akt5", eNeutronsGen_akt5, &b_eNeutronsGen_akt5);
   fChain->SetBranchAddress("eK0LGen_akt5", eK0LGen_akt5, &b_eK0LGen_akt5);
   fChain->SetBranchAddress("eK0SGen_akt5", eK0SGen_akt5, &b_eK0SGen_akt5);
   fChain->SetBranchAddress("eLambdasGen_akt5", eLambdasGen_akt5, &b_eLambdasGen_akt5);
   fChain->SetBranchAddress("eCsiGen_akt5", eCsiGen_akt5, &b_eCsiGen_akt5);
   fChain->SetBranchAddress("eOtherNeutralHadronsGen_akt5", eOtherNeutralHadronsGen_akt5, &b_eOtherNeutralHadronsGen_akt5);
   fChain->SetBranchAddress("ptMuonsGen_akt5", ptMuonsGen_akt5, &b_ptMuonsGen_akt5);
   fChain->SetBranchAddress("ptElectronsGen_akt5", ptElectronsGen_akt5, &b_ptElectronsGen_akt5);
   fChain->SetBranchAddress("ptPhotonsGen_akt5", ptPhotonsGen_akt5, &b_ptPhotonsGen_akt5);
   fChain->SetBranchAddress("ptTracksGen_akt5", ptTracksGen_akt5, &b_ptTracksGen_akt5);
   fChain->SetBranchAddress("ptNeutralHadronsGen_akt5", ptNeutralHadronsGen_akt5, &b_ptNeutralHadronsGen_akt5);
   fChain->SetBranchAddress("ptHFHadronsGen_akt5", ptHFHadronsGen_akt5, &b_ptHFHadronsGen_akt5);
   fChain->SetBranchAddress("ptHFEMGen_akt5", ptHFEMGen_akt5, &b_ptHFEMGen_akt5);
   fChain->SetBranchAddress("phiMuonsGen_akt5", phiMuonsGen_akt5, &b_phiMuonsGen_akt5);
   fChain->SetBranchAddress("phiElectronsGen_akt5", phiElectronsGen_akt5, &b_phiElectronsGen_akt5);
   fChain->SetBranchAddress("phiPhotonsGen_akt5", phiPhotonsGen_akt5, &b_phiPhotonsGen_akt5);
   fChain->SetBranchAddress("phiTracksGen_akt5", phiTracksGen_akt5, &b_phiTracksGen_akt5);
   fChain->SetBranchAddress("phiNeutralHadronsGen_akt5", phiNeutralHadronsGen_akt5, &b_phiNeutralHadronsGen_akt5);
   fChain->SetBranchAddress("phiHFHadronsGen_akt5", phiHFHadronsGen_akt5, &b_phiHFHadronsGen_akt5);
   fChain->SetBranchAddress("phiHFEMGen_akt5", phiHFEMGen_akt5, &b_phiHFEMGen_akt5);
   fChain->SetBranchAddress("etaMuonsGen_akt5", etaMuonsGen_akt5, &b_etaMuonsGen_akt5);
   fChain->SetBranchAddress("etaElectronsGen_akt5", etaElectronsGen_akt5, &b_etaElectronsGen_akt5);
   fChain->SetBranchAddress("etaPhotonsGen_akt5", etaPhotonsGen_akt5, &b_etaPhotonsGen_akt5);
   fChain->SetBranchAddress("etaTracksGen_akt5", etaTracksGen_akt5, &b_etaTracksGen_akt5);
   fChain->SetBranchAddress("etaNeutralHadronsGen_akt5", etaNeutralHadronsGen_akt5, &b_etaNeutralHadronsGen_akt5);
   fChain->SetBranchAddress("etaHFHadronsGen_akt5", etaHFHadronsGen_akt5, &b_etaHFHadronsGen_akt5);
   fChain->SetBranchAddress("etaHFEMGen_akt5", etaHFEMGen_akt5, &b_etaHFEMGen_akt5);
   fChain->SetBranchAddress("nJetGen_akt7", &nJetGen_akt7, &b_nJetGen_akt7);
   fChain->SetBranchAddress("ptJetGen_akt7 ", ptJetGen_akt7 , &b_ptJetGen_akt7 );
   fChain->SetBranchAddress("eJetGen_akt7  ", eJetGen_akt7  , &b_eJetGen_akt7  );
   fChain->SetBranchAddress("etaJetGen_akt7", etaJetGen_akt7, &b_etaJetGen_akt7);
   fChain->SetBranchAddress("phiJetGen_akt7", phiJetGen_akt7, &b_phiJetGen_akt7);
   fChain->SetBranchAddress("nJetGen_kt4", &nJetGen_kt4, &b_nJetGen_kt4);
   fChain->SetBranchAddress("ptJetGen_kt4 ", ptJetGen_kt4 , &b_ptJetGen_kt4 );
   fChain->SetBranchAddress("eJetGen_kt4  ", eJetGen_kt4  , &b_eJetGen_kt4  );
   fChain->SetBranchAddress("etaJetGen_kt4", etaJetGen_kt4, &b_etaJetGen_kt4);
   fChain->SetBranchAddress("phiJetGen_kt4", phiJetGen_kt4, &b_phiJetGen_kt4);
   fChain->SetBranchAddress("nJetGen_kt6", &nJetGen_kt6, &b_nJetGen_kt6);
   fChain->SetBranchAddress("ptJetGen_kt6 ", ptJetGen_kt6 , &b_ptJetGen_kt6 );
   fChain->SetBranchAddress("eJetGen_kt6  ", eJetGen_kt6  , &b_eJetGen_kt6  );
   fChain->SetBranchAddress("etaJetGen_kt6", etaJetGen_kt6, &b_etaJetGen_kt6);
   fChain->SetBranchAddress("phiJetGen_kt6", phiJetGen_kt6, &b_phiJetGen_kt6);
   fChain->SetBranchAddress("nJetGen_sis5", &nJetGen_sis5, &b_nJetGen_sis5);
   fChain->SetBranchAddress("ptJetGen_sis5", &ptJetGen_sis5, &b_ptJetGen_sis5);
   fChain->SetBranchAddress("eJetGen_sis5  ", &eJetGen_sis5  , &b_eJetGen_sis5  );
   fChain->SetBranchAddress("etaJetGen_sis5", &etaJetGen_sis5, &b_etaJetGen_sis5);
   fChain->SetBranchAddress("phiJetGen_sis5", &phiJetGen_sis5, &b_phiJetGen_sis5);
   fChain->SetBranchAddress("nJetGen_sis7", &nJetGen_sis7, &b_nJetGen_sis7);
   fChain->SetBranchAddress("ptJetGen_sis7 ", &ptJetGen_sis7 , &b_ptJetGen_sis7 );
   fChain->SetBranchAddress("eJetGen_sis7  ", &eJetGen_sis7  , &b_eJetGen_sis7  );
   fChain->SetBranchAddress("etaJetGen_sis7", &etaJetGen_sis7, &b_etaJetGen_sis7);
   fChain->SetBranchAddress("phiJetGen_sis7", &phiJetGen_sis7, &b_phiJetGen_sis7);
   fChain->SetBranchAddress("sMet  ", &sMet  , &b_sMet);
   fChain->SetBranchAddress("eMet  ", &eMet  , &b_eMet);
   fChain->SetBranchAddress("phiMet", &phiMet, &b_phiMet);
   fChain->SetBranchAddress("stcMet  ", &stcMet  , &b_stcMet);
   fChain->SetBranchAddress("etcMet  ", &etcMet  , &b_etcMet);
   fChain->SetBranchAddress("phitcMet", &phitcMet, &b_phitcMet);
   fChain->SetBranchAddress("spfMet  ", &spfMet  , &b_spfMet);
   fChain->SetBranchAddress("epfMet  ", &epfMet  , &b_epfMet);
   fChain->SetBranchAddress("phipfMet", &phipfMet, &b_phipfMet);
   fChain->SetBranchAddress("sMetGen  ", &sMetGen  , &b_sMetGen);
   fChain->SetBranchAddress("eMetGen  ", &eMetGen  , &b_eMetGen);
   fChain->SetBranchAddress("phiMetGen", &phiMetGen, &b_phiMetGen);
   fChain->SetBranchAddress("sMetGen2  ", &sMetGen2  , &b_sMetGen2);
   fChain->SetBranchAddress("eMetGen2  ", &eMetGen2  , &b_eMetGen2);
   fChain->SetBranchAddress("phiMetGen2", &phiMetGen2, &b_phiMetGen2);
   fChain->SetBranchAddress("nvertex", &nvertex, &b_nvertex);
   fChain->SetBranchAddress("vxMC", &vxMC, &b_vxMC);
   fChain->SetBranchAddress("vyMC", &vyMC, &b_vyMC);
   fChain->SetBranchAddress("vzMC", &vzMC, &b_vzMC);
   fChain->SetBranchAddress("vx", &vx, &b_vx);
   fChain->SetBranchAddress("vy", &vy, &b_vy);
   fChain->SetBranchAddress("vz", &vz, &b_vz);
   fChain->SetBranchAddress("vntracks", &vntracks, &b_vntracks);
   fChain->SetBranchAddress("vchi2", &vchi2, &b_vchi2);
   fChain->SetBranchAddress("vndof", &vndof, &b_vndof);
   fChain->SetBranchAddress("hltPass", &hltPass, &b_hltPass);
   fChain->SetBranchAddress("nHLT", &nHLT, &b_nHLT);
   fChain->SetBranchAddress("hltNamesLen", &hltNamesLen, &b_hltNamesLen);
   fChain->SetBranchAddress("HLTNames", &HLTNames, &b_HLTNames);
   fChain->SetBranchAddress("HLTResults", HLTResults, &b_HLTResults);
   Notify();
}

Bool_t tree_reader_V1::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void tree_reader_V1::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t tree_reader_V1::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}


void tree_reader_V1::Loop()
{
//   In a ROOT session, you can do:
//      Root > .L tree_reader_V1.C
//      Root > tree_reader_V1 t
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
