#define tree_reader_V8_cxx
#include "tree_reader_V8.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>

void tree_reader_V8::Loop()
{
//   In a ROOT session, you can do:
//      Root > .L tree_reader_V8.C
//      Root > tree_reader_V8 t
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

#ifdef tree_reader_V8_cxx
tree_reader_V8::tree_reader_V8(TTree *tree) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("root://eoscms//eos/cms/store/group/phys_higgs/meridian/MC/52xv1/GluGluToHToGG_M-125_8TeV-powheg-pythia6_Summer12-PU_S7_START52_V9-v1/output_1_1_3xh.root");
      if (!f || !f->IsOpen()) {
         f = new TFile("root://eoscms//eos/cms/store/group/phys_higgs/meridian/MC/52xv1/GluGluToHToGG_M-125_8TeV-powheg-pythia6_Summer12-PU_S7_START52_V9-v1/output_1_1_3xh.root");
      }
      TDirectory * dir = (TDirectory*)f->Get("root://eoscms//eos/cms/store/group/phys_higgs/meridian/MC/52xv1/GluGluToHToGG_M-125_8TeV-powheg-pythia6_Summer12-PU_S7_START52_V9-v1/output_1_1_3xh.root:/myanalysis");
      dir->GetObject("pippo",tree);

   }
   Init(tree);
}

tree_reader_V8::~tree_reader_V8()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t tree_reader_V8::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t tree_reader_V8::LoadTree(Long64_t entry)
{
// Set the environment to read one entry
   if (!fChain) return -5;
   Long64_t centry = fChain->LoadTree(entry);
   if (centry < 0) return centry;
   if (fChain->GetTreeNumber() != fCurrent) {
      fCurrent = fChain->GetTreeNumber();
      Notify();
   }
   return centry;
}

void tree_reader_V8::Init(TTree *tree)
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
   fChain->SetBranchAddress("genProcessId", &genProcessId, &b_genProcessId);
   fChain->SetBranchAddress("genQScale", &genQScale, &b_genQScale);
   fChain->SetBranchAddress("qPDF", &qPDF, &b_qPDF);
   fChain->SetBranchAddress("x1PDF", &x1PDF, &b_x1PDF);
   fChain->SetBranchAddress("x2PDF", &x2PDF, &b_x2PDF);
   fChain->SetBranchAddress("id1PDF", &id1PDF, &b_id1PDF);
   fChain->SetBranchAddress("id2PDF", &id2PDF, &b_id2PDF);
   fChain->SetBranchAddress("nWeightsPDF", nWeightsPDF, &b_nWeightsPDF);
   fChain->SetBranchAddress("pdfWeight", pdfWeight, &b_pdfWeight);
   fChain->SetBranchAddress("isMC", &isMC, &b_isMC);
   fChain->SetBranchAddress("store", &store, &b_store);
   fChain->SetBranchAddress("lbn", &lbn, &b_lbn);
   fChain->SetBranchAddress("bx", &bx, &b_bx);
   fChain->SetBranchAddress("orbit", &orbit, &b_orbit);
   fChain->SetBranchAddress("run", &run, &b_run);
   fChain->SetBranchAddress("event", &event, &b_event);
   fChain->SetBranchAddress("rhoPF", &rhoPF, &b_rhoPF);
   fChain->SetBranchAddress("rhoCalo", &rhoCalo, &b_rhoCalo);
   fChain->SetBranchAddress("rhoAllJets", &rhoAllJets, &b_rhoAllJets);
   fChain->SetBranchAddress("nMC", &nMC, &b_nMC);
   fChain->SetBranchAddress("pdgIdMC", pdgIdMC, &b_pdgIdMC);
   fChain->SetBranchAddress("statusMC", statusMC, &b_statusMC);
   fChain->SetBranchAddress("motherIDMC", motherIDMC, &b_motherIDMC);
   fChain->SetBranchAddress("ptMC ", ptMC , &b_ptMC );
   fChain->SetBranchAddress("eMC  ", eMC  , &b_eMC  );
   fChain->SetBranchAddress("etaMC", etaMC, &b_etaMC);
   fChain->SetBranchAddress("phiMC", phiMC, &b_phiMC);
   fChain->SetBranchAddress("pu_n", &pu_n, &b_pu_n);
   fChain->SetBranchAddress("pu_true_n", &pu_true_n, &b_pu_true_n);
   fChain->SetBranchAddress("pu_zpos", pu_zpos, &b_pu_zpos);
   fChain->SetBranchAddress("pu_sumpt_lowpt", pu_sumpt_lowpt, &b_pu_sumpt_lowpt);
   fChain->SetBranchAddress("pu_sumpt_highpt", pu_sumpt_highpt, &b_pu_sumpt_highpt);
   fChain->SetBranchAddress("pu_ntrks_lowpt", pu_ntrks_lowpt, &b_pu_ntrks_lowpt);
   fChain->SetBranchAddress("pu_ntrks_highpt", pu_ntrks_highpt, &b_pu_ntrks_highpt);
   fChain->SetBranchAddress("nPhot", &nPhot, &b_nPhot);
   fChain->SetBranchAddress("ptPhot ", ptPhot , &b_ptPhot );
   fChain->SetBranchAddress("ePhot  ", ePhot  , &b_ePhot  );
   fChain->SetBranchAddress("escPhot  ", escPhot  , &b_escPhot  );
   fChain->SetBranchAddress("escRegrPhot  ", escRegrPhot  , &b_escRegrPhot  );
   fChain->SetBranchAddress("escRegrPhotError  ", escRegrPhotError  , &b_escRegrPhotError  );
   fChain->SetBranchAddress("escPhFixPhot  ", escPhFixPhot  , &b_escPhFixPhot  );
   fChain->SetBranchAddress("escPhFixPhotError  ", escPhFixPhotError  , &b_escPhFixPhotError  );
   fChain->SetBranchAddress("escRawPhot  ", escRawPhot  , &b_escRawPhot  );
   fChain->SetBranchAddress("etascPhot  ", etascPhot  , &b_etascPhot  );
   fChain->SetBranchAddress("phiscPhot  ", phiscPhot  , &b_phiscPhot  );
   fChain->SetBranchAddress("xscPhot  ", xscPhot  , &b_xscPhot  );
   fChain->SetBranchAddress("yscPhot  ", yscPhot  , &b_yscPhot  );
   fChain->SetBranchAddress("zscPhot  ", zscPhot  , &b_zscPhot  );
   fChain->SetBranchAddress("xcaloPhot  ", xcaloPhot  , &b_xcaloPhot  );
   fChain->SetBranchAddress("ycaloPhot  ", ycaloPhot  , &b_ycaloPhot  );
   fChain->SetBranchAddress("zcaloPhot  ", zcaloPhot  , &b_zcaloPhot  );
   fChain->SetBranchAddress("eseedPhot  ", eseedPhot  , &b_eseedPhot  );
   fChain->SetBranchAddress("etaPhot", etaPhot, &b_etaPhot);
   fChain->SetBranchAddress("phiPhot", phiPhot, &b_phiPhot);
   fChain->SetBranchAddress("timePhot", timePhot, &b_timePhot);
   fChain->SetBranchAddress("e4SwissCrossPhot", e4SwissCrossPhot, &b_e4SwissCrossPhot);
   fChain->SetBranchAddress("hasPixelSeedPhot", hasPixelSeedPhot, &b_hasPixelSeedPhot);
   fChain->SetBranchAddress("hasMatchedPromptElePhot", hasMatchedPromptElePhot, &b_hasMatchedPromptElePhot);
   fChain->SetBranchAddress("hasMatchedConvPhot", hasMatchedConvPhot, &b_hasMatchedConvPhot);
   fChain->SetBranchAddress("isEBPhot", isEBPhot, &b_isEBPhot);
   fChain->SetBranchAddress("isEEPhot", isEEPhot, &b_isEEPhot);
   fChain->SetBranchAddress("isEBEEGapPhot", isEBEEGapPhot, &b_isEBEEGapPhot);
   fChain->SetBranchAddress("ntracksConvPhot", ntracksConvPhot, &b_ntracksConvPhot);
   fChain->SetBranchAddress("isValidVtxConvPhot", isValidVtxConvPhot, &b_isValidVtxConvPhot);
   fChain->SetBranchAddress("pairInvmassConvPhot", pairInvmassConvPhot, &b_pairInvmassConvPhot);
   fChain->SetBranchAddress("pairCotThetaSeperationConvPhot", pairCotThetaSeperationConvPhot, &b_pairCotThetaSeperationConvPhot);
   fChain->SetBranchAddress("pairmomentumXConvPhot", pairmomentumXConvPhot, &b_pairmomentumXConvPhot);
   fChain->SetBranchAddress("pairmomentumYConvPhot", pairmomentumYConvPhot, &b_pairmomentumYConvPhot);
   fChain->SetBranchAddress("pairmomentumZConvPhot", pairmomentumZConvPhot, &b_pairmomentumZConvPhot);
   fChain->SetBranchAddress("chi2ConvPhot", chi2ConvPhot, &b_chi2ConvPhot);
   fChain->SetBranchAddress("nDofConvPhot", nDofConvPhot, &b_nDofConvPhot);
   fChain->SetBranchAddress("eOverPConvPhot", eOverPConvPhot, &b_eOverPConvPhot);
   fChain->SetBranchAddress("convVxConvPhot", convVxConvPhot, &b_convVxConvPhot);
   fChain->SetBranchAddress("convVyConvPhot", convVyConvPhot, &b_convVyConvPhot);
   fChain->SetBranchAddress("convVzConvPhot", convVzConvPhot, &b_convVzConvPhot);
   fChain->SetBranchAddress("distOfMinimumApproachConvPhot", distOfMinimumApproachConvPhot, &b_distOfMinimumApproachConvPhot);
   fChain->SetBranchAddress("dPhiTracksAtVtxConvPhot", dPhiTracksAtVtxConvPhot, &b_dPhiTracksAtVtxConvPhot);
   fChain->SetBranchAddress("pid_isEM", pid_isEM, &b_pid_isEM);
   fChain->SetBranchAddress("pid_isLoose", pid_isLoose, &b_pid_isLoose);
   fChain->SetBranchAddress("pid_isTight", pid_isTight, &b_pid_isTight);
   fChain->SetBranchAddress("pid_jurECAL", pid_jurECAL, &b_pid_jurECAL);
   fChain->SetBranchAddress("pid_twrHCAL", pid_twrHCAL, &b_pid_twrHCAL);
   fChain->SetBranchAddress("pid_HoverE", pid_HoverE, &b_pid_HoverE);
   fChain->SetBranchAddress("pid_hlwTrack", pid_hlwTrack, &b_pid_hlwTrack);
   fChain->SetBranchAddress("pid_hlwTrackNoDz", pid_hlwTrackNoDz, &b_pid_hlwTrackNoDz);
   fChain->SetBranchAddress("pid_hlwTrackForCiC", pid_hlwTrackForCiC, &b_pid_hlwTrackBestRank);
   fChain->SetBranchAddress("pid_etawid", pid_etawid, &b_pid_etawid);
   fChain->SetBranchAddress("pid_jurECAL03", pid_jurECAL03, &b_pid_jurECAL03);
   fChain->SetBranchAddress("pid_twrHCAL03", pid_twrHCAL03, &b_pid_twrHCAL03);
   fChain->SetBranchAddress("pid_hlwTrack03", pid_hlwTrack03, &b_pid_hlwTrack03);
   fChain->SetBranchAddress("pid_hlwTrack03NoDz", pid_hlwTrack03NoDz, &b_pid_hlwTrack03NoDz);
   fChain->SetBranchAddress("pid_hlwTrack03ForCiC", pid_hlwTrack03ForCiC, &b_pid_hlwTrack03ForCiC);
   fChain->SetBranchAddress("pid_pfIsoCharged01ForCiC", pid_pfIsoCharged01ForCiC, &b_pid_pfIsoCharged01ForCiC);
   fChain->SetBranchAddress("pid_pfIsoCharged02ForCiC", pid_pfIsoCharged02ForCiC, &b_pid_pfIsoCharged02ForCiC);
   fChain->SetBranchAddress("pid_pfIsoCharged03ForCiC", pid_pfIsoCharged03ForCiC, &b_pid_pfIsoCharged03ForCiC);
   fChain->SetBranchAddress("pid_pfIsoCharged04ForCiC", pid_pfIsoCharged04ForCiC, &b_pid_pfIsoCharged04ForCiC);
   fChain->SetBranchAddress("pid_pfIsoCharged05ForCiC", pid_pfIsoCharged05ForCiC, &b_pid_pfIsoCharged05ForCiC);
   fChain->SetBranchAddress("pid_pfIsoCharged06ForCiC", pid_pfIsoCharged06ForCiC, &b_pid_pfIsoCharged06ForCiC);
   fChain->SetBranchAddress("pid_pfIsoPhotons01ForCiC", pid_pfIsoPhotons01ForCiC, &b_pid_pfIsoPhotons01ForCiC);
   fChain->SetBranchAddress("pid_pfIsoPhotons02ForCiC", pid_pfIsoPhotons02ForCiC, &b_pid_pfIsoPhotons02ForCiC);
   fChain->SetBranchAddress("pid_pfIsoPhotons03ForCiC", pid_pfIsoPhotons03ForCiC, &b_pid_pfIsoPhotons03ForCiC);
   fChain->SetBranchAddress("pid_pfIsoPhotons04ForCiC", pid_pfIsoPhotons04ForCiC, &b_pid_pfIsoPhotons04ForCiC);
   fChain->SetBranchAddress("pid_pfIsoPhotons05ForCiC", pid_pfIsoPhotons05ForCiC, &b_pid_pfIsoPhotons05ForCiC);
   fChain->SetBranchAddress("pid_pfIsoPhotons06ForCiC", pid_pfIsoPhotons06ForCiC, &b_pid_pfIsoPhotons06ForCiC);
   fChain->SetBranchAddress("pid_pfIsoNeutrals01ForCiC", pid_pfIsoNeutrals01ForCiC, &b_pid_pfIsoNeutrals01ForCiC);
   fChain->SetBranchAddress("pid_pfIsoNeutrals02ForCiC", pid_pfIsoNeutrals02ForCiC, &b_pid_pfIsoNeutrals02ForCiC);
   fChain->SetBranchAddress("pid_pfIsoNeutrals03ForCiC", pid_pfIsoNeutrals03ForCiC, &b_pid_pfIsoNeutrals03ForCiC);
   fChain->SetBranchAddress("pid_pfIsoNeutrals04ForCiC", pid_pfIsoNeutrals04ForCiC, &b_pid_pfIsoNeutrals04ForCiC);
   fChain->SetBranchAddress("pid_pfIsoNeutrals05ForCiC", pid_pfIsoNeutrals05ForCiC, &b_pid_pfIsoNeutrals05ForCiC);
   fChain->SetBranchAddress("pid_pfIsoNeutrals06ForCiC", pid_pfIsoNeutrals06ForCiC, &b_pid_pfIsoNeutrals06ForCiC);
   fChain->SetBranchAddress("ptiso004Phot", ptiso004Phot, &b_ptiso004Phot);
   fChain->SetBranchAddress("ntrkiso004Phot", ntrkiso004Phot, &b_ntrkiso004Phot);
   fChain->SetBranchAddress("ptiso035Phot", ptiso035Phot, &b_ptiso035Phot);
   fChain->SetBranchAddress("ntrkiso035Phot", ntrkiso035Phot, &b_ntrkiso035Phot);
   fChain->SetBranchAddress("ptiso04Phot", ptiso04Phot, &b_ptiso04Phot);
   fChain->SetBranchAddress("ntrkiso04Phot", ntrkiso04Phot, &b_ntrkiso04Phot);
   fChain->SetBranchAddress("hcalovecal04Phot", hcalovecal04Phot, &b_hcalovecal04Phot);
   fChain->SetBranchAddress("ecaliso04Phot", ecaliso04Phot, &b_ecaliso04Phot);
   fChain->SetBranchAddress("pid_scetawid", pid_scetawid, &b_pid_scetawid);
   fChain->SetBranchAddress("pid_scphiwid", pid_scphiwid, &b_pid_scphiwid);
   fChain->SetBranchAddress("pid_lambdaRatio", pid_lambdaRatio, &b_pid_lambdaRatio);
   fChain->SetBranchAddress("pid_esXwidth", pid_esXwidth, &b_pid_esXwidth);
   fChain->SetBranchAddress("pid_esYwidth", pid_esYwidth, &b_pid_esYwidth);
   fChain->SetBranchAddress("sMajMajPhot", sMajMajPhot, &b_sMajMajPhot);
   fChain->SetBranchAddress("sMinMinPhot", sMinMinPhot, &b_sMinMinPhot);
   fChain->SetBranchAddress("alphaPhot", alphaPhot, &b_alphaPhot);
   fChain->SetBranchAddress("sEtaEtaPhot", sEtaEtaPhot, &b_sEtaEtaPhot);
   fChain->SetBranchAddress("sEtaPhiPhot", sEtaPhiPhot, &b_sEtaPhiPhot);
   fChain->SetBranchAddress("sPhiPhiPhot", sPhiPhiPhot, &b_sPhiPhiPhot);
   fChain->SetBranchAddress("E1Phot", E1Phot, &b_E1Phot);
   fChain->SetBranchAddress("E2OverE9Phot", E2OverE9Phot, &b_E2OverE9Phot);
   fChain->SetBranchAddress("E4Phot", E4Phot, &b_E4Phot);
   fChain->SetBranchAddress("E9Phot", E9Phot, &b_E9Phot);
   fChain->SetBranchAddress("E25Phot", E25Phot, &b_E25Phot);
   fChain->SetBranchAddress("ieleassocPhot", ieleassocPhot, &b_ieleassocPhot);
   fChain->SetBranchAddress("pid_deltaRToTrackPhot", pid_deltaRToTrackPhot, &b_pid_deltaRToTrackPhot);
   fChain->SetBranchAddress("nElePhot", &nElePhot, &b_nElePhot);
   fChain->SetBranchAddress("pid_jurECALElePhot ", pid_jurECALElePhot , &b_pid_jurECALElePhot );
   fChain->SetBranchAddress("pid_twrHCALElePhot ", pid_twrHCALElePhot , &b_pid_twrHCALElePhot );
   fChain->SetBranchAddress("pid_HoverEElePhot ", pid_HoverEElePhot , &b_pid_HoverEElePhot );
   fChain->SetBranchAddress("pid_hlwTrackElePhot ", pid_hlwTrackElePhot , &b_pid_hlwTrackElePhot );
   fChain->SetBranchAddress("pid_etawidElePhot ", pid_etawidElePhot , &b_pid_etawidElePhot );
   fChain->SetBranchAddress("pid_dphivtxElePhot ", pid_dphivtxElePhot , &b_pid_dphivtxElePhot );
   fChain->SetBranchAddress("pid_detavtxElePhot ", pid_detavtxElePhot , &b_pid_detavtxElePhot );
   fChain->SetBranchAddress("pid_mishitsElePhot ", pid_mishitsElePhot , &b_pid_mishitsElePhot );
   fChain->SetBranchAddress("pid_distElePhot ", pid_distElePhot , &b_pid_distElePhot );
   fChain->SetBranchAddress("pid_dcotElePhot ", pid_dcotElePhot , &b_pid_dcotElePhot );
   fChain->SetBranchAddress("pid_ptElePhot ", pid_ptElePhot , &b_pid_ptElePhot );
   fChain->SetBranchAddress("nJet_akt5", &nJet_akt5, &b_nJet_akt5);
   fChain->SetBranchAddress("ptJet_akt5 ", ptJet_akt5 , &b_ptJet_akt5 );
   fChain->SetBranchAddress("ptCorrJet_akt5 ", ptCorrJet_akt5 , &b_ptCorrJet_akt5 );
   fChain->SetBranchAddress("eJet_akt5  ", eJet_akt5  , &b_eJet_akt5  );
   fChain->SetBranchAddress("etaJet_akt5", etaJet_akt5, &b_etaJet_akt5);
   fChain->SetBranchAddress("phiJet_akt5", phiJet_akt5, &b_phiJet_akt5);
   fChain->SetBranchAddress("emfJet_akt5", emfJet_akt5, &b_emfJet_akt5);
   fChain->SetBranchAddress("n90Jet_akt5", n90Jet_akt5, &b_n90Jet_akt5);
   fChain->SetBranchAddress("n90HitsJet_akt5", n90HitsJet_akt5, &b_n90HitsJet_akt5);
   fChain->SetBranchAddress("fHPDJet_akt5", fHPDJet_akt5, &b_fHPDJet_akt5);
   fChain->SetBranchAddress("fRBXJet_akt5", fRBXJet_akt5, &b_fRBXJet_akt5);
   fChain->SetBranchAddress("nJet_akt7", &nJet_akt7, &b_nJet_akt7);
   fChain->SetBranchAddress("ptJet_akt7 ", ptJet_akt7 , &b_ptJet_akt7 );
   fChain->SetBranchAddress("ptCorrJet_akt7 ", ptCorrJet_akt7 , &b_ptCorrJet_akt7 );
   fChain->SetBranchAddress("eJet_akt7  ", eJet_akt7  , &b_eJet_akt7  );
   fChain->SetBranchAddress("etaJet_akt7", etaJet_akt7, &b_etaJet_akt7);
   fChain->SetBranchAddress("phiJet_akt7", phiJet_akt7, &b_phiJet_akt7);
   fChain->SetBranchAddress("emfJet_akt7", emfJet_akt7, &b_emfJet_akt7);
   fChain->SetBranchAddress("n90Jet_akt7", n90Jet_akt7, &b_n90Jet_akt7);
   fChain->SetBranchAddress("n90HitsJet_akt7", n90HitsJet_akt7, &b_n90HitsJet_akt7);
   fChain->SetBranchAddress("fHPDJet_akt7", fHPDJet_akt7, &b_fHPDJet_akt7);
   fChain->SetBranchAddress("fRBXJet_akt7", fRBXJet_akt7, &b_fRBXJet_akt7);
   fChain->SetBranchAddress("nJet_pfakt5", &nJet_pfakt5, &b_nJet_pfakt5);
   fChain->SetBranchAddress("ptJet_pfakt5 ", ptJet_pfakt5 , &b_ptJet_pfakt5 );
   fChain->SetBranchAddress("ptCorrJet_pfakt5 ", ptCorrJet_pfakt5 , &b_ptCorrJet_pfakt5 );
   fChain->SetBranchAddress("eJet_pfakt5  ", eJet_pfakt5  , &b_eJet_pfakt5  );
   fChain->SetBranchAddress("etaJet_pfakt5", etaJet_pfakt5, &b_etaJet_pfakt5);
   fChain->SetBranchAddress("phiJet_pfakt5", phiJet_pfakt5, &b_phiJet_pfakt5);
   fChain->SetBranchAddress("ptDJet_pfakt5", ptDJet_pfakt5, &b_ptDJet_pfakt5);
   fChain->SetBranchAddress("ptD_QCJet_pfakt5", ptD_QCJet_pfakt5, &b_ptD_QCJet_pfakt5);
   fChain->SetBranchAddress("axis2_QCJet_pfakt5", axis2_QCJet_pfakt5, &b_axis2_QCJet_pfakt5);
   fChain->SetBranchAddress("rmsCandJet_pfakt5", rmsCandJet_pfakt5, &b_rmsCandJet_pfakt5);
   fChain->SetBranchAddress("jetId_dRMean_pfakt5", jetId_dRMean_pfakt5, &b_jetId_dRMean_pfakt5);
   fChain->SetBranchAddress("jetId_frac01_pfakt5", jetId_frac01_pfakt5, &b_jetId_frac01_pfakt5);
   fChain->SetBranchAddress("jetId_frac02_pfakt5", jetId_frac02_pfakt5, &b_jetId_frac02_pfakt5);
   fChain->SetBranchAddress("jetId_frac03_pfakt5", jetId_frac03_pfakt5, &b_jetId_frac03_pfakt5);
   fChain->SetBranchAddress("jetId_frac04_pfakt5", jetId_frac04_pfakt5, &b_jetId_frac04_pfakt5);
   fChain->SetBranchAddress("jetId_frac05_pfakt5", jetId_frac05_pfakt5, &b_jetId_frac05_pfakt5);
   fChain->SetBranchAddress("jetId_nNeutrals_pfakt5", jetId_nNeutrals_pfakt5, &b_jetId_nNeutrals_pfakt5);
   fChain->SetBranchAddress("jetId_beta_pfakt5", jetId_beta_pfakt5, &b_jetId_beta_pfakt5);
   fChain->SetBranchAddress("jetId_betaStar_pfakt5", jetId_betaStar_pfakt5, &b_jetId_betaStar_pfakt5);
   fChain->SetBranchAddress("jetId_dZ_pfakt5", jetId_dZ_pfakt5, &b_jetId_dZ_pfakt5);
   fChain->SetBranchAddress("jetId_nCharged_pfakt5", jetId_nCharged_pfakt5, &b_jetId_nCharged_pfakt5);
   fChain->SetBranchAddress("jetId_dR2Mean_pfakt5", jetId_dR2Mean_pfakt5, &b_jetId_dR2Mean_pfakt5);
   fChain->SetBranchAddress("jetId_betaStarClassic_pfakt5", jetId_betaStarClassic_pfakt5, &b_jetId_betaStarClassic_pfakt5);
   fChain->SetBranchAddress("jetIdSimple_mva_pfakt5", jetIdSimple_mva_pfakt5, &b_jetIdSimple_mva_pfakt5);
   fChain->SetBranchAddress("jetIdSimple_wp_pfakt5", jetIdSimple_wp_pfakt5, &b_jetIdSimple_wp_pfakt5);
   fChain->SetBranchAddress("jetIdFull_mva_pfakt5", jetIdFull_mva_pfakt5, &b_jetIdFull_mva_pfakt5);
   fChain->SetBranchAddress("jetIdFull_wp_pfakt5", jetIdFull_wp_pfakt5, &b_jetIdFull_wp_pfakt5);
   fChain->SetBranchAddress("jetIdCutBased_mva_pfakt5", jetIdCutBased_mva_pfakt5, &b_jetIdCutBased_mva_pfakt5);
   fChain->SetBranchAddress("jetIdCutBased_wp_pfakt5", jetIdCutBased_wp_pfakt5, &b_jetIdCutBased_wp_pfakt5);
   fChain->SetBranchAddress("beta_pfakt5", beta_pfakt5, &b_beta_pfakt5);
   fChain->SetBranchAddress("betaStar_pfakt5", betaStar_pfakt5, &b_betaStar_pfakt5);

   fChain->SetBranchAddress("nChargedHadronsgoodvtx_pfakt5", nChargedHadronsgoodvtx_pfakt5, &b_nChargedHadronsgoodvtx_pfakt5);
   fChain->SetBranchAddress("eChargedHadronsgoodvtx_pfakt5", eChargedHadronsgoodvtx_pfakt5, &b_eChargedHadronsgoodvtx_pfakt5);
   fChain->SetBranchAddress("ptChargedHadronsgoodvtx_pfakt5", ptChargedHadronsgoodvtx_pfakt5, &b_ptChargedHadronsgoodvtx_pfakt5);
   fChain->SetBranchAddress("ptChargedHadrons_pfakt5", ptChargedHadrons_pfakt5, &b_ptChargedHadrons_pfakt5);
   fChain->SetBranchAddress("ptPhotons_pfakt5", ptPhotons_pfakt5, &b_ptPhotons_pfakt5);
   fChain->SetBranchAddress("ptMuons_pfakt5", ptMuons_pfakt5, &b_ptMuons_pfakt5);
   fChain->SetBranchAddress("ptElectrons_pfakt5", ptElectrons_pfakt5, &b_ptElectrons_pfakt5);
   fChain->SetBranchAddress("ptNeutralHadrons_pfakt5", ptNeutralHadrons_pfakt5, &b_ptNeutralHadrons_pfakt5);
   fChain->SetBranchAddress("ptHFHadrons_pfakt5", ptHFHadrons_pfakt5, &b_ptHFHadrons_pfakt5);
   fChain->SetBranchAddress("ptHFEM_pfakt5", ptHFEM_pfakt5, &b_ptHFEM_pfakt5);
   fChain->SetBranchAddress("etaChargedHadronsgoodvtx_pfakt5", etaChargedHadronsgoodvtx_pfakt5, &b_etaChargedHadronsgoodvtx_pfakt5);
   fChain->SetBranchAddress("etaChargedHadrons_pfakt5", etaChargedHadrons_pfakt5, &b_etaChargedHadrons_pfakt5);
   fChain->SetBranchAddress("etaPhotons_pfakt5", etaPhotons_pfakt5, &b_etaPhotons_pfakt5);
   fChain->SetBranchAddress("etaMuons_pfakt5", etaMuons_pfakt5, &b_etaMuons_pfakt5);
   fChain->SetBranchAddress("etaElectrons_pfakt5", etaElectrons_pfakt5, &b_etaElectrons_pfakt5);
   fChain->SetBranchAddress("etaNeutralHadrons_pfakt5", etaNeutralHadrons_pfakt5, &b_etaNeutralHadrons_pfakt5);
   fChain->SetBranchAddress("etaHFHadrons_pfakt5", etaHFHadrons_pfakt5, &b_etaHFHadrons_pfakt5);
   fChain->SetBranchAddress("etaHFEM_pfakt5", etaHFEM_pfakt5, &b_etaHFEM_pfakt5);
   fChain->SetBranchAddress("phiChargedHadrons_pfakt5", phiChargedHadrons_pfakt5, &b_phiChargedHadrons_pfakt5);
   fChain->SetBranchAddress("phiChargedHadronsgoodvtx_pfakt5", phiChargedHadronsgoodvtx_pfakt5, &b_phiChargedHadronsgoodvtx_pfakt5);
   fChain->SetBranchAddress("phiPhotons_pfakt5", phiPhotons_pfakt5, &b_phiPhotons_pfakt5);
   fChain->SetBranchAddress("phiMuons_pfakt5", phiMuons_pfakt5, &b_phiMuons_pfakt5);
   fChain->SetBranchAddress("phiElectrons_pfakt5", phiElectrons_pfakt5, &b_phiElectrons_pfakt5);
   fChain->SetBranchAddress("phiNeutralHadrons_pfakt5", phiNeutralHadrons_pfakt5, &b_phiNeutralHadrons_pfakt5);
   fChain->SetBranchAddress("phiHFHadrons_pfakt5", phiHFHadrons_pfakt5, &b_phiHFHadrons_pfakt5);
   fChain->SetBranchAddress("phiHFEM_pfakt5", phiHFEM_pfakt5, &b_phiHFEM_pfakt5);
   fChain->SetBranchAddress("sumptChargedHadronsgoodvtx_pfakt5", sumptChargedHadronsgoodvtx_pfakt5, &b_sumptChargedHadronsgoodvtx_pfakt5);
   fChain->SetBranchAddress("sumptChargedHadrons_pfakt5", sumptChargedHadrons_pfakt5, &b_sumptChargedHadrons_pfakt5);
   fChain->SetBranchAddress("sumptPhotons_pfakt5", sumptPhotons_pfakt5, &b_sumptPhotons_pfakt5);
   fChain->SetBranchAddress("sumptMuons_pfakt5", sumptMuons_pfakt5, &b_sumptMuons_pfakt5);
   fChain->SetBranchAddress("sumptElectrons_pfakt5", sumptElectrons_pfakt5, &b_sumptElectrons_pfakt5);
   fChain->SetBranchAddress("sumptNeutralHadrons_pfakt5", sumptNeutralHadrons_pfakt5, &b_sumptNeutralHadrons_pfakt5);
   fChain->SetBranchAddress("sumptHFHadrons_pfakt5", sumptHFHadrons_pfakt5, &b_sumptHFHadrons_pfakt5);
   fChain->SetBranchAddress("sumptHFEM_pfakt5", sumptHFEM_pfakt5, &b_sumptHFEM_pfakt5);

   fChain->SetBranchAddress("npfcand_all", &npfcand_all, &b_npfcand_all);
   fChain->SetBranchAddress("nChargedHadrons_uncl", &nChargedHadrons_uncl, &b_nChargedHadrons_uncl);
   fChain->SetBranchAddress("nChargedHadronsgoodvtx_uncl", &nChargedHadronsgoodvtx_uncl, &b_nChargedHadronsgoodvtx_uncl);
   fChain->SetBranchAddress("nPhotons_uncl", &nPhotons_uncl, &b_nPhotons_uncl);
   fChain->SetBranchAddress("nMuons_uncl", &nMuons_uncl, &b_nMuons_uncl);
   fChain->SetBranchAddress("nElectrons_uncl", &nElectrons_uncl, &b_nElectrons_uncl);
   fChain->SetBranchAddress("nNeutralHadrons_uncl", &nNeutralHadrons_uncl, &b_nNeutralHadrons_uncl);
   fChain->SetBranchAddress("nHFHadrons_uncl", &nHFHadrons_uncl, &b_nHFHadrons_uncl);
   fChain->SetBranchAddress("nHFEM_uncl", &nHFEM_uncl, &b_nHFEM_uncl);
   fChain->SetBranchAddress("epfcand_all", &epfcand_all, &b_epfcand_all);
   fChain->SetBranchAddress("eChargedHadrons_uncl", &eChargedHadrons_uncl, &b_eChargedHadrons_uncl);
   fChain->SetBranchAddress("eChargedHadronsgoodvtx_uncl", &eChargedHadronsgoodvtx_uncl, &b_eChargedHadronsgoodvtx_uncl);
   fChain->SetBranchAddress("ePhotons_uncl", &ePhotons_uncl, &b_ePhotons_uncl);

   fChain->SetBranchAddress("eMuons_uncl", &eMuons_uncl, &b_eMuons_uncl);
   fChain->SetBranchAddress("eElectrons_uncl", &eElectrons_uncl, &b_eElectrons_uncl);
   fChain->SetBranchAddress("eNeutralHadrons_uncl", &eNeutralHadrons_uncl, &b_eNeutralHadrons_uncl);
   fChain->SetBranchAddress("eHFHadrons_uncl", &eHFHadrons_uncl, &b_eHFHadrons_uncl);
   fChain->SetBranchAddress("eHFEM_uncl", &eHFEM_uncl, &b_eHFEM_uncl);
   fChain->SetBranchAddress("ptpfcand_all", &ptpfcand_all, &b_ptpfcand_all);
   fChain->SetBranchAddress("ptChargedHadrons_uncl", &ptChargedHadrons_uncl, &b_ptChargedHadrons_uncl);
   fChain->SetBranchAddress("ptChargedHadronsgoodvtx_uncl", &ptChargedHadronsgoodvtx_uncl, &b_ptChargedHadronsgoodvtx_uncl);
   fChain->SetBranchAddress("ptPhotons_uncl", &ptPhotons_uncl, &b_ptPhotons_uncl);
   fChain->SetBranchAddress("ptMuons_uncl", &ptMuons_uncl, &b_ptMuons_uncl);
   fChain->SetBranchAddress("ptElectrons_uncl", &ptElectrons_uncl, &b_ptElectrons_uncl);
   fChain->SetBranchAddress("ptNeutralHadrons_uncl", &ptNeutralHadrons_uncl, &b_ptNeutralHadrons_uncl);
   fChain->SetBranchAddress("ptHFHadrons_uncl", &ptHFHadrons_uncl, &b_ptHFHadrons_uncl);
   fChain->SetBranchAddress("ptHFEM_uncl", &ptHFEM_uncl, &b_ptHFEM_uncl);
   fChain->SetBranchAddress("etapfcand_all", &etapfcand_all, &b_etapfcand_all);
   fChain->SetBranchAddress("etaChargedHadrons_uncl", &etaChargedHadrons_uncl, &b_etaChargedHadrons_uncl);
   fChain->SetBranchAddress("etaChargedHadronsgoodvtx_uncl", &etaChargedHadronsgoodvtx_uncl, &b_etaChargedHadronsgoodvtx_uncl);
   fChain->SetBranchAddress("etaPhotons_uncl", &etaPhotons_uncl, &b_etaPhotons_uncl);
   fChain->SetBranchAddress("etaMuons_uncl", &etaMuons_uncl, &b_etaMuons_uncl);
   fChain->SetBranchAddress("etaElectrons_uncl", &etaElectrons_uncl, &b_etaElectrons_uncl);
   fChain->SetBranchAddress("etaNeutralHadrons_uncl", &etaNeutralHadrons_uncl, &b_etaNeutralHadrons_uncl);
   fChain->SetBranchAddress("etaHFHadrons_uncl", &etaHFHadrons_uncl, &b_etaHFHadrons_uncl);
   fChain->SetBranchAddress("etaHFEM_uncl", &etaHFEM_uncl, &b_etaHFEM_uncl);
   fChain->SetBranchAddress("phipfcand_all", &phipfcand_all, &b_phipfcand_all);
   fChain->SetBranchAddress("phiChargedHadrons_uncl", &phiChargedHadrons_uncl, &b_phiChargedHadrons_uncl);
   fChain->SetBranchAddress("phiChargedHadronsgoodvtx_uncl", &phiChargedHadronsgoodvtx_uncl, &b_phiChargedHadronsgoodvtx_uncl);
   fChain->SetBranchAddress("phiPhotons_uncl", &phiPhotons_uncl, &b_phiPhotons_uncl);
   fChain->SetBranchAddress("phiMuons_uncl", &phiMuons_uncl, &b_phiMuons_uncl);
   fChain->SetBranchAddress("phiElectrons_uncl", &phiElectrons_uncl, &b_phiElectrons_uncl);
   fChain->SetBranchAddress("phiNeutralHadrons_uncl", &phiNeutralHadrons_uncl, &b_phiNeutralHadrons_uncl);
   fChain->SetBranchAddress("phiHFHadrons_uncl", &phiHFHadrons_uncl, &b_phiHFHadrons_uncl);
   fChain->SetBranchAddress("phiHFEM_uncl", &phiHFEM_uncl, &b_phiHFEM_uncl);
   fChain->SetBranchAddress("sumptpfcand_all", &sumptpfcand_all, &b_sumptpfcand_all);
   fChain->SetBranchAddress("sumptChargedHadrons_uncl", &sumptChargedHadrons_uncl, &b_sumptChargedHadrons_uncl);
   fChain->SetBranchAddress("sumptChargedHadronsgoodvtx_uncl", &sumptChargedHadronsgoodvtx_uncl, &b_sumptChargedHadronsgoodvtx_uncl);
   fChain->SetBranchAddress("sumptPhotons_uncl", &sumptPhotons_uncl, &b_sumptPhotons_uncl);
   fChain->SetBranchAddress("sumptMuons_uncl", &sumptMuons_uncl, &b_sumptMuons_uncl);
   fChain->SetBranchAddress("sumptElectrons_uncl", &sumptElectrons_uncl, &b_sumptElectrons_uncl);
   fChain->SetBranchAddress("sumptNeutralHadrons_uncl", &sumptNeutralHadrons_uncl, &b_sumptNeutralHadrons_uncl);
   fChain->SetBranchAddress("sumptHFHadrons_uncl", &sumptHFHadrons_uncl, &b_sumptHFHadrons_uncl);
   fChain->SetBranchAddress("sumptHFEM_uncl", &sumptHFEM_uncl, &b_sumptHFEM_uncl);

   fChain->SetBranchAddress("combinedSecondaryVertexBJetTags", combinedSecondaryVertexBJetTags, &b_combinedSecondaryVertexBJetTags);
   fChain->SetBranchAddress("combinedSecondaryVertexMVABJetTags", combinedSecondaryVertexMVABJetTags, &b_combinedSecondaryVertexMVABJetTags);
   fChain->SetBranchAddress("jetBProbabilityBJetTags", jetBProbabilityBJetTags, &b_jetBProbabilityBJetTags);
   fChain->SetBranchAddress("jetProbabilityBJetTags", jetProbabilityBJetTags, &b_jetProbabilityBJetTags);
   fChain->SetBranchAddress("simpleSecondaryVertexHighEffBJetTags", simpleSecondaryVertexHighEffBJetTags, &b_simpleSecondaryVertexHighEffBJetTags);
   fChain->SetBranchAddress("simpleSecondaryVertexHighPurBJetTags", simpleSecondaryVertexHighPurBJetTags, &b_simpleSecondaryVertexHighPurBJetTags);
   fChain->SetBranchAddress("softMuonBJetTags", softMuonBJetTags, &b_softMuonBJetTags);
   fChain->SetBranchAddress("softMuonByIP3dBJetTags", softMuonByIP3dBJetTags, &b_softMuonByIP3dBJetTags);
   fChain->SetBranchAddress("softMuonByPtBJetTags", softMuonByPtBJetTags, &b_softMuonByPtBJetTags);
   fChain->SetBranchAddress("softElectronByIP3dBJetTags", softElectronByIP3dBJetTags, &b_softElectronByIP3dBJetTags);
   fChain->SetBranchAddress("softElectronByPtBJetTags", softElectronByPtBJetTags, &b_softElectronByPtBJetTags);
   fChain->SetBranchAddress("trackCountingHighPurBJetTags", trackCountingHighPurBJetTags, &b_trackCountingHighPurBJetTags);
   fChain->SetBranchAddress("trackCountingHighEffBJetTags", trackCountingHighEffBJetTags, &b_trackCountingHighEffBJetTags);
   fChain->SetBranchAddress("nChg_QCJet_pfakt5", nChg_QC_pfakt5, &b_nChg_QC_pfakt5);
   fChain->SetBranchAddress("nNeutral_ptCutJet_pfakt5", nNeutral_ptCut_pfakt5, &b_nNeutral_ptCut_pfakt5);
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
   fChain->SetBranchAddress("nJetGen_akt7", &nJetGen_akt7, &b_nJetGen_akt7);
   fChain->SetBranchAddress("ptJetGen_akt7 ", ptJetGen_akt7 , &b_ptJetGen_akt7 );
   fChain->SetBranchAddress("eJetGen_akt7  ", eJetGen_akt7  , &b_eJetGen_akt7  );
   fChain->SetBranchAddress("etaJetGen_akt7", etaJetGen_akt7, &b_etaJetGen_akt7);
   fChain->SetBranchAddress("phiJetGen_akt7", phiJetGen_akt7, &b_phiJetGen_akt7);
   fChain->SetBranchAddress("nvertex", &nvertex, &b_nvertex);
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
   fChain->SetBranchAddress("stcMet", &stcMet  , &b_stcMet);
   fChain->SetBranchAddress("etcMet", &etcMet  , &b_etcMet);
   fChain->SetBranchAddress("phitcMet", &phitcMet, &b_phitcMet);
   fChain->SetBranchAddress("signiftcMet", &signiftcMet, &b_signiftcMet);
   fChain->SetBranchAddress("sglobalPfMet", &sglobalPfMet, &b_sglobalPfMet);
   fChain->SetBranchAddress("eglobalPfMet", &eglobalPfMet, &b_eglobalPfMet);
   fChain->SetBranchAddress("phiglobalPfMet", &phiglobalPfMet, &b_phiglobalPfMet);
   fChain->SetBranchAddress("signifglobalPfMet", &signifglobalPfMet, &b_signifglobalPfMet);
   fChain->SetBranchAddress("scentralPfMet", &scentralPfMet, &b_scentralPfMet);
   fChain->SetBranchAddress("ecentralPfMet", &ecentralPfMet, &b_ecentralPfMet);
   fChain->SetBranchAddress("phicentralPfMet", &phicentralPfMet, &b_phicentralPfMet);
   fChain->SetBranchAddress("signifcentralPfMet", &signifcentralPfMet, &b_signifcentralPfMet);
   fChain->SetBranchAddress("eassocPfMet", eassocPfMet, &b_eassocPfMet);
   fChain->SetBranchAddress("phiassocPfMet", phiassocPfMet, &b_phiassocPfMet);
   fChain->SetBranchAddress("signifassocPfMet", signifassocPfMet, &b_signifassocPfMet);
   fChain->SetBranchAddress("eassocOtherVtxPfMet", eassocOtherVtxPfMet, &b_eassocOtherVtxPfMet);
   fChain->SetBranchAddress("phiassocOtherVtxPfMet", phiassocOtherVtxPfMet, &b_phiassocOtherVtxPfMet);
   fChain->SetBranchAddress("signifassocOtherVtxPfMet", signifassocOtherVtxPfMet, &b_signifassocOtherVtxPfMet);
   fChain->SetBranchAddress("etrkPfMet", etrkPfMet, &b_etrkPfMet);
   fChain->SetBranchAddress("phitrkPfMet", phitrkPfMet, &b_phitrkPfMet);
   fChain->SetBranchAddress("signiftrkPfMet", signiftrkPfMet, &b_signiftrkPfMet);
   fChain->SetBranchAddress("ecleanPfMet", ecleanPfMet, &b_ecleanPfMet);
   fChain->SetBranchAddress("phicleanPfMet", phicleanPfMet, &b_phicleanPfMet);
   fChain->SetBranchAddress("signifcleanPfMet", signifcleanPfMet, &b_signifcleanPfMet);
   fChain->SetBranchAddress("ecleanedSaclayPfMet", ecleanedSaclayPfMet, &b_ecleanedSaclayPfMet);
   fChain->SetBranchAddress("phicleanedSaclayPfMet", phicleanedSaclayPfMet, &b_phicleanedSaclayPfMet);
   fChain->SetBranchAddress("signifcleanedSaclayPfMet", signifcleanedSaclayPfMet, &b_signifcleanedSaclayPfMet);
   fChain->SetBranchAddress("eminTypeICleanSaclayPfMet", eminTypeICleanSaclayPfMet, &b_eminTypeICleanSaclayPfMet);
   fChain->SetBranchAddress("phiminTypeICleanSaclayPfMet", phiminTypeICleanSaclayPfMet, &b_phiminTypeICleanSaclayPfMet);
   fChain->SetBranchAddress("signifminTypeICleanSaclayPfMet", signifminTypeICleanSaclayPfMet, &b_signifminTypeICleanSaclayPfMet);
   fChain->SetBranchAddress("globalPfSums", globalPfSums, &b_globalPfSums);
   fChain->SetBranchAddress("spfMet", &spfMet, &b_spfMet);
   fChain->SetBranchAddress("epfMet", &epfMet, &b_epfMet);
   fChain->SetBranchAddress("phipfMet", &phipfMet, &b_phipfMet);
   fChain->SetBranchAddress("signifpfMet", &signifpfMet, &b_signifpfMet);
   fChain->SetBranchAddress("spfMetType1", &spfMetType1, &b_spfMetType1);
   fChain->SetBranchAddress("epfMetType1", &epfMetType1, &b_epfMetType1);
   fChain->SetBranchAddress("phipfMetType1", &phipfMetType1, &b_phipfMetType1);
   fChain->SetBranchAddress("signifpfMetType1", &signifpfMetType1, &b_signifpfMetType1);
   fChain->SetBranchAddress("sMetGen", &sMetGen, &b_sMetGen);
   fChain->SetBranchAddress("eMetGen", &eMetGen, &b_eMetGen);
   fChain->SetBranchAddress("phiMetGen", &phiMetGen, &b_phiMetGen);
   fChain->SetBranchAddress("signifMetGen", &signifMetGen, &b_signifMetGen);
   fChain->SetBranchAddress("sMetGen2", &sMetGen2, &b_sMetGen2);
   fChain->SetBranchAddress("eMetGen2", &eMetGen2, &b_eMetGen2);
   fChain->SetBranchAddress("phiMetGen2", &phiMetGen2, &b_phiMetGen2);
   fChain->SetBranchAddress("vxMC", &vxMC, &b_vxMC);
   fChain->SetBranchAddress("vyMC", &vyMC, &b_vyMC);
   fChain->SetBranchAddress("vzMC", &vzMC, &b_vzMC);
   fChain->SetBranchAddress("vx", vx, &b_vx);
   fChain->SetBranchAddress("vy", vy, &b_vy);
   fChain->SetBranchAddress("vz", vz, &b_vz);
   fChain->SetBranchAddress("vntracks", vntracks, &b_vntracks);
   fChain->SetBranchAddress("vchi2", vchi2, &b_vchi2);
   fChain->SetBranchAddress("vndof", vndof, &b_vndof);
   fChain->SetBranchAddress("vlogsumpt2", vlogsumpt2, &b_vlogsumpt2);
   fChain->SetBranchAddress("nPreselPhotonPairs", &nPreselPhotonPairs, &b_nPreselPhotonPairs);
   fChain->SetBranchAddress("indexPreselPhot1", indexPreselPhot1, &b_indexPreselPhot1);
   fChain->SetBranchAddress("indexPreselPhot2", indexPreselPhot2, &b_indexPreselPhot2);
   fChain->SetBranchAddress("vrankPhotonPairs", vrankPhotonPairs, &b_vrankPhotonPairs);
   fChain->SetBranchAddress("vevtMvaPhotonPairs", vevtMvaPhotonPairs, &b_vevtMvaPhotonPairs);
   fChain->SetBranchAddress("vevtProbPhotonPairs", vevtProbPhotonPairs, &b_vevtProbPhotonPairs);
   fChain->SetBranchAddress("vptbalPhotonPairs", vptbalPhotonPairs, &b_vptbalPhotonPairs);
   fChain->SetBranchAddress("vptasymPhotonPairs", vptasymPhotonPairs, &b_vptasymPhotonPairs);
   fChain->SetBranchAddress("nHLT", &nHLT, &b_nHLT);
   fChain->SetBranchAddress("hltNamesLen", &hltNamesLen, &b_hltNamesLen);
   fChain->SetBranchAddress("HLTNames", &HLTNames, &b_HLTNames);
   fChain->SetBranchAddress("HLTResults", &HLTResults, &b_HLTResults);
   fChain->SetBranchAddress("nEle", &nEle, &b_nEle);
   fChain->SetBranchAddress("electron_px", electron_px, &b_electron_px);
   fChain->SetBranchAddress("electron_py", electron_py, &b_electron_py);
   fChain->SetBranchAddress("electron_pz", electron_pz, &b_electron_pz);
   fChain->SetBranchAddress("electron_vx", electron_vx, &b_electron_vx);
   fChain->SetBranchAddress("electron_vy", electron_vy, &b_electron_vy);
   fChain->SetBranchAddress("electron_vz", electron_vz, &b_electron_vz);
   fChain->SetBranchAddress("electron_pt", electron_pt, &b_electron_pt);
   fChain->SetBranchAddress("electron_eta", electron_eta, &b_electron_eta);
   fChain->SetBranchAddress("electron_phi", electron_phi, &b_electron_phi);
   fChain->SetBranchAddress("electron_energy", electron_energy, &b_electron_energy);
   fChain->SetBranchAddress("electron_charge", electron_charge, &b_electron_charge);
   fChain->SetBranchAddress("electron_fBrem", electron_fBrem, &b_electron_fBrem);
   fChain->SetBranchAddress("electron_dist", electron_dist, &b_electron_dist);
   fChain->SetBranchAddress("electron_dcot", electron_dcot, &b_electron_dcot);
   fChain->SetBranchAddress("electron_matchedConv", electron_matchedConv, &b_electron_matchedConv);
   fChain->SetBranchAddress("electron_misHits", electron_misHits, &b_electron_misHits);
   fChain->SetBranchAddress("electron_seedType", electron_seedType, &b_electron_seedType);
   fChain->SetBranchAddress("electron_EoP", electron_EoP, &b_electron_EoP);
   fChain->SetBranchAddress("electron_OneOverEMinusOneOverP", electron_OneOverEMinusOneOverP, &b_electron_OneOverEMinusOneOverP);
   fChain->SetBranchAddress("electron_r9", electron_r9, &b_electron_r9);
   fChain->SetBranchAddress("electron_nSubClusters", electron_nSubClusters, &b_electron_nSubClusters);
   fChain->SetBranchAddress("electron_trkIso", electron_trkIso, &b_electron_trkIso);
   fChain->SetBranchAddress("electron_ecalIso", electron_ecalIso, &b_electron_ecalIso);
   fChain->SetBranchAddress("electron_hcalIso", electron_hcalIso, &b_electron_hcalIso);
   fChain->SetBranchAddress("electron_trkIso03", electron_trkIso03, &b_electron_trkIso03);
   fChain->SetBranchAddress("electron_ecalIso03", electron_ecalIso03, &b_electron_ecalIso03);
   fChain->SetBranchAddress("electron_hcalIso03", electron_hcalIso03, &b_electron_hcalIso03);
   fChain->SetBranchAddress("electron_SigmaIetaIeta", electron_SigmaIetaIeta, &b_electron_SigmaIetaIeta);
   fChain->SetBranchAddress("electron_SigmaIphiIphi", electron_SigmaIphiIphi, &b_electron_SigmaIphiIphi);
   fChain->SetBranchAddress("electron_dEtaIn", electron_dEtaIn, &b_electron_dEtaIn);
   fChain->SetBranchAddress("electron_dPhiIn", electron_dPhiIn, &b_electron_dPhiIn);
   fChain->SetBranchAddress("electron_HoE", electron_HoE, &b_electron_HoE);
   fChain->SetBranchAddress("electron_pFlowMVA", electron_pFlowMVA, &b_electron_pFlowMVA);
   fChain->SetBranchAddress("electron_sc_energy", electron_sc_energy, &b_electron_sc_energy);
   fChain->SetBranchAddress("electron_sc_eta", electron_sc_eta, &b_electron_sc_eta);
   fChain->SetBranchAddress("electron_sc_phi", electron_sc_phi, &b_electron_sc_phi);
   fChain->SetBranchAddress("electron_ecalEnergy", electron_ecalEnergy, &b_electron_ecalEnergy);
   fChain->SetBranchAddress("electron_trackPatVtx", electron_trackPatVtx, &b_electron_trackPatVtx);
   fChain->SetBranchAddress("electron_mvaNonTrig", electron_mvaNonTrig, &b_electron_mvaNonTrig);
   fChain->SetBranchAddress("electron_mvaTrig", electron_mvaTrig, &b_electron_mvaTrig);
   fChain->SetBranchAddress("electron_chHad03Iso", electron_chHad03Iso, &b_electron_chHad03Iso);
   fChain->SetBranchAddress("electron_nHad03Iso",  electron_nHad03Iso,  &b_electron_nHad03Iso);
   fChain->SetBranchAddress("electron_phot03Iso",  electron_phot03Iso,  &b_electron_phot03Iso);
   fChain->SetBranchAddress("electron_chHad04Iso", electron_chHad04Iso, &b_electron_chHad04Iso);
   fChain->SetBranchAddress("electron_nHad04Iso",  electron_nHad04Iso,  &b_electron_nHad04Iso);
   fChain->SetBranchAddress("electron_phot04Iso",  electron_phot04Iso,  &b_electron_phot04Iso);
   fChain->SetBranchAddress("electron_chHad05Iso", electron_chHad05Iso, &b_electron_chHad05Iso);
   fChain->SetBranchAddress("electron_nHad05Iso",  electron_nHad05Iso,  &b_electron_nHad05Iso);
   fChain->SetBranchAddress("electron_phot05Iso",  electron_phot05Iso,  &b_electron_phot05Iso);
   fChain->SetBranchAddress("isBeamHaloGlobalLoosePass", &isBeamHaloGlobalLoosePass, &b_isBeamHaloGlobalLoosePass);
   fChain->SetBranchAddress("isBeamHaloGlobalTightPass", &isBeamHaloGlobalTightPass, &b_isBeamHaloGloablTightPass);
   fChain->SetBranchAddress("isBeamHaloHcalLoosePass", &isBeamHaloHcalLoosePass, &b_isBeamHaloHcalLoosePass);
   fChain->SetBranchAddress("isBeamHaloHcalTightPass", &isBeamHaloHcalTightPass, &b_isBeamHaloHcalTightPass);
   fChain->SetBranchAddress("isBeamHaloCSCLoosePass", &isBeamHaloCSCLoosePass, &b_isBeamHaloCSCLoosePass);
   fChain->SetBranchAddress("isBeamHaloCSCTightPass", &isBeamHaloCSCTightPass, &b_isBeamHaloCSCTightPass);
   fChain->SetBranchAddress("isBeamHaloEcalLoosePass", &isBeamHaloEcalLoosePass, &b_isBeamHaloEcalLoosePass);
   fChain->SetBranchAddress("isBeamHaloEcalTightPass", &isBeamHaloEcalTightPass, &b_isBeamHaloEcalTightPass);
   fChain->SetBranchAddress("isBeamHaloIDTightPass", &isBeamHaloIDTightPass, &b_isBeamHaloIDTightPass);
   fChain->SetBranchAddress("isBeamHaloIDLoosePass", &isBeamHaloIDLoosePass, &b_isBeamHaloIDLoosePass);
   fChain->SetBranchAddress("isSmellsLikeHalo_Tag", &isSmellsLikeHalo_Tag, &b_isSmellsLikeHalo_Tag);
   fChain->SetBranchAddress("isLooseHalo_Tag", &isLooseHalo_Tag, &b_isLooseHalo_Tag);
   fChain->SetBranchAddress("isTightHalo_Tag", &isTightHalo_Tag, &b_isTightHalo_Tag);
   fChain->SetBranchAddress("isExtremeTightHalo_Tag", &isExtremeTightHalo_Tag, &b_isExtremeTightHalo_Tag);
   fChain->SetBranchAddress("nMuons", &nMuons, &b_nMuons);
   fChain->SetBranchAddress("Muon_px", Muon_px, &b_Muon_px);
   fChain->SetBranchAddress("Muon_py", Muon_py, &b_Muon_py);
   fChain->SetBranchAddress("Muon_pz", Muon_pz, &b_Muon_pz);
   fChain->SetBranchAddress("Muon_vx", Muon_vx, &b_Muon_vx);
   fChain->SetBranchAddress("Muon_vy", Muon_vy, &b_Muon_vy);
   fChain->SetBranchAddress("Muon_vz", Muon_vz, &b_Muon_vz);
   fChain->SetBranchAddress("Muon_pt", Muon_pt, &b_Muon_pt);
   fChain->SetBranchAddress("Muon_eta", Muon_eta, &b_Muon_eta);
   fChain->SetBranchAddress("Muon_phi", Muon_phi, &b_Muon_phi);
   fChain->SetBranchAddress("Muon_energy", Muon_energy, &b_Muon_energy);
   fChain->SetBranchAddress("Muon_charge", Muon_charge, &b_Muon_charge);
   fChain->SetBranchAddress("Muon_isGlobalMuon", Muon_isGlobalMuon, &b_Muon_isGlobalMuon);
   fChain->SetBranchAddress("Muon_isTrackerMuon", Muon_isTrackerMuon, &b_Muon_isTrackerMuon);
   fChain->SetBranchAddress("Muon_isStandAloneMuon", Muon_isStandAloneMuon, &b_Muon_isStandAloneMuon);
   fChain->SetBranchAddress("Muon_InnerTrack_isNonnull", Muon_InnerTrack_isNonnull, &b_Muon_InnerTrack_isNonnull);
   fChain->SetBranchAddress("Muon_OuterTrack_isNonnull", Muon_OuterTrack_isNonnull, &b_Muon_OuterTrack_isNonnull);
   fChain->SetBranchAddress("Muon_OuterPoint_x", Muon_OuterPoint_x, &b_Muon_OuterPoint_x);
   fChain->SetBranchAddress("Muon_OuterPoint_y", Muon_OuterPoint_y, &b_Muon_OuterPoint_y);
   fChain->SetBranchAddress("Muon_OuterPoint_z", Muon_OuterPoint_z, &b_Muon_OuterPoint_z);
   fChain->SetBranchAddress("Muon_InnerPoint_x", Muon_InnerPoint_x, &b_Muon_InnerPoint_x);
   fChain->SetBranchAddress("Muon_InnerPoint_y", Muon_InnerPoint_y, &b_Muon_InnerPoint_y);
   fChain->SetBranchAddress("Muon_InnerPoint_z", Muon_InnerPoint_z, &b_Muon_InnerPoint_z);
   fChain->SetBranchAddress("Muon_trackIso", Muon_trackIso, &b_Muon_trackIso);
   fChain->SetBranchAddress("Muon_ecalIso", Muon_ecalIso, &b_Muon_ecalIso);
   fChain->SetBranchAddress("Muon_hcalIso", Muon_hcalIso, &b_Muon_hcalIso);
   fChain->SetBranchAddress("Muon_relIso", Muon_relIso, &b_Muon_relIso);
   fChain->SetBranchAddress("Muon_normChi2", Muon_normChi2, &b_Muon_normChi2);
   fChain->SetBranchAddress("Muon_validHits", Muon_validHits, &b_Muon_validHits);
   fChain->SetBranchAddress("Muon_tkHits", Muon_tkHits, &b_Muon_tkHits);
   fChain->SetBranchAddress("Muon_pixHits", Muon_pixHits, &b_Muon_pixHits);
   fChain->SetBranchAddress("Muon_numberOfMatches", Muon_numberOfMatches, &b_Muon_numberOfMatches);
   fChain->SetBranchAddress("Muon_pfiso04_chHad",Muon_pfiso04_chHad,&b_Muon_pfiso04_chHad);
   fChain->SetBranchAddress("Muon_pfiso04_chPar",Muon_pfiso04_chPar,&b_Muon_pfiso04_chPar);
   fChain->SetBranchAddress("Muon_pfiso04_nHad", Muon_pfiso04_nHad, &b_Muon_pfiso04_nHad);
   fChain->SetBranchAddress("Muon_pfiso04_Phot", Muon_pfiso04_Phot, &b_Muon_pfiso04_Phot);
   fChain->SetBranchAddress("Muon_pfiso04_PUPt", Muon_pfiso04_PUPt, &b_Muon_pfiso04_PUPt);
   fChain->SetBranchAddress("Muon_pfiso03_chHad",Muon_pfiso03_chHad,&b_Muon_pfiso03_chHad);
   fChain->SetBranchAddress("Muon_pfiso03_chPar",Muon_pfiso03_chPar,&b_Muon_pfiso03_chPar);
   fChain->SetBranchAddress("Muon_pfiso03_nHad", Muon_pfiso03_nHad, &b_Muon_pfiso03_nHad);
   fChain->SetBranchAddress("Muon_pfiso03_Phot", Muon_pfiso03_Phot, &b_Muon_pfiso03_Phot);
   fChain->SetBranchAddress("Muon_pfiso03_PUPt", Muon_pfiso03_PUPt, &b_Muon_pfiso03_PUPt);
   fChain->SetBranchAddress("Muon_isPFMuon", Muon_isPFMuon, &b_Muon_isPFMuon); 
   fChain->SetBranchAddress("Muon_trkLayerWithMeas", Muon_trkLayerWithMeas, &b_Muon_trkLayerWithMeas);
   fChain->SetBranchAddress("nCosmicMuons", &nCosmicMuons, &b_nCosmicMuons);
   fChain->SetBranchAddress("CosmicMuon_px", CosmicMuon_px, &b_CosmicMuon_px);
   fChain->SetBranchAddress("CosmicMuon_py", CosmicMuon_py, &b_CosmicMuon_py);
   fChain->SetBranchAddress("CosmicMuon_pz", CosmicMuon_pz, &b_CosmicMuon_pz);
   fChain->SetBranchAddress("CosmicMuon_pt", CosmicMuon_pt, &b_CosmicMuon_pt);
   fChain->SetBranchAddress("CosmicMuon_eta", CosmicMuon_eta, &b_CosmicMuon_eta);
   fChain->SetBranchAddress("CosmicMuon_phi", CosmicMuon_phi, &b_CosmicMuon_phi);
   fChain->SetBranchAddress("CosmicMuon_energy", CosmicMuon_energy, &b_CosmicMuon_energy);
   fChain->SetBranchAddress("CosmicMuon_charge", CosmicMuon_charge, &b_CosmicMuon_charge);
   fChain->SetBranchAddress("CosmicMuon_isGlobalMuon", CosmicMuon_isGlobalMuon, &b_CosmicMuon_isGlobalMuon);
   fChain->SetBranchAddress("CosmicMuon_isTrackerMuon", CosmicMuon_isTrackerMuon, &b_CosmicMuon_isTrackerMuon);
   fChain->SetBranchAddress("CosmicMuon_isStandAloneMuon", CosmicMuon_isStandAloneMuon, &b_CosmicMuon_isStandAloneMuon);
   fChain->SetBranchAddress("CosmicMuon_InnerTrack_isNonnull", CosmicMuon_InnerTrack_isNonnull, &b_CosmicMuon_InnerTrack_isNonnull);
   fChain->SetBranchAddress("CosmicMuon_OuterTrack_isNonnull", CosmicMuon_OuterTrack_isNonnull, &b_CosmicMuon_OuterTrack_isNonnull);
   fChain->SetBranchAddress("CosmicMuon_OuterPoint_x", CosmicMuon_OuterPoint_x, &b_CosmicMuon_OuterPoint_x);
   fChain->SetBranchAddress("CosmicMuon_OuterPoint_y", CosmicMuon_OuterPoint_y, &b_CosmicMuon_OuterPoint_y);
   fChain->SetBranchAddress("CosmicMuon_OuterPoint_z", CosmicMuon_OuterPoint_z, &b_CosmicMuon_OuterPoint_z);
   fChain->SetBranchAddress("Xsec", &Xsec, &b_Xsec);
   Notify();
}

Bool_t tree_reader_V8::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void tree_reader_V8::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t tree_reader_V8::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef tree_reader_V8_cxx
