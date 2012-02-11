#include "RedNtpTree.h"
#include "JSON.h"
#include <TH2.h>
#include <TStyle.h>
#include <TString.h>
#include <TObjArray.h>
#include <TCanvas.h>
#include <TMath.h>
#include <TStopwatch.h>
#include <iostream>
#include <vector>
#include <TLorentzVector.h>

#ifdef SMALL_VERTEX_VECTOR
#define MAX_PU_REWEIGHT 30
#else
#define MAX_PU_REWEIGHT 40
#endif

//#define DEBUG
//#define DEBUG1

using std::cout;
using std::endl;


RedNtpTree::RedNtpTree(TTree *tree, const TString& outname) : tree_reader_V7(tree), jsonFile(0) , ptweights_(0), scaleCorrections_(0)
{  
  hOutputFile   = TFile::Open(outname, "RECREATE" ) ;
  // must be set by the user 
  EquivLumi = -1.;
  xsection = -1.;
  NtotEvents = -1;
  SampleID = -1;
  gen_=new TRandom3(0);
  doPDFweight = 0;
  jetsyst_ = 0;
}



inline double delta_phi(double phi1, double phi2) {

  double dphi = TMath::Abs(phi1 - phi2);
  return (dphi <= TMath::Pi())? dphi : TMath::TwoPi() - dphi;

}

inline double delta_eta(double eta1, double eta2) {

  return (eta2 >= 0 ? eta1 - eta2 : eta2 - eta1);
}

RedNtpTree::~RedNtpTree()
{
   hOutputFile->Write() ;
   hOutputFile->Close() ;
   hOutputFile->Delete();
}



vector<int>  RedNtpTree::firstfour(Float_t *vec, vector<bool> *asso){

    // double max(-999); int idmax(-999);
    // double secondmax(-999); int idsecondmax(-999);
    // 
    // for (int i=0; i<int(asso->size()); i++) {
    // 
    //   if ( vec[i] > max && asso->at(i)) {
    //     max = vec[i];
    //     idmax = i;
    //   }
    // 
    // }
    // for (int i=0; i<int(asso->size()); i++) {
    // 
    //   if ( vec[i] > secondmax && asso->at(i) && i!= idmax) {
    //     secondmax = vec[i];
    //     idsecondmax = i;
    //   }
    // 
    // }
  
    vector<int> themax;
  
    for(int j=0; j<4; j++)
    {
        double maxtemp(-999); 
        int idmaxtemp(-999);
 
        for (int i=0; i<int(asso->size()); i++) 
        {
            bool skip(0);
            for(int ss=0; ss<j; ss++) 
            {
	            if ( i == themax.at(ss) )   
	                skip = 1;
            }
            if ( vec[i] > maxtemp && asso->at(i) && !skip) 
            {
	            maxtemp = vec[i];
	            idmaxtemp = i;
            }
        }
        themax.push_back(idmaxtemp);
    }
    return themax;
}




bool RedNtpTree::mcID(int i) 
{
    bool assoc(0);
    for(int j=0; j<nMC; j++)
    {
        double DR, DE;
    
        if(pdgIdMC[j] == 22 && statusMC[j] == 3)
        {
            DR = sqrt(delta_eta(etaPhot[i],etaMC[j])*delta_eta(etaPhot[i],etaMC[j]) + 
		        delta_phi(phiPhot[i],phiMC[j])*delta_phi(phiPhot[i],phiMC[j]) ) ;
            DE = TMath::Abs(ePhot[i]-eMC[j])/ePhot[i];
            if(DR < .1 && DE < .2) assoc = 1; 
        }
    }
    return assoc;
}


void RedNtpTree::Loop(int isgjetqcd, char* selection)
{
    if (fChain == 0) return;
    
    Long64_t nentries = fChain->GetEntriesFast();
    //   Long64_t nentries = 10000;
    
    Long64_t nbytes = 0, nb = 0;
    
    TStopwatch timer;
    
    //   JSON myjson("Cert_160404-163869_7TeV_May10ReReco_Collisions11_CMSSWConfig.txt");
    JSON* myjson=0;
    if (jsonFile)
    {
        std::cout << "Reading JSON" << jsonFile << std::endl;
        myjson=new JSON(jsonFile);
    }
    
    // hOutputFile = new TFile("output.root" , "RECREATE" ) ;
    
    hOutputFile->cd();   

    /********************************************************
     *                                                      *
     *                      HISTO INIT                      *
     *                                                      *
     ********************************************************/

    TH1D Dvz("Dvz","Dvz", 200, -10.,10.);
    TH1D Dvzbest("Dvzbest","Dvzbest", 200, -10.,10.);
    TH2D JECunc("JECunc","JECunc", 100, 0.,200.,100,0.,0.2);
    TH1D JECresovbf("JECresovbf","JECresovbf", 100, -0.5,0.5);
    TH1D JECresovh("JECresovh","JECresovh", 100, -0.5,0.5);
    
    TH1D nPDFweight1("nPDFweight1","nPDFweight1", 150, 0.,150.);
    TH1D nPDFweight2("nPDFweight2","nPDFweight2", 150, 0.,150.);
    TH1D nPDFweight3("nPDFweight3","nPDFweight3", 150, 0.,150.);
    TH1D nPDFweight4("nPDFweight4","nPDFweight4", 150, 0.,150.);
    TH1D nPDFweight5("nPDFweight5","nPDFweight5", 150, 0.,150.);
    TH1D nPDFweight6("nPDFweight6","nPDFweight6", 150, 0.,150.);
    TH1D nPDFweight7("nPDFweight7","nPDFweight7", 150, 0.,150.);
    TH1D nPDFweight8("nPDFweight8","nPDFweight8", 150, 0.,150.);
    TH1D nPDFweight9("nPDFweight9","nPDFweight9", 150, 0.,150.);
    TH1D nPDFweight10("nPDFweight10","nPDFweight10",150, 0.,150.);
    
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
    TH1D invmassjetgen("invmassjetgen","invmassjetgen", 100, 0.,170.);
    TH1D npunorew("npunorew","npunorew", 25, -0.5,24.5);
    TH1D npurew("npurew","npurew", 25, -0.5,24.5);
    TH1D nvtxnorew("nvtxnorew","nvtxnorew", 25, -0.5,24.5);
    TH1D nvtxrew("nvtxrew","nvtxrew", 25, -0.5,24.5);
    
    // isolation variables
    TH1D hcalisoassphot_EB("hcalisoassphot_EB","hcalisoassphot_EB",100,0.,1.);
    TH1D ecalisoassphot_EB("ecalisoassphot_EB","ecalisoassphot_EB",100,0.,1.);
    TH1D ptisoassphot_EB("ptisoassphot_EB","ptisoassphot_EB",100,0.,1.);
    TH1D ntrkisoassphot_EB("ntrkisoassphot_EB","ntrkisoassphot_EB",10,0.,10);
    TH1D sminminclusassphot_EB("sminminclusassphot_EB","sminminclusassphot_EB",100,0.,1.);
    TH1D smaxmaxclusassphot_EB("smaxmaxclusassphot_EB","smaxmaxclusassphot_EB",100,0.,1.);
    TH1D alphaclusassphot_EB("alphaclusassphot_EB","alphaclusassphot_EB",50,-1.57,1.57);
    TH1D hcalisoassjet_EB("hcalisoassjet_EB","hcalisoassjet_EB",100,0.,1.);
    TH1D ecalisoassjet_EB("ecalisoassjet_EB","ecalisoassjet_EB",100,0.,1.);
    TH1D ptisoassjet_EB("ptisoassjet_EB","ptisoassjet_EB",100,0.,1.);
    TH1D ntrkisoassjet_EB("ntrkisoassjet_EB","ntrkisoassjet_EB",10,0.,10);
    TH1D sminminclusassjet_EB("sminminclusassjet_EB","sminminclusassjet_EB",100,0.,1.);
    TH1D smaxmaxclusassjet_EB("smaxmaxclusassjet_EB","smaxmaxclusassjet_EB",100,0.,1.);
    TH1D alphaclusassjet_EB("alphaclusassjet_EB","alphaclusassjet_EB",50,-1.57,1.57);
    TH1D hcalisoassphot_EE("hcalisoassphot_EE","hcalisoassphot_EE",100,0.,1.);
    TH1D ecalisoassphot_EE("ecalisoassphot_EE","ecalisoassphot_EE",100,0.,1.);
    TH1D ptisoassphot_EE("ptisoassphot_EE","ptisoassphot_EE",100,0.,1.);
    TH1D ntrkisoassphot_EE("ntrkisoassphot_EE","ntrkisoassphot_EE",10,0.,10);
    TH1D sminminclusassphot_EE("sminminclusassphot_EE","sminminclusassphot_EE",100,0.,1.);
    TH1D smaxmaxclusassphot_EE("smaxmaxclusassphot_EE","smaxmaxclusassphot_EE",100,0.,1.);
    TH1D alphaclusassphot_EE("alphaclusassphot_EE","alphaclusassphot_EE",50,-1.57,1.57);
    TH1D hcalisoassjet_EE("hcalisoassjet_EE","hcalisoassjet_EE",100,0.,1.);
    TH1D ecalisoassjet_EE("ecalisoassjet_EE","ecalisoassjet_EE",100,0.,1.);
    TH1D ptisoassjet_EE("ptisoassjet_EE","ptisoassjet_EE",100,0.,1.);
    TH1D ntrkisoassjet_EE("ntrkisoassjet_EE","ntrkisoassjet_EE",10,0.,10);
    TH1D sminminclusassjet_EE("sminminclusassjet_EE","sminminclusassjet_EE",100,0.,1.);
    TH1D smaxmaxclusassjet_EE("smaxmaxclusassjet_EE","smaxmaxclusassjet_EE",100,0.,1.);
    TH1D alphaclusassjet_EE("alphaclusassjet_EE","alphaclusassjet_EE",50,-1.57,1.57);
      

    /********************************************************
     *                                                      *
     *                      TREE INIT                       *
     *                                                      *
     ********************************************************/


    ana_tree = new TTree ("AnaTree","Reduced tree for final analysis") ;
    ana_tree->Branch("run",&runRN,"run/I");
    ana_tree->Branch("event",&eventRN,"event/I");
    ana_tree->Branch("lumi",&lumi,"lumi/I");
    ana_tree->Branch("rhoPF",&rhoPFRN,"rhoPF/F");
    ana_tree->Branch("massgg",&massgg,"massgg/F");
    ana_tree->Branch("ptgg",&ptgg,"ptgg/F");
    ana_tree->Branch("ptggnewvtx",&ptggnewvtx,"ptggnewvtx/F");
    ana_tree->Branch("massggnewvtx",&massggnewvtx,"massggnewvtx/F");
    ana_tree->Branch("ptphot1",&ptphot1,"ptphot1/F");
    ana_tree->Branch("ptphot2",&ptphot2,"ptphot2/F");
    ana_tree->Branch("timephot1",&timephot1,"timephot1/F"); 
    ana_tree->Branch("timephot2",&timephot2,"timephot2/F"); 
    ana_tree->Branch("etaphot1",&etaphot1,"etaphot1/F");
    ana_tree->Branch("etaphot2",&etaphot2,"etaphot2/F");
    ana_tree->Branch("phiphot1",&phiphot1,"phiphot1/F");
    ana_tree->Branch("phiphot2",&phiphot2,"phiphot2/F");
    ana_tree->Branch("etascphot1",&etascphot1,"etascphot1/F");
    ana_tree->Branch("etascphot2",&etascphot2,"etascphot2/F");
    ana_tree->Branch("phiscphot1",&phiscphot1,"phiscphot1/F");
    ana_tree->Branch("phiscphot2",&phiscphot2,"phiscphot2/F");
    ana_tree->Branch("E1phot1",&E1phot1,"E1phot1/F");
    ana_tree->Branch("E1phot2",&E1phot2,"E1phot2/F");
    ana_tree->Branch("E9phot1",&E9phot1,"E9phot1/F");
    ana_tree->Branch("E9phot2",&E9phot2,"E9phot2/F");
    ana_tree->Branch("r9phot1",&r9phot1,"r9phot1/F");
    ana_tree->Branch("r9phot2",&r9phot2,"r9phot2/F");
    ana_tree->Branch("isemEGphot1",&isemEGphot1,"isemEGphot1/I");
    ana_tree->Branch("isemEGphot2",&isemEGphot2,"isemEGphot2/I");
    ana_tree->Branch("idloosenewEGphot1",&idloosenewEGphot1,"idloosenewEGphot1/I");
    ana_tree->Branch("idloosenewEGphot2",&idloosenewEGphot2,"idloosenewEGphot2/I");
    ana_tree->Branch("idloose006newEGphot1",&idloose006newEGphot1,"idloose006newEGphot1/I");
    ana_tree->Branch("idloose006newEGphot2",&idloose006newEGphot2,"idloose006newEGphot2/I");
    ana_tree->Branch("idtightnewEGphot1",&idtightnewEGphot1,"idtightnewEGphot1/I");
    ana_tree->Branch("idtightnewEGphot2",&idtightnewEGphot2,"idtightnewEGphot2/I");
    ana_tree->Branch("idhggtightnewEGphot1",&idhggtightnewEGphot1,"idhggtightnewEGphot1/I");
    ana_tree->Branch("idhggtightnewEGphot2",&idhggtightnewEGphot2,"idhggtightnewEGphot2/I");
    ana_tree->Branch("idloosenewpuEGphot1",&idloosenewpuEGphot1,"idloosenewpuEGphot1/I");
    ana_tree->Branch("idloosenewpuEGphot2",&idloosenewpuEGphot2,"idloosenewpuEGphot2/I");
    ana_tree->Branch("idtightnewpuEGphot1",&idtightnewpuEGphot1,"idtightnewpuEGphot1/I");
    ana_tree->Branch("idtightnewpuEGphot2",&idtightnewpuEGphot2,"idtightnewpuEGphot2/I");
    ana_tree->Branch("idhggtightnewpuEGphot1",&idhggtightnewpuEGphot1,"idhggtightnewpuEGphot1/I");
    ana_tree->Branch("idhggtightnewpuEGphot2",&idhggtightnewpuEGphot2,"idhggtightnewpuEGphot2/I");
    ana_tree->Branch("idcicphot1",&idcicphot1,"idcicphot1/I");
    ana_tree->Branch("idcicphot2",&idcicphot2,"idcicphot2/I");
    ana_tree->Branch("idcicnoelvetophot1",&idcicnoelvetophot1,"idcicnoelvetophot1/I");
    ana_tree->Branch("idcicnoelvetophot2",&idcicnoelvetophot2,"idcicnoelvetophot2/I");
    ana_tree->Branch("idlooseEGphot1",&idlooseEGphot1,"idlooseEGphot1/I");
    ana_tree->Branch("idlooseEGphot2",&idlooseEGphot2,"idlooseEGphot2/I");
    ana_tree->Branch("idtightEGphot1",&idtightEGphot1,"idtightEGphot1/I");
    ana_tree->Branch("idtightEGphot2",&idtightEGphot2,"idtightEGphot2/I");
    ana_tree->Branch("idloosephot1",&idloosephot1,"idloosephot1/I");
    ana_tree->Branch("idloosephot2",&idloosephot2,"idloosephot2/I");
    ana_tree->Branch("idmediumphot1",&idmediumphot1,"idmediumphot1/I");
    ana_tree->Branch("idmediumphot2",&idmediumphot2,"idmediumphot2/I");
    ana_tree->Branch("idloosecsphot1",&idloosecsphot1,"idloosecsphot1/I"); 
    ana_tree->Branch("idloosecsphot2",&idloosecsphot2,"idloosecsphot2/I"); 
    ana_tree->Branch("idmediumcsphot1",&idmediumcsphot1,"idmediumcsphot1/I"); 
    ana_tree->Branch("idmediumcsphot2",&idmediumcsphot2,"idmediumcsphot2/I"); 
    ana_tree->Branch("idelephot1",&idelephot1,"idelephot1/I");
    ana_tree->Branch("idelephot2",&idelephot2,"idelephot2/I");
    
    ana_tree->Branch("pid_isEMphot1",&pid_isEMphot1,"pid_isEMphot1/I");
    ana_tree->Branch("pid_isEMphot2",&pid_isEMphot2,"pid_isEMphot2/I");
    
    ana_tree->Branch("pid_haspixelseedphot1",&pid_haspixelseedphot1,"pid_haspixelseedphot1/I");
    ana_tree->Branch("pid_haspixelseedphot2",&pid_haspixelseedphot2,"pid_haspixelseedphot2/I");
    ana_tree->Branch("pid_jurECALphot1",&pid_jurECALphot1,"pid_jurECALphot1/F"); 
    ana_tree->Branch("pid_jurECALphot2",&pid_jurECALphot2,"pid_jurECALphot2/F"); 
    ana_tree->Branch("pid_twrHCALphot1",&pid_twrHCALphot1,"pid_twrHCALphot1/F");
    ana_tree->Branch("pid_twrHCALphot2",&pid_twrHCALphot2,"pid_twrHCALphot2/F");
    ana_tree->Branch("pid_HoverEphot1",&pid_HoverEphot1,"pid_HoverEphot1/F");
    ana_tree->Branch("pid_HoverEphot2",&pid_HoverEphot2,"pid_HoverEphot2/F");
    ana_tree->Branch("pid_hlwTrackphot1",&pid_hlwTrackphot1,"pid_hlwTrackphot1/F");
    ana_tree->Branch("pid_hlwTrackphot2",&pid_hlwTrackphot2,"pid_hlwTrackphot2/F");
    ana_tree->Branch("pid_etawidphot1",&pid_etawidphot1,"pid_etawidphot1/F");
    ana_tree->Branch("pid_etawidphot2",&pid_etawidphot2,"pid_etawidphot2/F");
    ana_tree->Branch("pid_hlwTrackNoDzphot1",&pid_hlwTrackNoDzphot1,"pid_hlwTrackNoDzphot1/F");
    ana_tree->Branch("pid_hlwTrackNoDzphot2",&pid_hlwTrackNoDzphot2,"pid_hlwTrackNoDzphot2/F");
    ana_tree->Branch("pid_hasMatchedConvphot1",&pid_hasMatchedConvphot1,"pid_hasMatchedConvphot1/I");
    ana_tree->Branch("pid_hasMatchedConvphot2",&pid_hasMatchedConvphot2,"pid_hasMatchedConvphot2/I");
    ana_tree->Branch("pid_hasMatchedPromptElephot1",&pid_hasMatchedPromptElephot1,"pid_hasMatchedPromptElephot1/I");
    ana_tree->Branch("pid_hasMatchedPromptElephot2",&pid_hasMatchedPromptElephot2,"pid_hasMatchedPromptElephot2/I");
    
    ana_tree->Branch("pid_sminphot1",&pid_sminphot1,"pid_sminphot1/F");
    ana_tree->Branch("pid_sminphot2",&pid_sminphot2,"pid_sminphot2/F");
    ana_tree->Branch("pid_smajphot1",&pid_smajphot1,"pid_smajphot1/F");
    ana_tree->Branch("pid_smajphot2",&pid_smajphot2,"pid_smajphot2/F");
    ana_tree->Branch("pid_ntrkphot1",&pid_ntrkphot1,"pid_ntrkphot1/I");
    ana_tree->Branch("pid_ntrkphot2",&pid_ntrkphot2,"pid_ntrkphot2/I");
    ana_tree->Branch("pid_ptisophot1",&pid_ptisophot1,"pid_ptisophot1/F");
    ana_tree->Branch("pid_ptisophot2",&pid_ptisophot2,"pid_ptisophot2/F");
    ana_tree->Branch("pid_ntrkcsphot1",&pid_ntrkcsphot1,"pid_ntrkcsphot1/I"); 
    ana_tree->Branch("pid_ntrkcsphot2",&pid_ntrkcsphot2,"pid_ntrkcsphot2/I"); 
    ana_tree->Branch("pid_ptisocsphot1",&pid_ptisocsphot1,"pid_ptisocsphot1/F"); 
    ana_tree->Branch("pid_ptisocsphot2",&pid_ptisocsphot2,"pid_ptisocsphot2/F"); 
    ana_tree->Branch("pid_ecalisophot1",&pid_ecalisophot1,"pid_ecalisophot1/F");
    ana_tree->Branch("pid_ecalisophot2",&pid_ecalisophot2,"pid_ecalisophot2/F");
    ana_tree->Branch("pid_hcalisophot1",&pid_hcalisophot1,"pid_hcalisophot1/F");
    ana_tree->Branch("pid_hcalisophot2",&pid_hcalisophot2,"pid_hcalisophot2/F");
    
    ana_tree->Branch("ptjet1",&ptjet1,"ptjet1/F");
    ana_tree->Branch("ptjet2",&ptjet2,"ptjet2/F");
    ana_tree->Branch("ptjet3",&ptjet3,"ptjet3/F");
    ana_tree->Branch("ptjet4",&ptjet4,"ptjet4/F");
    ana_tree->Branch("ptcorrjet1",&ptcorrjet1,"ptcorrjet1/F");
    ana_tree->Branch("ptcorrjet2",&ptcorrjet2,"ptcorrjet2/F");
    ana_tree->Branch("ptcorrjet3",&ptcorrjet3,"ptcorrjet3/F");
    ana_tree->Branch("ptcorrjet4",&ptcorrjet4,"ptcorrjet4/F");
    ana_tree->Branch("etajet1",&etajet1,"etajet1/F");
    ana_tree->Branch("etajet2",&etajet2,"etajet2/F");
    ana_tree->Branch("etajet3",&etajet3,"etajet3/F");
    ana_tree->Branch("etajet4",&etajet4,"etajet4/F");
    ana_tree->Branch("phijet1",&phijet1,"phijet1/F");
    ana_tree->Branch("phijet2",&phijet2,"phijet2/F");
    ana_tree->Branch("phijet3",&phijet3,"phijet3/F");
    ana_tree->Branch("phijet4",&phijet4,"phijet4/F");
    ana_tree->Branch("betajet1",&betajet1,"betajet1/F");
    ana_tree->Branch("betajet2",&betajet2,"betajet2/F");
    ana_tree->Branch("betastarjet1",&betastarjet1,"betastarjet1/F");
    ana_tree->Branch("betastarjet2",&betastarjet2,"betastarjet2/F");
    ana_tree->Branch("assjet1",&assjet1,"assjet1/I");
    ana_tree->Branch("assjet2",&assjet2,"assjet2/I");
    ana_tree->Branch("deltaeta",&deltaeta,"deltaeta/F");
    ana_tree->Branch("zeppenjet",&zeppenjet,"zeppenjet/F");
    ana_tree->Branch("deltaphi",&deltaphi,"deltaphi/F");
    ana_tree->Branch("deltaphinewvtx",&deltaphinewvtx,"deltaphinewvtx/F");
    ana_tree->Branch("invmassjet",&invmassjet,"invmassjet/F");
    ana_tree->Branch("invmass2g1j",&invmass2g1j,"invmass2g1j/F");
    ana_tree->Branch("invmass2g2j",&invmass2g2j,"invmass2g2j/F");
    ana_tree->Branch("nvtx",&nvtx,"nvtx/F");
    
    // ana_tree->Branch("met",&met,"met/F");
    // ana_tree->Branch("phimet",&phimet,"phimet/F");
    
    ana_tree->Branch("sMet", &sMet_, "sMet/F")  ;
    ana_tree->Branch("eMet", &eMet_, "eMet/F")  ;
    ana_tree->Branch("phiMet", &phiMet_, "phiMet/F");
    ana_tree->Branch("signifMet", &signifMet_, "signifMet/F");
    ana_tree->Branch("sCorrMet", &sCorrMet_, "sCorrMet/F")  ;
    ana_tree->Branch("eCorrMet", &eCorrMet_, "eCorrMet/F")  ;
    ana_tree->Branch("phiCorrMet", &phiCorrMet_, "phiCorrMet/F");
    ana_tree->Branch("signifCorrMet", &signifCorrMet_, "signifCorrMet/F");
    ana_tree->Branch("smuCorrMet", &smuCorrMet_, "smuCorrMet/F")  ;
    ana_tree->Branch("emuCorrMet", &emuCorrMet_, "emuCorrMet/F")  ;
    ana_tree->Branch("phimuCorrMet", &phimuCorrMet_, "phimuCorrMet/F");
    ana_tree->Branch("signifmuCorrMet", &signifmuCorrMet_, "signifmuCorrMet/F");
    ana_tree->Branch("sNoHFMet", &sNoHFMet_, "sNoHFMet/F")  ;
    ana_tree->Branch("eNoHFMet", &eNoHFMet_, "eNoHFMet/F")  ;
    ana_tree->Branch("phiNoHFMet", &phiNoHFMet_, "phiNoHFMet/F");
    ana_tree->Branch("signifNoHFMet", &signifNoHFMet_, "signifNoHFMet/F");
    ana_tree->Branch("stcMet", &stcMet_, "stcMet/F")  ;
    ana_tree->Branch("etcMet", &etcMet_, "etcMet/F")  ;
    ana_tree->Branch("phitcMet", &phitcMet_, "phitcMet/F");
    ana_tree->Branch("signiftcMet", &signiftcMet_, "signiftcMet/F");
    ana_tree->Branch("sglobalPfMet", &sglobalPfMet_, "sglobalPfMet/F");
    ana_tree->Branch("eglobalPfMet", &eglobalPfMet_, "eglobalPfMet/F");
    ana_tree->Branch("phiglobalPfMet", &phiglobalPfMet_, "phiglobalPfMet/F");
    ana_tree->Branch("signifglobalPfMet", &signifglobalPfMet_, "signifglobalPfMet/F");
    ana_tree->Branch("scentralPfMet", &scentralPfMet_, "scentralPfMet/F");
    ana_tree->Branch("ecentralPfMet", &ecentralPfMet_, "ecentralPfMet/F");
    ana_tree->Branch("phicentralPfMet", &phicentralPfMet_, "phicentralPfMet/F");
    ana_tree->Branch("signifcentralPfMet", &signifcentralPfMet_, "signifcentralPfMet/F");
    ana_tree->Branch("eassocPfMet", &eassocPfMet_, "eassocPfMet/F");   //[nvertex]
    ana_tree->Branch("phiassocPfMet", &phiassocPfMet_, "phiassocPfMet/F");   //[nvertex]
    ana_tree->Branch("signifassocPfMet", &signifassocPfMet_, "signifassocPfMet/F");   //[nvertex]
    ana_tree->Branch("eassocOtherVtxPfMet", &eassocOtherVtxPfMet_, "eassocOtherVtxPfMet/F");   //[nvertex]
    ana_tree->Branch("phiassocOtherVtxPfMet", &phiassocOtherVtxPfMet_, "phiassocOtherVtxPfMet/F");   //[nvertex]
    ana_tree->Branch("signifassocOtherVtxPfMet", &signifassocOtherVtxPfMet_, "signifassocOtherVtxPfMet/F");   //[nvertex]
    ana_tree->Branch("etrkPfMet", &etrkPfMet_, "etrkPfMet/F");   //[nvertex]
    ana_tree->Branch("phitrkPfMet", &phitrkPfMet_, "phitrkPfMet/F");   //[nvertex]
    ana_tree->Branch("signiftrkPfMet", &signiftrkPfMet_, "signiftrkPfMet/F");   //[nvertex]
    ana_tree->Branch("ecleanPfMet", &ecleanPfMet_, "ecleanPfMet/F");   //[nvertex]
    ana_tree->Branch("phicleanPfMet", &phicleanPfMet_, "phicleanPfMet/F");   //[nvertex]
    ana_tree->Branch("signifcleanPfMet", &signifcleanPfMet_, "signifcleanPfMet/F");   //[nvertex]
    ana_tree->Branch("ecleanedSaclayPfMet", &ecleanedSaclayPfMet_, "ecleanedSaclayPfMet/F");   //[nvertex]
    ana_tree->Branch("phicleanedSaclayPfMet", &phicleanedSaclayPfMet_, "phicleanedSaclayPfMet/F");   //[nvertex]
    ana_tree->Branch("signifcleanedSaclayPfMet", &signifcleanedSaclayPfMet_, "signifcleanedSaclayPfMet/F");   //[nvertex]
    ana_tree->Branch("eminTypeICleanSaclayPfMet", &eminTypeICleanSaclayPfMet_, "eminTypeICleanSaclayPfMet/F");   //[nvertex]
    ana_tree->Branch("phiminTypeICleanSaclayPfMet", &phiminTypeICleanSaclayPfMet_, "phiminTypeICleanSaclayPfMet/F");   //[nvertex]
    ana_tree->Branch("signifminTypeICleanSaclayPfMet", &signifminTypeICleanSaclayPfMet_, "signifminTypeICleanSaclayPfMet/F");   //[nvertex]
    ana_tree->Branch("globalPfSums", &globalPfSums_, "globalPfSums/F");
    ana_tree->Branch("spfMet", &spfMet_, "spfMet/F")  ;
    ana_tree->Branch("epfMet", &epfMet_, "epfMet/F")  ;
    ana_tree->Branch("phipfMet", &phipfMet_, "phipfMet/F");
    ana_tree->Branch("signifpfMet", &signifpfMet_, "signifpfMet/F");
    ana_tree->Branch("spfMetType1", &spfMetType1_, "spfMetType1/F");
    ana_tree->Branch("epfMetType1", &epfMetType1_, "epfMetType1/F");
    ana_tree->Branch("phipfMetType1", &phipfMetType1_, "phipfMetType1/F");
    ana_tree->Branch("signifpfMetType1", &signifpfMetType1_, "signifpfMetType1/F");
    ana_tree->Branch("sMetGen", &sMetGen_, "sMetGen/F")  ;
    ana_tree->Branch("eMetGen", &eMetGen_, "eMetGen/F")  ;
    ana_tree->Branch("phiMetGen", &phiMetGen_, "phiMetGen/F");
    ana_tree->Branch("signifMetGen", &signifMetGen_, "signifMetGen/F");
    ana_tree->Branch("sMetGen2", &sMetGen2_, "sMetGen2/F")  ;
    ana_tree->Branch("eMetGen2", &eMetGen2_, "eMetGen2/F")  ;
    ana_tree->Branch("phiMetGen2", &phiMetGen2_, "phiMetGen2/F");
    
    ana_tree->Branch("npu",&npu,"npu/I");
    ana_tree->Branch("NtotEvents",&NtotEvents,"NtotEvents/I");
    ana_tree->Branch("xsection",&xsection,"xsection/F");
    ana_tree->Branch("EquivLumi",&EquivLumi,"EquivLumi/F");
    ana_tree->Branch("SampleID",&SampleID,"SampleID/I");
    
    ana_tree->Branch("pu_weight",&pu_weight,"pu_weight/F");
    ana_tree->Branch("pt_weight",&pt_weight,"pt_weight/F");
    
    
    ana_tree->Branch("gen_custom_processId" , &gen_custom_processId, "gen_custom_processId/I");

    ana_tree->Branch("gen_pt_gamma1", &gen_pt_gamma1, "gen_pt_gamma1/F");
    ana_tree->Branch("gen_pt_gamma2", &gen_pt_gamma2, "gen_pt_gamma2/F");
    ana_tree->Branch("gen_eta_gamma1", &gen_eta_gamma1, "gen_eta_gamma1/F");
    ana_tree->Branch("gen_eta_gamma2", &gen_eta_gamma2, "gen_eta_gamma2/F");
    ana_tree->Branch("gen_phi_gamma1", &gen_phi_gamma1, "gen_phi_gamma1/F");
    ana_tree->Branch("gen_phi_gamma2", &gen_phi_gamma2, "gen_phi_gamma2/F");
    
    ana_tree->Branch("gen_pt_genjet1",      &gen_pt_genjet1,      "gen_pt_genjet1/F");         
    ana_tree->Branch("gen_pt_genjet2",      &gen_pt_genjet2,      "gen_pt_genjet2/F");       
    ana_tree->Branch("gen_eta_genjet1",     &gen_eta_genjet1,     "gen_eta_genjet1/F");         
    ana_tree->Branch("gen_eta_genjet2",     &gen_eta_genjet2,     "gen_eta_genjet2/F");        
    ana_tree->Branch("gen_phi_genjet1",     &gen_phi_genjet1,     "gen_phi_genjet1/F");        
    ana_tree->Branch("gen_phi_genjet2",     &gen_phi_genjet2,     "gen_phi_genjet2/F");         
    // ana_tree->Branch("gen_pt_VectorBoson",  &gen_pt_VectorBoson,  "gen_pt_VectorBoson/F");         
    // ana_tree->Branch("gen_phi_VectorBoson", &gen_phi_VectorBoson, "gen_phi_VectorBoson/F");         
    // ana_tree->Branch("gen_eta_VectorBoson", &gen_eta_VectorBoson, "gen_eta_VectorBoson/F");         
    ana_tree->Branch("gen_mass_diphoton",   &gen_mass_diphoton,   "gen_mass_diphoton/F");         
    ana_tree->Branch("gen_pt_diphoton",     &gen_pt_diphoton,     "gen_pt_diphoton/F");         
    ana_tree->Branch("gen_eta_diphoton",    &gen_eta_diphoton,    "gen_eta_diphoton/F");         
    ana_tree->Branch("gen_phi_diphoton",    &gen_phi_diphoton,    "gen_phi_diphoton/F");        
    ana_tree->Branch("gen_mass_dijet",      &gen_mass_dijet,      "gen_mass_dijet/F");         
    ana_tree->Branch("gen_pt_dijet",        &gen_pt_dijet,        "gen_pt_dijet/F");         
    ana_tree->Branch("gen_eta_dijet",       &gen_eta_dijet,       "gen_eta_dijet/F");         
    ana_tree->Branch("gen_phi_dijet",       &gen_phi_dijet,       "gen_phi_dijet/F");         
    ana_tree->Branch("gen_zeppenfeld",      &gen_zeppenfeld,      "gen_zeppenfeld/F");         


    if(doPDFweight){
        ana_tree->Branch("nWeightsPDF1",&nWeightsPDF1,"nWeightsPDF1/I");
        ana_tree->Branch("nWeightsPDF2",&nWeightsPDF2,"nWeightsPDF2/I");
        ana_tree->Branch("nWeightsPDF3",&nWeightsPDF3,"nWeightsPDF3/I");
        ana_tree->Branch("nWeightsPDF4",&nWeightsPDF4,"nWeightsPDF4/I");
        ana_tree->Branch("nWeightsPDF5",&nWeightsPDF5,"nWeightsPDF5/I");
        ana_tree->Branch("nWeightsPDF6",&nWeightsPDF6,"nWeightsPDF6/I");
        ana_tree->Branch("nWeightsPDF7",&nWeightsPDF7,"nWeightsPDF7/I");
        ana_tree->Branch("nWeightsPDF8",&nWeightsPDF8,"nWeightsPDF8/I");
        ana_tree->Branch("nWeightsPDF9",&nWeightsPDF9,"nWeightsPDF9/I");
        ana_tree->Branch("nWeightsPDF10",&nWeightsPDF10,"nWeightsPDF10/I");
        ana_tree->Branch("PDFweight1",&PDFweight1,"PDFweight1[nWeightsPDF1]/F");
        ana_tree->Branch("PDFweight2",&PDFweight2,"PDFweight2[nWeightsPDF2]/F");
        ana_tree->Branch("PDFweight3",&PDFweight3,"PDFweight3[nWeightsPDF3]/F");
        ana_tree->Branch("PDFweight4",&PDFweight4,"PDFweight4[nWeightsPDF4]/F");
        ana_tree->Branch("PDFweight5",&PDFweight5,"PDFweight5[nWeightsPDF5]/F");
        ana_tree->Branch("PDFweight6",&PDFweight6,"PDFweight6[nWeightsPDF6]/F");
        ana_tree->Branch("PDFweight7",&PDFweight7,"PDFweight7[nWeightsPDF7]/F");
        ana_tree->Branch("PDFweight8",&PDFweight8,"PDFweight8[nWeightsPDF8]/F");
        ana_tree->Branch("PDFweight9",&PDFweight9,"PDFweight9[nWeightsPDF9]/F");
        ana_tree->Branch("PDFweight10",&PDFweight10,"PDFweight10[nWeightsPDF10]/F");
    }


    /********************************************************
     *                                                      *
     *           SETTING PHOTON ID PARAMETERS               *
     *                                                      *
     ********************************************************/

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
    
    photonidcuts preselid; 
    preselid.hcaliso_rel=         1000; 
    preselid.hcaliso_abs=         8.; 
    preselid.ecaliso_rel=         1000; 
    preselid.ecaliso_abs=         10; 
    preselid.tracknb=             1000.; 
    preselid.trackiso_rel=        1000.; 
    preselid.sminmin=             0.60; 
    preselid.sminmin_min=         0.10; 
    preselid.smajmaj=             1000.;  
    
    photonidelecuts WP95id;
    WP95id.hovereisoEB=           0.15;
    WP95id.hcaliso_relEB=         0.12;
    WP95id.ecaliso_relEB=         2.0;
    WP95id.trackiso_relEB=        0.15;
    WP95id.setaetaEB=             0.01;
    WP95id.detaEB     =           0.007;
    WP95id.dphiEB     =           1000.;
    WP95id.minhitsEB  =           2.;
    WP95id.dcotEB     =           -1000.;
    WP95id.distEB     =           -1000.;
    WP95id.hovereisoEE=           0.07;
    WP95id.hcaliso_relEE=         0.05;
    WP95id.ecaliso_relEE=         0.06;
    WP95id.trackiso_relEE=        0.08;
    WP95id.setaetaEE=             0.03;
    WP95id.detaEE     =           0.01;
    WP95id.dphiEE     =           1000.;
    WP95id.minhitsEE  =           2.;
    WP95id.dcotEE     =           -1000;
    WP95id.distEE     =           -1000;
    
    photonidelecuts WP80id;
    WP80id.hovereisoEB=           0.04;
    WP80id.hcaliso_relEB=         0.10;
    WP80id.ecaliso_relEB=         0.07;
    WP80id.trackiso_relEB=        0.09;
    WP80id.setaetaEB=             0.01;
    WP80id.detaEB     =           0.004;
    WP80id.dphiEB     =           0.06;
    WP80id.minhitsEB  =           1.;
    WP80id.dcotEB     =           0.02;
    WP80id.distEB     =           0.02;
    WP80id.hovereisoEE=           0.025;
    WP80id.hcaliso_relEE=         0.025;
    WP80id.ecaliso_relEE=         0.05;
    WP80id.trackiso_relEE=        0.04;
    WP80id.setaetaEE=             0.03;
    WP80id.detaEE     =           0.007;
    WP80id.dphiEE     =           0.03;
    WP80id.minhitsEE  =           1.;
    WP80id.dcotEE     =           0.02;
    WP80id.distEE     =           0.02;
    
    photonidegcuts preselegid;
    preselegid.hovereiso=           0.15;
    preselegid.hcaliso_rel=         0.005;
    preselegid.hcaliso_abs=         10.;
    preselegid.ecaliso_rel=         0.012;
    preselegid.ecaliso_abs=         10.;
    preselegid.trackiso_rel=        0.002;
    preselegid.trackiso_abs=        10.;
    preselegid.setaetaEB=           0.017;
    preselegid.setaetaEE=           0.04;
    
    photonidegcuts looseegid;
    looseegid.hovereiso=           0.05;
    looseegid.hcaliso_rel=         0.0025;
    looseegid.hcaliso_abs=         2.2;
    looseegid.ecaliso_rel=         0.006;
    looseegid.ecaliso_abs=         4.2;
    looseegid.trackiso_rel=        0.001;
    looseegid.trackiso_abs=        3.5;
    looseegid.setaetaEB=           1000.;
    looseegid.setaetaEE=           1000.;
    
    photonidegcuts loose006egid;
    loose006egid.hovereiso=           0.05;
    loose006egid.hcaliso_rel=         0.0025;
    loose006egid.hcaliso_abs=         2.2;
    loose006egid.ecaliso_rel=         0.006;
    loose006egid.ecaliso_abs=         4.2;
    loose006egid.trackiso_rel=        0.001;
    loose006egid.trackiso_abs=        2.;
    loose006egid.setaetaEB=           0.0105;
    loose006egid.setaetaEE=           0.030;
    
    photonidegcuts tightegid;
    tightegid.hovereiso=           0.05;
    tightegid.hcaliso_rel=         0.0025;
    tightegid.hcaliso_abs=         2.2;
    tightegid.ecaliso_rel=         0.006;
    tightegid.ecaliso_abs=         4.2;
    tightegid.trackiso_rel=        0.001;
    tightegid.trackiso_abs=        2.;
    tightegid.setaetaEB=           0.013;
    tightegid.setaetaEE=           0.030;
    
    photonidegcuts hggtightid;
    hggtightid.hovereiso=           0.02;
    hggtightid.hcaliso_rel=         0.0025;
    hggtightid.hcaliso_abs=         2.;
    hggtightid.ecaliso_rel=         0.006;
    hggtightid.ecaliso_abs=         2.;
    hggtightid.trackiso_rel=        0.001;
    hggtightid.trackiso_abs=        1.5;
    hggtightid.setaetaEB=           0.010;
    hggtightid.setaetaEE=           0.028;
    
    photonidegcuts isemid;
    isemid.hovereiso=           1000.;
    isemid.hcaliso_rel=         0.0025;
    isemid.hcaliso_abs=         2.2;
    isemid.ecaliso_rel=         0.006;
    isemid.ecaliso_abs=         4.2;
    isemid.trackiso_rel=        1000.;
    isemid.trackiso_abs=        1000.;
    isemid.setaetaEB=             1000.;
    isemid.setaetaEE=             1000.;
  
   /********************************************************
    *                                                      *
    *            APPLYING CUTS IN CATEGORIES               *
    *                                                      *
    ********************************************************/
   
    for(int iLevel=0; iLevel<phoNCUTLEVELS; ++iLevel) 
    {
        float cic6_cuts_lead[phoNCUTS][phoCiC6NCATEGORIES];
        float cic6_cuts_sublead[phoNCUTS][phoCiC6NCATEGORIES];
        float cic4_cuts_lead[phoNCUTS][phoCiC4NCATEGORIES];
        float cic4_cuts_sublead[phoNCUTS][phoCiC4NCATEGORIES];
        SetPhotonCutsInCategories((phoCiCIDLevel)iLevel, &cic6_cuts_lead[0][0], &cic6_cuts_sublead[0][0], &cic4_cuts_lead[0][0], &cic4_cuts_sublead[0][0] );
     
        float * cic6_cuts_arrays_lead[phoNCUTS] = {
            &cic6_cut_lead_isosumoet[0][0], &cic6_cut_lead_isosumoetbad[0][0], &cic6_cut_lead_trkisooet[0][0], &cic6_cut_lead_sieie[0][0],
            &cic6_cut_lead_hovere[0][0], &cic6_cut_lead_r9[0][0], &cic6_cut_lead_drtotk_25_99[0][0], &cic6_cut_lead_pixel[0][0] 
        };
     
        float * cic6_cuts_arrays_sublead[phoNCUTS] = {
            &cic6_cut_sublead_isosumoet[0][0], &cic6_cut_sublead_isosumoetbad[0][0], &cic6_cut_sublead_trkisooet[0][0], 
            &cic6_cut_sublead_sieie[0][0], &cic6_cut_sublead_hovere[0][0], &cic6_cut_sublead_r9[0][0],
            &cic6_cut_sublead_drtotk_25_99[0][0], &cic6_cut_sublead_pixel[0][0]
        };
     
        float * cic4_cuts_arrays_lead[phoNCUTS] = {
            &cic4_cut_lead_isosumoet[0][0], &cic4_cut_lead_isosumoetbad[0][0], &cic4_cut_lead_trkisooet[0][0], &cic4_cut_lead_sieie[0][0],
            &cic4_cut_lead_hovere[0][0], &cic4_cut_lead_r9[0][0], &cic4_cut_lead_drtotk_25_99[0][0], &cic4_cut_lead_pixel[0][0] 
        } ;
     
        float * cic4_cuts_arrays_sublead[phoNCUTS] = {
            &cic4_cut_sublead_isosumoet[0][0], &cic4_cut_sublead_isosumoetbad[0][0], &cic4_cut_sublead_trkisooet[0][0], 
            &cic4_cut_sublead_sieie[0][0], &cic4_cut_sublead_hovere[0][0], &cic4_cut_sublead_r9[0][0],
            &cic4_cut_sublead_drtotk_25_99[0][0], &cic4_cut_sublead_pixel[0][0]
        };
     
        for(int iCut=0; iCut<phoNCUTS; ++iCut) {
            for(int iCat=0; iCat<phoCiC6NCATEGORIES; ++iCat) {
                cic6_cuts_arrays_lead[iCut][iLevel*phoCiC6NCATEGORIES+iCat] = cic6_cuts_lead[iCut][iCat];
                cic6_cuts_arrays_sublead[iCut][iLevel*phoCiC6NCATEGORIES+iCat] = cic6_cuts_sublead[iCut][iCat];
            }
            for(int iCat=0; iCat<phoCiC4NCATEGORIES; ++iCat) {
                cic4_cuts_arrays_lead[iCut][iLevel*phoCiC4NCATEGORIES+iCat] = cic4_cuts_lead[iCut][iCat];
                cic4_cuts_arrays_sublead[iCut][iLevel*phoCiC4NCATEGORIES+iCat] = cic4_cuts_sublead[iCut][iCat];
            }
        }
    }

   /********************************************************
    *                                                      *
    *                      CiC PLOTS                       *
    *                                                      *
    ********************************************************/

    for (int icat=0;icat<phoCiC4NCATEGORIES;++icat)
    {
        TString catName="cat";
        catName+=icat;
        catName+="_";
        
        cic4_cut_isosumoet[icat]=new TH1F("isosumoet_"+catName,"isosumoet_"+catName,100,-1.,25.);
        cic4_cut_isosumoetbad[icat]=new TH1F("isosumoetbad_"+catName,"isosumoetbad_"+catName,200,-1.,50.);
        cic4_cut_trkisooet[icat]=new TH1F("trkisooet_"+catName,"trkisooet_"+catName,100,0.,10.);
        cic4_cut_sieie[icat]=new TH1F("sieie_"+catName,"sieie_"+catName,200,0.,0.05);
        cic4_cut_hovere[icat]=new TH1F("hovere_"+catName,"hovere_"+catName,100,0.,0.1);
        cic4_cut_r9[icat]=new TH1F("r9_"+catName,"r9_"+catName,110,0.,1.1);
        cic4_cut_drtotk_25_99[icat]=new TH1F("drtotk_25_99_"+catName,"drtotk_25_99_"+catName,100,0.,0.5);
        cic4_cut_pixel[icat]=new TH1F("pixel_"+catName,"pixel_"+catName,20,-0.25,9.75);
    }


   /********************************************************
    *                                                      *
    *                   SETTING CUTS                       *
    *                                                      *
    ********************************************************/

    //event based cuts
    double ptphot1cut = 50;
    double ptphot2cut = 30;
    double ptjet1cut = 20;
    double ptjet2cut = 15;
    double deltaetacut = 2.5;
    double zeppencut = 2.5;
    double dijetmasscut = 300;
    double deltaphicut = 2.;

   // temp varables to ckeep track of the file being processed
   TString foldname("");
   TString currfilename("");
   int ifile(0);
   int nfiles = ((TChain*)fChain)->GetListOfFiles()->GetEntries();

   int nprocessed = 0;
   int nredntp = 0;
   timer.Start();


   /********************************************************
    *                                                      *
    *                       LOOP                           *
    *                                                      *
    ********************************************************/


    for (Long64_t jentry=0; jentry<nentries;jentry++) 
    {
        Long64_t ientry = LoadTree(jentry);
        if (ientry < 0) break;
        nb = fChain->GetEntry(jentry);   nbytes += nb;

        // if (Cut(ientry) < 0) continue;
        // json file event removal

        if (myjson && nMC<=0) 
	        if (!myjson->isGoodLS(run,lbn))
	        {
	            //	    std::cout << "Event skipped " << run << " " << lbn << std::endl;
	            continue;
	        }
    
        nprocessed++;
	
#ifdef DEBUG
        cout << "[DEBUG]" << endl;
        cout << "[DEBUG]" << endl;
#endif


        if (nprocessed%1000 == 0) cout << "Events " << nprocessed << " processed; Run " << run << " LS " << lbn << endl;
      
        if (scaleCorrections_)
	        correctPhotons(true);
      
        if (jetsyst_ && typejetsyst_>0 && typejetsyst_<5)
        {
            if(typejetsyst_ == 1) correctJets(1,0);
            if(typejetsyst_ == 2) correctJets(-1,0);
            if(typejetsyst_ == 3) correctJets(0,0.1);
            if(typejetsyst_ == 4) correctJets(0,-0.1);
        }

        // print name of crrent file
        currfilename = TString(fChain->GetCurrentFile()->GetName());
        if(currfilename != foldname) {
           ifile++;
           cout << "Opening file " << ifile << " of "  << nfiles << "\n"
                << currfilename  << "\n"
                << "------------------------------"
                << endl;
           foldname = currfilename;
        }
      
      // fill histos for PDF studies
        for(int iy=0; iy<nWeightsPDF[0] ; iy++)
            nPDFweight1.Fill(iy,pdfWeight[0][iy]);
        for(int iy=0; iy<nWeightsPDF[1] ; iy++)
            nPDFweight2.Fill(iy,pdfWeight[1][iy]);
        for(int iy=0; iy<nWeightsPDF[2] ; iy++)
            nPDFweight3.Fill(iy,pdfWeight[2][iy]);
        for(int iy=0; iy<nWeightsPDF[3] ; iy++)
            nPDFweight4.Fill(iy,pdfWeight[3][iy]);
        for(int iy=0; iy<nWeightsPDF[4] ; iy++)
            nPDFweight5.Fill(iy,pdfWeight[4][iy]);
        for(int iy=0; iy<nWeightsPDF[5] ; iy++)
            nPDFweight6.Fill(iy,pdfWeight[5][iy]);
        for(int iy=0; iy<nWeightsPDF[6] ; iy++)
            nPDFweight7.Fill(iy,pdfWeight[6][iy]);
        for(int iy=0; iy<nWeightsPDF[7] ; iy++)
            nPDFweight8.Fill(iy,pdfWeight[7][iy]);
        for(int iy=0; iy<nWeightsPDF[8] ; iy++)
            nPDFweight9.Fill(iy,pdfWeight[8][iy]);
        for(int iy=0; iy<nWeightsPDF[9] ; iy++)
            nPDFweight10.Fill(iy,pdfWeight[9][iy]);

        vector<bool> photassocMC, photassocMChiggs;
   
        int counter(0), countertt(0), ishiggsev(0);
        int isZH(0);
        int isWH(0);


   /********************************************************
    *                                                      *
    *                 LOOP :: GEN ANALYSIS                 *
    *                                                      *
    ********************************************************/

        /// init of mc related variables
        int higgsId=-1;
        int VHLeptonIndexes[2];
        int leptonCounter(0);
        int leptonIndex[10];
        double leptonPt[10];
        int iLep = 0;
        
        for(int i=0; i<nMC; i++)
        {      
            // cout << "pId:" << pdgIdMC[i] << "\tstatus:" << statusMC[i] << "\tmothIndex:" << motherIDMC[i];
            // if(motherIDMC[i]>=0 && motherIDMC[i]<nMC)
            //     cout << "\tmothId:" << pdgIdMC[motherIDMC[i]] << "\tmothStatus:"  << statusMC[motherIDMC[i]] << endl;
            // else
            //     cout << endl;

            if(pdgIdMC[i] == 25) 
            {
                ishiggsev=1;
                higgsId=i;
            }
            else if ( pdgIdMC[i] == 23 )
            {
                isZH = 1;
            }
            else if ( TMath::Abs(pdgIdMC[i])==24 )
            {
                isWH = 1;
            }

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
            
            if(TMath::Abs(pdgIdMC[i]) == 6 && TMath::Abs(pdgIdMC[motherIDMC[i]])<23)
                countertt++;
            
            /// if there are photons coming not from a photon or a higgs
            if(pdgIdMC[i] == 22 && statusMC[i] == 1 && TMath::Abs(pdgIdMC[motherIDMC[i]])<21)
                counter++;
        
            // leptonic VH events
            if( TMath::Abs(pdgIdMC[i]) <= 16 && TMath::Abs(pdgIdMC[i]) >= 11  && statusMC[i] == 3  &&
               (pdgIdMC[motherIDMC[i]] == 23 || TMath::Abs(pdgIdMC[motherIDMC[i]])==24 ) )
            {
                //cout << "p:" << pdgIdMC[i] << " status:" << statusMC[i] <<" pt:" << ptMC[i] << " eta:" << etaMC[i] << endl;
                //cout << "p:" << pdgIdMC[i] << " status:" << statusMC[i] <<" mothId:" << pdgIdMC[motherIDMC[i]] << " mothStatus:"  << statusMC[motherIDMC[i]] << endl;
                VHLeptonIndexes[leptonCounter] = i;
                leptonCounter++;
                // cout << endl <<"mother id:" << pdgIdMC[motherIDMC[i]] << " pt:" << ptMC[motherIDMC[i]] << " phi:" << phiMC[motherIDMC[i]] << endl;
            }

            // considering leptons not coming from the higgs
            // if( TMath::Abs(pdgIdMC[i]) <= 16 && TMath::Abs(pdgIdMC[i]) >= 11  && statusMC[i] == 3  &&
            //     TMath::Abs(pdgIdMC[i]) != 15 )
            // {
            //     if(iLep >=10 )
            //     {
            //         cout << "There are more than 10 leptons in the event. Skipping the others..." << endl;
            //         continue;
            //     }
            //     leptonIndex[iLep] = i;
            //     leptonPt[iLep] = ptMC[i];
            //     iLep++;
            // }
        }
   

       /***************************************************
        *                                                 *
        *           IDENTIFYING PHYSICS PROCESS           *
        *                                                 *
        ***************************************************/

        if(genProcessId == 10012) // GGF 
            gen_custom_processId = 1001;

        else if(genProcessId == 10001) // VBF
            gen_custom_processId = 2011;

        else if( ishiggsev && isZH  && genProcessId == 24)
        {
            if(leptonCounter==2)
            {
                if(TMath::Abs(pdgIdMC[VHLeptonIndexes[0]]) == 11) // electron
                    gen_custom_processId = 4101;

                else if(TMath::Abs(pdgIdMC[VHLeptonIndexes[0]]) == 13) // muon
                    gen_custom_processId = 4201;

                else if(TMath::Abs(pdgIdMC[VHLeptonIndexes[0]]) == 15) // tau
                    gen_custom_processId = 4301;

                else if(TMath::Abs(pdgIdMC[VHLeptonIndexes[0]]) == 12 || TMath::Abs(pdgIdMC[VHLeptonIndexes[0]]) == 14 || TMath::Abs(pdgIdMC[VHLeptonIndexes[0]]) == 16 ) // nu
                    gen_custom_processId = 4501;
            }
            else
                gen_custom_processId = 4401;
        }
        else if (ishiggsev && isWH && genProcessId == 26)
        {
            if(leptonCounter==2)
            {
                if(TMath::Abs(pdgIdMC[VHLeptonIndexes[0]]) == 11 || TMath::Abs(pdgIdMC[1]) == 11 ) // electron
                    gen_custom_processId = 3101;

                else if(TMath::Abs(pdgIdMC[VHLeptonIndexes[0]]) == 13 || TMath::Abs(pdgIdMC[1]) == 13 ) // muon
                    gen_custom_processId = 3201;

                else if(TMath::Abs(pdgIdMC[VHLeptonIndexes[0]]) == 15 || TMath::Abs(pdgIdMC[1]) == 15 ) // tau
                    gen_custom_processId = 3301;
            }
            else
                gen_custom_processId = 3401;
        }
        else
                gen_custom_processId = 9999;

        // cout << "gen_custom_processId : " << gen_custom_processId << endl; 
        // cout << endl << endl;


        if(isgjetqcd && counter > 1) continue; 
        //      To be used only when ttH is not produced separately  
        //      if(ishiggsev && countertt>0) continue; 
        
        vector<int> firstfourgenphot = firstfour(ptMC,&photassocMC);
        vector<int> firstfourhiggsgenphot = firstfour(ptMC,&photassocMChiggs);
        

#ifdef DEBUG
        cout << "[DEBUG] photassocMC.size() = " << photassocMC.size() << endl;
        cout << "[DEBUG] firstfourgenphot.size() = " <<  firstfourgenphot.size() << endl;
        cout << "[DEBUG] firstfourhiggsgenphot.size() = " <<  firstfourhiggsgenphot.size() << endl;
        cout << "[DEBUG] firstfourgenphot.at(0) = " << firstfourgenphot.at(0) << "  - pt = " << ptMC[firstfourgenphot.at(0) ] << endl;
        cout << "[DEBUG] firstfourgenphot.at(1) = " << firstfourgenphot.at(1) << "  - pt = " << ptMC[firstfourgenphot.at(1) ] << endl;
#endif


        npu = pu_n;
        if(npu<MAX_PU_REWEIGHT && puweights_.size()>0 && nMC>0) 
	        pu_weight = puweights_[npu];
        else
	        pu_weight = 1;

        //Pt Reweighting
        if (genProcessId==10012 && ptweights_!=0 && higgsId!=-1)
	    {
            //calculate bin size
            double binsize = (ptweights_->GetXaxis()->GetXmax()-ptweights_->GetXaxis()->GetXmin())/ptweights_->GetNbinsX();
            double higgspt = ptMC[higgsId];
            int bin = 0;

            // underflow protection: use underflow entry
            if(higgspt >= ptweights_->GetXaxis()->GetXmin()){
                bin = Int_t((higgspt-ptweights_->GetXaxis()->GetXmin())/binsize) + 1;
            }

            // overflow protection: use overflow entry
            // FIXME weights overflow bin seems to be 0. Not really using overflow but weight of last available bin
            if(bin > ptweights_->GetNbinsX()) bin=ptweights_->GetNbinsX();

            // std::cout <<" Bin Size "<< binsize <<std::endl;
            // std::cout <<" Higgs Pt "<< higgspt <<std::endl;
            // std::cout <<" Bin  "<< bin <<std::endl;
            // std::cout <<" KFactor "<<   ptweights_->GetBinContent(bin) <<std::endl;

	        // get KFactor
	        pt_weight=  ptweights_->GetBinContent(bin);
	        if (pt_weight==0)
	        {
	            std::cout <<"PTWEIGHT=0. THIS SHOULD NOT HAPPEN!!" << std::endl;
	            std::cout <<"Bin Size "<< binsize <<std::endl;
	            std::cout <<"Higgs Pt "<< higgspt <<std::endl;
	            std::cout <<"Bin  "<< bin <<std::endl;
	            std::cout <<"KFactor "<<   ptweights_->GetBinContent(bin) <<std::endl;
	        }
	    }
         else
	    {
	        pt_weight=1.;
	    }

        weight=pu_weight*pt_weight;
        
        npunorew.Fill(npu);
        npurew.Fill(npu,weight);
        nvtxnorew.Fill(nvertex);
        nvtxrew.Fill(nvertex,weight);
        
        ptphotgen1.Fill(ptMC[firstfourgenphot.at(0)],weight);
        ptphotgen2.Fill(ptMC[firstfourgenphot.at(1)],weight);
        etaphotgen1.Fill(etaMC[firstfourgenphot.at(0)],weight);
        etaphotgen2.Fill(etaMC[firstfourgenphot.at(1)],weight);
        
        ptphothiggsgen1.Fill(ptMC[firstfourhiggsgenphot.at(0)],weight);
        ptphothiggsgen2.Fill(ptMC[firstfourhiggsgenphot.at(1)],weight);
        etaphothiggsgen1.Fill(etaMC[firstfourhiggsgenphot.at(0)],weight);
        etaphothiggsgen2.Fill(etaMC[firstfourhiggsgenphot.at(1)],weight);
        
        
#ifdef DEBUG
        cout << "[DEBUG] before unprotected genjet" << endl;
#endif
        TLorentzVector jetgen1, jetgen2;	
        jetgen1.SetPtEtaPhiE(ptJetGen_akt5[0],etaJetGen_akt5[0],phiJetGen_akt5[0],eJetGen_akt5[0]);
        jetgen2.SetPtEtaPhiE(ptJetGen_akt5[1],etaJetGen_akt5[1],phiJetGen_akt5[1],eJetGen_akt5[1]);
        
        TLorentzVector sumgen = jetgen1 + jetgen2;
#ifdef DEBUG
        cout << "[DEBUG] after unprotected genjet" << endl;
#endif
        
        ptjetgen1.Fill(ptJetGen_akt5[0],weight);
        ptjetgen2.Fill(ptJetGen_akt5[1],weight);
        etajetgen1.Fill(etaJetGen_akt5[0],weight);
        etajetgen2.Fill(etaJetGen_akt5[1],weight);

        if(ptJetGen_akt5[0]>20 && ptJetGen_akt5[1]>20)
        {
            invmassjetgen.Fill(sumgen.M(),weight);
            deltaetajetgen.Fill(etaJetGen_akt5[0]-etaJetGen_akt5[1],weight);
        }

        if(etaJetGen_akt5[0]*etaJetGen_akt5[1]<0) 
            deltaetajetgencut.Fill(etaJetGen_akt5[0]-etaJetGen_akt5[1],weight);


       /***************************************************
        *                                                 *
        *   COMPUTING & FILLING TREE WITH GEN VARIABLES   *
        *                                                 *
        ***************************************************/


        /// sorting index arrays
        int ptJetGen_akt5_sortingIndex[NGENJETS];
        
        /// sort gen jets according to pt
        TMath::Sort(nJetGen_akt5, ptJetGen_akt5, ptJetGen_akt5_sortingIndex);
        
        /// check if analyzing signal or bkg MC sample
        vector<int>* genPhotPtr =  gen_custom_processId > 9000? &firstfourgenphot : &firstfourhiggsgenphot;
        int index_phot1 = genPhotPtr->at(0);
        int index_phot2 = genPhotPtr->at(1);
        
#ifdef DEBUG
        if( gen_custom_processId > 9000)
            cout << "[DEBUG] bkg process" << endl;
        else
            cout << "[DEBUG] sig process" << endl;
        

        cout << "[DEBUG] index_phot1 = " << index_phot1 << endl;
        cout << "[DEBUG] index_phot2 = " << index_phot2 << endl;
#endif


        bool genPreselection = index_phot1 >= 0 && index_phot2 >= 0 ;
        /// to avoid reading bad memory locations
        //genPreselection = genPreselection? ptMC[index_phot1] > 20. &&  ptMC[index_phot2] > 20. && TMath::Abs(etaMC[index_phot1]) < 3. && TMath::Abs(etaMC[index_phot2]) < 3 : 0; 
        genPreselection = genPreselection? ptMC[index_phot1] > 10. &&  ptMC[index_phot2] > 10. && TMath::Abs(etaMC[index_phot1]) < 4. && TMath::Abs(etaMC[index_phot2]) < 4 : 0; 


        
#ifdef DEBUG
        cout << "[DEBUG] genPreselection = " << genPreselection << endl;
#endif

        /// if there are good gen photons in the event
        if(!genPreselection)
        {
            SetAllGenVarToMinus999();

        }
        else
        {
            /// find isolated jets
            int foundJets = 0;
            int isoJetIndex[2];
            for(int ijet=0; ijet < nJetGen_akt5 && foundJets != 2 ; ++ijet)
            { 
                bool isIso(1);
                for(int kphot=0; kphot < 2 && foundJets != 2; ++kphot )
                    isIso &= ( sqrt( pow( delta_eta(etaJetGen_akt5[ptJetGen_akt5_sortingIndex[ijet]],etaMC[genPhotPtr->at(kphot)]),2 )  + 
                                     pow( delta_phi(phiJetGen_akt5[ptJetGen_akt5_sortingIndex[ijet]],phiMC[genPhotPtr->at(kphot)]),2 ) ) > 0.5 ) ;
            
                if(!isIso) continue;
                isoJetIndex[foundJets] = ptJetGen_akt5_sortingIndex[ijet];
                foundJets++;
            }
            
#ifdef DEBUG
            cout << "[DEBUG] before genPhot1/2 " << endl;
#endif
            TLorentzVector genPhot1, genPhot2;	

#ifdef DEBUG
            cout << "[DEBUG] genPhot1:: " << ptMC[index_phot1] << ", " << etaMC[index_phot1] << ", " << phiMC[index_phot1] << ", " << eMC[index_phot1] << endl;
            cout << "[DEBUG] genPhot2:: " << ptMC[index_phot2] << ", " << etaMC[index_phot2] << ", " << phiMC[index_phot2] << ", " << eMC[index_phot2] << endl;
#endif
            genPhot1.SetPtEtaPhiE( ptMC[index_phot1], etaMC[index_phot1], phiMC[index_phot1], eMC[index_phot1]);
            genPhot2.SetPtEtaPhiE( ptMC[index_phot2], etaMC[index_phot2], phiMC[index_phot2], eMC[index_phot2]);
            TLorentzVector diphot = genPhot1 + genPhot2;
            
#ifdef DEBUG
            cout << "[DEBUG] after genPhot1/2 " << endl;
            cout << "[DEBUG] before genJet1/2 " << endl;
#endif

            TLorentzVector genJet1, genJet2;
            genJet1.SetPtEtaPhiE( ptJetGen_akt5[isoJetIndex[0]],  etaJetGen_akt5[isoJetIndex[0]],  phiJetGen_akt5[isoJetIndex[0]],  eJetGen_akt5[isoJetIndex[0]]);
            
            genJet2.SetPtEtaPhiE( ptJetGen_akt5[isoJetIndex[1]], etaJetGen_akt5[isoJetIndex[1]], phiJetGen_akt5[isoJetIndex[1]], eJetGen_akt5[isoJetIndex[1]]);
            TLorentzVector dijet = genJet1 + genJet2;

#ifdef DEBUG
            cout << "[DEBUG] after genJet1/2 " << endl;
#endif
        
            double aveeta  = (etaJetGen_akt5[isoJetIndex[0]] + etaJetGen_akt5[isoJetIndex[1]])/2.;
            gen_zeppenfeld = diphot.Eta() - aveeta;

#ifdef DEBUG
            cout << "[DEBUG] after diphot.Eta() for zeppen computation " << endl;
#endif
            
            gen_pt_gamma1  = ptMC[index_phot1];
            gen_pt_gamma2  = ptMC[index_phot2];
            gen_eta_gamma1 = etaMC[index_phot1];
            gen_eta_gamma2 = etaMC[index_phot2] ;
            gen_phi_gamma1 = phiMC[index_phot1];
            gen_phi_gamma2 = phiMC[index_phot2] ;
            
            gen_pt_genjet1  =  ptJetGen_akt5[isoJetIndex[0]];
            gen_pt_genjet2  =  ptJetGen_akt5[isoJetIndex[1]];
            gen_eta_genjet1 =  etaJetGen_akt5[isoJetIndex[0]];
            gen_eta_genjet2 =  etaJetGen_akt5[isoJetIndex[1]];
            gen_phi_genjet1 =  phiJetGen_akt5[isoJetIndex[0]];
            gen_phi_genjet2 =  phiJetGen_akt5[isoJetIndex[1]];
            
            gen_mass_diphoton  = diphot.M();
            gen_pt_diphoton    = diphot.Pt();
            gen_eta_diphoton   = diphot.Eta();
            gen_phi_diphoton   = diphot.Phi();
                
            gen_mass_dijet = dijet.M();
            gen_pt_dijet   = dijet.Pt();
            gen_eta_dijet  = dijet.Eta();
            gen_phi_dijet  = dijet.Phi();
        }



        // skip events where the number of jets, photons, and vertexes is above the maximum allowed value
        if (nPhot>30) {
	        cout << "number of photons = " << nPhot << " and above threshold of 30; skipping" << endl;
	        continue;
        }

        if (nJet_akt5 > 200) {
	        cout << "number of nJet_akt5 = " << nJet_akt5 << " and above threshold of 200; skipping" << endl;
	        continue;
        }
        if (nJet_akt7 > 200) {
	        cout << "number of nJet_akt7 = " << nJet_akt7 << " and above threshold of 200; skipping" << endl;
	        continue;
        }
        if (nJet_pfakt5 > 200) {
	        cout << "number of nJet_pfakt5 = " << nJet_pfakt5 << " and above threshold of 200; skipping" << endl;
	        continue;
        }
        if (nJet_pfakt7 > 200) {
        	cout << "number of nJet_pfakt7 = " << nJet_pfakt7 << " and above threshold of 200; skipping" << endl;
	        continue;
        }
        if (nJetGen_akt5 > 200) {
	        cout << "number of nJetGen_akt5 = " << nJetGen_akt5 << " and above threshold of 200; skipping" << endl;
        	continue;
        }
        if (nJetGen_akt7 > 200) {
	        cout << "number of nJetGen_akt7 = " << nJetGen_akt7 << " and above threshold of 200; skipping" << endl;
        	continue;
        }
        if (nvertex > MAX_PU_REWEIGHT) {
        	cout << "number of nvertex = " << nvertex << " and above threshold of " << MAX_PU_REWEIGHT << "; skipping" << endl;
        	continue;
        }

        vector<bool> assophothiggs;
        vector<bool> assophot;
        vector<bool> jetphot;
        vector<bool> isophot;
        vector<bool> isophotele;
        vector<bool> isophotloose;
        vector<bool> isophotmedium;
        vector<bool> isophotloosecs; 
        vector<bool> isophotmediumcs; 
        vector<bool> isophotemeg;
        vector<bool> isophotlooseeg;
        vector<bool> isophotloose006eg;
        vector<bool> isophottighteg;
        vector<bool> isophothggtight;
        vector<bool> isophotloosepueg;
        vector<bool> isophottightpueg;
        vector<bool> isophothggtightpu;
        vector<int>  isocic;
        vector<int>  isocicnoelveto;
        
#ifdef DEBUG
        cout << "[DEBUG] before reco vector declaration" << endl;
#endif
        TLorentzVector thehiggs;
        TLorentzVector thehiggsnewvtx;
        TLorentzVector thejet1;
        TLorentzVector thejet2;
#ifdef DEBUG
        cout << "[DEBUG] after reco vector declaration" << endl;
#endif


       /***************************************************
        *                                                 *
        *                 RECO PHOTONS                    *
        *                                                 *
        ***************************************************/


#ifdef DEBUG
        cout << "[DEBUG] nPhot = " << nPhot << endl;
#endif

        for(int i=0; i<nPhot; i++)
        {
  
        bool assh(0);
        bool assp(0);
        bool assj(0);
        bool assjmc(0);
        
        // TEMP!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        //rhoPF = 0;
        // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        
            /// montecarlo association
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

            /// isolation from jets
            for(int j=0; j<nJetGen_akt5; j++)
            {
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
            vector<bool> idpassele(9);
            vector<bool> idpasseg(5);
            
            Dvz.Fill(vz[0]-vzMC);
            Dvzbest.Fill(vz[vrankPhotonPairs[0]]-vzMC);	
            
            if (ptPhot[i]>25. && assh)
              FillPhotonCiCSelectionVariable(i,vrankPhotonPairs[0]);

            // photon id used for preselection
            string finder(selection);
            bool preselection;
            
            if (finder == "superloose") preselection = cutID(i, superlooseid, &idpass);
            else if (finder == "loose") preselection = cutID(i, looseid, &idpass);
            else if (finder == "medium") preselection = cutID(i, mediumid, &idpass);
            else if (finder == "isem") preselection = cutIDEG(i, isemid, &idpasseg);
            else if (finder == "looseeg") preselection = cutIDEG(i, looseegid, &idpasseg);
            else if (finder == "tighteg") preselection = cutIDEG(i, tightegid, &idpasseg);
            else if (finder == "hggtighteg") preselection = cutIDEG(i, hggtightid, &idpasseg);
            else if (finder == "preselection") preselection = cutIDpresel(i, preselid, &idpass);
            else if (finder == "preselectionCS") preselection = cutIDEG(i, preselegid, &idpasseg);
            else if (finder == "looseegpu") preselection = cutIDEG(i, looseegid, &idpasseg,1);
            else if (finder == "tightegpu") preselection = cutIDEG(i, tightegid, &idpasseg,1);
            else if (finder == "hggtightegpu") preselection = cutIDEG(i, hggtightid, &idpasseg,1);
            else if (finder == "cicloose") preselection = PhotonCiCSelectionLevel(i,1,vrankPhotonPairs[0]) >= 1;	
            else if (finder == "cicloosenoeleveto") preselection = PhotonCiCSelectionLevel(i,0,vrankPhotonPairs[0]) >= 1;	
            else if (finder == "cicmedium") preselection = PhotonCiCSelectionLevel(i,1,vrankPhotonPairs[0]) >= 2;	
            else if (finder == "cictight") preselection = PhotonCiCSelectionLevel(i,1,vrankPhotonPairs[0]) >= 3;	
            else if (finder == "cicsuper") preselection = PhotonCiCSelectionLevel(i,1,vrankPhotonPairs[0]) >= 4;	
            else if (finder == "cichyper") preselection = PhotonCiCSelectionLevel(i,1,vrankPhotonPairs[0]) >= 5;	
            else if (finder == "mcass") preselection = mcID(i);
            else {
              cout << "NO SUCH " << selection << " PRESELECTION  AVAILABLE!!" << endl;
                  cout << "Good options are: superloose loose medium isem looseeg tighteg hggtighteg looseegpu tightegpu hggtightegpu preselection preselectionCS cicloose cicmedium cictight cicsuper cichyper mcass" << endl;
                  cout << "now exiting" << endl;
              exit(-1);
            }
            if(preselection) isophot.push_back(1); 
                else isophot.push_back(0);  
            
            if(cutIDele(i, WP80id, &idpassele)) isophotele.push_back(1); 
                else isophotele.push_back(0);  
            
            if(cutID(i, looseid, &idpass)) isophotloose.push_back(1); 
                else isophotloose.push_back(0);  
            
            if(cutID(i, mediumid, &idpass)) isophotmedium.push_back(1); 
                else isophotmedium.push_back(0);  
            
                if(cutIDcs(i, looseid, &idpass)) isophotloosecs.push_back(1);  
                else isophotloosecs.push_back(0);   
            
                if(cutIDcs(i, mediumid, &idpass)) isophotmediumcs.push_back(1);  
                else isophotmediumcs.push_back(0);   
            
            if(cutIDEG(i, isemid, &idpasseg)) isophotemeg.push_back(1); 
                else isophotemeg.push_back(0);  
            
            if(cutIDEG(i, looseegid, &idpasseg)) isophotlooseeg.push_back(1); 
                else isophotlooseeg.push_back(0);  
            
            if(cutIDEG(i, loose006egid, &idpasseg)) isophotloose006eg.push_back(1); 
                else isophotloose006eg.push_back(0);  
            
            if(cutIDEG(i, tightegid, &idpasseg)) isophottighteg.push_back(1); 
                else isophottighteg.push_back(0);  
            
            if(cutIDEG(i, hggtightid, &idpasseg)) isophothggtight.push_back(1); 
                else isophothggtight.push_back(0);  
            
            if(cutIDEG(i, looseegid, &idpasseg, 1)) isophotloosepueg.push_back(1); 
            else isophotloosepueg.push_back(0);  
            
            if(cutIDEG(i, tightegid, &idpasseg, 1)) isophottightpueg.push_back(1); 
            else isophottightpueg.push_back(0);  
            
            if(cutIDEG(i, hggtightid, &idpasseg, 1)) isophothggtightpu.push_back(1); 
            else isophothggtightpu.push_back(0);  
            
            isocic.push_back(PhotonCiCSelectionLevel(i,1,vrankPhotonPairs[0]));
            isocicnoelveto.push_back(PhotonCiCSelectionLevel(i,0,vrankPhotonPairs[0]));

            if( assp )	{
                ptphotassreco.Fill(ptPhot[i],weight);
                etaphotassreco.Fill(etaPhot[i],weight);
                if(ptPhot[i]>30) {
                    if(TMath::Abs(etascPhot[i])<1.47){
                        hcalisoassphot_EB.Fill(hcalovecal04Phot[i],weight);
                        ecalisoassphot_EB.Fill(ecaliso04Phot[i] / ePhot[i],weight);
                        ptisoassphot_EB.Fill(ptiso035Phot[i] / ptPhot[i],weight);
                        ntrkisoassphot_EB.Fill(ntrkiso035Phot[i],weight);
                        sminminclusassphot_EB.Fill(sMinMinPhot[i],weight);
                        smaxmaxclusassphot_EB.Fill(sMajMajPhot[i],weight);
                        alphaclusassphot_EB.Fill(alphaPhot[i],weight);
                    }
                    else if(TMath::Abs(etascPhot[i])<2.5){
                        hcalisoassphot_EE.Fill(hcalovecal04Phot[i],weight);
                        ecalisoassphot_EE.Fill(ecaliso04Phot[i] / ePhot[i],weight);
                        ptisoassphot_EE.Fill(ptiso035Phot[i] / ptPhot[i],weight);
                        ntrkisoassphot_EE.Fill(ntrkiso035Phot[i],weight);
                        sminminclusassphot_EE.Fill(sMinMinPhot[i],weight);
                        smaxmaxclusassphot_EE.Fill(sMajMajPhot[i],weight);	    
                        alphaclusassphot_EE.Fill(alphaPhot[i],weight);
                    }
                }
            }
            if( isophot.at(i) )	{
                ptphotisoreco.Fill(ptPhot[i],weight);
                etaphotisoreco.Fill(etaPhot[i],weight);
            }
            if( assp && isophot.at(i) )	{
                ptphotisoassreco.Fill(ptPhot[i],weight);
                etaphotisoassreco.Fill(etaPhot[i],weight);
            }
            if( !assp )	{
                ptphotnotassreco.Fill(ptPhot[i],weight);
                etaphotnotassreco.Fill(etaPhot[i],weight);
            }
	        if( !assp && isophot.at(i) )	{
	            ptphotisonotassreco.Fill(ptPhot[i],weight);
	            etaphotisonotassreco.Fill(etaPhot[i],weight);
	        }
	        if( !assp )	{
	            ptphotjetreco.Fill(ptPhot[i],weight);
	            etaphotjetreco.Fill(etaPhot[i],weight);
	            if(ptPhot[i]>30) {
	                if(TMath::Abs(etascPhot[i])<1.47){
	                    hcalisoassjet_EB.Fill(hcalovecal04Phot[i],weight);
	                    ecalisoassjet_EB.Fill(ecaliso04Phot[i] / ePhot[i],weight);
	                    ptisoassjet_EB.Fill(ptiso035Phot[i] / ptPhot[i],weight);
	                    ntrkisoassjet_EB.Fill(ntrkiso035Phot[i],weight);
	                    sminminclusassjet_EB.Fill(sMinMinPhot[i],weight);
	                    smaxmaxclusassjet_EB.Fill(sMajMajPhot[i],weight);
	                    alphaclusassjet_EB.Fill(alphaPhot[i],weight);
	                }
                    else if(TMath::Abs(etascPhot[i])<2.5){
	                    hcalisoassjet_EE.Fill(hcalovecal04Phot[i],weight);
	                    ecalisoassjet_EE.Fill(ecaliso04Phot[i] / ePhot[i],weight);
	                    ptisoassjet_EE.Fill(ptiso035Phot[i] / ptPhot[i],weight);
	                    ntrkisoassjet_EE.Fill(ntrkiso035Phot[i],weight);
	                    sminminclusassjet_EE.Fill(sMinMinPhot[i],weight);
	                    smaxmaxclusassjet_EE.Fill(sMajMajPhot[i],weight);	    
	                    alphaclusassjet_EE.Fill(alphaPhot[i],weight);
	                }
                }
	        }
	        if( assj && isophot.at(i) )	{
	            ptphotisojetreco.Fill(ptPhot[i],weight);
	            etaphotisojetreco.Fill(etaPhot[i],weight);
	        }
        }  

        vector<int> firstfourhiggsassphot = firstfour(ptPhot,&assophothiggs);
        vector<int> firstfourassphot = firstfour(ptPhot,&assophot);
        vector<int> firstfourisophot = firstfour(ptPhot,&isophot);      
        
        vector<bool> jetnohiggsphot;
        vector<bool> jetnoassphot;
        vector<bool> jetnoisophot;


       /***************************************************
        *                                                 *
        *                    RECO JETS                    *
        *                                                 *
        ***************************************************/


        for(int i=0; i<nJet_pfakt5; i++){

            bool assh(0);
            bool assp(0);
            bool assi(0);
            
            double DR;
            
            for(int k=0; k<2; k++){
              
                DR = sqrt(delta_eta(etaJet_pfakt5[i],etaPhot[firstfourhiggsassphot.at(k)])*delta_eta(etaJet_pfakt5[i],etaPhot[firstfourhiggsassphot.at(k)]) + 
                      delta_phi(phiJet_pfakt5[i],phiPhot[firstfourhiggsassphot.at(k)])*delta_phi(phiJet_pfakt5[i],phiPhot[firstfourhiggsassphot.at(k)]) ) ;
                if( DR < .5 ) assh = 1; 
                
                DR = sqrt(delta_eta(etaJet_pfakt5[i],etaPhot[firstfourassphot.at(k)])*delta_eta(etaJet_pfakt5[i],etaPhot[firstfourassphot.at(k)]) + 
                      delta_phi(phiJet_pfakt5[i],phiPhot[firstfourassphot.at(k)])*delta_phi(phiJet_pfakt5[i],phiPhot[firstfourassphot.at(k)]) ) ;
                if( DR < .5 ) assp = 1; 
                
                DR = sqrt(delta_eta(etaJet_pfakt5[i],etaPhot[firstfourisophot.at(k)])*delta_eta(etaJet_pfakt5[i],etaPhot[firstfourisophot.at(k)]) + 
                      delta_phi(phiJet_pfakt5[i],phiPhot[firstfourisophot.at(k)])*delta_phi(phiJet_pfakt5[i],phiPhot[firstfourisophot.at(k)]) ) ;
                if( DR < .5 ) assi = 1; 
            }
	
            bool goodetajet(1);
            
            if(TMath::Abs(etaJet_pfakt5[i]) > 4.7) goodetajet = 0;  
            if(TMath::Abs(etaJet_pfakt5[i]) < 2.0) 
              if( betaStar_pfakt5[i][vrankPhotonPairs[0]] > 0.2 * log( nvertex - 0.67 ) ) goodetajet = 0;
            
            if(!assh && goodetajet) jetnohiggsphot.push_back(1);
            else jetnohiggsphot.push_back(0); 
            
            if(!assp && goodetajet) jetnoassphot.push_back(1); 
            else jetnoassphot.push_back(0);  
            
            if(!assi && goodetajet) jetnoisophot.push_back(1); 
            else jetnoisophot.push_back(0);  
          
        }

        vector<int> firstfournohiggsjet = firstfour(ptCorrJet_pfakt5,&jetnohiggsphot);
        vector<int> firstfournoassjet = firstfour(ptCorrJet_pfakt5,&jetnoassphot);
        vector<int> firstfournoisojet = firstfour(ptCorrJet_pfakt5,&jetnoisophot);      
        
        if( firstfourhiggsassphot.at(0)>-1 && firstfourhiggsassphot.at(1)>-1 ) 
        { 
            TLorentzVector phot1, phot2;	
            phot1.SetPtEtaPhiE(ptPhot[firstfourhiggsassphot.at(0)],etaPhot[firstfourhiggsassphot.at(0)],phiPhot[firstfourhiggsassphot.at(0)],ePhot[firstfourhiggsassphot.at(0)]);
            phot2.SetPtEtaPhiE(ptPhot[firstfourhiggsassphot.at(1)],etaPhot[firstfourhiggsassphot.at(1)],phiPhot[firstfourhiggsassphot.at(1)],ePhot[firstfourhiggsassphot.at(1)]);
            
            TLorentzVector higgs = phot1 + phot2;
            
            higgsmasshiggsassreco.Fill(higgs.M(),weight);
            pthiggshiggsassreco.Fill(higgs.Pt(),weight);
            
            ptphothiggsassreco1.Fill(phot1.Pt(),weight);
            ptphothiggsassreco2.Fill(phot2.Pt(),weight);
            etaphothiggsassreco1.Fill(etaPhot[firstfourhiggsassphot.at(0)],weight);
            etaphothiggsassreco2.Fill(etaPhot[firstfourhiggsassphot.at(1)],weight);

        }

        double higgsmass(0), etahiggs(-999);

        if( firstfourassphot.at(0)>-1 && firstfourassphot.at(1)>-1 ) { 

            TLorentzVector phot1, phot2;	
            phot1.SetPtEtaPhiE(ptPhot[firstfourassphot.at(0)],etaPhot[firstfourassphot.at(0)],phiPhot[firstfourassphot.at(0)],ePhot[firstfourassphot.at(0)]);
            phot2.SetPtEtaPhiE(ptPhot[firstfourassphot.at(1)],etaPhot[firstfourassphot.at(1)],phiPhot[firstfourassphot.at(1)],ePhot[firstfourassphot.at(1)]);
            
            TLorentzVector higgs = phot1 + phot2;
            
            higgsmass = higgs.M();
            etahiggs = higgs.Eta();
            
            higgsmassassreco.Fill(higgs.M(),weight);
            pthiggsassreco.Fill(higgs.Pt(),weight);
            
            ptphotassreco1.Fill(phot1.Pt(),weight);
            ptphotassreco2.Fill(phot2.Pt(),weight);
            etaphotassreco1.Fill(etaPhot[firstfourassphot.at(0)],weight);
            etaphotassreco2.Fill(etaPhot[firstfourassphot.at(1)],weight);

        }	

        double higgsisomass(0), phihiggsiso(-999), etahiggsiso(-999), higgspt(-999.);
        double higgsisomassnewvtx(0), phihiggsisonewvtx(-999), etahiggsisonewvtx(-999), higgsptnewvtx(-999.);

        if( firstfourisophot.at(0)>-1 && firstfourisophot.at(1)>-1 ) 
        { 
#ifdef DEBUG
            cout << "[DEBUG]--- before reco photons" << endl;
#endif
            TLorentzVector phot1, phot2, phot1new, phot2new;	
            phot1.SetPtEtaPhiE(ptPhot[firstfourisophot.at(0)],etaPhot[firstfourisophot.at(0)],phiPhot[firstfourisophot.at(0)],ePhot[firstfourisophot.at(0)]);
            phot2.SetPtEtaPhiE(ptPhot[firstfourisophot.at(1)],etaPhot[firstfourisophot.at(1)],phiPhot[firstfourisophot.at(1)],ePhot[firstfourisophot.at(1)]);
            
            TLorentzVector higgs = phot1 + phot2;
            thehiggs = phot1 + phot2;
#ifdef DEBUG
            cout << "[DEBUG]--- after reco photons" << endl;
#endif
            
            higgsisomass = higgs.M();
            etahiggsiso = higgs.Eta();
            phihiggsiso = higgs.Phi();
                higgspt = higgs.Pt();
            
            higgsmassisoreco.Fill(higgs.M(),weight);
            higgsmassisorecofull.Fill(higgs.M(),weight);
            pthiggsisoreco.Fill(higgs.Pt(),weight);
            
            ptphotisoreco1.Fill(phot1.Pt(),weight);
            ptphotisoreco2.Fill(phot2.Pt(),weight);
            etaphotisoreco1.Fill(etaPhot[firstfourisophot.at(0)],weight);
            etaphotisoreco2.Fill(etaPhot[firstfourisophot.at(1)],weight);
            
            // recalculate photon kin with best vtx
            double xnew1 = xscPhot[firstfourisophot.at(0)] - vx[vrankPhotonPairs[0]];
            double ynew1 = yscPhot[firstfourisophot.at(0)] - vy[vrankPhotonPairs[0]];
            double znew1 = zscPhot[firstfourisophot.at(0)] - vz[vrankPhotonPairs[0]];
            double xnew2 = xscPhot[firstfourisophot.at(1)] - vx[vrankPhotonPairs[0]];
            double ynew2 = yscPhot[firstfourisophot.at(1)] - vy[vrankPhotonPairs[0]];
            double znew2 = zscPhot[firstfourisophot.at(1)] - vz[vrankPhotonPairs[0]];
            

#ifdef DEBUG
            cout << "[DEBUG]--- before newvtx photons" << endl;
#endif
            phot1new.SetX(xnew1); phot1new.SetY(ynew1); phot1new.SetZ(znew1);  phot1new.SetRho(ePhot[firstfourisophot.at(0)]);  phot1new.SetE(ePhot[firstfourisophot.at(0)]);
            phot2new.SetX(xnew2); phot2new.SetY(ynew2); phot2new.SetZ(znew2);  phot2new.SetRho(ePhot[firstfourisophot.at(1)]);  phot2new.SetE(ePhot[firstfourisophot.at(1)]);
            
            TLorentzVector higgsnew = phot1new + phot2new;
            thehiggsnewvtx = phot1new + phot2new;
#ifdef DEBUG
            cout << "[DEBUG]--- after newvtx photons" << endl;
#endif
            
            higgsisomassnewvtx = higgsnew.M();
            etahiggsisonewvtx = higgsnew.Eta();
            phihiggsisonewvtx = higgsnew.Phi();
            higgsptnewvtx = higgsnew.Pt();
        }	

        if( firstfourhiggsassphot.at(0)>-1 && firstfourhiggsassphot.at(1)>-1 ) 
        { 
            if( firstfournohiggsjet.at(0) > -1) {
              ptjethiggsassreco1.Fill(ptCorrJet_pfakt5[firstfournohiggsjet.at(0)],weight);
              etajethiggsassreco1.Fill(etaJet_pfakt5[firstfournohiggsjet.at(0)],weight);
            }
            if( firstfournohiggsjet.at(1) > -1) {
              ptjethiggsassreco2.Fill(ptCorrJet_pfakt5[firstfournohiggsjet.at(1)],weight);
              etajethiggsassreco2.Fill(etaJet_pfakt5[firstfournohiggsjet.at(1)],weight);
            }
            if( firstfournohiggsjet.at(0) > -1 && firstfournohiggsjet.at(1) > -1) {
              deltaetajethiggsassreco.Fill(etaJet_pfakt5[firstfournohiggsjet.at(0)]-etaJet_pfakt5[firstfournohiggsjet.at(1)]);
              double aveeta = (etaJet_pfakt5[firstfournohiggsjet.at(0)]+etaJet_pfakt5[firstfournohiggsjet.at(1)])/2;
              double zeppen1 = etaJet_pfakt5[firstfournohiggsjet.at(0)] - aveeta;
              double zeppen2 = etaJet_pfakt5[firstfournohiggsjet.at(1)] - aveeta;
              zeppenjethiggsassreco1.Fill(zeppen1,weight);
              zeppenjethiggsassreco2.Fill(zeppen2,weight);	
            }
        }
   
        double twojetsmass(0), etatwojets(-999), phitwojets(-999);

        if( firstfourassphot.at(0)>-1 && firstfourassphot.at(1)>-1 ) { 

            if( firstfournoassjet.at(0) > -1) {
              ptjetassreco1.Fill(ptCorrJet_pfakt5[firstfournoassjet.at(0)],weight);
              etajetassreco1.Fill(etaJet_pfakt5[firstfournoassjet.at(0)],weight);
            }
            if( firstfournoassjet.at(1) > -1) {
              ptjetassreco2.Fill(ptCorrJet_pfakt5[firstfournoassjet.at(1)],weight);
              etajetassreco2.Fill(etaJet_pfakt5[firstfournoassjet.at(1)]),weight;
            }
            if( firstfournoassjet.at(0) > -1 && firstfournoassjet.at(1) > -1) {
              TLorentzVector jet1, jet2;	
              jet1.SetPtEtaPhiE(ptCorrJet_pfakt5[firstfournoassjet.at(0)],etaJet_pfakt5[firstfournoassjet.at(0)],phiJet_pfakt5[firstfournoassjet.at(0)],eJet_pfakt5[firstfournoassjet.at(0)]/ptJet_pfakt5[firstfournoassjet.at(0)]*ptCorrJet_pfakt5[firstfournoassjet.at(0)]);
              jet2.SetPtEtaPhiE(ptCorrJet_pfakt5[firstfournoassjet.at(1)],etaJet_pfakt5[firstfournoassjet.at(1)],phiJet_pfakt5[firstfournoassjet.at(1)],eJet_pfakt5[firstfournoassjet.at(1)]/ptJet_pfakt5[firstfournoassjet.at(1)]*ptCorrJet_pfakt5[firstfournoassjet.at(1)]);
              
              TLorentzVector sum = jet1 + jet2;
              
              twojetsmass = sum.M();
              etatwojets = sum.Eta();
             
              deltaetajetassreco.Fill(etaJet_pfakt5[firstfournoassjet.at(0)]-etaJet_pfakt5[firstfournoassjet.at(1)],weight);
              double aveeta = (etaJet_pfakt5[firstfournoassjet.at(0)]+etaJet_pfakt5[firstfournoassjet.at(1)])/2;
              double zeppen1 = etaJet_pfakt5[firstfournoassjet.at(0)] - aveeta;
              double zeppen2 = etaJet_pfakt5[firstfournoassjet.at(1)] - aveeta;
              zeppenjetassreco1.Fill(zeppen1,weight);
              zeppenjetassreco2.Fill(zeppen2,weight);	
            }

            if(  ptCorrJet_pfakt5[firstfournoassjet.at(0)] > ptjet1cut && ptCorrJet_pfakt5[firstfournoassjet.at(1)] > ptjet2cut 
                && ptPhot[firstfourassphot.at(0)] > ptphot1cut && ptPhot[firstfourassphot.at(1)] > ptphot2cut )
            {
                if(TMath::Abs(etaJet_pfakt5[firstfournoassjet.at(0)]-etaJet_pfakt5[firstfournoassjet.at(1)])>deltaetacut)
                {
                    higgsmasscutreco.Fill(higgsmass,weight);
                    if(isophot.at(firstfourassphot.at(0)) && isophot.at(firstfourassphot.at(1)))  
                        higgsmassisorecocheck.Fill(higgsisomass,weight);	    

                    double aveeta = (etaJet_pfakt5[firstfournoassjet.at(0)]+etaJet_pfakt5[firstfournoassjet.at(1)])/2;
                    double zeppen = etahiggs - aveeta;
                    zeppenhiggsassreco.Fill(zeppen,weight);

                    if(TMath::Abs(zeppen)<zeppencut) 
                    {
                        higgsmasscutzeppreco.Fill(higgsmass,weight);
                        dijetmassassreco.Fill(twojetsmass,weight);
                        if(twojetsmass>250)
                            higgsmasscutzeppdijetreco.Fill(higgsmass,weight);
                    }
                }
            }

        }


       /***************************************************
        *                                                 *
        *        SAVING RECO VARIABLES IN TTREE           *
        *                                                 *
        ***************************************************/


      double twojetsmassiso(0), etatwojetsiso(-999), phitwojetsiso(-999);
   
      //bool recoPreselection = (firsttwoisophot.at(0)>-1 && firsttwoisophot.at(1)>-1 && ptPhot[firsttwoisophot.at(0)]>20 && ptPhot[firsttwoisophot.at(1)]>20);

      /// firstTWO --> firstFOUR
      bool recoPreselection = ( firstfourisophot.at(0)>-1 && firstfourisophot.at(1)>-1 && ptPhot[firstfourisophot.at(0)]>20 && ptPhot[firstfourisophot.at(1)]>20); 

#ifdef DEBUG
      cout << "[DEBUG] recoPreselection = " << recoPreselection << endl;
#endif
      
      if(!recoPreselection ) { 
        SetAllRecoVarToMinus999();
      }
      else {
	
	if( firstfournoisojet.at(0) > -1) {
	  ptjetisoreco1.Fill(ptCorrJet_pfakt5[firstfournoisojet.at(0)],weight);
	  etajetisoreco1.Fill(etaJet_pfakt5[firstfournoisojet.at(0)],weight);
	}
	if( firstfournoisojet.at(1) > -1) {
	  ptjetisoreco2.Fill(ptCorrJet_pfakt5[firstfournoisojet.at(1)],weight);
	  etajetisoreco2.Fill(etaJet_pfakt5[firstfournoisojet.at(1)],weight);
	}
	if( firstfournoisojet.at(0) > -1 && firstfournoisojet.at(1) > -1) {
	  

#ifdef DEBUG
        cout << "[DEBUG] before recojets" << endl;
#endif

	  TLorentzVector jet1, jet2;	
	  jet1.SetPtEtaPhiE(ptCorrJet_pfakt5[firstfournoisojet.at(0)],etaJet_pfakt5[firstfournoisojet.at(0)],phiJet_pfakt5[firstfournoisojet.at(0)],eJet_pfakt5[firstfournoisojet.at(0)]/ptJet_pfakt5[firstfournoisojet.at(0)]*ptCorrJet_pfakt5[firstfournoisojet.at(0)]);
	  jet2.SetPtEtaPhiE(ptCorrJet_pfakt5[firstfournoisojet.at(1)],etaJet_pfakt5[firstfournoisojet.at(1)],phiJet_pfakt5[firstfournoisojet.at(1)],eJet_pfakt5[firstfournoisojet.at(1)]/ptJet_pfakt5[firstfournoisojet.at(1)]*ptCorrJet_pfakt5[firstfournoisojet.at(1)]);
	  
	  TLorentzVector sum = jet1 + jet2;
	  thejet1 = jet1;
	  thejet2 = jet2;
#ifdef DEBUG
        cout << "[DEBUG] after recojets" << endl;
#endif

	  twojetsmassiso = sum.M();
	  etatwojetsiso = sum.Eta();
	  phitwojetsiso = sum.Phi();
	  
	  if(jetsyst_) 
	    JECunc.Fill(ptCorrJet_pfakt5[firstfournoisojet.at(0)],jetsyst_->getJESUncertainty(etaJet_pfakt5[firstfournoisojet.at(0)],ptCorrJet_pfakt5[firstfournoisojet.at(0)]));

	  int assjj(-999);
	  for(int j=0; j<nJetGen_akt5; j++){	
	    double DR = sqrt(delta_eta(etaJet_pfakt5[firstfournoisojet.at(0)],etaJetGen_akt5[j])*delta_eta(etaJet_pfakt5[firstfournoisojet.at(0)],etaJetGen_akt5[j]) + 
			     delta_phi(phiJet_pfakt5[firstfournoisojet.at(0)],phiJetGen_akt5[j])*delta_phi(phiJet_pfakt5[firstfournoisojet.at(0)],phiJetGen_akt5[j]) ) ;
	    if(DR < .1 && (TMath::Abs(ptCorrJet_pfakt5[firstfournoisojet.at(0)]-ptJetGen_akt5[j])/ptJetGen_akt5[j] < 0.5)) assjj = j; 
	  }
	  if(assjj>-1 && ptCorrJet_pfakt5[firstfournoisojet.at(0)]>30. && ptCorrJet_pfakt5[firstfournoisojet.at(1)]>20 ){
	    if(TMath::Abs(etaJet_pfakt5[firstfournoisojet.at(0)]-etaJet_pfakt5[firstfournoisojet.at(1)])>2.5)
	      JECresovbf.Fill((ptCorrJet_pfakt5[firstfournoisojet.at(0)]-ptJetGen_akt5[assjj])/ptJetGen_akt5[assjj]);
	    else
	      JECresovh.Fill((ptCorrJet_pfakt5[firstfournoisojet.at(0)]-ptJetGen_akt5[assjj])/ptJetGen_akt5[assjj]);
	  }
	  deltaetajetisoreco.Fill(etaJet_pfakt5[firstfournoisojet.at(0)]-etaJet_pfakt5[firstfournoisojet.at(1)],weight);
	  double aveeta = (etaJet_pfakt5[firstfournoisojet.at(0)]+etaJet_pfakt5[firstfournoisojet.at(1)])/2;
	  double zeppen1 = etaJet_pfakt5[firstfournoisojet.at(0)] - aveeta;
	  double zeppen2 = etaJet_pfakt5[firstfournoisojet.at(1)] - aveeta;
	  zeppenjetisoreco1.Fill(zeppen1,weight);
	  zeppenjetisoreco2.Fill(zeppen2,weight);	
	}

        nredntp++;

 	if(doPDFweight){
  	  nWeightsPDF1 = nWeightsPDF[0];
 	  nWeightsPDF2 = nWeightsPDF[1];
 	  nWeightsPDF3 = nWeightsPDF[2];
 	  nWeightsPDF4 = nWeightsPDF[3];
          nWeightsPDF5 = nWeightsPDF[4];
          nWeightsPDF6 = nWeightsPDF[5];
          nWeightsPDF7 = nWeightsPDF[6];
          nWeightsPDF8 = nWeightsPDF[7];
          nWeightsPDF9 = nWeightsPDF[8];
          nWeightsPDF10 = nWeightsPDF[9];
  	  for(int iy=0; iy<nWeightsPDF[0] ; iy++)
	    PDFweight1[iy] = pdfWeight[0][iy];
 	  for(int iy=0; iy<nWeightsPDF[1] ; iy++)
 	    PDFweight2[iy] = pdfWeight[1][iy];
 	  for(int iy=0; iy<nWeightsPDF[2] ; iy++)
 	    PDFweight3[iy] = pdfWeight[2][iy];		
          for(int iy=0; iy<nWeightsPDF[3] ; iy++)
            PDFweight4[iy] = pdfWeight[3][iy];
          for(int iy=0; iy<nWeightsPDF[4] ; iy++)
            PDFweight5[iy] = pdfWeight[4][iy];
          for(int iy=0; iy<nWeightsPDF[5] ; iy++)
            PDFweight6[iy] = pdfWeight[5][iy];
          for(int iy=0; iy<nWeightsPDF[6] ; iy++)
            PDFweight7[iy] = pdfWeight[6][iy];
          for(int iy=0; iy<nWeightsPDF[7] ; iy++)
            PDFweight8[iy] = pdfWeight[7][iy];
          for(int iy=0; iy<nWeightsPDF[8] ; iy++)
            PDFweight9[iy] = pdfWeight[8][iy];
          for(int iy=0; iy<nWeightsPDF[9] ; iy++)
            PDFweight10[iy] = pdfWeight[9][iy];
	}

	massgg = higgsisomass;
	ptgg = higgspt;
	massggnewvtx = higgsisomassnewvtx;
	ptggnewvtx = higgsptnewvtx;
        ptphot1 = ptPhot[firstfourisophot.at(0)]; 
        ptphot2 = ptPhot[firstfourisophot.at(1)];   
	timephot1 = timePhot[firstfourisophot.at(0)];
	timephot2 = timePhot[firstfourisophot.at(1)];  
	etaphot1 = etaPhot[firstfourisophot.at(0)];
	etaphot2 = etaPhot[firstfourisophot.at(1)];  
	phiphot1 = phiPhot[firstfourisophot.at(0)];
	phiphot2 = phiPhot[firstfourisophot.at(1)];  
	etascphot1 = etascPhot[firstfourisophot.at(0)];
	etascphot2 = etascPhot[firstfourisophot.at(1)];  
	phiscphot1 = phiscPhot[firstfourisophot.at(0)];
	phiscphot2 = phiscPhot[firstfourisophot.at(1)];  
	E1phot1 = E1Phot[firstfourisophot.at(0)];
	E1phot2 = E1Phot[firstfourisophot.at(1)];  
	E9phot1 = E9Phot[firstfourisophot.at(0)];
	E9phot2 = E9Phot[firstfourisophot.at(1)];  
        r9phot1 = E9Phot[firstfourisophot.at(0)]/escRawPhot[firstfourisophot.at(0)];
        r9phot2 = E9Phot[firstfourisophot.at(1)]/escRawPhot[firstfourisophot.at(1)];
	isemEGphot1 = isophotemeg.at(firstfourisophot.at(0));;
	isemEGphot2 = isophotemeg.at(firstfourisophot.at(1));;
	idloosenewEGphot1 = isophotlooseeg.at(firstfourisophot.at(0));
	idloosenewEGphot2 = isophotlooseeg.at(firstfourisophot.at(1));
	idloose006newEGphot1 = isophotloose006eg.at(firstfourisophot.at(0));
	idloose006newEGphot2 = isophotloose006eg.at(firstfourisophot.at(1));
	idtightnewEGphot1 = isophottighteg.at(firstfourisophot.at(0));
	idtightnewEGphot2 = isophottighteg.at(firstfourisophot.at(1));
	idhggtightnewEGphot1 = isophothggtight.at(firstfourisophot.at(0));
	idhggtightnewEGphot2 = isophothggtight.at(firstfourisophot.at(1));
	idloosenewpuEGphot1 = isophotloosepueg.at(firstfourisophot.at(0));
	idloosenewpuEGphot2 = isophotloosepueg.at(firstfourisophot.at(1));
	idtightnewpuEGphot1 = isophottightpueg.at(firstfourisophot.at(0));
	idtightnewpuEGphot2 = isophottightpueg.at(firstfourisophot.at(1));
	idhggtightnewpuEGphot1 = isophothggtightpu.at(firstfourisophot.at(0));
	idhggtightnewpuEGphot2 = isophothggtightpu.at(firstfourisophot.at(1));
	idcicphot1 = isocic.at(firstfourisophot.at(0));
	idcicphot2 = isocic.at(firstfourisophot.at(1));
	idcicnoelvetophot1 = isocicnoelveto.at(firstfourisophot.at(0));
	idcicnoelvetophot2 = isocicnoelveto.at(firstfourisophot.at(1));
	idlooseEGphot1 = pid_isLoose[firstfourisophot.at(0)];
	idlooseEGphot2 = pid_isLoose[firstfourisophot.at(1)];
	idtightEGphot1 = pid_isTight[firstfourisophot.at(0)];
	idtightEGphot2 = pid_isTight[firstfourisophot.at(1)];
	idloosephot1 = isophotloose.at(firstfourisophot.at(0));
	idloosephot2 = isophotloose.at(firstfourisophot.at(1));
	idmediumphot1 = isophotmedium.at(firstfourisophot.at(0));
	idmediumphot2 = isophotmedium.at(firstfourisophot.at(1));
        idloosecsphot1 = isophotloosecs.at(firstfourisophot.at(0)); 
        idloosecsphot2 = isophotloosecs.at(firstfourisophot.at(1)); 
        idmediumcsphot1 = isophotmediumcs.at(firstfourisophot.at(0)); 
        idmediumcsphot2 = isophotmediumcs.at(firstfourisophot.at(1)); 
	idelephot1 = isophotele.at(firstfourisophot.at(0));
	idelephot2 = isophotele.at(firstfourisophot.at(1));

        pid_haspixelseedphot1 =  hasPixelSeedPhot[firstfourisophot.at(0)]; 
        pid_haspixelseedphot2 =  hasPixelSeedPhot[firstfourisophot.at(1)]; 
        pid_isEMphot1 =  pid_isEM[firstfourisophot.at(0)];
        pid_isEMphot2 =  pid_isEM[firstfourisophot.at(1)];
        pid_jurECALphot1 =  pid_jurECAL[firstfourisophot.at(0)];
        pid_jurECALphot2 =  pid_jurECAL[firstfourisophot.at(1)];
        pid_twrHCALphot1 =  pid_twrHCAL[firstfourisophot.at(0)];
        pid_twrHCALphot2 =  pid_twrHCAL[firstfourisophot.at(1)];
        pid_HoverEphot1 =  pid_HoverE[firstfourisophot.at(0)];
        pid_HoverEphot2 =  pid_HoverE[firstfourisophot.at(1)];
        pid_hlwTrackphot1 =  pid_hlwTrack[firstfourisophot.at(0)];
        pid_hlwTrackphot2 =  pid_hlwTrack[firstfourisophot.at(1)];
        pid_etawidphot1 =  pid_etawid[firstfourisophot.at(0)];
        pid_etawidphot2 =  pid_etawid[firstfourisophot.at(1)];
        pid_hlwTrackNoDzphot1 =  pid_hlwTrackNoDz[firstfourisophot.at(0)];
        pid_hlwTrackNoDzphot2 =  pid_hlwTrackNoDz[firstfourisophot.at(1)];
        pid_hasMatchedConvphot1 =  hasMatchedConvPhot[firstfourisophot.at(0)]; 
        pid_hasMatchedConvphot2 =  hasMatchedConvPhot[firstfourisophot.at(1)]; 
        pid_hasMatchedPromptElephot1 =  hasMatchedPromptElePhot[firstfourisophot.at(0)]; 
        pid_hasMatchedPromptElephot2 =  hasMatchedPromptElePhot[firstfourisophot.at(1)]; 
	
	pid_sminphot1 =  sMinMinPhot[firstfourisophot.at(0)];
	pid_sminphot2 =  sMinMinPhot[firstfourisophot.at(1)];
	pid_smajphot1 =  sMajMajPhot[firstfourisophot.at(0)];
	pid_smajphot2 =  sMajMajPhot[firstfourisophot.at(1)];
	pid_ntrkphot1 =  ntrkiso035Phot[firstfourisophot.at(0)];
	pid_ntrkphot2 =  ntrkiso035Phot[firstfourisophot.at(1)];
	pid_ptisophot1 =  ptiso035Phot[firstfourisophot.at(0)];
	pid_ptisophot2 =  ptiso035Phot[firstfourisophot.at(1)];
	pid_ecalisophot1 =  ecaliso04Phot[firstfourisophot.at(0)];
	pid_ecalisophot2 =  ecaliso04Phot[firstfourisophot.at(1)];
	pid_hcalisophot1 =  hcalovecal04Phot[firstfourisophot.at(0)];
	pid_hcalisophot2 =  hcalovecal04Phot[firstfourisophot.at(1)];
	if(!ieleassocPhot[firstfourisophot.at(0)]){
	  pid_ntrkcsphot1 =  ntrkiso035Phot[firstfourisophot.at(0)]; 
	  pid_ptisocsphot1 =  ptiso035Phot[firstfourisophot.at(0)]; 
	}else{
          pid_ntrkcsphot1 =  ntrkiso035Phot[firstfourisophot.at(0)]-1;  
          pid_ptisocsphot1 =  ptiso035Phot[firstfourisophot.at(0)]-pid_ptElePhot[ieleassocPhot[firstfourisophot.at(0)]];  
	}
        if(!ieleassocPhot[firstfourisophot.at(1)]){ 
          pid_ntrkcsphot2 =  ntrkiso035Phot[firstfourisophot.at(1)];  
          pid_ptisocsphot2 =  ptiso035Phot[firstfourisophot.at(1)];  
        }else{ 
          pid_ntrkcsphot2 =  ntrkiso035Phot[firstfourisophot.at(1)]-1;   
          pid_ptisocsphot2 =  ptiso035Phot[firstfourisophot.at(1)]-pid_ptElePhot[ieleassocPhot[firstfourisophot.at(1)]];   
        } 

	if( firstfournoisojet.at(0) > -1) {
	  ptjet1 = ptJet_pfakt5[firstfournoisojet.at(0)];
	  ptcorrjet1 = ptCorrJet_pfakt5[firstfournoisojet.at(0)];	  
	  etajet1 = etaJet_pfakt5[firstfournoisojet.at(0)];
	  phijet1 = phiJet_pfakt5[firstfournoisojet.at(0)];
	  betajet1 = beta_pfakt5[firstfournoisojet.at(0)][vrankPhotonPairs[0]];
	  betastarjet1 = betaStar_pfakt5[firstfournoisojet.at(0)][vrankPhotonPairs[0]];
	  assjet1 = assoJet(firstfournoisojet.at(0));
	}else{
	  ptjet1 = -999;
	  ptcorrjet1 = -999;
	  etajet1 = -999;	 
	  phijet1 = -999;	 
	  betajet1 = -999.;
	  betastarjet1 = -999.;
	  assjet1 = -999.;
	}
	if( firstfournoisojet.at(1) > -1) {
	  ptjet2 = ptJet_pfakt5[firstfournoisojet.at(1)];
	  ptcorrjet2 = ptCorrJet_pfakt5[firstfournoisojet.at(1)];	  
	  etajet2 = etaJet_pfakt5[firstfournoisojet.at(1)];
	  phijet2 = phiJet_pfakt5[firstfournoisojet.at(1)];
	  betajet2 = beta_pfakt5[firstfournoisojet.at(1)][vrankPhotonPairs[0]];
	  betastarjet2 = betaStar_pfakt5[firstfournoisojet.at(1)][vrankPhotonPairs[0]];
	  assjet2 = assoJet(firstfournoisojet.at(1));
	}else{
	  ptjet2 = -999;
	  ptcorrjet2 = -999;
	  etajet2 = -999;	 
	  phijet2 = -999;	 
	  betajet2 = -999.;
	  betastarjet2 = -999.;
	  assjet2 = -999.;
	}
	if( firstfournoisojet.at(2) > -1) {
	  ptjet3 = ptJet_pfakt5[firstfournoisojet.at(2)];
	  ptcorrjet3 = ptCorrJet_pfakt5[firstfournoisojet.at(2)];	  
	  etajet3 = etaJet_pfakt5[firstfournoisojet.at(2)];
	  phijet3 = phiJet_pfakt5[firstfournoisojet.at(2)];
	}else{
	  ptjet3 = -999;
	  ptcorrjet3 = -999;
	  etajet3 = -999;	 
	  phijet3 = -999;	 
	}
	if( firstfournoisojet.at(3) > -1) {
	  ptjet4 = ptJet_pfakt5[firstfournoisojet.at(3)];
	  ptcorrjet4 = ptCorrJet_pfakt5[firstfournoisojet.at(3)];	  
	  etajet4 = etaJet_pfakt5[firstfournoisojet.at(3)];
	  phijet4 = phiJet_pfakt5[firstfournoisojet.at(3)];
	}else{
	  ptjet4 = -999;
	  ptcorrjet4 = -999;
	  etajet4 = -999;	 
	  phijet4 = -999;	 
	}
	if( firstfournoisojet.at(0) > -1 && firstfournoisojet.at(1) > -1) {
	  deltaeta = etaJet_pfakt5[firstfournoisojet.at(0)]-etaJet_pfakt5[firstfournoisojet.at(1)];
	  double aveeta = (etaJet_pfakt5[firstfournoisojet.at(0)]+etaJet_pfakt5[firstfournoisojet.at(1)])/2;
	  zeppenjet = etahiggsiso - aveeta;
	  invmassjet = twojetsmassiso;
	  deltaphi = delta_phi(phihiggsiso,phitwojetsiso);
	  deltaphinewvtx = delta_phi(phihiggsisonewvtx,phitwojetsiso);
	}else{
	  deltaeta = -999.;
	  zeppenjet = -999.;
	  invmassjet = -999.;
	  deltaphi = -999.;
	  deltaphinewvtx = -999.;
	}	  

    ///////////////////////////////////////////////////
	// met = epfMet;
	// phimet = phipfMet;

    sMet_ = sMet;
    eMet_ = eMet;
    phiMet_ = phiMet;
    signifMet_ = signifMet;
    sCorrMet_ = sCorrMet;
    eCorrMet_ = eCorrMet;
    phiCorrMet_ = phiCorrMet;
    signifCorrMet_ = signifCorrMet;
    smuCorrMet_ = smuCorrMet;
    emuCorrMet_ = emuCorrMet;
    phimuCorrMet_ = phimuCorrMet;
    signifmuCorrMet_ = signifmuCorrMet;
    sNoHFMet_ = sNoHFMet;
    eNoHFMet_ = eNoHFMet;
    phiNoHFMet_ = phiNoHFMet;
    signifNoHFMet_ = signifNoHFMet;
    stcMet_ = stcMet;
    etcMet_ = etcMet;
    phitcMet_ = phitcMet;
    signiftcMet_ = signiftcMet;
    sglobalPfMet_ = sglobalPfMet;
    eglobalPfMet_ = eglobalPfMet;
    phiglobalPfMet_ = phiglobalPfMet;
    signifglobalPfMet_ = signifglobalPfMet;
    scentralPfMet_ = scentralPfMet;
    ecentralPfMet_ = ecentralPfMet;
    phicentralPfMet_ = phicentralPfMet;
    signifcentralPfMet_ = signifcentralPfMet;
    ///vector
    eassocPfMet_ = eassocPfMet[vrankPhotonPairs[0]];
    phiassocPfMet_ = phiassocPfMet[vrankPhotonPairs[0]];
    signifassocPfMet_ = signifassocPfMet[vrankPhotonPairs[0]];
    eassocOtherVtxPfMet_ = eassocOtherVtxPfMet[vrankPhotonPairs[0]];
    phiassocOtherVtxPfMet_ = phiassocOtherVtxPfMet[vrankPhotonPairs[0]];
    signifassocOtherVtxPfMet_ = signifassocOtherVtxPfMet[vrankPhotonPairs[0]];
    etrkPfMet_ = etrkPfMet[vrankPhotonPairs[0]];
    phitrkPfMet_ = phitrkPfMet[vrankPhotonPairs[0]];
    signiftrkPfMet_ = signiftrkPfMet[vrankPhotonPairs[0]];
    ecleanPfMet_ = ecleanPfMet[vrankPhotonPairs[0]];
    phicleanPfMet_ = phicleanPfMet[vrankPhotonPairs[0]];
    signifcleanPfMet_ = signifcleanPfMet[vrankPhotonPairs[0]];
    ecleanedSaclayPfMet_ = ecleanedSaclayPfMet[vrankPhotonPairs[0]];
    phicleanedSaclayPfMet_ = phicleanedSaclayPfMet[vrankPhotonPairs[0]];
    signifcleanedSaclayPfMet_ = signifcleanedSaclayPfMet[vrankPhotonPairs[0]];
    eminTypeICleanSaclayPfMet_ = eminTypeICleanSaclayPfMet[vrankPhotonPairs[0]];
    phiminTypeICleanSaclayPfMet_ = phiminTypeICleanSaclayPfMet[vrankPhotonPairs[0]];
    signifminTypeICleanSaclayPfMet_ = signifminTypeICleanSaclayPfMet[vrankPhotonPairs[0]];
    globalPfSums_ = globalPfSums[vrankPhotonPairs[0]];
    /// end vector
    spfMet_ = spfMet;
    epfMet_ = epfMet;
    phipfMet_ = phipfMet;
    signifpfMet_ = signifpfMet;
    spfMetType1_ = spfMetType1;
    epfMetType1_ = epfMetType1;
    phipfMetType1_ = phipfMetType1;
    signifpfMetType1_ = signifpfMetType1;
    sMetGen_ = sMetGen;
    eMetGen_ = eMetGen;
    phiMetGen_ = phiMetGen;
    signifMetGen_ = signifMetGen;
    sMetGen2_ = sMetGen2;
    eMetGen2_ = eMetGen2;
    phiMetGen2_ = phiMetGen2;

/////// old
    /*
    sMet_  ;
    eMet_  ;
    phiMet_;
    signifMet_;
    sCorrMet_  ;
    eCorrMet_  ;
    phiCorrMet_;
    signifCorrMet_;
    smuCorrMet_  ;
    emuCorrMet_  ;
    phimuCorrMet_;
    signifmuCorrMet_;
    sNoHFMet_  ;
    eNoHFMet_  ;
    phiNoHFMet_;
    signifNoHFMet_;
    stcMet_  ;
    etcMet_  ;
    phitcMet_;
    signiftcMet_;
    sglobalPfMet_;
    eglobalPfMet_;
    phiglobalPfMet_;
    signifglobalPfMet_;
    scentralPfMet_;
    ecentralPfMet_;
    phicentralPfMet_;
    signifcentralPfMet_;
    eassocPfMet_;   //[nvertex]
    phiassocPfMet_;   //[nvertex]
    signifassocPfMet_;   //[nvertex]
    eassocOtherVtxPfMet_;   //[nvertex]
    phiassocOtherVtxPfMet_;   //[nvertex]
    signifassocOtherVtxPfMet_;   //[nvertex]
    etrkPfMet_;   //[nvertex]
    phitrkPfMet_;   //[nvertex]
    signiftrkPfMet_;   //[nvertex]
    ecleanPfMet_;   //[nvertex]
    phicleanPfMet_;   //[nvertex]
    signifcleanPfMet_;   //[nvertex]
    ecleanedSaclayPfMet_;   //[nvertex]
    phicleanedSaclayPfMet_;   //[nvertex]
    signifcleanedSaclayPfMet_;   //[nvertex]
    eminTypeICleanSaclayPfMet_;   //[nvertex]
    phiminTypeICleanSaclayPfMet_;   //[nvertex]
    signifminTypeICleanSaclayPfMet_;   //[nvertex]
    globalPfSums_;
    spfMet_  ;
    epfMet_  ;
    phipfMet_;
    signifpfMet_;
    spfMetType1_;
    epfMetType1_;
    phipfMetType1_;
    signifpfMetType1_;
    sMetGen_  ;
    eMetGen_  ;
    phiMetGen_;
    signifMetGen_;
    sMetGen2_  ;
    eMetGen2_  ;
    phiMetGen2_;
    */
    //////////////////////////////////////////////////

	nvtx = nvertex;

        runRN = run;
        eventRN = event;
        lumi = lbn;
        rhoPFRN = rhoPF;
	
#ifdef DEBUG
    cout << "[DEBUG] before gamma-jet combination" << endl;
#endif

	TLorentzVector twog1j = thehiggs + thejet1;
	TLorentzVector twog2j = thehiggs + thejet1 + thejet2; 
#ifdef DEBUG
    cout << "[DEBUG] after gamma-jet combination" << endl;
#endif
	invmass2g1j = twog1j.M();
	invmass2g2j = twog2j.M();
	

	if(ptPhot[firstfourisophot.at(0)] > ptphot1cut && ptPhot[firstfourisophot.at(1)] > ptphot2cut){
	  higgsmassjustisocutreco.Fill(higgsisomass,weight);
	  higgsmassjustisocutrecofull.Fill(higgsisomass,weight);
	}

	if( ptCorrJet_pfakt5[firstfournoisojet.at(0)] > ptjet1cut && ptCorrJet_pfakt5[firstfournoisojet.at(1)] > ptjet2cut 
	    && ptPhot[firstfourisophot.at(0)] > ptphot1cut && ptPhot[firstfourisophot.at(1)] > ptphot2cut){
	  higgsmassisojetptcutreco.Fill(higgsisomass,weight);
	  higgsmassisojetptcutrecofull.Fill(higgsisomass,weight);
	  deltaetajetreco.Fill(etaJet_pfakt5[firstfournoisojet.at(0)]-etaJet_pfakt5[firstfournoisojet.at(1)],weight);
	  //	  double zeppen = higgsreco_pt - aveeta;
	  if(TMath::Abs(etaJet_pfakt5[firstfournoisojet.at(0)]-etaJet_pfakt5[firstfournoisojet.at(1)])>deltaetacut){
	    higgsmassisocutreco.Fill(higgsisomass,weight);
	    higgsmassisocutrecofull.Fill(higgsisomass,weight);
	    double aveeta = (etaJet_pfakt5[firstfournoisojet.at(0)]+etaJet_pfakt5[firstfournoisojet.at(1)])/2;
	    double zeppen = etahiggsiso - aveeta;
	    zeppenhiggsisoreco.Fill(zeppen,weight);
	    if(TMath::Abs(zeppen)<zeppencut) {
	      higgsmassisocutzeppreco.Fill(higgsisomass,weight);
	      higgsmassisocutzepprecofull.Fill(higgsisomass,weight);
	      dijetmassisoreco.Fill(twojetsmassiso,weight);
	      if(twojetsmassiso>dijetmasscut){
		higgsmassisocutzeppdijetreco.Fill(higgsisomass,weight);
		higgsmassisocutzeppdijetrecofull.Fill(higgsisomass,weight);
	      }
	    }
	  }
	}

      } 
      

#ifdef DEBUG1
        if(recoPreselection && !genPreselection)
        {
            TLorentzVector genPhot1, genPhot2;	
            genPhot1.SetPtEtaPhiE( ptMC[index_phot1], etaMC[index_phot1], phiMC[index_phot1], eMC[index_phot1]);
            genPhot2.SetPtEtaPhiE( ptMC[index_phot2], etaMC[index_phot2], phiMC[index_phot2], eMC[index_phot2]);
            TLorentzVector diphot = genPhot1 + genPhot2;

            cout << endl;
            if( gen_custom_processId > 9000)
                cout << "[DEBUG] bkg process" << endl;
            else
                cout << "[DEBUG] sig process" << endl;

            cout << "[DEBUG] firstfourgenphot.at(0) = " << firstfourgenphot.at(0) << endl;
            cout << "[DEBUG] firstfourgenphot.at(1) = " << firstfourgenphot.at(1) << endl;
            cout << "[DEBUG] index_phot1 = " << index_phot1 << endl;
            cout << "[DEBUG] index_phot2 = " << index_phot2 << endl;
            cout << "[DEBUG] genPhot1:: " << ptMC[index_phot1] << ", " << etaMC[index_phot1] << ", " << phiMC[index_phot1] << ", " << eMC[index_phot1] << endl;
            cout << "[DEBUG] genPhot2:: " << ptMC[index_phot2] << ", " << etaMC[index_phot2] << ", " << phiMC[index_phot2] << ", " << eMC[index_phot2] << endl;
            cout << "[DEBUG] gen_mass_diphoton:: " << diphot.M() << endl;
            cout << "[DEBUG] recoPhot1:: " <<  ptPhot[firstfourisophot.at(0)] << ", " << etaPhot[firstfourisophot.at(0)] << ", " << phiPhot[firstfourisophot.at(0)] << endl; 
            cout << "[DEBUG] recoPhot2:: " <<  ptPhot[firstfourisophot.at(1)] << ", " << etaPhot[firstfourisophot.at(1)] << ", " << phiPhot[firstfourisophot.at(1)] << endl; 
            cout << "[DEBUG] reco  higgs mass:: " << higgsisomass << endl;
        }
#endif
    

      if(recoPreselection)
	    ana_tree->Fill();

	

   } /// loop over events

   timer.Stop();
   cout << "Original number of events: " << NtotEvents << endl;
   cout << "Processed events:          " << nprocessed << endl; 
   cout << "Processed events/s (CPU Time):          " << ((float)nprocessed)/timer.CpuTime() << endl; 
   cout << "Processed events/s (Real Time):          " << ((float)nprocessed)/timer.RealTime() << endl; 
   cout << "Events in reduced ntuple:  " << nredntp << endl; 
   
   hOutputFile->Write() ;
   if (myjson)
     delete myjson;
}

void RedNtpTree::SetPuWeights(std::string puWeightFile)
{
  if (puWeightFile == "")
    {
      std::cout << "you need a weights file to use this function" << std::endl;
       return;
    }
  
  std::cout << "PU REWEIGHTING:: Using file " << puWeightFile << std::endl;
  
  TFile *f_pu  = new TFile(puWeightFile.c_str(),"READ");

  f_pu->cd();

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
    weight=weights->GetBinContent(i+1);
    sumPuWeights+=weight;
    puweights_.push_back(weight);
  }
  
  //std::cout << "weights sum is " << sumPuWeights << std::endl;
}

// std::vector<std::string> tokenize_str(const std::string & str,
//                                       const std::string & delims=", \t")
// {
//   using namespace std;
//   // Skip delims at beginning, find start of first token
//   string::size_type lastPos = str.find_first_not_of(delims.c_str(), 0, delims.length());
//   // Find next delimiter @ end of token
//   string::size_type pos     = str.find_first_of(delims.c_str(), lastPos, delims.length());
 
//   // output vector
//   vector<string> tokens;
 
//   while (string::npos != pos || string::npos != lastPos)
//     {
//       // Found a token, add it to the vector.
//       tokens.push_back(str.substr(lastPos, pos - lastPos));
//       // Skip delims.  Note the "not_of". this is beginning of token
//       lastPos = str.find_first_not_of(delims.c_str(), pos , delims.length());
//       // Find next delimiter at end of token.
//       pos     = str.find_first_of(delims.c_str(), lastPos, delims.length());
//     }
 
//   return tokens;
// }
std::vector<std::string> tokenize_str(const std::string & str,
				      const std::string & delims=", \t")
{
  using namespace std;
  // Skip delims at beginning, find start of first token
  string::size_type lastPos = 0;
  
  // Find next delimiter @ end of token
  string::size_type pos = str.find(delims, 0);
  if (pos == string::npos)
    pos = str.length();
  
   // output vector
  vector<string> tokens;
  
  while (string::npos != lastPos)
    {
      // Found a token, add it to the vector.
      tokens.push_back(str.substr(lastPos, pos - lastPos));
      // Skip delims.  Note the "not_of". this is beginning of token
      lastPos = str.find(delims, pos+delims.length());
      if (lastPos == string::npos && pos!= str.length())
	lastPos=pos+delims.length();
      // Find next delimiter at end of token.
      pos     = str.find(delims, lastPos+1);
      if (pos == string::npos)
 	pos = str.length();
    }

   return tokens;
 }


void RedNtpTree::SetPtWeights(std::string ptWeightFile)
{
  if (ptWeightFile == "")
    {
      std::cout << "you need a weights file to use this function" << std::endl;
       return;
    }
  
  std::cout << "PT REWEIGHTING:: Using file " << ptWeightFile << std::endl;

  std::vector<std::string> tokens=tokenize_str(ptWeightFile,"Kfactors_");
  //std::vector<std::string> tokens=tokenize_str(ptWeightFile,"weight_ptH_");
  TString massValue;

//    for (int i=0;i<tokens.size();++i)
//     std::cout << tokens[i] << std::endl;

  if (tokens.size()>1)
    {
      std::vector<std::string> newTokens= tokenize_str(tokens[1],"_");
      if (newTokens.size()>0)
	massValue=TString(newTokens[0]);
    }
  
  std::cout << "PT REWEIGHTING:: mass values used " << massValue << std::endl;

  TFile *f_pt  = new TFile(ptWeightFile.c_str(),"READ");
  f_pt->cd();
  
  ptweights_ =(TH1D*)  f_pt->Get("kfactors/kfact_mh"+massValue+"_ren"+massValue+"_fac"+massValue);
  //f_pt->Get("powheg_weight/weight_hqt_fehipro_fit_120")->Draw();
  //ptweights_ =(TH1D*)  f_pt->Get("powheg_weight/weight_hqt_fehipro_fit_120");

  if (!ptweights_ )
    {
      std::cout << "weights histograms  not found in file " << ptWeightFile << std::endl;
      return;
    }
}

bool RedNtpTree::assoJet(int i){

  bool ass(0);
  for(int j=0; j<nJetGen_akt5; j++){	
    double DR = sqrt(delta_eta(etaJet_pfakt5[i],etaJetGen_akt5[j])*delta_eta(etaJet_pfakt5[i],etaJetGen_akt5[j]) + 
		     delta_phi(phiJet_pfakt5[i],phiJetGen_akt5[j])*delta_phi(phiJet_pfakt5[i],phiJetGen_akt5[j]) ) ;
    if(DR < .1 && TMath::Abs(ptCorrJet_pfakt5[i]-ptJetGen_akt5[j])/ptJetGen_akt5[j]  < 0.5) ass = 1; 
  }

  return ass;

}

void RedNtpTree::correctJets(int shift, float smear)
{

  for(int i=0; i<nJet_pfakt5; i++){
    
    double increase_endcap = 1;
    if(TMath::Abs(etaJet_pfakt5[i])>2.5 && TMath::Abs(etaJet_pfakt5[i])<3.4) increase_endcap = 2;  
    ptCorrJet_pfakt5[i] *= 1 + increase_endcap * shift * jetsyst_->getJESUncertainty(etaJet_pfakt5[i],ptCorrJet_pfakt5[i]);
    ptJet_pfakt5[i] *= 1 + increase_endcap * shift * jetsyst_->getJESUncertainty(etaJet_pfakt5[i],ptCorrJet_pfakt5[i]);
    eJet_pfakt5[i] *= 1 + increase_endcap * shift * jetsyst_->getJESUncertainty(etaJet_pfakt5[i],eJet_pfakt5[i]);
    
    if(smear){
      
      int ass(-999);
      for(int j=0; j<nJetGen_akt5; j++){	
	double DR = sqrt(delta_eta(etaJet_pfakt5[i],etaJetGen_akt5[j])*delta_eta(etaJet_pfakt5[i],etaJetGen_akt5[j]) + 
			 delta_phi(phiJet_pfakt5[i],phiJetGen_akt5[j])*delta_phi(phiJet_pfakt5[i],phiJetGen_akt5[j]) ) ;
	if(DR < .1 && TMath::Abs(ptCorrJet_pfakt5[i]-ptJetGen_akt5[j])/ptJetGen_akt5[j]  < 0.5) ass = j; 
      }
      
      if(ass>-1){
	double scaling = (ptJetGen_akt5[ass] + (1 + smear) * (ptCorrJet_pfakt5[i] - ptJetGen_akt5[ass]))/ptCorrJet_pfakt5[i];
	ptCorrJet_pfakt5[i] *= scaling;
	ptJet_pfakt5[i] *= scaling;
	eJet_pfakt5[i] *= scaling;
      }      

    }

  }

}

void RedNtpTree::correctPhotons(bool energyRegression)
{
  for (int iPho=0;iPho<nPhot;++iPho)
    {
      bool isEBPho=(fabs(etascPhot[iPho])<1.479);
      float R9Pho=E9Phot[iPho]/escRawPhot[iPho];
      float scaleCorrection=scaleCorrections_->getScaleOffset(run,isEBPho,R9Pho,fabs(etascPhot[iPho]));
      float smearing=scaleCorrections_->getSmearing(run,isEBPho,R9Pho,fabs(etascPhot[iPho]));
      //      std::cout << scaleCorrection << "," << smearing << " run " << run << " isEB " << isEBPho << " R9 " << R9Pho << std::endl;

      //In  MC apply smearing as energy correction
      if (nMC>0)
	scaleCorrection=gen_->Gaus(1.,smearing);

      //energies correction
      if (!energyRegression)
	{
	  ptPhot[iPho]=ptPhot[iPho]*scaleCorrection;   //[nPhot]
	  ePhot[iPho]=ePhot[iPho]*scaleCorrection;   //[nPhot]
	}
      else
	{
	  ptPhot[iPho]=escRegrPhot[iPho]/TMath::CosH(etaPhot[iPho])*scaleCorrection;   //[nPhot]
	  ePhot[iPho]=escRegrPhot[iPho]*scaleCorrection;   //[nPhot]
	}
      escPhot[iPho]=escPhot[iPho]*scaleCorrection;   //[nPhot]
      escRawPhot[iPho]=escRawPhot[iPho]*scaleCorrection;   //[nPhot]
      eseedPhot[iPho]=eseedPhot[iPho]*scaleCorrection;   //[nPhot]
      E1Phot[iPho]=E1Phot[iPho]*scaleCorrection;   //[nPhot]
      E9Phot[iPho]=E9Phot[iPho]*scaleCorrection;   //[nPhot]
      E25Phot[iPho]=E25Phot[iPho]*scaleCorrection;   //[nPhot]
    }
}

void RedNtpTree::DoPDFWeighting()
{
  
  doPDFweight = 1;
  std::cout << "writing weights for PDF systematics out " << std::endl;
  
}


void RedNtpTree::SetAllGenVarToMinus999()
{
    gen_pt_gamma1 = -999;
    gen_pt_gamma2 = -999;
    gen_eta_gamma1 = -999;
    gen_eta_gamma2 = -999;
    gen_phi_gamma1 = -999;
    gen_phi_gamma2 = -999;
    
    gen_pt_genjet1 = -999;
    gen_pt_genjet2 = -999;
    gen_eta_genjet1 = -999;
    gen_eta_genjet2 = -999;
    gen_phi_genjet1 = -999;
    gen_phi_genjet2 = -999;
    
    // gen_pt_VectorBoson = -999;
    // gen_phi_VectorBoson = -999;
    // gen_eta_VectorBoson = -999;
    
    gen_mass_diphoton = -999;
    gen_pt_diphoton = -999;
    gen_eta_diphoton = -999;
    gen_phi_diphoton = -999;
    
    gen_mass_dijet = -999;
    gen_pt_dijet = -999;
    gen_eta_dijet = -999;
    gen_phi_dijet = -999;
    
    gen_zeppenfeld = -999;
}


void RedNtpTree::SetAllRecoVarToMinus999()
{
massgg = -999;
ptgg = -999;
massggnewvtx = -999;
ptggnewvtx = -999;
ptphot1 = -999;
ptphot2 = -999;
etaphot1 = -999;
etaphot2 = -999;
phiphot1 = -999;
phiphot2 = -999;
timephot1 = -999; 
timephot2 = -999; 
E1phot1 = -999;
E1phot2 = -999;
E9phot1 = -999;
E9phot2 = -999;
ptjet1 = -999;
ptjet2 = -999;
ptjet3 = -999;
ptjet4 = -999;
ptcorrjet1 = -999;
ptcorrjet2 = -999;
ptcorrjet3 = -999;
ptcorrjet4 = -999;
etajet1 = -999;
etajet2 = -999;
etajet3 = -999;
etajet4 = -999;
phijet1 = -999;
phijet2 = -999;
phijet3 = -999;
phijet4 = -999;
betajet1 = -999;
betajet2 = -999;
betastarjet1 = -999;
betastarjet2 = -999;
assjet1 = -999;
assjet2 = -999;
deltaeta = -999;
zeppenjet = -999;
deltaphi = -999;
deltaphinewvtx = -999;
invmassjet = -999;
invmass2g1j = -999;
invmass2g2j = -999;
nvtx = -999;

   //////////////////////////////////////
    sMet_   = -999;
    eMet_   = -999;
    phiMet_ = -999;
    signifMet_ = -999;
    sCorrMet_   = -999;
    eCorrMet_   = -999;
    phiCorrMet_ = -999;
    signifCorrMet_ = -999;
    smuCorrMet_   = -999;
    emuCorrMet_   = -999;
    phimuCorrMet_ = -999;
    signifmuCorrMet_ = -999;
    sNoHFMet_   = -999;
    eNoHFMet_   = -999;
    phiNoHFMet_ = -999;
    signifNoHFMet_ = -999;
    stcMet_   = -999;
    etcMet_   = -999;
    phitcMet_ = -999;
    signiftcMet_ = -999;
    sglobalPfMet_ = -999;
    eglobalPfMet_ = -999;
    phiglobalPfMet_ = -999;
    signifglobalPfMet_ = -999;
    scentralPfMet_ = -999;
    ecentralPfMet_ = -999;
    phicentralPfMet_ = -999;
    signifcentralPfMet_ = -999;
    eassocPfMet_ = -999;   //[nvertex]
    phiassocPfMet_ = -999;   //[nvertex]
    signifassocPfMet_ = -999;   //[nvertex]
    eassocOtherVtxPfMet_ = -999;   //[nvertex]
    phiassocOtherVtxPfMet_ = -999;   //[nvertex]
    signifassocOtherVtxPfMet_ = -999;   //[nvertex]
    etrkPfMet_ = -999;   //[nvertex]
    phitrkPfMet_ = -999;   //[nvertex]
    signiftrkPfMet_ = -999;   //[nvertex]
    ecleanPfMet_ = -999;   //[nvertex]
    phicleanPfMet_ = -999;   //[nvertex]
    signifcleanPfMet_ = -999;   //[nvertex]
    ecleanedSaclayPfMet_ = -999;   //[nvertex]
    phicleanedSaclayPfMet_ = -999;   //[nvertex]
    signifcleanedSaclayPfMet_ = -999;   //[nvertex]
    eminTypeICleanSaclayPfMet_ = -999;   //[nvertex]
    phiminTypeICleanSaclayPfMet_ = -999;   //[nvertex]
    signifminTypeICleanSaclayPfMet_ = -999;   //[nvertex]
    globalPfSums_ = -999;
    spfMet_   = -999;
    epfMet_   = -999;
    phipfMet_ = -999;
    signifpfMet_ = -999;
    spfMetType1_ = -999;
    epfMetType1_ = -999;
    phipfMetType1_ = -999;
    signifpfMetType1_ = -999;
    sMetGen_   = -999;
    eMetGen_   = -999;
    phiMetGen_ = -999;
    signifMetGen_ = -999;
    sMetGen2_   = -999;
    eMetGen2_   = -999;
    phiMetGen2_ = -999;
   //////////////////////////////////////

npu = -999;
isemEGphot1 = -999;
isemEGphot2 = -999;
idloosenewEGphot1 = -999;
idloosenewEGphot2 = -999;
idloose006newEGphot1 = -999;
idloose006newEGphot2 = -999;
idtightnewEGphot1 = -999;
idtightnewEGphot2 = -999;
idhggtightnewEGphot1 = -999;
idhggtightnewEGphot2 = -999;
idloosenewpuEGphot1 = -999;
idloosenewpuEGphot2 = -999;
idtightnewpuEGphot1 = -999;
idtightnewpuEGphot2 = -999;
idhggtightnewpuEGphot1 = -999;
idhggtightnewpuEGphot2 = -999;
idcicphot1 = -999;
idcicphot2 = -999;
idcicnoelvetophot1 = -999;
idcicnoelvetophot2 = -999;
idlooseEGphot1 = -999;
idlooseEGphot2 = -999;
idtightEGphot1 = -999;
idtightEGphot2 = -999;
idloosephot1 = -999; 
idloosephot2 = -999; 
idmediumphot1 = -999; 
idmediumphot2 = -999; 
idloosecsphot1 = -999;
idloosecsphot2 = -999;
idmediumcsphot1 = -999;
idmediumcsphot2 = -999;
idelephot1 = -999;
idelephot2 = -999;
    pid_haspixelseedphot1 = -999; 
    pid_haspixelseedphot2 = -999; 
    pid_isEMphot1 = -999;
    pid_isEMphot2 = -999;
       pid_jurECALphot1 = -999;
       pid_jurECALphot2 = -999;
       pid_twrHCALphot1 = -999;
       pid_twrHCALphot2 = -999;
       pid_HoverEphot1 = -999;
       pid_HoverEphot2 = -999;
       pid_hlwTrackphot1 = -999;
       pid_hlwTrackphot2 = -999;
       pid_etawidphot1 = -999;
       pid_etawidphot2 = -999;
       pid_sminphot1 = -999;
       pid_sminphot2 = -999;
       pid_smajphot1 = -999;
       pid_smajphot2 = -999;
       pid_ntrkphot1 = -999;
       pid_ntrkphot2 = -999;
       pid_ptisophot1 = -999;
       pid_ptisophot2 = -999;
       pid_ntrkcsphot1 = -999; 
       pid_ntrkcsphot2 = -999; 
       pid_ptisocsphot1 = -999; 
       pid_ptisocsphot2 = -999; 
       pid_ecalisophot1 = -999;
       pid_ecalisophot2 = -999;
       pid_hcalisophot1 = -999;
       pid_hcalisophot2 = -999;
    runRN = -999;
    eventRN = -999;
    lumi = -999;
      rhoPFRN = -999;
      pid_hlwTrackNoDzphot1 = -999;
      pid_hlwTrackNoDzphot2 = -999;
      pid_hasMatchedConvphot1 = -999;
      pid_hasMatchedConvphot2 = -999;
      pid_hasMatchedPromptElephot1 = -999;
      pid_hasMatchedPromptElephot2 = -999;
      r9phot1 = -999;
      r9phot2 = -999;
      etascphot1 = -999;
      etascphot2 = -999;
      phiscphot1 = -999;
      phiscphot2 = -999;
      pu_weight = -999;
      pt_weight = -999;

   nWeightsPDF1 = -999;
   nWeightsPDF2 = -999;
   nWeightsPDF3 = -999;
   nWeightsPDF4 = -999;
   nWeightsPDF5 = -999;
   nWeightsPDF6 = -999;
   nWeightsPDF7 = -999;
   nWeightsPDF8 = -999;
   nWeightsPDF9 = -999;
   nWeightsPDF10 = -999;
   
   weight = -999;
}



/************************************
 *                                  *
 *                                  *
 *         Reco  Selection          *
 *                                  *
 *                                  *
 ************************************/



bool RedNtpTree::cutID(int i, photonidcuts const& pid, vector<bool> *vpass) {

  // Use photon supercluster energy (would be e5x5 if r9>0.93 otherwise)
  bool ntrkiso = ntrkiso035Phot[i] < pid.tracknb;
  bool ptiso = (ptiso035Phot[i] / ptPhot[i] < pid.trackiso_rel);
//   bool ecaliso = (ecaliso04Phot[i] / ePhot[i] < pid.ecaliso_rel ||
//                    ecaliso04Phot[i] < pid.ecaliso_abs);
  // in order to fix bug in endcap use equivalent egamma variables 
  bool ecaliso = (pid_jurECAL[i]*cosh(etaPhot[i]) / ePhot[i] < pid.ecaliso_rel/2. ||
                  pid_jurECAL[i]*cosh(etaPhot[i]) < pid.ecaliso_abs/2.);
//    double fhcal = hcalovecal04Phot[i];
//    bool hcaliso = (fhcal < pid.hcaliso_rel ||
// 		  fhcal*ptPhot[i] < pid.hcaliso_abs);
   // in order to fix bug in endcap use equivalent egamma variables
  double fhcal = pid_HoverE[i] + pid_twrHCAL[i] / ptPhot[i];
  bool hcaliso = (fhcal < pid.hcaliso_rel/2. ||
		  fhcal*ptPhot[i] < pid.hcaliso_abs/2.);
  bool smaj = sMajMajPhot[i] < pid.smajmaj;
  bool smin = sMinMinPhot[i] < pid.sminmin;
  bool smin_min = sMinMinPhot[i] > pid.sminmin_min;
  //bool eta = TMath::Abs(etaPhot[i]) < 2.5; 
  bool eta = true;


//   if(TMath::Abs(etaPhot[i]) > 1.4442) {
//     smaj = 1; smin = 1; smin_min = 1;
//   }
  
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


bool RedNtpTree::cutIDcs(int i, photonidcuts const& pid, vector<bool> *vpass) { 
 
  // Use photon supercluster energy (would be e5x5 if r9>0.93 otherwise) 
  bool ntrkiso = ntrkiso035Phot[i] < pid.tracknb;  
  bool ptiso = (ptiso035Phot[i])/ ptPhot[i] < pid.trackiso_rel; 

  if(ieleassocPhot[i] > -1){
    ntrkiso = ntrkiso035Phot[i]-1 < pid.tracknb;   
    ptiso = (ptiso035Phot[i]-pid_ptElePhot[ieleassocPhot[i]])/ ptPhot[i] < pid.trackiso_rel;  
  }
  bool ecaliso = (ecaliso04Phot[i] / ePhot[i] < pid.ecaliso_rel || 
                  ecaliso04Phot[i] < pid.ecaliso_abs); 
  double fhcal = hcalovecal04Phot[i]; 
  // in order to fix bug in endcap use equivalent egamma variables 
  //double fhcal = pid_HoverE[i] + pid_twrHCAL[i] / ptPhot[i]; 
  bool hcaliso = (fhcal < pid.hcaliso_rel || 
                  fhcal*ptPhot[i] < pid.hcaliso_abs); 
  bool smaj = sMajMajPhot[i] < pid.smajmaj; 
  bool smin = sMinMinPhot[i] < pid.sminmin; 
  bool smin_min = sMinMinPhot[i] > pid.sminmin_min; 
  //bool eta = TMath::Abs(etaPhot[i]) < 2.5;  
  bool eta = true; 
 
 
  //   if(TMath::Abs(etaPhot[i]) > 1.4442) { 
  //     smaj = 1; smin = 1; smin_min = 1; 
  //   } 
   
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
 


bool RedNtpTree::cutIDpresel(int i, photonidcuts const& pid, vector<bool> *vpass) {

  // Use photon supercluster energy (would be e5x5 if r9>0.93 otherwise)
  bool ntrkiso = ntrkiso035Phot[i] < pid.tracknb;
  bool ptiso = (ptiso035Phot[i] / ptPhot[i] < pid.trackiso_rel);
  bool ecaliso =  (ecaliso04Phot[i] / ePhot[i] < pid.ecaliso_rel ||
		   ecaliso04Phot[i] < pid.ecaliso_abs) ||
                  (pid_jurECAL[i] < ptPhot[i] * pid.ecaliso_rel + pid.ecaliso_abs);
  double fhcal = hcalovecal04Phot[i];
  // in order to fix bug in endcap use equivalent egamma variables
  //double fhcal = pid_HoverE[i] + pid_twrHCAL[i] / ptPhot[i];
  bool hcaliso = (fhcal < pid.hcaliso_rel ||
                  fhcal*ptPhot[i] < pid.hcaliso_abs) ||
                 (pid_twrHCAL[i] < ptPhot[i] * pid.hcaliso_rel + pid.hcaliso_abs);
  bool smaj = sMajMajPhot[i] < pid.smajmaj;
  bool smin = sMinMinPhot[i] < pid.sminmin;
  bool smin_min = sMinMinPhot[i] > pid.sminmin_min;
  //bool eta = TMath::Abs(etaPhot[i]) < 2.5; 
  bool eta = true;


//   if(TMath::Abs(etaPhot[i]) > 1.4442) {
//     smaj = 1; smin = 1; smin_min = 1;
//   }
  
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



bool RedNtpTree::cutIDele(int i, photonidelecuts const& pid, vector<bool> *vpass) {

  if(ieleassocPhot[i] < 0) return 0;

  // Use photon supercluster energy (would be e5x5 if r9>0.93 otherwise)
  bool ptiso,ecaliso, hcaliso, hoveiso, setaeta, deta, dphi, minhits, dcot, dist;
  if(TMath::Abs(etascPhot[i]) < 1.4442) {
    ptiso = pid_hlwTrackElePhot[ieleassocPhot[i]] < ptPhot[i] * pid.trackiso_relEB;
    ecaliso = pid_jurECALElePhot[ieleassocPhot[i]] < ptPhot[i] * pid.ecaliso_relEB;
    hcaliso = pid_twrHCALElePhot[ieleassocPhot[i]] < ptPhot[i] * pid.hcaliso_relEB;
    hoveiso = pid_HoverEElePhot[ieleassocPhot[i]] < pid.hovereisoEB;
    setaeta = pid_etawidElePhot[ieleassocPhot[i]] < pid.setaetaEB;
    deta = pid_detavtxElePhot[ieleassocPhot[i]] < pid.detaEB;
    dphi = pid_dphivtxElePhot[ieleassocPhot[i]] < pid.dphiEB;
    minhits = pid_mishitsElePhot[ieleassocPhot[i]] < pid.minhitsEB;
    dcot = TMath::Abs(pid_dcotElePhot[ieleassocPhot[i]]) > pid.dcotEB;      
    dist = TMath::Abs(pid_distElePhot[ieleassocPhot[i]]) > pid.distEB;      
  }else{
    ptiso = pid_hlwTrackElePhot[ieleassocPhot[i]] < ptPhot[i] * pid.trackiso_relEE;
    ecaliso = pid_jurECALElePhot[ieleassocPhot[i]] < ptPhot[i] * pid.ecaliso_relEE;
    hcaliso = pid_twrHCALElePhot[ieleassocPhot[i]] < ptPhot[i] * pid.hcaliso_relEE;
    hoveiso = pid_HoverEElePhot[ieleassocPhot[i]] < pid.hovereisoEE;
    setaeta = pid_etawidElePhot[ieleassocPhot[i]] < pid.setaetaEE;
    deta = pid_detavtxElePhot[ieleassocPhot[i]] < pid.detaEE;
    dphi = pid_dphivtxElePhot[ieleassocPhot[i]] < pid.dphiEE;
    minhits = pid_mishitsElePhot[ieleassocPhot[i]] < pid.minhitsEE;    
    dcot = TMath::Abs(pid_dcotElePhot[ieleassocPhot[i]]) > pid.dcotEE;     
    dist = TMath::Abs(pid_distElePhot[ieleassocPhot[i]]) > pid.distEE;     
  }

  if (vpass) {
    //assert((*vpass).size()==9);
    if((*vpass).size()!=9) { cout << "major failure! (*vpass).size()!=9.. die!" << endl; exit(0) ; }
    (*vpass)[0] = ptiso;
    (*vpass)[1] = ecaliso;
    (*vpass)[2] = hcaliso;
    (*vpass)[3] = hoveiso;
    (*vpass)[4] = setaeta;
    (*vpass)[5] = deta;
    (*vpass)[6] = dphi;
    (*vpass)[7] = minhits;
    (*vpass)[8] = dcot; 
    (*vpass)[9] = dist; 

  }

  return (ptiso && hcaliso && ecaliso && hoveiso && setaeta && deta && minhits && dcot && dist);
}



bool RedNtpTree::cutIDEG(int i, photonidegcuts const& pid, vector<bool> *vpass, bool PU) {

  // Use photon supercluster energy (would be e5x5 if r9>0.93 otherwise)
  bool ptiso = (pid_hlwTrack[i] < ptPhot[i] * pid.trackiso_rel + pid.trackiso_abs);
  bool ecaliso = (pid_jurECAL[i] < ptPhot[i] * pid.ecaliso_rel + pid.ecaliso_abs);
  bool hcaliso = (pid_twrHCAL[i] < ptPhot[i] * pid.hcaliso_rel + pid.hcaliso_abs);
  bool hoveiso = (pid_HoverE[i] < pid.hovereiso);
  bool setaeta = pid_etawid[i] < pid.setaetaEB;

  if(PU){
    if(TMath::Abs(etascPhot[i]) < 1.4442) {
      //    ptiso = (pid_hlwTrack[i] < ptPhot[i] * pid.trackiso_rel + 1.08998 + 8.86335e-02*rhoPF - 1.5 + pid.trackiso_abs);
      ptiso = (pid_hlwTrackNoDz[i] < ptPhot[i] * pid.trackiso_rel + 8.34071e-01 + 5.48136e-01*rhoPF - 1.5 + pid.trackiso_abs);
      ecaliso = (pid_jurECAL[i] < ptPhot[i] * pid.ecaliso_rel + 1.58995 + 2.98677e-01*rhoPF - 2.0 + pid.ecaliso_abs );
      hcaliso = (pid_twrHCAL[i] < ptPhot[i] * pid.hcaliso_rel + 1.49628 + 2.44899e-01*rhoPF - 2.0 + pid.hcaliso_abs );
      hoveiso = (pid_HoverE[i] < 1.96440e-02 + 1.00859e-03*rhoPF - 0.02 + pid.hovereiso);
    }else{
      //    ptiso = (pid_hlwTrack[i] < ptPhot[i] * pid.trackiso_rel + 1.24664 + 7.01932e-02*rhoPF - 1.5 + pid.trackiso_abs);
      ptiso = (pid_hlwTrackNoDz[i] < ptPhot[i] * pid.trackiso_rel + 8.86732e-01 + 5.25491e-01*rhoPF  - 1.5 + pid.trackiso_abs);
      ecaliso = (pid_jurECAL[i] < ptPhot[i] * pid.ecaliso_rel + 8.32333e-01 + 1.91840e-01*rhoPF - 2.0 + pid.ecaliso_abs );
      hcaliso = (pid_twrHCAL[i] < ptPhot[i] * pid.hcaliso_rel + 1.24901 + 2.74598e-01*rhoPF - 2.0 + pid.hcaliso_abs );
      hoveiso = (pid_HoverE[i] < 1.95369e-02 + 1.14826e-03*rhoPF - 0.02 + pid.hovereiso);
    }
  }

  if(TMath::Abs(etascPhot[i]) > 1.4442) {
    setaeta = pid_etawid[i] < pid.setaetaEE;
  }  

  if (vpass) {
    //assert((*vpass).size()==7);
    if((*vpass).size()!=5) { cout << "major failure! (*vpass).size()!=7.. die!" << endl; exit(0) ; }
    (*vpass)[0] = ptiso;
    (*vpass)[1] = ecaliso;
    (*vpass)[2] = hcaliso;
    (*vpass)[3] = hoveiso;
    (*vpass)[4] = setaeta;
  }

  return (ptiso && hcaliso && ecaliso && hoveiso && setaeta);
}






/************************************
 *                                  *
 *                                  *
 *       Cuts in Categories         *
 *                                  *
 *                                  *
 ************************************/




// CiC SELECTION CODE BEGIN - SSIMON
// ---------------------------------------------------------------------------------------------------------------------------------------------
void RedNtpTree::SetPhotonCutsInCategories(phoCiCIDLevel cutlevel, float * cic6_allcuts_lead, float * cic6_allcuts_sublead, float * cic4_allcuts_lead, float * cic4_allcuts_sublead) {

  //thresholds are in this order below
  // isosumoet[]
  // isosumoetbad[]
  // trkisooetom[]
  // sieie[]
  // hovere[]
  // r9[]
  // drtotk_25_99[]
  // pixel[]

// 6 categories
// phoNOCUTS      - thresholds so all photons pass
// phoLOOSE       - sob value=0.0002         - iteration 8 - eff=0.947448  fake=0.0783937
// phoMEDIUM      - sob value=0.0004         - iteration 6 - eff=0.928017  fake=0.0572683
// phoTIGHT       - sob value=0.0008         - iteration 6 - eff=0.895238  fake=0.0392572
// phoSUPERTIGHT  - sob value=0.0016         - iteration 6 - eff=0.849812  fake=0.0256949
// phoHYPERTIGHT1 - sob value=0.0032         - iteration 6 - eff=0.784283  fake=0.016346
// phoHYPERTIGHT4 - sob value=0.025          - iteration 6 - eff=0.41176   fake=0.00217666


// 4 categories
// phoNOCUTS      - thresholds so all photons pass
// phoLOOSE       - sob value=0.0002         - iteration 8 - eff=0.939229  fake=0.0815158
// phoMEDIUM      - sob value=0.0004         - iteration 6 - eff=0.91754   fake=0.0581047
// phoTIGHT       - sob value=0.0008         - iteration 6 - eff=0.886869  fake=0.041063
// phoSUPERTIGHT  - sob value=0.0016         - iteration 6 - eff=0.844314  fake=0.0286033
// phoHYPERTIGHT1 - sob value=0.0032         - iteration 6 - eff=0.774552  fake=0.0191603


  const unsigned int ncuts = 8;
  const unsigned int ncat_cic6 = 6;
  const unsigned int ncat_cic4 = 4;
  switch(cutlevel) {
    case(phoNOCUTS) : {
                        float cic6_allcuts_temp_lead[] = { 
                          1e+09,     1e+09,     1e+09,     1e+09,     1e+09,     1e+09,
                          1e+09,     1e+09,     1e+09,     1e+09,     1e+09,     1e+09,
                          1e+09,     1e+09,     1e+09,     1e+09,     1e+09,     1e+09,
                          1e+09,     1e+09,     1e+09,     1e+09,     1e+09,     1e+09,
                          1e+09,     1e+09,     1e+09,     1e+09,     1e+09,     1e+09,
                          -1e+09,    -1e+09,    -1e+09,    -1e+09,    -1e+09,    -1e+09,
                          -1e+09,    -1e+09,    -1e+09,    -1e+09,    -1e+09,    -1e+09,
                          1.5,         1.5,         1.5,         1.5,         1.5,         1.5 };
                        float cic6_allcuts_temp_sublead[] = { 
                          1e+09,     1e+09,     1e+09,     1e+09,     1e+09,     1e+09,
                          1e+09,     1e+09,     1e+09,     1e+09,     1e+09,     1e+09,
                          1e+09,     1e+09,     1e+09,     1e+09,     1e+09,     1e+09,
                          1e+09,     1e+09,     1e+09,     1e+09,     1e+09,     1e+09,
                          1e+09,     1e+09,     1e+09,     1e+09,     1e+09,     1e+09,
                          -1e+09,    -1e+09,    -1e+09,    -1e+09,    -1e+09,    -1e+09,
                          -1e+09,    -1e+09,    -1e+09,    -1e+09,    -1e+09,    -1e+09,
                          1.5,         1.5,         1.5,         1.5,         1.5,         1.5 };
                        float cic4_allcuts_temp_lead[] = { 
                          1e+09,     1e+09,     1e+09,     1e+09,
                          1e+09,     1e+09,     1e+09,     1e+09,
                          1e+09,     1e+09,     1e+09,     1e+09,
                          1e+09,     1e+09,     1e+09,     1e+09,
                          1e+09,     1e+09,     1e+09,     1e+09,
                          -1e+09,    -1e+09,    -1e+09,    -1e+09,
                          -1e+09,    -1e+09,    -1e+09,    -1e+09,
                          1.5,         1.5,         1.5,         1.5 };
                        float cic4_allcuts_temp_sublead[] = { 
                          1e+09,     1e+09,     1e+09,     1e+09,
                          1e+09,     1e+09,     1e+09,     1e+09,
                          1e+09,     1e+09,     1e+09,     1e+09,
                          1e+09,     1e+09,     1e+09,     1e+09,
                          1e+09,     1e+09,     1e+09,     1e+09,
                          -1e+09,    -1e+09,    -1e+09,    -1e+09,
                          -1e+09,    -1e+09,    -1e+09,    -1e+09,
                          1.5,         1.5,         1.5,         1.5 };
                        for(int i=0;i!=ncuts*ncat_cic6;++i) { cic6_allcuts_lead[i]=cic6_allcuts_temp_lead[i];
                          cic6_allcuts_sublead[i]=cic6_allcuts_temp_sublead[i]; }
                          for(int i=0;i!=ncuts*ncat_cic4;++i) { cic4_allcuts_lead[i]=cic4_allcuts_temp_lead[i];
                            cic4_allcuts_sublead[i]=cic4_allcuts_temp_sublead[i]; }
                      } break;
    case(phoLOOSE) : {
                       float cic6_allcuts_temp_lead[] = { 
                         14.1278,     11.7187,     9.78826,     10.9814,     9.21945,     8.89621,
                         72.5178,     59.1506,     85.1822,     93.8969,     74.2109,     14.4058,
                         7.89015,     5.61652,     4.45536,     5.87563,     4.24725,     2.96206,
                         0.0114196,   0.0109898,   0.0100549,    0.029265,   0.0290002,   0.0279397,
                         0.0907646,   0.0791189,   0.0835245,    0.102617,   0.0596196,    0.098899,
                         0.94,    0.899976,    0.262285,    0.94,    0.90,    0.276953,
                         12.0314,     98.0038,  0.00968623,  0.00636153,  0.00476398,  0.00610842,
                         1.5,         1.5,         1.5,         1.5,         1.5,         1.5 };
                       float cic6_allcuts_temp_sublead[] = { 
                         14.1278,     11.7187,     9.78826,     10.9814,     9.21945,     8.89621,
                         72.5178,     59.1506,     85.1822,     93.8969,     74.2109,     14.4058,
                         7.89015,     5.61652,     4.45536,     5.87563,     4.24725,     2.96206,
                         0.0114196,   0.0109898,   0.0100549,    0.029265,   0.0290002,   0.0279397,
                         0.0907646,   0.0791189,   0.0835245,    0.102617,   0.0596196,    0.098899,
                         0.94,    0.899976,    0.262285,    0.94,    0.90,    0.276953,
                         12.0314,     98.0038,  0.00968623,  0.00636153,  0.00476398,  0.00610842,
                         1.5,         1.5,         1.5,         1.5,         1.5,         1.5 };
//                        float cic4_allcuts_temp_lead[] = { 
//                          8.2468,     4.16322,     5.42672,     2.58495,
//                          66.7249,     69.6142,     84.8978,     7.20773,
//                          7.54866,     4.56393,     5.24503,     2.52207,
//                          0.0111626,   0.0102519,   0.0290281,   0.0282623,
//                          0.0903079,   0.0893613,    0.101156,   0.0730662,
//                          0.94,    0.307556,    0.921175,    0.287359,
//                          3.8468,     1.78477,  0.00621473,  0.00547537,
//                          1.5,         1.5,         1.5,         1.5 };
//                        float cic4_allcuts_temp_sublead[] = { 
//                          8.2468,     4.16322,     5.42672,     2.58495,
//                          66.7249,     69.6142,     84.8978,     7.20773,
//                          7.54866,     4.56393,     5.24503,     2.52207,
//                          0.0111626,   0.0102519,   0.0290281,   0.0282623,
//                          0.0903079,   0.0893613,    0.101156,   0.0730662,
//                          0.94,    0.307556,    0.921175,    0.287359,
//                          3.8468,     1.78477,  0.00621473,  0.00547537,
//                          1.5,         1.5,         1.5,         1.5 };
		       float cic4_allcuts_temp_lead[] = { 
			 8.2,       4.1,       5.4,       2.6,
			 67,        69,        85,       7.2,
			 7.5,       4.5,       5.2,       2.5,
			 0.0112,    0.0102,     0.029,     0.028,
			 0.09,     0.089,     0.101,     0.073, 
			 0.94,      0.31,      0.92,      0.29, 
			 0.26,     0.029,    0.0062,    0.0055,
			 1.5,         1.5,         1.5,         1.5
		       };

                       float cic4_allcuts_temp_sublead[] = {
			 8.2,       4.1,       5.4,       2.6,
			 67,        69,        85,       7.2,
			 7.5,       4.5,       5.2,       2.5,
			 0.0112,    0.0102,     0.029,     0.028,
			 0.09,     0.089,     0.101,     0.073, 
			 0.94,      0.31,      0.92,      0.29, 
			 0.26,     0.029,    0.0062,    0.0055, 
                         1.5,         1.5,         1.5,         1.5 
		       };

                       for(int i=0;i!=ncuts*ncat_cic6;++i) { cic6_allcuts_lead[i]=cic6_allcuts_temp_lead[i];
                         cic6_allcuts_sublead[i]=cic6_allcuts_temp_sublead[i]; }
                         for(int i=0;i!=ncuts*ncat_cic4;++i) { cic4_allcuts_lead[i]=cic4_allcuts_temp_lead[i];
                           cic4_allcuts_sublead[i]=cic4_allcuts_temp_sublead[i]; }
                     } break;
    case(phoMEDIUM) : {
                        float cic6_allcuts_temp_lead[] = {  
                          12.5084,      10.156,     9.23141,     10.0482,     8.34498,     8.73704,
                          70.9011,     50.0742,     21.9926,     24.2436,     18.7884,     12.6882,
                          6.58797,     4.68564,     4.38815,     5.67876,     2.41162,     2.19991,
                          0.0110266,   0.0106749,  0.00983011,   0.0287021,   0.0286817,   0.0272739,
                          0.0891215,   0.0763711,   0.0798623,   0.0911974,   0.0511163,   0.0627764,
                          0.94,    0.90,    0.274434,    0.94,    0.90,    0.276953,
                          96.5654,     98.9721,   0.0119942,   0.0111399,  0.00855448,    0.012159,
                          1.5,         1.5,         1.5,         1.5,         1.5,         1.5 };
                        float cic6_allcuts_temp_sublead[] = {  
                          12.5084,      10.156,     9.23141,     10.0482,     8.34498,     8.73704,
                          70.9011,     50.0742,     21.9926,     24.2436,     18.7884,     12.6882,
                          6.58797,     4.68564,     4.38815,     5.67876,     2.41162,     2.19991,
                          0.0110266,   0.0106749,  0.00983011,   0.0287021,   0.0286817,   0.0272739,
                          0.0891215,   0.0763711,   0.0798623,   0.0911974,   0.0511163,   0.0627764,
                          0.94,    0.90,    0.274434,    0.94,    0.90,    0.276953,
                          96.5654,     98.9721,   0.0119942,   0.0111399,  0.00855448,    0.012159,
                          1.5,         1.5,         1.5,         1.5,         1.5,         1.5 };
//                         float cic4_allcuts_temp_lead[] = {  
//                           6.48191,     3.25515,     3.44913,     2.15119,
//                           63.7204,     11.7537,     12.9501,     3.47679,
//                           6.43647,      3.5981,     3.78298,     2.10611,
//                           0.0108868,   0.0100955,   0.0286898,   0.0278915,
//                           0.0889167,   0.0888774,    0.090159,   0.0610415,
//                           0.94,    0.309849,     0.93949,    0.287359,
//                           96.3363,     1.78477,   0.0109295,   0.0110844,
//                           1.5,         1.5,         1.5,         1.5 };
//                         float cic4_allcuts_temp_sublead[] = {  
//                           6.48191,     3.25515,     3.44913,     2.15119,
//                           63.7204,     11.7537,     12.9501,     3.47679,
//                           6.43647,      3.5981,     3.78298,     2.10611,
//                           0.0108868,   0.0100955,   0.0286898,   0.0278915,
//                           0.0889167,   0.0888774,    0.090159,   0.0610415,
//                           0.94,    0.309849,     0.93949,    0.287359,
//                           96.3363,     1.78477,   0.0109295,   0.0110844,
//                           1.5,         1.5,         1.5,         1.5 };

                        float cic4_allcuts_temp_lead[] = {  
			  6.4,       3.2,       3.4,       2.2,
			  64,      10.8,        13,       3.5,
			  6.4,       3.4,       3.8,       2.1,
			  0.0109,      0.01,     0.029,     0.028 ,
			  0.089,     0.079,      0.09,     0.061 ,
			  0.94,      0.32,      0.94,      0.29,
			  0.98,     0.029,    0.0109,    0.0111,
                          1.5,         1.5,         1.5,         1.5 };
                        float cic4_allcuts_temp_sublead[] = {  
			  6.4,       3.2,       3.4,       2.2,
			  64,      10.8,        13,       3.5,
			  6.4,       3.4,       3.8,       2.1,
			  0.0109,      0.01,     0.029,     0.028 ,
			  0.089,     0.079,      0.09,     0.061 ,
			  0.94,      0.32,      0.94,      0.29,
			  0.98,     0.029,    0.0109,    0.0111,
                          1.5,         1.5,         1.5,         1.5 
			};
                        for(int i=0;i!=ncuts*ncat_cic6;++i) { cic6_allcuts_lead[i]=cic6_allcuts_temp_lead[i];
                          cic6_allcuts_sublead[i]=cic6_allcuts_temp_sublead[i]; }
                          for(int i=0;i!=ncuts*ncat_cic4;++i) { cic4_allcuts_lead[i]=cic4_allcuts_temp_lead[i];
                            cic4_allcuts_sublead[i]=cic4_allcuts_temp_sublead[i]; }
                      } break;
    case(phoTIGHT) : {
                       float cic6_allcuts_temp_lead[] = { 
                         11.1845,     9.28445,     8.98759,     9.19055,     7.94171,     8.16991,
                         70.7835,     16.7873,     13.7361,     15.6259,     13.2407,     10.3932,
                         5.76122,     3.97439,     2.89137,     4.62749,     2.34848,      1.9302,
                         0.010781,   0.0104673,  0.00965497,   0.0284936,    0.028082,   0.0270328,
                         0.0844869,   0.0703749,    0.060775,   0.0881813,   0.0502974,    0.060877,
                         0.94,    0.90,       0.321,    0.94,    0.90,    0.279956,
                         98.9318,     98.9992,   0.0146256,   0.0207672,     34.1809,   0.0261029,
                         1.5,         1.5,         1.5,         1.5,         1.5,         1.5 };
                       float cic6_allcuts_temp_sublead[] = { 
                         11.1845,     9.28445,     8.98759,     9.19055,     7.94171,     8.16991,
                         70.7835,     16.7873,     13.7361,     15.6259,     13.2407,     10.3932,
                         5.76122,     3.97439,     2.89137,     4.62749,     2.34848,      1.9302,
                         0.010781,   0.0104673,  0.00965497,   0.0284936,    0.028082,   0.0270328,
                         0.0844869,   0.0703749,    0.060775,   0.0881813,   0.0502974,    0.060877,
                         0.94,    0.90,       0.321,    0.94,    0.90,    0.279956,
                         98.9318,     98.9992,   0.0146256,   0.0207672,     34.1809,   0.0261029,
                         1.5,         1.5,         1.5,         1.5,         1.5,         1.5 };
//                        float cic4_allcuts_temp_lead[] = { 
//                          4.62434,     2.81143,     2.50887,     1.45731,
//                          62.3241,     5.45376,     7.32095,     2.48286,
//                          4.76431,     2.87144,     3.78298,      1.6321,
//                          0.0107094,  0.00995029,   0.0284777,    0.027055,
//                          0.0865544,   0.0875894,   0.0871544,   0.0499186,
//                          0.94,     0.33192,    0.94,    0.287359,
//                          98.9254,     1.78477,   0.0207514,   0.0278736,
//                          1.5,         1.5,         1.5,         1.5 };
//                        float cic4_allcuts_temp_sublead[] = { 
//                          4.62434,     2.81143,     2.50887,     1.45731,
//                          62.3241,     5.45376,     7.32095,     2.48286,
//                          4.76431,     2.87144,     3.78298,      1.6321,
//                          0.0107094,  0.00995029,   0.0284777,    0.027055,
//                          0.0865544,   0.0875894,   0.0871544,   0.0499186,
//                          0.94,     0.33192,    0.94,    0.287359,
//                          98.9254,     1.78477,   0.0207514,   0.0278736,
//                          1.5,         1.5,         1.5,         1.5 };

                       float cic4_allcuts_temp_lead[] = { 
			 4.7,       2.8,       2.5,      1.46,		
			 62,       5.2,       7.3,       2.5,
			 4.7,       2.9,       3.8,      1.63,
			 0.0107,    0.0099,     0.028,     0.027,
			 0.087,     0.065,     0.087,      0.05,
			 0.94,      0.34,      0.94,      0.29,
			 1,     0.029,     0.021,     0.028,
                         1.5,         1.5,         1.5,         1.5 };
                       float cic4_allcuts_temp_sublead[] = {
			 4.7,       2.8,       2.5,      1.46,		
			 62,       5.2,       7.3,       2.5,
			 4.7,       2.9,       3.8,      1.63,
			 0.0107,    0.0099,     0.028,     0.027,
			 0.087,     0.065,     0.087,      0.05,
			 0.94,      0.34,      0.94,      0.29,
			 1,     0.029,     0.021,     0.028, 
                         1.5,         1.5,         1.5,         1.5 };

                       for(int i=0;i!=ncuts*ncat_cic6;++i) { cic6_allcuts_lead[i]=cic6_allcuts_temp_lead[i];
                         cic6_allcuts_sublead[i]=cic6_allcuts_temp_sublead[i]; }
                         for(int i=0;i!=ncuts*ncat_cic4;++i) { cic4_allcuts_lead[i]=cic4_allcuts_temp_lead[i];
                           cic4_allcuts_sublead[i]=cic4_allcuts_temp_sublead[i]; }
                     } break;
    case(phoSUPERTIGHT) : {
                            float cic6_allcuts_temp_lead[] = { 
                              10.0171,     8.81037,     8.74909,     8.47393,     7.94171,     7.47883,
                              54.9366,     14.3545,     11.5208,      12.939,     10.2496,      9.7095,
                              4.11252,     3.35092,     2.49296,     2.05592,     1.67021,     1.66678,
                              0.0106315,   0.0101656,  0.00950936,   0.0283215,   0.0276216,   0.0263378,
                              0.0823828,   0.0598641,   0.0494497,   0.0706222,   0.0502974,    0.060877,
                              0.94,    0.90,       0.321,    0.94,    0.90,    0.282153,
                              98.9981,          99,   0.0216484,     96.2292,     97.1855,     96.2294,
                              1.5,         1.5,         1.5,         1.5,         1.5,         1.5 };
                            float cic6_allcuts_temp_sublead[] = { 
                              10.0171,     8.81037,     8.74909,     8.47393,     7.94171,     7.47883,
                              54.9366,     14.3545,     11.5208,      12.939,     10.2496,      9.7095,
                              4.11252,     3.35092,     2.49296,     2.05592,     1.67021,     1.66678,
                              0.0106315,   0.0101656,  0.00950936,   0.0283215,   0.0276216,   0.0263378,
                              0.0823828,   0.0598641,   0.0494497,   0.0706222,   0.0502974,    0.060877,
                              0.94,    0.90,       0.321,    0.94,    0.90,    0.282153,
                              98.9981,          99,   0.0216484,     96.2292,     97.1855,     96.2294,
                              1.5,         1.5,         1.5,         1.5,         1.5,         1.5 };
//                             float cic4_allcuts_temp_lead[] = { 
//                               3.77459,     2.18305,     1.76811,     1.30029,
//                               11.5519,     3.48306,     3.87394,     1.89038,
//                               3.47103,      2.1822,     2.26931,     1.43769,
//                               0.0105631,  0.00974116,   0.0282572,   0.0265096,
//                               0.0834648,   0.0821447,   0.0648449,   0.0476437,
//                               0.94,    0.355935,    0.94,    0.316358,
//                               98.9979,     1.94497,     96.2292,     96.2294,
//                               1.5,         1.5,         1.5,         1.5 };
//                             float cic4_allcuts_temp_sublead[] = { 
//                               3.77459,     2.18305,     1.76811,     1.30029,
//                               11.5519,     3.48306,     3.87394,     1.89038,
//                               3.47103,      2.1822,     2.26931,     1.43769,
//                               0.0105631,  0.00974116,   0.0282572,   0.0265096,
//                               0.0834648,   0.0821447,   0.0648449,   0.0476437,
//                               0.94,    0.355935,    0.94,    0.316358,
//                               98.9979,     1.94497,     96.2292,     96.2294,
//                               1.5,         1.5,         1.5,         1.5 };
			    //New agreed rounded numbers PM 2011.06.14
			    float cic4_allcuts_temp_lead[] = {
                              3.8,       2.2,      1.77,      1.29,
                              11.7,       3.4,       3.9,      1.84,
                              3.5,       2.2,       2.3,      1.45,
                              0.0106,    0.0097,     0.028,     0.027,
                              0.082,     0.062,     0.065,     0.048,
                              0.94,      0.36,      0.94,      0.32,
                              1,     0.062,      0.97,      0.97,
                              1.5,         1.5,         1.5,         1.5 };
                            float cic4_allcuts_temp_sublead[] = {
                              3.8,       2.2,      1.77,      1.29,
                              11.7,       3.4,       3.9,      1.84,
                              3.5,       2.2,       2.3,      1.45,
                              0.0106,    0.0097,     0.028,     0.027,
                              0.082,     0.062,     0.065,     0.048,
                              0.94,      0.36,      0.94,      0.32,
                              1,     0.062,      0.97,      0.97,
                              1.5,         1.5,         1.5,         1.5 };

                            for(int i=0;i!=ncuts*ncat_cic6;++i) { cic6_allcuts_lead[i]=cic6_allcuts_temp_lead[i];
                              cic6_allcuts_sublead[i]=cic6_allcuts_temp_sublead[i]; }
                              for(int i=0;i!=ncuts*ncat_cic4;++i) { cic4_allcuts_lead[i]=cic4_allcuts_temp_lead[i];
                                cic4_allcuts_sublead[i]=cic4_allcuts_temp_sublead[i]; }
                          } break;
    case(phoHYPERTIGHT1) : {
                             float cic6_allcuts_temp_lead[] = { 
                               9.14323,     8.13617,     7.43416,     7.97795,     5.88227,     6.60691,
                               16.4126,     10.7813,     10.1764,     11.3829,     8.63128,     8.75289,
                               3.49873,     2.93013,     2.00419,     1.60673,     1.36163,     1.36132,
                               0.0105033,  0.00999387,  0.00946607,   0.0282088,   0.0273334,   0.0256399,
                               0.0782034,   0.0598641,   0.0273668,   0.0553324,   0.0502974,   0.0465477,
                               0.94,    0.90,    0.347653,    0.94,    0.90,    0.301546,
                               98.9999,          99,     1.92089,     98.9224,     98.9492,     98.9224,
                               1.5,         1.5,         1.5,         1.5,         1.5,         1.5 };
                             float cic6_allcuts_temp_sublead[] = { 
                               9.14323,     8.13617,     7.43416,     7.97795,     5.88227,     6.60691,
                               16.4126,     10.7813,     10.1764,     11.3829,     8.63128,     8.75289,
                               3.49873,     2.93013,     2.00419,     1.60673,     1.36163,     1.36132,
                               0.0105033,  0.00999387,  0.00946607,   0.0282088,   0.0273334,   0.0256399,
                               0.0782034,   0.0598641,   0.0273668,   0.0553324,   0.0502974,   0.0465477,
                               0.94,    0.90,    0.347653,    0.94,    0.90,    0.301546,
                               98.9999,          99,     1.92089,     98.9224,     98.9492,     98.9224,
                               1.5,         1.5,         1.5,         1.5,         1.5,         1.5 };
//                              float cic4_allcuts_temp_lead[] = { 
//                                3.16373,     1.76295,     1.38862,     1.10001,
//                                6.07834,      2.7333,     2.84236,    0.769064,
//                                3.36688,     1.86825,     1.66995,     1.43769,
//                                0.0104201,  0.00944203,    0.027691,   0.0245544,
//                                0.0763125,   0.0299019,   0.0469716,   0.0433841,
//                                0.94,    0.409296,    0.94,    0.342725,
//                                98.9999,     96.2831,     98.9224,     98.9224,
//                                1.5,         1.5,         1.5,         1.5 };
//                              float cic4_allcuts_temp_sublead[] = { 
//                                3.16373,     1.76295,     1.38862,     1.10001,
//                                6.07834,      2.7333,     2.84236,    0.769064,
//                                3.36688,     1.86825,     1.66995,     1.43769,
//                                0.0104201,  0.00944203,    0.027691,   0.0245544,
//                                0.0763125,   0.0299019,   0.0469716,   0.0433841,
//                                0.94,    0.409296,    0.94,    0.342725,
//                                98.9999,     96.2831,     98.9224,     98.9224,
//                                1.5,         1.5,         1.5,         1.5 };
                             float cic4_allcuts_temp_lead[] = { 
			       3.2,      1.76,      1.39,      1.18,
			       6.1,       2.7,       2.8,      0.66,
			       3.4,      1.86,      1.67,      1.44,
			       0.0104,    0.0094,     0.028,     0.025,
			       0.076,      0.03,     0.047,     0.046,
			       0.94,      0.41,      0.94,      0.34,
			       1,      0.97,         1,         1 ,
                               1.5,         1.5,         1.5,         1.5 };
                             float cic4_allcuts_temp_sublead[] = { 
			       3.2,      1.76,      1.39,      1.18,
			       6.1,       2.7,       2.8,      0.66,
			       3.4,      1.86,      1.67,      1.44,
			       0.0104,    0.0094,     0.028,     0.025,
			       0.076,      0.03,     0.047,     0.046,
			       0.94,      0.41,      0.94,      0.34,
			       1,      0.97,         1,         1 ,
                               1.5,         1.5,         1.5,         1.5 };
                             for(int i=0;i!=ncuts*ncat_cic6;++i) { cic6_allcuts_lead[i]=cic6_allcuts_temp_lead[i];
                               cic6_allcuts_sublead[i]=cic6_allcuts_temp_sublead[i]; }
                               for(int i=0;i!=ncuts*ncat_cic4;++i) { cic4_allcuts_lead[i]=cic4_allcuts_temp_lead[i];
                                 cic4_allcuts_sublead[i]=cic4_allcuts_temp_sublead[i]; }
                           } break;
     case(phoHYPERTIGHT2) : {
                             float cic6_allcuts_temp_lead[] = { 
                               8.57184,     6.64014,     6.82022,     7.13109,     5.88011,      6.2565,
                               13.4065,     10.4316,     9.18551,     9.30193,     7.51729,     7.30382,
                               2.73319,     2.93013,     1.55723,     1.54876,     1.05254,     1.36132,
                               0.0103615,  0.00978982,  0.00940152,   0.0279141,   0.0260354,   0.0241246,
                               0.0572816,   0.0232443,   0.0173437,   0.0553324,   0.0365276,   0.0465477,
                               0.94,    0.90,    0.367082,    0.94,    0.90,    0.579434,
                               99,          99,     96.2824,     98.9978,     98.9986,     98.9978,
                               1.5,         1.5,         1.5,         1.5,         1.5,         1.5 };
                             float cic6_allcuts_temp_sublead[] = { 
                               8.57184,     6.64014,     6.82022,     7.13109,     5.88011,      6.2565,
                               13.4065,     10.4316,     9.18551,     9.30193,     7.51729,     7.30382,
                               2.73319,     2.93013,     1.55723,     1.54876,     1.05254,     1.36132,
                               0.0103615,  0.00978982,  0.00940152,   0.0279141,   0.0260354,   0.0241246,
                               0.0572816,   0.0232443,   0.0173437,   0.0553324,   0.0365276,   0.0465477,
                               0.94,    0.90,    0.367082,    0.94,    0.90,    0.579434,
                               99,          99,     96.2824,     98.9978,     98.9986,     98.9978,
                               1.5,         1.5,         1.5,         1.5,         1.5,         1.5 };
//                              float cic4_allcuts_temp_lead[] = { 
//                                2.62247,     1.39105,     1.32741,    0.901527,
//                                5.07275,     1.65372,     1.37996,  -0.0456764,
//                                2.71641,      1.6007,     1.55087,     1.43769,
//                                0.0101475,  0.00924974,   0.0273546,   0.0230396,
//                                0.0484683,   0.0189641,   0.0320885,  0.00121448,
//                                0.94,    0.483342,     0.94,    0.551242,
//                                99,     98.9239,     98.9978,     98.9978,
//                                1.5,         1.5,         1.5,         1.5 };
//                              float cic4_allcuts_temp_sublead[] = { 
//                                2.62247,     1.39105,     1.32741,    0.901527,
//                                5.07275,     1.65372,     1.37996,  -0.0456764,
//                                2.71641,      1.6007,     1.55087,     1.43769,
//                                0.0101475,  0.00924974,   0.0273546,   0.0230396,
//                                0.0484683,   0.0189641,   0.0320885,  0.00121448,
//                                0.94,    0.483342,     0.94,    0.551242,
//                                99,     98.9239,     98.9978,     98.9978,
//                                1.5,         1.5,         1.5,         1.5 };
                             float cic4_allcuts_temp_lead[] = { 
			       2.6,      1.31,      1.33,      0.82,
			       5.1,      1.62,      1.38,  -0.224864,
			       2.9,       1.6,      1.55,      1.44,
			       0.0101,    0.0093,     0.027,     0.023,
			       0.048,    0.0189,     0.032,    0.0085,
			       0.94,      0.47,      0.94,      0.52,
			       1,         1,         1,         1,
                               1.5,         1.5,         1.5,         1.5 };
                             float cic4_allcuts_temp_sublead[] = {
			       2.6,      1.31,      1.33,      0.82,
			       5.1,      1.62,      1.38,  -0.224864,
			       2.9,       1.6,      1.55,      1.44,
			       0.0101,    0.0093,     0.027,     0.023,
			       0.048,    0.0189,     0.032,    0.0085,
			       0.94,      0.47,      0.94,      0.52,
			       1,         1,         1,         1, 
                               1.5,         1.5,         1.5,         1.5 };
                             for(int i=0;i!=ncuts*ncat_cic6;++i) { cic6_allcuts_lead[i]=cic6_allcuts_temp_lead[i];
                               cic6_allcuts_sublead[i]=cic6_allcuts_temp_sublead[i]; }
                               for(int i=0;i!=ncuts*ncat_cic4;++i) { cic4_allcuts_lead[i]=cic4_allcuts_temp_lead[i];
                                 cic4_allcuts_sublead[i]=cic4_allcuts_temp_sublead[i]; }
                           } break;
    case(phoHYPERTIGHT3) : {
                             float cic6_allcuts_temp_lead[] = { 
                               7.97897,     6.64014,     6.60332,     5.14765,     5.02192,     5.72775,
                               11.3476,     8.93788,     8.36279,     7.88566,     5.83093,     6.66771,
                               2.348,     2.59173,     1.55158,     1.54876,     0.98618,     1.06927,
                               0.0100676,  0.00971589,  0.00932669,   0.0279141,    0.025781,   0.0229432,
                               0.0372854,   0.0215628,   0.0132992,   0.0412051,   0.0322458,   0.0465477,
                               0.94,    0.90,    0.375623,    0.94,    0.90,    0.579434,
                               99,          99,     98.9239,     98.9999,     98.9997,     98.9987,
                               1.5,         1.5,         1.5,         1.5,         1.5,         1.5 };
                             float cic6_allcuts_temp_sublead[] = { 
                               7.97897,     6.64014,     6.60332,     5.14765,     5.02192,     5.72775,
                               11.3476,     8.93788,     8.36279,     7.88566,     5.83093,     6.66771,
                               2.348,     2.59173,     1.55158,     1.54876,     0.98618,     1.06927,
                               0.0100676,  0.00971589,  0.00932669,   0.0279141,    0.025781,   0.0229432,
                               0.0372854,   0.0215628,   0.0132992,   0.0412051,   0.0322458,   0.0465477,
                               0.94,    0.90,    0.375623,    0.94,    0.90,    0.579434,
                               99,          99,     98.9239,     98.9999,     98.9997,     98.9987,
                               1.5,         1.5,         1.5,         1.5,         1.5,         1.5 };
//                              float cic4_allcuts_temp_lead[] = { 
//                                1.84128,     1.02327,      1.2096,   -0.247273,
//                                3.70571,    0.948406,     1.37996,   -0.756893,
//                                1.94586,     1.53214,      1.4844,    0.111795,
//                                0.00990065,  0.00914404,   0.0272978,   0.0229497,
//                                0.0422667,   0.0185691,    0.022612,  0.000437212,
//                                0.94,    0.692359,    0.965613,    0.551242,
//                                99,     98.9979,     98.9999,     98.9999,
//                                1.5,         1.5,         1.5,         1.5 };
//                              float cic4_allcuts_temp_sublead[] = { 
//                                1.84128,     1.02327,      1.2096,   -0.247273,
//                                3.70571,    0.948406,     1.37996,   -0.756893,
//                                1.94586,     1.53214,      1.4844,    0.111795,
//                                0.00990065,  0.00914404,   0.0272978,   0.0229497,
//                                0.0422667,   0.0185691,    0.022612,  0.000437212,
//                                0.94,    0.692359,    0.965613,    0.551242,
//                                99,     98.9979,     98.9999,     98.9999,
//                                1.5,         1.5,         1.5,         1.5 };

                             float cic4_allcuts_temp_lead[] = { 
			       1.85,      0.96,      1.21,  -0.028513,
			       3.7,      0.97,      1.38,  -0.880416,
			       1.93,       1.4,      1.48,     0.056,
			       0.0099,    0.0092,     0.027,     0.023,
			       0.042,    0.0173,     0.023,    0.0085,
			       0.94,      0.69,      0.97,      0.52,
			       1,         1,         1,         1,
                               1.5,         1.5,         1.5,         1.5 };
                             float cic4_allcuts_temp_sublead[] = { 
			       1.85,      0.96,      1.21,  -0.028513,
			       3.7,      0.97,      1.38,  -0.880416,
			       1.93,       1.4,      1.48,     0.056,
			       0.0099,    0.0092,     0.027,     0.023,
			       0.042,    0.0173,     0.023,    0.0085,
			       0.94,      0.69,      0.97,      0.52,
			       1,         1,         1,         1,
                               1.5,         1.5,         1.5,         1.5 };
                             for(int i=0;i!=ncuts*ncat_cic6;++i) { cic6_allcuts_lead[i]=cic6_allcuts_temp_lead[i];
                               cic6_allcuts_sublead[i]=cic6_allcuts_temp_sublead[i]; }
                               for(int i=0;i!=ncuts*ncat_cic4;++i) { cic4_allcuts_lead[i]=cic4_allcuts_temp_lead[i];
                                 cic4_allcuts_sublead[i]=cic4_allcuts_temp_sublead[i]; }
                           } break;
    case(phoHYPERTIGHT4) : {
                             float cic6_allcuts_temp_lead[] = { 
                               6.53539,     6.07874,     5.51521,     4.78731,     5.00511,     4.90969,
                               9.30747,      8.0574,     7.70153,     7.43339,     5.44326,     6.66771,
                               1.96543,     1.78829,    0.819072,     1.54876,     0.98618,    0.255192,
                               0.0100676,  0.00919753,  0.00911379,   0.0278098,   0.0249354,   0.0221531,
                               0.03099,   0.0153957,   0.0132992,   0.0214415,   0.0322458,   0.0138186,
                               0.94,    0.90,    0.397401,    0.94,    0.90,     0.68715,
                               99,          99,     98.9979,          99,          99,     98.9987,
                               1.5,         1.5,         1.5,         1.5,         1.5,         1.5 };
                             float cic6_allcuts_temp_sublead[] = { 
                               6.53539,     6.07874,     5.51521,     4.78731,     5.00511,     4.90969,
                               9.30747,      8.0574,     7.70153,     7.43339,     5.44326,     6.66771,
                               1.96543,     1.78829,    0.819072,     1.54876,     0.98618,    0.255192,
                               0.0100676,  0.00919753,  0.00911379,   0.0278098,   0.0249354,   0.0221531,
                               0.03099,   0.0153957,   0.0132992,   0.0214415,   0.0322458,   0.0138186,
                               0.94,    0.90,    0.397401,    0.94,    0.90,     0.68715,
                               99,          99,     98.9979,          99,          99,     98.9987,
                               1.5,         1.5,         1.5,         1.5,         1.5,         1.5 };
//                              float cic4_allcuts_temp_lead[] = { 
//                                1.30284,    0.283147,       1.155,   -0.247273,
//                                1.74788,    0.657085,      1.1369,    -1.01141,
//                                1.45901,    0.842083,     1.47572,    0.111795,
//                                0.00976527,   0.0089837,   0.0257602,   0.0229497,
//                                0.0365257,  0.000519817,   0.0197826,  1.22391e-05,
//                                0.94,    0.692359,    0.967428,    0.551242,
//                                99,     98.9999,          99,          99,
//                                1.5,         1.5,         1.5,         1.5 };
//                              float cic4_allcuts_temp_sublead[] = { 
//                                1.30284,    0.283147,       1.155,   -0.247273,
//                                1.74788,    0.657085,      1.1369,    -1.01141,
//                                1.45901,    0.842083,     1.47572,    0.111795,
//                                0.00976527,   0.0089837,   0.0257602,   0.0229497,
//                                0.0365257,  0.000519817,   0.0197826,  1.22391e-05,
//                                0.94,    0.692359,    0.967428,    0.551242,
//                                99,     98.9999,          99,          99,
//                                1.5,         1.5,         1.5,         1.5 };

                             float cic4_allcuts_temp_lead[] = { 
			       1.31,       0.3,      1.15,  -0.028513,
			       1.72,      0.69,      1.14,  -0.880416,
			       1.42,      0.76,      1.48,     0.056,
			       0.0098,     0.009,     0.026,     0.023,
			       0.037,   0.00049,    0.0198,   0.00024,
			       0.94,      0.69,      0.97,      0.73,
			       1,         1,         1,         1,
                               1.5,         1.5,         1.5,         1.5 };
                             float cic4_allcuts_temp_sublead[] = { 
			       1.31,       0.3,      1.15,  -0.028513,
			       1.72,      0.69,      1.14,  -0.880416,
			       1.42,      0.76,      1.48,     0.056,
			       0.0098,     0.009,     0.026,     0.023,
			       0.037,   0.00049,    0.0198,   0.00024,
			       0.94,      0.69,      0.97,      0.73,
			       1,         1,         1,         1,
                               1.5,         1.5,         1.5,         1.5 };
                             for(int i=0;i!=ncuts*ncat_cic6;++i) { cic6_allcuts_lead[i]=cic6_allcuts_temp_lead[i];
                               cic6_allcuts_sublead[i]=cic6_allcuts_temp_sublead[i]; }
                               for(int i=0;i!=ncuts*ncat_cic4;++i) { cic4_allcuts_lead[i]=cic4_allcuts_temp_lead[i];
                                 cic4_allcuts_sublead[i]=cic4_allcuts_temp_sublead[i]; }
                           } break;
  default:std::cout << "UNKNOWN phoCiCIDLevel: " << cutlevel << std::endl;

  }
}


void RedNtpTree::FillPhotonCiCSelectionVariable(int photon_index, int vtx_index)
{
  int photon_category = PhotonCategory(photon_index);

  float val_tkiso = pid_hlwTrack03ForCiC[photon_index][vtx_index];
  float val_ecaliso = pid_jurECAL03[photon_index];
  float val_hcaliso = pid_twrHCAL[photon_index];
  float val_ecalisobad = pid_jurECAL[photon_index];
  float val_hcalisobad = pid_twrHCAL[photon_index];
  float val_tkisobad = 0;
  for(int j=0;j<nvertex;j++)
    if(pid_hlwTrackForCiC[photon_index][j]>val_tkisobad) val_tkisobad = pid_hlwTrackForCiC[photon_index][j];
  float val_sieie = pid_etawid[photon_index];
  float val_hoe = pid_HoverE[photon_index];
  float val_r9 = E9Phot[photon_index]/escRawPhot[photon_index];
  float val_drtotk_25_99 = pid_deltaRToTrackPhot[photon_index];
  float val_pixel = (float)hasPixelSeedPhot[photon_index];

  float isosumconst = 0.;
  float isosumconstbad = 0.;

  float rhofacbad=0.52, rhofac=0.17;
  float val_isosumoet=(val_tkiso+val_ecaliso+val_hcaliso+isosumconst-rhoPF*rhofac)*50./ptPhot[photon_index];
  float val_isosumoetbad=(val_tkisobad+val_ecalisobad+val_hcalisobad+isosumconstbad-rhoPF*rhofacbad)*50./ptPhot[photon_index];
  float val_trkisooet=(val_tkiso)*50./ptPhot[photon_index];

  cic4_cut_isosumoet[photon_category]->Fill(val_isosumoet,weight);
  cic4_cut_isosumoetbad[photon_category]->Fill(val_isosumoetbad,weight);
  cic4_cut_trkisooet[photon_category]->Fill(val_trkisooet,weight);
  cic4_cut_sieie[photon_category]->Fill(val_sieie,weight);
  cic4_cut_hovere[photon_category]->Fill(val_hoe,weight);
  cic4_cut_r9[photon_category]->Fill(val_r9,weight);
  cic4_cut_drtotk_25_99[photon_category]->Fill(val_drtotk_25_99,weight);
  cic4_cut_pixel[photon_category]->Fill(val_pixel,weight);

}



int RedNtpTree::PhotonCiCSelectionLevel( int photon_index , bool electronVeto, int vtx_index) {

  int cutlevelpassed = -1;

  int photon_category = PhotonCategory(photon_index);

  float val_tkiso = pid_hlwTrack03ForCiC[photon_index][vtx_index];
  float val_ecaliso = pid_jurECAL03[photon_index];
  float val_hcaliso = pid_twrHCAL[photon_index];
  float val_ecalisobad = pid_jurECAL[photon_index];
  float val_hcalisobad = pid_twrHCAL[photon_index];
  float val_tkisobad = 0;
  for(int j=0;j<nvertex;j++)
    if(pid_hlwTrackForCiC[photon_index][j]>val_tkisobad) val_tkisobad = pid_hlwTrackForCiC[photon_index][j];
  float val_sieie = pid_etawid[photon_index];
  float val_hoe = pid_HoverE[photon_index];
  float val_r9 = E9Phot[photon_index]/escRawPhot[photon_index];
  float val_drtotk_25_99 = pid_deltaRToTrackPhot[photon_index];
  float val_pixel = (float)hasPixelSeedPhot[photon_index];

  float isosumconst = 0.;
  float isosumconstbad = 0.;

  //PM 2011.05.30 Changed according to new values
  float rhofacbad=0.52, rhofac=0.17;
  float val_isosumoet=(val_tkiso+val_ecaliso+val_hcaliso+isosumconst-rhoPF*rhofac)*50./ptPhot[photon_index];
  float val_isosumoetbad=(val_tkisobad+val_ecalisobad+val_hcalisobad+isosumconstbad-rhoPF*rhofacbad)*50./ptPhot[photon_index];
  float val_trkisooet=(val_tkiso)*50./ptPhot[photon_index];

  bool ph_passcut[phoNCUTLEVELS][8];
  //   float variable[8];
  //  float cut[8];
  for(int iCUTLEVEL=0;iCUTLEVEL!=(int)phoNCUTLEVELS;++iCUTLEVEL) {

    //     variable[0]=val_isosumoet;
    //     variable[1]=val_isosumoetbad;
    //     variable[2]=val_trkisooet;
    //     variable[3]=val_sieie;        
    //     variable[4]=val_hoe;          
    //     variable[5]=val_r9;           
    //     variable[6]=val_drtotk_25_99; 
    //     variable[7]=val_pixel;
    
    //     cut[0]=cic4_cut_lead_isosumoet[iCUTLEVEL][photon_category];
    //     cut[1]=cic4_cut_lead_isosumoetbad[iCUTLEVEL][photon_category];
    //     cut[2]=cic4_cut_lead_trkisooet[iCUTLEVEL][photon_category];
    //     cut[3]=cic4_cut_lead_sieie[iCUTLEVEL][photon_category];        
    //     cut[4]=cic4_cut_lead_hovere[iCUTLEVEL][photon_category];          
    //     cut[5]=cic4_cut_lead_r9[iCUTLEVEL][photon_category];           
    //     cut[6]=cic4_cut_lead_drtotk_25_99[iCUTLEVEL][photon_category]; 
    //     cut[7]=cic4_cut_lead_pixel[iCUTLEVEL][photon_category];        

    ph_passcut[iCUTLEVEL][0] = (val_isosumoet        <=   cic4_cut_lead_isosumoet[iCUTLEVEL][photon_category]     );
    ph_passcut[iCUTLEVEL][1] = (val_isosumoetbad     <=   cic4_cut_lead_isosumoetbad[iCUTLEVEL][photon_category]  );
    ph_passcut[iCUTLEVEL][2] = (val_trkisooet        <=   cic4_cut_lead_trkisooet[iCUTLEVEL][photon_category]     );
    ph_passcut[iCUTLEVEL][3] = (val_sieie            <=   cic4_cut_lead_sieie[iCUTLEVEL][photon_category]         );
    ph_passcut[iCUTLEVEL][4] = (val_hoe              <=   cic4_cut_lead_hovere[iCUTLEVEL][photon_category]        );
    ph_passcut[iCUTLEVEL][5] = (val_r9             >=     cic4_cut_lead_r9[iCUTLEVEL][photon_category]            );// gt cut
    ph_passcut[iCUTLEVEL][6] = electronVeto ? (val_drtotk_25_99   >=     cic4_cut_lead_drtotk_25_99[iCUTLEVEL][photon_category]  ) : true;// gt cut
    ph_passcut[iCUTLEVEL][7] = electronVeto ? (val_pixel            <=   cic4_cut_lead_pixel[iCUTLEVEL][photon_category]         ) : true;
    bool ph_passcut_all = true;
    for(int icut=0;icut!=8;++icut) {
      ph_passcut_all = ph_passcut_all && ph_passcut[iCUTLEVEL][icut];
      if (!ph_passcut[iCUTLEVEL][icut])
	break;
    }
    if(ph_passcut_all) {
      if( cutlevelpassed != iCUTLEVEL - 1 ) {
	std::cerr << "photon " << photon_index << " (category " << photon_category << ") in run/event " << run << "/" << event << " passed CiC cut level " 
		  << iCUTLEVEL << " but not "  << iCUTLEVEL - 1 << ". Did you load your cut values correctly? "<< std::endl;
	/// assert( 0 );
      }
      cutlevelpassed=iCUTLEVEL;
    }
  }

  return cutlevelpassed;

}


