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
#define MAX_PU_REWEIGHT 60
#else
#define MAX_PU_REWEIGHT 60   // for 2012. Was 40 at the end of 2011
#endif

#define MET_SHIFT_2012

#define LEPTONS_2011 0
#define LEPTONS_2012 1

//#define DEBUG
//#define DEBUG1

using std::cout;
using std::endl;

RedNtpTree::RedNtpTree(TTree *tree, const TString& outname) : tree_reader_V8(tree), jsonFile(0) , ptweights_(0), scaleCorrections_(0)
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
  
  tmvaReaderID_Single_Endcap=0;
  tmvaReaderID_Single_Barrel=0;
  // myTree = new TTree("cicTree_structure_","");
  // TString treeVariables = "runCIC/I:eventCIC/I:isosumoet/F:isoecalet/F:isohcalet/F:isotrackeret/F:isosumoetbad/F:isoecaletbad/F:isohcaletbad/F:isotrackeretbad/F:sieie/F:hoe/F:r9/F:drtotk_25_99/F:pixel/F";
  // myTree->Branch("cicTree_structure_",&(tree_.runCIC),treeVariables);
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



vector<int>  RedNtpTree::firstones(Float_t *vec, vector<bool> *asso, int number){

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
  
    for(int j=0; j<number; j++)
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
    //   Long64_t nentries = 1000;
    
    Long64_t nbytes = 0, nb = 0;
    
    TStopwatch timer;

    if (!tmvaReaderID_Single_Barrel || !tmvaReaderID_Single_Endcap)
      SetAllMVA();
    
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
    jetDR = new TH2D("jetDR","jetDR", 20, 0.,100.,100,0,1.);
    jetresp_vs_pt = new TH2D("jetresp_vs_pt","jetresp_vs_pt", 20, 20.,420.,100,-.35,.35);
    jetresp_vs_eta = new TH2D("jetresp_vs_eta","jetresp_vs_eta", 200, -5.,5.,100,-.35,.35);
    jetresp_vs_npu = new TH2D("jetresp_vs_npu","jetresp_vs_npu", 20, 10.,50.,100,-.35,.35);
    jetresp_vs_eta_50 = new TH2D("jetresp_vs_eta_50","jetresp_vs_eta_50", 200, -5.,5.,100,-.35,.35);
    jetresp_vs_eta_50_abs = new TH2D("jetresp_vs_eta_50_abs","jetresp_vs_eta_50_abs", 100, 0.,5.,100,-.35,.35);
    jetresp_vs_npu_50 = new TH2D("jetresp_vs_npu_50","jetresp_vs_npu_50", 20, 10.,50.,100,-.35,.35);
    jetresp_vs_eta_150 = new TH2D("jetresp_vs_eta_150","jetresp_vs_eta_150", 200, -5.,5.,100,-.35,.35);
    jetresp_vs_npu_150 = new TH2D("jetresp_vs_npu_150","jetresp_vs_npu_150", 20, 10.,50.,100,-.35,.35);
    jetresp_vs_pt_forward = new TH2D("jetresp_vs_pt_forward","jetresp_vs_pt_forward", 20, 20.,420.,100,-.35,.35);
    jetresp_vs_npu_forward = new TH2D("jetresp_vs_npu_forward","jetresp_vs_npu_forward", 20, 10.,50.,100,-.35,.35);
 
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
    TH1D ptphotgen1wl("ptphotgen1wl","ptphotgen1wl", 100, 0.,300.);
    TH1D ptphotgen1wh("ptphotgen1wh","ptphotgen1wh", 100, 0.,300.);
    TH1D ptphotgen1zl("ptphotgen1zl","ptphotgen1zl", 100, 0.,300.);
    TH1D ptphotgen1zh("ptphotgen1zh","ptphotgen1zh", 100, 0.,300.);
    TH1D ptphotgen1zn("ptphotgen1zn","ptphotgen1zn", 100, 0.,300.);
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
      
    TH1D* h1_deltaR_jetpart = new TH1D("deltaR_jetpart","",100, 0., 10.);
      

    /********************************************************
     *                                                      *
     *                      TREE INIT                       *
     *                                                      *
     ********************************************************/


    ana_tree = new TTree ("AnaTree","Reduced tree for final analysis") ;
    ana_tree->Branch("run",&runRN,"run/I");
    ana_tree->Branch("event",&eventRN,"event/I");
    ana_tree->Branch("lumi",&lumi,"lumi/I");
    ana_tree->Branch("H_event",&H_event,"H_event/O");
    ana_tree->Branch("V_event",&V_event,"V_event/O");
    ana_tree->Branch("Zbb_event",&Zbb_event,"Zbb_event/O");
    ana_tree->Branch("Vqq_event",&Vqq_event,"Vqq_event/O");
    ana_tree->Branch("rhoPF",&rhoPFRN,"rhoPF/F");
    ana_tree->Branch("massgg",&massgg,"massgg/F");
    ana_tree->Branch("ptgg",&ptgg,"ptgg/F");
    ana_tree->Branch("ptggnewvtx",&ptggnewvtx,"ptggnewvtx/F");
    ana_tree->Branch("phigg",&phigg,"phigg/F");
    ana_tree->Branch("etagg",&etagg,"etagg/F");
    ana_tree->Branch("massggnewvtx",&massggnewvtx,"massggnewvtx/F");
    ana_tree->Branch("ptphot1",&ptphot1,"ptphot1/F");
    ana_tree->Branch("ptphot2",&ptphot2,"ptphot2/F");
    ana_tree->Branch("deltaRToTrackphot1",&deltaRToTrackphot1,"deltaRToTrackphot1/F");
    ana_tree->Branch("deltaRToTrackphot2",&deltaRToTrackphot2,"deltaRToTrackphot2/F");
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
    ana_tree->Branch("energyErrphot1",&energyErrphot1,"energyErrphot1/F");
    ana_tree->Branch("energyErrphot2",&energyErrphot2,"energyErrphot2/F");
    ana_tree->Branch("energySmearingphot1",&energySmearingphot1,"energySmearingphot1/F");
    ana_tree->Branch("energySmearingphot2",&energySmearingphot2,"energySmearingphot2/F");
    ana_tree->Branch("r9phot1",&r9phot1,"r9phot1/F");
    ana_tree->Branch("r9phot2",&r9phot2,"r9phot2/F");
    ana_tree->Branch("isemEGphot1",&isemEGphot1,"isemEGphot1/I");
    ana_tree->Branch("isemEGphot2",&isemEGphot2,"isemEGphot2/I");
    ana_tree->Branch("promptGamma",&promptGamma,"promptGamma/I");
    ana_tree->Branch("LOGamma",    &LOGamma,    "LOGamma/I");
    ana_tree->Branch("ISRGamma",   &ISRGamma,   "ISRGamma/I");
    ana_tree->Branch("FSRGamma",   &FSRGamma,   "FSRGamma/I");

//     ana_tree->Branch("idloosenewEGphot1",&idloosenewEGphot1,"idloosenewEGphot1/I");
//     ana_tree->Branch("idloosenewEGphot2",&idloosenewEGphot2,"idloosenewEGphot2/I");
//     ana_tree->Branch("idloose006newEGphot1",&idloose006newEGphot1,"idloose006newEGphot1/I");
//     ana_tree->Branch("idloose006newEGphot2",&idloose006newEGphot2,"idloose006newEGphot2/I");
//     ana_tree->Branch("idtightnewEGphot1",&idtightnewEGphot1,"idtightnewEGphot1/I");
//     ana_tree->Branch("idtightnewEGphot2",&idtightnewEGphot2,"idtightnewEGphot2/I");
//     ana_tree->Branch("idhggtightnewEGphot1",&idhggtightnewEGphot1,"idhggtightnewEGphot1/I");
//     ana_tree->Branch("idhggtightnewEGphot2",&idhggtightnewEGphot2,"idhggtightnewEGphot2/I");
//     ana_tree->Branch("idloosenewpuEGphot1",&idloosenewpuEGphot1,"idloosenewpuEGphot1/I");
//     ana_tree->Branch("idloosenewpuEGphot2",&idloosenewpuEGphot2,"idloosenewpuEGphot2/I");
//     ana_tree->Branch("idtightnewpuEGphot1",&idtightnewpuEGphot1,"idtightnewpuEGphot1/I");
//     ana_tree->Branch("idtightnewpuEGphot2",&idtightnewpuEGphot2,"idtightnewpuEGphot2/I");
//     ana_tree->Branch("idhggtightnewpuEGphot1",&idhggtightnewpuEGphot1,"idhggtightnewpuEGphot1/I");
//     ana_tree->Branch("idhggtightnewpuEGphot2",&idhggtightnewpuEGphot2,"idhggtightnewpuEGphot2/I");
    ana_tree->Branch("idmvaphot1",&idmvaphot1,"idmvaphot1/F");
    ana_tree->Branch("idmvaphot2",&idmvaphot2,"idmvaphot2/F");
    ana_tree->Branch("idcicphot1",&idcicphot1,"idcicphot1/I");
    ana_tree->Branch("idcicphot2",&idcicphot2,"idcicphot2/I");
    ana_tree->Branch("idcicnoelvetophot1",&idcicnoelvetophot1,"idcicnoelvetophot1/I");
    ana_tree->Branch("idcicnoelvetophot2",&idcicnoelvetophot2,"idcicnoelvetophot2/I");
    ana_tree->Branch("idcicpfphot1",&idcicpfphot1,"idcicpfphot1/I");
    ana_tree->Branch("idcicpfphot2",&idcicpfphot2,"idcicpfphot2/I");
    ana_tree->Branch("idcicpfnoelvetophot1",&idcicpfnoelvetophot1,"idcicpfnoelvetophot1/I");
    ana_tree->Branch("idcicpfnoelvetophot2",&idcicpfnoelvetophot2,"idcicpfnoelvetophot2/I");
//     ana_tree->Branch("idlooseEGphot1",&idlooseEGphot1,"idlooseEGphot1/I");
//     ana_tree->Branch("idlooseEGphot2",&idlooseEGphot2,"idlooseEGphot2/I");
//     ana_tree->Branch("idtightEGphot1",&idtightEGphot1,"idtightEGphot1/I");
//     ana_tree->Branch("idtightEGphot2",&idtightEGphot2,"idtightEGphot2/I");
//     ana_tree->Branch("idloosephot1",&idloosephot1,"idloosephot1/I");
//     ana_tree->Branch("idloosephot2",&idloosephot2,"idloosephot2/I");
//     ana_tree->Branch("idmediumphot1",&idmediumphot1,"idmediumphot1/I");
//     ana_tree->Branch("idmediumphot2",&idmediumphot2,"idmediumphot2/I");
//     ana_tree->Branch("idloosecsphot1",&idloosecsphot1,"idloosecsphot1/I"); 
//     ana_tree->Branch("idloosecsphot2",&idloosecsphot2,"idloosecsphot2/I"); 
//     ana_tree->Branch("idmediumcsphot1",&idmediumcsphot1,"idmediumcsphot1/I"); 
//     ana_tree->Branch("idmediumcsphot2",&idmediumcsphot2,"idmediumcsphot2/I"); 
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
    
    ana_tree->Branch("njets", &njets, "njets/I");
    ana_tree->Branch("ecorrjet",  ecorrjet,  "ecorrjet[njets]/F");
    ana_tree->Branch("ptjet",  ptjet,  "ptjet[njets]/F");
    ana_tree->Branch("ptcorrjet",  ptcorrjet,  "ptcorrjet[njets]/F");
    ana_tree->Branch("etajet", etajet, "etajet[njets]/F");
    ana_tree->Branch("phijet", phijet, "phijet[njets]/F");
    ana_tree->Branch("betajet", betajet, "betajet[njets]/F");
    ana_tree->Branch("betastarjet", betastarjet, "betastarjet[njets]/F");
    ana_tree->Branch("btagvtxjet", btagvtxjet, "btagvtxjet[njets]/F");
    ana_tree->Branch("btagjprobjet", btagjprobjet, "btagjprobjet[njets]/F");
    ana_tree->Branch("ptDjet", ptDjet, "ptDjet[njets]/F");
    ana_tree->Branch("rmsjet", rmsjet, "rmsjet[njets]/F");
    ana_tree->Branch("ntrkjet", ntrkjet, "ntrkjet[njets]/F");
    ana_tree->Branch("nneutjet", nneutjet, "nneutjet[njets]/F");
    ana_tree->Branch("jetIdSimple_mvajet", jetIdSimple_mvajet, "jetIdSimple_mvajet[njets]/F");
    ana_tree->Branch("jetIdFull_mvajet", jetIdFull_mvajet, "jetIdFull_mvajet[njets]/F");
    ana_tree->Branch("jetId_dR2Meanjet", jetId_dR2Meanjet, "jetId_dR2Meanjet[njets]/F");
    ana_tree->Branch("jetId_betaStarClassicjet", jetId_betaStarClassicjet, "jetId_betaStarClassicjet[njets]/F");
    ana_tree->Branch("jetId_frac01jet", jetId_frac01jet, "jetId_frac01jet[njets]/F");
    ana_tree->Branch("jetId_frac02jet", jetId_frac02jet, "jetId_frac02jet[njets]/F");
    ana_tree->Branch("jetId_frac03jet", jetId_frac03jet, "jetId_frac03jet[njets]/F");
    ana_tree->Branch("jetId_frac04jet", jetId_frac04jet, "jetId_frac04jet[njets]/F");
    ana_tree->Branch("jetId_frac05jet", jetId_frac05jet, "jetId_frac05jet[njets]/F");
    ana_tree->Branch("jetId_betajet", jetId_betajet, "jetId_betajet[njets]/F");
    ana_tree->Branch("jetId_betaStarjet", jetId_betaStarjet, "jetId_betaStarjet[njets]/F");
    ana_tree->Branch("jetIdCutBased_wpjet", jetIdCutBased_wpjet, "jetIdCutBased_wpjet[njets]/F");
    ana_tree->Branch("jetIdSimple_wpjet", jetIdSimple_wpjet, "jetIdSimple_wpjet[njets]/F");
    ana_tree->Branch("jetIdFull_wpjet", jetIdFull_wpjet, "jetIdFull_wpjet[njets]/F");
    ana_tree->Branch("assjet",assjet,"assjet[njets]/I");
    ana_tree->Branch("partPdgIDjet",partPdgIDjet,"partPdgIDjet[njets]/I");
    ana_tree->Branch("partMomPdgIDjet",partMomPdgIDjet,"partMomPdgIDjet[njets]/I");

//  ana_tree->Branch("assjet2",&assjet2,"assjet2/I");
//  ana_tree->Branch("deltaeta",&deltaeta,"deltaeta/F");
//  ana_tree->Branch("zeppenjet",&zeppenjet,"zeppenjet/F");
//  ana_tree->Branch("deltaphi",&deltaphi,"deltaphi/F");
//  ana_tree->Branch("deltaphinewvtx",&deltaphinewvtx,"deltaphinewvtx/F");
//  ana_tree->Branch("deltaphigg",&deltaphigg,"deltaphigg/F");
//  ana_tree->Branch("invmassjet",&invmassjet,"invmassjet/F");
//  ana_tree->Branch("invmass2g1j",&invmass2g1j,"invmass2g1j/F");
//  ana_tree->Branch("invmass2g2j",&invmass2g2j,"invmass2g2j/F");
//  ana_tree->Branch("pt2g2j",&pt2g2j,"pt2g2j/F");
//  ana_tree->Branch("eta2j",&eta2j,"eta2j/F");
//  ana_tree->Branch("phi2j",&phi2j,"phi2j/F");
//  ana_tree->Branch("pt2j",&pt2j,"pt2j/F");

    ana_tree->Branch("nvtx",&nvtx,"nvtx/F");

    ana_tree->Branch("vtxId",&vtxId,"vtxId/I");
    ana_tree->Branch("vtxPos_x",&vtxPos_x,"vtxPos_x/F");
    ana_tree->Branch("vtxPos_y",&vtxPos_y,"vtxPos_y/F");
    ana_tree->Branch("vtxPos_z",&vtxPos_z,"vtxPos_z/F");
    ana_tree->Branch("vtxIdMVA",&vtxIdMVA,"vtxIdMVA/F");
    ana_tree->Branch("vtxIdEvtProb",&vtxIdEvtProb,"vtxIdEvtProb/F");

    
    // ana_tree->Branch("met",&met,"met/F");
    // ana_tree->Branch("phimet",&phimet,"phimet/F");
    
    ana_tree->Branch("sMet", &sMet_, "sMet/F")  ;
    ana_tree->Branch("eMet", &eMet_, "eMet/F")  ;
    ana_tree->Branch("phiMet", &phiMet_, "phiMet/F");
    ana_tree->Branch("signifMet", &signifMet_, "signifMet/F");
    ana_tree->Branch("eSmearedMet",&eSmearedMet_,"eSmearedMet/F");
    ana_tree->Branch("phiSmearedMet",&phiSmearedMet_,"phiSmearedMet/F");
    ana_tree->Branch("eShiftedMet",&eShiftedMet_,"eShiftedMet/F");
    ana_tree->Branch("phiShiftedMet",&phiShiftedMet_,"phiShiftedMet/F");
    ana_tree->Branch("eShiftedScaledMet",&eShiftedScaledMet_,"eShiftedScaledMet/F");
    ana_tree->Branch("phiShiftedScaledMet",&phiShiftedScaledMet_,"phiShiftedScaledMet/F");
    ana_tree->Branch("eSmearedShiftedMet",&eSmearedShiftedMet_,"eSmearedShiftedMet/F");
    ana_tree->Branch("phiSmearedShiftedMet",&phiSmearedShiftedMet_,"phiSmearedShiftedMet/F");
    ana_tree->Branch("eShiftedScaledMetPUcorr",&eShiftedScaledMetPUcorr_,"eShiftedScaledMetPUcorr/F");
    ana_tree->Branch("phiShiftedScaledMetPUcorr",&phiShiftedScaledMetPUcorr_,"phiShiftedScaledMetPUcorr/F");
    ana_tree->Branch("eSmearedShiftedMePUcorrt",&eSmearedShiftedMetPUcorr_,"eSmearedShiftedMetPUcorr/F");
    ana_tree->Branch("phiSmearedShiftedMetPUcorr",&phiSmearedShiftedMetPUcorr_,"phiSmearedShiftedMetPUcorr/F");
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
    
  //ana_tree->Branch("gen_pt_genjet1",      &gen_pt_genjet1,      "gen_pt_genjet1/F");         
  //ana_tree->Branch("gen_pt_genjet2",      &gen_pt_genjet2,      "gen_pt_genjet2/F");       
  //ana_tree->Branch("gen_eta_genjet1",     &gen_eta_genjet1,     "gen_eta_genjet1/F");         
  //ana_tree->Branch("gen_eta_genjet2",     &gen_eta_genjet2,     "gen_eta_genjet2/F");        
  //ana_tree->Branch("gen_phi_genjet1",     &gen_phi_genjet1,     "gen_phi_genjet1/F");        
  //ana_tree->Branch("gen_phi_genjet2",     &gen_phi_genjet2,     "gen_phi_genjet2/F");         
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
    ana_tree->Branch("gen_pt_lep1",      &gen_pt_lep1,      "gen_pt_lep1/F");         
    ana_tree->Branch("gen_pt_lep2",      &gen_pt_lep2,      "gen_pt_lep2/F");         
    ana_tree->Branch("gen_eta_lep1",     &gen_eta_lep1,     "gen_eta_lep1/F");         
    ana_tree->Branch("gen_eta_lep2",     &gen_eta_lep2,     "gen_eta_lep2/F");         
    ana_tree->Branch("gen_phi_lep1",     &gen_phi_lep1,     "gen_phi_lep1/F");         
    ana_tree->Branch("gen_phi_lep2",     &gen_phi_lep2,     "gen_phi_lep2/F");         
    ana_tree->Branch("gen_pid_lep1",     &gen_pid_lep1,     "gen_pid_lep1/I");         
    ana_tree->Branch("gen_pid_lep2",     &gen_pid_lep2,     "gen_pid_lep2/I");         


    // tight selected electrons
    ana_tree->Branch("ptele1",    &ptele1,    "ptele1/F");
    ana_tree->Branch("ptele2",    &ptele2,    "ptele2/F");
    ana_tree->Branch("etaele1",   &etaele1,   "etaele1/F");
    ana_tree->Branch("etaele2",   &etaele2,   "etaele2/F");
    ana_tree->Branch("phiele1",   &phiele1,   "phiele1/F");
    ana_tree->Branch("phiele2",   &phiele2,   "phiele2/F");
    ana_tree->Branch("eneele1",   &eneele1,   "eneele1/F");
    ana_tree->Branch("eneele2",   &eneele2,   "eneele2/F");
    ana_tree->Branch("sIeIeele1", &sIeIeele1, "sIeIeele1/F");
    ana_tree->Branch("sIeIeele2", &sIeIeele2, "sIeIeele2/F");
    ana_tree->Branch("dphiele1",  &dphiele1,  "dphiele1/F");
    ana_tree->Branch("dphiele2",  &dphiele2,  "dphiele2/F");
    ana_tree->Branch("detaele1",  &detaele1,  "detaele1/F");
    ana_tree->Branch("detaele2",  &detaele2,  "detaele2/F");
    ana_tree->Branch("mhitsele1", &mhitsele1, "mhitsele1/I");
    ana_tree->Branch("mhitsele2", &mhitsele2, "mhitsele2/I");
    ana_tree->Branch("d0ele1",    &d0ele1,    "d0ele1/F");
    ana_tree->Branch("d0ele2",    &d0ele2,    "d0ele2/F");
    ana_tree->Branch("dzele1",    &dzele1,    "dzele1/F");
    ana_tree->Branch("dzele2",    &dzele2,    "dzele2/F");
    ana_tree->Branch("invMassele1g1",   &invMassele1g1,   "invMassele1g1/F");
    ana_tree->Branch("invMassele1g2",   &invMassele1g2,   "invMassele1g2/F");
    ana_tree->Branch("invMassele2g1",   &invMassele2g1,   "invMassele2g1/F");
    ana_tree->Branch("invMassele2g2",   &invMassele2g2,   "invMassele2g2/F");

    if (LEPTONS_2011) {
      ana_tree->Branch("dcotele1",    &dcotele1,    "dcotele1/F");
      ana_tree->Branch("dcotele2",    &dcotele2,    "dcotele2/F");
      ana_tree->Branch("distele1",    &distele1,    "distele1/F");
      ana_tree->Branch("distele2",    &distele2,    "distele2/F");
      ana_tree->Branch("isoele1",     &isoele1,     "isoele1/F");
      ana_tree->Branch("isoele2",     &isoele2,     "isoele2/F");
      ana_tree->Branch("fullisoele1", &fullisoele1, "fullisoele1/F");
      ana_tree->Branch("fullisoele2", &fullisoele2, "fullisoele2/F");
    }
    if (LEPTONS_2012) {
      ana_tree->Branch("oEmoPele1",      &oEmoPele1,      "oEmoPele1/F");
      ana_tree->Branch("oEmoPele2",      &oEmoPele2,      "oEmoPele2/F");
      ana_tree->Branch("mvanotrigele1",  &mvanotrigele1,  "mvanotrigele1/F");
      ana_tree->Branch("mvanotrigele2",  &mvanotrigele2,  "mvanotrigele2/F");
      ana_tree->Branch("mvatrigele1",    &mvatrigele1,    "mvatrigele1/F");
      ana_tree->Branch("mvatrigele2",    &mvatrigele2,    "mvatrigele2/F");
      ana_tree->Branch("matchconvele1",  &matchconvele1,  "matchconvele1/I");
      ana_tree->Branch("matchconvele2",  &matchconvele2,  "matchconvele2/I");
      ana_tree->Branch("chHadIso03ele1", &chHadIso03ele1, "chHadIso03ele1/F");
      ana_tree->Branch("chHadIso03ele2", &chHadIso03ele2, "chHadIso03ele2/F");
      ana_tree->Branch("nHadIso03ele1",  &nHadIso03ele1,  "nHadIso03ele1/F");
      ana_tree->Branch("nHadIso03ele2",  &nHadIso03ele2,  "nHadIso03ele2/F");
      ana_tree->Branch("photIso03ele1",  &photIso03ele1,  "photIso03ele1/F");
      ana_tree->Branch("photIso03ele2",  &photIso03ele2,  "photIso03ele2/F");
    }


    // loose selected electrons
    ana_tree->Branch("pteleloose1",    &pteleloose1,    "pteleloose1/F");
    ana_tree->Branch("pteleloose2",    &pteleloose2,    "pteleloose2/F");
    ana_tree->Branch("etaeleloose1",   &etaeleloose1,   "etaeleloose1/F");
    ana_tree->Branch("etaeleloose2",   &etaeleloose2,   "etaeleloose2/F");
    ana_tree->Branch("phieleloose1",   &phieleloose1,   "phieleloose1/F");
    ana_tree->Branch("phieleloose2",   &phieleloose2,   "phieleloose2/F");
    ana_tree->Branch("eneeleloose1",   &eneeleloose1,   "eneeleloose1/F");
    ana_tree->Branch("eneeleloose2",   &eneeleloose2,   "eneeleloose2/F");
    ana_tree->Branch("sIeIeeleloose1", &sIeIeeleloose1, "sIeIeeleloose1/F");
    ana_tree->Branch("sIeIeeleloose2", &sIeIeeleloose2, "sIeIeeleloose2/F");
    ana_tree->Branch("dphieleloose1",  &dphieleloose1,  "dphieleloose1/F");
    ana_tree->Branch("dphieleloose2",  &dphieleloose2,  "dphieleloose2/F");
    ana_tree->Branch("detaeleloose1",  &detaeleloose1,  "detaeleloose1/F");
    ana_tree->Branch("detaeleloose2",  &detaeleloose2,  "detaeleloose2/F");
    ana_tree->Branch("mhitseleloose1", &mhitseleloose1, "mhitseleloose1/I");
    ana_tree->Branch("mhitseleloose2", &mhitseleloose2, "mhitseleloose2/I");
    ana_tree->Branch("d0eleloose1",    &d0eleloose1,    "d0eleloose1/F");
    ana_tree->Branch("d0eleloose2",    &d0eleloose2,    "d0eleloose2/F");
    ana_tree->Branch("dzeleloose1",    &dzeleloose1,    "dzeleloose1/F");
    ana_tree->Branch("dzeleloose2",    &dzeleloose2,    "dzeleloose2/F");
    ana_tree->Branch("invMasseleloose1g1",   &invMasseleloose1g1,   "invMasseleloose1g1/F");
    ana_tree->Branch("invMasseleloose1g2",   &invMasseleloose1g2,   "invMasseleloose1g2/F");
    ana_tree->Branch("invMasseleloose2g1",   &invMasseleloose2g1,   "invMasseleloose2g1/F");
    ana_tree->Branch("invMasseleloose2g2",   &invMasseleloose2g2,   "invMasseleloose2g2/F");
    if (LEPTONS_2011) {
      ana_tree->Branch("dcoteleloose1",    &dcoteleloose1,    "dcoteleloose1/F");
      ana_tree->Branch("dcoteleloose2",    &dcoteleloose2,    "dcoteleloose2/F");
      ana_tree->Branch("disteleloose1",    &disteleloose1,    "disteleloose1/F");
      ana_tree->Branch("disteleloose2",    &disteleloose2,    "disteleloose2/F");
      ana_tree->Branch("isoeleloose1",     &isoeleloose1,     "isoeleloose1/F");
      ana_tree->Branch("isoeleloose2",     &isoeleloose2,     "isoeleloose2/F");
      ana_tree->Branch("fullisoeleloose1", &fullisoeleloose1, "fullisoeleloose1/F");
      ana_tree->Branch("fullisoeleloose2", &fullisoeleloose2, "fullisoeleloose2/F");
    }
    if (LEPTONS_2012) {
      ana_tree->Branch("oEmoPeleloose1",      &oEmoPeleloose1,      "oEmoPeleloose1/F");
      ana_tree->Branch("oEmoPeleloose2",      &oEmoPeleloose2,      "oEmoPeleloose2/F");
      ana_tree->Branch("mvanotrigeleloose1",  &mvanotrigeleloose1,  "mvanotrigeleloose1/F");
      ana_tree->Branch("mvanotrigeleloose2",  &mvanotrigeleloose2,  "mvanotrigeleloose2/F");
      ana_tree->Branch("mvatrigeleloose1",    &mvatrigeleloose1,    "mvatrigeleloose1/F");
      ana_tree->Branch("mvatrigeleloose2",    &mvatrigeleloose2,    "mvatrigeleloose2/F");
      ana_tree->Branch("matchconveleloose1",  &matchconveleloose1,  "matchconveleloose1/I");
      ana_tree->Branch("matchconveleloose2",  &matchconveleloose2,  "matchconveleloose2/I");
      ana_tree->Branch("chHadIso03eleloose1", &chHadIso03eleloose1, "chHadIso03eleloose1/F");
      ana_tree->Branch("chHadIso03eleloose2", &chHadIso03eleloose2, "chHadIso03eleloose2/F");
      ana_tree->Branch("nHadIso03eleloose1",  &nHadIso03eleloose1,  "nHadIso03eleloose1/F");
      ana_tree->Branch("nHadIso03eleloose2",  &nHadIso03eleloose2,  "nHadIso03eleloose2/F");
      ana_tree->Branch("photIso03eleloose1",  &photIso03eleloose1,  "photIso03eleloose1/F");
      ana_tree->Branch("photIso03eleloose2",  &photIso03eleloose2,  "photIso03eleloose2/F");
    }


    // tight selected muons
    ana_tree->Branch("ptmu1",      &ptmu1,      "ptmu1/F");
    ana_tree->Branch("ptmu2",      &ptmu2,      "ptmu2/F");
    ana_tree->Branch("etamu1",     &etamu1,     "etamu1/F");
    ana_tree->Branch("etamu2",     &etamu2,     "etamu2/F");
    ana_tree->Branch("phimu1",     &phimu1,     "phimu1/F");
    ana_tree->Branch("phimu2",     &phimu2,     "phimu2/F");
    ana_tree->Branch("enemu1",     &enemu1,     "enemu1/F");
    ana_tree->Branch("enemu2",     &enemu2,     "enemu2/F");
    ana_tree->Branch("pixhitsmu1", &pixhitsmu1, "pixhitsmu1/I");
    ana_tree->Branch("pixhitsmu2", &pixhitsmu2, "pixhitsmu2/I");
    ana_tree->Branch("trkhitsmu1", &trkhitsmu1, "trkhitsmu1/I");
    ana_tree->Branch("trkhitsmu2", &trkhitsmu2, "trkhitsmu2/I");
    ana_tree->Branch("hitsmu1",    &hitsmu1,    "hitsmu1/I");
    ana_tree->Branch("hitsmu2",    &hitsmu2,    "hitsmu2/I");
    ana_tree->Branch("chi2mu1",    &chi2mu1,    "chi2mu1/F");
    ana_tree->Branch("chi2mu2",    &chi2mu2,    "chi2mu2/F");
    ana_tree->Branch("matchmu1",   &matchmu1,   "matchmu1/I");
    ana_tree->Branch("matchmu2",   &matchmu2,   "matchmu2/I");
    ana_tree->Branch("d0mu1",      &d0mu1,      "d0mu1/F");
    ana_tree->Branch("d0mu2",      &d0mu2,      "d0mu2/F");
    ana_tree->Branch("dzmu1",      &dzmu1,      "dzmu1/F");
    ana_tree->Branch("dzmu2",      &dzmu2,      "dzmu2/F");

    if (LEPTONS_2011) {
      ana_tree->Branch("isomu1",     &isomu1,     "isomu1/F");
      ana_tree->Branch("isomu2",     &isomu2,     "isomu2/F");
    }
    if (LEPTONS_2012) {
      ana_tree->Branch("chHadmu1",   &chHadmu1,   "chHadmu1/F");
      ana_tree->Branch("chHadmu2",   &chHadmu2,   "chHadmu2/F");
      ana_tree->Branch("nHadmu1",    &nHadmu1,    "nHadmu1/F");
      ana_tree->Branch("nHadmu2",    &nHadmu2,    "nHadmu2/F");
      ana_tree->Branch("photmu1",    &photmu1,    "photmu1/F");
      ana_tree->Branch("photmu2",    &photmu2,    "photmu2/F");
      ana_tree->Branch("puptmu1",    &puptmu1,    "puptmu1/F");
      ana_tree->Branch("puptmu2",    &puptmu2,    "puptmu2/F");
    }

    // loose selected muons
    ana_tree->Branch("ptmuloose1",      &ptmuloose1,      "ptmuloose1/F");
    ana_tree->Branch("ptmuloose2",      &ptmuloose2,      "ptmuloose2/F");
    ana_tree->Branch("etamuloose1",     &etamuloose1,     "etamuloose1/F");
    ana_tree->Branch("etamuloose2",     &etamuloose2,     "etamuloose2/F");
    ana_tree->Branch("phimuloose1",     &phimuloose1,     "phimuloose1/F");
    ana_tree->Branch("phimuloose2",     &phimuloose2,     "phimuloose2/F");
    ana_tree->Branch("enemuloose1",     &enemuloose1,     "enemuloose1/F");
    ana_tree->Branch("enemuloose2",     &enemuloose2,     "enemuloose2/F");
    ana_tree->Branch("pixhitsmuloose1", &pixhitsmuloose1, "pixhitsmuloose1/I");
    ana_tree->Branch("pixhitsmuloose2", &pixhitsmuloose2, "pixhitsmuloose2/I");
    ana_tree->Branch("trkhitsmuloose1", &trkhitsmuloose1, "trkhitsmuloose1/I");
    ana_tree->Branch("trkhitsmuloose2", &trkhitsmuloose2, "trkhitsmuloose2/I");
    ana_tree->Branch("hitsmuloose1",    &hitsmuloose1,    "hitsmuloose1/I");
    ana_tree->Branch("hitsmuloose2",    &hitsmuloose2,    "hitsmuloose2/I");
    ana_tree->Branch("chi2muloose1",    &chi2muloose1,    "chi2muloose1/F");
    ana_tree->Branch("chi2muloose2",    &chi2muloose2,    "chi2muloose2/F");
    ana_tree->Branch("matchmuloose1",   &matchmuloose1,   "matchmuloose1/I");
    ana_tree->Branch("matchmuloose2",   &matchmuloose2,   "matchmuloose2/I");
    ana_tree->Branch("d0muloose1",      &d0muloose1,      "d0muloose1/F");
    ana_tree->Branch("d0muloose2",      &d0muloose2,      "d0muloose2/F");
    ana_tree->Branch("dzmuloose1",      &dzmuloose1,      "dzmuloose1/F");
    ana_tree->Branch("dzmuloose2",      &dzmuloose2,      "dzmuloose2/F");

    if (LEPTONS_2011) {
      ana_tree->Branch("isomuloose1",     &isomuloose1,     "isomuloose1/F");
      ana_tree->Branch("isomuloose2",     &isomuloose2,     "isomuloose2/F");
    }
    if (LEPTONS_2012) {
      ana_tree->Branch("chHadmuloose1",   &chHadmuloose1,   "chHadmuloose1/F");
      ana_tree->Branch("chHadmuloose2",   &chHadmuloose2,   "chHadmuloose2/F");
      ana_tree->Branch("nHadmuloose1",    &nHadmuloose1,    "nHadmuloose1/F");
      ana_tree->Branch("nHadmuloose2",    &nHadmuloose2,    "nHadmuloose2/F");
      ana_tree->Branch("photmuloose1",    &photmuloose1,    "photmuloose1/F");
      ana_tree->Branch("photmuloose2",    &photmuloose2,    "photmuloose2/F");
      ana_tree->Branch("puptmuloose1",    &puptmuloose1,    "puptmuloose1/F");
      ana_tree->Branch("puptmuloose2",    &puptmuloose2,    "puptmuloose2/F");
    }

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
    isemid.setaetaEB=           1000.;
    isemid.setaetaEE=           1000.;
    
    // Lepton tag selection: electrons - Hgg 2011 analysis
    // this is WP85, 2011 cut based: https://twiki.cern.ch/twiki/bin/view/CMS/SimpleCutBasedEleID2011
    electronidcuts eletag2011;
    eletag2011.eta       = 2.5;
    eletag2011.crack1    = 1.4442;
    eletag2011.crack2    = 1.566;
    eletag2011.pt        = 5.;           
    eletag2011.setaetaEB = 0.01;
    eletag2011.setaetaEE = 0.031;
    eletag2011.dphiEB    = 0.039;
    eletag2011.dphiEE    = 0.028;
    eletag2011.detaEB    = 0.005;
    eletag2011.detaEE    = 0.007;
    eletag2011.minhitsEB = 0;
    eletag2011.minhitsEE = 0;
    eletag2011.dcotEB    = 0.02;
    eletag2011.dcotEE    = 0.02;
    eletag2011.distEB    = 0.02;
    eletag2011.distEE    = 0.02;
    eletag2011.d0EB      = 0.02;
    eletag2011.d0EE      = 0.02;
    eletag2011.dzEB      = 0.1;
    eletag2011.dzEE      = 0.1;
    eletag2011.iso_relEB = 0.053;
    eletag2011.iso_relEE = 0.042;

    // Lepton tag selection: electrons - this is WP95, 2011 cut based 
    // https://twiki.cern.ch/twiki/bin/view/CMS/SimpleCutBasedEleID2011
    electronidcuts eletagLoose2011;
    eletagLoose2011.eta       = 2.5;
    eletagLoose2011.crack1    = 1.4442;
    eletagLoose2011.crack2    = 1.566;
    eletagLoose2011.pt        = 5.;     
    eletagLoose2011.setaetaEB = 0.012;
    eletagLoose2011.setaetaEE = 0.031;
    eletagLoose2011.dphiEB    = 0.8;
    eletagLoose2011.dphiEE    = 0.7;
    eletagLoose2011.detaEB    = 0.007;
    eletagLoose2011.detaEE    = 0.011;
    eletagLoose2011.minhitsEB = 0;
    eletagLoose2011.minhitsEE = 0;
    eletagLoose2011.dcotEB    = 0.;
    eletagLoose2011.dcotEE    = 0.;
    eletagLoose2011.distEB    = 0.;
    eletagLoose2011.distEE    = 0.;
    eletagLoose2011.d0EB      = 0.02;
    eletagLoose2011.d0EE      = 0.02;
    eletagLoose2011.dzEB      = 0.1;
    eletagLoose2011.dzEE      = 0.1;
    eletagLoose2011.iso_relEB = 0.15;
    eletagLoose2011.iso_relEE = 0.10;

    // Lepton tag selection: muons
    // 2011 Hgg analysis, close to the tight ID in https://twiki.cern.ch/twiki/bin/view/CMSPublic/SWGuideMuonId
    muonidcuts mutag2011;
    mutag2011.eta     = 2.4;
    mutag2011.pt      = 5.;       
    mutag2011.pixhits = 0;
    mutag2011.tkhits  = 10;   
    mutag2011.hits    = 0;
    mutag2011.chi2    = 10;
    mutag2011.match   = 1;
    mutag2011.d0      = 0.02;
    mutag2011.dz      = 0.1;
    mutag2011.iso_rel = 0.1;

    // Lepton tag selection: muons
    // loose isolation from in https://twiki.cern.ch/twiki/bin/view/CMSPublic/SWGuideMuonId
    // + relaxed number of tracker hits
    muonidcuts mutagLoose2011;
    mutagLoose2011.eta     = 2.4;
    mutagLoose2011.pt      = 5.;       
    mutagLoose2011.pixhits = 0;
    mutagLoose2011.tkhits  = 8;   
    mutagLoose2011.hits    = 0;
    mutagLoose2011.chi2    = 10;
    mutagLoose2011.match   = 1;
    mutagLoose2011.d0      = 0.02;
    mutagLoose2011.dz      = 0.1;
    mutagLoose2011.iso_rel = 0.15;

    // Lepton tag selection 2012: electrons
    // this is the medium WP (80) 
    // in https://twiki.cern.ch/twiki/bin/view/CMS/EgammaCutBasedIdentification
    electronidcuts2012 eletag2012;
    eletag2012.eta       = 2.5;
    eletag2012.crack1    = 1.4442;
    eletag2012.crack2    = 1.566;
    eletag2012.pt        = 5.;
    eletag2012.setaetaEB = 0.01;
    eletag2012.setaetaEE = 0.03;
    eletag2012.dphiEB    = 0.06;   
    eletag2012.dphiEE    = 0.03;   
    eletag2012.detaEB    = 0.004;
    eletag2012.detaEE    = 0.007;
    eletag2012.hoeEB     = 0.12;
    eletag2012.hoeEE     = 0.10;
    eletag2012.oemopEB   = 0.05;
    eletag2012.oemopEE   = 0.05;
    eletag2012.d0EB      = 0.02;
    eletag2012.d0EE      = 0.02;
    eletag2012.dzEB      = 0.1;
    eletag2012.dzEE      = 0.1;
    eletag2012.minhitsEB = 1;
    eletag2012.minhitsEE = 1;
    eletag2012.iso_relEB = 0.15;  
    eletag2012.iso_relEE = 0.15;  

    // Lepton tag selection 2012: electrons
    // this is the loose WP in https://twiki.cern.ch/twiki/bin/view/CMS/EgammaCutBasedIdentification
    electronidcuts2012 eletagLoose2012;
    eletagLoose2012.eta       = 2.5;
    eletagLoose2012.crack1    = 1.4442;
    eletagLoose2012.crack2    = 1.566;
    eletagLoose2012.pt        = 5.;
    eletagLoose2012.setaetaEB = 0.01;
    eletagLoose2012.setaetaEE = 0.03;
    eletagLoose2012.dphiEB    = 0.15;   
    eletagLoose2012.dphiEE    = 0.10;   
    eletagLoose2012.detaEB    = 0.007;
    eletagLoose2012.detaEE    = 0.009;
    eletagLoose2012.hoeEB     = 0.12;
    eletagLoose2012.hoeEE     = 0.10;
    eletagLoose2012.oemopEB   = 0.05;
    eletagLoose2012.oemopEE   = 0.05;
    eletagLoose2012.d0EB      = 0.02;
    eletagLoose2012.d0EE      = 0.02;
    eletagLoose2012.dzEB      = 0.2;
    eletagLoose2012.dzEE      = 0.2;
    eletagLoose2012.minhitsEB = 1;
    eletagLoose2012.minhitsEE = 1;
    eletagLoose2012.iso_relEB = 0.15;  
    eletagLoose2012.iso_relEE = 0.15;  

    // Lepton tag selection 2012: tight muon selection as defined in 
    // https://twiki.cern.ch/twiki/bin/view/CMSPublic/SWGuideMuonId
    // and loose isolation cut from https://twiki.cern.ch/twiki/bin/view/CMSPublic/SWGuideMuonId#Muon_Isolation
    muonidcuts2012 mutagLoose2012;
    mutagLoose2012.eta     = 2.4;
    mutagLoose2012.pt      = 5.;
    mutagLoose2012.chi2    = 10;
    mutagLoose2012.hits    = 0;
    mutagLoose2012.match   = 1;
    mutagLoose2012.pixhits = 0;
    mutagLoose2012.withm   = 5;
    mutagLoose2012.d0      = 0.2;   
    mutagLoose2012.dz      = 0.5;   
    mutagLoose2012.iso_rel = 0.2;   

    // Lepton tag selection 2012: tight muon selection as defined in 
    // https://twiki.cern.ch/twiki/bin/view/CMSPublic/SWGuideMuonId 
    // + tight isolation cut from https://twiki.cern.ch/twiki/bin/view/CMSPublic/SWGuideMuonId#Muon_Isolation
    muonidcuts2012 mutag2012;
    mutag2012.eta     = 2.4;
    mutag2012.pt      = 5.;
    mutag2012.chi2    = 10;
    mutag2012.hits    = 0;
    mutag2012.match   = 1;
    mutag2012.pixhits = 0;
    mutag2012.withm   = 5;
    mutag2012.d0      = 0.2;   
    mutag2012.dz      = 0.5;   
    mutag2012.iso_rel = 0.12;   


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
        float cic4pf_cuts_lead[phoNCUTS][phoCiC4NCATEGORIES];
        float cic4pf_cuts_sublead[phoNCUTS][phoCiC4NCATEGORIES];
        SetPhotonCutsInCategories((phoCiCIDLevel)iLevel, &cic6_cuts_lead[0][0], &cic6_cuts_sublead[0][0], &cic4_cuts_lead[0][0], &cic4_cuts_sublead[0][0] ,
				  	            &cic4pf_cuts_lead[0][0], &cic4pf_cuts_sublead[0][0]);
     
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
     
	float * cic4pf_cuts_arrays_lead[phoNCUTS] = {
            &cic4pf_cut_lead_isosumoet[0][0], 
            &cic4pf_cut_lead_isosumoetbad[0][0], 
            &cic4pf_cut_lead_trkisooet[0][0], 
            &cic4pf_cut_lead_sieie[0][0],
            &cic4pf_cut_lead_hovere[0][0], 
            &cic4pf_cut_lead_r9[0][0], 
            &cic4pf_cut_lead_drtotk_25_99[0][0], 
            &cic4pf_cut_lead_pixel[0][0] 
        };
        
        float * cic4pf_cuts_arrays_sublead[phoNCUTS] = {
	  &cic4pf_cut_sublead_isosumoet[0][0], 
	  &cic4pf_cut_sublead_isosumoetbad[0][0], 
	  &cic4pf_cut_sublead_trkisooet[0][0], 
	  &cic4pf_cut_sublead_sieie[0][0], 
	  &cic4pf_cut_sublead_hovere[0][0], 
	  &cic4pf_cut_sublead_r9[0][0],
	  &cic4pf_cut_sublead_drtotk_25_99[0][0], 
	  &cic4pf_cut_sublead_pixel[0][0]
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
	  for(int iCat=0; iCat<phoCiC4NCATEGORIES; ++iCat) {
	    cic4pf_cuts_arrays_lead[iCut][iLevel*phoCiC4NCATEGORIES+iCat] = cic4pf_cuts_lead[iCut][iCat];
	    cic4pf_cuts_arrays_sublead[iCut][iLevel*phoCiC4NCATEGORIES+iCat] = cic4pf_cuts_sublead[iCut][iCat];
	  }
        }
    } // end of loop over all photon cut levels

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


        /// bug fix:
        /// when  nPreselPhotonPairs==0 the vrank variables are not initialized
        if(nPreselPhotonPairs==0)
        {
            indexPreselPhot1[0] = 0;   //[nPreselPhotonPairs]
            indexPreselPhot2[0] = 0;   //[nPreselPhotonPairs]
            vrankPhotonPairs[0] = 0;   //[nPreselPhotonPairs]
            vevtMvaPhotonPairs[0] = 0;   //[nPreselPhotonPairs]
            vevtProbPhotonPairs[0] = 0;   //[nPreselPhotonPairs]
            vptbalPhotonPairs[0] = 0;   //[nPreselPhotonPairs]
            vptasymPhotonPairs[0] = 0;   //[nPreselPhotonPairs]
        }

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
        
        vector<int> firstfourgenphot = firstones(ptMC,&photassocMC,4);
        vector<int> firstfourhiggsgenphot = firstones(ptMC,&photassocMChiggs,4);
        
	// gen level info for leptons  
	vector<int> genVHLepton;
	if (gen_custom_processId<3400 && gen_custom_processId>3100) {         // W leptonic
	  genVHLepton.push_back(VHLeptonIndexes[0]);
	  genVHLepton.push_back(-999);
	} else if (gen_custom_processId<4400 && gen_custom_processId>4100) {  // Z leptonic 
	  genVHLepton.push_back(VHLeptonIndexes[0]);
	  genVHLepton.push_back(VHLeptonIndexes[1]);
	} else {
	  genVHLepton.push_back(-999);
	  genVHLepton.push_back(-999);
	}

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
	if (gen_custom_processId==3101 || gen_custom_processId==3201 || gen_custom_processId==3301) ptphotgen1wl.Fill(ptMC[firstfourgenphot.at(0)],weight);
	if (gen_custom_processId==4101 || gen_custom_processId==4201 || gen_custom_processId==4301) ptphotgen1zl.Fill(ptMC[firstfourgenphot.at(0)],weight);
	if (gen_custom_processId==3401) ptphotgen1wh.Fill(ptMC[firstfourgenphot.at(0)],weight);
	if (gen_custom_processId==4401) ptphotgen1zh.Fill(ptMC[firstfourgenphot.at(0)],weight);
	if (gen_custom_processId==4501) ptphotgen1zn.Fill(ptMC[firstfourgenphot.at(0)],weight);
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

        // define what event it is:
        H_event=false;
        V_event=false;
        Vqq_event=false;
        Zbb_event=false;

        for(Int_t iPartMC=0; iPartMC<nMC; ++iPartMC) {
        
          if( statusMC[iPartMC]!=3 ) continue;
        
          if( pdgIdMC[iPartMC]==25 ) H_event = true; 
          if( pdgIdMC[iPartMC]==23 || abs(pdgIdMC[iPartMC])==24 ) V_event = true; 
          if( abs(pdgIdMC[iPartMC])==5 && pdgIdMC[motherIDMC[iPartMC]]==23 ) Zbb_event = true; 
          if( abs(pdgIdMC[iPartMC])<=5 && (pdgIdMC[motherIDMC[iPartMC]]==23 || abs(pdgIdMC[motherIDMC[iPartMC]])==24) ) Vqq_event = true; 
        
        } //for MC particles


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

         // TLorentzVector genJet1, genJet2;
         // genJet1.SetPtEtaPhiE( ptJetGen_akt5[isoJetIndex[0]],  etaJetGen_akt5[isoJetIndex[0]],  phiJetGen_akt5[isoJetIndex[0]],  eJetGen_akt5[isoJetIndex[0]]);
         // 
         // genJet2.SetPtEtaPhiE( ptJetGen_akt5[isoJetIndex[1]], etaJetGen_akt5[isoJetIndex[1]], phiJetGen_akt5[isoJetIndex[1]], eJetGen_akt5[isoJetIndex[1]]);
         // TLorentzVector dijet = genJet1 + genJet2;

#ifdef DEBUG
            cout << "[DEBUG] after genJet1/2 " << endl;
#endif
        
         // double aveeta  = (etaJetGen_akt5[isoJetIndex[0]] + etaJetGen_akt5[isoJetIndex[1]])/2.;
         // gen_zeppenfeld = diphot.Eta() - aveeta;

#ifdef DEBUG
            cout << "[DEBUG] after diphot.Eta() for zeppen computation " << endl;
#endif
            
            gen_pt_gamma1  = ptMC[index_phot1];
            gen_pt_gamma2  = ptMC[index_phot2];
            gen_eta_gamma1 = etaMC[index_phot1];
            gen_eta_gamma2 = etaMC[index_phot2] ;
            gen_phi_gamma1 = phiMC[index_phot1];
            gen_phi_gamma2 = phiMC[index_phot2] ;
            
         // gen_pt_genjet1  =  ptJetGen_akt5[isoJetIndex[0]];
         // gen_pt_genjet2  =  ptJetGen_akt5[isoJetIndex[1]];
         // gen_eta_genjet1 =  etaJetGen_akt5[isoJetIndex[0]];
         // gen_eta_genjet2 =  etaJetGen_akt5[isoJetIndex[1]];
         // gen_phi_genjet1 =  phiJetGen_akt5[isoJetIndex[0]];
         // gen_phi_genjet2 =  phiJetGen_akt5[isoJetIndex[1]];
            
            gen_mass_diphoton  = diphot.M();
            gen_pt_diphoton    = diphot.Pt();
            gen_eta_diphoton   = diphot.Eta();
            gen_phi_diphoton   = diphot.Phi();
                
         // gen_mass_dijet = dijet.M();
         // gen_pt_dijet   = dijet.Pt();
         // gen_eta_dijet  = dijet.Eta();
         // gen_phi_dijet  = dijet.Phi();
        }

	// gen. level variables for lepton tag
	int index_lep1 = genVHLepton.at(0);
	int index_lep2 = genVHLepton.at(1);
	if (index_lep1>-1 && index_lep2>-1) {
	  if(ptMC[index_lep1]>ptMC[index_lep2]) {
	    gen_pt_lep1  = ptMC[index_lep1];
	    gen_pt_lep2  = ptMC[index_lep2];
	    gen_eta_lep1 = etaMC[index_lep1];
	    gen_eta_lep2 = etaMC[index_lep2];
	    gen_phi_lep1 = phiMC[index_lep1];
	    gen_phi_lep2 = phiMC[index_lep2];
	    gen_pid_lep1 = pdgIdMC[index_lep1];
	    gen_pid_lep2 = pdgIdMC[index_lep2];
	  } else {
	    gen_pt_lep1  = ptMC[index_lep2];
	    gen_pt_lep2  = ptMC[index_lep1];
	    gen_eta_lep1 = etaMC[index_lep2];
	    gen_eta_lep2 = etaMC[index_lep1];
	    gen_phi_lep1 = phiMC[index_lep2];
	    gen_phi_lep2 = phiMC[index_lep1];
	    gen_pid_lep1 = pdgIdMC[index_lep2];
	    gen_pid_lep2 = pdgIdMC[index_lep1];
	  } 

	} else if (index_lep1>-1 && index_lep2<0) {
	  gen_pt_lep1  = ptMC[index_lep1];
	  gen_eta_lep1 = etaMC[index_lep1];
	  gen_phi_lep1 = phiMC[index_lep1];
	  gen_pid_lep1 = pdgIdMC[index_lep1];
	  gen_pt_lep2  = -500.;
	  gen_eta_lep2 = -500.;
	  gen_phi_lep2 = -500.;
	  gen_pid_lep2 = -500;	  

	} else if (index_lep2>-1 && index_lep1<0) {
	  gen_pt_lep1  = ptMC[index_lep2];
	  gen_eta_lep1 = etaMC[index_lep2];
	  gen_phi_lep1 = phiMC[index_lep2];
	  gen_pid_lep1 = pdgIdMC[index_lep2];
	  gen_pt_lep2  = -500.;
	  gen_eta_lep2 = -500.;
	  gen_phi_lep2 = -500.;
	  gen_pid_lep2 = -500;
	  
	} else if (index_lep2<0 && index_lep1<0) {
	  gen_pt_lep1  = -500.;
	  gen_eta_lep1 = -500.;
	  gen_phi_lep1 = -500.;
	  gen_pid_lep1 = -500;
	  gen_pt_lep2  = -500.;
	  gen_eta_lep2 = -500.;
	  gen_phi_lep2 = -500.;
	  gen_pid_lep2 = -500;
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
        vector<float>  isomva;
        vector<int>  isocic;
        vector<int>  isocicnoelveto;
        vector<int>  isocicpf;
        vector<int>  isocicpfnoelveto;
        
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
            
	    vtxId=vrankPhotonPairs[0];
	    vtxPos_x=vx[vrankPhotonPairs[0]];
	    vtxPos_y=vy[vrankPhotonPairs[0]];
	    vtxPos_z=vz[vrankPhotonPairs[0]];
	    vtxIdMVA=vevtMvaPhotonPairs[0];
	    vtxIdEvtProb=vevtProbPhotonPairs[0];

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
            else if (finder == "preselectionMVA") preselection = PhotonMITPreSelection(i,vrankPhotonPairs[0],1);	
            else if (finder == "preselectionMVAnoeleveto") preselection = PhotonMITPreSelection(i,vrankPhotonPairs[0],0);	
            else if (finder == "cicloose") preselection = PhotonCiCSelectionLevel(i,1,vrankPhotonPairs[0],0) >= 1;	
            else if (finder == "cicloosenoeleveto") preselection = PhotonCiCSelectionLevel(i,0,vrankPhotonPairs[0],0) >= 1;	
            else if (finder == "cicmedium") preselection = PhotonCiCSelectionLevel(i,1,vrankPhotonPairs[0],0) >= 2;	
            else if (finder == "cictight") preselection = PhotonCiCSelectionLevel(i,1,vrankPhotonPairs[0],0) >= 3;	
            else if (finder == "cicsuper") preselection = PhotonCiCSelectionLevel(i,1,vrankPhotonPairs[0],0) >= 4;	
            else if (finder == "cichyper") preselection = PhotonCiCSelectionLevel(i,1,vrankPhotonPairs[0],0) >= 5;	
            else if (finder == "cicpfloose") preselection = PhotonCiCSelectionLevel(i,1,vrankPhotonPairs[0],1) >= 1;	
            else if (finder == "cicpfloosenoeleveto") preselection = PhotonCiCSelectionLevel(i,0,vrankPhotonPairs[0],1) >= 1;	
            else if (finder == "cicpfmedium") preselection = PhotonCiCSelectionLevel(i,1,vrankPhotonPairs[0],1) >= 2;	
            else if (finder == "cicpftight") preselection = PhotonCiCSelectionLevel(i,1,vrankPhotonPairs[0],1) >= 3;	
            else if (finder == "cicpfsuper") preselection = PhotonCiCSelectionLevel(i,1,vrankPhotonPairs[0],1) >= 4;	
            else if (finder == "cicpfhyper") preselection = PhotonCiCSelectionLevel(i,1,vrankPhotonPairs[0],1) >= 5;	
            else if (finder == "mcass") preselection = mcID(i);
            else {
              cout << "NO SUCH " << selection << " PRESELECTION  AVAILABLE!!" << endl;
                  cout << "Good options are: superloose loose medium isem looseeg tighteg hggtighteg looseegpu tightegpu hggtightegpu preselection preselectionCS preselectionMVA preselectionMVAnoeleveto cicloose cicloosenoeleveto cicmedium cictight cicsuper cichyper mcass cicpfloose cicpfloosenoeleveto cicpfmedium cicpftight cicpfsuper cicpfhyper mcass" << endl;
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
            
            isomva.push_back(PhotonIDMVANew(i,vrankPhotonPairs[0]));
            isocic.push_back(PhotonCiCSelectionLevel(i,1,vrankPhotonPairs[0],0));
            isocicnoelveto.push_back(PhotonCiCSelectionLevel(i,0,vrankPhotonPairs[0],0));
            isocicpf.push_back(PhotonCiCSelectionLevel(i,1,vrankPhotonPairs[0],1));
            isocicpfnoelveto.push_back(PhotonCiCSelectionLevel(i,0,vrankPhotonPairs[0],1));

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
       
        vector<int> firstfourhiggsassphot = firstones(ptPhot,&assophothiggs,4);
        vector<int> firstfourassphot = firstones(ptPhot,&assophot,4);
        vector<int> firstfourisophot = firstones(ptPhot,&isophot,4);      
        
        vector<bool> jetnohiggsphot;
        vector<bool> jetnoassphot;
	vector<bool> jetgoodnoisophot;
	jetnoisophot.clear();

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

	    if(TMath::Abs(etaJet_pfakt5[i]) < 2.5) {
	      if(betaStar_pfakt5[i][vrankPhotonPairs[0]] > 0.2 * log( nvertex - 0.67 ) ) goodetajet = 0;
	      if(rmsCandJet_pfakt5[i] > 0.06) goodetajet = 0;
// 	    } else if(TMath::Abs(etaJet_pfakt5[i]) < 2.75){
// 	      if(betaStar_pfakt5[i][vrankPhotonPairs[0]] > 0.3 * log( nvertex - 0.67 ) ) goodetajet = 0;
//	      if(rmsCandJet_pfakt5[i] > 0.06) goodetajet = 0;	     
	    } else if(TMath::Abs(etaJet_pfakt5[i]) < 3){
 	      if(rmsCandJet_pfakt5[i] > 0.05) goodetajet = 0;
 	    } else {
	      if(rmsCandJet_pfakt5[i] > 0.055) goodetajet = 0;
 	    }
		            
            if(!assh && goodetajet) jetnohiggsphot.push_back(1);
            else jetnohiggsphot.push_back(0); 
            
            if(!assp && goodetajet) jetnoassphot.push_back(1); 
            else jetnoassphot.push_back(0);  
            
            if(!assi && goodetajet) jetgoodnoisophot.push_back(1); 
            else jetgoodnoisophot.push_back(0);  
 
	    if(!assi) jetnoisophot.push_back(1); 
            else jetnoisophot.push_back(0);  
          
	    int ass_here(-999);
	    double DRmin_here(999.);
	    for(int j=0; j<nJetGen_akt5; j++){
	      double DR = sqrt(delta_eta(etaJet_pfakt5[i],etaJetGen_akt5[j])*delta_eta(etaJet_pfakt5[i],etaJetGen_akt5[j]) +
			       delta_phi(phiJet_pfakt5[i],phiJetGen_akt5[j])*delta_phi(phiJet_pfakt5[i],phiJetGen_akt5[j]) ) ;
	      double expres = ErrEt(ptCorrJet_pfakt5[i],etaJet_pfakt5[i]);
	      //      if(DR < DRmin && (ptCorrJet_pfakt5[i]-ptJetGen_akt5[j])/ptCorrJet_pfakt5[i] < 5. * expres) {
	      if(DR < DRmin_here && TMath::Abs(ptCorrJet_pfakt5[i]-ptJetGen_akt5[j])/ptJetGen_akt5[j] < .75) {
		ass_here = j;
		DRmin_here = DR;
	      }
	    }
	    
	    if(DRmin_here > 0.1 )  ass_here = -999;
	    
	    if(!assi && ass_here>-1) {
	      jetDR->Fill(ptJetGen_akt5[ass_here],DRmin_here);
	      jetresp_vs_pt->Fill(ptJetGen_akt5[ass_here],(ptCorrJet_pfakt5[i]-ptJetGen_akt5[ass_here])/ptJetGen_akt5[ass_here]);
	      if(ptJetGen_akt5[ass_here]>20 && ptJetGen_akt5[ass_here]<50) {
		jetresp_vs_eta->Fill(etaJetGen_akt5[ass_here],(ptCorrJet_pfakt5[i]-ptJetGen_akt5[ass_here])/ptJetGen_akt5[ass_here]);
		jetresp_vs_npu->Fill(npu,(ptCorrJet_pfakt5[i]-ptJetGen_akt5[ass_here])/ptJetGen_akt5[ass_here]);
	      }
	      if(ptJetGen_akt5[ass_here]>50 && ptJetGen_akt5[ass_here]<150) {
		jetresp_vs_eta_50->Fill(etaJetGen_akt5[ass_here],(ptCorrJet_pfakt5[i]-ptJetGen_akt5[ass_here])/ptJetGen_akt5[ass_here]);
		jetresp_vs_eta_50_abs->Fill(TMath::Abs(etaJetGen_akt5[ass_here]),(ptCorrJet_pfakt5[i]-ptJetGen_akt5[ass_here])/ptJetGen_akt5[ass_here]);
		jetresp_vs_npu_50->Fill(npu,(ptCorrJet_pfakt5[i]-ptJetGen_akt5[ass_here])/ptJetGen_akt5[ass_here]);
	      }
	      if(ptJetGen_akt5[ass_here]>150) {
		jetresp_vs_eta_150->Fill(etaJetGen_akt5[ass_here],(ptCorrJet_pfakt5[i]-ptJetGen_akt5[ass_here])/ptJetGen_akt5[ass_here]);
		jetresp_vs_npu_150->Fill(npu,(ptCorrJet_pfakt5[i]-ptJetGen_akt5[ass_here])/ptJetGen_akt5[ass_here]);
	      }
	      if(ptJetGen_akt5[ass_here]>20 && TMath::Abs(etaJet_akt5[i])>3.) {
		jetresp_vs_npu_forward->Fill(npu,(ptCorrJet_pfakt5[i]-ptJetGen_akt5[ass_here])/ptJetGen_akt5[ass_here]);
		jetresp_vs_pt_forward->Fill(ptJetGen_akt5[ass_here],(ptCorrJet_pfakt5[i]-ptJetGen_akt5[ass_here])/ptJetGen_akt5[ass_here]);		
	      }
	      
	    }

	}

        vector<int> firsttennohiggsjet = firstones(ptCorrJet_pfakt5,&jetnohiggsphot,10);
        vector<int> firsttennoassjet = firstones(ptCorrJet_pfakt5,&jetnoassphot,10);
        vector<int> firsttennoisojet = firstones(ptCorrJet_pfakt5,&jetgoodnoisophot,10);      
        
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


        // jets
        for( unsigned ijet=0; ijet<10; ++ijet ) {
          ptjet[ijet] = -999.;
          ptcorrjet[ijet] = -999.;
          ecorrjet[ijet] = -999.;
          etajet[ijet] = -999.;
          phijet[ijet] = -999.;
        }

        njets = 0;
        
        for( unsigned ijet=0; ijet<firsttennoisojet.size(); ++ijet ) {

          if( firsttennoisojet.at(ijet)>=0 ) {
           
            if( ptCorrJet_pfakt5[firsttennoisojet.at(ijet)] < 20. ) continue;
            if( njets >=10 ) continue;

            ptjet[njets] = ptJet_pfakt5[firsttennoisojet.at(ijet)];
            ptcorrjet[njets] = ptCorrJet_pfakt5[firsttennoisojet.at(ijet)];	  
            ecorrjet[njets] = ptcorrjet[njets]/ptjet[njets]*eJet_pfakt5[firsttennoisojet.at(ijet)];
            etajet[njets] = etaJet_pfakt5[firsttennoisojet.at(ijet)];
            phijet[njets] = phiJet_pfakt5[firsttennoisojet.at(ijet)];
            betajet[njets] = beta_pfakt5[firsttennoisojet.at(ijet)][vrankPhotonPairs[0]];
            betastarjet[njets] = betaStar_pfakt5[firsttennoisojet.at(ijet)][vrankPhotonPairs[0]];
            assjet[njets] = assoJet(firsttennoisojet.at(ijet));
            btagvtxjet[njets] = simpleSecondaryVertexHighEffBJetTags[firsttennoisojet.at(ijet)];
            btagtrkjet[njets] = trackCountingHighEffBJetTags[firsttennoisojet.at(ijet)];	  
            btagjprobjet[njets] = jetProbabilityBJetTags[firsttennoisojet.at(ijet)];	  
            ptDjet[njets] = ptDJet_pfakt5[firsttennoisojet.at(ijet)];
            rmsjet[njets] = rmsCandJet_pfakt5[firsttennoisojet.at(ijet)];
            ntrkjet[njets] = nChargedHadrons_pfakt5[firsttennoisojet.at(ijet)];
            nneutjet[njets] = nPhotons_pfakt5[firsttennoisojet.at(ijet)] + nNeutralHadrons_pfakt5[firsttennoisojet.at(ijet)] + nHFHadrons_pfakt5[firsttennoisojet.at(ijet)] + nHFEM_pfakt5[firsttennoisojet.at(ijet)];
            jetIdSimple_mvajet[njets] = jetIdSimple_mva_pfakt5[firsttennoisojet.at(ijet)];
            jetIdFull_mvajet[njets] = jetIdFull_mva_pfakt5[firsttennoisojet.at(ijet)];
            jetId_dR2Meanjet[njets] = jetId_dR2Mean_pfakt5[firsttennoisojet.at(ijet)];
            jetId_betaStarClassicjet[njets] = jetId_betaStarClassic_pfakt5[firsttennoisojet.at(ijet)];
            jetIdCutBased_wpjet[njets] = jetIdCutBased_wp_pfakt5[firsttennoisojet.at(ijet)];
            jetIdSimple_wpjet[njets] = jetIdSimple_wp_pfakt5[firsttennoisojet.at(ijet)];	  
            jetIdFull_wpjet[njets] = jetIdFull_wp_pfakt5[firsttennoisojet.at(ijet)];	  
            jetId_frac01jet[njets] = jetId_frac01_pfakt5[firsttennoisojet.at(ijet)];
            jetId_frac02jet[njets] = jetId_frac02_pfakt5[firsttennoisojet.at(ijet)];
            jetId_frac03jet[njets] = jetId_frac03_pfakt5[firsttennoisojet.at(ijet)];
            jetId_frac04jet[njets] = jetId_frac04_pfakt5[firsttennoisojet.at(ijet)];
            jetId_frac05jet[njets] = jetId_frac05_pfakt5[firsttennoisojet.at(ijet)];
            jetId_betajet[njets] = jetId_beta_pfakt5[firsttennoisojet.at(ijet)];
            jetId_betaStarjet[njets] = jetId_betaStar_pfakt5[firsttennoisojet.at(ijet)];

            // match to parton
            Float_t deltaRMCmin = 999.;
            Int_t pdgIdPart_found = 0;
            Int_t pdgIdMomPart_found = 0;
         
            for(Int_t iPartMC=0; iPartMC<nMC; ++iPartMC) {
         
              if( statusMC[iPartMC]!=3 ) continue;
         
              if( ptMC[iPartMC]<0.1 ) continue;
         
              TLorentzVector jet;
              jet.SetPtEtaPhiE( ptcorrjet[njets], etajet[njets], phijet[njets], ecorrjet[njets] );
              TLorentzVector parton;
              parton.SetPtEtaPhiE( ptMC[iPartMC], etaMC[iPartMC], phiMC[iPartMC], eMC[iPartMC] );
         
              Int_t pdgId = pdgIdMC[iPartMC];
            
              Float_t deltaRMC = jet.DeltaR(parton);
         
              bool goodPdgId = ( (fabs(pdgId)<=9) || (fabs(pdgId)==21) );
              if( !goodPdgId ) continue;
            
              if( (deltaRMC < deltaRMCmin) && goodPdgId ) {
                deltaRMCmin = deltaRMC;
                pdgIdPart_found = pdgIdMC[iPartMC];
                pdgIdMomPart_found = pdgIdMC[motherIDMC[iPartMC]];
              }
         
            } //for MC particles
        
        
            partPdgIDjet[njets] = ( deltaRMCmin<0.5 ) ? pdgIdPart_found : -999; 
            partMomPdgIDjet[njets] = ( deltaRMCmin<0.5 ) ? pdgIdMomPart_found : -999; 

            h1_deltaR_jetpart->Fill( deltaRMCmin );

            njets++;

          } // if good index 

        } // for firsttennoisojet



        if( firstfourhiggsassphot.at(0)>-1 && firstfourhiggsassphot.at(1)>-1 ) 
        { 
            if( firsttennohiggsjet.at(0) > -1) {
              ptjethiggsassreco1.Fill(ptCorrJet_pfakt5[firsttennohiggsjet.at(0)],weight);
              etajethiggsassreco1.Fill(etaJet_pfakt5[firsttennohiggsjet.at(0)],weight);
            }
            if( firsttennohiggsjet.at(1) > -1) {
              ptjethiggsassreco2.Fill(ptCorrJet_pfakt5[firsttennohiggsjet.at(1)],weight);
              etajethiggsassreco2.Fill(etaJet_pfakt5[firsttennohiggsjet.at(1)],weight);
            }
            if( firsttennohiggsjet.at(0) > -1 && firsttennohiggsjet.at(1) > -1) {
              deltaetajethiggsassreco.Fill(etaJet_pfakt5[firsttennohiggsjet.at(0)]-etaJet_pfakt5[firsttennohiggsjet.at(1)]);
              double aveeta = (etaJet_pfakt5[firsttennohiggsjet.at(0)]+etaJet_pfakt5[firsttennohiggsjet.at(1)])/2;
              double zeppen1 = etaJet_pfakt5[firsttennohiggsjet.at(0)] - aveeta;
              double zeppen2 = etaJet_pfakt5[firsttennohiggsjet.at(1)] - aveeta;
              zeppenjethiggsassreco1.Fill(zeppen1,weight);
              zeppenjethiggsassreco2.Fill(zeppen2,weight);	
            }
        }
   
        double twojetsmass(0), etatwojets(-999), phitwojets(-999), pttwojets(-999);

        if( firstfourassphot.at(0)>-1 && firstfourassphot.at(1)>-1 ) { 

            if( firsttennoassjet.at(0) > -1) {
              ptjetassreco1.Fill(ptCorrJet_pfakt5[firsttennoassjet.at(0)],weight);
              etajetassreco1.Fill(etaJet_pfakt5[firsttennoassjet.at(0)],weight);
            }
            if( firsttennoassjet.at(1) > -1) {
              ptjetassreco2.Fill(ptCorrJet_pfakt5[firsttennoassjet.at(1)],weight);
              etajetassreco2.Fill(etaJet_pfakt5[firsttennoassjet.at(1)]),weight;
            }
            if( firsttennoassjet.at(0) > -1 && firsttennoassjet.at(1) > -1) {
              TLorentzVector jet1, jet2;	
              jet1.SetPtEtaPhiE(ptCorrJet_pfakt5[firsttennoassjet.at(0)],etaJet_pfakt5[firsttennoassjet.at(0)],phiJet_pfakt5[firsttennoassjet.at(0)],eJet_pfakt5[firsttennoassjet.at(0)]/ptJet_pfakt5[firsttennoassjet.at(0)]*ptCorrJet_pfakt5[firsttennoassjet.at(0)]);
              jet2.SetPtEtaPhiE(ptCorrJet_pfakt5[firsttennoassjet.at(1)],etaJet_pfakt5[firsttennoassjet.at(1)],phiJet_pfakt5[firsttennoassjet.at(1)],eJet_pfakt5[firsttennoassjet.at(1)]/ptJet_pfakt5[firsttennoassjet.at(1)]*ptCorrJet_pfakt5[firsttennoassjet.at(1)]);
              
              TLorentzVector sum = jet1 + jet2;
              
              twojetsmass = sum.M();
              etatwojets = sum.Eta();
	      phitwojets = sum.Phi();
	      pttwojets = sum.Pt();
             
              deltaetajetassreco.Fill(etaJet_pfakt5[firsttennoassjet.at(0)]-etaJet_pfakt5[firsttennoassjet.at(1)],weight);
              double aveeta = (etaJet_pfakt5[firsttennoassjet.at(0)]+etaJet_pfakt5[firsttennoassjet.at(1)])/2;
              double zeppen1 = etaJet_pfakt5[firsttennoassjet.at(0)] - aveeta;
              double zeppen2 = etaJet_pfakt5[firsttennoassjet.at(1)] - aveeta;
              zeppenjetassreco1.Fill(zeppen1,weight);
              zeppenjetassreco2.Fill(zeppen2,weight);	
            }

            if(  ptCorrJet_pfakt5[firsttennoassjet.at(0)] > ptjet1cut && ptCorrJet_pfakt5[firsttennoassjet.at(1)] > ptjet2cut 
                && ptPhot[firstfourassphot.at(0)] > ptphot1cut && ptPhot[firstfourassphot.at(1)] > ptphot2cut )
            {
                if(TMath::Abs(etaJet_pfakt5[firsttennoassjet.at(0)]-etaJet_pfakt5[firsttennoassjet.at(1)])>deltaetacut)
                {
                    higgsmasscutreco.Fill(higgsmass,weight);
                    if(isophot.at(firstfourassphot.at(0)) && isophot.at(firstfourassphot.at(1)))  
                        higgsmassisorecocheck.Fill(higgsisomass,weight);	    

                    double aveeta = (etaJet_pfakt5[firsttennoassjet.at(0)]+etaJet_pfakt5[firsttennoassjet.at(1)])/2;
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
    
	/****************************************************
	 *                                                  *
	 *        LEPTON TAG                                *
	 *                                                  *
	 ****************************************************/

	// electron tag 
	vector<bool> idpasseletag2011(11);
	vector<bool> idpasseletag2012(13);
	vector<bool> idpasseletagLoose2011(11);
	vector<bool> idpasseletagLoose2012(13);

	// tight selection
	int firstEle       = -999;
	int secondEle      = -999;
	double firstElePt  = -998.;
	double secondElePt = -999.;

        for(int iEle=0; iEle<nEle; iEle++){
	  
	  if (LEPTONS_2011 && !leptonCutsEle2011(iEle, eletag2011, &idpasseletag2011)) continue; 
	  if (LEPTONS_2012 && !leptonCutsEle2012(iEle, eletag2012, &idpasseletag2012)) continue; 
	  
	  if (electron_pt[iEle]>=secondElePt && electron_pt[iEle]<firstElePt) {
	    secondEle=iEle;
	    secondElePt=electron_pt[iEle];
	  } else if (electron_pt[iEle]>=firstElePt && electron_pt[iEle]>=secondElePt) {
	    secondEle=firstEle;
	    secondElePt=firstElePt;
	    firstEle=iEle;
	    firstElePt=electron_pt[iEle];
	  }
	}

	// loose selection
	int firstEleLoose       = -999;
	int secondEleLoose      = -999;
	double firstEleLoosePt  = -998.;
	double secondEleLoosePt = -999.;

        for(int iEle=0; iEle<nEle; iEle++){
	  
	  if (LEPTONS_2011 && !leptonCutsEle2011(iEle, eletagLoose2011, &idpasseletagLoose2011)) continue; 
	  if (LEPTONS_2012 && !leptonCutsEle2012(iEle, eletagLoose2012, &idpasseletagLoose2012)) continue;   
	  
	  if (electron_pt[iEle]>=secondEleLoosePt && electron_pt[iEle]<firstEleLoosePt) {
	    secondEleLoose=iEle;
	    secondEleLoosePt=electron_pt[iEle];
	  } else if (electron_pt[iEle]>=firstEleLoosePt && electron_pt[iEle]>=secondEleLoosePt) {
	    secondEleLoose=firstEleLoose;
	    secondEleLoosePt=firstEleLoosePt;
	    firstEleLoose=iEle;
	    firstEleLoosePt=electron_pt[iEle];
	  }
	}

	// muon tag 
	vector<bool> idpassmutag2011(11);
	vector<bool> idpassmutag2012(12);
	vector<bool> idpassmutagLoose2011(11);
	vector<bool> idpassmutagLoose2012(12);

	// tight selection
	int firstMu       = -999;
	int secondMu      = -999;
	double firstMuPt  = -998.;
	double secondMuPt = -999.;
        for(int iMu=0; iMu<nMuons; iMu++){

	  if (LEPTONS_2011 && !leptonCutsMu2011(iMu, mutag2011, &idpassmutag2011)) continue; 
	  if (LEPTONS_2012 && !leptonCutsMu2012(iMu, mutag2012, &idpassmutag2012)) continue; 

	  if (Muon_pt[iMu]>=secondMuPt && Muon_pt[iMu]<firstMuPt) {
	    secondMu=iMu;
	    secondMuPt=Muon_pt[iMu];
	  } else if (Muon_pt[iMu]>=firstMuPt && Muon_pt[iMu]>=secondMuPt) {
	    secondMu=firstMu;
	    secondMuPt=firstMuPt;
	    firstMu=iMu;
	    firstMuPt=Muon_pt[iMu];
	  }
	}
	
	// loose selection
	int firstMuLoose       = -999;
	int secondMuLoose      = -999;
	double firstMuLoosePt  = -998.;
	double secondMuLoosePt = -999.;
        for(int iMu=0; iMu<nMuons; iMu++){

	  if (LEPTONS_2011 && !leptonCutsMu2011(iMu, mutagLoose2011, &idpassmutagLoose2011)) continue; 
	  if (LEPTONS_2012 && !leptonCutsMu2012(iMu, mutagLoose2012, &idpassmutagLoose2012)) continue;  

	  if (Muon_pt[iMu]>=secondMuLoosePt && Muon_pt[iMu]<firstMuLoosePt) {
	    secondMuLoose=iMu;
	    secondMuLoosePt=Muon_pt[iMu];
	  } else if (Muon_pt[iMu]>=firstMuLoosePt && Muon_pt[iMu]>=secondMuLoosePt) {
	    secondMuLoose=firstMuLoose;
	    secondMuLoosePt=firstMuLoosePt;
	    firstMuLoose=iMu;
	    firstMuLoosePt=Muon_pt[iMu];
	  }
	}
	
	// filling variables for the tree - lepton tag

	// tight electron selection
	if (firstEle>=0) {
	  ptele1    = electron_pt[firstEle];
	  etaele1   = electron_sc_eta[firstEle];
	  phiele1   = electron_phi[firstEle];
	  eneele1   = electron_energy[firstEle];
	  sIeIeele1 = electron_SigmaIetaIeta[firstEle];
	  dphiele1  = electron_dPhiIn[firstEle];
	  detaele1  = electron_dEtaIn[firstEle];
	  mhitsele1 = electron_misHits[firstEle];
	  d0ele1    = eleDxyPV(firstEle,vrankPhotonPairs[0]); 
	  dzele1    = eleDzPV(firstEle,vrankPhotonPairs[0]); 

	  TVector3 t3ele1, t3phot1, t3phot2;
	  t3ele1.SetPtEtaPhi(ptele1,etaele1,phiele1);
	  t3phot1.SetPtEtaPhi(ptPhot[firstfourisophot.at(0)],etaPhot[firstfourisophot.at(0)],phiPhot[firstfourisophot.at(0)]);
	  t3phot2.SetPtEtaPhi(ptPhot[firstfourisophot.at(1)],etaPhot[firstfourisophot.at(1)],phiPhot[firstfourisophot.at(1)]);
	  float eneEle1  = ptele1/(fabs(sin(t3ele1.Theta())));
	  float enePhot1 = ptPhot[firstfourisophot.at(0)]/(fabs(sin(t3phot1.Theta())));
	  float enePhot2 = ptPhot[firstfourisophot.at(1)]/(fabs(sin(t3phot2.Theta())));
	  TLorentzVector t4ele1, t4phot1, t4phot2;
	  t4ele1.SetPtEtaPhiE(ptele1,etaele1,phiele1,eneEle1);
	  t4phot1.SetPtEtaPhiE(ptPhot[firstfourisophot.at(0)],etaPhot[firstfourisophot.at(0)],phiPhot[firstfourisophot.at(0)],enePhot1);
	  t4phot2.SetPtEtaPhiE(ptPhot[firstfourisophot.at(1)],etaPhot[firstfourisophot.at(1)],phiPhot[firstfourisophot.at(1)],enePhot2);
	  invMassele1g1 = (t4phot1 + t4ele1).M();
	  invMassele1g2 = (t4phot2 + t4ele1).M();
	  
	  if (LEPTONS_2012) {
	    oEmoPele1      = fabs(1./electron_ecalEnergy[firstEle] - 1./electron_trackPatVtx[firstEle]);
	    mvanotrigele1  = electron_mvaNonTrig[firstEle];    
	    mvatrigele1    = electron_mvaTrig[firstEle];       
	    matchconvele1  = electron_matchedConv[firstEle];   
	    chHadIso03ele1 = electron_chHad03Iso[firstEle];    
	    nHadIso03ele1  = electron_nHad03Iso[firstEle];     
	    photIso03ele1  = electron_phot03Iso[firstEle];     
	  }
	  if (LEPTONS_2011) {
	    dcotele1 = electron_dcot[firstEle];
	    distele1 = electron_dist[firstEle];
	    float fullHcal = electron_hcalIso03[firstEle] + electron_HoE[firstEle]*electron_sc_energy[firstEle]/cosh(electron_sc_eta[firstEle]);
	    if (fabs(electron_sc_eta[firstEle])<1.4442) {
	      isoele1 = electron_trkIso03[firstEle] + std::max(0.,(electron_ecalIso03[firstEle]-1.)) + electron_hcalIso03[firstEle] - rhoPF*TMath::Pi()*0.3*0.3; 
	      fullisoele1 = electron_trkIso03[firstEle] + std::max(0.,(electron_ecalIso03[firstEle]-1.)) + fullHcal - rhoPF*TMath::Pi()*0.3*0.3; 
	    } else {
	      isoele1 = electron_trkIso03[firstEle] + electron_ecalIso03[firstEle] + electron_hcalIso03[firstEle] - rhoPF*TMath::Pi()*0.3*0.3; 
	      fullisoele1 = electron_trkIso03[firstEle] + electron_ecalIso03[firstEle] + fullHcal - rhoPF*TMath::Pi()*0.3*0.3; 
	    }
	  }
	  
	} else {
	  ptele1    = -500.;
	  etaele1   = -500.;
	  phiele1   = -500.;
	  eneele1   = -500.;
	  sIeIeele1 = -500.;
	  dphiele1  = -500.;
	  detaele1  = -500.;
	  mhitsele1 = -500;
	  d0ele1    = -500.;
	  dzele1    = -500.;
	  invMassele1g1 = -500.;
	  invMassele1g2 = -500.;

	  if (LEPTONS_2012) {
	    oEmoPele1      = -500.;  
	    mvanotrigele1  = -500.;
	    mvatrigele1    = -500.;
	    matchconvele1  = -500;
	    chHadIso03ele1 = -500.;
	    nHadIso03ele1  = -500.;
	    photIso03ele1  = -500.;
	  }
	  if (LEPTONS_2011) {
	    dcotele1    = -500.;
	    distele1    = -500.;
	    fullisoele1 = -500.; 
	    isoele1     = -500.;
	  }
	} 	    

	if (secondEle>=0) {
	  ptele2    = electron_pt[secondEle];
	  etaele2   = electron_sc_eta[secondEle];
	  phiele2   = electron_phi[secondEle];
	  eneele2   = electron_energy[secondEle];
	  sIeIeele2 = electron_SigmaIetaIeta[secondEle];
	  dphiele2  = electron_dPhiIn[secondEle];
	  detaele2  = electron_dEtaIn[secondEle];
	  mhitsele2 = electron_misHits[secondEle];
	  d0ele2    = eleDxyPV(secondEle,vrankPhotonPairs[0]); 
	  dzele2    = eleDzPV(secondEle,vrankPhotonPairs[0]); 

	  TVector3 t3ele2, t3phot1, t3phot2;
	  t3ele2.SetPtEtaPhi(ptele2,etaele2,phiele2);
	  t3phot1.SetPtEtaPhi(ptPhot[firstfourisophot.at(0)],etaPhot[firstfourisophot.at(0)],phiPhot[firstfourisophot.at(0)]);
	  t3phot2.SetPtEtaPhi(ptPhot[firstfourisophot.at(1)],etaPhot[firstfourisophot.at(1)],phiPhot[firstfourisophot.at(1)]);
	  float eneEle2  = ptele2/(fabs(sin(t3ele2.Theta())));
	  float enePhot1 = ptPhot[firstfourisophot.at(0)]/(fabs(sin(t3phot1.Theta())));
	  float enePhot2 = ptPhot[firstfourisophot.at(1)]/(fabs(sin(t3phot2.Theta())));
	  TLorentzVector t4ele2, t4phot1, t4phot2;
	  t4ele2.SetPtEtaPhiE(ptele2,etaele2,phiele2,eneEle2);
	  t4phot1.SetPtEtaPhiE(ptPhot[firstfourisophot.at(0)],etaPhot[firstfourisophot.at(0)],phiPhot[firstfourisophot.at(0)],enePhot1);
	  t4phot2.SetPtEtaPhiE(ptPhot[firstfourisophot.at(1)],etaPhot[firstfourisophot.at(1)],phiPhot[firstfourisophot.at(1)],enePhot2);
	  invMassele2g1 = (t4phot1 + t4ele2).M();
	  invMassele2g2 = (t4phot2 + t4ele2).M();

	  if (LEPTONS_2012) {
	    oEmoPele2      = fabs(1./electron_ecalEnergy[secondEle] - 1./electron_trackPatVtx[secondEle]);
	    mvanotrigele2  = electron_mvaNonTrig[secondEle];    
	    mvatrigele2    = electron_mvaTrig[secondEle];       
	    matchconvele2  = electron_matchedConv[secondEle];   
	    chHadIso03ele2 = electron_chHad03Iso[secondEle];    
	    nHadIso03ele2  = electron_nHad03Iso[secondEle];     
	    photIso03ele2  = electron_phot03Iso[secondEle];     
	  }
	  if (LEPTONS_2011) {
	    dcotele2  = electron_dcot[secondEle];
	    distele2  = electron_dist[secondEle];
	    float fullHcal = electron_hcalIso03[secondEle] + electron_HoE[secondEle]*electron_sc_energy[secondEle]/cosh(electron_sc_eta[secondEle]);	  
	    if (fabs(electron_sc_eta[secondEle])<1.4442) {
	      isoele2 = electron_trkIso03[secondEle] + std::max(0.,(electron_ecalIso03[secondEle]-1.)) + electron_hcalIso03[secondEle] - rhoPF*TMath::Pi()*0.3*0.3; 
	      fullisoele2 = electron_trkIso03[secondEle] + std::max(0.,(electron_ecalIso03[secondEle]-1.)) + fullHcal - rhoPF*TMath::Pi()*0.3*0.3; 
	    } else {
	      isoele2     = electron_trkIso03[secondEle] + electron_ecalIso03[secondEle] + electron_hcalIso03[secondEle] - rhoPF*TMath::Pi()*0.3*0.3; 
	      fullisoele2 = electron_trkIso03[secondEle] + electron_ecalIso03[secondEle] + fullHcal - rhoPF*TMath::Pi()*0.3*0.3; 
	    }
	  }
	  
	} else {
	  ptele2    = -500.;
	  etaele2   = -500.;
	  phiele2   = -500.;
	  eneele2   = -500.;
	  sIeIeele2 = -500.;
	  dphiele2  = -500.;
	  detaele2  = -500.;
	  mhitsele2 = -500;
	  d0ele2    = -500.;
	  dzele2    = -500.;
	  invMassele2g1 = -500.;
	  invMassele2g2 = -500.;

	  if(LEPTONS_2012) {
	    oEmoPele2      = -500.; 
	    mvanotrigele2  = -500.;
	    mvatrigele2    = -500.;
	    matchconvele2  = -500;
	    chHadIso03ele2 = -500.;
	    nHadIso03ele2  = -500.;
	    photIso03ele2  = -500.;
	  }
	  if(LEPTONS_2011) {
	    dcotele2    = -500.;
	    distele2    = -500.;
	    isoele2     = -500.;
	    fullisoele2 = -500.;
	  }
	}	 

	// loose electron selection
	if (firstEleLoose>=0) {
	  pteleloose1    = electron_pt[firstEleLoose];
	  etaeleloose1   = electron_sc_eta[firstEleLoose];
	  phieleloose1   = electron_phi[firstEleLoose];
	  eneeleloose1   = electron_energy[firstEleLoose];
	  sIeIeeleloose1 = electron_SigmaIetaIeta[firstEleLoose];
	  dphieleloose1  = electron_dPhiIn[firstEleLoose];
	  detaeleloose1  = electron_dEtaIn[firstEleLoose];
	  mhitseleloose1 = electron_misHits[firstEleLoose];
	  d0eleloose1    = eleDxyPV(firstEleLoose,vrankPhotonPairs[0]); 
	  dzeleloose1    = eleDzPV(firstEleLoose,vrankPhotonPairs[0]); 

	  TVector3 t3eleloose1, t3phot1, t3phot2;
	  t3eleloose1.SetPtEtaPhi(pteleloose1,etaeleloose1,phieleloose1);
	  t3phot1.SetPtEtaPhi(ptPhot[firstfourisophot.at(0)],etaPhot[firstfourisophot.at(0)],phiPhot[firstfourisophot.at(0)]);
	  t3phot2.SetPtEtaPhi(ptPhot[firstfourisophot.at(1)],etaPhot[firstfourisophot.at(1)],phiPhot[firstfourisophot.at(1)]);
	  float eneEleloose1  = pteleloose1/(fabs(sin(t3eleloose1.Theta())));
	  float enePhot1 = ptPhot[firstfourisophot.at(0)]/(fabs(sin(t3phot1.Theta())));
	  float enePhot2 = ptPhot[firstfourisophot.at(1)]/(fabs(sin(t3phot2.Theta())));
	  TLorentzVector t4eleloose1, t4phot1, t4phot2;
	  t4eleloose1.SetPtEtaPhiE(pteleloose1,etaeleloose1,phieleloose1,eneEleloose1);
	  t4phot1.SetPtEtaPhiE(ptPhot[firstfourisophot.at(0)],etaPhot[firstfourisophot.at(0)],phiPhot[firstfourisophot.at(0)],enePhot1);
	  t4phot2.SetPtEtaPhiE(ptPhot[firstfourisophot.at(1)],etaPhot[firstfourisophot.at(1)],phiPhot[firstfourisophot.at(1)],enePhot2);
	  invMasseleloose1g1 = (t4phot1 + t4eleloose1).M();
	  invMasseleloose1g2 = (t4phot2 + t4eleloose1).M();

	  if (LEPTONS_2012) {
	    oEmoPeleloose1      = fabs(1./electron_ecalEnergy[firstEleLoose] - 1./electron_trackPatVtx[firstEleLoose]);
	    mvanotrigeleloose1  = electron_mvaNonTrig[firstEleLoose];    
	    mvatrigeleloose1    = electron_mvaTrig[firstEleLoose];       
	    matchconveleloose1  = electron_matchedConv[firstEleLoose];   
	    chHadIso03eleloose1 = electron_chHad03Iso[firstEleLoose];    
	    nHadIso03eleloose1  = electron_nHad03Iso[firstEleLoose];     
	    photIso03eleloose1  = electron_phot03Iso[firstEleLoose];     
	  }
	  if (LEPTONS_2011) {
	    dcoteleloose1 = electron_dcot[firstEleLoose];
	    disteleloose1 = electron_dist[firstEleLoose];
	    float fullHcal = electron_hcalIso03[firstEleLoose] + electron_HoE[firstEleLoose]*electron_sc_energy[firstEleLoose]/cosh(electron_sc_eta[firstEleLoose]);
	    if (fabs(electron_sc_eta[firstEleLoose])<1.4442) {
	      isoeleloose1 = electron_trkIso03[firstEleLoose] + std::max(0.,(electron_ecalIso03[firstEleLoose]-1.)) + electron_hcalIso03[firstEleLoose] - rhoPF*TMath::Pi()*0.3*0.3; 
	      fullisoeleloose1 = electron_trkIso03[firstEleLoose] + std::max(0.,(electron_ecalIso03[firstEleLoose]-1.)) + fullHcal - rhoPF*TMath::Pi()*0.3*0.3; 
	    } else {
	      isoeleloose1 = electron_trkIso03[firstEleLoose] + electron_ecalIso03[firstEleLoose] + electron_hcalIso03[firstEleLoose] - rhoPF*TMath::Pi()*0.3*0.3; 
	      fullisoeleloose1 = electron_trkIso03[firstEleLoose] + electron_ecalIso03[firstEleLoose] + fullHcal - rhoPF*TMath::Pi()*0.3*0.3; 
	    }
	  }
	  
	} else {
	  pteleloose1    = -500.;
	  etaeleloose1   = -500.;
	  phieleloose1   = -500.;
	  eneeleloose1   = -500.;
	  sIeIeeleloose1 = -500.;
	  dphieleloose1  = -500.;
	  detaeleloose1  = -500.;
	  mhitseleloose1 = -500;
	  d0eleloose1    = -500.;
	  dzeleloose1    = -500.;
	  invMasseleloose1g1 = -500.;
	  invMasseleloose1g2 = -500.;
	  if (LEPTONS_2012) {
	    oEmoPeleloose1      = -500.;  
	    mvanotrigeleloose1  = -500.;
	    mvatrigeleloose1    = -500.;
	    matchconveleloose1  = -500;
	    chHadIso03eleloose1 = -500.;
	    nHadIso03eleloose1  = -500.;
	    photIso03eleloose1  = -500.;
	  }
	  if (LEPTONS_2011) {
	    dcoteleloose1    = -500.;
	    disteleloose1    = -500.;
	    fullisoeleloose1 = -500.; 
	    isoeleloose1     = -500.;
	  }
	} 	    

	if (secondEleLoose>=0) {
	  pteleloose2    = electron_pt[secondEleLoose];
	  etaeleloose2   = electron_sc_eta[secondEleLoose];
	  phieleloose2   = electron_phi[secondEleLoose];
	  eneeleloose2   = electron_energy[secondEleLoose];
	  sIeIeeleloose2 = electron_SigmaIetaIeta[secondEleLoose];
	  dphieleloose2  = electron_dPhiIn[secondEleLoose];
	  detaeleloose2  = electron_dEtaIn[secondEleLoose];
	  mhitseleloose2 = electron_misHits[secondEleLoose];
	  d0eleloose2    = eleDxyPV(secondEleLoose,vrankPhotonPairs[0]); 
	  dzeleloose2    = eleDzPV(secondEleLoose,vrankPhotonPairs[0]); 

	  TVector3 t3eleloose2, t3phot1, t3phot2;
	  t3eleloose2.SetPtEtaPhi(pteleloose2,etaeleloose2,phieleloose2);
	  t3phot1.SetPtEtaPhi(ptPhot[firstfourisophot.at(0)],etaPhot[firstfourisophot.at(0)],phiPhot[firstfourisophot.at(0)]);
	  t3phot2.SetPtEtaPhi(ptPhot[firstfourisophot.at(1)],etaPhot[firstfourisophot.at(1)],phiPhot[firstfourisophot.at(1)]);
	  float eneEleloose2  = pteleloose2/(fabs(sin(t3eleloose2.Theta())));
	  float enePhot1 = ptPhot[firstfourisophot.at(0)]/(fabs(sin(t3phot1.Theta())));
	  float enePhot2 = ptPhot[firstfourisophot.at(1)]/(fabs(sin(t3phot2.Theta())));
	  TLorentzVector t4eleloose2, t4phot1, t4phot2;
	  t4eleloose2.SetPtEtaPhiE(pteleloose2,etaeleloose2,phieleloose2,eneEleloose2);
	  t4phot1.SetPtEtaPhiE(ptPhot[firstfourisophot.at(0)],etaPhot[firstfourisophot.at(0)],phiPhot[firstfourisophot.at(0)],enePhot1);
	  t4phot2.SetPtEtaPhiE(ptPhot[firstfourisophot.at(1)],etaPhot[firstfourisophot.at(1)],phiPhot[firstfourisophot.at(1)],enePhot2);
	  invMasseleloose2g1 = (t4phot1 + t4eleloose2).M();
	  invMasseleloose2g2 = (t4phot2 + t4eleloose2).M();

	  if (LEPTONS_2012) {
	    oEmoPeleloose2      = fabs(1./electron_ecalEnergy[secondEleLoose] - 1./electron_trackPatVtx[secondEleLoose]);
	    mvanotrigeleloose2  = electron_mvaNonTrig[secondEleLoose];    
	    mvatrigeleloose2    = electron_mvaTrig[secondEleLoose];       
	    matchconveleloose2  = electron_matchedConv[secondEleLoose];   
	    chHadIso03eleloose2 = electron_chHad03Iso[secondEleLoose];    
	    nHadIso03eleloose2  = electron_nHad03Iso[secondEleLoose];     
	    photIso03eleloose2  = electron_phot03Iso[secondEleLoose];     
	  }
	  if (LEPTONS_2011) {
	    dcoteleloose2  = electron_dcot[secondEleLoose];
	    disteleloose2  = electron_dist[secondEleLoose];
	    float fullHcal = electron_hcalIso03[secondEleLoose] + electron_HoE[secondEleLoose]*electron_sc_energy[secondEleLoose]/cosh(electron_sc_eta[secondEleLoose]);	  
	    if (fabs(electron_sc_eta[secondEleLoose])<1.4442) {
	      isoeleloose2 = electron_trkIso03[secondEleLoose] + std::max(0.,(electron_ecalIso03[secondEleLoose]-1.)) + electron_hcalIso03[secondEleLoose] - rhoPF*TMath::Pi()*0.3*0.3; 
	      fullisoeleloose2 = electron_trkIso03[secondEleLoose] + std::max(0.,(electron_ecalIso03[secondEleLoose]-1.)) + fullHcal - rhoPF*TMath::Pi()*0.3*0.3; 
	    } else {
	      isoeleloose2     = electron_trkIso03[secondEleLoose] + electron_ecalIso03[secondEleLoose] + electron_hcalIso03[secondEleLoose] - rhoPF*TMath::Pi()*0.3*0.3; 
	      fullisoeleloose2 = electron_trkIso03[secondEleLoose] + electron_ecalIso03[secondEleLoose] + fullHcal - rhoPF*TMath::Pi()*0.3*0.3; 
	    }
	  }
	  
	} else {
	  pteleloose2    = -500.;
	  etaeleloose2   = -500.;
	  phieleloose2   = -500.;
	  eneeleloose2   = -500.;
	  sIeIeeleloose2 = -500.;
	  dphieleloose2  = -500.;
	  detaeleloose2  = -500.;
	  mhitseleloose2 = -500;
	  d0eleloose2    = -500.;
	  dzeleloose2    = -500.;
	  invMasseleloose2g1 = -500.;
	  invMasseleloose2g2 = -500.;
	  if(LEPTONS_2012) {
	    oEmoPeleloose2      = -500.; 
	    mvanotrigeleloose2  = -500.;
	    mvatrigeleloose2    = -500.;
	    matchconveleloose2  = -500;
	    chHadIso03eleloose2 = -500.;
	    nHadIso03eleloose2  = -500.;
	    photIso03eleloose2  = -500.;
	  }
	  if(LEPTONS_2011) {
	    dcoteleloose2    = -500.;
	    disteleloose2    = -500.;
	    isoeleloose2     = -500.;
	    fullisoeleloose2 = -500.;
	  }
	}	 

	// muons: tight selection
	if (firstMu>=0) {
	  ptmu1      = Muon_pt[firstMu];
	  etamu1     = Muon_eta[firstMu];
	  phimu1     = Muon_phi[firstMu];
	  enemu1     = Muon_energy[firstMu];
	  pixhitsmu1 = Muon_pixHits[firstMu];
	  trkhitsmu1 = Muon_tkHits[firstMu];
	  hitsmu1    = Muon_validHits[firstMu];
	  chi2mu1    = Muon_normChi2[firstMu];
	  matchmu1   = Muon_numberOfMatches[firstMu];
	  d0mu1      = muonDxyPV(firstMu,vrankPhotonPairs[0]);
	  dzmu1      = muonDzPV(firstMu,vrankPhotonPairs[0]);

	  if(LEPTONS_2012) {
	    chHadmu1 = Muon_pfiso03_chHad[firstMu];
	    nHadmu1  = Muon_pfiso03_nHad[firstMu];
	    photmu1  = Muon_pfiso03_Phot[firstMu];
	    puptmu1  = Muon_pfiso03_PUPt[firstMu];
	  }
	  if(LEPTONS_2011) {
	    isomu1 = Muon_trackIso[firstMu] + Muon_ecalIso[firstMu] + Muon_hcalIso[firstMu] - rhoPF*TMath::Pi()*0.3*0.3; 
	  }
	} else {
	  ptmu1      = -500.;
	  etamu1     = -500.;
	  phimu1     = -500.;
	  enemu1     = -500.;
	  pixhitsmu1 = -500;
	  trkhitsmu1 = -500;
	  hitsmu1    = -500;
	  chi2mu1    = -500.;
	  matchmu1   = -500;
	  d0mu1      = -500.;
	  dzmu1      = -500.;

	  if(LEPTONS_2012) {
	    chHadmu1 = -500.;
	    nHadmu1  = -500.;
	    photmu1  = -500.;
	    puptmu1  = -500.;
	  }
	  if(LEPTONS_2011) {
	    isomu1 = -500.;
	  }
	}

	if (secondMu>=0) {
	  ptmu2      = Muon_pt[secondMu];
	  etamu2     = Muon_eta[secondMu];
	  phimu2     = Muon_phi[secondMu];
	  enemu2     = Muon_energy[secondMu];
	  pixhitsmu2 = Muon_pixHits[secondMu];
	  trkhitsmu2 = Muon_tkHits[secondMu];
	  hitsmu2    = Muon_validHits[secondMu];
	  chi2mu2    = Muon_normChi2[secondMu];
	  matchmu2   = Muon_numberOfMatches[secondMu];
	  d0mu2      = muonDxyPV(secondMu,vrankPhotonPairs[0]);
	  dzmu2      = muonDzPV(secondMu,vrankPhotonPairs[0]);

	  if (LEPTONS_2012) {
	    chHadmu2 = Muon_pfiso03_chHad[secondMu];
	    nHadmu2  = Muon_pfiso03_nHad[secondMu];
	    photmu2  = Muon_pfiso03_Phot[secondMu];
	    puptmu2  = Muon_pfiso03_PUPt[secondMu];
	  }
	  if (LEPTONS_2011) {
	    isomu2 = Muon_trackIso[secondMu] + Muon_ecalIso[secondMu] + Muon_hcalIso[secondMu] - rhoPF*TMath::Pi()*0.3*0.3;  
	  }
	} else {
	  ptmu2      = -500.;
	  etamu2     = -500.;
	  phimu2     = -500.;
	  enemu2     = -500.;
	  pixhitsmu2 = -500;
	  trkhitsmu2 = -500;
	  hitsmu2    = -500;
	  chi2mu2    = -500.;
	  matchmu2   = -500;
	  d0mu2      = -500.;
	  dzmu2      = -500.;

	  if (LEPTONS_2012) {
	    chHadmu2   = -500.;
	    nHadmu2    = -500.;
	    photmu2    = -500.;
	    puptmu2    = -500.;
	  }
	  if (LEPTONS_2011) {
	    isomu2 = -500.;
	  }
	}


	// muons: loose selection
	if (firstMuLoose>=0) {
	  ptmuloose1      = Muon_pt[firstMuLoose];
	  etamuloose1     = Muon_eta[firstMuLoose];
	  phimuloose1     = Muon_phi[firstMuLoose];
	  enemuloose1     = Muon_energy[firstMuLoose];
	  pixhitsmuloose1 = Muon_pixHits[firstMuLoose];
	  trkhitsmuloose1 = Muon_tkHits[firstMuLoose];
	  hitsmuloose1    = Muon_validHits[firstMuLoose];
	  chi2muloose1    = Muon_normChi2[firstMuLoose];
	  matchmuloose1   = Muon_numberOfMatches[firstMuLoose];
	  d0muloose1      = muonDxyPV(firstMuLoose,vrankPhotonPairs[0]);
	  dzmuloose1      = muonDzPV(firstMuLoose,vrankPhotonPairs[0]);

	  if(LEPTONS_2012) {
	    chHadmuloose1 = Muon_pfiso03_chHad[firstMuLoose];
	    nHadmuloose1  = Muon_pfiso03_nHad[firstMuLoose];
	    photmuloose1  = Muon_pfiso03_Phot[firstMuLoose];
	    puptmuloose1  = Muon_pfiso03_PUPt[firstMuLoose];
	  }
	  if(LEPTONS_2011) {
	    isomuloose1 = Muon_trackIso[firstMuLoose] + Muon_ecalIso[firstMuLoose] + Muon_hcalIso[firstMuLoose] - rhoPF*TMath::Pi()*0.3*0.3; 
	  }
	} else {
	  ptmuloose1      = -500.;
	  etamuloose1     = -500.;
	  phimuloose1     = -500.;
	  enemuloose1     = -500.;
	  pixhitsmuloose1 = -500;
	  trkhitsmuloose1 = -500;
	  hitsmuloose1    = -500;
	  chi2muloose1    = -500.;
	  matchmuloose1   = -500;
	  d0muloose1      = -500.;
	  dzmuloose1      = -500.;

	  if(LEPTONS_2012) {
	    chHadmuloose1 = -500.;
	    nHadmuloose1  = -500.;
	    photmuloose1  = -500.;
	    puptmuloose1  = -500.;
	  }
	  if(LEPTONS_2011) {
	    isomuloose1 = -500.;
	  }
	}

	if (secondMuLoose>=0) {
	  ptmuloose2      = Muon_pt[secondMuLoose];
	  etamuloose2     = Muon_eta[secondMuLoose];
	  phimuloose2     = Muon_phi[secondMuLoose];
	  enemuloose2     = Muon_energy[secondMuLoose];
	  pixhitsmuloose2 = Muon_pixHits[secondMuLoose];
	  trkhitsmuloose2 = Muon_tkHits[secondMuLoose];
	  hitsmuloose2    = Muon_validHits[secondMuLoose];
	  chi2muloose2    = Muon_normChi2[secondMuLoose];
	  matchmuloose2   = Muon_numberOfMatches[secondMuLoose];
	  d0muloose2      = muonDxyPV(secondMuLoose,vrankPhotonPairs[0]);
	  dzmuloose2      = muonDzPV(secondMuLoose,vrankPhotonPairs[0]);

	  if (LEPTONS_2012) {
	    chHadmuloose2 = Muon_pfiso03_chHad[secondMuLoose];
	    nHadmuloose2  = Muon_pfiso03_nHad[secondMuLoose];
	    photmuloose2  = Muon_pfiso03_Phot[secondMuLoose];
	    puptmuloose2  = Muon_pfiso03_PUPt[secondMuLoose];
	  }
	  if (LEPTONS_2011) {
	    isomuloose2 = Muon_trackIso[secondMuLoose] + Muon_ecalIso[secondMuLoose] + Muon_hcalIso[secondMuLoose] - rhoPF*TMath::Pi()*0.3*0.3;  
	  }
	} else {
	  ptmuloose2      = -500.;
	  etamuloose2     = -500.;
	  phimuloose2     = -500.;
	  enemuloose2     = -500.;
	  pixhitsmuloose2 = -500;
	  trkhitsmuloose2 = -500;
	  hitsmuloose2    = -500;
	  chi2muloose2    = -500.;
	  matchmuloose2   = -500;
	  d0muloose2      = -500.;
	  dzmuloose2      = -500.;

	  if (LEPTONS_2012) {
	    chHadmuloose2   = -500.;
	    nHadmuloose2    = -500.;
	    photmuloose2    = -500.;
	    puptmuloose2    = -500.;
	  }
	  if (LEPTONS_2011) {
	    isomuloose2 = -500.;
	  }
	}
	  
	/***************************************************
        *                                                 *
        *        SAVING RECO VARIABLES IN TTREE           *
        *                                                 *
        ***************************************************/


      double twojetsmassiso(0), etatwojetsiso(-999), phitwojetsiso(-999), pttwojetsiso(-999);
   
      //bool recoPreselection = (firsttwoisophot.at(0)>-1 && firsttwoisophot.at(1)>-1 && ptPhot[firsttwoisophot.at(0)]>20 && ptPhot[firsttwoisophot.at(1)]>20);

      /// firstTWO --> firstten
      bool recoPreselection = ( firstfourisophot.at(0)>-1 && firstfourisophot.at(1)>-1 && ptPhot[firstfourisophot.at(0)]>20 && ptPhot[firstfourisophot.at(1)]>20); 

#ifdef DEBUG
      cout << "[DEBUG] recoPreselection = " << recoPreselection << endl;
#endif
      
      if(!recoPreselection ) { 
        SetAllRecoVarToMinus999();
      }
      else {
	
	if( firsttennoisojet.at(0) > -1) {
	  ptjetisoreco1.Fill(ptCorrJet_pfakt5[firsttennoisojet.at(0)],weight);
	  etajetisoreco1.Fill(etaJet_pfakt5[firsttennoisojet.at(0)],weight);
	}
	if( firsttennoisojet.at(1) > -1) {
	  ptjetisoreco2.Fill(ptCorrJet_pfakt5[firsttennoisojet.at(1)],weight);
	  etajetisoreco2.Fill(etaJet_pfakt5[firsttennoisojet.at(1)],weight);
	}
	if( firsttennoisojet.at(0) > -1 && firsttennoisojet.at(1) > -1) {
	  

#ifdef DEBUG
        cout << "[DEBUG] before recojets" << endl;
#endif

	  TLorentzVector jet1, jet2;	
	  jet1.SetPtEtaPhiE(ptCorrJet_pfakt5[firsttennoisojet.at(0)],etaJet_pfakt5[firsttennoisojet.at(0)],phiJet_pfakt5[firsttennoisojet.at(0)],eJet_pfakt5[firsttennoisojet.at(0)]/ptJet_pfakt5[firsttennoisojet.at(0)]*ptCorrJet_pfakt5[firsttennoisojet.at(0)]);
	  jet2.SetPtEtaPhiE(ptCorrJet_pfakt5[firsttennoisojet.at(1)],etaJet_pfakt5[firsttennoisojet.at(1)],phiJet_pfakt5[firsttennoisojet.at(1)],eJet_pfakt5[firsttennoisojet.at(1)]/ptJet_pfakt5[firsttennoisojet.at(1)]*ptCorrJet_pfakt5[firsttennoisojet.at(1)]);
	  
	  TLorentzVector sum = jet1 + jet2;
	  thejet1 = jet1;
	  thejet2 = jet2;
#ifdef DEBUG
        cout << "[DEBUG] after recojets" << endl;
#endif

	  twojetsmassiso = sum.M();
	  etatwojetsiso = sum.Eta();
	  phitwojetsiso = sum.Phi();
	  pttwojetsiso = sum.Pt();
	  
	  if(jetsyst_) 
	    JECunc.Fill(ptCorrJet_pfakt5[firsttennoisojet.at(0)],jetsyst_->getJESUncertainty(etaJet_pfakt5[firsttennoisojet.at(0)],ptCorrJet_pfakt5[firsttennoisojet.at(0)]));

	  int assjj(-999);
	  for(int j=0; j<nJetGen_akt5; j++){	
	    double DR = sqrt(delta_eta(etaJet_pfakt5[firsttennoisojet.at(0)],etaJetGen_akt5[j])*delta_eta(etaJet_pfakt5[firsttennoisojet.at(0)],etaJetGen_akt5[j]) + 
			     delta_phi(phiJet_pfakt5[firsttennoisojet.at(0)],phiJetGen_akt5[j])*delta_phi(phiJet_pfakt5[firsttennoisojet.at(0)],phiJetGen_akt5[j]) ) ;
	    if(DR < .1 && (TMath::Abs(ptCorrJet_pfakt5[firsttennoisojet.at(0)]-ptJetGen_akt5[j])/ptJetGen_akt5[j] < 0.5)) assjj = j; 
	  }
	  if(assjj>-1 && ptCorrJet_pfakt5[firsttennoisojet.at(0)]>30. && ptCorrJet_pfakt5[firsttennoisojet.at(1)]>20 ){
	    if(TMath::Abs(etaJet_pfakt5[firsttennoisojet.at(0)]-etaJet_pfakt5[firsttennoisojet.at(1)])>2.5)
	      JECresovbf.Fill((ptCorrJet_pfakt5[firsttennoisojet.at(0)]-ptJetGen_akt5[assjj])/ptJetGen_akt5[assjj]);
	    else
	      JECresovh.Fill((ptCorrJet_pfakt5[firsttennoisojet.at(0)]-ptJetGen_akt5[assjj])/ptJetGen_akt5[assjj]);
	  }
	  deltaetajetisoreco.Fill(etaJet_pfakt5[firsttennoisojet.at(0)]-etaJet_pfakt5[firsttennoisojet.at(1)],weight);
	  double aveeta = (etaJet_pfakt5[firsttennoisojet.at(0)]+etaJet_pfakt5[firsttennoisojet.at(1)])/2;
	  double zeppen1 = etaJet_pfakt5[firsttennoisojet.at(0)] - aveeta;
	  double zeppen2 = etaJet_pfakt5[firsttennoisojet.at(1)] - aveeta;
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
 	phigg = phihiggsisonewvtx;
	etagg = etahiggsisonewvtx;
	deltaphigg = delta_phi(phiPhot[firstfourisophot.at(0)],phiPhot[firstfourisophot.at(1)]);
	ptphot1 = ptPhot[firstfourisophot.at(0)]; 
        ptphot2 = ptPhot[firstfourisophot.at(1)];   
	deltaRToTrackphot1 = pid_deltaRToTrackPhot[firstfourisophot.at(0)];
	deltaRToTrackphot2 = pid_deltaRToTrackPhot[firstfourisophot.at(1)];
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
	energyErrphot1 = escRegrPhotError[firstfourisophot.at(0)];
	energyErrphot2 = escRegrPhotError[firstfourisophot.at(1)];  
	energySmearingphot1 = smearEnePhot[firstfourisophot.at(0)];
	energySmearingphot2 = smearEnePhot[firstfourisophot.at(1)];  
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
	idmvaphot1 = isomva.at(firstfourisophot.at(0));
	idmvaphot2 = isomva.at(firstfourisophot.at(1));
	idcicphot1 = isocic.at(firstfourisophot.at(0));
	idcicphot2 = isocic.at(firstfourisophot.at(1));
	idcicnoelvetophot1 = isocicnoelveto.at(firstfourisophot.at(0));
	idcicnoelvetophot2 = isocicnoelveto.at(firstfourisophot.at(1));
	idcicpfphot1 = isocicpf.at(firstfourisophot.at(0));
	idcicpfphot2 = isocicpf.at(firstfourisophot.at(1));
	idcicpfnoelvetophot1 = isocicpfnoelveto.at(firstfourisophot.at(0));
	idcicpfnoelvetophot2 = isocicpfnoelveto.at(firstfourisophot.at(1));
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

/*
	if( firsttennoisojet.at(0) > -1) {
	  ptjet1 = ptJet_pfakt5[firsttennoisojet.at(0)];
	  ptcorrjet1 = ptCorrJet_pfakt5[firsttennoisojet.at(0)];	  
	  ecorrjet1 = ptcorrjet1/ptjet1*eJet_pfakt5[firsttennoisojet.at(0)];
	  etajet1 = etaJet_pfakt5[firsttennoisojet.at(0)];
	  phijet1 = phiJet_pfakt5[firsttennoisojet.at(0)];
	  betajet1 = beta_pfakt5[firsttennoisojet.at(0)][vrankPhotonPairs[0]];
	  betastarjet1 = betaStar_pfakt5[firsttennoisojet.at(0)][vrankPhotonPairs[0]];
	  assjet1 = assoJet(firsttennoisojet.at(0));
 	  btagvtxjet1 = simpleSecondaryVertexHighEffBJetTags[firsttennoisojet.at(0)];
 	  btagtrkjet1 = trackCountingHighEffBJetTags[firsttennoisojet.at(0)];	  
 	  btagjprobjet1 = jetProbabilityBJetTags[firsttennoisojet.at(0)];	  
 	  ptDjet1 = ptDJet_pfakt5[firsttennoisojet.at(0)];
	  rmsjet1 = rmsCandJet_pfakt5[firsttennoisojet.at(0)];
 	  ntrkjet1 = nChargedHadrons_pfakt5[firsttennoisojet.at(0)];
 	  nneutjet1 = nPhotons_pfakt5[firsttennoisojet.at(0)] + nNeutralHadrons_pfakt5[firsttennoisojet.at(0)] + nHFHadrons_pfakt5[firsttennoisojet.at(0)] + nHFEM_pfakt5[firsttennoisojet.at(0)];
 	  jetIdSimple_mvajet1 = jetIdSimple_mva_pfakt5[firsttennoisojet.at(0)];
 	  jetIdFull_mvajet1 = jetIdFull_mva_pfakt5[firsttennoisojet.at(0)];
 	  jetId_dR2Meanjet1 = jetId_dR2Mean_pfakt5[firsttennoisojet.at(0)];
 	  jetId_betaStarClassicjet1 = jetId_betaStarClassic_pfakt5[firsttennoisojet.at(0)];
 	  jetIdCutBased_wpjet1 = jetIdCutBased_wp_pfakt5[firsttennoisojet.at(0)];
 	  jetIdSimple_wpjet1 = jetIdSimple_wp_pfakt5[firsttennoisojet.at(0)];	  
 	  jetIdFull_wpjet1 = jetIdFull_wp_pfakt5[firsttennoisojet.at(0)];	  
	  jetId_frac01jet1 = jetId_frac01_pfakt5[firsttennoisojet.at(0)];
	  jetId_frac02jet1 = jetId_frac02_pfakt5[firsttennoisojet.at(0)];
	  jetId_frac03jet1 = jetId_frac03_pfakt5[firsttennoisojet.at(0)];
	  jetId_frac04jet1 = jetId_frac04_pfakt5[firsttennoisojet.at(0)];
	  jetId_frac05jet1 = jetId_frac05_pfakt5[firsttennoisojet.at(0)];
          jetId_betajet1 = jetId_beta_pfakt5[firsttennoisojet.at(0)];
	  jetId_betaStarjet1 = jetId_betaStar_pfakt5[firsttennoisojet.at(0)];


     // match to parton
     Float_t deltaRMCmin = 999.;
     Int_t pdgIdPart_found = 0;
     Int_t pdgIdMomPart_found = 0;

     for(Int_t iPartMC=0; iPartMC<nMC; ++iPartMC) {

       if( statusMC[iPartMC]!=3 ) continue;

       if( ptMC[iPartMC]<0.1 ) continue;

       TLorentzVector jet;
       jet.SetPtEtaPhiE( ptcorrjet1, etajet1, phijet1, ecorrjet1 );
       TLorentzVector parton;
       parton.SetPtEtaPhiE( ptMC[iPartMC], etaMC[iPartMC], phiMC[iPartMC], eMC[iPartMC] );

       Int_t pdgId = pdgIdMC[iPartMC];
     
       Float_t deltaRMC = jet.DeltaR(parton);

       bool goodPdgId = ( (fabs(pdgId)<=9) || (fabs(pdgId)==21) );
       if( !goodPdgId ) continue;
     
       if( (deltaRMC < deltaRMCmin) && goodPdgId ) {
         deltaRMCmin = deltaRMC;
         pdgIdPart_found = pdgIdMC[iPartMC];
         pdgIdMomPart_found = pdgIdMC[motherIDMC[iPartMC]];
       }

     } //for MC particles


     partPdgIDjet1 = ( deltaRMCmin<0.5 ) ? pdgIdPart_found : -999; 
     partMomPdgIDjet1 = ( deltaRMCmin<0.5 ) ? pdgIdMomPart_found : -999; 




	}else{
	  ptjet1 = -999;
	  ptcorrjet1 = -999;
	  ecorrjet1 = -999;
	  etajet1 = -999;	 
	  phijet1 = -999;	 
	  betajet1 = -999.;
	  betastarjet1 = -999.;
	  assjet1 = -999.;
 	  btagvtxjet1 = -999.;
 	  btagtrkjet1 = -999.;
 	  ptDjet1 = -999.;
	  rmsjet1 = -999.;
 	  ntrkjet1 = -999.;
 	  nneutjet1 = -999.; 
 	  jetIdSimple_mvajet1 = -999.; 
 	  jetIdFull_mvajet1   = -999.; 
 	  jetId_dR2Meanjet1   = -999.; 
 	  jetId_betaStarClassicjet1 = -999.; 
 	  jetIdCutBased_wpjet1 = -999.; 
 	  jetIdSimple_wpjet1 = -999.; 
 	  jetIdFull_wpjet1 = -999.; 
 	  jetId_frac01jet1 = -999.; 
 	  jetId_frac02jet1 = -999.; 
 	  jetId_frac03jet1 = -999.; 
 	  jetId_frac04jet1 = -999.; 
 	  jetId_frac05jet1 = -999.; 
 	  jetId_betajet1 = -999.; 
 	  jetId_betaStarjet1 = -999.; 
        partPdgIDjet1 = -999;
        partMomPdgIDjet1 = -999;
	}
	if( firsttennoisojet.at(1) > -1) {
	  ptjet2 = ptJet_pfakt5[firsttennoisojet.at(1)];
	  ptcorrjet2 = ptCorrJet_pfakt5[firsttennoisojet.at(1)];	  
	  ecorrjet2 = ptcorrjet2/ptjet2*eJet_pfakt5[firsttennoisojet.at(1)];
	  etajet2 = etaJet_pfakt5[firsttennoisojet.at(1)];
	  phijet2 = phiJet_pfakt5[firsttennoisojet.at(1)];
	  betajet2 = beta_pfakt5[firsttennoisojet.at(1)][vrankPhotonPairs[0]];
	  betastarjet2 = betaStar_pfakt5[firsttennoisojet.at(1)][vrankPhotonPairs[0]];
	  assjet2 = assoJet(firsttennoisojet.at(1));
 	  btagvtxjet2 = simpleSecondaryVertexHighEffBJetTags[firsttennoisojet.at(1)];
 	  btagtrkjet2 = trackCountingHighEffBJetTags[firsttennoisojet.at(1)];	  
 	  btagjprobjet2 = jetProbabilityBJetTags[firsttennoisojet.at(1)];	  
 	  ptDjet2 = ptDJet_pfakt5[firsttennoisojet.at(1)];
	  rmsjet2 = rmsCandJet_pfakt5[firsttennoisojet.at(1)];
 	  ntrkjet2 = nChargedHadrons_pfakt5[firsttennoisojet.at(1)];
	  // 	  nneutjet2 = nNeutralHadrons_pfakt5[firsttennoisojet.at(1)];
 	  nneutjet2 = nPhotons_pfakt5[firsttennoisojet.at(1)] + nNeutralHadrons_pfakt5[firsttennoisojet.at(1)] + nHFHadrons_pfakt5[firsttennoisojet.at(1)] + nHFEM_pfakt5[firsttennoisojet.at(1)];
	  jetIdSimple_mvajet2 = jetIdSimple_mva_pfakt5[firsttennoisojet.at(1)];
 	  jetIdFull_mvajet2 = jetIdFull_mva_pfakt5[firsttennoisojet.at(1)];
 	  jetId_dR2Meanjet2 = jetId_dR2Mean_pfakt5[firsttennoisojet.at(1)];
 	  jetId_betaStarClassicjet2 = jetId_betaStarClassic_pfakt5[firsttennoisojet.at(1)];
 	  jetIdCutBased_wpjet2 = jetIdCutBased_wp_pfakt5[firsttennoisojet.at(1)];
 	  jetIdSimple_wpjet2 = jetIdSimple_wp_pfakt5[firsttennoisojet.at(1)];	  
 	  jetIdFull_wpjet2 = jetIdFull_wp_pfakt5[firsttennoisojet.at(1)];	  
	  jetId_frac01jet2 = jetId_frac01_pfakt5[firsttennoisojet.at(1)];
	  jetId_frac02jet2 = jetId_frac02_pfakt5[firsttennoisojet.at(1)];
	  jetId_frac03jet2 = jetId_frac03_pfakt5[firsttennoisojet.at(1)];
	  jetId_frac04jet2 = jetId_frac04_pfakt5[firsttennoisojet.at(1)];
	  jetId_frac05jet2 = jetId_frac05_pfakt5[firsttennoisojet.at(1)];
          jetId_betajet2 = jetId_beta_pfakt5[firsttennoisojet.at(1)];
	  jetId_betaStarjet2 = jetId_betaStar_pfakt5[firsttennoisojet.at(1)];
     
     
     // match to parton
     Float_t deltaRMCmin = 999.;
     Int_t pdgIdPart_found = 0;
     Int_t pdgIdMomPart_found = 0;

     for(Int_t iPartMC=0; iPartMC<nMC; ++iPartMC) {

       if( statusMC[iPartMC]!=3 ) continue;

       if( ptMC[iPartMC]<0.1 ) continue;

       TLorentzVector jet;
       jet.SetPtEtaPhiE( ptcorrjet2, etajet2, phijet2, ecorrjet2 );
       TLorentzVector parton;
       parton.SetPtEtaPhiE( ptMC[iPartMC], etaMC[iPartMC], phiMC[iPartMC], eMC[iPartMC] );

       Int_t pdgId = pdgIdMC[iPartMC];
     
       Float_t deltaRMC = jet.DeltaR(parton);

       bool goodPdgId = ( (fabs(pdgId)<=9) || (fabs(pdgId)==21) );
       if( !goodPdgId ) continue;
     
       if( (deltaRMC < deltaRMCmin) && goodPdgId ) {
         deltaRMCmin = deltaRMC;
         pdgIdPart_found = pdgIdMC[iPartMC];
         pdgIdMomPart_found = pdgIdMC[motherIDMC[iPartMC]];
       }

     } //for MC particles


     partPdgIDjet2 = ( deltaRMCmin<0.5 ) ? pdgIdPart_found : -999; 
     partMomPdgIDjet2 = ( deltaRMCmin<0.5 ) ? pdgIdMomPart_found : -999; 


	}else{
	  ptjet2 = -999;
	  ptcorrjet2 = -999;
	  etajet2 = -999;	 
	  phijet2 = -999;	 
	  betajet2 = -999.;
	  betastarjet2 = -999.;
	  assjet2 = -999.;
 	  btagvtxjet2 = -999.;
 	  btagtrkjet2 = -999.;
 	  ptDjet2 = -999.;
 	  rmsjet2 = -999.;
	  ntrkjet2 = -999.;
 	  nneutjet2 = -999.; 
 	  jetIdSimple_mvajet2 = -999.; 
 	  jetIdFull_mvajet2   = -999.; 
 	  jetId_dR2Meanjet2   = -999.; 
 	  jetId_betaStarClassicjet2 = -999.; 
 	  jetIdCutBased_wpjet2 = -999.; 
 	  jetIdSimple_wpjet2 = -999.; 
 	  jetIdFull_wpjet2 = -999.; 
	  jetId_frac01jet2 = -999.; 
 	  jetId_frac02jet2 = -999.; 
 	  jetId_frac03jet2 = -999.; 
 	  jetId_frac04jet2 = -999.; 
 	  jetId_frac05jet2 = -999.; 
 	  jetId_betajet2 = -999.; 
 	  jetId_betaStarjet2 = -999.; 
        partPdgIDjet2 = -999;
        partMomPdgIDjet2 = -999;
	}
	if( firsttennoisojet.at(2) > -1) {
	  ptjet3 = ptJet_pfakt5[firsttennoisojet.at(2)];
	  ptcorrjet3 = ptCorrJet_pfakt5[firsttennoisojet.at(2)];	  
	  ecorrjet3 = ptcorrjet3/ptjet3*eJet_pfakt5[firsttennoisojet.at(2)];
	  etajet3 = etaJet_pfakt5[firsttennoisojet.at(2)];
	  phijet3 = phiJet_pfakt5[firsttennoisojet.at(2)];
 	  btagvtxjet3 = simpleSecondaryVertexHighEffBJetTags[firsttennoisojet.at(2)];
 	  btagtrkjet3 = trackCountingHighEffBJetTags[firsttennoisojet.at(2)];	  
 	  btagjprobjet3 = jetProbabilityBJetTags[firsttennoisojet.at(2)];	  
 	  ptDjet3 = ptDJet_pfakt5[firsttennoisojet.at(2)];
	  rmsjet3 = rmsCandJet_pfakt5[firsttennoisojet.at(2)];
 	  ntrkjet3 = nChargedHadrons_pfakt5[firsttennoisojet.at(2)];
 	  nneutjet3 = nPhotons_pfakt5[firsttennoisojet.at(2)] + nNeutralHadrons_pfakt5[firsttennoisojet.at(2)] + nHFHadrons_pfakt5[firsttennoisojet.at(2)] + nHFEM_pfakt5[firsttennoisojet.at(2)];
     
     
     // match to parton
     Float_t deltaRMCmin = 999.;
     Int_t pdgIdPart_found = 0;
     Int_t pdgIdMomPart_found = 0;

     for(Int_t iPartMC=0; iPartMC<nMC; ++iPartMC) {

       if( statusMC[iPartMC]!=3 ) continue;

       if( ptMC[iPartMC]<0.1 ) continue;

       TLorentzVector jet;
       jet.SetPtEtaPhiE( ptcorrjet3, etajet3, phijet3, ecorrjet3 );
       TLorentzVector parton;
       parton.SetPtEtaPhiE( ptMC[iPartMC], etaMC[iPartMC], phiMC[iPartMC], eMC[iPartMC] );

       Int_t pdgId = pdgIdMC[iPartMC];
     
       Float_t deltaRMC = jet.DeltaR(parton);

       bool goodPdgId = ( (fabs(pdgId)<=9) || (fabs(pdgId)==21) );
       if( !goodPdgId ) continue;
     
       if( (deltaRMC < deltaRMCmin) && goodPdgId ) {
         deltaRMCmin = deltaRMC;
         pdgIdPart_found = pdgIdMC[iPartMC];
         pdgIdMomPart_found = pdgIdMC[motherIDMC[iPartMC]];
       }

     } //for MC particles


     partPdgIDjet3 = ( deltaRMCmin<0.5 ) ? pdgIdPart_found : -999; 
     partMomPdgIDjet3 = ( deltaRMCmin<0.5 ) ? pdgIdMomPart_found : -999; 



	}else{
	  ptjet3 = -999;
	  ptcorrjet3 = -999;
	  etajet3 = -999;	 
	  phijet3 = -999;	 
 	  btagvtxjet3 = -999;
 	  btagtrkjet3 = -999;
 	  btagjprobjet3 = -999;
 	  ptDjet3 = -999;
	  rmsjet3 = -999;
 	  ntrkjet3 = -999;
 	  nneutjet3 = -999;
        partPdgIDjet3 = -999;
        partMomPdgIDjet3 = -999;
	}
	if( firsttennoisojet.at(3) > -1) {
	  ptjet4 = ptJet_pfakt5[firsttennoisojet.at(3)];
	  ptcorrjet4 = ptCorrJet_pfakt5[firsttennoisojet.at(3)];	  
	  ecorrjet4 = ptcorrjet4/ptjet4*eJet_pfakt5[firsttennoisojet.at(3)];
	  etajet4 = etaJet_pfakt5[firsttennoisojet.at(3)];
	  phijet4 = phiJet_pfakt5[firsttennoisojet.at(3)];
 	  btagvtxjet4 = simpleSecondaryVertexHighEffBJetTags[firsttennoisojet.at(3)];
 	  btagtrkjet4 = trackCountingHighEffBJetTags[firsttennoisojet.at(3)];	  
 	  btagjprobjet4 = jetProbabilityBJetTags[firsttennoisojet.at(3)];	  
 	  ptDjet4 = ptDJet_pfakt5[firsttennoisojet.at(3)];
	  rmsjet4 = rmsCandJet_pfakt5[firsttennoisojet.at(3)];
 	  ntrkjet4 = nChargedHadrons_pfakt5[firsttennoisojet.at(3)];
 	  nneutjet4 = nPhotons_pfakt5[firsttennoisojet.at(3)] + nNeutralHadrons_pfakt5[firsttennoisojet.at(3)] + nHFHadrons_pfakt5[firsttennoisojet.at(3)] + nHFEM_pfakt5[firsttennoisojet.at(3)];
     
     
     // match to parton
     Float_t deltaRMCmin = 999.;
     Int_t pdgIdPart_found = 0;
     Int_t pdgIdMomPart_found = 0;

     for(Int_t iPartMC=0; iPartMC<nMC; ++iPartMC) {

       if( statusMC[iPartMC]!=3 ) continue;

       if( ptMC[iPartMC]<0.1 ) continue;

       TLorentzVector jet;
       jet.SetPtEtaPhiE( ptcorrjet4, etajet4, phijet4, ecorrjet4 );
       TLorentzVector parton;
       parton.SetPtEtaPhiE( ptMC[iPartMC], etaMC[iPartMC], phiMC[iPartMC], eMC[iPartMC] );

       Int_t pdgId = pdgIdMC[iPartMC];
     
       Float_t deltaRMC = jet.DeltaR(parton);

       bool goodPdgId = ( (fabs(pdgId)<=9) || (fabs(pdgId)==21) );
       if( !goodPdgId ) continue;
     
       if( (deltaRMC < deltaRMCmin) && goodPdgId ) {
         deltaRMCmin = deltaRMC;
         pdgIdPart_found = pdgIdMC[iPartMC];
         pdgIdMomPart_found = pdgIdMC[motherIDMC[iPartMC]];
       }

     } //for MC particles


     partPdgIDjet4 = ( deltaRMCmin<0.5 ) ? pdgIdPart_found : -999; 
     partMomPdgIDjet4 = ( deltaRMCmin<0.5 ) ? pdgIdMomPart_found : -999; 



	}else{
	  ptjet4 = -999;
	  ptcorrjet4 = -999;
	  etajet4 = -999;	 
	  phijet4 = -999;	 
 	  btagvtxjet4 = -999;
 	  btagtrkjet4 = -999;
 	  btagjprobjet4 = -999;
 	  ptDjet4 = -999;
	  rmsjet4 = -999;
 	  ntrkjet4 = -999;
 	  nneutjet4 = -999;
        partPdgIDjet4 = -999;
        partMomPdgIDjet4 = -999;
	}
	if( firsttennoisojet.at(0) > -1 && firsttennoisojet.at(1) > -1) {
	  deltaeta = etaJet_pfakt5[firsttennoisojet.at(0)]-etaJet_pfakt5[firsttennoisojet.at(1)];
	  double aveeta = (etaJet_pfakt5[firsttennoisojet.at(0)]+etaJet_pfakt5[firsttennoisojet.at(1)])/2;
	  zeppenjet = etahiggsiso - aveeta;
	  invmassjet = twojetsmassiso;
	  eta2j = etatwojetsiso;
	  phi2j = phitwojetsiso;
	  pt2j = pttwojetsiso;
	  deltaphi = delta_phi(phihiggsiso,phitwojetsiso);
	  deltaphinewvtx = delta_phi(phihiggsisonewvtx,phitwojetsiso);	  
	}else{
	  deltaeta = -999.;
	  zeppenjet = -999.;
	  invmassjet = -999.;
	  eta2j = -999.;
	  phi2j = -999.;
	  pt2j = -999.;
	  deltaphi = -999.;
	  deltaphinewvtx = -999.;
	}	  
*/


    ///////////////////////////////////////////////////
	// met = epfMet;
	// phimet = phipfMet;

    sMet_ = sMet;
    eMet_ = eMet;
    phiMet_ = phiMet;
    TLorentzVector tlvPFmet;
    tlvPFmet.SetPtEtaPhiE(epfMet,0,phipfMet,epfMet);
    TLorentzVector theSmearedMet = correctMet(tlvPFmet);
    TLorentzVector theSmearedMetPUcorr = correctMet(tlvPFmet,1,0,1);
    TLorentzVector theShiftedMet = shiftMet(tlvPFmet);
    TLorentzVector theShiftedScaledMet = correctMet(theShiftedMet,0,1);
    TLorentzVector theShiftedScaledMetPUcorr = correctMet(theShiftedMet,0,1,1);
    TLorentzVector theSmearedShiftedMet = shiftMet(theSmearedMet);
    TLorentzVector theSmearedShiftedMetPUcorr = shiftMet(theSmearedMetPUcorr);
    eSmearedMet_   = theSmearedMet.Pt();
    phiSmearedMet_ = theSmearedMet.Phi();
    eShiftedMet_   = theShiftedMet.Pt();
    phiShiftedMet_ = theShiftedMet.Phi();
    eShiftedScaledMet_   = theShiftedScaledMet.Pt();
    phiShiftedScaledMet_ = theShiftedScaledMet.Phi();
    eSmearedShiftedMet_   = theSmearedShiftedMet.Pt();
    phiSmearedShiftedMet_ = theSmearedShiftedMet.Phi();
    eShiftedScaledMetPUcorr_   = theShiftedScaledMetPUcorr.Pt();
    phiShiftedScaledMetPUcorr_ = theShiftedScaledMetPUcorr.Phi();
    eSmearedShiftedMetPUcorr_   = theSmearedShiftedMetPUcorr.Pt();
    phiSmearedShiftedMetPUcorr_ = theSmearedShiftedMetPUcorr.Phi();
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

	LOGamma  = countLOGenGamma();
	ISRGamma = countISRGenGamma();
	FSRGamma = countFSRGenGamma();
	promptGamma = LOGamma + ISRGamma + FSRGamma;
	
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
	pt2g2j = twog2j.Pt();
	

	if(ptPhot[firstfourisophot.at(0)] > ptphot1cut && ptPhot[firstfourisophot.at(1)] > ptphot2cut){
	  higgsmassjustisocutreco.Fill(higgsisomass,weight);
	  higgsmassjustisocutrecofull.Fill(higgsisomass,weight);
	}

	if( ptCorrJet_pfakt5[firsttennoisojet.at(0)] > ptjet1cut && ptCorrJet_pfakt5[firsttennoisojet.at(1)] > ptjet2cut 
	    && ptPhot[firstfourisophot.at(0)] > ptphot1cut && ptPhot[firstfourisophot.at(1)] > ptphot2cut){
	  higgsmassisojetptcutreco.Fill(higgsisomass,weight);
	  higgsmassisojetptcutrecofull.Fill(higgsisomass,weight);
	  deltaetajetreco.Fill(etaJet_pfakt5[firsttennoisojet.at(0)]-etaJet_pfakt5[firsttennoisojet.at(1)],weight);
	  //	  double zeppen = higgsreco_pt - aveeta;
	  if(TMath::Abs(etaJet_pfakt5[firsttennoisojet.at(0)]-etaJet_pfakt5[firsttennoisojet.at(1)])>deltaetacut){
	    higgsmassisocutreco.Fill(higgsisomass,weight);
	    higgsmassisocutrecofull.Fill(higgsisomass,weight);
	    double aveeta = (etaJet_pfakt5[firsttennoisojet.at(0)]+etaJet_pfakt5[firsttennoisojet.at(1)])/2;
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
  double cut, cutptlow, cutpthigh;
  for(int j=0; j<nJetGen_akt5; j++){	
    double DR = sqrt(delta_eta(etaJet_pfakt5[i],etaJetGen_akt5[j])*delta_eta(etaJet_pfakt5[i],etaJetGen_akt5[j]) + 
		     delta_phi(phiJet_pfakt5[i],phiJetGen_akt5[j])*delta_phi(phiJet_pfakt5[i],phiJetGen_akt5[j]) ) ;
    //    if(DR < .1 && TMath::Abs(ptCorrJet_pfakt5[i]-ptJetGen_akt5[j])/ptJetGen_akt5[j]  < 0.5) ass = 1; 
    if(TMath::Abs(etaJetGen_akt5[j])<2.5) {cut = 0.1; cutptlow = -0.75; cutpthigh = 1.;}
    else if(TMath::Abs(etaJetGen_akt5[j])<3.5) {cut = 0.15; cutptlow = -0.75; cutpthigh = 1.;}
    else if(TMath::Abs(etaJetGen_akt5[j])<4.5) {cut = 0.12; cutptlow = -0.75; cutpthigh = 1.;}
    else {cut = 0.2; cutptlow = -0.75; cutpthigh = 1.;}
    if( (DR < cut + 0.3 * exp(-0.015*(ptJetGen_akt5[j]-10))) &&  
	((ptCorrJet_pfakt5[i]-ptJetGen_akt5[j])/ptJetGen_akt5[j]  > cutptlow)  && 
        ((ptCorrJet_pfakt5[i]-ptJetGen_akt5[j])/ptJetGen_akt5[j]  < cutpthigh)
	)  ass = 1;
  }

  return ass;

}

// pfjet resolutions. taken from AN-2010-371
double RedNtpTree::ErrEt( double Et, double Eta) {
  
  double InvPerr2;

  double N, S, C, m;
  if(fabs(Eta) < 0.5 ) {
    N = 3.96859;
    S = 0.18348;
    C = 0.;
    m = 0.62627;
  } else if( fabs(Eta) < 1. ) {
    N = 3.55226;
    S = 0.24026;
    C = 0.;
    m = 0.52571;
  } else if( fabs(Eta) < 1.5 ) {
    N = 4.54826;
    S = 0.22652;
    C = 0.;
    m = 0.58963;
  } else if( fabs(Eta) < 2. ) {
    N = 4.62622;
    S = 0.23664;
    C = 0.;
    m = 0.48738;
  } else if( fabs(Eta) < 3. ) {
    N = 2.53324;
    S = 0.34306;
    C = 0.;
    m = 0.28662;
//   } else if( fabs(Eta) < 3. ) {
//     N = -3.33814;
//     S = 0.73360;
//     C = 0.;
//     m = 0.08264;
  } else if( fabs(Eta) < 5. ) {
    N = 2.95397;
    S = 0.11619;
    C = 0.;
    m = 0.96086;
  }

  // this is the absolute resolution (squared), not sigma(pt)/pt
  // so have to multiply by pt^2, thats why m+1 instead of m-1
  InvPerr2 =  (N * fabs(N) ) + (S * S) * pow(Et, m+1) + (C * C) * Et * Et ;


  return sqrt(InvPerr2)/Et;

}

void RedNtpTree::correctJets(int shift, float smear)
{

  for(int i=0; i<nJet_pfakt5; i++){
    
    double increase_endcap = 1;
    //    if(TMath::Abs(etaJet_pfakt5[i])>2.5 && TMath::Abs(etaJet_pfakt5[i])<3.4) increase_endcap = 2;  
    ptCorrJet_pfakt5[i] *= 1 + increase_endcap * shift * jetsyst_->getJESUncertainty(etaJet_pfakt5[i],ptCorrJet_pfakt5[i]);
    ptJet_pfakt5[i] *= 1 + increase_endcap * shift * jetsyst_->getJESUncertainty(etaJet_pfakt5[i],ptCorrJet_pfakt5[i]);
    eJet_pfakt5[i] *= 1 + increase_endcap * shift * jetsyst_->getJESUncertainty(etaJet_pfakt5[i],eJet_pfakt5[i]);
    
    if(smear){
      
      int ass(-999);
      for(int j=0; j<nJetGen_akt5; j++){	
	double DR = sqrt(delta_eta(etaJet_pfakt5[i],etaJetGen_akt5[j])*delta_eta(etaJet_pfakt5[i],etaJetGen_akt5[j]) + 
			 delta_phi(phiJet_pfakt5[i],phiJetGen_akt5[j])*delta_phi(phiJet_pfakt5[i],phiJetGen_akt5[j]) ) ;
	if(DR < .1 + 0.3 * exp(-0.015*(ptJetGen_akt5[j]-10)) && TMath::Abs(ptCorrJet_pfakt5[i]-ptJetGen_akt5[j])/ptJetGen_akt5[j]  < 0.5) ass = j; 
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

TLorentzVector RedNtpTree::correctMet(TLorentzVector uncormet, bool smearing, bool scale, bool PUremoval) {
  
  TLorentzVector jetSumSmeared;
  jetSumSmeared.SetXYZT(0.,0.,0.,0);
  
  TLorentzVector jetSumUnsmeared;
  jetSumUnsmeared.SetXYZT(0.,0.,0.,0);

  // associating reco - gen met                                                                                                            
  for(int i=0; i<nJet_pfakt5; i++){
    
    // remove identified photons
    if(!jetnoisophot.at(i)) continue;
    
    bool goodetajet(1);
    if(PUremoval){
      if(TMath::Abs(etaJet_pfakt5[i]) < 2.5) {
	if(betaStar_pfakt5[i][vrankPhotonPairs[0]] > 0.2 * log( nvertex - 0.67 ) ) goodetajet = 0;
	if(rmsCandJet_pfakt5[i] > 0.06) goodetajet = 0;
      } else if(TMath::Abs(etaJet_pfakt5[i]) < 3){
	if(rmsCandJet_pfakt5[i] > 0.05) goodetajet = 0;
      } else {
	if(rmsCandJet_pfakt5[i] > 0.055) goodetajet = 0;
      }
    }

    // smearing via association with genjets
    int ass(-999);
    double DRmin(999.);
    for(int j=0; j<nJetGen_akt5; j++){
      double DR = sqrt(delta_eta(etaJet_pfakt5[i],etaJetGen_akt5[j])*delta_eta(etaJet_pfakt5[i],etaJetGen_akt5[j]) +
		       delta_phi(phiJet_pfakt5[i],phiJetGen_akt5[j])*delta_phi(phiJet_pfakt5[i],phiJetGen_akt5[j]) ) ;
      double expres = ErrEt(ptCorrJet_pfakt5[i],etaJet_pfakt5[i]);
      if(DR < DRmin && (ptCorrJet_pfakt5[i]-ptJetGen_akt5[j])/ptCorrJet_pfakt5[i] < 5. * expres) {
	ass = j;
	DRmin = DR;
      }
    }
    
    if(DRmin > 0.1 + 0.3 * exp(-0.05*(ptJetGen_akt5[ass]-10)))  ass = -999;

//     if(ass>-1) jetDR->Fill(ptJetGen_akt5[ass],DRmin);
//     if(ass>-1) jetresp->Fill(ptJetGen_akt5[ass],(ptCorrJet_pfakt5[i]-ptJetGen_akt5[ass])/ptJetGen_akt5[ass]);
    
    // smearing for non-associated jets, using expected resolutions
    float smear = -999.;
    if (fabs(etaJet_pfakt5[i])<=1.1)                               smear = 1.06177;
    if (fabs(etaJet_pfakt5[i])<=1.7 && fabs(etaJet_pfakt5[i])>1.1) smear = 1.08352;
    if (fabs(etaJet_pfakt5[i])<=2.3 && fabs(etaJet_pfakt5[i])>1.7) smear = 1.02911;
    if (fabs(etaJet_pfakt5[i])>2.3)                                smear = 1.15288;
    
    double shift(0);
    if(ass>-1)
      shift = (smear-1) * (ptCorrJet_pfakt5[i] - ptJetGen_akt5[ass])/ptCorrJet_pfakt5[i];
    else {
      double expres = ErrEt(ptJet_pfakt5[i],etaJet_pfakt5[i]);
      double relsmear = expres * sqrt(smear*smear-1);
      TRandom3 gen(int(eventRN+ptJet_pfakt5[i]*1000));
      shift = gen.Gaus(0.,relsmear);
    }

    float ptSmeared  = ptJet_pfakt5[i];
    float eneSmeared = eJet_pfakt5[i];

    if(smearing && shift>-1 && shift < 2) {
      ptSmeared  *= 1 + shift;
      eneSmeared *= 1 + shift;
    }

    // JEC scaling to correct for residual jet corrections
    if(scale) {
      double factor(1);
      if(TMath::Abs(etaJet_pfakt5[i])<1.5) factor = 1.015;
      else if(TMath::Abs(etaJet_pfakt5[i])<3) factor = 1.04;
      else factor = 1.15;
      ptSmeared  *= factor;
      eneSmeared *= factor;
    }

    TLorentzVector thisJetSmeared;
    thisJetSmeared.SetPtEtaPhiE(ptSmeared,etaJet_pfakt5[i],phiJet_pfakt5[i],eneSmeared);
    
    TLorentzVector thisJetUnsmeared;

    thisJetUnsmeared.SetPtEtaPhiE(ptJet_pfakt5[i],etaJet_pfakt5[i],phiJet_pfakt5[i],eJet_pfakt5[i]);
    
    //    if(!PUremoval || ptJet_pfakt5[i]>50) ass=0;
    if (ptJet_pfakt5[i]>10 && TMath::Abs(etaJet_pfakt5[i])<4.7) {
      //      if(ass>-1) jetSumSmeared   += thisJetSmeared;
      jetSumSmeared   += thisJetSmeared;
      jetSumUnsmeared += thisJetUnsmeared;
    }

  }

  TLorentzVector correctedMet;
  correctedMet = uncormet + jetSumUnsmeared - jetSumSmeared;

  return correctedMet;
}

TLorentzVector RedNtpTree::shiftMet(TLorentzVector uncormet) {

  TLorentzVector correctedMet;
  
  // correction for METx, METy bias
  double px(0), py(0), e(0);

#ifdef MET_SHIFT_2012
  if(nMC==0){
    px = uncormet.Pt()*cos(uncormet.Phi())-0.006239*spfMet+0.662;
    py = uncormet.Pt()*sin(uncormet.Phi())+0.004613*spfMet-0.673;
  // MC
  }else{
    px = uncormet.Pt()*cos(uncormet.Phi())+0.00135*spfMet-0.021;
    py = uncormet.Pt()*sin(uncormet.Phi())+0.00371*spfMet-0.826;
   }
#else
// data
   if(nMC==0){
     px = uncormet.Pt()*cos(uncormet.Phi())-0.00563109*spfMet+0.959742;
     py = uncormet.Pt()*sin(uncormet.Phi())+0.00586162*spfMet-0.540137;
   // MC
   }else{
     px = uncormet.Pt()*cos(uncormet.Phi())-0.00069992*spfMet+0.430059;
     py = uncormet.Pt()*sin(uncormet.Phi())+0.00262869*spfMet+0.210784;
   }
#endif

  e = sqrt(px*px+py*py);
  
  correctedMet.SetPxPyPzE(px,py,0,e);

  return correctedMet;
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
      smearEnePhot[iPho]=smearing;
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

    gen_pt_lep1  = -999.;
    gen_pt_lep2  = -999.;
    gen_eta_lep1 = -999.;
    gen_eta_lep2 = -999.;
    gen_phi_lep1 = -999.;
    gen_phi_lep2 = -999.;
    gen_pid_lep1 = -999;
    gen_pid_lep2 = -999;
}


void RedNtpTree::SetAllRecoVarToMinus999()
{
massgg = -999;
ptgg = -999;
massggnewvtx = -999;
ptggnewvtx = -999;
ptphot1 = -999;
ptphot2 = -999;
deltaRToTrackphot1 = -999;
deltaRToTrackphot2 = -999;
etaphot1 = -999;
etaphot2 = -999;
phiphot1 = -999;
phiphot2 = -999;
// timephot1 = -999; 
// timephot2 = -999; 
E1phot1 = -999;
E1phot2 = -999;
E9phot1 = -999;
E9phot2 = -999;
energyErrphot1 = -999;
energyErrphot2 = -999;
energySmearingphot1 = -999;
energySmearingphot2 = -999;
for( unsigned i=0; i<10; ++i ) {
  ptjet[i] = -999;
  ptcorrjet[i] = -999;
  etajet[i] = -999;
  phijet[i] = -999;
  betajet[i] = -999;
  betastarjet[i] = -999;
  assjet[i] = -999;
}
deltaeta = -999;
zeppenjet = -999;
deltaphi = -999;
deltaphinewvtx = -999;
deltaphigg = -999;
invmassjet = -999;
invmass2g1j = -999;
invmass2g2j = -999;
 nvtx = -999;
 vtxId = -999;
 vtxPos_x = -999;
 vtxPos_y = -999;
 vtxPos_z = -999;
 vtxIdMVA = -999;
 vtxIdEvtProb = -999;
   //////////////////////////////////////
    sMet_   = -999;
    eMet_   = -999;
    phiMet_ = -999;
    signifMet_ = -999;
    eSmearedMet_ = -999;
    phiSmearedMet_ = -999;
    eShiftedMet_ = -999;
    phiShiftedMet_ = -999;
    eShiftedScaledMet_ = -999;
    phiShiftedScaledMet_ = -999;
    eSmearedShiftedMet_ = -999;
    phiSmearedShiftedMet_ = -999;
    eShiftedScaledMetPUcorr_ = -999;
    phiShiftedScaledMetPUcorr_ = -999;
    eSmearedShiftedMetPUcorr_ = -999;
    phiSmearedShiftedMetPUcorr_ = -999;
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
//     ecleanedSaclayPfMet_ = -999;   //[nvertex]
//     phicleanedSaclayPfMet_ = -999;   //[nvertex]
//     signifcleanedSaclayPfMet_ = -999;   //[nvertex]
//     eminTypeICleanSaclayPfMet_ = -999;   //[nvertex]
//     phiminTypeICleanSaclayPfMet_ = -999;   //[nvertex]
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

   ptele1    = -999.;
   ptele2    = -999.;
   etaele1   = -999.;
   etaele2   = -999.;
   phiele1   = -999.;
   phiele2   = -999.;
   eneele1   = -999.;
   eneele2   = -999.;
   sIeIeele1 = -999.;
   sIeIeele2 = -999.;
   dphiele1  = -999.;
   dphiele2  = -999.;
   detaele1  = -999.;
   detaele2  = -999.;
   mhitsele1 = -999;
   mhitsele2 = -999;
   d0ele1    = -999.;
   d0ele2    = -999.;
   dzele1    = -999.;
   dzele2    = -999.;
   invMassele1g1 = -999.;
   invMassele1g2 = -999.;
   invMassele2g1 = -999.;
   invMassele2g2 = -999.;

   if (LEPTONS_2012) {
     oEmoPele1      = -999.;
     oEmoPele2      = -999.;
     mvanotrigele1  = -999.;
     mvanotrigele2  = -999.;
     mvatrigele1    = -999.;
     mvatrigele2    = -999.;
     matchconvele1  = -999;
     matchconvele2  = -999;
     chHadIso03ele1 = -999.;
     chHadIso03ele2 = -999.;
     nHadIso03ele1  = -999.;
     nHadIso03ele2  = -999.;
     photIso03ele1  = -999.;
     photIso03ele2  = -999.;
   }
   if (LEPTONS_2011) {
     dcotele1  = -999.;
     dcotele2  = -999.;
     distele1  = -999.;
     distele2  = -999.;
     isoele1   = -999.;
     isoele2   = -999.;
     fullisoele1 = -999.;
     fullisoele2 = -999.;
   }

   pteleloose1    = -999.;
   pteleloose2    = -999.;
   etaeleloose1   = -999.;
   etaeleloose2   = -999.;
   phieleloose1   = -999.;
   phieleloose2   = -999.;
   eneeleloose1   = -999.;
   eneeleloose2   = -999.;
   sIeIeeleloose1 = -999.;
   sIeIeeleloose2 = -999.;
   dphieleloose1  = -999.;
   dphieleloose2  = -999.;
   detaeleloose1  = -999.;
   detaeleloose2  = -999.;
   mhitseleloose1 = -999;
   mhitseleloose2 = -999;
   d0eleloose1    = -999.;
   d0eleloose2    = -999.;
   dzeleloose1    = -999.;
   dzeleloose2    = -999.;
   invMasseleloose1g1 = -999.;
   invMasseleloose1g2 = -999.;
   invMasseleloose2g1 = -999.;
   invMasseleloose2g2 = -999.;
   if (LEPTONS_2012) {
     oEmoPeleloose1      = -999.;
     oEmoPeleloose2      = -999.;
     mvanotrigeleloose1  = -999.;
     mvanotrigeleloose2  = -999.;
     mvatrigeleloose1    = -999.;
     mvatrigeleloose2    = -999.;
     matchconveleloose1  = -999;
     matchconveleloose2  = -999;
     chHadIso03eleloose1 = -999.;
     chHadIso03eleloose2 = -999.;
     nHadIso03eleloose1  = -999.;
     nHadIso03eleloose2  = -999.;
     photIso03eleloose1  = -999.;
     photIso03eleloose2  = -999.;
   }
   if (LEPTONS_2011) {
     dcoteleloose1  = -999.;
     dcoteleloose2  = -999.;
     disteleloose1  = -999.;
     disteleloose2  = -999.;
     isoeleloose1   = -999.;
     isoeleloose2   = -999.;
     fullisoeleloose1 = -999.;
     fullisoeleloose2 = -999.;
   }

   ptmu1     = -999.;
   ptmu2     = -999.;
   etamu1    = -999.;
   etamu2    = -999.;
   phimu1    = -999.;
   phimu2    = -999.;
   enemu1    = -999.;
   enemu2    = -999.;
   pixhitsmu1 = -999;
   pixhitsmu2 = -999;
   trkhitsmu1 = -999;
   trkhitsmu2 = -999;
   hitsmu1 = -999;
   hitsmu2 = -999;
   chi2mu1 = -999.;   
   chi2mu2 = -999.;   
   matchmu1 = -999;
   matchmu2 = -999;
   d0mu1 = -999.;
   d0mu2 = -999.;
   dzmu1 = -999.;
   dzmu2 = -999.;

   if(LEPTONS_2012) {
     chHadmu1 = -999.;
     chHadmu2 = -999.;
     nHadmu1  = -999.;
     nHadmu2  = -999.;
     photmu1  = -999.;
     photmu2  = -999.;
     puptmu1  = -999.;
     puptmu2  = -999.;
   }
   if (LEPTONS_2011) {
     isomu1 = -999.;
     isomu2 = -999.;
   }

   ptmuloose1      = -999.;
   ptmuloose2      = -999.;
   etamuloose1     = -999.;
   etamuloose2     = -999.;
   phimuloose1     = -999.;
   phimuloose2     = -999.;
   enemuloose1     = -999.;
   enemuloose2     = -999.;
   pixhitsmuloose1 = -999;
   pixhitsmuloose2 = -999;
   trkhitsmuloose1 = -999;
   trkhitsmuloose2 = -999;
   hitsmuloose1  = -999;
   hitsmuloose2  = -999;
   chi2muloose1  = -999.;   
   chi2muloose2  = -999.;   
   matchmuloose1 = -999;
   matchmuloose2 = -999;
   d0muloose1 = -999.;
   d0muloose2 = -999.;
   dzmuloose1 = -999.;
   dzmuloose2 = -999.;
   if(LEPTONS_2012) {
     chHadmuloose1 = -999.;
     chHadmuloose2 = -999.;
     nHadmuloose1  = -999.;
     nHadmuloose2  = -999.;
     photmuloose1  = -999.;
     photmuloose2  = -999.;
     puptmuloose1  = -999.;
     puptmuloose2  = -999.;
   }
   if (LEPTONS_2011) {
     isomuloose1 = -999.;
     isomuloose2 = -999.;
   }

   promptGamma = -999;
   LOGamma     = -999;
   ISRGamma    = -999;
   FSRGamma    = -999;
 
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

bool RedNtpTree::leptonCutsEle2011(int iEle, electronidcuts const& pid, vector<bool> *vpass) {

  bool pt, eta, crack;
  bool setaeta, deta, dphi;
  bool minhits, dconv;
  bool d0, dz;
  bool isol;

  // acceptance
  pt    = electron_pt[iEle] > pid.pt;      
  eta   = fabs(electron_sc_eta[iEle]) < pid.eta; 
  crack = fabs(electron_sc_eta[iEle]) < pid.crack1 || fabs(electron_sc_eta[iEle]) > pid.crack2;

  // electronId + conv.rejection + impact parameters + isolation
  float d0Ele = eleDxyPV(iEle,vrankPhotonPairs[0]);   
  float dzEle = eleDzPV(iEle,vrankPhotonPairs[0]);   

  float fullHcal = electron_hcalIso03[iEle] + fabs(electron_HoE[iEle]*electron_sc_energy[iEle]/cosh(electron_sc_eta[iEle]));
  float electronIsoEB = electron_trkIso03[iEle] + std::max(0.,(electron_ecalIso03[iEle]-1.)) + fullHcal - rhoPF*TMath::Pi()*0.3*0.3; 
  float electronIsoEE = electron_trkIso03[iEle] + electron_ecalIso03[iEle] + fullHcal - rhoPF*TMath::Pi()*0.3*0.3; 
  
  if (fabs(electron_sc_eta[iEle])<1.4442) {
    setaeta = electron_SigmaIetaIeta[iEle] < pid.setaetaEB;
    deta    = fabs(electron_dEtaIn[iEle]) < pid.detaEB;
    dphi    = fabs(electron_dPhiIn[iEle]) < pid.dphiEB;  
    minhits = electron_misHits[iEle] <= pid.minhitsEB;
    dconv   = (fabs(electron_dcot[iEle]) > pid.dcotEB) || (fabs(electron_dist[iEle]) > pid.distEB);  
    d0      = fabs(d0Ele) < pid.d0EB;
    dz      = fabs(dzEle) < pid.dzEB;
    isol    = electronIsoEB < electron_pt[iEle]* pid.iso_relEB;
  } else {
    setaeta = electron_SigmaIetaIeta[iEle] < pid.setaetaEE;
    deta    = fabs(electron_dEtaIn[iEle]) < pid.detaEE;
    dphi    = fabs(electron_dPhiIn[iEle]) < pid.dphiEE;
    minhits = electron_misHits[iEle] <= pid.minhitsEE;
    dconv   = (fabs(electron_dcot[iEle]) > pid.dcotEE) || (fabs(electron_dist[iEle]) > pid.distEE); 
    d0      = fabs(d0Ele) < pid.d0EE;
    dz      = fabs(dzEle) < pid.dzEE;
    isol    = electronIsoEE < electron_pt[iEle]* pid.iso_relEE;
  }

  if (vpass) {
    if((*vpass).size()!=11) { cout << "major failure in LeptonCutsEle2011! (*vpass).size()!=11.. die!" << endl; exit(0) ; }
    (*vpass)[0]  = pt;
    (*vpass)[1]  = eta;
    (*vpass)[2]  = crack;
    (*vpass)[3]  = setaeta;
    (*vpass)[4]  = deta;
    (*vpass)[5]  = dphi;
    (*vpass)[6]  = minhits;
    (*vpass)[7]  = dconv;  
    (*vpass)[8]  = d0;
    (*vpass)[9]  = dz;
    (*vpass)[10] = isol;
  }

  return (pt && eta && crack && setaeta && deta && dphi && minhits && dconv && d0 && dz && isol);
}

bool RedNtpTree::leptonCutsEle2012(int iEle, electronidcuts2012 const& pid, vector<bool> *vpass) {

  bool pt, eta, crack;
  bool setaeta, deta, dphi, hoe, oeMop;
  bool minhits, matchconv;
  bool d0, dz;
  bool isol;

  // acceptance                                                                                                                            
  pt    = electron_pt[iEle] > pid.pt;
  eta   = fabs(electron_sc_eta[iEle]) < pid.eta;
  crack = fabs(electron_sc_eta[iEle]) < pid.crack1 || fabs(electron_sc_eta[iEle]) > pid.crack2;

  // impact parameters                                                                                                                     
  float d0Ele = eleDxyPV(iEle,vrankPhotonPairs[0]);
  float dzEle = eleDzPV(iEle,vrankPhotonPairs[0]);

  // effective areas                                                                                                                       
  float abseta = fabs(electron_sc_eta[iEle]);
  ElectronEffectiveArea::ElectronEffectiveAreaTarget effAreaTarget_   = ElectronEffectiveArea::kEleEAData2012;
  ElectronEffectiveArea::ElectronEffectiveAreaType effAreaGamma_      = ElectronEffectiveArea::kEleGammaIso03;
  ElectronEffectiveArea::ElectronEffectiveAreaType effAreaNeutralHad_ = ElectronEffectiveArea::kEleNeutralHadronIso03;
  ElectronEffectiveArea::ElectronEffectiveAreaType effAreaGammaAndNeutralHad_ = ElectronEffectiveArea::kEleGammaAndNeutralHadronIso03;
  float eff_area_ga   = ElectronEffectiveArea::GetElectronEffectiveArea(effAreaGamma_, abseta, effAreaTarget_);
  float eff_area_nh   = ElectronEffectiveArea::GetElectronEffectiveArea(effAreaNeutralHad_, abseta, effAreaTarget_);
  float eff_area_ganh = ElectronEffectiveArea::GetElectronEffectiveArea(effAreaGammaAndNeutralHad_, abseta, effAreaTarget_);
  float eff_area_sum  = eff_area_ga + eff_area_nh;
  // isolation                                                                                                                             
  float theIsolation = electron_chHad03Iso[iEle];
  theIsolation += max<float>(0.,electron_nHad03Iso[iEle]+electron_phot03Iso[iEle]-eff_area_sum*rhoAllJets);

  // full selection                                                                                                                        
  if (abseta<1.4442) {
    setaeta   = electron_SigmaIetaIeta[iEle] < pid.setaetaEB;
    deta      = fabs(electron_dEtaIn[iEle]) < pid.detaEB;
    dphi      = fabs(electron_dPhiIn[iEle]) < pid.dphiEB;
    hoe       = electron_HoE[iEle] < pid.hoeEB;
    oeMop     = fabs(1./electron_ecalEnergy[iEle] - 1./electron_trackPatVtx[iEle]) < pid.oemopEB ;
    d0        = fabs(d0Ele) < pid.d0EB;
    dz        = fabs(dzEle) < pid.dzEB;
    minhits   = electron_misHits[iEle] <= pid.minhitsEB;
    matchconv = electron_matchedConv[iEle]==0;
    isol      = theIsolation < electron_pt[iEle]* pid.iso_relEB;
  } else {
    setaeta = electron_SigmaIetaIeta[iEle] < pid.setaetaEE;
    deta    = fabs(electron_dEtaIn[iEle]) < pid.detaEE;
    dphi    = fabs(electron_dPhiIn[iEle]) < pid.dphiEE;
    hoe     = electron_HoE[iEle] < pid.hoeEE;
    oeMop   = fabs(1./electron_ecalEnergy[iEle] - 1./electron_trackPatVtx[iEle]) < pid.oemopEE;
    d0      = fabs(d0Ele) < pid.d0EE;
    dz      = fabs(dzEle) < pid.dzEE;
    minhits = electron_misHits[iEle] <= pid.minhitsEE;
    matchconv = electron_matchedConv[iEle]==0;
    isol    = theIsolation < electron_pt[iEle]* pid.iso_relEE;
  }

  if (vpass) {
    if((*vpass).size()!=13) { cout << "major failure in LeptonCutsEle2012! (*vpass).size()!=13.. die!" << endl; exit(0) ; }
    (*vpass)[0]  = pt;
    (*vpass)[1]  = eta;
    (*vpass)[2]  = crack;
    (*vpass)[3]  = setaeta;
    (*vpass)[4]  = deta;
    (*vpass)[5]  = dphi;
    (*vpass)[6]  = hoe;
    (*vpass)[7]  = oeMop;
    (*vpass)[8]  = d0;
    (*vpass)[9]  = dz;
    (*vpass)[10] = minhits;
    (*vpass)[11] = matchconv;
    (*vpass)[12] = isol;
  }

  return (pt && eta && crack && setaeta && deta && dphi && hoe && oeMop && d0 && dz && minhits && matchconv && isol);
}

bool RedNtpTree::leptonCutsMu2011(int iMu, muonidcuts const& pid, vector<bool> *vpass) {

  bool pt, eta;
  bool pixhits, tkhits, globalhits, chi2, match, globAndTrk;
  bool d0, dz;
  bool isol;

  // acceptance
  pt  = Muon_pt[iMu] > pid.pt;      
  eta = fabs(Muon_eta[iMu]) < pid.eta; 
	   
  // muonId 
  pixhits    = Muon_pixHits[iMu] > pid.pixhits;          
  tkhits     = Muon_tkHits[iMu] > pid.tkhits;            
                                                         
  globalhits = Muon_validHits[iMu] > pid.hits;           
  chi2       = Muon_normChi2[iMu] < pid.chi2;            
  match      = Muon_numberOfMatches[iMu] > pid.match;
  
  globAndTrk = Muon_isGlobalMuon[iMu] && Muon_isTrackerMuon[iMu];

  // impact parameter 
  float d0Muon = muonDxyPV(iMu,vrankPhotonPairs[0]);
  float dzMuon = muonDzPV(iMu,vrankPhotonPairs[0]);
  d0 = fabs(d0Muon) < pid.d0;
  dz = fabs(dzMuon) < pid.dz;
  
  // isolation 
  float muonIso  = Muon_trackIso[iMu] + Muon_ecalIso[iMu] + Muon_hcalIso[iMu] - rhoPF*TMath::Pi()*0.3*0.3;	   
  float relMuIso = muonIso/Muon_pt[iMu];
  isol = relMuIso < pid.iso_rel;

  if (vpass) {
    if((*vpass).size()!=11) { cout << "major failure! (*vpass).size()!=10.. die!" << endl; exit(0) ; }
    (*vpass)[0] = pt;
    (*vpass)[1] = eta;
    (*vpass)[2] = pixhits;
    (*vpass)[3] = tkhits;
    (*vpass)[4] = globalhits;
    (*vpass)[5] = chi2;
    (*vpass)[6] = match;
    (*vpass)[7] = d0;
    (*vpass)[8] = dz;
    (*vpass)[9] = isol;
    (*vpass)[10] = globAndTrk;
  }

  return (pt && eta && pixhits && tkhits && globalhits && chi2 && match && d0 && dz && isol && globAndTrk);
}

bool RedNtpTree::leptonCutsMu2012(int iMu, muonidcuts2012 const& pid, vector<bool> *vpass) {

  bool pt, eta;
  bool globalmu, pfmu;
  bool pixhits, globalhits, chi2, match, withm;
  bool d0, dz;
  bool isol;

  // acceptance                                                                                                                            
  pt  = Muon_pt[iMu] > pid.pt;
  eta = fabs(Muon_eta[iMu]) < pid.eta;

  // muonId                                                                                                                                
  globalmu   = Muon_isGlobalMuon[iMu];       
  pfmu       = Muon_isPFMuon[iMu];           
  chi2       = Muon_normChi2[iMu] < pid.chi2;
  globalhits = Muon_validHits[iMu] > pid.hits;
  match      = Muon_numberOfMatches[iMu] > pid.match;
  pixhits    = Muon_pixHits[iMu] > pid.pixhits;
  withm      = Muon_trkLayerWithMeas[iMu] > pid.withm;   

  // impact parameter                                                                                                                      
  float d0Muon = muonDxyPV(iMu,vrankPhotonPairs[0]);
  float dzMuon = muonDzPV(iMu,vrankPhotonPairs[0]);
  d0 = fabs(d0Muon) < pid.d0;
  dz = fabs(dzMuon) < pid.dz;

  // isolation                                                                                                                             
  float muonIso  = Muon_pfiso04_chHad[iMu] + max(0., Muon_pfiso04_nHad[iMu]+Muon_pfiso04_Phot[iMu]-0.5*Muon_pfiso04_PUPt[iMu]);
  float relMuIso = muonIso/Muon_pt[iMu];
  isol = relMuIso < pid.iso_rel;

  if (vpass) {
    if((*vpass).size()!=12) { cout << "major failure! (*vpass).size()!=12.. die!" << endl; exit(0) ; }
    (*vpass)[0]  = pt;
    (*vpass)[1]  = eta;
    (*vpass)[2]  = globalmu;
    (*vpass)[3]  = pfmu;
    (*vpass)[4]  = chi2;
    (*vpass)[5]  = globalhits;
    (*vpass)[6]  = match;
    (*vpass)[7]  = pixhits;
    (*vpass)[8]  = withm;
    (*vpass)[9]  = d0;
    (*vpass)[10] = dz;
    (*vpass)[11] = isol;
  }

  return (pt && eta && globalmu && pfmu && chi2 && globalhits && match && pixhits && withm && d0 && dz && isol);
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

double RedNtpTree::eleDzPV(int iele, int iPV) {
  TVector3 PVPos(vx[iPV],vy[iPV],vz[iPV]);
  TVector3 lepVPos(electron_vx[iele],electron_vy[iele],electron_vz[iele]);  
  TVector3 lepMom(electron_px[iele],electron_py[iele],electron_pz[iele]);
  return trackDzPV(PVPos,lepVPos,lepMom);
}

double RedNtpTree::eleDxyPV(int iele, int iPV) {
  TVector3 PVPos(vx[iPV],vy[iPV],vz[iPV]);
  TVector3 lepVPos(electron_vx[iele],electron_vy[iele],electron_vz[iele]);
  TVector3 lepMom(electron_px[iele],electron_py[iele],electron_pz[iele]);
  return trackDxyPV(PVPos,lepVPos,lepMom);
}

double RedNtpTree::muonDzPV(int imu, int iPV) {
  TVector3 PVPos(vx[iPV],vy[iPV],vz[iPV]);
  TVector3 lepVPos(Muon_vx[imu],Muon_vy[imu],Muon_vz[imu]);
  TVector3 lepMom(Muon_px[imu],Muon_py[imu],Muon_pz[imu]);
  return trackDzPV(PVPos,lepVPos,lepMom);
}

double RedNtpTree::muonDxyPV(int imu, int iPV) {
  TVector3 PVPos(vx[iPV],vy[iPV],vz[iPV]);
  TVector3 lepVPos(Muon_vx[imu],Muon_vy[imu],Muon_vz[imu]);
  TVector3 lepMom(Muon_px[imu],Muon_py[imu],Muon_pz[imu]);
  return trackDxyPV(PVPos,lepVPos,lepMom);
}

double RedNtpTree::trackDzPV(TVector3 PVPos, TVector3 trackVPos, TVector3 trackMom) {
  float trackPt = trackMom.Pt();
  return (trackVPos.Z()-PVPos.Z()) - ((trackVPos.X()-PVPos.X())*trackMom.X()+(trackVPos.Y()-PVPos.Y())*trackMom.Y())/trackPt *trackMom.Pz()/trackPt;
}

double RedNtpTree::trackDxyPV(TVector3 PVPos, TVector3 trackVPos, TVector3 trackMom) {
  return ( - (trackVPos.X()-PVPos.X())*trackMom.Y() + (trackVPos.Y()-PVPos.Y())*trackMom.X() ) / trackMom.Pt();
}

int RedNtpTree::countLOGenGamma(){
  
  int totLO = 0;
  for (int ii=0; ii<nMC; ii++) {
    int myStatus = statusMC[ii];
    int myId     = pdgIdMC[ii];
    if (myStatus==3 && myId==22) {
      int myMoth   = motherIDMC[ii];
      int myMothId = abs(pdgIdMC[myMoth]);
      if (myMothId<=25) totLO++;   // quarks, gluons, W, Z and ZHiggs as mothers                  
    }
  }
  return totLO;
}

int RedNtpTree::countISRGenGamma(){
  
  int totISR = 0;
  for (int ii=0; ii<nMC; ii++) {
    int myStatus = statusMC[ii];
    int myId     = pdgIdMC[ii];
    if (myStatus==1 && myId==22) {
      int myMoth   = motherIDMC[ii];
      int myMothId = abs(pdgIdMC[myMoth]);
      if (myMothId<11 || myMothId==21) totISR++;   // quarks and gluons as mothers                  
    }
  }
  return totISR;
}

int RedNtpTree::countFSRGenGamma(){
  
  int totFSR = 0;
  for (int ii=0; ii<nMC; ii++) {
    int myStatus = statusMC[ii];
    int myId     = pdgIdMC[ii];
    if (myStatus==1 && myId==22) {
      int myMoth   = motherIDMC[ii];
      int myMothId = abs(pdgIdMC[myMoth]);
      if (myMothId>10 && myMothId<21) totFSR++;   // leptons as mothers                  
    }
  }
  return totFSR;
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
void RedNtpTree::SetPhotonCutsInCategories(phoCiCIDLevel cutlevel, float * cic6_allcuts_lead, float * cic6_allcuts_sublead, 
					float * cic4_allcuts_lead, float * cic4_allcuts_sublead,
					float * cic4pf_allcuts_lead, float * cic4pf_allcuts_sublead) {

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
// phoHYPERTIGHT2 - sob value=0.00625        - iteration 6 - eff=0.699128  fake=0.00991561
// phoHYPERTIGHT3 - sob value=0.0125         - iteration 6 - eff=0.573171  fake=0.00520159
// phoHYPERTIGHT4 - sob value=0.025          - iteration 6 - eff=0.41176   fake=0.00217666


// 4 categories
// phoNOCUTS      - thresholds so all photons pass
// phoLOOSE       - sob value=0.0002         - iteration 8 - eff=0.939229  fake=0.0815158
// phoMEDIUM      - sob value=0.0004         - iteration 6 - eff=0.91754   fake=0.0581047
// phoTIGHT       - sob value=0.0008         - iteration 6 - eff=0.886869  fake=0.041063
// phoSUPERTIGHT  - sob value=0.0016         - iteration 6 - eff=0.844314  fake=0.0286033
// phoHYPERTIGHT1 - sob value=0.0032         - iteration 6 - eff=0.774552  fake=0.0191603
// phoHYPERTIGHT2 - sob value=0.00625        - iteration 6 - eff=0.67859   fake=0.0121262
// phoHYPERTIGHT3 - sob value=0.0125         - iteration 6 - eff=0.521328  fake=0.00626966
// phoHYPERTIGHT4 - sob value=0.025          - iteration 6 - eff=0.381192  fake=0.00357649



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

        float cic4pf_allcuts_temp_lead[] = { 
            1e+09,     1e+09,     1e+09,     1e+09,
            1e+09,     1e+09,     1e+09,     1e+09,
            1e+09,     1e+09,     1e+09,     1e+09,
            1e+09,     1e+09,     1e+09,     1e+09,
            1e+09,     1e+09,     1e+09,     1e+09,
            -1e+09,    -1e+09,    -1e+09,    -1e+09,
            -1e+09,    -1e+09,    -1e+09,    -1e+09,
            1.5,         1.5,         1.5,         1.5 };
        float cic4pf_allcuts_temp_sublead[] = { 
            1e+09,     1e+09,     1e+09,     1e+09,
            1e+09,     1e+09,     1e+09,     1e+09,
            1e+09,     1e+09,     1e+09,     1e+09,
            1e+09,     1e+09,     1e+09,     1e+09,
            1e+09,     1e+09,     1e+09,     1e+09,
            -1e+09,    -1e+09,    -1e+09,    -1e+09,
            -1e+09,    -1e+09,    -1e+09,    -1e+09,
            1.5,         1.5,         1.5,         1.5 };
        for(int i=0;i!=ncuts*ncat_cic4;++i) { 
            cic4pf_allcuts_lead[i]    = cic4pf_allcuts_temp_lead[i];
            cic4pf_allcuts_sublead[i] = cic4pf_allcuts_temp_sublead[i]; 
        }
          
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
        float cic4_allcuts_temp_lead[] = { 
            8.2,       4.1,       5.4,       2.6,
            67,        69,        85,       7.2,
            7.5,       4.5,       5.2,       2.5,
            0.0112,    0.0102,     0.029,     0.028,
            0.09,     0.089,     0.101,     0.073,
            0.94,      0.31,      0.92,      0.29,
            0.26,     0.029,    0.0062,    0.0055,
            1.5,         1.5,         1.5,         1.5};
        float cic4_allcuts_temp_sublead[] = { 
            8.2,       4.1,       5.4,       2.6,
            67,        69,        85,       7.2,
            7.5,       4.5,       5.2,       2.5,
            0.0112,    0.0102,     0.029,     0.028,
            0.09,     0.089,     0.101,     0.073,
            0.94,      0.31,      0.92,      0.29,
            0.26,     0.029,    0.0062,    0.0055,
            1.5,         1.5,         1.5,         1.5};
        for(int i=0;i!=ncuts*ncat_cic6;++i) { 
            cic6_allcuts_lead[i]=cic6_allcuts_temp_lead[i];
            cic6_allcuts_sublead[i]=cic6_allcuts_temp_sublead[i]; }
        for(int i=0;i!=ncuts*ncat_cic4;++i) { cic4_allcuts_lead[i]=cic4_allcuts_temp_lead[i];
                           cic4_allcuts_sublead[i]=cic4_allcuts_temp_sublead[i]; }

        float cic4pf_allcuts_temp_lead[] = { 
            8.9,       6.3,       9.8,       6.8,
            43,      19.4,        24,       7.9,
            6.2,       4.3,         5,       4.3,
            0.0117,    0.0105,     0.031,     0.031,
            0.137,      0.14,     0.145,     0.143,
            0.94,      0.25,      0.93,      0.24,
            1,    0.0136,    0.0138,    0.0122,
            1.5,         1.5,         1.5,         1.5};
        float cic4pf_allcuts_temp_sublead[] = { 
            8.9,       6.3,       9.8,       6.8,
            43,      19.4,        24,       7.9,
            6.2,       4.3,         5,       4.3,
            0.0117,    0.0105,     0.031,     0.031,
            0.137,      0.14,     0.145,     0.143,
            0.94,      0.25,      0.93,      0.24,
            1,    0.0136,    0.0138,    0.0122,
            1.5,         1.5,         1.5,         1.5};
        for(int i=0;i!=ncuts*ncat_cic4;++i) { 
            cic4pf_allcuts_lead[i]    = cic4pf_allcuts_temp_lead[i];
            cic4pf_allcuts_sublead[i] = cic4pf_allcuts_temp_sublead[i]; 
        }
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
        float cic4_allcuts_temp_lead[] = {  
            6.4,       3.2,       3.4,       2.2,
            64,      10.8,        13,       3.5,
            6.4,       3.4,       3.8,       2.1,
            0.0109,      0.01,     0.029,     0.028,
            0.089,     0.079,      0.09,     0.061,
            0.94,      0.32,      0.94,      0.29,
            0.98,     0.029,    0.0109,    0.0111,
            1.5,         1.5,         1.5,         1.5 };
        float cic4_allcuts_temp_sublead[] = {  
            6.4,       3.2,       3.4,       2.2,
            64,      10.8,        13,       3.5,
            6.4,       3.4,       3.8,       2.1,
            0.0109,      0.01,     0.029,     0.028,
            0.089,     0.079,      0.09,     0.061,
            0.94,      0.32,      0.94,      0.29,
            0.98,     0.029,    0.0109,    0.0111,
            1.5,         1.5,         1.5,         1.5 };
        for(int i=0;i!=ncuts*ncat_cic6;++i) { 
            cic6_allcuts_lead[i]=cic6_allcuts_temp_lead[i];
            cic6_allcuts_sublead[i]=cic6_allcuts_temp_sublead[i]; }
        for(int i=0;i!=ncuts*ncat_cic4;++i) { 
            cic4_allcuts_lead[i]=cic4_allcuts_temp_lead[i];
            cic4_allcuts_sublead[i]=cic4_allcuts_temp_sublead[i]; }

        float cic4pf_allcuts_temp_lead[] = { 
            8.4,       5.4,       6.1,       6.2,
            43,       8.6,       8.9,       6.1,
            5.5,       3.8,       3.8,       3.4,
            0.0116,    0.0104,     0.031,     0.029,
            0.137,     0.103,     0.145,     0.128,
            0.94,      0.25,      0.94,      0.24,
            1,     0.021,    0.0138,    0.0162, 
            1.5,         1.5,         1.5,         1.5};
        float cic4pf_allcuts_temp_sublead[] = { 
            8.4,       5.4,       6.1,       6.2,
            43,       8.6,       8.9,       6.1,
            5.5,       3.8,       3.8,       3.4,
            0.0116,    0.0104,     0.031,     0.029,
            0.137,     0.103,     0.145,     0.128,
            0.94,      0.25,      0.94,      0.24,
            1,     0.021,    0.0138,    0.0162, 
            1.5,         1.5,         1.5,         1.5};

        for(int i=0;i!=ncuts*ncat_cic4;++i) { 
            cic4pf_allcuts_lead[i]    = cic4pf_allcuts_temp_lead[i];
            cic4pf_allcuts_sublead[i] = cic4pf_allcuts_temp_sublead[i]; 
        }
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
        for(int i=0;i!=ncuts*ncat_cic6;++i) { 
            cic6_allcuts_lead[i]=cic6_allcuts_temp_lead[i];
            cic6_allcuts_sublead[i]=cic6_allcuts_temp_sublead[i]; }
        for(int i=0;i!=ncuts*ncat_cic4;++i) { 
            cic4_allcuts_lead[i]=cic4_allcuts_temp_lead[i];
            cic4_allcuts_sublead[i]=cic4_allcuts_temp_sublead[i]; }
           
        float cic4pf_allcuts_temp_lead[] = {
            6.6,       4.9,       5.6,       4.2,
            11.5,       7.2,       8.4,       5.5,
            4.4,       2.7,       3.3,       2.5,
            0.0111,    0.0103,     0.028,     0.028,
            0.124,     0.103,     0.142,     0.128,
            0.94,      0.25,      0.94,      0.24,
            1,     0.032,     0.024,    0.0173,
            1.5,         1.5,         1.5,         1.5};
        float cic4pf_allcuts_temp_sublead[] = { 
            6.6,       4.9,       5.6,       4.2,
            11.5,       7.2,       8.4,       5.5,
            4.4,       2.7,       3.3,       2.5,
            0.0111,    0.0103,     0.028,     0.028,
            0.124,     0.103,     0.142,     0.128,
            0.94,      0.25,      0.94,      0.24,
            1,     0.032,     0.024,    0.0173,
            1.5,         1.5,         1.5,         1.5};

        for(int i=0;i!=ncuts*ncat_cic4;++i) { 
            cic4pf_allcuts_lead[i]    = cic4pf_allcuts_temp_lead[i];
            cic4pf_allcuts_sublead[i] = cic4pf_allcuts_temp_sublead[i]; 
        }
           
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
        float cic4_allcuts_temp_lead[] = { 

                         // category from PhotonCategory(..)
                         // cat 0      cat 1    cat 2    cat 3 
                         // high R9    low R9   high R9  low R9
                         // barrel     barrel   endcap   endcap

            3.8,     2.2,     1.77,    1.29,   // isosumoet       sum of isolation cone energies divided by Et (?)   (upper cut, maximum value)
            11.7,    3.4,     3.9,     1.84,   // isosumoetbad    same as isosumoet but for 'bad' (worst ?) vertex   (upper cut, maximum value)
            3.5,     2.2,     2.3,     1.45,   // trkisooetom     tracker isolation cone only (divided by Et)        (upper cut, maximum value)
            0.0106,  0.0097,  0.028,   0.027,  // sieie           sigma(ieta,ieta)                                   (upper cut, maximum value)
            0.082,   0.062,   0.065,   0.048,  // hovere          H/E                                                (upper cut, maximum value)
            0.94,    0.36,    0.94,    0.32,   // r9              R9                                               (lower cut, minimum value) 
            1.,      0.062,   0.97,    0.97,   // drtotk_25_99    Delta R to track ?                               (lower cut, minimum value)
            1.5,     1.5,     1.5,     1.5 };  // pixel           has pixel seed ?                                   (upper cut, maximum value)
        float cic4_allcuts_temp_sublead[] = { 
            3.8,     2.2,     1.77,    1.29,
            11.7,    3.4,     3.9,     1.84,
            3.5,     2.2,     2.3,     1.45,
            0.0106,  0.0097,  0.028,   0.027,
            0.082,   0.062,   0.065,   0.048,
            0.94,    0.36,    0.94,    0.32,
            1.,      0.062,   0.97,    0.97,
            1.5,     1.5,     1.5,     1.5 };
        for(int i=0;i!=ncuts*ncat_cic6;++i) { 
            cic6_allcuts_lead[i]=cic6_allcuts_temp_lead[i];
            cic6_allcuts_sublead[i]=cic6_allcuts_temp_sublead[i]; }
        for(int i=0;i!=ncuts*ncat_cic4;++i) { 
            cic4_allcuts_lead[i]=cic4_allcuts_temp_lead[i];
            cic4_allcuts_sublead[i]=cic4_allcuts_temp_sublead[i]; }

        float cic4pf_allcuts_temp_lead[] = { 
            6,       4.7,       5.6,       3.6,
            10,       6.5,       5.6,       4.4,
            3.8,       2.5,       3.1,       2.2,
            0.0108,    0.0102,     0.028,     0.028,
            0.124,     0.092,     0.142,     0.063,
            0.94,      0.28,      0.94,      0.24,
            1,      0.99,      0.99,     0.028,
            1.5,         1.5,         1.5,         1.5};

        float cic4pf_allcuts_temp_sublead[] = { 
            6,       4.7,       5.6,       3.6,
            10,       6.5,       5.6,       4.4,
            3.8,       2.5,       3.1,       2.2,
            0.0108,    0.0102,     0.028,     0.028,
            0.124,     0.092,     0.142,     0.063,
            0.94,      0.28,      0.94,      0.24,
            1,      0.99,      0.99,     0.028,
            1.5,         1.5,         1.5,         1.5};

	float cic4pf8tev_allcuts_temp_lead[] = {     
	  6.3,       5.6,       5.8,       5.1,
	  18.9,         8,        10,       6.2,
	  4.5,       2.8,         4,      1.62,
	  0.0125,    0.0103,     0.029,     0.028,
	  0.141,     0.138,      0.12,     0.091,
	  0.94,      0.33,      0.94,      0.37,
	  1,     0.051,     0.054,     0.064,
	  1.5,         1.5,         1.5,         1.5};
	
        float cic4pf8tev_allcuts_temp_sublead[] = {  
	  6.3,       5.6,       5.8,       5.1,
	  18.9,         8,        10,       6.2,
	  4.5,       2.8,         4,      1.62,
	  0.0125,    0.0103,     0.029,     0.028,
	  0.141,     0.138,      0.12,     0.091,
	  0.94,      0.33,      0.94,      0.37,
	  1,     0.051,     0.054,     0.064,
	  1.5,         1.5,         1.5,         1.5};
	
        for(int i=0;i!=ncuts*ncat_cic4;++i) { 
	  if (cicVersion == "7TeV") {
	    cic4pf_allcuts_lead[i]    = cic4pf_allcuts_temp_lead[i];
	    cic4pf_allcuts_sublead[i] = cic4pf_allcuts_temp_sublead[i]; 
	  } else {
	    cic4pf_allcuts_lead[i]    = cic4pf8tev_allcuts_temp_lead[i];
	    cic4pf_allcuts_sublead[i] = cic4pf8tev_allcuts_temp_sublead[i]; 
	  }
        }
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
        float cic4_allcuts_temp_lead[] = { 
            3.2,      1.76,      1.39,      1.18,
            6.1,       2.7,       2.8,      0.66,
            3.4,      1.86,      1.67,      1.44,
            0.0104,    0.0094,     0.028,     0.025,
            0.076,      0.03,     0.047,     0.046,
            0.94,      0.41,      0.94,      0.34,
            1,      0.97,         1,         1,
            1.5,         1.5,         1.5,         1.5 };
        float cic4_allcuts_temp_sublead[] = { 
            3.2,      1.76,      1.39,      1.18,
            6.1,       2.7,       2.8,      0.66,
            3.4,      1.86,      1.67,      1.44,
            0.0104,    0.0094,     0.028,     0.025,
            0.076,      0.03,     0.047,     0.046,
            0.94,      0.41,      0.94,      0.34,
            1,      0.97,         1,         1,
            1.5,         1.5,         1.5,         1.5 };
        for(int i=0;i!=ncuts*ncat_cic6;++i) { 
            cic6_allcuts_lead[i]=cic6_allcuts_temp_lead[i];
            cic6_allcuts_sublead[i]=cic6_allcuts_temp_sublead[i]; }
        for(int i=0;i!=ncuts*ncat_cic4;++i) { 
            cic4_allcuts_lead[i]=cic4_allcuts_temp_lead[i];
            cic4_allcuts_sublead[i]=cic4_allcuts_temp_sublead[i]; }

        float cic4pf_allcuts_temp_lead[] = { 
            5.6,       4.3,       5.5,       2.8,
            9.5,         6,       5.1,       4.3,
            3.3,       2.3,       2.2,       1.2,
            0.0107,    0.0101,     0.028,     0.028,
            0.124,     0.092,     0.116,     0.059,
            0.94,      0.38,      0.94,      0.24,
            1,      0.99,      0.99,     0.028,
            1.5,         1.5,         1.5,         1.5};

        float cic4pf_allcuts_temp_sublead[] = { 
            5.6,       4.3,       5.5,       2.8,
            9.5,         6,       5.1,       4.3,
            3.3,       2.3,       2.2,       1.2,
            0.0107,    0.0101,     0.028,     0.028,
            0.124,     0.092,     0.116,     0.059,
            0.94,      0.38,      0.94,      0.24,
            1,      0.99,      0.99,     0.028,
            1.5,         1.5,         1.5,         1.5};

        float cic4pf8tev_allcuts_temp_lead[] = { 
	       5.9,       4.8,       5.3,       3.9,
	      11.2,         8,        10,       5.8,
	       4.2,       2.8,       2.8,      1.41,
	       9.4,       8.3,       5.9,       4.5,
	    0.0125,    0.0101,     0.028,     0.028,
	     0.141,     0.138,      0.12,     0.058,
	      0.94,      0.33,      0.94,      0.39,
	         1,     0.095,      0.77,       0.1,
	       1.5,         1.5,         1.5,         1.5};

        float cic4pf8tev_allcuts_temp_sublead[] = {   
	  5.9,       4.8,       5.3,       3.9,
	  11.2,         8,        10,       5.8,
	  4.2,       2.8,       2.8,      1.41,
	  9.4,       8.3,       5.9,       4.5,
	  0.0125,    0.0101,     0.028,     0.028,
	  0.141,     0.138,      0.12,     0.058,
	  0.94,      0.33,      0.94,      0.39,
	  1,     0.095,      0.77,       0.1,
	  1.5,         1.5,         1.5,         1.5};

        for(int i=0;i!=ncuts*ncat_cic4;++i) { 
	  if (cicVersion == "7TeV") {
            cic4pf_allcuts_lead[i]    = cic4pf_allcuts_temp_lead[i];
            cic4pf_allcuts_sublead[i] = cic4pf_allcuts_temp_sublead[i]; 
	  } else {
	    cic4pf_allcuts_lead[i]    = cic4pf8tev_allcuts_temp_lead[i];
            cic4pf_allcuts_sublead[i] = cic4pf8tev_allcuts_temp_sublead[i]; 
	  }
        }
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
        for(int i=0;i!=ncuts*ncat_cic6;++i) { 
            cic6_allcuts_lead[i]=cic6_allcuts_temp_lead[i];
            cic6_allcuts_sublead[i]=cic6_allcuts_temp_sublead[i]; }
        for(int i=0;i!=ncuts*ncat_cic4;++i) { 
            cic4_allcuts_lead[i]=cic4_allcuts_temp_lead[i];
            cic4_allcuts_sublead[i]=cic4_allcuts_temp_sublead[i]; }

        float cic4pf_allcuts_temp_lead[] = { 
            5.1,       3.8,       3.3,       2.6,
            7.5,         5,       4.7,         4,
            3,       2.2,      1.35,      0.38,
            0.0107,      0.01,     0.028,     0.027,
            0.124,     0.092,     0.116,     0.059,
            0.94,      0.39,      0.94,      0.24,
            1,         1,      0.99,      0.06,
            1.5,         1.5,         1.5,         1.5};
        float cic4pf_allcuts_temp_sublead[] = { 
            5.1,       3.8,       3.3,       2.6,
            7.5,         5,       4.7,         4,
            3,       2.2,      1.35,      0.38,
            0.0107,      0.01,     0.028,     0.027,
            0.124,     0.092,     0.116,     0.059,
            0.94,      0.39,      0.94,      0.24,
            1,         1,      0.99,      0.06,
            1.5,         1.5,         1.5,         1.5};
        for(int i=0;i!=ncuts*ncat_cic4;++i) { 
            cic4pf_allcuts_lead[i]    = cic4pf_allcuts_temp_lead[i];
            cic4pf_allcuts_sublead[i] = cic4pf_allcuts_temp_sublead[i]; 
        }
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
        for(int i=0;i!=ncuts*ncat_cic6;++i) { 
            cic6_allcuts_lead[i]=cic6_allcuts_temp_lead[i];
            cic6_allcuts_sublead[i]=cic6_allcuts_temp_sublead[i]; }
        for(int i=0;i!=ncuts*ncat_cic4;++i) { 
            cic4_allcuts_lead[i]=cic4_allcuts_temp_lead[i];
            cic4_allcuts_sublead[i]=cic4_allcuts_temp_sublead[i]; }

        float cic4pf_allcuts_temp_lead[] = { 
            4.5,       2.7,       2.6,       2.4,
            5.9,       4.8,       4.6,       3.7,
            2.5,      1.48,      0.87,      0.38,
            0.0106,      0.01,     0.028,     0.027,
            0.099,     0.092,     0.104,     0.059,
            0.94,      0.39,      0.94,      0.24,
            1,         1,         1,      0.78,
            1.5,         1.5,         1.5,         1.5};
        float cic4pf_allcuts_temp_sublead[] = { 
            4.5,       2.7,       2.6,       2.4,
            5.9,       4.8,       4.6,       3.7,
            2.5,      1.48,      0.87,      0.38,
            0.0106,      0.01,     0.028,     0.027,
            0.099,     0.092,     0.104,     0.059,
            0.94,      0.39,      0.94,      0.24,
            1,         1,         1,      0.78,
            1.5,         1.5,         1.5,         1.5};
        for(int i=0;i!=ncuts*ncat_cic4;++i) { 
            cic4pf_allcuts_lead[i]    = cic4pf_allcuts_temp_lead[i];
            cic4pf_allcuts_sublead[i] = cic4pf_allcuts_temp_sublead[i]; 
        }
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
        for(int i=0;i!=ncuts*ncat_cic6;++i) { 
            cic6_allcuts_lead[i]=cic6_allcuts_temp_lead[i];
            cic6_allcuts_sublead[i]=cic6_allcuts_temp_sublead[i]; }
        for(int i=0;i!=ncuts*ncat_cic4;++i) { 
            cic4_allcuts_lead[i]=cic4_allcuts_temp_lead[i];
            cic4_allcuts_sublead[i]=cic4_allcuts_temp_sublead[i]; }

        float cic4pf_allcuts_temp_lead[] = { 
            2.7,       2.4,       2.4,       2.3,
            5.9,       3.7,       4.6,       2.7,
            1.28,      1.15,      0.65,      0.38,
            .0106,    0.0099,     0.028,     0.027,
            0.095,     0.092,     0.065,     0.059,
            0.94,      0.39,      0.94,      0.24,
            1,         1,         1,      0.97,
            1.5,         1.5,         1.5,         1.5} ;

        float cic4pf_allcuts_temp_sublead[] = { 
            2.7,       2.4,       2.4,       2.3,
            5.9,       3.7,       4.6,       2.7,
            1.28,      1.15,      0.65,      0.38,
            .0106,     0.0099,     0.028,     0.027,
            0.095,     0.092,     0.065,     0.059,
            0.94,      0.39,      0.94,      0.24,
            1,         1,         1,      0.97,
            1.5,         1.5,         1.5,         1.5} ;

        for(int i=0;i!=ncuts*ncat_cic4;++i) { 
            cic4pf_allcuts_lead[i]    = cic4pf_allcuts_temp_lead[i];
            cic4pf_allcuts_sublead[i] = cic4pf_allcuts_temp_sublead[i]; 
        }
    } break;
    
    case(phoHYPERTIGHT5) : {
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
        for(int i=0;i!=ncuts*ncat_cic6;++i) { 
            cic6_allcuts_lead[i]=cic6_allcuts_temp_lead[i];
            cic6_allcuts_sublead[i]=cic6_allcuts_temp_sublead[i]; }
        for(int i=0;i!=ncuts*ncat_cic4;++i) { 
            cic4_allcuts_lead[i]=cic4_allcuts_temp_lead[i];
            cic4_allcuts_sublead[i]=cic4_allcuts_temp_sublead[i]; }


        float cic4pf_allcuts_temp_lead[] = { 
            2.5,       2.2,       2.3,       2.2,
            5.9,         3,       3.2,       2.7,
            1.05,      1.15,    0.0065,    0.0038,
            0.0106,    0.0099,     0.027,     0.027,
            0.086,     0.055,     0.065,    0.0032,
            0.94,      0.39,      0.95,      0.78,
            1,         1,         1,      0.97,
            1.5,         1.5,         1.5,         1.5} ;

        float cic4pf_allcuts_temp_sublead[] = { 
            2.5,       2.2,       2.3,       2.2,
            5.9,         3,       3.2,       2.7,
            1.05,      1.15,    0.0065,    0.0038,
            0.0106,    0.0099,     0.027,     0.027,
            0.086,     0.055,     0.065,    0.0032,
            0.94,      0.39,      0.95,      0.78,
            1,         1,         1,      0.97,
            1.5,         1.5,         1.5,         1.5} ;

        for(int i=0;i!=ncuts*ncat_cic4;++i) { 
            cic4pf_allcuts_lead[i]    = cic4pf_allcuts_temp_lead[i];
            cic4pf_allcuts_sublead[i] = cic4pf_allcuts_temp_sublead[i]; 
        }
    } break;
    
    case(phoHYPERTIGHT6) : {
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
        for(int i=0;i!=ncuts*ncat_cic6;++i) { 
            cic6_allcuts_lead[i]=cic6_allcuts_temp_lead[i];
            cic6_allcuts_sublead[i]=cic6_allcuts_temp_sublead[i]; }
        for(int i=0;i!=ncuts*ncat_cic4;++i) { 
            cic4_allcuts_lead[i]=cic4_allcuts_temp_lead[i];
            cic4_allcuts_sublead[i]=cic4_allcuts_temp_sublead[i]; }

        float cic4pf_allcuts_temp_lead[] = { 
            2.3,         2,       2.2,      1.88,
            5.5,       2.2,       2.1,      1.74,
            0.96,      0.99,   6.6e-05,   3.9e-05,
            0.0106,    0.0098,     0.026,     0.027,
            0.032,     0.035,      0.03,   3.2e-05,
            0.94,      0.39,      0.96,      0.86,
            1,         1,         1,      0.97,
            1.5,         1.5,         1.5,         1.5} ;

        float cic4pf_allcuts_temp_sublead[] = { 
            2.3,         2,       2.2,      1.88,
            5.5,       2.2,       2.1,      1.74,
            0.96,      0.99,   6.6e-05,   3.9e-05,
            0.0106,    0.0098,     0.026,     0.027,
            0.032,     0.035,      0.03,   3.2e-05,
            0.94,      0.39,      0.96,      0.86,
            1,         1,         1,      0.97,
            1.5,         1.5,         1.5,         1.5} ;

        for(int i=0;i!=ncuts*ncat_cic4;++i) { 
            cic4pf_allcuts_lead[i]    = cic4pf_allcuts_temp_lead[i];
            cic4pf_allcuts_sublead[i] = cic4pf_allcuts_temp_sublead[i]; 
        }
    } break;
    
    case(phoHYPERTIGHT7) : {
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
        for(int i=0;i!=ncuts*ncat_cic6;++i) { 
            cic6_allcuts_lead[i]=cic6_allcuts_temp_lead[i];
            cic6_allcuts_sublead[i]=cic6_allcuts_temp_sublead[i]; }
        for(int i=0;i!=ncuts*ncat_cic4;++i) { 
            cic4_allcuts_lead[i]=cic4_allcuts_temp_lead[i];
            cic4_allcuts_sublead[i]=cic4_allcuts_temp_sublead[i]; }

        float cic4pf_allcuts_temp_lead[] = { 
            2.2,      1.64,         2,      1.88,
            4.2,       2.2,       2.1,      1.39,
            0.96,      0.01,     1e-06,         0,
            0.0101,     0.009,     0.026,     0.027,
            0.017,   0.00156,     0.023,   3.2e-05,
            0.94,      0.68,      0.96,      0.86,
            1,         1,         1,      0.99,
            1.5,         1.5,         1.5,         1.5} ;

        float cic4pf_allcuts_temp_sublead[] = { 
            2.2,      1.64,         2,      1.88,
            4.2,       2.2,       2.1,      1.39,
            0.96,      0.01,     1e-06,         0,
            0.0101,     0.009,     0.026,     0.027,
            0.017,   0.00156,     0.023,   3.2e-05,
            0.94,      0.68,      0.96,      0.86,
            1,         1,         1,      0.99,
            1.5,         1.5,         1.5,         1.5} ;

        for(int i=0;i!=ncuts*ncat_cic4;++i) { 
            cic4pf_allcuts_lead[i]    = cic4pf_allcuts_temp_lead[i];
            cic4pf_allcuts_sublead[i] = cic4pf_allcuts_temp_sublead[i]; 
        }
    } break;
    
    case(phoHYPERTIGHT8) : {
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
        for(int i=0;i!=ncuts*ncat_cic6;++i) { 
            cic6_allcuts_lead[i]=cic6_allcuts_temp_lead[i];
            cic6_allcuts_sublead[i]=cic6_allcuts_temp_sublead[i]; }
        for(int i=0;i!=ncuts*ncat_cic4;++i) { 
            cic4_allcuts_lead[i]=cic4_allcuts_temp_lead[i];
            cic4_allcuts_sublead[i]=cic4_allcuts_temp_sublead[i]; }

        float cic4pf_allcuts_temp_lead[] = { 
            2,      1.62,      1.87,      1.88,
            4.2,      1.57,      1.38,      1.39,
            0.0097,  0.000101,         0,         0,
            0.0098,    0.0088,     0.025,     0.027,
            0.0164,   0.00156,   0.00023,   3.2e-05,
            0.95,      0.68,      0.97,      0.86,
            1,         1,         1,         1,
            1.5,         1.5,         1.5,         1.5} ;

        float cic4pf_allcuts_temp_sublead[] = { 
            2,      1.62,      1.87,      1.88,
            4.2,      1.57,      1.38,      1.39,
            0.0097,  0.000101,         0,         0,
            0.0098,    0.0088,     0.025,     0.027,
            0.0164,   0.00156,   0.00023,   3.2e-05,
            0.95,      0.68,      0.97,      0.86,
            1,         1,         1,         1,
            1.5,         1.5,         1.5,         1.5} ;

        for(int i=0;i!=ncuts*ncat_cic4;++i) { 
            cic4pf_allcuts_lead[i]    = cic4pf_allcuts_temp_lead[i];
            cic4pf_allcuts_sublead[i] = cic4pf_allcuts_temp_sublead[i]; 
        }
    } break;
    
    case(phoHYPERTIGHT9) : {
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
        for(int i=0;i!=ncuts*ncat_cic6;++i) { 
            cic6_allcuts_lead[i]=cic6_allcuts_temp_lead[i];
            cic6_allcuts_sublead[i]=cic6_allcuts_temp_sublead[i]; }
        for(int i=0;i!=ncuts*ncat_cic4;++i) { 
            cic4_allcuts_lead[i]=cic4_allcuts_temp_lead[i];
            cic4_allcuts_sublead[i]=cic4_allcuts_temp_sublead[i]; }

        float cic4pf_allcuts_temp_lead[] = { 
            2,      1.47,      1.87,      1.88,
            3.1,      1.37,      1.38,      1.39,
            8e-05,     1e-06,         0,         0,
            0.0098,    0.0087,     0.024,     0.027,
            0.000165,   0.00093,   0.00023,   3.2e-05,
            0.95,      0.68,      0.97,      0.86,
            1,         1,         1,         1,
            1.5,         1.5,         1.5,         1.5} ;

        float cic4pf_allcuts_temp_sublead[] = { 
            2,      1.47,      1.87,      1.88,
            3.1,      1.37,      1.38,      1.39,
            8e-05,     1e-06,         0,         0,
            0.0098,    0.0087,     0.024,     0.027,
            0.000165,   0.00093,   0.00023,   3.2e-05,
            0.95,      0.68,      0.97,      0.86,
            1,         1,         1,         1,
            1.5,         1.5,         1.5,         1.5} ;

        for(int i=0;i!=ncuts*ncat_cic4;++i) { 
            cic4pf_allcuts_lead[i]    = cic4pf_allcuts_temp_lead[i];
            cic4pf_allcuts_sublead[i] = cic4pf_allcuts_temp_sublead[i]; 
        }
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

int RedNtpTree::PhotonCiCSelectionLevel( int photon_index , bool electronVeto, int vtx_index, bool usePF ) 
{

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
  
  if (usePF)
    {
      val_tkisobad = -99;
      for(int iv=0; iv < nvertex; iv++) {
	if( pid_pfIsoCharged04ForCiC[photon_index][iv] > val_tkisobad) {
	  val_tkisobad = pid_pfIsoCharged04ForCiC[photon_index][iv];
	}
      }
      
      val_tkiso        = pid_pfIsoCharged03ForCiC[photon_index][vtx_index];
      val_ecaliso      = pid_pfIsoPhotons03ForCiC[photon_index];
      val_ecalisobad   = pid_pfIsoPhotons04ForCiC[photon_index];
      val_hcaliso = 0.;
      val_hcalisobad = 0.;
    }
  
  float val_sieie = pid_etawid[photon_index];
  float val_hoe = pid_HoverE[photon_index];
  float val_r9 = E9Phot[photon_index]/escRawPhot[photon_index];
  float val_drtotk_25_99 = pid_deltaRToTrackPhot[photon_index];
  float val_pixel = (float)hasPixelSeedPhot[photon_index];
  bool  val_conv = !hasMatchedPromptElePhot[photon_index];
  
  float isosumconst = 0.;
  float isosumconstbad = 0.;
  
  //PM 2011.05.30 Changed according to new values
  float rhofacbad=0.52, rhofac=0.17;
  if (usePF)
    {
      rhofacbad=0.23;
      rhofac=0.09;
    }

  float rho=rhoPF;
  if (usePF)
    rho=rhoAllJets;

  float pfIsoOffset=0;
  if (usePF)
    pfIsoOffset=2.5;
  float val_isosumoet=(val_tkiso+val_ecaliso+val_hcaliso+pfIsoOffset+isosumconst-rho*rhofac)*50./ptPhot[photon_index];
  float val_isosumoetbad=(val_tkisobad+val_ecalisobad+val_hcalisobad+pfIsoOffset+isosumconstbad-rho*rhofacbad)*50./ptPhot[photon_index];
  float val_trkisooet=(val_tkiso)*50./ptPhot[photon_index];

  /*
  tree_.runCIC=run;
  tree_.eventCIC=event;
  tree_.isosumoet=val_isosumoet;
  tree_.isoecalet=val_ecaliso; 
  tree_.isohcalet=val_hcaliso; 
  tree_.isotrackeret=val_tkiso; 
  tree_.isosumoetbad=val_isosumoetbad;
  tree_.isoecaletbad=val_ecalisobad;
  tree_.isohcaletbad=val_hcalisobad; 
  tree_.isotrackeretbad=val_tkisobad; 
  tree_.sieie=val_sieie; 
  tree_.hoe=val_hoe; 
  tree_.r9=val_r9; 
  tree_.drtotk_25_99=val_drtotk_25_99; 
  tree_.pixel=val_pixel; 
  myTree->Fill();
  */


  std::vector<std::vector<bool> > ph_passcut;
  ph_passcut.resize(phoNCUTLEVELS,std::vector<bool>(8,true) );

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
    
    if (!usePF)
      {
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
    else
      {

	ph_passcut[iCUTLEVEL][0] = (val_isosumoet        <=   cic4pf_cut_lead_isosumoet[iCUTLEVEL][photon_category]     );
	ph_passcut[iCUTLEVEL][1] = (val_isosumoetbad     <=   cic4pf_cut_lead_isosumoetbad[iCUTLEVEL][photon_category]  );
	ph_passcut[iCUTLEVEL][2] = (val_trkisooet        <=   cic4pf_cut_lead_trkisooet[iCUTLEVEL][photon_category]     );
	ph_passcut[iCUTLEVEL][3] = (val_sieie            <=   cic4pf_cut_lead_sieie[iCUTLEVEL][photon_category]         );
	ph_passcut[iCUTLEVEL][4] = (val_hoe              <=   cic4pf_cut_lead_hovere[iCUTLEVEL][photon_category]        );
	ph_passcut[iCUTLEVEL][5] = (val_r9               >=   cic4pf_cut_lead_r9[iCUTLEVEL][photon_category]            );// gt cut
	ph_passcut[iCUTLEVEL][6] = electronVeto ? (val_conv) : true; // gt cut
        
	bool ph_passcut_all = true;
	for(int icut=0;icut!=8;++icut) {
	  ph_passcut_all = ph_passcut_all && ph_passcut[iCUTLEVEL][icut];
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
  }
  
  return cutlevelpassed;

}

bool RedNtpTree::PhotonMITPreSelection( int photon_index, int vertex_index, bool electronVeto) {

  int photon_category = PhotonCategory(photon_index);

  float mitCuts_hoe[4]                 = {0.082,0.075,0.075,0.075};                                        
  float mitCuts_sieie[4]               = {0.014,0.014,0.034,0.034};                                        
  float mitCuts_ecaliso[4]             = {50,4,50,4};                                                      
  float mitCuts_hcaliso[4]             = {50,4,50,4};                                                      
  float mitCuts_trkiso[4]              = {50,4,50,4};                                                      
  //float mitCuts_hcalecal[4]            = {3,3,3,3};                                                        
  //float mitCuts_abstrkiso[4]           = {2.8,2.8,2.8,2.8};                                                
  //float mitCuts_trkiso_hollow03[4]     = {4,4,4,4};                                                       
  //float mitCuts_drtotk_25_99[4]	= {0.26,0.029,0.0062,0.0055};
  float mitCuts_pfiso[4]               = {4,4,4,4};
  
  float val_hoe        = pid_HoverE[photon_index];
  float val_sieie      = pid_etawid[photon_index];                                                          
  float val_ecaliso = pid_jurECAL03[photon_index] - 0.012*ptPhot[photon_index];                              
  float val_hcaliso = pid_twrHCAL03[photon_index] - 0.005*ptPhot[photon_index]; 
  float val_trkiso  = pid_hlwTrack03[photon_index] - 0.002*ptPhot[photon_index]; 
  
  //float val_hcalecal   = (pho_ecalsumetconedr03[photon_index]+pho_hcalsumetconedr03[photon_index]-rho_algo1*rhofac);                                             
  //float val_abstrkiso  = (*pho_tkiso_recvtx_030_002_0000_10_01)[photon_index][vertex_index];                
  //float val_trkiso_hollow03 = pho_trksumpthollowconedr03[photon_index];                                    
  //float val_drtotk_25_99 = pho_drtotk_25_99[photon_index];
  int val_pho_isconv = !hasMatchedPromptElePhot[photon_index];
  float val_pfiso02 = pid_pfIsoCharged02ForCiC[photon_index][vertex_index];
  
  if (val_hoe             >= mitCuts_hoe[photon_category]         ) return false;                                           
  if (val_sieie           >= mitCuts_sieie[photon_category]       ) return false;
  if (val_ecaliso         >= mitCuts_ecaliso[photon_category]     ) return false;
  if (val_hcaliso         >= mitCuts_hcaliso[photon_category]     ) return false;                                           
  if (val_trkiso          >= mitCuts_trkiso[photon_category]      ) return false;
  //if (val_hcalecal        >= mitCuts_hcalecal[photon_category]    ) return false;
  //if (val_abstrkiso       >= mitCuts_abstrkiso[photon_category]   ) return false;                   
  // if (val_drtotk_25_99    <  mitCuts_drtotk_25_99[photon_category]   ) return false; // Electron Rejection based on CiC for now
  if ((!val_pho_isconv && electronVeto) ) return false; // Electron Rejection based Conversion Safe Veto
  //if (val_trkiso_hollow03 >= mitCuts_trkiso_hollow03[photon_category]) return false;                                        
  if (val_pfiso02 >= mitCuts_pfiso[photon_category]) return false;            
  
  return true;
}

Float_t RedNtpTree::PhotonIDMVANew(Int_t iPhoton, Int_t vtx)  
{
  Float_t mva = 999.;

  double pfchargedisobad03=0.;
  for(int ivtx=0; ivtx<nvertex; ivtx++) {
    pfchargedisobad03=pid_pfIsoCharged03ForCiC[iPhoton][ivtx]>pfchargedisobad03?pid_pfIsoCharged03ForCiC[iPhoton][ivtx]:pfchargedisobad03;
  }

  tmva_photonid_pfchargedisogood03 = pid_pfIsoCharged03ForCiC[iPhoton][vtx];
  tmva_photonid_pfchargedisobad03  = pfchargedisobad03;
  tmva_photonid_pfphotoniso03      = pid_pfIsoPhotons03ForCiC[iPhoton];
  tmva_photonid_pfneutraliso03     = pid_pfIsoNeutrals03ForCiC[iPhoton]; 
  
  tmva_photonid_sieie        = pid_etawid[iPhoton];
  tmva_photonid_sieip        = sEtaPhiPhot[iPhoton];
  tmva_photonid_etawidth     = pid_scetawid[iPhoton];
  tmva_photonid_phiwidth     = pid_scphiwid[iPhoton];
  tmva_photonid_r9           = E9Phot[iPhoton]/escRawPhot[iPhoton];
  tmva_photonid_lambdaratio  = pid_lambdaRatio[iPhoton];
  
  tmva_photonid_s4ratio  = E4Phot[iPhoton]/E25Phot[iPhoton];
  tmva_photonid_eventrho = rhoAllJets;
  tmva_photonid_sceta    = etascPhot[iPhoton];
  
  float rr2=pid_esXwidth[iPhoton]*pid_esXwidth[iPhoton]+pid_esYwidth[iPhoton]*pid_esYwidth[iPhoton];
  tmva_photonid_ESEffSigmaRR = 0.0; 
  if(rr2>0. && rr2<999999.) 
    tmva_photonid_ESEffSigmaRR = sqrt(rr2);
    
  if (isEBPhot[iPhoton])
    mva = tmvaReaderID_Single_Barrel->EvaluateMVA("AdaBoost");
  else
    mva = tmvaReaderID_Single_Endcap->EvaluateMVA("AdaBoost");

  return mva;
}

void RedNtpTree::SetAllMVA() {
  tmvaReaderID_Single_Barrel = new TMVA::Reader("!Color:Silent");
  tmvaReaderID_Single_Barrel->AddVariable("myphoton_pfchargedisogood03",   &tmva_photonid_pfchargedisogood03 );
  tmvaReaderID_Single_Barrel->AddVariable("myphoton_pfchargedisobad03",   &tmva_photonid_pfchargedisobad03 );
  tmvaReaderID_Single_Barrel->AddVariable("myphoton_pfphotoniso03",   &tmva_photonid_pfphotoniso03 );
  tmvaReaderID_Single_Barrel->AddVariable("myphoton_sieie",   &tmva_photonid_sieie );
  tmvaReaderID_Single_Barrel->AddVariable("myphoton_sieip",   &tmva_photonid_sieip );
  tmvaReaderID_Single_Barrel->AddVariable("myphoton_etawidth",   &tmva_photonid_etawidth );
  tmvaReaderID_Single_Barrel->AddVariable("myphoton_phiwidth",   &tmva_photonid_phiwidth );
  tmvaReaderID_Single_Barrel->AddVariable("myphoton_r9",   &tmva_photonid_r9 );
  tmvaReaderID_Single_Barrel->AddVariable("myphoton_s4ratio",   &tmva_photonid_s4ratio );
  tmvaReaderID_Single_Barrel->AddVariable("myphoton_SCeta",   &tmva_photonid_sceta );
  tmvaReaderID_Single_Barrel->AddVariable("event_rho",   &tmva_photonid_eventrho );
  
  tmvaReaderID_Single_Endcap = new TMVA::Reader("!Color:Silent");
  tmvaReaderID_Single_Endcap->AddVariable("myphoton_pfchargedisogood03",   &tmva_photonid_pfchargedisogood03 );
  tmvaReaderID_Single_Endcap->AddVariable("myphoton_pfchargedisobad03",   &tmva_photonid_pfchargedisobad03 );
  tmvaReaderID_Single_Endcap->AddVariable("myphoton_pfphotoniso03",   &tmva_photonid_pfphotoniso03 );
  tmvaReaderID_Single_Endcap->AddVariable("myphoton_sieie",   &tmva_photonid_sieie );
  tmvaReaderID_Single_Endcap->AddVariable("myphoton_sieip",   &tmva_photonid_sieip );
  tmvaReaderID_Single_Endcap->AddVariable("myphoton_etawidth",   &tmva_photonid_etawidth );
  tmvaReaderID_Single_Endcap->AddVariable("myphoton_phiwidth",   &tmva_photonid_phiwidth );
  tmvaReaderID_Single_Endcap->AddVariable("myphoton_r9",   &tmva_photonid_r9 );
  tmvaReaderID_Single_Endcap->AddVariable("myphoton_s4ratio",   &tmva_photonid_s4ratio );
  tmvaReaderID_Single_Endcap->AddVariable("myphoton_SCeta",   &tmva_photonid_sceta );
  tmvaReaderID_Single_Endcap->AddVariable("event_rho",   &tmva_photonid_eventrho );
  tmvaReaderID_Single_Endcap->AddVariable("myphoton_ESEffSigmaRR",   &tmva_photonid_ESEffSigmaRR );


  tmvaReaderID_Single_Barrel->BookMVA("AdaBoost",photonLevelNewIDMVA_EB.c_str());
  tmvaReaderID_Single_Endcap->BookMVA("AdaBoost",photonLevelNewIDMVA_EE.c_str());
}
