#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <string>

#include <TTree.h>
#include <TChain.h>
#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TH2.h>
#include <TH1F.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TVector2.h>
#include <TVector3.h>
#include <TAxis.h>
#include <TString.h>

#include "RedNtpTree.h"
#include "IsGJet.h"
#include "CrossSection.h"

#include "EnergyScaleCorrection.h"
//#include "EnergyScaleCorrectionSet.h"

using namespace std;

int main(int argc, char* argv[]) {

      //================ Parameters 
      if(argc < 6 || argc>14 ) {
        cout << "Usage:  ./tmp/redntpApp  listfile   outputfile   selection idmvaEB_file idmvaEE_file diPhotMVA_file jsonfile(optional) puweight(optional) ptweight(optional) scaleCorrections(optional) doPDFrew(optional) systJets(optional)\n" 
             << "    listfile:    list of root files incusing protocol eg dcap:/// .....\n"
             << "    outputfile:  name of output root file  eg output.root\n"
             << "    selection:   selection for preselecting events"  
             << "       options: superloose loose medium isem looseeg tighteg hggtighteg looseegpu tightegpu hggtightegpu preselection preselectionCS preselectionMVA preselectionMVAnoeleveto cicloose cicloosenoeleveto cicmedium cictight cicsuper cichyper cicpfloose cicpfloosenoeleveto cicpfmedium cicpftight cicpfsuper cicpfhyper mcass\n"
	     << "   weights for EB photonID MVA"
	     << "   weights for EE photonID MVA"
             << "   jsonfile: jsonfile used to select RUN/LS when looping over data. -1 if not used"
             << "   puweight: puweight for MC nPU reweighting. -1 if not used"
             << "   ptweight: ptweight for MC HiggsPt reweighting for GluGlu. -1 if not used"
             << "   scalCorrection: ...."
             << "   doPDFrew: if 1 it dumps info for PDF systematics "
             << "   systJet: if !=-1 file to perform JEC uncerntainties "
             << "   typesystJet: type of jet syst to perform 1=+1sigma JEC, 2=-1sigma JEC, 3=+10% JER,  4=-10% JER  "
             << endl;
        exit(-1);
      }


      //  1st option: nome del file contenete lista di root file
     
      // Input list
      char listName[500];
      sprintf(listName,argv[1]); 

      // Output filename (.root)  
      TString OutputFileName(argv[2]);
      
      // Name of input tree objects in (.root) files 
      char treeName[100] = "myanalysis/pippo";
      //sprintf(treeName,argv[2]);

      // fai TChain
      TChain *chain = new TChain(treeName);
      char pName[500];
      ifstream is(listName);
      if(! is.good()) {
         cout << "int main() >> ERROR : file " << listName << " not read" << endl;
         is.close();
         exit(-1);
      }
      cout << "Reading list : " << listName << " ......." << endl;
  
      while( is.getline(pName, 500, '\n') ) {
	 if (pName[0] == '#') continue;
	   //cout << "   Add: " << pName << endl;
	   chain->Add(pName); 
      }
      is.close();


      //4th option:  name of flat file with cuts
      char  selection[100];
      sprintf(selection,argv[3]);
      string finder(selection);
      if(finder == "") sprintf(selection,"looseeg");
      cout << "Photon selection is : " << selection << endl;

       // find cross section for this list
       float myxsec = CrossSection(listName);

       // filter for 2gam + jets. this is included in GJets samples but we use dedicated DiPhotonjets-madgraph
       int isGJetQCD = IsGJet(listName);

       // compute equivalent luminosity
       //Long64_t  ntot = chain->GetEntries();
       Long64_t  ntot = 1;
       double lumi = ntot/myxsec;
       cout << "#events: " << ntot << "\t xsec: " << myxsec << " pb\t equiv. lumi: " 
            << lumi/1000. << " fb-1"
            << endl;

       // run analysis code
       RedNtpTree tool(chain, OutputFileName);
       tool.cicVersion="7TeV";
       tool.SetNtotXsection( ntot, myxsec );
       tool.photonLevelNewIDMVA_EB=std::string(argv[4]);
       tool.photonLevelNewIDMVA_EE=std::string(argv[5]);
       tool.diPhotonMVAweights=std::string(argv[6]);

       if (argc>7 && std::string(argv[7]) != "-1")
 	 tool.SetJsonFile(argv[7]);

       if (argc>8 && std::string(argv[8]) != "-1")
	 tool.SetPuWeights(std::string(argv[8]));

       if (argc>9 && std::string(argv[9]) != "-1")
	 tool.SetPtWeights(std::string(argv[9]));

       if (argc>10 && std::string(argv[10]) != "-1")
	 {
	   //	   EnergyScaleCorrection::energyScaleParameters scaleCorrections;
	   //scaleCorrections.parameterSetName="";
	   TString scaleCorrectionFile(argv[10]);
	   //	   fillCorrections(scaleCorrectionSet,scaleCorrections); 
	   if (scaleCorrectionFile!="")
	     tool.setEnergyScaleCorrections(scaleCorrectionFile,"Hgg_eta_R9");
	 }

       if (argc>11 && std::string(argv[11]) != "-1")
	 tool.DoPDFWeighting();

       if (argc>12 && std::string(argv[12]) != "-1" && argc>13 && std::string(argv[13]) != "-1"){
	 TString scaleJetSysFile(argv[12]);
	 TString stringtypesyst(argv[13]);
	 char* endptr;
	 int typesyst = strtol (stringtypesyst, &endptr, 0);
	 tool.setJetSystematics(scaleJetSysFile,typesyst);;
       }
 
       std::cout << "DONE with settings starting loop" << std::endl;


       tool.Loop(isGJetQCD, selection);
}
