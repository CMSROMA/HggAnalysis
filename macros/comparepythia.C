void comparepythia(char* default, char* check, char* name1, char* name2, char *var, int cuts,int bin, double min, double max, char* axis)
{

  gROOT->SetStyle("Plain");
   
  gStyle->SetPalette(1);
  gStyle->SetOptStat(1111111);  // Show overflow, underflow + SumOfWeights 
  gStyle->SetOptFit(111110); 
  gStyle->SetOptFile(1); 
  
  gStyle->SetMarkerStyle(20);
  gStyle->SetMarkerSize(.3);
  gStyle->SetMarkerColor(1);

  TCanvas* c0 = new TCanvas("c0"," ",200,10,500,500);
  c0->Clear();

//   gROOT->ProcessLine(".L fillPlot.C++");

  TFile *g1 = TFile::Open(check);
  TFile *g2 = TFile::Open(default);
//   TFile *g1 = TFile::Open(" /castor/cern.ch/user/d/delre/Higgs/reduced/redntp.42xv3.cicloose.-1.1/testVBF-pythia.root");
//   TFile *g2 = TFile::Open(" /castor/cern.ch/user/d/delre/Higgs/reduced/redntp.42xv3.cicloose.-1.1/testVBF-powheg.root");
  //  TFile *g2 = TFile::Open("root://pccmsrm23.cern.ch:1094//u2/xrootd/delre/Higgs/reduced/redntp.42xv3.preselectionCS.May10PromptV4.v3/merged/redntp_VBF_HToGG_M-120_7TeV-powheg-pythia6.root");
  fillPlot pythia((TTree*)g1->Get("AnaTree"), 1);
  fillPlot powheg((TTree*)g2->Get("AnaTree"), 1);

  // normalization
  int n_pythia = ((TH1D*)g1->Get("ptphotgen1"))->GetEntries();
  int n_powheg = ((TH1D*)g2->Get("ptphotgen1"))->GetEntries();

  // cuts
  float ptleadvbf=60.;
  //  float ptleadvbf=55.;
  float ptsubleadvbf=25.;
  float ptjetleadvbf=30.;
  //  float ptjetsubleadvbf=20.;
  float ptjetsubleadvbf=30.;  
  float deltaetavbf=3.;
  //float deltaetavbf=3.5;
  float zeppenvbf=2.5;
  float mjjvbf=500.;
  //float mjjvbf=350.;
  float deltaphi=2.6;
// 
//   float zeppenvbf=2.;
//   float mjjvbf=400.;
//   float deltaphi=0;

  pythia.DoPuReweight();
  //  pythia.DoPtReweight();
  //  powheg.DoPtReweight();
  powheg.DoPuReweight();
//   powheg.setCic(4);
//   pythia.setCic(4);
  powheg.setCic(15);
  pythia.setCic(15);

//   //photon cuts
//   pythia.Setcuts(ptleadvbf,ptsubleadvbf,-10000, -10000, -10000,-10000,0,0,-10000,-1,-1);
//   powheg.Setcuts(ptleadvbf,ptsubleadvbf,-10000, -10000, -10000,-10000,0,0,-10000,-1,-1);

  //no cut
  if(cuts == 0){
    pythia.Setcuts(-10000,-10000,-10000, -10000, -10000,-10000,0,0,-10000,0,-1,-1);
    powheg.Setcuts(-10000,-10000,-10000, -10000, -10000,-10000,0,0,-10000,0,-1,-1);
  }

  //pt phot cuts
  if(cuts == 1){
    pythia.Setcuts(ptleadvbf,ptsubleadvbf,-10000, -10000, -10000,-10000,0,0,-10000,0,-1,-1);
    powheg.Setcuts(ptleadvbf,ptsubleadvbf,-10000, -10000, -10000,-10000,0,0,-10000,0,-1,-1);
  }
  //+pt jet cuts
  if(cuts == 2){
    pythia.Setcuts(ptleadvbf,ptsubleadvbf,-10000, -10000, ptjetleadvbf,ptjetsubleadvbf,0,0,-10000,0,-1,-1);
    powheg.Setcuts(ptleadvbf,ptsubleadvbf,-10000, -10000, ptjetleadvbf,ptjetsubleadvbf,0,0,-10000,0,-1,-1);
  }
  //+deltaeta cut
  if(cuts == 3){
    pythia.Setcuts(ptleadvbf,ptsubleadvbf,-10000, -10000, ptjetleadvbf,ptjetsubleadvbf,deltaetavbf,0,-10000,0,-1,-1);
    powheg.Setcuts(ptleadvbf,ptsubleadvbf,-10000, -10000, ptjetleadvbf,ptjetsubleadvbf,deltaetavbf,0,-10000,0,-1,-1);
  }
  //+zeppend cut
  if(cuts == 4){
    pythia.Setcuts(ptleadvbf,ptsubleadvbf,-10000, -10000, ptjetleadvbf,ptjetsubleadvbf,deltaetavbf,zeppenvbf,-10000,0,-1,-1);
    powheg.Setcuts(ptleadvbf,ptsubleadvbf,-10000, -10000, ptjetleadvbf,ptjetsubleadvbf,deltaetavbf,zeppenvbf,-10000,0,-1,-1);
  }
  //+mjj
  if(cuts == 5){
    pythia.Setcuts(ptleadvbf,ptsubleadvbf,-10000, -10000, ptjetleadvbf,ptjetsubleadvbf,deltaetavbf,zeppenvbf,mjjvbf,0,-1,-1);
    powheg.Setcuts(ptleadvbf,ptsubleadvbf,-10000, -10000, ptjetleadvbf,ptjetsubleadvbf,deltaetavbf,zeppenvbf,mjjvbf,0,-1,-1);
  }
  //allcuts
  if(cuts == 6){
    pythia.Setcuts(ptleadvbf,ptsubleadvbf,-10000, -10000, ptjetleadvbf,ptjetsubleadvbf,deltaetavbf,zeppenvbf,mjjvbf,deltaphi,-1,-1,0);
    powheg.Setcuts(ptleadvbf,ptsubleadvbf,-10000, -10000, ptjetleadvbf,ptjetsubleadvbf,deltaetavbf,zeppenvbf,mjjvbf,deltaphi,-1,-1,0);
  }

  TH1D *pythiaplot = pythia.Plot(var,name2, bin, min, max);
  TH1D *powhegplot = powheg.Plot(var,name1, bin, min, max);

//   double effpythia = pythiaplot->GetEntries()/n_pythia;
//   double effpowheg = powhegplot->GetEntries()/n_powheg;
  double effpythia = pythiaplot->Integral(0,bin+1)/n_pythia;
  double effpowheg = powhegplot->Integral(0,bin+1)/n_powheg;
  double erreffpythia = sqrt(effpythia*(1-effpythia)/n_pythia);
  double erreffpowheg = sqrt(effpowheg*(1-effpowheg)/n_powheg);

  cout << "eff pythia = " << effpythia << " +/- " << erreffpythia << endl;
  cout << "eff powheg = " << effpowheg << " +/- " << erreffpowheg << endl;
  double ratioeff = effpythia/effpowheg;
  double errratioeff = sqrt ( erreffpythia*erreffpythia/(effpythia*effpythia) + erreffpowheg*erreffpowheg*effpythia*effpythia/pow(effpowheg,4) );
  cout << "ratio eff pythia/eff powheg " <<  ratioeff << " +/- " << errratioeff << endl;


  //legenda
  TLegendEntry *legge;
  TLegend *leg;
  leg = new TLegend(0.6,0.6,0.85,0.85);
  leg->SetFillStyle(0); leg->SetBorderSize(0); leg->SetTextSize(0.05);
  leg->SetFillColor(0);
  legge = leg->AddEntry(pythiaplot, name2, "p");
  legge = leg->AddEntry(powhegplot, name1, "f");

  powhegplot->SetStats(0);
  powhegplot->Draw();
  powhegplot->SetTitle("");
  powhegplot->SetXTitle(axis);
  pythiaplot->Scale(double(n_powheg)/n_pythia);
  pythiaplot->SetMarkerStyle(8);
  pythiaplot->SetMarkerSize(.9);
  pythiaplot->SetMarkerColor(kRed);
  pythiaplot->Draw("pesame");

  leg->Draw();

  char name[100];
  sprintf(name,"%s%s%i%s",var,"_",cuts,".gif");
  c0->SaveAs(name);

  pythiaplot->SetNormFactor(powhegplot->GetEntries());
  double maximumpow = powhegplot->GetMaximum();
  double maximumpyt = pythiaplot->GetMaximum();
  if(maximumpyt>maximumpow) powhegplot->SetMaximum(maximumpyt*1.1);
  else powhegplot->SetMaximum(maximumpow*1.1);
  powhegplot->Draw();
  pythiaplot->Draw("pesame");

  leg->Draw();

  sprintf(name,"%s%s%i%s",var,"_",cuts,"_norm.gif");
  c0->SaveAs(name);
 
 
}
