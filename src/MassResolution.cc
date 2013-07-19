#include "MassResolution.h"
#include "RedNtpTree.h"
//----------------------------------------------------------//
// Project:   MassResolution
// Author:    Matt Kenzie (matthew.william.kenzie@cern.ch)
// Modified:  25/08/2011
// Admins:    Matth Kenzie (matthew.william.kenzie@cern.ch)
//---------------------------------------------------------//

/*
See MassResolution.h for instructions
*/

MassResolution::MassResolution(){}

void MassResolution::Setup(EnergyScaleCorrection* escale,RedNtpTree* event, int phot1, int phot2, int vtx_index, double beamspotSigma_in)
{
  escale_=escale;
  event_=event;
  beamspotSigma= beamspotSigma_in;
  
  leadPhoton= phot1;
  subleadPhoton= phot2;

  vertex = TVector3(event_->vx[vtx_index],event_->vy[vtx_index],event_->vz[vtx_index]);

  //lead_sc_pos    = leadPhoton->caloPosition();
  //sublead_sc_pos = subleadPhoton->caloPosition();

  lead_Eres    = event_->escRegrPhotError[leadPhoton];
  sublead_Eres = event_->escRegrPhotError[subleadPhoton];

  lead_r9      = event_->E9Phot[leadPhoton]/event_->escRawPhot[leadPhoton];
  sublead_r9   = event_->E9Phot[subleadPhoton]/event_->escRawPhot[subleadPhoton];

  lead_iDet    = (fabs(event->etascPhot[leadPhoton])<1.479);
  sublead_iDet = (fabs(event->etascPhot[subleadPhoton])<1.479);

  lead_p4=event->p4Phot(leadPhoton,vtx_index);
  sublead_p4=event->p4Phot(subleadPhoton,vtx_index);

  higgsMass=(lead_p4+sublead_p4).M();

}


// return the mass resolution given correct vertex
double MassResolution::massResolutionCorrVtx(){

  double lead_E = lead_p4.E();
  double sublead_E = sublead_p4.E();
  double alpha = lead_p4.Angle(sublead_p4.Vect());
  double lead_sig = leadPhotonResolution();
  double sublead_sig = subleadPhotonResolution();
  double alpha_sig = angleResolutionCorrVtx();
  
  return 0.5*higgsMass*TMath::Sqrt(((lead_sig*lead_sig)/(lead_E*lead_E))+((sublead_sig*sublead_sig)/(sublead_E*sublead_E))+((alpha_sig*alpha_sig)*(TMath::Sin(alpha)/(1.-TMath::Cos(alpha)))*(TMath::Sin(alpha)/(1.-TMath::Cos(alpha)))));

}
// return the mass resolution wrong vertex
double MassResolution::massResolutionWrongVtx(){
  
//  TLorentzVector lead_p4=leadPhoton->p4(vertex->X(),vertex->Y(),vertex->Z());
//  TLorentzVector sublead_p4=subleadPhoton->p4(vertex->X(),vertex->Y(),vertex->Z());

//  double lead_E = lead_p4.E();
//  double sublead_E = sublead_p4.E();
//  double alpha = lead_p4.Angle(sublead_p4.Vect());
//  double lead_sig = leadPhotonResolution();
//  double sublead_sig = subleadPhotonResolution();
  double alpha_sig = higgsMass*0.5*angleResolutionWrongVtx();
  
   double sigmaM = massResolutionEonly();
//  return 0.5*higgsMass*TMath::Sqrt(((lead_sig*lead_sig)/(lead_E*lead_E))+((sublead_sig*sublead_sig)/(sublead_E*sublead_E))+((alpha_sig*alpha_sig)));
  return TMath::Sqrt((sigmaM*sigmaM)+(alpha_sig*alpha_sig));

}

// return energy contribution to mass resolution only
double MassResolution::massResolutionEonly() {

  double lead_E = lead_p4.E(); 
  double sublead_E = sublead_p4.E(); 
  double lead_sig = leadPhotonResolution();
  double sublead_sig = subleadPhotonResolution();

  return 0.5*higgsMass*TMath::Sqrt((lead_sig*lead_sig)/(lead_E*lead_E)+(sublead_sig*sublead_sig)/(sublead_E*sublead_E));
}

double MassResolution::massResolutionAonly(){
	double aRes = angleResolution();
	return 0.5*higgsMass*aRes;
}


// return angle resolution given the vertex choice is correct
double MassResolution::angleResolutionCorrVtx() {
  return propagateDz(dzResolutionCorrVtx());
}

// return angle resolution given the vertex choice is wrong
double MassResolution::angleResolutionWrongVtx() {
  return propagateDz(dzResolutionWrongVtx());
}

// return angle resolution given a convolution of correct/wrong vertex as func of higgsPt
double MassResolution::angleResolution() {
  return propagateDz(dzResolution());
}
// return lead photon resolution without smearing
double MassResolution::leadPhotonResolutionNoSmear() {
  return lead_Eres;
}
// return sublead photon resolution without smearing
double MassResolution::subleadPhotonResolutionNoSmear() {
  return sublead_Eres;
}
// return lead photon resolution 
double MassResolution::leadPhotonResolution() {
//   TLorentzVector lead_p4=leadPhoton->p4(vertex->X(),vertex->Y(),vertex->Z());
//   bool sphericalLeadPhoton_=leadPhoton->isSphericalPhoton();
//  std::cout << " MassResolution -- Lead special ? " << sphericalLeadPhoton_ <<std::endl;
  return getPhotonResolution(lead_p4.E(),lead_Eres,lead_r9, event_->etascPhot[leadPhoton], lead_iDet);
}
// return sublead photon resolution
double MassResolution::subleadPhotonResolution() {
//   TLorentzVector sublead_p4=subleadPhoton->p4(vertex->X(),vertex->Y(),vertex->Z());
//   bool sphericalSubleadPhoton_=subleadPhoton->isSphericalPhoton();
//  std::cout << " MassResolution -- SubLead special ? " << sphericalSubleadPhoton_ <<std::endl;
  return getPhotonResolution(sublead_p4.E(),sublead_Eres,sublead_r9, event_->etascPhot[subleadPhoton],sublead_iDet);
}

// Actually compute resolution given a photon
double MassResolution::getPhotonResolution(double photonEnergy, double photonResolution, double r9,  double scEta, bool iDet) {





  /// Moved to config file PM 4/4/12 
  /// double categoryResolution = ispherical ? 0.0067*photonEnergy : _eSmearPars.smearing_sigma[myCategory]*photonEnergy;	
  double categoryResolution = escale_->getSmearing(event_->run,iDet,r9,fabs(scEta))*photonEnergy;	
  return TMath::Sqrt(categoryResolution*categoryResolution + photonResolution*photonResolution);

}

//return dz resolution given correct vertex (used 10mm)
double MassResolution::dzResolutionCorrVtx() {
  return 0.1;
}
//return dz resolution given wrong vertex (using sqrt(2)*5.8cm)
double MassResolution::dzResolutionWrongVtx() {
  return TMath::Sqrt(2.)*beamspotSigma;
}
  
//return dz resolution from dz wrong and dz right (stored in TGraph as func of higgsPt)
double MassResolution::dzResolution() { 
  return dz;
}

// propagate error on z to error on angle
double MassResolution::propagateDz(double dz){

//  TLorentzVector lead_p4=leadPhoton->p4(vertex->X(),vertex->Y(),vertex->Z());
//  TLorentzVector sublead_p4=subleadPhoton->p4(vertex->X(),vertex->Y(),vertex->Z());

//  double alpha = //lead_p4.Angle(sublead_p4.Vect());
//  if (alpha!= sublead_p4.Angle(lead_p4.Vect())) std::cout << "Error: Angle between photons not consistent" << std::endl;
  
  TVector3 LeadPosition = (TVector3(event_->xscPhot[leadPhoton],event_->yscPhot[leadPhoton],event_->zscPhot[leadPhoton])) - vertex;
  TVector3 SubLeadPosition = (TVector3(event_->xscPhot[subleadPhoton],event_->yscPhot[subleadPhoton],event_->zscPhot[subleadPhoton])) - vertex;

  /*
  double x1 = leadPhoton->caloPosition().X();
  double y1 = leadPhoton->caloPosition().Y();
  double z1 = leadPhoton->caloPosition().Z();

  double x2 = subleadPhoton->caloPosition().X();
  double y2 = subleadPhoton->caloPosition().Y();
  double z2 = subleadPhoton->caloPosition().Z();
*/
  
//   double x1 = leadPhoton->caloPosition().X()-vertex->X();
//   double y1 = leadPhoton->caloPosition().Y()-vertex->Y();
//   double z1 = leadPhoton->caloPosition().Z()-vertex->Z();
 
//   double x2 = subleadPhoton->caloPosition().X()-vertex->X();
//   double y2 = subleadPhoton->caloPosition().Y()-vertex->Y();
//   double z2 = subleadPhoton->caloPosition().Z()-vertex->Z();

  
  double r1 = LeadPosition.Mag();
  double r2 = SubLeadPosition.Mag();

  double cos_term = TMath::Cos(LeadPosition.Phi()-SubLeadPosition.Phi());
  double sech1 = SecH(LeadPosition.Eta());
  double sech2 = SecH(SubLeadPosition.Eta());
  double tanh1 = TanH(LeadPosition.Eta());
  double tanh2 = TanH(SubLeadPosition.Eta());

  double numerator1 = sech1*(sech1*tanh2-tanh1*sech2*cos_term);
  double numerator2 = sech2*(sech2*tanh1-tanh2*sech1*cos_term);
  double denominator = 1. - tanh1*tanh2 - sech1*sech2*cos_term;

  double ResTerm = (-1.*dz/denominator)*(numerator1/r1 + numerator2/r2);

  //double angleResolution = ResTerm*(1.-TMath::Cos(alpha))/TMath::Sin(alpha);
  double angleResolution = ResTerm;

  return angleResolution;

}

// utility functions
double MassResolution::SecH(double x){
  return 1.0/TMath::CosH(x);
}

double MassResolution::TanH(double x){
  return TMath::TanH(x);
}

