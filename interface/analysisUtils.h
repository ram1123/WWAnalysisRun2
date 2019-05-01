#ifndef analysisUtils_h
#define analysisUtils_h

#include<iostream>
#include "TMath.h"
#include "TLorentzVector.h"
#include "TF1.h"
#include "TH1F.h"

#include "CondFormats/JetMETObjects/interface/JetCorrectorParameters.h"
#include "CondFormats/JetMETObjects/interface/FactorizedJetCorrector.h"
#include "CondFormats/JetMETObjects/interface/JetCorrectionUncertainty.h"

double getDeltaPhi(double phi1, double phi2 );

double deltaPhi(const double& phi1, const double& phi2);

double deltaEta(const double& eta1, const double& eta2);

double deltaR(const double& eta1, const double& phi1,
              const double& eta2, const double& phi2);

void computeAngles(TLorentzVector thep4H, TLorentzVector thep4Z1, TLorentzVector thep4M11, TLorentzVector thep4M12, 
		  TLorentzVector thep4Z2, TLorentzVector thep4M21, TLorentzVector thep4M22,
		  double& costheta1,  double& costheta2,  double& Phi,
		  double& costhetastar, double& Phi1);

double GetSFs_Lepton(double xAxis, double yAxis, TH1F *h1, TString posX="pt", TString posY="eta" );

double GetMin(double x, double y);
double GetMax(double x, double y);
//float getPUPPIweight(float puppipt, float puppieta );
float getPUPPIweight(TF1* puppisd_corrGEN, TF1* puppisd_corrRECO_cen, TF1* puppisd_corrRECO_for, float puppipt, float puppieta );
double func( double pt, double eta, JetCorrectionUncertainty *fJetUnc); 
double GetPt_MET(double pfMET,  double phi, double pz);
double GetEta_MET(double pfMET, double phi, double pz);
#endif
