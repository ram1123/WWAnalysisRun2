#ifndef UTILS_HH
#define UTILS_HH

#include <TH1F.h>
#include "WWAnalysis/WWAnalysisRun2/interface/WVJJData.hh"
#include "WWAnalysis/WWAnalysisRun2/interface/ElectronID.hh"
#include "WWAnalysis/WWAnalysisRun2/interface/MuonID.hh"
#include "WWAnalysis/WWAnalysisRun2/interface/JetID.hh"

#include "CondFormats/JetMETObjects/interface/JetCorrectorParameters.h"
#include "CondFormats/JetMETObjects/interface/FactorizedJetCorrector.h"
#include "CondFormats/JetMETObjects/interface/JetCorrectionUncertainty.h"

float deltaR(float eta1, float phi1, float eta2, float phi2) {

  float dPhi = fabs(phi1-phi2);
  if (dPhi>6.283185308) dPhi -= 6.283185308;
  if (dPhi>3.141592654) dPhi = 3.141592654 - dPhi;

  float dEta = fabs(eta1-eta2);

  return sqrt( dPhi*dPhi + dEta*dEta );

}

float GetSFs_Lepton(double pt, double eta, TH1F *h1){
  if (pt > h1->GetYaxis()->GetXmax())
    pt = h1->GetYaxis()->GetXmax() - 1.0;
  if (pt < h1->GetYaxis()->GetXmin())
    pt = h1->GetYaxis()->GetXmin() + 1.0;
  
  return h1->GetBinContent(h1->GetXaxis()->FindFixBin(eta), h1->GetYaxis()->FindFixBin(pt));
}

double GetJECunc( double pt, double eta, JetCorrectionUncertainty *fJetUnc) { 
  fJetUnc->setJetPt ( pt  );
  fJetUnc->setJetEta( eta );
  return fJetUnc->getUncertainty(true);
}

#endif
