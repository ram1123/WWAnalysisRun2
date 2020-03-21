#ifndef SCALE_FACTORS_HH
#define SCALE_FACTORS_HH

#include "TFile.h"
#include "TF1.h"
#include "TH1.h"
#include "TH2.h"

#include "CondFormats/JetMETObjects/interface/JetCorrectorParameters.h"
#include "CondFormats/JetMETObjects/interface/FactorizedJetCorrector.h"
#include "CondFormats/JetMETObjects/interface/JetCorrectionUncertainty.h"

#include "WWAnalysis/WWAnalysisRun2/interface/Utils.hh"
//#include "WWAnalysis/WWAnalysisRun2/interface/analysisUtils.h"
//#include "CondFormats/BTauObjects/interface/BTagCalibration.h"
//#include "CondTools/BTau/interface/BTagCalibrationReader.h"
//#include "WVJJAna/Selection/interface/BtagUnc.hh"

class ScaleFactors {

public:
  ScaleFactors(int era) {
    if (era!=2016 && era!=2017) {
      std::cout << "Invalid Era! Using 2016." << std::endl;
      era_ = 2016;
    }
    else {
      era_=era;
    }
    setupSFs();
  };

  ~ScaleFactors() {};

  void setupSFs();

  float GetPUWeight(int nPUmean, int var) {
    if (var==1) return puWeightsUp->GetBinContent(nPUmean);
    else if (var==-1) return puWeightsDown->GetBinContent(nPUmean);

    return puWeights->GetBinContent(nPUmean);
  };

  float GetLeptonWeights(double lepPt, double lepEta, int pid) {
    float totWeight=1.0;
    if (era_==2016 && pid==11) {
      totWeight*=GetSFs_Lepton(lepPt, abs(lepEta), hIDIsoEle);
      totWeight*=GetSFs_Lepton(lepPt, abs(lepEta), hGSFCorrEle);
      totWeight/=GetSFs_Lepton(lepPt, abs(lepEta), hTriggerEle);
    }
    if (era_==2016 && pid==13) {
      totWeight*=GetSFs_Lepton(lepPt, abs(lepEta), hIDMuA);
      totWeight*=GetSFs_Lepton(lepPt, abs(lepEta), hIsoMuA);
      totWeight*=GetSFs_Lepton(lepPt, abs(lepEta), hTriggerMuA);
    }

    return totWeight;
  };

  TString GetBtagCSV() { return btagCSV; };

  TString GetAK4JecUnc() { return jecUncAK4; };

  TString GetAK8JecUnc() { return jecUncAK8; };

  float GetL1PFWeightJet(float jetPt, float jetEta) {
    return 1.0 - hL1PF_jetpt->GetBinContent(hL1PF_jetpt->FindBin(jetEta, jetPt));
  };

  float GetL1PFWeightPhot(float photPt, float photEta) {
    return 1.0 - hL1PF_photpt->GetBinContent(hL1PF_photpt->FindBin(photEta, photPt));
  }

private:

  int era_;

  // pileup

  //int nBinsPU;

  TFile* pileupFileMC;
  TH1D* puWeights;
  TH1D* puWeightsUp;
  TH1D* puWeightsDown;

  // electrons

  TFile* IDIsoEle;
  TH1F *hIDIsoEle;

  TFile* GSFCorrEle;
  TH1F *hGSFCorrEle;

  TFile* TriggerEle;
  TH1F* hTriggerEle;

  // muons

  TFile* IDMuA;
  TH1F *hIDMuA;
  TFile* IDMuB;
  TH1F *hIDMuB;

  TFile* IsoMuA;
  TH1F *hIsoMuA;
  TFile* IsoMuB;
  TH1F *hIsoMuB;

  TFile* TriggerMuA;
  TH1F* hTriggerMuA;
  TFile* TriggerMuB;
  TH1F* hTriggerMuB;

  // puppi

  TFile *fPuppi;
  TF1* puppisd_corrGEN;
  TF1* puppisd_corrRECO_cen;
  TF1* puppisd_corrRECO_for;

  // b-tagging

  TString btagCSV;
  //BTagCalibration calib;
  //BTagCalibrationReader bTagReader;
  
  // jec

  TString jecUncAK4;
  TString jecUncAK8;

  //JetCorrectorParameters paramAK4chs;
  //JetCorrectorParameters paramAK8puppi;

  //JetCorrectionUncertainty *fJetUnc_AK4chs;
  //JetCorrectionUncertainty *fJetUnc_AK8puppi;

  // L1 prefiring

  TFile *L1PF_jetpt;
  TH2F *hL1PF_jetpt;

  TFile *L1PF_photpt;
  TH2F *hL1PF_photpt;

  std::vector<TString> eleTriggerNames;
  std::vector<TString> muTriggerNames;
  

};


#endif
