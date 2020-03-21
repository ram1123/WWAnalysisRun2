#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <cmath>
#include <algorithm>
#include <sstream>
#include <cstdio>
#include <cstdlib>
#include <stdint.h>
#include <iomanip>
#include <ctime>
#include <map>
#include <math.h>

#include "TROOT.h"
#include "TFile.h"
#include "TSystem.h"
#include "TStyle.h"
#include "TChain.h"
#include "TProfile.h"
#include "TMath.h"
#include "TString.h"
#include "TClass.h"
#include "TApplication.h"
#include "TLorentzVector.h"
#include "TF1.h"
#include "TH2.h"
#include <TClonesArray.h>

#include "CondFormats/JetMETObjects/interface/JetCorrectorParameters.h"
#include "CondFormats/JetMETObjects/interface/FactorizedJetCorrector.h"
#include "CondFormats/JetMETObjects/interface/JetCorrectionUncertainty.h"

#include "CondFormats/BTauObjects/interface/BTagCalibration.h"
#include "CondTools/BTau/interface/BTagCalibrationReader.h"
//#include "../BtagUnc.hh"

#include "BaconAna/DataFormats/interface/BaconAnaDefs.hh"
#include "BaconAna/DataFormats/interface/TEventInfo.hh"
#include "BaconAna/DataFormats/interface/TGenEventInfo.hh"
#include "BaconAna/DataFormats/interface/TGenParticle.hh"
#include "BaconAna/DataFormats/interface/TMuon.hh"
#include "BaconAna/DataFormats/interface/TElectron.hh"
#include "BaconAna/DataFormats/interface/TVertex.hh"
#include "BaconAna/DataFormats/interface/TJet.hh"
#include "BaconAna/DataFormats/interface/TAddJet.hh"
#include "BaconAna/Utils/interface/TTrigger.hh"
#include "BaconAna/DataFormats/interface/TLHEWeight.hh"

#include "WWAnalysis/WWAnalysisRun2/interface/ScaleFactors.hh"
#include "WWAnalysis/WWAnalysisRun2/interface/Utils.hh"
#include "WWAnalysis/WWAnalysisRun2/interface/WVJJData.hh"
#include "WWAnalysis/WWAnalysisRun2/interface/METzCalculator.h"

int main (int ac, char** av) {

  std::string inputFile = av[1];
  std::string outputFile = av[2];
  //int isMC = atoi(av[3]);
  int era = atoi(av[4]);

  const float MUON_MASS = 0.1056583745;
  const float ELE_MASS  = 0.000511;
  const float W_MASS = 80.385;
  //const float Z_MASS = 91.1876;

  //lepton cuts
  const float LEP_PT_VETO_CUT = 20;
  const float EL_PT_CUT = 35;
  const float EL_ETA_CUT = 2.5;
  const float MU_PT_CUT = 35;
  const float MU_ETA_CUT = 2.4;

  //ak8 jet cuts
  const float AK8_MIN_PT = 200;
  const float AK8_MAX_ETA = 2.4;
  const float AK8_MIN_SDM = 40;
  const float AK8_MAX_SDM = 150;

  //ak4 jet cuts
  //const float AK4_PT_VETO_CUT = 20;
  const float AK4_ETA_CUT = 2.4;
  const float AK4_PT_CUT = 30;
  const float AK4_JJ_MIN_M = 40.0;
  const float AK4_JJ_MAX_M = 150.0;
  const float VBF_MJJ_CUT= 500;

  //cleaning cuts
  const float AK8_LEP_DR_CUT = 1.0;
  const float AK4_AK8_DR_CUT = 0.8;
  const float AK4_DR_CUT = 0.3;

  //2016 csv tag thresholds
  const float CSV_LOOSE_2016 = 0.5426;
  const float CSV_MEDIUM_2016 = 0.8484;
  const float CSV_TIGHT_2016 = 0.9535;

  //2017 csv tag thresholds
  const float CSV_LOOSE_2017 = 0.5803;
  const float CSV_MEDIUM_2017 = 0.8838;
  const float CSV_TIGHT_2017 = 0.9693;
  //https://twiki.cern.ch/twiki/bin/viewauth/CMS/BtagRecommendation94X

  ScaleFactors scaleFactor(era);

  //misc
  //const int RUN_BOUND = 278820;

  std::vector<TLorentzVector> tightMuon;
  std::vector<TLorentzVector> tightEle;
  std::vector<TLorentzVector> goodAK4Jets;

  std::string iHLTFile="${CMSSW_BASE}/src/BaconAna/DataFormats/data/HLTFile_25ns";
  const std::string cmssw_base = getenv("CMSSW_BASE");
  std::string cmssw_base_env = "${CMSSW_BASE}";
  size_t start_pos = iHLTFile.find(cmssw_base_env);
  if(start_pos != std::string::npos) {
    iHLTFile.replace(start_pos, cmssw_base_env.length(), cmssw_base);
  }

  const baconhep::TTrigger triggerMenu(iHLTFile);

  //PILEUP WEIGHT FACTORS
  //int nbinsPU=75;
  //int nbinsPU=100;
  //TFile* pileupFileMC = TFile::Open("data/puWeights_80x_37ifb.root");
  //TFile* pileupFileMC = TFile::Open("data/puWeights_90x_2017.root");
  //TH1D* puWeights = (TH1D*)pileupFileMC->Get("puWeights");
  //TH1D* puWeightsUp = (TH1D*)pileupFileMC->Get("puWeightsUp");
  //TH1D* puWeightsDown = (TH1D*)pileupFileMC->Get("puWeightsDown");
  //puWeights->SetBins(nbinsPU,0,nbinsPU);
  //puWeightsUp->SetBins(nbinsPU,0,nbinsPU);
  //puWeightsDown->SetBins(nbinsPU,0,nbinsPU);

  //ELECTRON CORRECTIONS
  //TFile* IDIsoEle = TFile::Open("data/SF2016/egammaEffi_EGM2D_TightCutBasedIDSF.root","READ");
  //TFile* IDIsoEle = TFile::Open("data/SF2017/egammaEffi.txt_EGM2D_runBCDEF_passingTight94X.root","READ");
  //TH1F *hIDIsoEle = (TH1F*)IDIsoEle->Get("EGamma_SF2D");

  //TFile* GSFCorrEle = TFile::Open("data/SF2016/egammaEffi_SF2D_GSF_tracking.root","READ");
  //TFile* GSFCorrEle = TFile::Open("data/SF2017/egammaEffi.txt_EGM2D_runBCDEF_passingRECO.root","READ");
  //TH1F *hGSFCorrEle = (TH1F*)GSFCorrEle->Get("EGamma_SF2D");

  //TFile* TriggerEle = TFile::Open("data/SF2016/ElectronTrigger_SF.root","READ");
  //TH1F* hTriggerEle = (TH1F*)TriggerEle->Get("HLT_Ele27");
  
  //MUONS CORRECTIONS
  //not using currently
  //TFile* IDMuA = TFile::Open("data/SF2016/MuonID_RunBCDEF_23SepReReco_19p72fb.root","READ");
  //TH1F *hIDMuA = (TH1F*)IDMuA->Get("MC_NUM_TightID_DEN_genTracks_PAR_pt_eta/abseta_pt_ratio");
  //
  //TFile* IDMuB = TFile::Open("data/SF2016/MuonID_RunGH_23SepReReco_16p146fb.root","READ");
  //TH1F *hIDMuB = (TH1F*)IDMuB->Get("MC_NUM_TightID_DEN_genTracks_PAR_pt_eta/abseta_pt_ratio");
  //
  //TFile* IsoMuA = TFile::Open("data/SF2016/MuonIso_RunBCDEF_23SepReReco_19p72fb.root","READ");
  //TH1F *hIsoMuA = (TH1F*)IsoMuA->Get("TightISO_TightID_pt_eta/abseta_pt_ratio");
  //
  //TFile* IsoMuB = TFile::Open("data/SF2016/MuonIso_RunGH_23SepReReco_16p146fb.root","READ");
  //TH1F *hIsoMuB = (TH1F*)IsoMuB->Get("TightISO_TightID_pt_eta/abseta_pt_ratio");
  //
  //TFile* TriggerMuA = TFile::Open("data/SF2016/MuonTrigger_RunBCDEF_23SepReReco_19p72fb.root","READ");
  //TH1F* hTriggerMuA = (TH1F*)TriggerMuA->Get("IsoMu24_OR_IsoTkMu24_PtEtaBins/abseta_pt_ratio");
  //
  //TFile* TriggerMuB = TFile::Open("data/SF2016/MuonTrigger_RunGH_23SepReReco_16p146fb.root","READ");
  //TH1F* hTriggerMuB = (TH1F*)TriggerMuB->Get("IsoMu24_OR_IsoTkMu24_PtEtaBins/abseta_pt_ratio");

  //PUPPI CORRECTIONS
  //2016?
  TFile* file = TFile::Open("data/puppiCorr.root","READ");
  TF1* puppisd_corrGEN      = (TF1*)file->Get("puppiJECcorr_gen");
  TF1* puppisd_corrRECO_cen = (TF1*)file->Get("puppiJECcorr_reco_0eta1v3");
  TF1* puppisd_corrRECO_for = (TF1*)file->Get("puppiJECcorr_reco_1v3eta2v5");

  //B TAGGING
  //BTagCalibration calib("csvv2", "data/CSVv2_Moriond17_B_H.csv");//2016
  BTagCalibration calib("csvv2", scaleFactor.GetBtagCSV().Data());
  BTagCalibrationReader bTagReader(BTagEntry::OP_LOOSE,  // working point: can be OP_LOOSE, OP_MEDIUM, OP_TIGHT 
                                   "central",             // label for the central value (see the scale factor file)
                                   {"up","down"});        // vector of labels for systematics
  bTagReader.load(calib, BTagEntry::FLAV_B, "comb");      // use the "comb" measurements for b-jets
  bTagReader.load(calib, BTagEntry::FLAV_C, "comb");      // use the "comb" measurements for c-jets
  bTagReader.load(calib, BTagEntry::FLAV_UDSG, "incl");   // use the "incl" measurements for light jets

  //JET CORRECTIONS
  //2016
  //JetCorrectorParameters paramAK4chs("data/Summer16_23Sep2016V4_MC/Summer16_23Sep2016V4_MC_Uncertainty_AK4PFchs.txt");
  //JetCorrectorParameters paramAK8puppi("data/Summer16_23Sep2016V4_MC/Summer16_23Sep2016V4_MC_Uncertainty_AK8PFPuppi.txt");

  //2017
  JetCorrectorParameters paramAK4chs(scaleFactor.GetAK4JecUnc().Data());
  JetCorrectorParameters paramAK8puppi(scaleFactor.GetAK8JecUnc().Data());
  JetCorrectionUncertainty *fJetUnc_AK4chs = new JetCorrectionUncertainty(paramAK4chs);
  JetCorrectionUncertainty *fJetUnc_AK8puppi = new JetCorrectionUncertainty(paramAK8puppi);

  //L1 PREFIRING
  //TFile *L1PF_jetpt = TFile::Open("data/SF2017/L1prefiring_jetpt_2017BtoF.root");
  //TH2F *hL1PF_jetpt = (TH2F*) L1PF_jetpt->Get("L1prefiring_jetpt_2017BtoF");

  //
  //
  //   INPUT/OUTPUT
  //
  //

  TFile *of = new TFile(outputFile.c_str(),"RECREATE");
  TTree *ot = new TTree("Events","Events");
  WVJJData* WVJJTree = new WVJJData(ot);  
  TH1F *totalEvents = new TH1F("TotalEvents","TotalEvents",2,-1,1);

  baconhep::TEventInfo *info   = new baconhep::TEventInfo();
  baconhep::TGenEventInfo *gen = new baconhep::TGenEventInfo();
  TClonesArray *vertexArr   = new TClonesArray("baconhep::TVertex");
  TClonesArray *muonArr     = new TClonesArray("baconhep::TMuon");
  TClonesArray *electronArr = new TClonesArray("baconhep::TElectron");
  TClonesArray *AK4Arr      = new TClonesArray("baconhep::TJet");
  TClonesArray *PuppiAK8Arr = new TClonesArray("baconhep::TJet");
  TClonesArray *PuppiAK8AddArr = new TClonesArray("baconhep::TAddJet");
  TClonesArray *lheWgtArr   = new TClonesArray("baconhep::TLHEWeight");

  TFile *f=0;
  TTree *t=0;

  std::ifstream ifs;
  ifs.open(inputFile.data());
  assert(ifs.is_open());
  std::string line;
  while (getline(ifs,line)) {
    std::stringstream ss(line);
    std::string filetoopen;
    ss >> filetoopen;

    f = TFile::Open(TString("root://cmseos.fnal.gov/")+TString(filetoopen),"read");
    t = (TTree *)f->Get("Events");

    if (t==NULL) continue;
    
    t->SetBranchAddress("Info", &info);            TBranch *infoBr = t->GetBranch("Info");
    t->SetBranchAddress("PV",   &vertexArr);       TBranch *vertexBr = t->GetBranch("PV");
    t->SetBranchAddress("Muon", &muonArr);         TBranch *muonBr = t->GetBranch("Muon");
    t->SetBranchAddress("Electron", &electronArr); TBranch *electronBr = t->GetBranch("Electron");
    t->SetBranchAddress("AK4CHS", &AK4Arr);        TBranch *AK4Br = t->GetBranch("AK4CHS");
    t->SetBranchAddress("AK8Puppi", &PuppiAK8Arr); TBranch *puppiAK8Br = t->GetBranch("AK8Puppi");
    t->SetBranchAddress("AddAK8Puppi", &PuppiAK8AddArr); TBranch *puppiAK8AddBr = t->GetBranch("AddAK8Puppi");
    
    TBranch *genBr=0;
    TBranch *lheBr=0;
    
    if (t->GetListOfBranches()->FindObject("GenEvtInfo")) {
      t->SetBranchAddress("GenEvtInfo", &gen); genBr = t->GetBranch("GenEvtInfo");
    }
    if (t->GetListOfBranches()->FindObject("LHEWeight")) {
      t->SetBranchAddress("LHEWeight",&lheWgtArr); lheBr = t->GetBranch("LHEWeight");
    }
    
    totalEvents->SetBinContent(2,totalEvents->GetBinContent(2)+t->GetEntries());
    
    for (uint i=0; i < t->GetEntries(); i++) {
      //for (uint i=0; i < 200; i++) {
      WVJJTree->clearVars();
      infoBr->GetEntry(i);
      
      if (genBr) {
	genBr->GetEntry(i);
	if (gen->weight<0) totalEvents->Fill(-0.5);
      }

      if (era==2016) {
	if(!(triggerMenu.pass("HLT_IsoMu24_v*",info->triggerBits) || triggerMenu.pass("HLT_IsoTkMu24_v*",info->triggerBits) ||  
	     triggerMenu.pass("HLT_Ele27_WPTight_Gsf_v*",info->triggerBits))) continue;
      }
      if (era==2017) {
	if(!(triggerMenu.pass("HLT_IsoMu27_v*",info->triggerBits) || triggerMenu.pass("HLT_Ele32_WPTight_Gsf_v*",info->triggerBits) || 
	     triggerMenu.pass("HLT_Ele35_WPTight_Gsf_v*", info->triggerBits))) continue;
      }

      tightMuon.clear();
      tightEle.clear();
      
      WVJJTree->run = info->runNum;
      WVJJTree->evt = info->evtNum;
      WVJJTree->ls = info->lumiSec;
      
      vertexArr->Clear();
      vertexBr->GetEntry(i);
      
      WVJJTree->nPV = vertexArr->GetEntries();
      WVJJTree->nPU_mean = info->nPUmean;

      WVJJTree->puWeight = scaleFactor.GetPUWeight(info->nPUmean, 0);

      // LEPTON SELECTION
      
      muonArr->Clear();
      muonBr->GetEntry(i);
      
      int nTightEle=0;
      int nTightMu=0;
      int nVetoEle=0;
      int nVetoMu=0;
      
      for (int j=0; j < muonArr->GetEntries(); j++) {
	const baconhep::TMuon *mu = (baconhep::TMuon*)((*muonArr)[j]);
	
	if ( abs(mu->eta) > MU_ETA_CUT ) continue;
	if ( mu->pt < LEP_PT_VETO_CUT ) continue;
	
	if (!passMuonLoose(mu)) continue;
	nVetoMu++;
	
	if (!passMuonTight(mu)) continue;
	if ( mu->pt < MU_PT_CUT ) continue;
	nTightMu++;
	tightMuon.push_back(TLorentzVector(0,0,0,0));
	tightMuon.back().SetPtEtaPhiM(mu->pt, mu->eta, mu->phi, MUON_MASS);
	
	if ( mu->pt > WVJJTree->lep1_pt ) {
	  
	  WVJJTree->lep2_pt = WVJJTree->lep1_pt;
	  WVJJTree->lep2_eta = WVJJTree->lep1_eta;
	  WVJJTree->lep2_phi = WVJJTree->lep1_phi;
	  WVJJTree->lep2_m = WVJJTree->lep1_m;
	  WVJJTree->lep2_iso = WVJJTree->lep1_iso;
	  WVJJTree->lep2_q = WVJJTree->lep1_q;
	  
	  WVJJTree->lep1_pt = mu->pt;
	  WVJJTree->lep1_eta = mu->eta;
	  WVJJTree->lep1_phi = mu->phi;
	  WVJJTree->lep1_m = MUON_MASS;
	  WVJJTree->lep1_iso = (mu->chHadIso + TMath::Max(mu->neuHadIso + mu->gammaIso - 0.5*(mu->puIso), 0.0))/(mu->pt);
	  WVJJTree->lep1_q = mu->q;
	  
	}
	else if ( mu->pt > WVJJTree->lep2_pt ) {
	  
	  WVJJTree->lep2_pt = mu->pt;
	  WVJJTree->lep2_eta = mu->eta;
	  WVJJTree->lep2_phi = mu->phi;
	  WVJJTree->lep2_m = MUON_MASS;
	  WVJJTree->lep2_iso = (mu->chHadIso + TMath::Max(mu->neuHadIso + mu->gammaIso - 0.5*(mu->puIso), 0.0))/(mu->pt);
	  WVJJTree->lep2_q = mu->q;
	}
      }
      
      electronArr->Clear();
      electronBr->GetEntry(i);
      
      for (int j=0; j < electronArr->GetEntries(); j++) {
	const baconhep::TElectron *ele = (baconhep::TElectron*)((*electronArr)[j]);
	
	if ( abs(ele->eta) > EL_ETA_CUT ) continue;
	if ( ele->pt < LEP_PT_VETO_CUT ) continue;
	float eleIso = ele->chHadIso + TMath::Max(ele->gammaIso + ele->neuHadIso - info->rhoIso * eleEffArea(ele->eta,era), float(0.0));
	
	if (!passEleLoose(ele, eleIso, era)) continue;
	nVetoEle++;
	
	if (!passEleTight(ele, eleIso, era)) continue;
	nTightEle++;
	
	tightEle.push_back(TLorentzVector(0,0,0,0));
	tightEle.back().SetPtEtaPhiM(ele->pt,ele->eta,ele->phi,ELE_MASS);
	
	if ( ele->pt < EL_PT_CUT ) continue;
	
	//don't try to select electrons unless we don't already
	//have muons
	if (WVJJTree->lep1_m == MUON_MASS) continue;
	
	if ( ele->pt > WVJJTree->lep1_pt ) {
	  
	  WVJJTree->lep2_pt = WVJJTree->lep1_pt;
	  WVJJTree->lep2_eta = WVJJTree->lep1_eta;
	  WVJJTree->lep2_phi = WVJJTree->lep1_phi;
	  WVJJTree->lep2_m = WVJJTree->lep1_m;
	  WVJJTree->lep2_iso = WVJJTree->lep1_iso;
	  WVJJTree->lep2_q = WVJJTree->lep1_q;
	  
	  WVJJTree->lep1_pt = ele->pt;
	  WVJJTree->lep1_eta = ele->eta;
	  WVJJTree->lep1_phi = ele->phi;
	  WVJJTree->lep1_m = ELE_MASS;
	  WVJJTree->lep1_iso = eleIso;
	  WVJJTree->lep1_q = ele->q;
	  
	}
	else if ( ele->pt > WVJJTree->lep2_pt ) {
	  
	  WVJJTree->lep2_pt = ele->pt;
	  WVJJTree->lep2_eta = ele->eta;
	  WVJJTree->lep2_phi = ele->phi;
	  WVJJTree->lep2_m = ELE_MASS;
	  WVJJTree->lep2_iso = eleIso;
	  WVJJTree->lep2_q = ele->q;
	  
	}
      }
      
      //check conditions
      if(!(WVJJTree->lep1_pt>0)) continue;
      if ((nTightMu+nTightEle)==0) continue; //no leptons with required ID
      if((nVetoEle+nVetoMu)>2) continue;
      if(nTightMu>0 && nVetoEle>0) continue;
      if(nTightEle>0 && nVetoMu>0) continue;
      if(nTightMu==1 && nVetoMu>1) continue;
      if(nTightEle==1 && nVetoEle>1) continue;
      
      //if (WVJJTree->lep1_pt < 0) continue; // no lepton canddiates
      //if (nVetoLeps>2) continue; // too many leptons
      //if ((nTightMu==1||nTightEle==1) && nVetoLeps>1) continue;
      
      //muon scale variations
      //if (WVJJTree->lep1_m == MU_MASS) {}
      
      //electron scale variations
      if (WVJJTree->lep1_m == ELE_MASS) {
	WVJJTree->lep1_pt_scaleUp = 1.01 * WVJJTree->lep1_pt;
	WVJJTree->lep1_pt_scaleDn = 0.99 * WVJJTree->lep1_pt;
	if (WVJJTree->lep2_pt>0) {
	  WVJJTree->lep2_pt_scaleUp = 1.01 * WVJJTree->lep2_pt;
	  WVJJTree->lep2_pt_scaleDn = 0.99 * WVJJTree->lep2_pt;
	}
      }
      
      if (WVJJTree->lep1_pt > 0 && WVJJTree->lep2_pt > 0) {
	
	TLorentzVector lep1(0,0,0,0);
	lep1.SetPtEtaPhiM( WVJJTree->lep1_pt, WVJJTree->lep1_eta, WVJJTree->lep1_phi, WVJJTree->lep1_m );
	
	TLorentzVector lep2(0,0,0,0);
	lep2.SetPtEtaPhiM( WVJJTree->lep2_pt, WVJJTree->lep2_eta, WVJJTree->lep2_phi, WVJJTree->lep2_m );
	
	TLorentzVector dilep = lep1+lep2;
	
	WVJJTree->dilep_m = dilep.M();
	WVJJTree->dilep_pt = dilep.Pt();
	WVJJTree->dilep_eta = dilep.Eta();
	WVJJTree->dilep_phi = dilep.Phi();
	
	//dilepton scale variations
	
      }
      
      //lepton ID/iso/trigger efficiencies
      //if (WVJJTree->lep1_m == ELE_MASS) {
      WVJJTree->lep1_idEffWeight = scaleFactor.GetLeptonWeights(WVJJTree->lep1_pt, WVJJTree->lep1_eta, WVJJTree->lep1_m == ELE_MASS ? 11 : 13);
	//WVJJTree->lep1_idEffWeight = GetSFs_Lepton(WVJJTree->lep1_pt, WVJJTree->lep1_eta, hIDIsoEle);
	//WVJJTree->lep1_idEffWeight *= GetSFs_Lepton(WVJJTree->lep1_pt,WVJJTree->lep1_eta, hGSFCorrEle);
	//WVJJTree->lep1_idEffWeight *= GetSFs_Lepton(WVJJTree->lep1_pt,WVJJTree->lep1_eta, hTriggerEle);
	
      if (WVJJTree->lep2_pt>0) {
	WVJJTree->lep2_idEffWeight = scaleFactor.GetLeptonWeights(WVJJTree->lep2_pt, WVJJTree->lep2_eta, WVJJTree->lep1_m == ELE_MASS ? 11 : 13);
	  //WVJJTree->lep2_idEffWeight = GetSFs_Lepton(WVJJTree->lep2_pt, WVJJTree->lep2_eta, hIDIsoEle);
	  //WVJJTree->lep2_idEffWeight *= GetSFs_Lepton(WVJJTree->lep2_pt,WVJJTree->lep2_eta, hGSFCorrEle);
	  //WVJJTree->lep2_idEffWeight *= GetSFs_Lepton(WVJJTree->lep2_pt,WVJJTree->lep2_eta, hTriggerEle);
	  //do we even want the trigger eff for the subleading lepton?
      }

      //else if (WVJJTree->lep1_m == MUON_MASS) {
      //	WVJJTree->lep1_idEffWeight = ScaleFactors.GetLeptonWeights(WVJJTree->lep1_pt, WVJJTree->lep1_eta, 13);
      //	//WVJJTree->lep1_idEffWeight = 1.0; // not implemented yet
      //	//WVJJTree->lep2_idEffWeight = 1.0; // not implemented yet
      //
      //	if (WVJJTree->lep2_pt>0) {
      //	  WVJJTree->lep2_idEffWeight = ScaleFactors.GetLeptonWeights(WVJJTree->lep2_pt, WVJJTree->lep2_eta, 13);
      //	}
      //}
      
      // MET
      
      if (info->pfMETC < 0) continue;
      WVJJTree->MET_2017raw = info->pfMETC;
      WVJJTree->MET_phi = info->pfMETCphi;
      
      TLorentzVector tempMet(0,0,0,0);
      tempMet.SetPxPyPzE(WVJJTree->MET_2017raw*TMath::Cos(WVJJTree->MET_phi), 
			 WVJJTree->MET_2017raw*TMath::Sin(WVJJTree->MET_phi), 
			 0.0, WVJJTree->MET_2017raw);
      
      TLorentzVector tempMet_Up = tempMet;
      TLorentzVector tempMet_Dn = tempMet;
      
      AK4Arr->Clear();
      AK4Br->GetEntry(i);
      
      //2017 only correction for ECAL noise
      for (int j=0; j<AK4Arr->GetEntries(); j++) {
	const baconhep::TJet *ak4jet = (baconhep::TJet*)((*AK4Arr)[j]);
	//if(ak4jet->pt < AK4_PT_VETO_CUT) continue;
	
	float jecUnc = GetJECunc(ak4jet->pt, ak4jet->eta, fJetUnc_AK4chs);
	
	//HARDCODED
	if (era==2017 && ak4jet->ptRaw < 50 && abs(ak4jet->eta)>2.65 && abs(ak4jet->eta)<3.139) {
	  TLorentzVector tempRawJet(0,0,0,0);
	  tempRawJet.SetPtEtaPhiM(ak4jet->ptRaw, ak4jet->eta, ak4jet->phi, ak4jet->mass);
	  
	  tempMet+=tempRawJet;
	  tempMet_Up+=tempRawJet;
	  tempMet_Dn+=tempRawJet;
	}
	
	else {
	  TLorentzVector tempJet(0,0,0,0);
	  tempJet.SetPtEtaPhiM(ak4jet->pt, ak4jet->eta, ak4jet->phi, ak4jet->mass);
	  
	  TLorentzVector tempJetVar(0,0,0,0);
	  tempJetVar.SetPtEtaPhiM(ak4jet->pt*(1.0+jecUnc), ak4jet->eta, ak4jet->phi, ak4jet->mass);
	  tempMet_Up += tempJetVar - tempJet;
	  
	  tempJetVar.SetPtEtaPhiM(ak4jet->pt*(1.0-jecUnc), ak4jet->eta, ak4jet->phi, ak4jet->mass);
	  tempMet_Dn += tempJetVar - tempJet;
	}
      }
      
      WVJJTree->MET = tempMet.Pt();
      WVJJTree->MET_scaleUp = tempMet_Up.Pt();
      WVJJTree->MET_scaleDn = tempMet_Dn.Pt();
      
      if (WVJJTree->lep2_pt<0) {
	
	TLorentzVector lep1(0,0,0,0);
	lep1.SetPtEtaPhiM( WVJJTree->lep1_pt, WVJJTree->lep1_eta, WVJJTree->lep1_phi, WVJJTree->lep1_m );
	
	METzCalculator NeutrinoPz_type0;
	NeutrinoPz_type0.SetMET(tempMet);
	NeutrinoPz_type0.SetLepton(lep1);
	NeutrinoPz_type0.SetLeptonType(WVJJTree->lep1_m == ELE_MASS ? "el" : "mu");
	
	WVJJTree->neu_pz_type0 = NeutrinoPz_type0.Calculate();
	
	NeutrinoPz_type0.SetMET(tempMet_Up);
	WVJJTree->neu_pz_type0_scaleUp = NeutrinoPz_type0.Calculate();
	
	NeutrinoPz_type0.SetMET(tempMet_Dn);
	WVJJTree->neu_pz_type0_scaleDn = NeutrinoPz_type0.Calculate();
      
	TLorentzVector neutrino(0,0,0,0);
	neutrino.SetPxPyPzE(tempMet.Px(), tempMet.Py(), WVJJTree->neu_pz_type0, 
			    sqrt(tempMet.Px()*tempMet.Px()+tempMet.Py()*tempMet.Py()+ WVJJTree->neu_pz_type0* WVJJTree->neu_pz_type0));
	
	TLorentzVector bosLep = lep1+neutrino;
	
	WVJJTree->dilep_m = bosLep.M();
	WVJJTree->dilep_pt = bosLep.Pt();
	WVJJTree->dilep_eta = bosLep.Eta();
	WVJJTree->dilep_phi = bosLep.Phi();
	
	//dilepton scale variations
	
      }
      
      // AK8
      
      PuppiAK8Arr->Clear();
      puppiAK8Br->GetEntry(i);
      PuppiAK8AddArr->Clear();
      puppiAK8AddBr->GetEntry(i);
      
      float dmW = 3000.0;
      int nGoodFatJet=0;
      
      for (int j=0; j<PuppiAK8Arr->GetEntries(); j++) {
	const baconhep::TJet *ak8jet = (baconhep::TJet*)((*PuppiAK8Arr)[j]);
	const baconhep::TAddJet *ak8addjet = (baconhep::TAddJet*)((*PuppiAK8AddArr)[j]);
	
	if ( ak8jet->pt < AK8_MIN_PT ||  fabs(ak8jet->eta) > AK8_MAX_ETA ) continue;
	
	float jecUnc = GetJECunc(ak8jet->pt, ak8jet->eta, fJetUnc_AK8puppi);
	
	if ( ak8addjet->mass_sd0 < AK8_MIN_SDM || ak8addjet->mass_sd0 > AK8_MAX_SDM ) continue;
	if ( fabs(ak8addjet->mass_sd0 - W_MASS) > dmW ) continue;
	
	bool isClean=true;
	//lepton cleaning
	for ( std::size_t k=0; k<tightEle.size(); k++) {
	  if (deltaR(tightEle.at(k).Eta(), tightEle.at(k).Phi(),
		     ak8jet->eta, ak8jet->phi) < AK8_LEP_DR_CUT)
	    isClean = false;
	}
	for ( std::size_t k=0; k<tightMuon.size(); k++) {
	  if (deltaR(tightMuon.at(k).Eta(), tightMuon.at(k).Phi(),
		     ak8jet->eta, ak8jet->phi) < AK8_LEP_DR_CUT)
	    isClean = false;
	}
	if ( isClean == false ) continue;
	
	WVJJTree->bos_PuppiAK8_m_sd0 = ak8addjet->mass_sd0;
	WVJJTree->bos_PuppiAK8_m_sd0_corr = 
	  ak8addjet->mass_sd0 * getPUPPIweight(puppisd_corrGEN, puppisd_corrRECO_cen, puppisd_corrRECO_for,
					       ak8jet->pt, ak8jet->eta);
	WVJJTree->bos_PuppiAK8_tau2tau1 = ak8addjet->tau2/ak8addjet->tau1;
	WVJJTree->bos_PuppiAK8_pt = ak8jet->pt;
	WVJJTree->bos_PuppiAK8_eta = ak8jet->eta;
	WVJJTree->bos_PuppiAK8_phi = ak8jet->phi;
	
	WVJJTree->bos_PuppiAK8_m_sd0_corr_scaleUp =  WVJJTree->bos_PuppiAK8_m_sd0_corr*(1.0+jecUnc);
	WVJJTree->bos_PuppiAK8_m_sd0_corr_scaleDn =  WVJJTree->bos_PuppiAK8_m_sd0_corr*(1.0-jecUnc);
	WVJJTree->bos_PuppiAK8_pt_scaleUp =  WVJJTree->bos_PuppiAK8_pt*(1.0+jecUnc);
	WVJJTree->bos_PuppiAK8_pt_scaleDn =  WVJJTree->bos_PuppiAK8_pt*(1.0-jecUnc);

	WVJJTree->bos_PuppiAK8_e2_sdb1 = ak8addjet->e2_sdb1;
	WVJJTree->bos_PuppiAK8_e3_sdb1 = ak8addjet->e3_sdb1;
	WVJJTree->bos_PuppiAK8_e3_v1_sdb1 = ak8addjet->e3_v1_sdb1;
	WVJJTree->bos_PuppiAK8_e3_v2_sdb1 = ak8addjet->e3_v2_sdb1;
	WVJJTree->bos_PuppiAK8_e4_v1_sdb1 = ak8addjet->e4_v1_sdb1;
	WVJJTree->bos_PuppiAK8_e4_v2_sdb1 = ak8addjet->e4_v2_sdb1;
	//WVJJTree->bos_PuppiAK8_e2_sdb2 = ak8addjet->e2_sdb2;
	//WVJJTree->bos_PuppiAK8_e3_sdb2 = ak8addjet->e3_sdb2;
	//WVJJTree->bos_PuppiAK8_e3_v1_sdb2 = ak8addjet->e3_v1_sdb2;
	//WVJJTree->bos_PuppiAK8_e3_v2_sdb2 = ak8addjet->e3_v2_sdb2;
	//WVJJTree->bos_PuppiAK8_e4_v1_sdb2 = ak8addjet->e4_v1_sdb2;
	//WVJJTree->bos_PuppiAK8_e4_v2_sdb2 = ak8addjet->e4_v2_sdb2;
	
	dmW = fabs(ak8addjet->mass_sd0 - W_MASS);
	nGoodFatJet++;
      }
      
      goodAK4Jets.clear();
      AK4Arr->Clear();
      AK4Br->GetEntry(i);
      
      for (int j=0; j<AK4Arr->GetEntries(); j++) {
	const baconhep::TJet *ak4jet = (baconhep::TJet*)((*AK4Arr)[j]);
	//if(ak4jet->pt < AK4_PT_VETO_CUT) continue;
	
	float jecUnc = GetJECunc(ak4jet->pt, ak4jet->eta, fJetUnc_AK4chs);
	
	//jet energy scale variations
	if ( ak4jet->pt < AK4_PT_CUT && ak4jet->pt*(1.0+jecUnc) < AK4_PT_CUT && ak4jet->pt*(1.0-jecUnc) < AK4_PT_CUT) continue;
	if (!passAK4JetLoose(ak4jet,era)) continue;
	
	//2017 only
	if (ak4jet->ptRaw < 50 && abs(ak4jet->eta)>2.65 && abs(ak4jet->eta)<3.139) continue;
	
	if (abs(ak4jet->eta)<2.4 && ak4jet->pt>30) {
	  if (ak4jet->csv > (era==2016) ? CSV_LOOSE_2016 : CSV_LOOSE_2017) WVJJTree->nBtag_loose++;
	  if (ak4jet->csv > (era==2016) ? CSV_MEDIUM_2016 : CSV_MEDIUM_2017) WVJJTree->nBtag_medium++;
	  if (ak4jet->csv > (era==2016) ? CSV_TIGHT_2016 : CSV_TIGHT_2017) WVJJTree->nBtag_tight++;
	}
	
	bool isClean=true;
	// object cleaning
	
	if (nGoodFatJet>0) {
	  if (deltaR(WVJJTree->bos_PuppiAK8_eta, WVJJTree->bos_PuppiAK8_phi,
		     ak4jet->eta, ak4jet->phi) < AK4_AK8_DR_CUT) {
	    isClean = false;
	  }
	}
	
	for ( std::size_t k=0; k<goodAK4Jets.size(); k++) {
	  if (deltaR(goodAK4Jets.at(k).Eta(), goodAK4Jets.at(k).Phi(),
		     ak4jet->eta, ak4jet->phi) < AK4_DR_CUT) {
	    isClean = false;
	  }
	}
	for ( std::size_t k=0; k<tightEle.size(); k++) {
	  if (deltaR(tightEle.at(k).Eta(), tightEle.at(k).Phi(),
		     ak4jet->eta,   ak4jet->phi) < AK4_DR_CUT) {
	    isClean = false;
	  }
	}
	for ( std::size_t k=0; k<tightMuon.size(); k++) {
	  if (deltaR(tightMuon.at(k).Eta(), tightMuon.at(k).Phi(),
		     ak4jet->eta,   ak4jet->phi) < AK4_DR_CUT) {
	    isClean = false;
	  }
	}
	
	if ( isClean == false ) continue;
	
	goodAK4Jets.push_back(TLorentzVector(0,0,0,0));
	goodAK4Jets.back().SetPtEtaPhiM(ak4jet->pt, ak4jet->eta, ak4jet->phi, ak4jet->mass);
	
      }
      
      int nGoodDijet=0;
      
      uint sel1=1000, sel2=1000;
      if (nGoodFatJet==0) {
	TLorentzVector tmpV1, tmpV2;
	for (uint j=0; j<goodAK4Jets.size(); j++) {
	  if ( fabs(goodAK4Jets.at(j).Eta()) < AK4_ETA_CUT ) continue;
	  for(uint k=j+1; k<goodAK4Jets.size(); k++) {
	    if ( fabs(goodAK4Jets.at(k).Eta()) < AK4_ETA_CUT ) continue;
	    TLorentzVector tmpV=goodAK4Jets.at(j)+goodAK4Jets.at(k);
	    
	    if (tmpV.M()>=AK4_JJ_MIN_M && tmpV.M()<=AK4_JJ_MAX_M && tmpV.Pt()<AK8_MIN_PT) {
	      
	      WVJJTree->bos_j1_AK4_pt =  goodAK4Jets.at(j).Pt();
	      WVJJTree->bos_j1_AK4_eta = goodAK4Jets.at(j).Eta();
	      WVJJTree->bos_j1_AK4_phi = goodAK4Jets.at(j).Phi();
	      WVJJTree->bos_j1_AK4_m =   goodAK4Jets.at(j).M();
	      
	      WVJJTree->bos_j2_AK4_pt =  goodAK4Jets.at(k).Pt();
	      WVJJTree->bos_j2_AK4_eta = goodAK4Jets.at(k).Eta();
	      WVJJTree->bos_j2_AK4_phi = goodAK4Jets.at(k).Phi();
	      WVJJTree->bos_j2_AK4_m =   goodAK4Jets.at(k).M();
	      
	      WVJJTree->bos_AK4AK4_pt =  tmpV.Pt();
	      WVJJTree->bos_AK4AK4_eta = tmpV.Eta();
	      WVJJTree->bos_AK4AK4_phi = tmpV.Phi();
	      WVJJTree->bos_AK4AK4_m =   tmpV.M();
	      
	      sel1=j; sel2=k;
	      nGoodDijet++;
	      
	    }
	    if (nGoodDijet>0) break;
	  }
	  if (nGoodDijet>0) break;
	}
	
	if (nGoodDijet==0) continue;
	
	float jecUnc1 = GetJECunc(WVJJTree->bos_j1_AK4_pt, WVJJTree->bos_j1_AK4_eta, fJetUnc_AK4chs);
	float jecUnc2 = GetJECunc(WVJJTree->bos_j2_AK4_pt, WVJJTree->bos_j2_AK4_eta, fJetUnc_AK4chs);
	
	WVJJTree->bos_j1_AK4_pt_scaleUp = WVJJTree->bos_j1_AK4_pt*(1.0+jecUnc1);
	WVJJTree->bos_j1_AK4_pt_scaleDn = WVJJTree->bos_j1_AK4_pt*(1.0-jecUnc1);
	WVJJTree->bos_j1_AK4_m_scaleUp = WVJJTree->bos_j1_AK4_m*(1.0+jecUnc1);
	WVJJTree->bos_j1_AK4_m_scaleDn = WVJJTree->bos_j1_AK4_m*(1.0-jecUnc1);
	
	WVJJTree->bos_j2_AK4_pt_scaleUp = WVJJTree->bos_j2_AK4_pt*(1.0+jecUnc2);
	WVJJTree->bos_j2_AK4_pt_scaleDn = WVJJTree->bos_j2_AK4_pt*(1.0-jecUnc2);
	WVJJTree->bos_j2_AK4_m_scaleUp = WVJJTree->bos_j2_AK4_m*(1.0+jecUnc2);
	WVJJTree->bos_j2_AK4_m_scaleDn = WVJJTree->bos_j2_AK4_m*(1.0-jecUnc2);
	
	TLorentzVector tempBos1(0,0,0,0);
	TLorentzVector tempBos2(0,0,0,0);
	
	tempBos1.SetPtEtaPhiM(WVJJTree->bos_j1_AK4_pt_scaleUp, WVJJTree->bos_j1_AK4_eta,
			      WVJJTree->bos_j1_AK4_phi, WVJJTree->bos_j1_AK4_m_scaleUp);
	tempBos2.SetPtEtaPhiM(WVJJTree->bos_j2_AK4_pt_scaleUp, WVJJTree->bos_j2_AK4_eta,
			      WVJJTree->bos_j2_AK4_phi, WVJJTree->bos_j2_AK4_m_scaleUp);
	
	WVJJTree->bos_AK4AK4_pt_scaleUp = (tempBos1+tempBos2).Pt();
	WVJJTree->bos_AK4AK4_m_scaleUp = (tempBos1+tempBos2).M();
	
	tempBos1.SetPtEtaPhiM(WVJJTree->bos_j1_AK4_pt_scaleDn, WVJJTree->bos_j1_AK4_eta,
			      WVJJTree->bos_j1_AK4_phi, WVJJTree->bos_j1_AK4_m_scaleDn);
	tempBos2.SetPtEtaPhiM(WVJJTree->bos_j2_AK4_pt_scaleDn, WVJJTree->bos_j2_AK4_eta,
			      WVJJTree->bos_j2_AK4_phi, WVJJTree->bos_j2_AK4_m_scaleDn);
	
	WVJJTree->bos_AK4AK4_pt_scaleDn = (tempBos1+tempBos2).Pt();
	WVJJTree->bos_AK4AK4_m_scaleDn = (tempBos1+tempBos2).M();
	
      }
      
      //check we have a hadronic boson candidate
      if ( nGoodFatJet == 0 && nGoodDijet == 0 ) continue;
      
      float tmpMassMax = 0.0;
      int vbf1=-1, vbf2=-1;
      
      for (uint j=0; j<goodAK4Jets.size(); j++) {
	if (j==sel1 || j==sel2) continue;
	for(uint k=j+1; k<goodAK4Jets.size(); k++) {
	  if (k==sel1 || k==sel2) continue;
	  TLorentzVector tempVBF = goodAK4Jets.at(j) + goodAK4Jets.at(k);
	  if ( tempVBF.M() < VBF_MJJ_CUT ) continue;
	  if ( tempVBF.M() < tmpMassMax ) continue;
	  tmpMassMax = tempVBF.M();
	  vbf1=j; vbf2=k;
	}
      }
      
      if (vbf1==-1 && vbf2==-1) continue;
      
      TLorentzVector tempVBF = goodAK4Jets.at(vbf1) + goodAK4Jets.at(vbf2);
      
      WVJJTree->vbf1_AK4_pt = goodAK4Jets.at(vbf1).Pt();
      WVJJTree->vbf1_AK4_eta = goodAK4Jets.at(vbf1).Eta();
      WVJJTree->vbf1_AK4_phi = goodAK4Jets.at(vbf1).Phi();
      WVJJTree->vbf1_AK4_m = goodAK4Jets.at(vbf1).M();
      
      WVJJTree->vbf2_AK4_pt = goodAK4Jets.at(vbf2).Pt();
      WVJJTree->vbf2_AK4_eta = goodAK4Jets.at(vbf2).Eta();
      WVJJTree->vbf2_AK4_phi = goodAK4Jets.at(vbf2).Phi();
      WVJJTree->vbf2_AK4_m = goodAK4Jets.at(vbf2).M();
      
      WVJJTree->vbf_pt = tempVBF.Pt();
      WVJJTree->vbf_eta = tempVBF.Eta();
      WVJJTree->vbf_phi = tempVBF.Phi();
      WVJJTree->vbf_m = tempVBF.M();
      
      WVJJTree->vbf_deta = abs( WVJJTree->vbf2_AK4_eta - WVJJTree->vbf1_AK4_eta );
      
      TLorentzVector tempVBF1(0,0,0,0);
      TLorentzVector tempVBF2(0,0,0,0);
      
      float jecUnc1 = GetJECunc(WVJJTree->vbf1_AK4_pt, WVJJTree->vbf1_AK4_eta, fJetUnc_AK4chs);
      float jecUnc2 = GetJECunc(WVJJTree->vbf2_AK4_pt, WVJJTree->vbf2_AK4_eta, fJetUnc_AK4chs);
      
      WVJJTree->vbf1_AK4_pt_scaleUp = WVJJTree->vbf1_AK4_pt*(1.0+jecUnc1);
      WVJJTree->vbf1_AK4_pt_scaleDn = WVJJTree->vbf1_AK4_pt*(1.0-jecUnc1);
      WVJJTree->vbf1_AK4_m_scaleUp = WVJJTree->vbf1_AK4_m*(1.0+jecUnc1);
      WVJJTree->vbf1_AK4_m_scaleDn = WVJJTree->vbf1_AK4_m*(1.0-jecUnc1);
      
      WVJJTree->vbf2_AK4_pt_scaleUp = WVJJTree->vbf2_AK4_pt*(1.0+jecUnc2);
      WVJJTree->vbf2_AK4_pt_scaleDn = WVJJTree->vbf2_AK4_pt*(1.0-jecUnc2);
      WVJJTree->vbf2_AK4_m_scaleUp = WVJJTree->vbf2_AK4_m*(1.0+jecUnc2);
      WVJJTree->vbf2_AK4_m_scaleDn = WVJJTree->vbf2_AK4_m*(1.0-jecUnc2);
      
      tempVBF1.SetPtEtaPhiM(WVJJTree->vbf1_AK4_pt_scaleUp, WVJJTree->vbf1_AK4_eta,
			    WVJJTree->vbf1_AK4_phi, WVJJTree->vbf1_AK4_m_scaleUp);
      tempVBF2.SetPtEtaPhiM(WVJJTree->vbf2_AK4_pt_scaleUp, WVJJTree->vbf2_AK4_eta,
			    WVJJTree->vbf2_AK4_phi, WVJJTree->vbf2_AK4_m_scaleUp);
      
      WVJJTree->vbf_pt_scaleUp = (tempVBF1+tempVBF2).Pt();
      WVJJTree->vbf_m_scaleUp = (tempVBF1+tempVBF2).M();
      
      tempVBF1.SetPtEtaPhiM(WVJJTree->vbf1_AK4_pt_scaleDn, WVJJTree->vbf1_AK4_eta,
			    WVJJTree->vbf1_AK4_phi, WVJJTree->vbf1_AK4_m_scaleDn);
      tempVBF2.SetPtEtaPhiM(WVJJTree->vbf2_AK4_pt_scaleDn, WVJJTree->vbf2_AK4_eta,
			    WVJJTree->vbf2_AK4_phi, WVJJTree->vbf2_AK4_m_scaleDn);
      
      WVJJTree->vbf_pt_scaleDn = (tempVBF1+tempVBF2).Pt();
      WVJJTree->vbf_m_scaleDn = (tempVBF1+tempVBF2).M();
      
      TLorentzVector bosHad(0,0,0,0), bosHad_up(0,0,0,0), bosHad_dn(0,0,0,0);
      TLorentzVector bosLep(0,0,0,0), bosLep_up(0,0,0,0), bosLep_dn(0,0,0,0);
      
      //boosted event
      if (WVJJTree->bos_PuppiAK8_m_sd0_corr > 0) {
	bosHad.SetPtEtaPhiM(WVJJTree->bos_PuppiAK8_pt, WVJJTree->bos_PuppiAK8_eta,
			    WVJJTree->bos_PuppiAK8_phi, WVJJTree->bos_PuppiAK8_m_sd0_corr);
	bosHad_up.SetPtEtaPhiM(WVJJTree->bos_PuppiAK8_pt_scaleUp, WVJJTree->bos_PuppiAK8_eta,
			       WVJJTree->bos_PuppiAK8_phi, WVJJTree->bos_PuppiAK8_m_sd0_corr_scaleUp);
	bosHad_dn.SetPtEtaPhiM(WVJJTree->bos_PuppiAK8_pt_scaleDn, WVJJTree->bos_PuppiAK8_eta,
			       WVJJTree->bos_PuppiAK8_phi, WVJJTree->bos_PuppiAK8_m_sd0_corr_scaleDn);
      }
    //resolved event
      else if (WVJJTree->bos_AK4AK4_m > 0) {
	bosHad.SetPtEtaPhiM(WVJJTree->bos_AK4AK4_pt, WVJJTree->bos_AK4AK4_eta, 
			    WVJJTree->bos_AK4AK4_phi, WVJJTree->bos_AK4AK4_m);
	bosHad_up.SetPtEtaPhiM(WVJJTree->bos_AK4AK4_pt_scaleUp, WVJJTree->bos_AK4AK4_eta, 
			       WVJJTree->bos_AK4AK4_phi, WVJJTree->bos_AK4AK4_m_scaleUp);
	bosHad_dn.SetPtEtaPhiM(WVJJTree->bos_AK4AK4_pt_scaleDn, WVJJTree->bos_AK4AK4_eta, 
			       WVJJTree->bos_AK4AK4_phi, WVJJTree->bos_AK4AK4_m_scaleDn);
      }
      
      bosLep.SetPtEtaPhiM(WVJJTree->dilep_pt, WVJJTree->dilep_eta, WVJJTree->dilep_phi, WVJJTree->dilep_m);
      //lepton variations not implemented
      bosLep_up.SetPtEtaPhiM(WVJJTree->dilep_pt, WVJJTree->dilep_eta, WVJJTree->dilep_phi, WVJJTree->dilep_m);
      bosLep_dn.SetPtEtaPhiM(WVJJTree->dilep_pt, WVJJTree->dilep_eta, WVJJTree->dilep_phi, WVJJTree->dilep_m);
      
      TLorentzVector diBos = bosHad+bosLep;
      TLorentzVector diBos_up = bosHad_up+bosLep_up;
      TLorentzVector diBos_dn = bosHad_dn+bosLep_dn;
      
      WVJJTree->dibos_m = diBos.M();
      WVJJTree->dibos_pt = diBos.Pt();
      WVJJTree->dibos_eta = diBos.Eta();
      WVJJTree->dibos_phi = diBos.Phi();
      
      WVJJTree->dibos_m_scaleUp = diBos_up.M();
      WVJJTree->dibos_m_scaleDn = diBos_dn.M();
      WVJJTree->dibos_pt_scaleUp = diBos_up.Pt();
      WVJJTree->dibos_pt_scaleDn = diBos_dn.Pt();
      
      if (WVJJTree->lep2_pt < 0) {
	WVJJTree->bosCent = std::min( std::min(bosHad.Eta(), bosLep.Eta()) - std::min(WVJJTree->vbf1_AK4_eta, WVJJTree->vbf2_AK4_eta), 
				      std::max(WVJJTree->vbf1_AK4_eta, WVJJTree->vbf2_AK4_eta) - std::max(bosHad.Eta(), bosLep.Eta()) );
      }
      
      WVJJTree->zeppLep = bosLep.Eta() - 0.5*(WVJJTree->vbf1_AK4_eta + WVJJTree->vbf2_AK4_eta);
      WVJJTree->zeppHad = bosHad.Eta() - 0.5*(WVJJTree->vbf1_AK4_eta + WVJJTree->vbf2_AK4_eta);
      
      if(lheBr) {
	lheWgtArr->Clear();
	for (int j=0; j<lheWgtArr->GetEntries(); j++) {
	  const baconhep::TLHEWeight *lhe = (baconhep::TLHEWeight*)((*lheWgtArr)[j]);
	  WVJJTree->LHEWeight[i] = lhe->weight;
	}
      }
      
      //need to add photon part
      for (int j=0; j<AK4Arr->GetEntries(); j++) {
	const baconhep::TJet *ak4jet = (baconhep::TJet*)((*AK4Arr)[j]);
	if (abs(ak4jet->eta) > 2.0 && abs(ak4jet->eta)<3.0 && ak4jet->pt>30 && ak4jet->pt<500) {
	  //float pfRate = hL1PF_jetpt->GetBinContent(hL1PF_jetpt->FindBin(ak4jet->eta, ak4jet->pt));
	  WVJJTree->L1PFWeight *= scaleFactor.GetL1PFWeightJet(ak4jet->pt, ak4jet->eta);
	}
      }
      
      ot->Fill();
    }
    delete t; delete f;
    t=0; f=0;
  }
  of->Write();
  of->Close();
  
  delete f;
  f=0; t=0;
  
  return 0;

}
