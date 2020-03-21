#include "WWAnalysis/WWAnalysisRun2/interface/ScaleFactors.hh"

void ScaleFactors::setupSFs() {

  if (era_==2016) {
    
    // pileup
    
    //nBinsPU=75;
    pileupFileMC = TFile::Open("data/puWeights_80x_37ifb.root");

    // electrons

    IDIsoEle = TFile::Open("data/SF2016/egammaEffi_EGM2D_TightCutBasedIDSF.root","READ");
    hIDIsoEle = (TH1F*)IDIsoEle->Get("EGamma_SF2D");

    GSFCorrEle = TFile::Open("data/SF2016/egammaEffi_SF2D_GSF_tracking.root","READ");
    hGSFCorrEle = (TH1F*)GSFCorrEle->Get("EGamma_SF2D");

    TriggerEle = TFile::Open("data/SF2016/ElectronTrigger_SF.root","READ");
    hTriggerEle = (TH1F*)TriggerEle->Get("HLT_Ele27");

    // muons

    IDMuA = TFile::Open("data/SF2016/MuonID_RunBCDEF_23SepReReco_19p72fb.root","READ");
    hIDMuA = (TH1F*)IDMuA->Get("MC_NUM_TightID_DEN_genTracks_PAR_pt_eta/abseta_pt_ratio");
    IDMuB = TFile::Open("data/SF2016/MuonID_RunGH_23SepReReco_16p146fb.root","READ");
    hIDMuB = (TH1F*)IDMuB->Get("MC_NUM_TightID_DEN_genTracks_PAR_pt_eta/abseta_pt_ratio");

    IsoMuA = TFile::Open("data/SF2016/MuonIso_RunBCDEF_23SepReReco_19p72fb.root","READ");
    hIsoMuA = (TH1F*)IsoMuA->Get("TightISO_TightID_pt_eta/abseta_pt_ratio");
    IsoMuB = TFile::Open("data/SF2016/MuonIso_RunGH_23SepReReco_16p146fb.root","READ");
    hIsoMuB = (TH1F*)IsoMuB->Get("TightISO_TightID_pt_eta/abseta_pt_ratio");

    TriggerMuA = TFile::Open("data/SF2016/MuonTrigger_RunBCDEF_23SepReReco_19p72fb.root","READ");
    hTriggerMuA = (TH1F*)TriggerMuA->Get("IsoMu24_OR_IsoTkMu24_PtEtaBins/abseta_pt_ratio");
    TriggerMuB = TFile::Open("data/SF2016/MuonTrigger_RunGH_23SepReReco_16p146fb.root","READ");
    hTriggerMuB = (TH1F*)TriggerMuB->Get("IsoMu24_OR_IsoTkMu24_PtEtaBins/abseta_pt_ratio");

    btagCSV="data/CSVv2_Moriond17_B_H.csv";
    //calib(TString("csvv2"), TString("data/CSVv2_Moriond17_B_H.csv"));

    jecUncAK4="data/Summer16_23Sep2016V4_MC/Summer16_23Sep2016V4_MC_Uncertainty_AK4PFchs.txt";
    jecUncAK8="data/Summer16_23Sep2016V4_MC/Summer16_23Sep2016V4_MC_Uncertainty_AK8PFPuppi.txt";

    //paramAK4chs("data/Summer16_23Sep2016V4_MC/Summer16_23Sep2016V4_MC_Uncertainty_AK4PFchs.txt");
    //paramAK8puppi("data/Summer16_23Sep2016V4_MC/Summer16_23Sep2016V4_MC_Uncertainty_AK8PFPuppi.txt");

    L1PF_jetpt = TFile::Open("data/SF2016/L1prefiring_jetpt_2016BtoH.root");
    hL1PF_jetpt = (TH2F*) L1PF_jetpt->Get("L1prefiring_jetpt_2016BtoH");

    L1PF_photpt = TFile::Open("data/SF2016/L1prefiring_photonpt_2016BtoH.root");
    hL1PF_photpt = (TH2F*) L1PF_photpt->Get("L1prefiring_photonpt_2016BtoH");

    eleTriggerNames.push_back("HLT_Ele27_WPTight_Gsf_v*");
    muTriggerNames.push_back("HLT_IsoMu24_v*");
    muTriggerNames.push_back("HLT_IsoTkMu24_v*");

  } else if (era_==2017) {
    
    // pileup
    
    //nBinsPU=100;
    pileupFileMC = TFile::Open("data/puWeights_90x_2017.root");

    // electrons

    IDIsoEle = TFile::Open("data/SF2017/egammaEffi.txt_EGM2D_runBCDEF_passingTight94X.root","READ");
    hIDIsoEle = (TH1F*)IDIsoEle->Get("EGamma_SF2D");

    GSFCorrEle = TFile::Open("data/SF2017/egammaEffi.txt_EGM2D_runBCDEF_passingRECO.root","READ");
    hGSFCorrEle = (TH1F*)GSFCorrEle->Get("EGamma_SF2D");

    //TriggerEle = TFile::Open("data/SF2016/ElectronTrigger_SF.root","READ");
    //hTriggerEle = (TH1F*)TriggerEle->Get("HLT_Ele27");

    // muons

    //IDMuA = TFile::Open("data/SF2016/MuonID_RunBCDEF_23SepReReco_19p72fb.root","READ");
    //hIDMuA = (TH1F*)IDMuA->Get("MC_NUM_TightID_DEN_genTracks_PAR_pt_eta/abseta_pt_ratio");
    //IDMuB = TFile::Open("data/SF2016/MuonID_RunGH_23SepReReco_16p146fb.root","READ");
    //hIDMuB = (TH1F*)IDMuB->Get("MC_NUM_TightID_DEN_genTracks_PAR_pt_eta/abseta_pt_ratio");

    //IsoMuA = TFile::Open("data/SF2016/MuonIso_RunBCDEF_23SepReReco_19p72fb.root","READ");
    //hIsoMuA = (TH1F*)IsoMuA->Get("TightISO_TightID_pt_eta/abseta_pt_ratio");
    //IsoMuB = TFile::Open("data/SF2016/MuonIso_RunGH_23SepReReco_16p146fb.root","READ");
    //hIsoMuB = (TH1F*)IsoMuB->Get("TightISO_TightID_pt_eta/abseta_pt_ratio");

    //TriggerMuA = TFile::Open("data/SF2016/MuonTrigger_RunBCDEF_23SepReReco_19p72fb.root","READ");
    //hTriggerMuA = (TH1F*)TriggerMuA->Get("IsoMu24_OR_IsoTkMu24_PtEtaBins/abseta_pt_ratio");
    //TriggerMuB = TFile::Open("data/SF2016/MuonTrigger_RunGH_23SepReReco_16p146fb.root","READ");
    //hTriggerMuB = (TH1F*)TriggerMuB->Get("IsoMu24_OR_IsoTkMu24_PtEtaBins/abseta_pt_ratio");

    // b-tagging

    btagCSV="data/CSVv2_94XSF_V2_B_F.csv";
    //calib("csvv2", "data/CSVv2_94XSF_V2_B_F.csv");

    // jec

    jecUncAK4="data/Fall17_17Nov2017_V32_MC/Fall17_17Nov2017_V32_MC_Uncertainty_AK4PFchs.txt";
    jecUncAK8="data/Fall17_17Nov2017_V32_MC/Fall17_17Nov2017_V32_MC_Uncertainty_AK8PFPuppi.txt";

    //paramAK4chs("data/Fall17_17Nov2017_V32_MC/Fall17_17Nov2017_V32_MC_Uncertainty_AK4PFchs.txt");
    //paramAK8puppi("data/Fall17_17Nov2017_V32_MC/Fall17_17Nov2017_V32_MC_Uncertainty_AK8PFPuppi.txt");

    L1PF_jetpt = TFile::Open("data/SF2017/L1prefiring_jetpt_2017BtoF.root");
    hL1PF_jetpt = (TH2F*) L1PF_jetpt->Get("L1prefiring_jetpt_2017BtoF");

    L1PF_photpt = TFile::Open("data/SF2017/L1prefiring_photonpt_2017BtoF.root");
    hL1PF_photpt = (TH2F*) L1PF_photpt->Get("L1prefiring_photonpt_2017BtoF");

    eleTriggerNames.push_back("HLT_Ele32_WPTight_Gsf_v*");
    muTriggerNames.push_back("HLT_IsoMu27_v*");

  }

  // pileup
  
  puWeights = (TH1D*)pileupFileMC->Get("puWeights");
  puWeightsUp = (TH1D*)pileupFileMC->Get("puWeightsUp");
  puWeightsDown = (TH1D*)pileupFileMC->Get("puWeightsDown");
  
  //puWeights->SetBins(nbinsPU,0,nbinsPU);
  //puWeightsUp->SetBins(nbinsPU,0,nbinsPU);
  //puWeightsDown->SetBins(nbinsPU,0,nbinsPU);
  
  // puppi

  fPuppi = TFile::Open("data/puppiCorr.root","READ");
  puppisd_corrGEN      = (TF1*)fPuppi->Get("puppiJECcorr_gen");
  puppisd_corrRECO_cen = (TF1*)fPuppi->Get("puppiJECcorr_reco_0eta1v3");
  puppisd_corrRECO_for = (TF1*)fPuppi->Get("puppiJECcorr_reco_1v3eta2v5");
  
  
  // b-tagging
  
  //bTagReader(BTagEntry::OP_LOOSE,  // working point: can be OP_LOOSE, OP_MEDIUM, OP_TIGHT                              
  //	     "central",             // label for the central value (see the scale factor file)                         
  //	     {"up","down"});        // vector of labels for systematics                                                
  //
  //bTagReader.load(calib, BTagEntry::FLAV_B, "comb");      // use the "comb" measurements for b-jets                    
  //bTagReader.load(calib, BTagEntry::FLAV_C, "comb");      // use the "comb" measurements for c-jets                    
  //bTagReader.load(calib, BTagEntry::FLAV_UDSG, "incl");   // use the "incl" measurements for light jets
  
  // jec
  
  //fJetUnc_AK4chs = new JetCorrectionUncertainty(paramAK4chs);
  //fJetUnc_AK8puppi = new JetCorrectionUncertainty(paramAK8puppi);


}
