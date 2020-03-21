#include "WWAnalysis/WWAnalysisRun2/interface/WVJJData.hh"

void WVJJData::init() {

  clearVars();

  //metadata and weights
  fTree->Branch("run",&run,"run/I");
  fTree->Branch("ls",&ls,"ls/I");
  fTree->Branch("evt",&evt,"evt/I");
  fTree->Branch("nPV",&nPV,"nPV/F");
  fTree->Branch("nPU_mean",&nPU_mean,"nPU_mean/F");
  fTree->Branch("genWeight",&genWeight,"genWeight/F");
  fTree->Branch("puWeight",&puWeight,"puWeight/F");
  fTree->Branch("puWeight_Up",&puWeight_Up,"puWeight_Up/F");
  fTree->Branch("puWeight_Dn",&puWeight_Dn,"puWeight_Dn/F");
  fTree->Branch("L1PFWeight",&L1PFWeight,"L1PFWeight/F");
  fTree->Branch("LHEWeight",&LHEWeight[0],"LHEWeight[1164]/F");
  //btag counters
  fTree->Branch("nBtag_loose",&nBtag_loose,"nBtag_loose/I");
  fTree->Branch("nBtag_medium",&nBtag_medium,"nBtag_medium/I");
  fTree->Branch("nBtag_tight",&nBtag_tight,"nBtag_tight/I");
  //lepton 1
  fTree->Branch("lep1_pt",&lep1_pt,"lep1_pt/F");
  fTree->Branch("lep1_eta",&lep1_eta,"lep1_eta/F");
  fTree->Branch("lep1_phi",&lep1_phi,"lep1_phi/F");
  fTree->Branch("lep1_m",&lep1_m,"lep1_m/F");
  fTree->Branch("lep1_q",&lep1_q,"lep1_q/F");
  fTree->Branch("lep1_iso",&lep1_iso,"lep1_iso/F");
  fTree->Branch("lep1_idEffWeight",&lep1_idEffWeight,"lep1_idEffWeight/F");
  //lepton 1 scale variations
  fTree->Branch("lep1_pt_scaleUp",&lep1_pt_scaleUp,"lep1_pt_scaleUp/F");
  fTree->Branch("lep1_pt_scaleDn",&lep1_pt_scaleDn,"lep1_pt_scaleDn/F");
  //lepton 2
  fTree->Branch("lep2_pt",&lep2_pt,"lep2_pt/F");
  fTree->Branch("lep2_eta",&lep2_eta,"lep2_eta/F");
  fTree->Branch("lep2_phi",&lep2_phi,"lep2_phi/F");
  fTree->Branch("lep2_m",&lep2_m,"lep2_m/F");
  fTree->Branch("lep2_q",&lep2_q,"lep2_q/F");
  fTree->Branch("lep2_iso",&lep2_iso,"lep2_iso/F");
  fTree->Branch("lep2_idEffWeight",&lep2_idEffWeight,"lep2_idEffWeight/F");  
  //lepton 2 scale variations
  fTree->Branch("lep2_pt_scaleUp",&lep2_pt_scaleUp,"lep2_pt_scaleUp/F");
  fTree->Branch("lep2_pt_scaleDn",&lep2_pt_scaleDn,"lep2_pt_scaleDn/F");
  //dilepton final state
  fTree->Branch("dilep_m",&dilep_m,"dilep_m/F");
  fTree->Branch("dilep_pt",&dilep_pt,"dilep_pt/F");
  fTree->Branch("dilep_eta",&dilep_eta,"dilep_eta/F");
  fTree->Branch("dilep_phi",&dilep_phi,"dilep_phi/F");
  //dilepton scale variations
  fTree->Branch("dilep_m_scaleUp",&dilep_m_scaleUp,"dilep_m_scaleUp/F");
  fTree->Branch("dilep_m_scaleDn",&dilep_m_scaleDn,"dilep_m_scaleDn/F");
  fTree->Branch("dilep_pt_scaleUp",&dilep_pt_scaleUp,"dilep_pt_scaleUp/F");
  fTree->Branch("dilep_pt_scaleDn",&dilep_pt_scaleDn,"dilep_pt_scaleDn/F");
  //MET
  fTree->Branch("MET",&MET,"MET/F");
  fTree->Branch("MET_phi",&MET_phi,"MET_phi/F");
  fTree->Branch("MET_2017raw",&MET_2017raw,"MET_2017raw/F");
  fTree->Branch("MET_scaleUp",&MET_scaleUp,"MET_scaleUp/F");
  fTree->Branch("MET_scaleDn",&MET_scaleDn,"MET_scaleDn/F");
  //W neutrino pZ
  fTree->Branch("neu_pz_type0",&neu_pz_type0,"neu_pz_type0/F");
  fTree->Branch("neu_pz_type0_scaleUp",&neu_pz_type0_scaleUp,"neu_pz_type0_scaleUp/F");
  fTree->Branch("neu_pz_type0_scaleDn",&neu_pz_type0_scaleDn,"neu_pz_type0_scaleDn/F");
  //VBF jet 1
  fTree->Branch("vbf1_AK4_pt",&vbf1_AK4_pt,"vbf1_AK4_pt/F");
  fTree->Branch("vbf1_AK4_eta",&vbf1_AK4_eta,"vbf1_AK4_eta/F");
  fTree->Branch("vbf1_AK4_phi",&vbf1_AK4_phi,"vbf1_AK4_phi/F");
  fTree->Branch("vbf1_AK4_m",&vbf1_AK4_m,"vbf1_AK4_m/F");
  fTree->Branch("vbf1_AK4_gqid",&vbf1_AK4_qgid,"vbf1_AK4_qgid/F");
  fTree->Branch("vbf1_AK4_axis2",&vbf1_AK4_axis2,"vbf1_AK4_axis2/F");
  fTree->Branch("vbf1_AK4_ptD",&vbf1_AK4_ptD,"vbf1_AK4_ptD/F");
  //VBF jet 1 variations
  fTree->Branch("vbf1_AK4_pt_scaleUp",&vbf1_AK4_pt_scaleUp,"vbf1_AK4_pt_scaleUp/F");
  fTree->Branch("vbf1_AK4_pt_scaleDn",&vbf1_AK4_pt_scaleDn,"vbf1_AK4_pt_scaleDn/F");
  fTree->Branch("vbf1_AK4_m_scaleUp",&vbf1_AK4_m_scaleUp,"vbf1_AK4_m_scaleUp/F");
  fTree->Branch("vbf1_AK4_m_scaleDn",&vbf1_AK4_m_scaleDn,"vbf1_AK4_m_scaleDn/F");
  //VBF jet 2
  fTree->Branch("vbf2_AK4_pt",&vbf2_AK4_pt,"vbf2_AK4_pt/F");
  fTree->Branch("vbf2_AK4_eta",&vbf2_AK4_eta,"vbf2_AK4_eta/F");
  fTree->Branch("vbf2_AK4_phi",&vbf2_AK4_phi,"vbf2_AK4_phi/F");
  fTree->Branch("vbf2_AK4_m",&vbf2_AK4_m,"vbf2_AK4_m/F");
  fTree->Branch("vbf2_AK4_gqid",&vbf2_AK4_qgid,"vbf2_AK4_qgid/F");
  fTree->Branch("vbf2_AK4_axis2",&vbf2_AK4_axis2,"vbf2_AK4_axis2/F");
  fTree->Branch("vbf2_AK4_ptD",&vbf2_AK4_ptD,"vbf2_AK4_ptD/F");
  //VBF jet 2 variations
  fTree->Branch("vbf2_AK4_pt_scaleUp",&vbf2_AK4_pt_scaleUp,"vbf2_AK4_pt_scaleUp/F");
  fTree->Branch("vbf2_AK4_pt_scaleDn",&vbf2_AK4_pt_scaleDn,"vbf2_AK4_pt_scaleDn/F");
  fTree->Branch("vbf2_AK4_m_scaleUp",&vbf2_AK4_m_scaleUp,"vbf2_AK4_m_scaleUp/F");
  fTree->Branch("vbf2_AK4_m_scaleDn",&vbf2_AK4_m_scaleDn,"vbf2_AK4_m_scaleDn/F");
  //VBF dijet object
  fTree->Branch("vbf_pt", &vbf_pt, "vbf_pt/F");
  fTree->Branch("vbf_eta",&vbf_eta,"vbf_eta/F");
  fTree->Branch("vbf_phi",&vbf_phi,"vbf_phi/F");
  fTree->Branch("vbf_m",  &vbf_m,  "vbf_m/F");
  //VBF dijet variations
  fTree->Branch("vbf_pt_scaleUp", &vbf_pt_scaleUp, "vbf_pt_scaleUp/F");
  fTree->Branch("vbf_pt_scaleDn", &vbf_pt_scaleDn, "vbf_pt_scaleDn/F");
  fTree->Branch("vbf_m_scaleUp", &vbf_m_scaleUp, "vbf_m_scaleUp/F");
  fTree->Branch("vbf_m_scaleDn", &vbf_m_scaleDn, "vbf_m_scaleDn/F");
  //Boson AK8 jet
  fTree->Branch("bos_PuppiAK8_m_sd0",&bos_PuppiAK8_m_sd0,"bos_PuppiAK8_m_sd0/F");
  fTree->Branch("bos_PuppiAK8_m_sd0_corr",&bos_PuppiAK8_m_sd0_corr,"bos_PuppiAK8_m_sd0_corr/F");
  fTree->Branch("bos_PuppiAK8_pt",&bos_PuppiAK8_pt,"bos_PuppiAK8_pt/F");
  fTree->Branch("bos_PuppiAK8_eta",&bos_PuppiAK8_eta,"bos_PuppiAK8_eta/F");
  fTree->Branch("bos_PuppiAK8_phi",&bos_PuppiAK8_phi,"bos_PuppiAK8_phi/F");
  fTree->Branch("bos_PuppiAK8_tau2tau1",&bos_PuppiAK8_tau2tau1,"bos_PuppiAK8_tau2tau1/F");
  //Boson AK8 jet variations
  fTree->Branch("bos_PuppiAK8_m_sd0_corr_scaleUp",&bos_PuppiAK8_m_sd0_corr_scaleUp,"bos_PuppiAK8_m_sd0_corr_scaleUp/F");
  fTree->Branch("bos_PuppiAK8_m_sd0_corr_scaleDn",&bos_PuppiAK8_m_sd0_corr_scaleDn,"bos_PuppiAK8_m_sd0_corr_scaleDn/F");
  fTree->Branch("bos_PuppiAK8_pt_scaleUp",&bos_PuppiAK8_pt_scaleUp,"bos_PuppiAK8_pt_scaleUp/F");
  fTree->Branch("bos_PuppiAK8_pt_scaleDn",&bos_PuppiAK8_pt_scaleDn,"bos_PuppiAK8_pt_scaleDn/F");
  //Boson AK8 correlation variations
  fTree->Branch("bos_PuppiAK8_e2_sdb1", &bos_PuppiAK8_e2_sdb1, "bos_PuppiAK8_e2_sdb1/F");
  fTree->Branch("bos_PuppiAK8_e3_sdb1", &bos_PuppiAK8_e3_sdb1, "bos_PuppiAK8_e3_sdb1/F");
  fTree->Branch("bos_PuppiAK8_e3_v1_sdb1", &bos_PuppiAK8_e3_v1_sdb1, "bos_PuppiAK8_e3_v1_sdb1/F");
  fTree->Branch("bos_PuppiAK8_e3_v2_sdb1", &bos_PuppiAK8_e3_v2_sdb1, "bos_PuppiAK8_e3_v2_sdb1/F");
  fTree->Branch("bos_PuppiAK8_e4_v1_sdb1", &bos_PuppiAK8_e4_v1_sdb1, "bos_PuppiAK8_e4_v1_sdb1/F");
  fTree->Branch("bos_PuppiAK8_e4_v2_sdb1", &bos_PuppiAK8_e4_v2_sdb1, "bos_PuppiAK8_e4_v2_sdb1/F");
  fTree->Branch("bos_PuppiAK8_e2_sdb2", &bos_PuppiAK8_e2_sdb2, "bos_PuppiAK8_e2_sdb2/F");
  fTree->Branch("bos_PuppiAK8_e3_sdb2", &bos_PuppiAK8_e3_sdb2, "bos_PuppiAK8_e3_sdb2/F");
  fTree->Branch("bos_PuppiAK8_e3_v1_sdb2", &bos_PuppiAK8_e3_v1_sdb2, "bos_PuppiAK8_e3_v1_sdb2/F");
  fTree->Branch("bos_PuppiAK8_e3_v2_sdb2", &bos_PuppiAK8_e3_v2_sdb2, "bos_PuppiAK8_e3_v2_sdb2/F");
  fTree->Branch("bos_PuppiAK8_e4_v1_sdb2", &bos_PuppiAK8_e4_v1_sdb2, "bos_PuppiAK8_e4_v1_sdb2/F");
  fTree->Branch("bos_PuppiAK8_e4_v2_sdb2", &bos_PuppiAK8_e4_v2_sdb2, "bos_PuppiAK8_e4_v2_sdb2/F");
  //Boson AK4 jet 1
  fTree->Branch("bos_j1_AK4_pt",&bos_j1_AK4_pt,"bos_j1_AK4_pt/F");
  fTree->Branch("bos_j1_AK4_eta",&bos_j1_AK4_eta,"bos_j1_AK4_eta/F");
  fTree->Branch("bos_j1_AK4_phi",&bos_j1_AK4_phi,"bos_j1_AK4_phi/F");
  fTree->Branch("bos_j1_AK4_m",&bos_j1_AK4_m,"bos_j1_AK4_m/F");
  //Boson AK4 jet 1 variations
  fTree->Branch("bos_j1_AK4_pt_scaleUp",&bos_j1_AK4_pt_scaleUp,"bos_j1_AK4_pt_scaleUp/F");
  fTree->Branch("bos_j1_AK4_pt_scaleDn",&bos_j1_AK4_pt_scaleDn,"bos_j1_AK4_pt_scaleDn/F");
  fTree->Branch("bos_j1_AK4_m_scaleUp",&bos_j1_AK4_m_scaleUp,"bos_j1_AK4_m_scaleUp/F");
  fTree->Branch("bos_j1_AK4_m_scaleDn",&bos_j1_AK4_m_scaleDn,"bos_j1_AK4_m_scaleDn/F");
  //Boson AK4 jet 2
  fTree->Branch("bos_j2_AK4_pt",&bos_j2_AK4_pt,"bos_j2_AK4_pt/F");
  fTree->Branch("bos_j2_AK4_eta",&bos_j2_AK4_eta,"bos_j2_AK4_eta/F");
  fTree->Branch("bos_j2_AK4_phi",&bos_j2_AK4_phi,"bos_j2_AK4_phi/F");
  fTree->Branch("bos_j2_AK4_m",&bos_j2_AK4_m,"bos_j2_AK4_m/F");
  //Boson AK4 jet 2 variations
  fTree->Branch("bos_j2_AK4_pt_scaleUp",&bos_j2_AK4_pt_scaleUp,"bos_j2_AK4_pt_scaleUp/F");
  fTree->Branch("bos_j2_AK4_pt_scaleDn",&bos_j2_AK4_pt_scaleDn,"bos_j2_AK4_pt_scaleDn/F");
  fTree->Branch("bos_j2_AK4_m_scaleUp",&bos_j2_AK4_m_scaleUp,"bos_j2_AK4_m_scaleUp/F");
  fTree->Branch("bos_j2_AK4_m_scaleDn",&bos_j2_AK4_m_scaleDn,"bos_j2_AK4_m_scaleDn/F");
  //Boson dijet object
  fTree->Branch("bos_AK4AK4_pt", &bos_AK4AK4_pt, "bos_AK4AK4_pt/F");
  fTree->Branch("bos_AK4AK4_eta",&bos_AK4AK4_eta,"bos_AK4AK4_eta/F");
  fTree->Branch("bos_AK4AK4_phi",&bos_AK4AK4_phi,"bos_AK4AK4_phi/F");
  fTree->Branch("bos_AK4AK4_m",  &bos_AK4AK4_m,  "bos_AK4AK4_m/F");
  //Boson dijet variations
  fTree->Branch("bos_AK4AK4_pt_scaleUp", &bos_AK4AK4_pt_scaleUp, "bos_AK4AK4_pt_scaleUp/F");
  fTree->Branch("bos_AK4AK4_pt_scaleDn", &bos_AK4AK4_pt_scaleDn, "bos_AK4AK4_pt_scaleDn/F");
  fTree->Branch("bos_AK4AK4_m_scaleUp", &bos_AK4AK4_m_scaleUp, "bos_AK4AK4_m_scaleUp/F");
  fTree->Branch("bos_AK4AK4_m_scaleDn", &bos_AK4AK4_m_scaleDn, "bos_AK4AK4_m_scaleDn/F");
  //final state variables
  fTree->Branch("dibos_m", &dibos_m, "dibos_m/F");
  fTree->Branch("dibos_pt", &dibos_pt, "dibos_pt/F");
  fTree->Branch("dibos_eta", &dibos_eta, "dibos_eta/F");
  fTree->Branch("dibos_phi", &dibos_phi, "dibos_phi/F");

  fTree->Branch("dibos_m_scaleUp", &dibos_m_scaleUp, "dibos_m_scaleUp/F");
  fTree->Branch("dibos_m_scaleDn", &dibos_m_scaleDn, "dibos_m_scaleDn/F");
  fTree->Branch("dibos_pt_scaleUp", &dibos_pt_scaleUp, "dibos_pt_scaleUp/F");
  fTree->Branch("dibos_pt_scaleDn", &dibos_pt_scaleDn, "dibos_pt_scaleDn/F");
  
  fTree->Branch("bosCent", &bosCent, "bosCent/F");
  fTree->Branch("zeppLep", &zeppLep, "zeppLep/F");
  fTree->Branch("zeppHad", &zeppHad, "zeppHad/F");

};

void WVJJData::clearVars() {

  //------------------------------------//
  //       METADATA AND EVENT WEIGHTS   //
  //------------------------------------//
  
  run = 0;
  ls = 0;
  evt = 0;

  nPV = 0;
  nPU_mean = 0;
  
  genWeight = 1.0;
  puWeight = 1.0;
  puWeight_Up = 1.0;
  puWeight_Dn = 1.0;

  L1PFWeight = 1.0;

  LHEWeight[1164] = {};

  nBtag_loose = 0;
  nBtag_medium = 0;
  nBtag_tight = 0;
  
  //------------------------------------//
  //       LEPTONS                      //
  //------------------------------------//

  //lepton 1
  lep1_pt = -999.0;
  lep1_eta = -999.0;
  lep1_phi = -999.0;
  lep1_m = -999.0;
  lep1_q = -999.0;
  lep1_iso = -999.0;
  lep1_idEffWeight = 1.0;

  //lepton 1 scale variations
  lep1_pt_scaleUp = -999.0;
  lep1_pt_scaleDn = -999.0;
  
  //lepton 2
  lep2_pt = -999.0;
  lep2_eta = -999.0;
  lep2_phi = -999.0;
  lep2_m = -999.0;
  lep2_q = -999.0;
  lep2_iso = -999.0;
  lep2_idEffWeight = 1.0;

  //lepton 2 scale variations
  lep2_pt_scaleUp = -999.0;
  lep2_pt_scaleDn = -999.0;

  //dilepton final state
  dilep_m = -999.0;
  dilep_pt = -999.0;
  dilep_eta = -999.0;
  dilep_phi = -999.0;

  //dilepton scale variations
  dilep_m_scaleUp = -999.0;
  dilep_m_scaleDn = -999.0;
  dilep_pt_scaleUp = -999.0;
  dilep_pt_scaleDn = -999.0;

  //------------------------------------//
  //       MET                          //
  //------------------------------------//

  MET = -999.0;
  MET_phi = -999.0;

  MET_2017raw = -999.0;

  MET_scaleUp = -999.0;
  MET_scaleDn = -999.0;

  neu_pz_type0 = -999.0;
  neu_pz_type0_scaleUp = -999.0;
  neu_pz_type0_scaleDn = -999.0;

  //------------------------------------//
  //       VBF/TAGGING JETS             //
  //------------------------------------//

  //VBF jet 1
  vbf1_AK4_pt = -999.0;
  vbf1_AK4_eta = -999.0;
  vbf1_AK4_phi = -999.0;
  vbf1_AK4_m = -999.0;
  vbf1_AK4_qgid = -999.0;
  vbf1_AK4_axis2 = -999.0;
  vbf1_AK4_ptD = -999.0;

  //VBF jet 1 variations
  vbf1_AK4_pt_scaleUp = -999.0;
  vbf1_AK4_pt_scaleDn = -999.0;
  vbf1_AK4_m_scaleUp = -999.0;
  vbf1_AK4_m_scaleDn = -999.0;

  //VBF jet 2
  vbf2_AK4_pt = -999.0;
  vbf2_AK4_eta = -999.0;
  vbf2_AK4_phi = -999.0;
  vbf2_AK4_m = -999.0;
  vbf2_AK4_qgid = -999.0;
  vbf2_AK4_axis2 = -999.0;
  vbf2_AK4_ptD = -999.0;

  //VBF jet 2 variations
  vbf2_AK4_pt_scaleUp = -999.0;
  vbf2_AK4_pt_scaleDn = -999.0;
  vbf2_AK4_m_scaleUp = -999.0;
  vbf2_AK4_m_scaleDn = -999.0;

  //VBF dijet object
  vbf_pt = -999.0;
  vbf_eta = -999.0;
  vbf_phi = -999.0;
  vbf_m = -999.0;

  //VBF dijet variations
  vbf_pt_scaleUp = -999.0;
  vbf_pt_scaleDn = -999.0;
  vbf_m_scaleUp = -999.0;
  vbf_m_scaleDn = -999.0;

  //------------------------------------//
  //       HADRONIC BOOSTED OBJECTS     //
  //------------------------------------//

  //Boson AK8 jet
  bos_PuppiAK8_m_sd0 = -999.0;
  bos_PuppiAK8_m_sd0_corr = -999.0;
  bos_PuppiAK8_pt = -999.0;
  bos_PuppiAK8_eta = -999.0;
  bos_PuppiAK8_phi = -999.0;
  bos_PuppiAK8_tau2tau1 = -999.0;

  //Boson AK8 jet variations
  bos_PuppiAK8_m_sd0_corr_scaleUp = -999.0;
  bos_PuppiAK8_m_sd0_corr_scaleDn = -999.0;
  bos_PuppiAK8_pt_scaleUp = -999.0;
  bos_PuppiAK8_pt_scaleDn = -999.0;

  bos_PuppiAK8_e2_sdb1 = -999.0;
  bos_PuppiAK8_e3_sdb1 = -999.0;
  bos_PuppiAK8_e3_v1_sdb1 = -999.0;
  bos_PuppiAK8_e3_v2_sdb1 = -999.0;
  bos_PuppiAK8_e4_v1_sdb1 = -999.0;
  bos_PuppiAK8_e4_v2_sdb1 = -999.0;

  bos_PuppiAK8_e2_sdb2 = -999.0;
  bos_PuppiAK8_e3_sdb2 = -999.0;
  bos_PuppiAK8_e3_v1_sdb2 = -999.0;
  bos_PuppiAK8_e3_v2_sdb2 = -999.0;
  bos_PuppiAK8_e4_v1_sdb2 = -999.0;
  bos_PuppiAK8_e4_v2_sdb2 = -999.0;

  //------------------------------------//
  //       HADRONIC RESOLVED OBJECTS    //
  //------------------------------------//

  //Boson AK4 jet 1
  bos_j1_AK4_pt = -999.0;
  bos_j1_AK4_eta = -999.0;
  bos_j1_AK4_phi = -999.0;
  bos_j1_AK4_m = -999.0;

  //Boson AK4 jet 1 variations
  bos_j1_AK4_pt_scaleUp = -999.0;
  bos_j1_AK4_pt_scaleDn = -999.0;
  bos_j1_AK4_m_scaleUp = -999.0;
  bos_j1_AK4_m_scaleDn = -999.0;

  //Boson AK4 jet 2
  bos_j2_AK4_pt = -999.0;
  bos_j2_AK4_eta = -999.0;
  bos_j2_AK4_phi = -999.0;
  bos_j2_AK4_m = -999.0;

  //Boson AK4 jet 2 variations
  bos_j2_AK4_pt_scaleUp = -999.0;
  bos_j2_AK4_pt_scaleDn = -999.0;
  bos_j2_AK4_m_scaleUp = -999.0;
  bos_j2_AK4_m_scaleDn = -999.0;

  //Boson dijet object
  bos_AK4AK4_pt = -999.0;
  bos_AK4AK4_eta = -999.0;
  bos_AK4AK4_phi = -999.0;
  bos_AK4AK4_m = -999.0;

  dibos_m = -999.0;
  dibos_pt = -999.0;
  dibos_eta = -999.0;
  dibos_phi = -999.0;

  dibos_m_scaleUp = -999.0;
  dibos_m_scaleDn = -999.0;
  dibos_pt_scaleUp = -999.0;
  dibos_pt_scaleDn = -999.0;

  bosCent = -999.0;
  zeppLep = -999.0;
  zeppHad = -999.0;

};

