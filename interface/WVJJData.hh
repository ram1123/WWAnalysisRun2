#ifndef WVJJ_DATA_HH
#define WVJJ_DATA_HH

#include "TTree.h"

class WVJJData {

public:

  WVJJData(TTree *t) {
    fTree = t;
    init();
  };
  
  ~WVJJData() {
    delete fTree;
  };

  void clearVars();
  void init();  
  
  TTree *fTree;

  //------------------------------------//
  //       METADATA AND EVENT WEIGHTS   //
  //------------------------------------//
  
  uint run;
  uint ls;
  uint evt;
  float nPV;
  float nPU_mean;
  float genWeight;
  float puWeight;
  float puWeight_Up;
  float puWeight_Dn;
  float L1PFWeight;
  float LHEWeight[1164] = {};

  //btag counters
  uint nBtag_loose;
  uint nBtag_medium;
  uint nBtag_tight;
  
  //------------------------------------//
  //       LEPTONS                      //
  //------------------------------------//

  //lepton 1
  float lep1_pt;
  float lep1_eta;
  float lep1_phi;
  float lep1_m;
  float lep1_q;
  float lep1_iso;
  float lep1_idEffWeight;

  //lepton 1 scale variations
  float lep1_pt_scaleUp;
  float lep1_pt_scaleDn;
  
  //lepton 2
  float lep2_pt;
  float lep2_eta;
  float lep2_phi;
  float lep2_m;
  float lep2_q;
  float lep2_iso;
  float lep2_idEffWeight;

  //lepton 2 scale variations
  float lep2_pt_scaleUp;
  float lep2_pt_scaleDn;

  //dilepton final state
  float dilep_m;
  float dilep_pt;
  float dilep_eta;
  float dilep_phi;

  //dilepton scale variations
  float dilep_m_scaleUp;
  float dilep_m_scaleDn;
  float dilep_pt_scaleUp;
  float dilep_pt_scaleDn;

  //------------------------------------//
  //       MET                          //
  //------------------------------------//

  float MET;
  float MET_phi;
  float MET_2017raw;
  float MET_scaleUp;
  float MET_scaleDn;

  //W neutrino pZ
  float neu_pz_type0;
  float neu_pz_type0_scaleUp;
  float neu_pz_type0_scaleDn;

  //------------------------------------//
  //       VBF/TAGGING JETS             //
  //------------------------------------//

  //VBF jet 1
  float vbf1_AK4_pt;
  float vbf1_AK4_eta;
  float vbf1_AK4_phi;
  float vbf1_AK4_m;
  //float vbf1_AK4_e;
  float vbf1_AK4_qgid;
  float vbf1_AK4_axis2;
  float vbf1_AK4_ptD;

  //VBF jet 1 variations
  float vbf1_AK4_pt_scaleUp;
  float vbf1_AK4_pt_scaleDn;
  float vbf1_AK4_m_scaleUp;
  float vbf1_AK4_m_scaleDn;
  //float vbf1_AK4_e_scaleUp;
  //float vbf1_AK4_e_scaleDn;

  //VBF jet 2
  float vbf2_AK4_pt;
  float vbf2_AK4_eta;
  float vbf2_AK4_phi;
  float vbf2_AK4_m;
  //float vbf2_AK4_e;
  float vbf2_AK4_qgid;
  float vbf2_AK4_axis2;
  float vbf2_AK4_ptD;

  //VBF jet 2 variations
  float vbf2_AK4_pt_scaleUp;
  float vbf2_AK4_pt_scaleDn;
  float vbf2_AK4_m_scaleUp;
  float vbf2_AK4_m_scaleDn;
  //float vbf2_AK4_e_scaleUp;
  //float vbf2_AK4_e_scaleDn;

  //VBF dijet object
  float vbf_pt;
  float vbf_eta;
  float vbf_phi;
  float vbf_m;
  float vbf_deta;

  //VBF dijet variations
  float vbf_pt_scaleUp;
  float vbf_pt_scaleDn;
  float vbf_m_scaleUp;
  float vbf_m_scaleDn;

  //------------------------------------//
  //       HADRONIC BOOSTED OBJECTS     //
  //------------------------------------//

  //Boson AK8 jet
  float bos_PuppiAK8_m_sd0;
  float bos_PuppiAK8_m_sd0_corr;
  float bos_PuppiAK8_pt;
  float bos_PuppiAK8_eta;
  float bos_PuppiAK8_phi;
  //float bos_PuppiAK8_e_ungroomed;
  float bos_PuppiAK8_tau2tau1;

  //Boson AK8 jet variations
  float bos_PuppiAK8_m_sd0_corr_scaleUp;
  float bos_PuppiAK8_m_sd0_corr_scaleDn;
  float bos_PuppiAK8_pt_scaleUp;
  float bos_PuppiAK8_pt_scaleDn;
  //float bos_PuppiAK8_e_ungroomed_scaleUp;
  //float bos_PuppiAK8_e_ungroomed_scaleDn;

  float bos_PuppiAK8_e2_sdb1; 
  float bos_PuppiAK8_e3_sdb1; 
  float bos_PuppiAK8_e3_v1_sdb1; 
  float bos_PuppiAK8_e3_v2_sdb1; 
  float bos_PuppiAK8_e4_v1_sdb1; 
  float bos_PuppiAK8_e4_v2_sdb1;
  float bos_PuppiAK8_e2_sdb2; 
  float bos_PuppiAK8_e3_sdb2; 
  float bos_PuppiAK8_e3_v1_sdb2; 
  float bos_PuppiAK8_e3_v2_sdb2; 
  float bos_PuppiAK8_e4_v1_sdb2; 
  float bos_PuppiAK8_e4_v2_sdb2; 

  //------------------------------------//
  //       HADRONIC RESOLVED OBJECTS    //
  //------------------------------------//

  //Boson AK4 jet 1
  float bos_j1_AK4_pt;
  float bos_j1_AK4_eta;
  float bos_j1_AK4_phi;
  float bos_j1_AK4_m;
  //float bos_j1_AK4_e;

  //Boson AK4 jet 1 variations
  float bos_j1_AK4_pt_scaleUp;
  float bos_j1_AK4_pt_scaleDn;
  float bos_j1_AK4_m_scaleUp;
  float bos_j1_AK4_m_scaleDn;
  //float bos_j1_AK4_e_scaleUp;
  //float bos_j1_AK4_e_scaleDn;

  //Boson AK4 jet 2
  float bos_j2_AK4_pt;
  float bos_j2_AK4_eta;
  float bos_j2_AK4_phi;
  float bos_j2_AK4_m;
  //float bos_j2_AK4_e;

  //Boson AK4 jet 2 variations
  float bos_j2_AK4_pt_scaleUp;
  float bos_j2_AK4_pt_scaleDn;
  float bos_j2_AK4_m_scaleUp;
  float bos_j2_AK4_m_scaleDn;
  //float bos_j2_AK4_e_scaleUp;
  //float bos_j2_AK4_e_scaleDn;

  //Boson dijet object
  float bos_AK4AK4_pt;
  float bos_AK4AK4_eta;
  float bos_AK4AK4_phi;
  float bos_AK4AK4_m;

  //Boson dijet variations
  float bos_AK4AK4_pt_scaleUp;
  float bos_AK4AK4_pt_scaleDn;
  float bos_AK4AK4_m_scaleUp;
  float bos_AK4AK4_m_scaleDn;

  //------------------------------------//
  //       FINAL STATE VARIABLES        //
  //------------------------------------//

  float dibos_m;
  float dibos_pt;
  float dibos_eta;
  float dibos_phi;

  float dibos_m_scaleUp;
  float dibos_m_scaleDn;
  float dibos_pt_scaleUp;
  float dibos_pt_scaleDn;
  
  float bosCent;
  float zeppLep;
  float zeppHad;

};

#endif
