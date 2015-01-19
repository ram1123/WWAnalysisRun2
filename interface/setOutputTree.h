#ifndef __setOutputTree__
#define __setOutputTree__

#include "TTree.h"
#include "TChain.h"

// Declaration of leaf types
extern int run;
extern int event;
extern int nJets;
extern int nVtx;
extern float met;
extern float met_px;
extern float met_py;
extern float met_pz_type0;
extern float met_pz_type2;
extern float leptonPt;
extern float leptonEta;
extern float leptonPhi;
extern float leptonE;
extern float AK8jetPt[10];
extern float AK8jetEta[10];
extern float AK8jetPhi[10];
extern float AK8jetE[10];
extern float AK8jetPrunedMass[10];
extern float AK8jetTrimmedMass[10];
extern float AK8jetFilteredMass[10];
extern float AK8jetTau1[10];
extern float AK8jetTau2[10];
extern float AK8jetTau3[10];
extern float jetPt[10];
extern float jetEta[10];
extern float jetPhi[10];
extern float jetE[10];
extern float jet_bDiscr[10];
extern int genBosonPdgId[10];
extern float genBosonPt[10];
extern float genBosonEta[10];
extern float genBosonPhi[10];
extern float genBosonE[10];
extern int genLeptonPdgId[10];
extern float genLeptonPt[10];
extern float genLeptonEta[10];
extern float genLeptonPhi[10];
extern float genLeptonE[10];
extern float genNuPt[10];
extern float genNuEta[10];
extern float genNuPhi[10];
extern float genNuE[10];
extern float deltaR_lak8jet;
extern float deltaphi_METak8jet;
extern float deltaphi_Vak8jet;
extern float W_pt;
extern float W_eta;
extern float W_phi;
extern float W_E;
extern float W_mt;

// List of branches
extern TBranch *b_run;
extern TBranch *b_event;
extern TBranch *b_nJets;
extern TBranch *b_nVtx;
extern TBranch *b_met;
extern TBranch *b_met_px;
extern TBranch *b_met_py;
extern TBranch *b_met_pz_type0;
extern TBranch *b_met_pz_type2;
extern TBranch *b_leptonPt;
extern TBranch *b_leptonEta;
extern TBranch *b_leptonPhi;
extern TBranch *b_leptonE;
extern TBranch *b_AK8jetPt;
extern TBranch *b_AK8jetEta;
extern TBranch *b_AK8jetPhi;
extern TBranch *b_AK8jetE;
extern TBranch *b_AK8jetPrunedMass;
extern TBranch *b_AK8jetTrimmedMass;
extern TBranch *b_AK8jetFilteredMass;
extern TBranch *b_AK8jetTau1;
extern TBranch *b_AK8jetTau2;
extern TBranch *b_AK8jetTau3;
extern TBranch *b_jetPt;
extern TBranch *b_jetEta;
extern TBranch *b_jetPhi;
extern TBranch *b_jetE;
extern TBranch *b_jet_bDiscr;
extern TBranch *b_genBosonPdgId;
extern TBranch *b_genBosonPt;
extern TBranch *b_genBosonEta;
extern TBranch *b_genBosonPhi;
extern TBranch *b_genBosonE;
extern TBranch *b_genLeptonPdgId;
extern TBranch *b_genLeptonPt;
extern TBranch *b_genLeptonEta;
extern TBranch *b_genLeptonPhi;
extern TBranch *b_genLeptonE;
extern TBranch *b_genNuPt;
extern TBranch *b_genNuEta;
extern TBranch *b_genNuPhi;
extern TBranch *b_genNuE;
extern TBranch *b_deltaR_lak8jet;
extern TBranch *b_deltaphi_METak8jet;
extern TBranch *b_deltaphi_Vak8jet;
extern TBranch *b_W_pt;
extern TBranch *b_W_eta;
extern TBranch *b_W_phi;
extern TBranch *b_W_E;
extern TBranch *b_W_mt;

void InitRecoTree(TTree* nt);

void init();

void SetOutTree(TTree* outTree);

#endif
