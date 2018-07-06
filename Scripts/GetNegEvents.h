//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Mon Nov 20 17:43:00 2017 by ROOT version 6.06/01
// from TTree Events/Events
// found on file: /eos/uscms/store/user/lnujj/WpWm_aQGC_Ntuples_Ram/FirstStepOutput/BaconNtuples/WminusToLNuZTo2JJJ_EWK_LO_aQGC_MJJ100PTJ10_TuneCUETP8M1_13TeV-madgraph-pythia8/Output_1.root
//////////////////////////////////////////////////////////

#ifndef GetNegEvents_h
#define GetNegEvents_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.
#include "TClonesArray.h"

class GetNegEvents {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

// Fixed size dimensions of array or collections stored in the TTree if any.
   static const Int_t kMaxGenParticle = 94;
   static const Int_t kMaxLHEWeight = 1164;
   static const Int_t kMaxElectron = 5;
   static const Int_t kMaxMuon = 36;
   static const Int_t kMaxTau = 10;
   static const Int_t kMaxPhoton = 8;
   static const Int_t kMaxPV = 87;
   static const Int_t kMaxAK4CHS = 18;
   static const Int_t kMaxAK8CHS = 7;
   static const Int_t kMaxAddAK8CHS = 7;
   static const Int_t kMaxAK4Puppi = 11;
   static const Int_t kMaxAK8Puppi = 6;
   static const Int_t kMaxAddAK8Puppi = 6;

   // Declaration of leaf types
 //baconhep::TEventInfo *Info;
   UInt_t          runNum;
   UInt_t          evtNum;
   UInt_t          lumiSec;
   UInt_t          metFilterFailBits;
   UInt_t          nPU;
   UInt_t          nPUm;
   UInt_t          nPUp;
   Float_t         nPUmean;
   Float_t         nPUmeanm;
   Float_t         nPUmeanp;
   Float_t         pvx;
   Float_t         pvy;
   Float_t         pvz;
   Float_t         bsx;
   Float_t         bsy;
   Float_t         bsz;
   Float_t         caloMET;
   Float_t         caloMETphi;
   Float_t         pfMET;
   Float_t         pfMETphi;
   Float_t         pfMETCov00;
   Float_t         pfMETCov01;
   Float_t         pfMETCov11;
   Float_t         pfMETC;
   Float_t         pfMETCphi;
   Float_t         pfMETCCov00;
   Float_t         pfMETCCov01;
   Float_t         pfMETCCov11;
   Float_t         pfMETCjerup;
   Float_t         pfMETCjerdn;
   Float_t         pfMETCjenup;
   Float_t         pfMETCjendn;
   Float_t         pfMETCuncup;
   Float_t         pfMETCuncdn;
   Float_t         pfMETCjrsup;
   Float_t         pfMETCjrsdn;
   Float_t         pfMETCphijerup;
   Float_t         pfMETCphijerdn;
   Float_t         pfMETCphijenup;
   Float_t         pfMETCphijendn;
   Float_t         pfMETCphiuncup;
   Float_t         pfMETCphiuncdn;
   Float_t         pfMETCphijrsup;
   Float_t         pfMETCphijrsdn;
   Float_t         puppET;
   Float_t         puppETphi;
   Float_t         puppETCov00;
   Float_t         puppETCov01;
   Float_t         puppETCov11;
   Float_t         puppETC;
   Float_t         puppETCphi;
   Float_t         puppETCCov00;
   Float_t         puppETCCov01;
   Float_t         puppETCCov11;
   Float_t         alpacaMET;
   Float_t         alpacaMETphi;
   Float_t         pcpMET;
   Float_t         pcpMETphi;
   Float_t         trkMET;
   Float_t         trkMETphi;
   Float_t         rhoIso;
   Float_t         rhoJet;
 //bitset<256>     triggerBits;
   Bool_t          hasGoodPV;
 //baconhep::TGenEventInfo *GenEvtInfo;
   Int_t           id_1;
   Int_t           id_2;
   Float_t         x_1;
   Float_t         x_2;
   Float_t         scalePDF;
   Float_t         xs;
   Float_t         weight;
   Int_t           GenParticle_;
   Int_t           GenParticle_parent[kMaxGenParticle];   //[GenParticle_]
   Int_t           GenParticle_pdgId[kMaxGenParticle];   //[GenParticle_]
   Int_t           GenParticle_status[kMaxGenParticle];   //[GenParticle_]
   Float_t         GenParticle_pt[kMaxGenParticle];   //[GenParticle_]
   Float_t         GenParticle_eta[kMaxGenParticle];   //[GenParticle_]
   Float_t         GenParticle_phi[kMaxGenParticle];   //[GenParticle_]
   Float_t         GenParticle_mass[kMaxGenParticle];   //[GenParticle_]
   Float_t         GenParticle_y[kMaxGenParticle];   //[GenParticle_]
   //Int_t           LHEWeight_;
   //Int_t           LHEWeight_id[kMaxLHEWeight];   //[LHEWeight_]
   //Float_t         LHEWeight_weight[kMaxLHEWeight];   //[LHEWeight_]
   Int_t           Electron_;
   Float_t         Electron_pt[kMaxElectron];   //[Electron_]
   Float_t         Electron_eta[kMaxElectron];   //[Electron_]
   Float_t         Electron_phi[kMaxElectron];   //[Electron_]
   Float_t         Electron_scEt[kMaxElectron];   //[Electron_]
   Float_t         Electron_scEta[kMaxElectron];   //[Electron_]
   Float_t         Electron_scPhi[kMaxElectron];   //[Electron_]
   Float_t         Electron_ecalEnergy[kMaxElectron];   //[Electron_]
   Float_t         Electron_pfPt[kMaxElectron];   //[Electron_]
   Float_t         Electron_pfEta[kMaxElectron];   //[Electron_]
   Float_t         Electron_pfPhi[kMaxElectron];   //[Electron_]
   Float_t         Electron_trkIso[kMaxElectron];   //[Electron_]
   Float_t         Electron_ecalIso[kMaxElectron];   //[Electron_]
   Float_t         Electron_hcalIso[kMaxElectron];   //[Electron_]
   Float_t         Electron_hcalDepth1Iso[kMaxElectron];   //[Electron_]
   Float_t         Electron_chHadIso[kMaxElectron];   //[Electron_]
   Float_t         Electron_gammaIso[kMaxElectron];   //[Electron_]
   Float_t         Electron_neuHadIso[kMaxElectron];   //[Electron_]
   Float_t         Electron_puIso[kMaxElectron];   //[Electron_]
   Float_t         Electron_ecalPFClusIso[kMaxElectron];   //[Electron_]
   Float_t         Electron_hcalPFClusIso[kMaxElectron];   //[Electron_]
   Float_t         Electron_puppiChHadIso[kMaxElectron];   //[Electron_]
   Float_t         Electron_puppiGammaIso[kMaxElectron];   //[Electron_]
   Float_t         Electron_puppiNeuHadIso[kMaxElectron];   //[Electron_]
   Float_t         Electron_puppiChHadIsoNoLep[kMaxElectron];   //[Electron_]
   Float_t         Electron_puppiGammaIsoNoLep[kMaxElectron];   //[Electron_]
   Float_t         Electron_puppiNeuHadIsoNoLep[kMaxElectron];   //[Electron_]
   Float_t         Electron_d0[kMaxElectron];   //[Electron_]
   Float_t         Electron_dz[kMaxElectron];   //[Electron_]
   Float_t         Electron_sip3d[kMaxElectron];   //[Electron_]
   Float_t         Electron_sieie[kMaxElectron];   //[Electron_]
   Float_t         Electron_e1x5[kMaxElectron];   //[Electron_]
   Float_t         Electron_e2x5[kMaxElectron];   //[Electron_]
   Float_t         Electron_e5x5[kMaxElectron];   //[Electron_]
   Float_t         Electron_r9[kMaxElectron];   //[Electron_]
   Float_t         Electron_eoverp[kMaxElectron];   //[Electron_]
   Float_t         Electron_hovere[kMaxElectron];   //[Electron_]
   Float_t         Electron_fbrem[kMaxElectron];   //[Electron_]
   Float_t         Electron_dEtaInSeed[kMaxElectron];   //[Electron_]
   Float_t         Electron_dEtaIn[kMaxElectron];   //[Electron_]
   Float_t         Electron_dPhiIn[kMaxElectron];   //[Electron_]
   Float_t         Electron_mva[kMaxElectron];   //[Electron_]
   Int_t           Electron_q[kMaxElectron];   //[Electron_]
   Bool_t          Electron_isConv[kMaxElectron];   //[Electron_]
   UInt_t          Electron_nMissingHits[kMaxElectron];   //[Electron_]
   UInt_t          Electron_typeBits[kMaxElectron];   //[Electron_]
   UInt_t          Electron_fiducialBits[kMaxElectron];   //[Electron_]
   Int_t           Electron_classification[kMaxElectron];   //[Electron_]
   Int_t           Electron_scID[kMaxElectron];   //[Electron_]
   Int_t           Electron_trkID[kMaxElectron];   //[Electron_]
 //bitset<256>     Electron_hltMatchBits[kMaxElectron];
   Int_t           Muon_;
   Float_t         Muon_pt[kMaxMuon];   //[Muon_]
   Float_t         Muon_eta[kMaxMuon];   //[Muon_]
   Float_t         Muon_phi[kMaxMuon];   //[Muon_]
   Float_t         Muon_ptErr[kMaxMuon];   //[Muon_]
   Float_t         Muon_staPt[kMaxMuon];   //[Muon_]
   Float_t         Muon_staEta[kMaxMuon];   //[Muon_]
   Float_t         Muon_staPhi[kMaxMuon];   //[Muon_]
   Float_t         Muon_pfPt[kMaxMuon];   //[Muon_]
   Float_t         Muon_pfEta[kMaxMuon];   //[Muon_]
   Float_t         Muon_pfPhi[kMaxMuon];   //[Muon_]
   Float_t         Muon_trkIso[kMaxMuon];   //[Muon_]
   Float_t         Muon_ecalIso[kMaxMuon];   //[Muon_]
   Float_t         Muon_hcalIso[kMaxMuon];   //[Muon_]
   Float_t         Muon_chHadIso[kMaxMuon];   //[Muon_]
   Float_t         Muon_gammaIso[kMaxMuon];   //[Muon_]
   Float_t         Muon_neuHadIso[kMaxMuon];   //[Muon_]
   Float_t         Muon_puIso[kMaxMuon];   //[Muon_]
   Float_t         Muon_puppiChHadIso[kMaxMuon];   //[Muon_]
   Float_t         Muon_puppiGammaIso[kMaxMuon];   //[Muon_]
   Float_t         Muon_puppiNeuHadIso[kMaxMuon];   //[Muon_]
   Float_t         Muon_puppiChHadIsoNoLep[kMaxMuon];   //[Muon_]
   Float_t         Muon_puppiGammaIsoNoLep[kMaxMuon];   //[Muon_]
   Float_t         Muon_puppiNeuHadIsoNoLep[kMaxMuon];   //[Muon_]
   Float_t         Muon_d0[kMaxMuon];   //[Muon_]
   Float_t         Muon_dz[kMaxMuon];   //[Muon_]
   Float_t         Muon_sip3d[kMaxMuon];   //[Muon_]
   Float_t         Muon_tkNchi2[kMaxMuon];   //[Muon_]
   Float_t         Muon_muNchi2[kMaxMuon];   //[Muon_]
   Float_t         Muon_trkKink[kMaxMuon];   //[Muon_]
   Float_t         Muon_glbKink[kMaxMuon];   //[Muon_]
   Float_t         Muon_trkHitFrac[kMaxMuon];   //[Muon_]
   Float_t         Muon_chi2LocPos[kMaxMuon];   //[Muon_]
   Float_t         Muon_segComp[kMaxMuon];   //[Muon_]
   Float_t         Muon_caloComp[kMaxMuon];   //[Muon_]
   Int_t           Muon_q[kMaxMuon];   //[Muon_]
   Int_t           Muon_nValidHits[kMaxMuon];   //[Muon_]
   UInt_t          Muon_typeBits[kMaxMuon];   //[Muon_]
   UInt_t          Muon_selectorBits[kMaxMuon];   //[Muon_]
   UInt_t          Muon_pogIDBits[kMaxMuon];   //[Muon_]
   UInt_t          Muon_nTkHits[kMaxMuon];   //[Muon_]
   UInt_t          Muon_nPixHits[kMaxMuon];   //[Muon_]
   UInt_t          Muon_nTkLayers[kMaxMuon];   //[Muon_]
   UInt_t          Muon_nPixLayers[kMaxMuon];   //[Muon_]
   UInt_t          Muon_nMatchStn[kMaxMuon];   //[Muon_]
   Int_t           Muon_trkID[kMaxMuon];   //[Muon_]
 //bitset<256>     Muon_hltMatchBits[kMaxMuon];
   Int_t           Tau_;
   Float_t         Tau_pt[kMaxTau];   //[Tau_]
   Float_t         Tau_eta[kMaxTau];   //[Tau_]
   Float_t         Tau_phi[kMaxTau];   //[Tau_]
   Float_t         Tau_m[kMaxTau];   //[Tau_]
   Float_t         Tau_e[kMaxTau];   //[Tau_]
   Int_t           Tau_q[kMaxTau];   //[Tau_]
   Float_t         Tau_dzLeadChHad[kMaxTau];   //[Tau_]
   Float_t         Tau_d0LeadChHad[kMaxTau];   //[Tau_]
   UInt_t          Tau_nSignalChHad[kMaxTau];   //[Tau_]
   UInt_t          Tau_nSignalGamma[kMaxTau];   //[Tau_]
   Int_t           Tau_decaymode[kMaxTau];   //[Tau_]
   Float_t         Tau_antiEleMVA6[kMaxTau];   //[Tau_]
   Float_t         Tau_antiEleMVA6Cat[kMaxTau];   //[Tau_]
   Float_t         Tau_rawMuonRejection[kMaxTau];   //[Tau_]
   Float_t         Tau_rawIso3Hits[kMaxTau];   //[Tau_]
   Float_t         Tau_rawIsoMVA3oldDMwoLT[kMaxTau];   //[Tau_]
   Float_t         Tau_rawIsoMVA3oldDMwLT[kMaxTau];   //[Tau_]
   Float_t         Tau_rawIsoMVA3newDMwoLT[kMaxTau];   //[Tau_]
   Float_t         Tau_rawIsoMVA3newDMwLT[kMaxTau];   //[Tau_]
   Float_t         Tau_puppiChHadIso[kMaxTau];   //[Tau_]
   Float_t         Tau_puppiGammaIso[kMaxTau];   //[Tau_]
   Float_t         Tau_puppiNeuHadIso[kMaxTau];   //[Tau_]
   Float_t         Tau_puppiChHadIsoNoLep[kMaxTau];   //[Tau_]
   Float_t         Tau_puppiGammaIsoNoLep[kMaxTau];   //[Tau_]
   Float_t         Tau_puppiNeuHadIsoNoLep[kMaxTau];   //[Tau_]
   ULong_t         Tau_hpsDisc[kMaxTau];   //[Tau_]
 //bitset<256>     Tau_hltMatchBits[kMaxTau];
   Int_t           Photon_;
   Float_t         Photon_pt[kMaxPhoton];   //[Photon_]
   Float_t         Photon_eta[kMaxPhoton];   //[Photon_]
   Float_t         Photon_phi[kMaxPhoton];   //[Photon_]
   Float_t         Photon_scEt[kMaxPhoton];   //[Photon_]
   Float_t         Photon_scEta[kMaxPhoton];   //[Photon_]
   Float_t         Photon_scPhi[kMaxPhoton];   //[Photon_]
   Float_t         Photon_trkIso[kMaxPhoton];   //[Photon_]
   Float_t         Photon_ecalIso[kMaxPhoton];   //[Photon_]
   Float_t         Photon_hcalIso[kMaxPhoton];   //[Photon_]
   Float_t         Photon_chHadIso[kMaxPhoton];   //[Photon_]
   Float_t         Photon_gammaIso[kMaxPhoton];   //[Photon_]
   Float_t         Photon_neuHadIso[kMaxPhoton];   //[Photon_]
   Float_t         Photon_mva[kMaxPhoton];   //[Photon_]
   Float_t         Photon_hovere[kMaxPhoton];   //[Photon_]
   Float_t         Photon_sthovere[kMaxPhoton];   //[Photon_]
   Float_t         Photon_sieie[kMaxPhoton];   //[Photon_]
   Float_t         Photon_sipip[kMaxPhoton];   //[Photon_]
   Float_t         Photon_r9[kMaxPhoton];   //[Photon_]
   UInt_t          Photon_fiducialBits[kMaxPhoton];   //[Photon_]
   UInt_t          Photon_typeBits[kMaxPhoton];   //[Photon_]
   Int_t           Photon_scID[kMaxPhoton];   //[Photon_]
   Bool_t          Photon_hasPixelSeed[kMaxPhoton];   //[Photon_]
   Bool_t          Photon_passElectronVeto[kMaxPhoton];   //[Photon_]
   Bool_t          Photon_isConv[kMaxPhoton];   //[Photon_]
 //bitset<256>     Photon_hltMatchBits[kMaxPhoton];
   Int_t           PV_;
   UInt_t          PV_nTracksFit[kMaxPV];   //[PV_]
   Float_t         PV_ndof[kMaxPV];   //[PV_]
   Float_t         PV_chi2[kMaxPV];   //[PV_]
   Float_t         PV_x[kMaxPV];   //[PV_]
   Float_t         PV_y[kMaxPV];   //[PV_]
   Float_t         PV_z[kMaxPV];   //[PV_]
   Int_t           AK4CHS_;
   Float_t         AK4CHS_pt[kMaxAK4CHS];   //[AK4CHS_]
   Float_t         AK4CHS_eta[kMaxAK4CHS];   //[AK4CHS_]
   Float_t         AK4CHS_phi[kMaxAK4CHS];   //[AK4CHS_]
   Float_t         AK4CHS_mass[kMaxAK4CHS];   //[AK4CHS_]
   Float_t         AK4CHS_ptRaw[kMaxAK4CHS];   //[AK4CHS_]
   Float_t         AK4CHS_unc[kMaxAK4CHS];   //[AK4CHS_]
   Float_t         AK4CHS_area[kMaxAK4CHS];   //[AK4CHS_]
   Float_t         AK4CHS_d0[kMaxAK4CHS];   //[AK4CHS_]
   Float_t         AK4CHS_dz[kMaxAK4CHS];   //[AK4CHS_]
   Float_t         AK4CHS_csv[kMaxAK4CHS];   //[AK4CHS_]
   Float_t         AK4CHS_bmva[kMaxAK4CHS];   //[AK4CHS_]
   Float_t         AK4CHS_cvb[kMaxAK4CHS];   //[AK4CHS_]
   Float_t         AK4CHS_cvl[kMaxAK4CHS];   //[AK4CHS_]
   Float_t         AK4CHS_qgid[kMaxAK4CHS];   //[AK4CHS_]
   Float_t         AK4CHS_axis2[kMaxAK4CHS];   //[AK4CHS_]
   Float_t         AK4CHS_ptD[kMaxAK4CHS];   //[AK4CHS_]
   Int_t           AK4CHS_mult[kMaxAK4CHS];   //[AK4CHS_]
   Float_t         AK4CHS_q[kMaxAK4CHS];   //[AK4CHS_]
   Float_t         AK4CHS_mva[kMaxAK4CHS];   //[AK4CHS_]
   Float_t         AK4CHS_beta[kMaxAK4CHS];   //[AK4CHS_]
   Float_t         AK4CHS_betaStar[kMaxAK4CHS];   //[AK4CHS_]
   Float_t         AK4CHS_dR2Mean[kMaxAK4CHS];   //[AK4CHS_]
   Float_t         AK4CHS_pullY[kMaxAK4CHS];   //[AK4CHS_]
   Float_t         AK4CHS_pullPhi[kMaxAK4CHS];   //[AK4CHS_]
   Float_t         AK4CHS_chPullY[kMaxAK4CHS];   //[AK4CHS_]
   Float_t         AK4CHS_chPullPhi[kMaxAK4CHS];   //[AK4CHS_]
   Float_t         AK4CHS_neuPullY[kMaxAK4CHS];   //[AK4CHS_]
   Float_t         AK4CHS_neuPullPhi[kMaxAK4CHS];   //[AK4CHS_]
   Float_t         AK4CHS_chEmFrac[kMaxAK4CHS];   //[AK4CHS_]
   Float_t         AK4CHS_neuEmFrac[kMaxAK4CHS];   //[AK4CHS_]
   Float_t         AK4CHS_chHadFrac[kMaxAK4CHS];   //[AK4CHS_]
   Float_t         AK4CHS_neuHadFrac[kMaxAK4CHS];   //[AK4CHS_]
   Float_t         AK4CHS_muonFrac[kMaxAK4CHS];   //[AK4CHS_]
   Float_t         AK4CHS_genpt[kMaxAK4CHS];   //[AK4CHS_]
   Float_t         AK4CHS_geneta[kMaxAK4CHS];   //[AK4CHS_]
   Float_t         AK4CHS_genphi[kMaxAK4CHS];   //[AK4CHS_]
   Float_t         AK4CHS_genm[kMaxAK4CHS];   //[AK4CHS_]
   Int_t           AK4CHS_partonFlavor[kMaxAK4CHS];   //[AK4CHS_]
   Int_t           AK4CHS_hadronFlavor[kMaxAK4CHS];   //[AK4CHS_]
   UInt_t          AK4CHS_nCharged[kMaxAK4CHS];   //[AK4CHS_]
   UInt_t          AK4CHS_nNeutrals[kMaxAK4CHS];   //[AK4CHS_]
   UInt_t          AK4CHS_nParticles[kMaxAK4CHS];   //[AK4CHS_]
 //bitset<256>     AK4CHS_hltMatchBits[kMaxAK4CHS];
   Int_t           AK8CHS_;
   Float_t         AK8CHS_pt[kMaxAK8CHS];   //[AK8CHS_]
   Float_t         AK8CHS_eta[kMaxAK8CHS];   //[AK8CHS_]
   Float_t         AK8CHS_phi[kMaxAK8CHS];   //[AK8CHS_]
   Float_t         AK8CHS_mass[kMaxAK8CHS];   //[AK8CHS_]
   Float_t         AK8CHS_ptRaw[kMaxAK8CHS];   //[AK8CHS_]
   Float_t         AK8CHS_unc[kMaxAK8CHS];   //[AK8CHS_]
   Float_t         AK8CHS_area[kMaxAK8CHS];   //[AK8CHS_]
   Float_t         AK8CHS_d0[kMaxAK8CHS];   //[AK8CHS_]
   Float_t         AK8CHS_dz[kMaxAK8CHS];   //[AK8CHS_]
   Float_t         AK8CHS_csv[kMaxAK8CHS];   //[AK8CHS_]
   Float_t         AK8CHS_bmva[kMaxAK8CHS];   //[AK8CHS_]
   Float_t         AK8CHS_cvb[kMaxAK8CHS];   //[AK8CHS_]
   Float_t         AK8CHS_cvl[kMaxAK8CHS];   //[AK8CHS_]
   Float_t         AK8CHS_qgid[kMaxAK8CHS];   //[AK8CHS_]
   Float_t         AK8CHS_axis2[kMaxAK8CHS];   //[AK8CHS_]
   Float_t         AK8CHS_ptD[kMaxAK8CHS];   //[AK8CHS_]
   Int_t           AK8CHS_mult[kMaxAK8CHS];   //[AK8CHS_]
   Float_t         AK8CHS_q[kMaxAK8CHS];   //[AK8CHS_]
   Float_t         AK8CHS_mva[kMaxAK8CHS];   //[AK8CHS_]
   Float_t         AK8CHS_beta[kMaxAK8CHS];   //[AK8CHS_]
   Float_t         AK8CHS_betaStar[kMaxAK8CHS];   //[AK8CHS_]
   Float_t         AK8CHS_dR2Mean[kMaxAK8CHS];   //[AK8CHS_]
   Float_t         AK8CHS_pullY[kMaxAK8CHS];   //[AK8CHS_]
   Float_t         AK8CHS_pullPhi[kMaxAK8CHS];   //[AK8CHS_]
   Float_t         AK8CHS_chPullY[kMaxAK8CHS];   //[AK8CHS_]
   Float_t         AK8CHS_chPullPhi[kMaxAK8CHS];   //[AK8CHS_]
   Float_t         AK8CHS_neuPullY[kMaxAK8CHS];   //[AK8CHS_]
   Float_t         AK8CHS_neuPullPhi[kMaxAK8CHS];   //[AK8CHS_]
   Float_t         AK8CHS_chEmFrac[kMaxAK8CHS];   //[AK8CHS_]
   Float_t         AK8CHS_neuEmFrac[kMaxAK8CHS];   //[AK8CHS_]
   Float_t         AK8CHS_chHadFrac[kMaxAK8CHS];   //[AK8CHS_]
   Float_t         AK8CHS_neuHadFrac[kMaxAK8CHS];   //[AK8CHS_]
   Float_t         AK8CHS_muonFrac[kMaxAK8CHS];   //[AK8CHS_]
   Float_t         AK8CHS_genpt[kMaxAK8CHS];   //[AK8CHS_]
   Float_t         AK8CHS_geneta[kMaxAK8CHS];   //[AK8CHS_]
   Float_t         AK8CHS_genphi[kMaxAK8CHS];   //[AK8CHS_]
   Float_t         AK8CHS_genm[kMaxAK8CHS];   //[AK8CHS_]
   Int_t           AK8CHS_partonFlavor[kMaxAK8CHS];   //[AK8CHS_]
   Int_t           AK8CHS_hadronFlavor[kMaxAK8CHS];   //[AK8CHS_]
   UInt_t          AK8CHS_nCharged[kMaxAK8CHS];   //[AK8CHS_]
   UInt_t          AK8CHS_nNeutrals[kMaxAK8CHS];   //[AK8CHS_]
   UInt_t          AK8CHS_nParticles[kMaxAK8CHS];   //[AK8CHS_]
 //bitset<256>     AK8CHS_hltMatchBits[kMaxAK8CHS];
   Int_t           AddAK8CHS_;
   UInt_t          AddAK8CHS_index[kMaxAddAK8CHS];   //[AddAK8CHS_]
   Float_t         AddAK8CHS_mass_prun[kMaxAddAK8CHS];   //[AddAK8CHS_]
   Float_t         AddAK8CHS_mass_trim[kMaxAddAK8CHS];   //[AddAK8CHS_]
   Float_t         AddAK8CHS_mass_sd0[kMaxAddAK8CHS];   //[AddAK8CHS_]
   Float_t         AddAK8CHS_c2_0[kMaxAddAK8CHS];   //[AddAK8CHS_]
   Float_t         AddAK8CHS_c2_0P2[kMaxAddAK8CHS];   //[AddAK8CHS_]
   Float_t         AddAK8CHS_c2_0P5[kMaxAddAK8CHS];   //[AddAK8CHS_]
   Float_t         AddAK8CHS_c2_1P0[kMaxAddAK8CHS];   //[AddAK8CHS_]
   Float_t         AddAK8CHS_c2_2P0[kMaxAddAK8CHS];   //[AddAK8CHS_]
   Float_t         AddAK8CHS_e2_b1[kMaxAddAK8CHS];   //[AddAK8CHS_]
   Float_t         AddAK8CHS_e3_b1[kMaxAddAK8CHS];   //[AddAK8CHS_]
   Float_t         AddAK8CHS_e3_v1_b1[kMaxAddAK8CHS];   //[AddAK8CHS_]
   Float_t         AddAK8CHS_e3_v2_b1[kMaxAddAK8CHS];   //[AddAK8CHS_]
   Float_t         AddAK8CHS_e4_v1_b1[kMaxAddAK8CHS];   //[AddAK8CHS_]
   Float_t         AddAK8CHS_e4_v2_b1[kMaxAddAK8CHS];   //[AddAK8CHS_]
   Float_t         AddAK8CHS_e2_b2[kMaxAddAK8CHS];   //[AddAK8CHS_]
   Float_t         AddAK8CHS_e3_b2[kMaxAddAK8CHS];   //[AddAK8CHS_]
   Float_t         AddAK8CHS_e3_v1_b2[kMaxAddAK8CHS];   //[AddAK8CHS_]
   Float_t         AddAK8CHS_e3_v2_b2[kMaxAddAK8CHS];   //[AddAK8CHS_]
   Float_t         AddAK8CHS_e4_v1_b2[kMaxAddAK8CHS];   //[AddAK8CHS_]
   Float_t         AddAK8CHS_e4_v2_b2[kMaxAddAK8CHS];   //[AddAK8CHS_]
   Float_t         AddAK8CHS_e2_sdb1[kMaxAddAK8CHS];   //[AddAK8CHS_]
   Float_t         AddAK8CHS_e3_sdb1[kMaxAddAK8CHS];   //[AddAK8CHS_]
   Float_t         AddAK8CHS_e3_v1_sdb1[kMaxAddAK8CHS];   //[AddAK8CHS_]
   Float_t         AddAK8CHS_e3_v2_sdb1[kMaxAddAK8CHS];   //[AddAK8CHS_]
   Float_t         AddAK8CHS_e4_v1_sdb1[kMaxAddAK8CHS];   //[AddAK8CHS_]
   Float_t         AddAK8CHS_e4_v2_sdb1[kMaxAddAK8CHS];   //[AddAK8CHS_]
   Float_t         AddAK8CHS_e2_sdb2[kMaxAddAK8CHS];   //[AddAK8CHS_]
   Float_t         AddAK8CHS_e3_sdb2[kMaxAddAK8CHS];   //[AddAK8CHS_]
   Float_t         AddAK8CHS_e3_v1_sdb2[kMaxAddAK8CHS];   //[AddAK8CHS_]
   Float_t         AddAK8CHS_e3_v2_sdb2[kMaxAddAK8CHS];   //[AddAK8CHS_]
   Float_t         AddAK8CHS_e4_v1_sdb2[kMaxAddAK8CHS];   //[AddAK8CHS_]
   Float_t         AddAK8CHS_e4_v2_sdb2[kMaxAddAK8CHS];   //[AddAK8CHS_]
   Float_t         AddAK8CHS_qjet[kMaxAddAK8CHS];   //[AddAK8CHS_]
   Float_t         AddAK8CHS_tau1[kMaxAddAK8CHS];   //[AddAK8CHS_]
   Float_t         AddAK8CHS_tau2[kMaxAddAK8CHS];   //[AddAK8CHS_]
   Float_t         AddAK8CHS_tau3[kMaxAddAK8CHS];   //[AddAK8CHS_]
   Float_t         AddAK8CHS_tau4[kMaxAddAK8CHS];   //[AddAK8CHS_]
   Float_t         AddAK8CHS_doublecsv[kMaxAddAK8CHS];   //[AddAK8CHS_]
   Float_t         AddAK8CHS_Double_sub[kMaxAddAK8CHS];   //[AddAK8CHS_]
   Float_t         AddAK8CHS_sj1_pt[kMaxAddAK8CHS];   //[AddAK8CHS_]
   Float_t         AddAK8CHS_sj1_eta[kMaxAddAK8CHS];   //[AddAK8CHS_]
   Float_t         AddAK8CHS_sj1_phi[kMaxAddAK8CHS];   //[AddAK8CHS_]
   Float_t         AddAK8CHS_sj1_m[kMaxAddAK8CHS];   //[AddAK8CHS_]
   Float_t         AddAK8CHS_sj1_csv[kMaxAddAK8CHS];   //[AddAK8CHS_]
   Float_t         AddAK8CHS_sj1_qgid[kMaxAddAK8CHS];   //[AddAK8CHS_]
   Float_t         AddAK8CHS_sj1_q[kMaxAddAK8CHS];   //[AddAK8CHS_]
   Float_t         AddAK8CHS_sj2_pt[kMaxAddAK8CHS];   //[AddAK8CHS_]
   Float_t         AddAK8CHS_sj2_eta[kMaxAddAK8CHS];   //[AddAK8CHS_]
   Float_t         AddAK8CHS_sj2_phi[kMaxAddAK8CHS];   //[AddAK8CHS_]
   Float_t         AddAK8CHS_sj2_m[kMaxAddAK8CHS];   //[AddAK8CHS_]
   Float_t         AddAK8CHS_sj2_csv[kMaxAddAK8CHS];   //[AddAK8CHS_]
   Float_t         AddAK8CHS_sj2_qgid[kMaxAddAK8CHS];   //[AddAK8CHS_]
   Float_t         AddAK8CHS_sj2_q[kMaxAddAK8CHS];   //[AddAK8CHS_]
   Float_t         AddAK8CHS_sj3_pt[kMaxAddAK8CHS];   //[AddAK8CHS_]
   Float_t         AddAK8CHS_sj3_eta[kMaxAddAK8CHS];   //[AddAK8CHS_]
   Float_t         AddAK8CHS_sj3_phi[kMaxAddAK8CHS];   //[AddAK8CHS_]
   Float_t         AddAK8CHS_sj3_m[kMaxAddAK8CHS];   //[AddAK8CHS_]
   Float_t         AddAK8CHS_sj3_csv[kMaxAddAK8CHS];   //[AddAK8CHS_]
   Float_t         AddAK8CHS_sj3_qgid[kMaxAddAK8CHS];   //[AddAK8CHS_]
   Float_t         AddAK8CHS_sj3_q[kMaxAddAK8CHS];   //[AddAK8CHS_]
   Float_t         AddAK8CHS_sj4_pt[kMaxAddAK8CHS];   //[AddAK8CHS_]
   Float_t         AddAK8CHS_sj4_eta[kMaxAddAK8CHS];   //[AddAK8CHS_]
   Float_t         AddAK8CHS_sj4_phi[kMaxAddAK8CHS];   //[AddAK8CHS_]
   Float_t         AddAK8CHS_sj4_m[kMaxAddAK8CHS];   //[AddAK8CHS_]
   Float_t         AddAK8CHS_sj4_csv[kMaxAddAK8CHS];   //[AddAK8CHS_]
   Float_t         AddAK8CHS_sj4_qgid[kMaxAddAK8CHS];   //[AddAK8CHS_]
   Float_t         AddAK8CHS_sj4_q[kMaxAddAK8CHS];   //[AddAK8CHS_]
   Float_t         AddAK8CHS_pullAngle[kMaxAddAK8CHS];   //[AddAK8CHS_]
   UInt_t          AddAK8CHS_topTagType[kMaxAddAK8CHS];   //[AddAK8CHS_]
   UInt_t          AddAK8CHS_top_n_subjets[kMaxAddAK8CHS];   //[AddAK8CHS_]
   Float_t         AddAK8CHS_top_m_min[kMaxAddAK8CHS];   //[AddAK8CHS_]
   Float_t         AddAK8CHS_top_m_123[kMaxAddAK8CHS];   //[AddAK8CHS_]
   Float_t         AddAK8CHS_top_fRec[kMaxAddAK8CHS];   //[AddAK8CHS_]
   Float_t         AddAK8CHS_topchi2[kMaxAddAK8CHS];   //[AddAK8CHS_]
   Int_t           AK4Puppi_;
   Float_t         AK4Puppi_pt[kMaxAK4Puppi];   //[AK4Puppi_]
   Float_t         AK4Puppi_eta[kMaxAK4Puppi];   //[AK4Puppi_]
   Float_t         AK4Puppi_phi[kMaxAK4Puppi];   //[AK4Puppi_]
   Float_t         AK4Puppi_mass[kMaxAK4Puppi];   //[AK4Puppi_]
   Float_t         AK4Puppi_ptRaw[kMaxAK4Puppi];   //[AK4Puppi_]
   Float_t         AK4Puppi_unc[kMaxAK4Puppi];   //[AK4Puppi_]
   Float_t         AK4Puppi_area[kMaxAK4Puppi];   //[AK4Puppi_]
   Float_t         AK4Puppi_d0[kMaxAK4Puppi];   //[AK4Puppi_]
   Float_t         AK4Puppi_dz[kMaxAK4Puppi];   //[AK4Puppi_]
   Float_t         AK4Puppi_csv[kMaxAK4Puppi];   //[AK4Puppi_]
   Float_t         AK4Puppi_bmva[kMaxAK4Puppi];   //[AK4Puppi_]
   Float_t         AK4Puppi_cvb[kMaxAK4Puppi];   //[AK4Puppi_]
   Float_t         AK4Puppi_cvl[kMaxAK4Puppi];   //[AK4Puppi_]
   Float_t         AK4Puppi_qgid[kMaxAK4Puppi];   //[AK4Puppi_]
   Float_t         AK4Puppi_axis2[kMaxAK4Puppi];   //[AK4Puppi_]
   Float_t         AK4Puppi_ptD[kMaxAK4Puppi];   //[AK4Puppi_]
   Int_t           AK4Puppi_mult[kMaxAK4Puppi];   //[AK4Puppi_]
   Float_t         AK4Puppi_q[kMaxAK4Puppi];   //[AK4Puppi_]
   Float_t         AK4Puppi_mva[kMaxAK4Puppi];   //[AK4Puppi_]
   Float_t         AK4Puppi_beta[kMaxAK4Puppi];   //[AK4Puppi_]
   Float_t         AK4Puppi_betaStar[kMaxAK4Puppi];   //[AK4Puppi_]
   Float_t         AK4Puppi_dR2Mean[kMaxAK4Puppi];   //[AK4Puppi_]
   Float_t         AK4Puppi_pullY[kMaxAK4Puppi];   //[AK4Puppi_]
   Float_t         AK4Puppi_pullPhi[kMaxAK4Puppi];   //[AK4Puppi_]
   Float_t         AK4Puppi_chPullY[kMaxAK4Puppi];   //[AK4Puppi_]
   Float_t         AK4Puppi_chPullPhi[kMaxAK4Puppi];   //[AK4Puppi_]
   Float_t         AK4Puppi_neuPullY[kMaxAK4Puppi];   //[AK4Puppi_]
   Float_t         AK4Puppi_neuPullPhi[kMaxAK4Puppi];   //[AK4Puppi_]
   Float_t         AK4Puppi_chEmFrac[kMaxAK4Puppi];   //[AK4Puppi_]
   Float_t         AK4Puppi_neuEmFrac[kMaxAK4Puppi];   //[AK4Puppi_]
   Float_t         AK4Puppi_chHadFrac[kMaxAK4Puppi];   //[AK4Puppi_]
   Float_t         AK4Puppi_neuHadFrac[kMaxAK4Puppi];   //[AK4Puppi_]
   Float_t         AK4Puppi_muonFrac[kMaxAK4Puppi];   //[AK4Puppi_]
   Float_t         AK4Puppi_genpt[kMaxAK4Puppi];   //[AK4Puppi_]
   Float_t         AK4Puppi_geneta[kMaxAK4Puppi];   //[AK4Puppi_]
   Float_t         AK4Puppi_genphi[kMaxAK4Puppi];   //[AK4Puppi_]
   Float_t         AK4Puppi_genm[kMaxAK4Puppi];   //[AK4Puppi_]
   Int_t           AK4Puppi_partonFlavor[kMaxAK4Puppi];   //[AK4Puppi_]
   Int_t           AK4Puppi_hadronFlavor[kMaxAK4Puppi];   //[AK4Puppi_]
   UInt_t          AK4Puppi_nCharged[kMaxAK4Puppi];   //[AK4Puppi_]
   UInt_t          AK4Puppi_nNeutrals[kMaxAK4Puppi];   //[AK4Puppi_]
   UInt_t          AK4Puppi_nParticles[kMaxAK4Puppi];   //[AK4Puppi_]
 //bitset<256>     AK4Puppi_hltMatchBits[kMaxAK4Puppi];
   Int_t           AK8Puppi_;
   Float_t         AK8Puppi_pt[kMaxAK8Puppi];   //[AK8Puppi_]
   Float_t         AK8Puppi_eta[kMaxAK8Puppi];   //[AK8Puppi_]
   Float_t         AK8Puppi_phi[kMaxAK8Puppi];   //[AK8Puppi_]
   Float_t         AK8Puppi_mass[kMaxAK8Puppi];   //[AK8Puppi_]
   Float_t         AK8Puppi_ptRaw[kMaxAK8Puppi];   //[AK8Puppi_]
   Float_t         AK8Puppi_unc[kMaxAK8Puppi];   //[AK8Puppi_]
   Float_t         AK8Puppi_area[kMaxAK8Puppi];   //[AK8Puppi_]
   Float_t         AK8Puppi_d0[kMaxAK8Puppi];   //[AK8Puppi_]
   Float_t         AK8Puppi_dz[kMaxAK8Puppi];   //[AK8Puppi_]
   Float_t         AK8Puppi_csv[kMaxAK8Puppi];   //[AK8Puppi_]
   Float_t         AK8Puppi_bmva[kMaxAK8Puppi];   //[AK8Puppi_]
   Float_t         AK8Puppi_cvb[kMaxAK8Puppi];   //[AK8Puppi_]
   Float_t         AK8Puppi_cvl[kMaxAK8Puppi];   //[AK8Puppi_]
   Float_t         AK8Puppi_qgid[kMaxAK8Puppi];   //[AK8Puppi_]
   Float_t         AK8Puppi_axis2[kMaxAK8Puppi];   //[AK8Puppi_]
   Float_t         AK8Puppi_ptD[kMaxAK8Puppi];   //[AK8Puppi_]
   Int_t           AK8Puppi_mult[kMaxAK8Puppi];   //[AK8Puppi_]
   Float_t         AK8Puppi_q[kMaxAK8Puppi];   //[AK8Puppi_]
   Float_t         AK8Puppi_mva[kMaxAK8Puppi];   //[AK8Puppi_]
   Float_t         AK8Puppi_beta[kMaxAK8Puppi];   //[AK8Puppi_]
   Float_t         AK8Puppi_betaStar[kMaxAK8Puppi];   //[AK8Puppi_]
   Float_t         AK8Puppi_dR2Mean[kMaxAK8Puppi];   //[AK8Puppi_]
   Float_t         AK8Puppi_pullY[kMaxAK8Puppi];   //[AK8Puppi_]
   Float_t         AK8Puppi_pullPhi[kMaxAK8Puppi];   //[AK8Puppi_]
   Float_t         AK8Puppi_chPullY[kMaxAK8Puppi];   //[AK8Puppi_]
   Float_t         AK8Puppi_chPullPhi[kMaxAK8Puppi];   //[AK8Puppi_]
   Float_t         AK8Puppi_neuPullY[kMaxAK8Puppi];   //[AK8Puppi_]
   Float_t         AK8Puppi_neuPullPhi[kMaxAK8Puppi];   //[AK8Puppi_]
   Float_t         AK8Puppi_chEmFrac[kMaxAK8Puppi];   //[AK8Puppi_]
   Float_t         AK8Puppi_neuEmFrac[kMaxAK8Puppi];   //[AK8Puppi_]
   Float_t         AK8Puppi_chHadFrac[kMaxAK8Puppi];   //[AK8Puppi_]
   Float_t         AK8Puppi_neuHadFrac[kMaxAK8Puppi];   //[AK8Puppi_]
   Float_t         AK8Puppi_muonFrac[kMaxAK8Puppi];   //[AK8Puppi_]
   Float_t         AK8Puppi_genpt[kMaxAK8Puppi];   //[AK8Puppi_]
   Float_t         AK8Puppi_geneta[kMaxAK8Puppi];   //[AK8Puppi_]
   Float_t         AK8Puppi_genphi[kMaxAK8Puppi];   //[AK8Puppi_]
   Float_t         AK8Puppi_genm[kMaxAK8Puppi];   //[AK8Puppi_]
   Int_t           AK8Puppi_partonFlavor[kMaxAK8Puppi];   //[AK8Puppi_]
   Int_t           AK8Puppi_hadronFlavor[kMaxAK8Puppi];   //[AK8Puppi_]
   UInt_t          AK8Puppi_nCharged[kMaxAK8Puppi];   //[AK8Puppi_]
   UInt_t          AK8Puppi_nNeutrals[kMaxAK8Puppi];   //[AK8Puppi_]
   UInt_t          AK8Puppi_nParticles[kMaxAK8Puppi];   //[AK8Puppi_]
 //bitset<256>     AK8Puppi_hltMatchBits[kMaxAK8Puppi];
   Int_t           AddAK8Puppi_;
   UInt_t          AddAK8Puppi_index[kMaxAddAK8Puppi];   //[AddAK8Puppi_]
   Float_t         AddAK8Puppi_mass_prun[kMaxAddAK8Puppi];   //[AddAK8Puppi_]
   Float_t         AddAK8Puppi_mass_trim[kMaxAddAK8Puppi];   //[AddAK8Puppi_]
   Float_t         AddAK8Puppi_mass_sd0[kMaxAddAK8Puppi];   //[AddAK8Puppi_]
   Float_t         AddAK8Puppi_c2_0[kMaxAddAK8Puppi];   //[AddAK8Puppi_]
   Float_t         AddAK8Puppi_c2_0P2[kMaxAddAK8Puppi];   //[AddAK8Puppi_]
   Float_t         AddAK8Puppi_c2_0P5[kMaxAddAK8Puppi];   //[AddAK8Puppi_]
   Float_t         AddAK8Puppi_c2_1P0[kMaxAddAK8Puppi];   //[AddAK8Puppi_]
   Float_t         AddAK8Puppi_c2_2P0[kMaxAddAK8Puppi];   //[AddAK8Puppi_]
   Float_t         AddAK8Puppi_e2_b1[kMaxAddAK8Puppi];   //[AddAK8Puppi_]
   Float_t         AddAK8Puppi_e3_b1[kMaxAddAK8Puppi];   //[AddAK8Puppi_]
   Float_t         AddAK8Puppi_e3_v1_b1[kMaxAddAK8Puppi];   //[AddAK8Puppi_]
   Float_t         AddAK8Puppi_e3_v2_b1[kMaxAddAK8Puppi];   //[AddAK8Puppi_]
   Float_t         AddAK8Puppi_e4_v1_b1[kMaxAddAK8Puppi];   //[AddAK8Puppi_]
   Float_t         AddAK8Puppi_e4_v2_b1[kMaxAddAK8Puppi];   //[AddAK8Puppi_]
   Float_t         AddAK8Puppi_e2_b2[kMaxAddAK8Puppi];   //[AddAK8Puppi_]
   Float_t         AddAK8Puppi_e3_b2[kMaxAddAK8Puppi];   //[AddAK8Puppi_]
   Float_t         AddAK8Puppi_e3_v1_b2[kMaxAddAK8Puppi];   //[AddAK8Puppi_]
   Float_t         AddAK8Puppi_e3_v2_b2[kMaxAddAK8Puppi];   //[AddAK8Puppi_]
   Float_t         AddAK8Puppi_e4_v1_b2[kMaxAddAK8Puppi];   //[AddAK8Puppi_]
   Float_t         AddAK8Puppi_e4_v2_b2[kMaxAddAK8Puppi];   //[AddAK8Puppi_]
   Float_t         AddAK8Puppi_e2_sdb1[kMaxAddAK8Puppi];   //[AddAK8Puppi_]
   Float_t         AddAK8Puppi_e3_sdb1[kMaxAddAK8Puppi];   //[AddAK8Puppi_]
   Float_t         AddAK8Puppi_e3_v1_sdb1[kMaxAddAK8Puppi];   //[AddAK8Puppi_]
   Float_t         AddAK8Puppi_e3_v2_sdb1[kMaxAddAK8Puppi];   //[AddAK8Puppi_]
   Float_t         AddAK8Puppi_e4_v1_sdb1[kMaxAddAK8Puppi];   //[AddAK8Puppi_]
   Float_t         AddAK8Puppi_e4_v2_sdb1[kMaxAddAK8Puppi];   //[AddAK8Puppi_]
   Float_t         AddAK8Puppi_e2_sdb2[kMaxAddAK8Puppi];   //[AddAK8Puppi_]
   Float_t         AddAK8Puppi_e3_sdb2[kMaxAddAK8Puppi];   //[AddAK8Puppi_]
   Float_t         AddAK8Puppi_e3_v1_sdb2[kMaxAddAK8Puppi];   //[AddAK8Puppi_]
   Float_t         AddAK8Puppi_e3_v2_sdb2[kMaxAddAK8Puppi];   //[AddAK8Puppi_]
   Float_t         AddAK8Puppi_e4_v1_sdb2[kMaxAddAK8Puppi];   //[AddAK8Puppi_]
   Float_t         AddAK8Puppi_e4_v2_sdb2[kMaxAddAK8Puppi];   //[AddAK8Puppi_]
   Float_t         AddAK8Puppi_qjet[kMaxAddAK8Puppi];   //[AddAK8Puppi_]
   Float_t         AddAK8Puppi_tau1[kMaxAddAK8Puppi];   //[AddAK8Puppi_]
   Float_t         AddAK8Puppi_tau2[kMaxAddAK8Puppi];   //[AddAK8Puppi_]
   Float_t         AddAK8Puppi_tau3[kMaxAddAK8Puppi];   //[AddAK8Puppi_]
   Float_t         AddAK8Puppi_tau4[kMaxAddAK8Puppi];   //[AddAK8Puppi_]
   Float_t         AddAK8Puppi_doublecsv[kMaxAddAK8Puppi];   //[AddAK8Puppi_]
   Float_t         AddAK8Puppi_Double_sub[kMaxAddAK8Puppi];   //[AddAK8Puppi_]
   Float_t         AddAK8Puppi_sj1_pt[kMaxAddAK8Puppi];   //[AddAK8Puppi_]
   Float_t         AddAK8Puppi_sj1_eta[kMaxAddAK8Puppi];   //[AddAK8Puppi_]
   Float_t         AddAK8Puppi_sj1_phi[kMaxAddAK8Puppi];   //[AddAK8Puppi_]
   Float_t         AddAK8Puppi_sj1_m[kMaxAddAK8Puppi];   //[AddAK8Puppi_]
   Float_t         AddAK8Puppi_sj1_csv[kMaxAddAK8Puppi];   //[AddAK8Puppi_]
   Float_t         AddAK8Puppi_sj1_qgid[kMaxAddAK8Puppi];   //[AddAK8Puppi_]
   Float_t         AddAK8Puppi_sj1_q[kMaxAddAK8Puppi];   //[AddAK8Puppi_]
   Float_t         AddAK8Puppi_sj2_pt[kMaxAddAK8Puppi];   //[AddAK8Puppi_]
   Float_t         AddAK8Puppi_sj2_eta[kMaxAddAK8Puppi];   //[AddAK8Puppi_]
   Float_t         AddAK8Puppi_sj2_phi[kMaxAddAK8Puppi];   //[AddAK8Puppi_]
   Float_t         AddAK8Puppi_sj2_m[kMaxAddAK8Puppi];   //[AddAK8Puppi_]
   Float_t         AddAK8Puppi_sj2_csv[kMaxAddAK8Puppi];   //[AddAK8Puppi_]
   Float_t         AddAK8Puppi_sj2_qgid[kMaxAddAK8Puppi];   //[AddAK8Puppi_]
   Float_t         AddAK8Puppi_sj2_q[kMaxAddAK8Puppi];   //[AddAK8Puppi_]
   Float_t         AddAK8Puppi_sj3_pt[kMaxAddAK8Puppi];   //[AddAK8Puppi_]
   Float_t         AddAK8Puppi_sj3_eta[kMaxAddAK8Puppi];   //[AddAK8Puppi_]
   Float_t         AddAK8Puppi_sj3_phi[kMaxAddAK8Puppi];   //[AddAK8Puppi_]
   Float_t         AddAK8Puppi_sj3_m[kMaxAddAK8Puppi];   //[AddAK8Puppi_]
   Float_t         AddAK8Puppi_sj3_csv[kMaxAddAK8Puppi];   //[AddAK8Puppi_]
   Float_t         AddAK8Puppi_sj3_qgid[kMaxAddAK8Puppi];   //[AddAK8Puppi_]
   Float_t         AddAK8Puppi_sj3_q[kMaxAddAK8Puppi];   //[AddAK8Puppi_]
   Float_t         AddAK8Puppi_sj4_pt[kMaxAddAK8Puppi];   //[AddAK8Puppi_]
   Float_t         AddAK8Puppi_sj4_eta[kMaxAddAK8Puppi];   //[AddAK8Puppi_]
   Float_t         AddAK8Puppi_sj4_phi[kMaxAddAK8Puppi];   //[AddAK8Puppi_]
   Float_t         AddAK8Puppi_sj4_m[kMaxAddAK8Puppi];   //[AddAK8Puppi_]
   Float_t         AddAK8Puppi_sj4_csv[kMaxAddAK8Puppi];   //[AddAK8Puppi_]
   Float_t         AddAK8Puppi_sj4_qgid[kMaxAddAK8Puppi];   //[AddAK8Puppi_]
   Float_t         AddAK8Puppi_sj4_q[kMaxAddAK8Puppi];   //[AddAK8Puppi_]
   Float_t         AddAK8Puppi_pullAngle[kMaxAddAK8Puppi];   //[AddAK8Puppi_]
   UInt_t          AddAK8Puppi_topTagType[kMaxAddAK8Puppi];   //[AddAK8Puppi_]
   UInt_t          AddAK8Puppi_top_n_subjets[kMaxAddAK8Puppi];   //[AddAK8Puppi_]
   Float_t         AddAK8Puppi_top_m_min[kMaxAddAK8Puppi];   //[AddAK8Puppi_]
   Float_t         AddAK8Puppi_top_m_123[kMaxAddAK8Puppi];   //[AddAK8Puppi_]
   Float_t         AddAK8Puppi_top_fRec[kMaxAddAK8Puppi];   //[AddAK8Puppi_]
   Float_t         AddAK8Puppi_topchi2[kMaxAddAK8Puppi];   //[AddAK8Puppi_]

   // List of branches
   TBranch        *b_Info_runNum;   //!
   TBranch        *b_Info_evtNum;   //!
   TBranch        *b_Info_lumiSec;   //!
   TBranch        *b_Info_metFilterFailBits;   //!
   TBranch        *b_Info_nPU;   //!
   TBranch        *b_Info_nPUm;   //!
   TBranch        *b_Info_nPUp;   //!
   TBranch        *b_Info_nPUmean;   //!
   TBranch        *b_Info_nPUmeanm;   //!
   TBranch        *b_Info_nPUmeanp;   //!
   TBranch        *b_Info_pvx;   //!
   TBranch        *b_Info_pvy;   //!
   TBranch        *b_Info_pvz;   //!
   TBranch        *b_Info_bsx;   //!
   TBranch        *b_Info_bsy;   //!
   TBranch        *b_Info_bsz;   //!
   TBranch        *b_Info_caloMET;   //!
   TBranch        *b_Info_caloMETphi;   //!
   TBranch        *b_Info_pfMET;   //!
   TBranch        *b_Info_pfMETphi;   //!
   TBranch        *b_Info_pfMETCov00;   //!
   TBranch        *b_Info_pfMETCov01;   //!
   TBranch        *b_Info_pfMETCov11;   //!
   TBranch        *b_Info_pfMETC;   //!
   TBranch        *b_Info_pfMETCphi;   //!
   TBranch        *b_Info_pfMETCCov00;   //!
   TBranch        *b_Info_pfMETCCov01;   //!
   TBranch        *b_Info_pfMETCCov11;   //!
   TBranch        *b_Info_pfMETCjerup;   //!
   TBranch        *b_Info_pfMETCjerdn;   //!
   TBranch        *b_Info_pfMETCjenup;   //!
   TBranch        *b_Info_pfMETCjendn;   //!
   TBranch        *b_Info_pfMETCuncup;   //!
   TBranch        *b_Info_pfMETCuncdn;   //!
   TBranch        *b_Info_pfMETCjrsup;   //!
   TBranch        *b_Info_pfMETCjrsdn;   //!
   TBranch        *b_Info_pfMETCphijerup;   //!
   TBranch        *b_Info_pfMETCphijerdn;   //!
   TBranch        *b_Info_pfMETCphijenup;   //!
   TBranch        *b_Info_pfMETCphijendn;   //!
   TBranch        *b_Info_pfMETCphiuncup;   //!
   TBranch        *b_Info_pfMETCphiuncdn;   //!
   TBranch        *b_Info_pfMETCphijrsup;   //!
   TBranch        *b_Info_pfMETCphijrsdn;   //!
   TBranch        *b_Info_puppET;   //!
   TBranch        *b_Info_puppETphi;   //!
   TBranch        *b_Info_puppETCov00;   //!
   TBranch        *b_Info_puppETCov01;   //!
   TBranch        *b_Info_puppETCov11;   //!
   TBranch        *b_Info_puppETC;   //!
   TBranch        *b_Info_puppETCphi;   //!
   TBranch        *b_Info_puppETCCov00;   //!
   TBranch        *b_Info_puppETCCov01;   //!
   TBranch        *b_Info_puppETCCov11;   //!
   TBranch        *b_Info_alpacaMET;   //!
   TBranch        *b_Info_alpacaMETphi;   //!
   TBranch        *b_Info_pcpMET;   //!
   TBranch        *b_Info_pcpMETphi;   //!
   TBranch        *b_Info_trkMET;   //!
   TBranch        *b_Info_trkMETphi;   //!
   TBranch        *b_Info_rhoIso;   //!
   TBranch        *b_Info_rhoJet;   //!
   TBranch        *b_Info_hasGoodPV;   //!
   TBranch        *b_GenEvtInfo_id_1;   //!
   TBranch        *b_GenEvtInfo_id_2;   //!
   TBranch        *b_GenEvtInfo_x_1;   //!
   TBranch        *b_GenEvtInfo_x_2;   //!
   TBranch        *b_GenEvtInfo_scalePDF;   //!
   TBranch        *b_GenEvtInfo_xs;   //!
   TBranch        *b_GenEvtInfo_weight;   //!
   TBranch        *b_GenParticle_;   //!
   TBranch        *b_GenParticle_parent;   //!
   TBranch        *b_GenParticle_pdgId;   //!
   TBranch        *b_GenParticle_status;   //!
   TBranch        *b_GenParticle_pt;   //!
   TBranch        *b_GenParticle_eta;   //!
   TBranch        *b_GenParticle_phi;   //!
   TBranch        *b_GenParticle_mass;   //!
   TBranch        *b_GenParticle_y;   //!
   //TBranch        *b_LHEWeight_;   //!
   //TBranch        *b_LHEWeight_id;   //!
   //TBranch        *b_LHEWeight_weight;   //!
   TBranch        *b_Electron_;   //!
   TBranch        *b_Electron_pt;   //!
   TBranch        *b_Electron_eta;   //!
   TBranch        *b_Electron_phi;   //!
   TBranch        *b_Electron_scEt;   //!
   TBranch        *b_Electron_scEta;   //!
   TBranch        *b_Electron_scPhi;   //!
   TBranch        *b_Electron_ecalEnergy;   //!
   TBranch        *b_Electron_pfPt;   //!
   TBranch        *b_Electron_pfEta;   //!
   TBranch        *b_Electron_pfPhi;   //!
   TBranch        *b_Electron_trkIso;   //!
   TBranch        *b_Electron_ecalIso;   //!
   TBranch        *b_Electron_hcalIso;   //!
   TBranch        *b_Electron_hcalDepth1Iso;   //!
   TBranch        *b_Electron_chHadIso;   //!
   TBranch        *b_Electron_gammaIso;   //!
   TBranch        *b_Electron_neuHadIso;   //!
   TBranch        *b_Electron_puIso;   //!
   TBranch        *b_Electron_ecalPFClusIso;   //!
   TBranch        *b_Electron_hcalPFClusIso;   //!
   TBranch        *b_Electron_puppiChHadIso;   //!
   TBranch        *b_Electron_puppiGammaIso;   //!
   TBranch        *b_Electron_puppiNeuHadIso;   //!
   TBranch        *b_Electron_puppiChHadIsoNoLep;   //!
   TBranch        *b_Electron_puppiGammaIsoNoLep;   //!
   TBranch        *b_Electron_puppiNeuHadIsoNoLep;   //!
   TBranch        *b_Electron_d0;   //!
   TBranch        *b_Electron_dz;   //!
   TBranch        *b_Electron_sip3d;   //!
   TBranch        *b_Electron_sieie;   //!
   TBranch        *b_Electron_e1x5;   //!
   TBranch        *b_Electron_e2x5;   //!
   TBranch        *b_Electron_e5x5;   //!
   TBranch        *b_Electron_r9;   //!
   TBranch        *b_Electron_eoverp;   //!
   TBranch        *b_Electron_hovere;   //!
   TBranch        *b_Electron_fbrem;   //!
   TBranch        *b_Electron_dEtaInSeed;   //!
   TBranch        *b_Electron_dEtaIn;   //!
   TBranch        *b_Electron_dPhiIn;   //!
   TBranch        *b_Electron_mva;   //!
   TBranch        *b_Electron_q;   //!
   TBranch        *b_Electron_isConv;   //!
   TBranch        *b_Electron_nMissingHits;   //!
   TBranch        *b_Electron_typeBits;   //!
   TBranch        *b_Electron_fiducialBits;   //!
   TBranch        *b_Electron_classification;   //!
   TBranch        *b_Electron_scID;   //!
   TBranch        *b_Electron_trkID;   //!
   TBranch        *b_Muon_;   //!
   TBranch        *b_Muon_pt;   //!
   TBranch        *b_Muon_eta;   //!
   TBranch        *b_Muon_phi;   //!
   TBranch        *b_Muon_ptErr;   //!
   TBranch        *b_Muon_staPt;   //!
   TBranch        *b_Muon_staEta;   //!
   TBranch        *b_Muon_staPhi;   //!
   TBranch        *b_Muon_pfPt;   //!
   TBranch        *b_Muon_pfEta;   //!
   TBranch        *b_Muon_pfPhi;   //!
   TBranch        *b_Muon_trkIso;   //!
   TBranch        *b_Muon_ecalIso;   //!
   TBranch        *b_Muon_hcalIso;   //!
   TBranch        *b_Muon_chHadIso;   //!
   TBranch        *b_Muon_gammaIso;   //!
   TBranch        *b_Muon_neuHadIso;   //!
   TBranch        *b_Muon_puIso;   //!
   TBranch        *b_Muon_puppiChHadIso;   //!
   TBranch        *b_Muon_puppiGammaIso;   //!
   TBranch        *b_Muon_puppiNeuHadIso;   //!
   TBranch        *b_Muon_puppiChHadIsoNoLep;   //!
   TBranch        *b_Muon_puppiGammaIsoNoLep;   //!
   TBranch        *b_Muon_puppiNeuHadIsoNoLep;   //!
   TBranch        *b_Muon_d0;   //!
   TBranch        *b_Muon_dz;   //!
   TBranch        *b_Muon_sip3d;   //!
   TBranch        *b_Muon_tkNchi2;   //!
   TBranch        *b_Muon_muNchi2;   //!
   TBranch        *b_Muon_trkKink;   //!
   TBranch        *b_Muon_glbKink;   //!
   TBranch        *b_Muon_trkHitFrac;   //!
   TBranch        *b_Muon_chi2LocPos;   //!
   TBranch        *b_Muon_segComp;   //!
   TBranch        *b_Muon_caloComp;   //!
   TBranch        *b_Muon_q;   //!
   TBranch        *b_Muon_nValidHits;   //!
   TBranch        *b_Muon_typeBits;   //!
   TBranch        *b_Muon_selectorBits;   //!
   TBranch        *b_Muon_pogIDBits;   //!
   TBranch        *b_Muon_nTkHits;   //!
   TBranch        *b_Muon_nPixHits;   //!
   TBranch        *b_Muon_nTkLayers;   //!
   TBranch        *b_Muon_nPixLayers;   //!
   TBranch        *b_Muon_nMatchStn;   //!
   TBranch        *b_Muon_trkID;   //!
   TBranch        *b_Tau_;   //!
   TBranch        *b_Tau_pt;   //!
   TBranch        *b_Tau_eta;   //!
   TBranch        *b_Tau_phi;   //!
   TBranch        *b_Tau_m;   //!
   TBranch        *b_Tau_e;   //!
   TBranch        *b_Tau_q;   //!
   TBranch        *b_Tau_dzLeadChHad;   //!
   TBranch        *b_Tau_d0LeadChHad;   //!
   TBranch        *b_Tau_nSignalChHad;   //!
   TBranch        *b_Tau_nSignalGamma;   //!
   TBranch        *b_Tau_decaymode;   //!
   TBranch        *b_Tau_antiEleMVA6;   //!
   TBranch        *b_Tau_antiEleMVA6Cat;   //!
   TBranch        *b_Tau_rawMuonRejection;   //!
   TBranch        *b_Tau_rawIso3Hits;   //!
   TBranch        *b_Tau_rawIsoMVA3oldDMwoLT;   //!
   TBranch        *b_Tau_rawIsoMVA3oldDMwLT;   //!
   TBranch        *b_Tau_rawIsoMVA3newDMwoLT;   //!
   TBranch        *b_Tau_rawIsoMVA3newDMwLT;   //!
   TBranch        *b_Tau_puppiChHadIso;   //!
   TBranch        *b_Tau_puppiGammaIso;   //!
   TBranch        *b_Tau_puppiNeuHadIso;   //!
   TBranch        *b_Tau_puppiChHadIsoNoLep;   //!
   TBranch        *b_Tau_puppiGammaIsoNoLep;   //!
   TBranch        *b_Tau_puppiNeuHadIsoNoLep;   //!
   TBranch        *b_Tau_hpsDisc;   //!
   TBranch        *b_Photon_;   //!
   TBranch        *b_Photon_pt;   //!
   TBranch        *b_Photon_eta;   //!
   TBranch        *b_Photon_phi;   //!
   TBranch        *b_Photon_scEt;   //!
   TBranch        *b_Photon_scEta;   //!
   TBranch        *b_Photon_scPhi;   //!
   TBranch        *b_Photon_trkIso;   //!
   TBranch        *b_Photon_ecalIso;   //!
   TBranch        *b_Photon_hcalIso;   //!
   TBranch        *b_Photon_chHadIso;   //!
   TBranch        *b_Photon_gammaIso;   //!
   TBranch        *b_Photon_neuHadIso;   //!
   TBranch        *b_Photon_mva;   //!
   TBranch        *b_Photon_hovere;   //!
   TBranch        *b_Photon_sthovere;   //!
   TBranch        *b_Photon_sieie;   //!
   TBranch        *b_Photon_sipip;   //!
   TBranch        *b_Photon_r9;   //!
   TBranch        *b_Photon_fiducialBits;   //!
   TBranch        *b_Photon_typeBits;   //!
   TBranch        *b_Photon_scID;   //!
   TBranch        *b_Photon_hasPixelSeed;   //!
   TBranch        *b_Photon_passElectronVeto;   //!
   TBranch        *b_Photon_isConv;   //!
   TBranch        *b_PV_;   //!
   TBranch        *b_PV_nTracksFit;   //!
   TBranch        *b_PV_ndof;   //!
   TBranch        *b_PV_chi2;   //!
   TBranch        *b_PV_x;   //!
   TBranch        *b_PV_y;   //!
   TBranch        *b_PV_z;   //!
   TBranch        *b_AK4CHS_;   //!
   TBranch        *b_AK4CHS_pt;   //!
   TBranch        *b_AK4CHS_eta;   //!
   TBranch        *b_AK4CHS_phi;   //!
   TBranch        *b_AK4CHS_mass;   //!
   TBranch        *b_AK4CHS_ptRaw;   //!
   TBranch        *b_AK4CHS_unc;   //!
   TBranch        *b_AK4CHS_area;   //!
   TBranch        *b_AK4CHS_d0;   //!
   TBranch        *b_AK4CHS_dz;   //!
   TBranch        *b_AK4CHS_csv;   //!
   TBranch        *b_AK4CHS_bmva;   //!
   TBranch        *b_AK4CHS_cvb;   //!
   TBranch        *b_AK4CHS_cvl;   //!
   TBranch        *b_AK4CHS_qgid;   //!
   TBranch        *b_AK4CHS_axis2;   //!
   TBranch        *b_AK4CHS_ptD;   //!
   TBranch        *b_AK4CHS_mult;   //!
   TBranch        *b_AK4CHS_q;   //!
   TBranch        *b_AK4CHS_mva;   //!
   TBranch        *b_AK4CHS_beta;   //!
   TBranch        *b_AK4CHS_betaStar;   //!
   TBranch        *b_AK4CHS_dR2Mean;   //!
   TBranch        *b_AK4CHS_pullY;   //!
   TBranch        *b_AK4CHS_pullPhi;   //!
   TBranch        *b_AK4CHS_chPullY;   //!
   TBranch        *b_AK4CHS_chPullPhi;   //!
   TBranch        *b_AK4CHS_neuPullY;   //!
   TBranch        *b_AK4CHS_neuPullPhi;   //!
   TBranch        *b_AK4CHS_chEmFrac;   //!
   TBranch        *b_AK4CHS_neuEmFrac;   //!
   TBranch        *b_AK4CHS_chHadFrac;   //!
   TBranch        *b_AK4CHS_neuHadFrac;   //!
   TBranch        *b_AK4CHS_muonFrac;   //!
   TBranch        *b_AK4CHS_genpt;   //!
   TBranch        *b_AK4CHS_geneta;   //!
   TBranch        *b_AK4CHS_genphi;   //!
   TBranch        *b_AK4CHS_genm;   //!
   TBranch        *b_AK4CHS_partonFlavor;   //!
   TBranch        *b_AK4CHS_hadronFlavor;   //!
   TBranch        *b_AK4CHS_nCharged;   //!
   TBranch        *b_AK4CHS_nNeutrals;   //!
   TBranch        *b_AK4CHS_nParticles;   //!
   TBranch        *b_AK8CHS_;   //!
   TBranch        *b_AK8CHS_pt;   //!
   TBranch        *b_AK8CHS_eta;   //!
   TBranch        *b_AK8CHS_phi;   //!
   TBranch        *b_AK8CHS_mass;   //!
   TBranch        *b_AK8CHS_ptRaw;   //!
   TBranch        *b_AK8CHS_unc;   //!
   TBranch        *b_AK8CHS_area;   //!
   TBranch        *b_AK8CHS_d0;   //!
   TBranch        *b_AK8CHS_dz;   //!
   TBranch        *b_AK8CHS_csv;   //!
   TBranch        *b_AK8CHS_bmva;   //!
   TBranch        *b_AK8CHS_cvb;   //!
   TBranch        *b_AK8CHS_cvl;   //!
   TBranch        *b_AK8CHS_qgid;   //!
   TBranch        *b_AK8CHS_axis2;   //!
   TBranch        *b_AK8CHS_ptD;   //!
   TBranch        *b_AK8CHS_mult;   //!
   TBranch        *b_AK8CHS_q;   //!
   TBranch        *b_AK8CHS_mva;   //!
   TBranch        *b_AK8CHS_beta;   //!
   TBranch        *b_AK8CHS_betaStar;   //!
   TBranch        *b_AK8CHS_dR2Mean;   //!
   TBranch        *b_AK8CHS_pullY;   //!
   TBranch        *b_AK8CHS_pullPhi;   //!
   TBranch        *b_AK8CHS_chPullY;   //!
   TBranch        *b_AK8CHS_chPullPhi;   //!
   TBranch        *b_AK8CHS_neuPullY;   //!
   TBranch        *b_AK8CHS_neuPullPhi;   //!
   TBranch        *b_AK8CHS_chEmFrac;   //!
   TBranch        *b_AK8CHS_neuEmFrac;   //!
   TBranch        *b_AK8CHS_chHadFrac;   //!
   TBranch        *b_AK8CHS_neuHadFrac;   //!
   TBranch        *b_AK8CHS_muonFrac;   //!
   TBranch        *b_AK8CHS_genpt;   //!
   TBranch        *b_AK8CHS_geneta;   //!
   TBranch        *b_AK8CHS_genphi;   //!
   TBranch        *b_AK8CHS_genm;   //!
   TBranch        *b_AK8CHS_partonFlavor;   //!
   TBranch        *b_AK8CHS_hadronFlavor;   //!
   TBranch        *b_AK8CHS_nCharged;   //!
   TBranch        *b_AK8CHS_nNeutrals;   //!
   TBranch        *b_AK8CHS_nParticles;   //!
   TBranch        *b_AddAK8CHS_;   //!
   TBranch        *b_AddAK8CHS_index;   //!
   TBranch        *b_AddAK8CHS_mass_prun;   //!
   TBranch        *b_AddAK8CHS_mass_trim;   //!
   TBranch        *b_AddAK8CHS_mass_sd0;   //!
   TBranch        *b_AddAK8CHS_c2_0;   //!
   TBranch        *b_AddAK8CHS_c2_0P2;   //!
   TBranch        *b_AddAK8CHS_c2_0P5;   //!
   TBranch        *b_AddAK8CHS_c2_1P0;   //!
   TBranch        *b_AddAK8CHS_c2_2P0;   //!
   TBranch        *b_AddAK8CHS_e2_b1;   //!
   TBranch        *b_AddAK8CHS_e3_b1;   //!
   TBranch        *b_AddAK8CHS_e3_v1_b1;   //!
   TBranch        *b_AddAK8CHS_e3_v2_b1;   //!
   TBranch        *b_AddAK8CHS_e4_v1_b1;   //!
   TBranch        *b_AddAK8CHS_e4_v2_b1;   //!
   TBranch        *b_AddAK8CHS_e2_b2;   //!
   TBranch        *b_AddAK8CHS_e3_b2;   //!
   TBranch        *b_AddAK8CHS_e3_v1_b2;   //!
   TBranch        *b_AddAK8CHS_e3_v2_b2;   //!
   TBranch        *b_AddAK8CHS_e4_v1_b2;   //!
   TBranch        *b_AddAK8CHS_e4_v2_b2;   //!
   TBranch        *b_AddAK8CHS_e2_sdb1;   //!
   TBranch        *b_AddAK8CHS_e3_sdb1;   //!
   TBranch        *b_AddAK8CHS_e3_v1_sdb1;   //!
   TBranch        *b_AddAK8CHS_e3_v2_sdb1;   //!
   TBranch        *b_AddAK8CHS_e4_v1_sdb1;   //!
   TBranch        *b_AddAK8CHS_e4_v2_sdb1;   //!
   TBranch        *b_AddAK8CHS_e2_sdb2;   //!
   TBranch        *b_AddAK8CHS_e3_sdb2;   //!
   TBranch        *b_AddAK8CHS_e3_v1_sdb2;   //!
   TBranch        *b_AddAK8CHS_e3_v2_sdb2;   //!
   TBranch        *b_AddAK8CHS_e4_v1_sdb2;   //!
   TBranch        *b_AddAK8CHS_e4_v2_sdb2;   //!
   TBranch        *b_AddAK8CHS_qjet;   //!
   TBranch        *b_AddAK8CHS_tau1;   //!
   TBranch        *b_AddAK8CHS_tau2;   //!
   TBranch        *b_AddAK8CHS_tau3;   //!
   TBranch        *b_AddAK8CHS_tau4;   //!
   TBranch        *b_AddAK8CHS_doublecsv;   //!
   TBranch        *b_AddAK8CHS_Double_sub;   //!
   TBranch        *b_AddAK8CHS_sj1_pt;   //!
   TBranch        *b_AddAK8CHS_sj1_eta;   //!
   TBranch        *b_AddAK8CHS_sj1_phi;   //!
   TBranch        *b_AddAK8CHS_sj1_m;   //!
   TBranch        *b_AddAK8CHS_sj1_csv;   //!
   TBranch        *b_AddAK8CHS_sj1_qgid;   //!
   TBranch        *b_AddAK8CHS_sj1_q;   //!
   TBranch        *b_AddAK8CHS_sj2_pt;   //!
   TBranch        *b_AddAK8CHS_sj2_eta;   //!
   TBranch        *b_AddAK8CHS_sj2_phi;   //!
   TBranch        *b_AddAK8CHS_sj2_m;   //!
   TBranch        *b_AddAK8CHS_sj2_csv;   //!
   TBranch        *b_AddAK8CHS_sj2_qgid;   //!
   TBranch        *b_AddAK8CHS_sj2_q;   //!
   TBranch        *b_AddAK8CHS_sj3_pt;   //!
   TBranch        *b_AddAK8CHS_sj3_eta;   //!
   TBranch        *b_AddAK8CHS_sj3_phi;   //!
   TBranch        *b_AddAK8CHS_sj3_m;   //!
   TBranch        *b_AddAK8CHS_sj3_csv;   //!
   TBranch        *b_AddAK8CHS_sj3_qgid;   //!
   TBranch        *b_AddAK8CHS_sj3_q;   //!
   TBranch        *b_AddAK8CHS_sj4_pt;   //!
   TBranch        *b_AddAK8CHS_sj4_eta;   //!
   TBranch        *b_AddAK8CHS_sj4_phi;   //!
   TBranch        *b_AddAK8CHS_sj4_m;   //!
   TBranch        *b_AddAK8CHS_sj4_csv;   //!
   TBranch        *b_AddAK8CHS_sj4_qgid;   //!
   TBranch        *b_AddAK8CHS_sj4_q;   //!
   TBranch        *b_AddAK8CHS_pullAngle;   //!
   TBranch        *b_AddAK8CHS_topTagType;   //!
   TBranch        *b_AddAK8CHS_top_n_subjets;   //!
   TBranch        *b_AddAK8CHS_top_m_min;   //!
   TBranch        *b_AddAK8CHS_top_m_123;   //!
   TBranch        *b_AddAK8CHS_top_fRec;   //!
   TBranch        *b_AddAK8CHS_topchi2;   //!
   TBranch        *b_AK4Puppi_;   //!
   TBranch        *b_AK4Puppi_pt;   //!
   TBranch        *b_AK4Puppi_eta;   //!
   TBranch        *b_AK4Puppi_phi;   //!
   TBranch        *b_AK4Puppi_mass;   //!
   TBranch        *b_AK4Puppi_ptRaw;   //!
   TBranch        *b_AK4Puppi_unc;   //!
   TBranch        *b_AK4Puppi_area;   //!
   TBranch        *b_AK4Puppi_d0;   //!
   TBranch        *b_AK4Puppi_dz;   //!
   TBranch        *b_AK4Puppi_csv;   //!
   TBranch        *b_AK4Puppi_bmva;   //!
   TBranch        *b_AK4Puppi_cvb;   //!
   TBranch        *b_AK4Puppi_cvl;   //!
   TBranch        *b_AK4Puppi_qgid;   //!
   TBranch        *b_AK4Puppi_axis2;   //!
   TBranch        *b_AK4Puppi_ptD;   //!
   TBranch        *b_AK4Puppi_mult;   //!
   TBranch        *b_AK4Puppi_q;   //!
   TBranch        *b_AK4Puppi_mva;   //!
   TBranch        *b_AK4Puppi_beta;   //!
   TBranch        *b_AK4Puppi_betaStar;   //!
   TBranch        *b_AK4Puppi_dR2Mean;   //!
   TBranch        *b_AK4Puppi_pullY;   //!
   TBranch        *b_AK4Puppi_pullPhi;   //!
   TBranch        *b_AK4Puppi_chPullY;   //!
   TBranch        *b_AK4Puppi_chPullPhi;   //!
   TBranch        *b_AK4Puppi_neuPullY;   //!
   TBranch        *b_AK4Puppi_neuPullPhi;   //!
   TBranch        *b_AK4Puppi_chEmFrac;   //!
   TBranch        *b_AK4Puppi_neuEmFrac;   //!
   TBranch        *b_AK4Puppi_chHadFrac;   //!
   TBranch        *b_AK4Puppi_neuHadFrac;   //!
   TBranch        *b_AK4Puppi_muonFrac;   //!
   TBranch        *b_AK4Puppi_genpt;   //!
   TBranch        *b_AK4Puppi_geneta;   //!
   TBranch        *b_AK4Puppi_genphi;   //!
   TBranch        *b_AK4Puppi_genm;   //!
   TBranch        *b_AK4Puppi_partonFlavor;   //!
   TBranch        *b_AK4Puppi_hadronFlavor;   //!
   TBranch        *b_AK4Puppi_nCharged;   //!
   TBranch        *b_AK4Puppi_nNeutrals;   //!
   TBranch        *b_AK4Puppi_nParticles;   //!
   TBranch        *b_AK8Puppi_;   //!
   TBranch        *b_AK8Puppi_pt;   //!
   TBranch        *b_AK8Puppi_eta;   //!
   TBranch        *b_AK8Puppi_phi;   //!
   TBranch        *b_AK8Puppi_mass;   //!
   TBranch        *b_AK8Puppi_ptRaw;   //!
   TBranch        *b_AK8Puppi_unc;   //!
   TBranch        *b_AK8Puppi_area;   //!
   TBranch        *b_AK8Puppi_d0;   //!
   TBranch        *b_AK8Puppi_dz;   //!
   TBranch        *b_AK8Puppi_csv;   //!
   TBranch        *b_AK8Puppi_bmva;   //!
   TBranch        *b_AK8Puppi_cvb;   //!
   TBranch        *b_AK8Puppi_cvl;   //!
   TBranch        *b_AK8Puppi_qgid;   //!
   TBranch        *b_AK8Puppi_axis2;   //!
   TBranch        *b_AK8Puppi_ptD;   //!
   TBranch        *b_AK8Puppi_mult;   //!
   TBranch        *b_AK8Puppi_q;   //!
   TBranch        *b_AK8Puppi_mva;   //!
   TBranch        *b_AK8Puppi_beta;   //!
   TBranch        *b_AK8Puppi_betaStar;   //!
   TBranch        *b_AK8Puppi_dR2Mean;   //!
   TBranch        *b_AK8Puppi_pullY;   //!
   TBranch        *b_AK8Puppi_pullPhi;   //!
   TBranch        *b_AK8Puppi_chPullY;   //!
   TBranch        *b_AK8Puppi_chPullPhi;   //!
   TBranch        *b_AK8Puppi_neuPullY;   //!
   TBranch        *b_AK8Puppi_neuPullPhi;   //!
   TBranch        *b_AK8Puppi_chEmFrac;   //!
   TBranch        *b_AK8Puppi_neuEmFrac;   //!
   TBranch        *b_AK8Puppi_chHadFrac;   //!
   TBranch        *b_AK8Puppi_neuHadFrac;   //!
   TBranch        *b_AK8Puppi_muonFrac;   //!
   TBranch        *b_AK8Puppi_genpt;   //!
   TBranch        *b_AK8Puppi_geneta;   //!
   TBranch        *b_AK8Puppi_genphi;   //!
   TBranch        *b_AK8Puppi_genm;   //!
   TBranch        *b_AK8Puppi_partonFlavor;   //!
   TBranch        *b_AK8Puppi_hadronFlavor;   //!
   TBranch        *b_AK8Puppi_nCharged;   //!
   TBranch        *b_AK8Puppi_nNeutrals;   //!
   TBranch        *b_AK8Puppi_nParticles;   //!
   TBranch        *b_AddAK8Puppi_;   //!
   TBranch        *b_AddAK8Puppi_index;   //!
   TBranch        *b_AddAK8Puppi_mass_prun;   //!
   TBranch        *b_AddAK8Puppi_mass_trim;   //!
   TBranch        *b_AddAK8Puppi_mass_sd0;   //!
   TBranch        *b_AddAK8Puppi_c2_0;   //!
   TBranch        *b_AddAK8Puppi_c2_0P2;   //!
   TBranch        *b_AddAK8Puppi_c2_0P5;   //!
   TBranch        *b_AddAK8Puppi_c2_1P0;   //!
   TBranch        *b_AddAK8Puppi_c2_2P0;   //!
   TBranch        *b_AddAK8Puppi_e2_b1;   //!
   TBranch        *b_AddAK8Puppi_e3_b1;   //!
   TBranch        *b_AddAK8Puppi_e3_v1_b1;   //!
   TBranch        *b_AddAK8Puppi_e3_v2_b1;   //!
   TBranch        *b_AddAK8Puppi_e4_v1_b1;   //!
   TBranch        *b_AddAK8Puppi_e4_v2_b1;   //!
   TBranch        *b_AddAK8Puppi_e2_b2;   //!
   TBranch        *b_AddAK8Puppi_e3_b2;   //!
   TBranch        *b_AddAK8Puppi_e3_v1_b2;   //!
   TBranch        *b_AddAK8Puppi_e3_v2_b2;   //!
   TBranch        *b_AddAK8Puppi_e4_v1_b2;   //!
   TBranch        *b_AddAK8Puppi_e4_v2_b2;   //!
   TBranch        *b_AddAK8Puppi_e2_sdb1;   //!
   TBranch        *b_AddAK8Puppi_e3_sdb1;   //!
   TBranch        *b_AddAK8Puppi_e3_v1_sdb1;   //!
   TBranch        *b_AddAK8Puppi_e3_v2_sdb1;   //!
   TBranch        *b_AddAK8Puppi_e4_v1_sdb1;   //!
   TBranch        *b_AddAK8Puppi_e4_v2_sdb1;   //!
   TBranch        *b_AddAK8Puppi_e2_sdb2;   //!
   TBranch        *b_AddAK8Puppi_e3_sdb2;   //!
   TBranch        *b_AddAK8Puppi_e3_v1_sdb2;   //!
   TBranch        *b_AddAK8Puppi_e3_v2_sdb2;   //!
   TBranch        *b_AddAK8Puppi_e4_v1_sdb2;   //!
   TBranch        *b_AddAK8Puppi_e4_v2_sdb2;   //!
   TBranch        *b_AddAK8Puppi_qjet;   //!
   TBranch        *b_AddAK8Puppi_tau1;   //!
   TBranch        *b_AddAK8Puppi_tau2;   //!
   TBranch        *b_AddAK8Puppi_tau3;   //!
   TBranch        *b_AddAK8Puppi_tau4;   //!
   TBranch        *b_AddAK8Puppi_doublecsv;   //!
   TBranch        *b_AddAK8Puppi_Double_sub;   //!
   TBranch        *b_AddAK8Puppi_sj1_pt;   //!
   TBranch        *b_AddAK8Puppi_sj1_eta;   //!
   TBranch        *b_AddAK8Puppi_sj1_phi;   //!
   TBranch        *b_AddAK8Puppi_sj1_m;   //!
   TBranch        *b_AddAK8Puppi_sj1_csv;   //!
   TBranch        *b_AddAK8Puppi_sj1_qgid;   //!
   TBranch        *b_AddAK8Puppi_sj1_q;   //!
   TBranch        *b_AddAK8Puppi_sj2_pt;   //!
   TBranch        *b_AddAK8Puppi_sj2_eta;   //!
   TBranch        *b_AddAK8Puppi_sj2_phi;   //!
   TBranch        *b_AddAK8Puppi_sj2_m;   //!
   TBranch        *b_AddAK8Puppi_sj2_csv;   //!
   TBranch        *b_AddAK8Puppi_sj2_qgid;   //!
   TBranch        *b_AddAK8Puppi_sj2_q;   //!
   TBranch        *b_AddAK8Puppi_sj3_pt;   //!
   TBranch        *b_AddAK8Puppi_sj3_eta;   //!
   TBranch        *b_AddAK8Puppi_sj3_phi;   //!
   TBranch        *b_AddAK8Puppi_sj3_m;   //!
   TBranch        *b_AddAK8Puppi_sj3_csv;   //!
   TBranch        *b_AddAK8Puppi_sj3_qgid;   //!
   TBranch        *b_AddAK8Puppi_sj3_q;   //!
   TBranch        *b_AddAK8Puppi_sj4_pt;   //!
   TBranch        *b_AddAK8Puppi_sj4_eta;   //!
   TBranch        *b_AddAK8Puppi_sj4_phi;   //!
   TBranch        *b_AddAK8Puppi_sj4_m;   //!
   TBranch        *b_AddAK8Puppi_sj4_csv;   //!
   TBranch        *b_AddAK8Puppi_sj4_qgid;   //!
   TBranch        *b_AddAK8Puppi_sj4_q;   //!
   TBranch        *b_AddAK8Puppi_pullAngle;   //!
   TBranch        *b_AddAK8Puppi_topTagType;   //!
   TBranch        *b_AddAK8Puppi_top_n_subjets;   //!
   TBranch        *b_AddAK8Puppi_top_m_min;   //!
   TBranch        *b_AddAK8Puppi_top_m_123;   //!
   TBranch        *b_AddAK8Puppi_top_fRec;   //!
   TBranch        *b_AddAK8Puppi_topchi2;   //!

   GetNegEvents(TTree *tree=0);
   virtual ~GetNegEvents();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop(std::string outFileName);
   //virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef GetNegEvents_cxx
GetNegEvents::GetNegEvents(TTree *tree) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("/eos/uscms/store/user/lnujj/WpWm_aQGC_Ntuples_Ram/FirstStepOutput/BaconNtuples/WminusToLNuZTo2JJJ_EWK_LO_aQGC_MJJ100PTJ10_TuneCUETP8M1_13TeV-madgraph-pythia8/Output_1.root");
      if (!f || !f->IsOpen()) {
         f = new TFile("/eos/uscms/store/user/lnujj/WpWm_aQGC_Ntuples_Ram/FirstStepOutput/BaconNtuples/WminusToLNuZTo2JJJ_EWK_LO_aQGC_MJJ100PTJ10_TuneCUETP8M1_13TeV-madgraph-pythia8/Output_1.root");
      }
      f->GetObject("Events",tree);

   }
   Init(tree);
}

GetNegEvents::~GetNegEvents()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t GetNegEvents::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t GetNegEvents::LoadTree(Long64_t entry)
{
// Set the environment to read one entry
   if (!fChain) return -5;
   Long64_t centry = fChain->LoadTree(entry);
   if (centry < 0) return centry;
   if (fChain->GetTreeNumber() != fCurrent) {
      fCurrent = fChain->GetTreeNumber();
      Notify();
   }
   return centry;
}

void GetNegEvents::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("runNum", &runNum, &b_Info_runNum);
   fChain->SetBranchAddress("evtNum", &evtNum, &b_Info_evtNum);
   fChain->SetBranchAddress("lumiSec", &lumiSec, &b_Info_lumiSec);
   fChain->SetBranchAddress("metFilterFailBits", &metFilterFailBits, &b_Info_metFilterFailBits);
   fChain->SetBranchAddress("nPU", &nPU, &b_Info_nPU);
   fChain->SetBranchAddress("nPUm", &nPUm, &b_Info_nPUm);
   fChain->SetBranchAddress("nPUp", &nPUp, &b_Info_nPUp);
   fChain->SetBranchAddress("nPUmean", &nPUmean, &b_Info_nPUmean);
   fChain->SetBranchAddress("nPUmeanm", &nPUmeanm, &b_Info_nPUmeanm);
   fChain->SetBranchAddress("nPUmeanp", &nPUmeanp, &b_Info_nPUmeanp);
   fChain->SetBranchAddress("pvx", &pvx, &b_Info_pvx);
   fChain->SetBranchAddress("pvy", &pvy, &b_Info_pvy);
   fChain->SetBranchAddress("pvz", &pvz, &b_Info_pvz);
   fChain->SetBranchAddress("bsx", &bsx, &b_Info_bsx);
   fChain->SetBranchAddress("bsy", &bsy, &b_Info_bsy);
   fChain->SetBranchAddress("bsz", &bsz, &b_Info_bsz);
   fChain->SetBranchAddress("caloMET", &caloMET, &b_Info_caloMET);
   fChain->SetBranchAddress("caloMETphi", &caloMETphi, &b_Info_caloMETphi);
   fChain->SetBranchAddress("pfMET", &pfMET, &b_Info_pfMET);
   fChain->SetBranchAddress("pfMETphi", &pfMETphi, &b_Info_pfMETphi);
   fChain->SetBranchAddress("pfMETCov00", &pfMETCov00, &b_Info_pfMETCov00);
   fChain->SetBranchAddress("pfMETCov01", &pfMETCov01, &b_Info_pfMETCov01);
   fChain->SetBranchAddress("pfMETCov11", &pfMETCov11, &b_Info_pfMETCov11);
   fChain->SetBranchAddress("pfMETC", &pfMETC, &b_Info_pfMETC);
   fChain->SetBranchAddress("pfMETCphi", &pfMETCphi, &b_Info_pfMETCphi);
   fChain->SetBranchAddress("pfMETCCov00", &pfMETCCov00, &b_Info_pfMETCCov00);
   fChain->SetBranchAddress("pfMETCCov01", &pfMETCCov01, &b_Info_pfMETCCov01);
   fChain->SetBranchAddress("pfMETCCov11", &pfMETCCov11, &b_Info_pfMETCCov11);
   fChain->SetBranchAddress("pfMETCjerup", &pfMETCjerup, &b_Info_pfMETCjerup);
   fChain->SetBranchAddress("pfMETCjerdn", &pfMETCjerdn, &b_Info_pfMETCjerdn);
   fChain->SetBranchAddress("pfMETCjenup", &pfMETCjenup, &b_Info_pfMETCjenup);
   fChain->SetBranchAddress("pfMETCjendn", &pfMETCjendn, &b_Info_pfMETCjendn);
   fChain->SetBranchAddress("pfMETCuncup", &pfMETCuncup, &b_Info_pfMETCuncup);
   fChain->SetBranchAddress("pfMETCuncdn", &pfMETCuncdn, &b_Info_pfMETCuncdn);
   fChain->SetBranchAddress("pfMETCjrsup", &pfMETCjrsup, &b_Info_pfMETCjrsup);
   fChain->SetBranchAddress("pfMETCjrsdn", &pfMETCjrsdn, &b_Info_pfMETCjrsdn);
   fChain->SetBranchAddress("pfMETCphijerup", &pfMETCphijerup, &b_Info_pfMETCphijerup);
   fChain->SetBranchAddress("pfMETCphijerdn", &pfMETCphijerdn, &b_Info_pfMETCphijerdn);
   fChain->SetBranchAddress("pfMETCphijenup", &pfMETCphijenup, &b_Info_pfMETCphijenup);
   fChain->SetBranchAddress("pfMETCphijendn", &pfMETCphijendn, &b_Info_pfMETCphijendn);
   fChain->SetBranchAddress("pfMETCphiuncup", &pfMETCphiuncup, &b_Info_pfMETCphiuncup);
   fChain->SetBranchAddress("pfMETCphiuncdn", &pfMETCphiuncdn, &b_Info_pfMETCphiuncdn);
   fChain->SetBranchAddress("pfMETCphijrsup", &pfMETCphijrsup, &b_Info_pfMETCphijrsup);
   fChain->SetBranchAddress("pfMETCphijrsdn", &pfMETCphijrsdn, &b_Info_pfMETCphijrsdn);
   fChain->SetBranchAddress("puppET", &puppET, &b_Info_puppET);
   fChain->SetBranchAddress("puppETphi", &puppETphi, &b_Info_puppETphi);
   fChain->SetBranchAddress("puppETCov00", &puppETCov00, &b_Info_puppETCov00);
   fChain->SetBranchAddress("puppETCov01", &puppETCov01, &b_Info_puppETCov01);
   fChain->SetBranchAddress("puppETCov11", &puppETCov11, &b_Info_puppETCov11);
   fChain->SetBranchAddress("puppETC", &puppETC, &b_Info_puppETC);
   fChain->SetBranchAddress("puppETCphi", &puppETCphi, &b_Info_puppETCphi);
   fChain->SetBranchAddress("puppETCCov00", &puppETCCov00, &b_Info_puppETCCov00);
   fChain->SetBranchAddress("puppETCCov01", &puppETCCov01, &b_Info_puppETCCov01);
   fChain->SetBranchAddress("puppETCCov11", &puppETCCov11, &b_Info_puppETCCov11);
   fChain->SetBranchAddress("alpacaMET", &alpacaMET, &b_Info_alpacaMET);
   fChain->SetBranchAddress("alpacaMETphi", &alpacaMETphi, &b_Info_alpacaMETphi);
   fChain->SetBranchAddress("pcpMET", &pcpMET, &b_Info_pcpMET);
   fChain->SetBranchAddress("pcpMETphi", &pcpMETphi, &b_Info_pcpMETphi);
   fChain->SetBranchAddress("trkMET", &trkMET, &b_Info_trkMET);
   fChain->SetBranchAddress("trkMETphi", &trkMETphi, &b_Info_trkMETphi);
   fChain->SetBranchAddress("rhoIso", &rhoIso, &b_Info_rhoIso);
   fChain->SetBranchAddress("rhoJet", &rhoJet, &b_Info_rhoJet);
   fChain->SetBranchAddress("hasGoodPV", &hasGoodPV, &b_Info_hasGoodPV);
   fChain->SetBranchAddress("id_1", &id_1, &b_GenEvtInfo_id_1);
   fChain->SetBranchAddress("id_2", &id_2, &b_GenEvtInfo_id_2);
   fChain->SetBranchAddress("x_1", &x_1, &b_GenEvtInfo_x_1);
   fChain->SetBranchAddress("x_2", &x_2, &b_GenEvtInfo_x_2);
   fChain->SetBranchAddress("scalePDF", &scalePDF, &b_GenEvtInfo_scalePDF);
   fChain->SetBranchAddress("xs", &xs, &b_GenEvtInfo_xs);
   fChain->SetBranchAddress("weight", &weight, &b_GenEvtInfo_weight);
   fChain->SetBranchAddress("GenParticle", &GenParticle_, &b_GenParticle_);
   fChain->SetBranchAddress("GenParticle.parent", GenParticle_parent, &b_GenParticle_parent);
   fChain->SetBranchAddress("GenParticle.pdgId", GenParticle_pdgId, &b_GenParticle_pdgId);
   fChain->SetBranchAddress("GenParticle.status", GenParticle_status, &b_GenParticle_status);
   fChain->SetBranchAddress("GenParticle.pt", GenParticle_pt, &b_GenParticle_pt);
   fChain->SetBranchAddress("GenParticle.eta", GenParticle_eta, &b_GenParticle_eta);
   fChain->SetBranchAddress("GenParticle.phi", GenParticle_phi, &b_GenParticle_phi);
   fChain->SetBranchAddress("GenParticle.mass", GenParticle_mass, &b_GenParticle_mass);
   fChain->SetBranchAddress("GenParticle.y", GenParticle_y, &b_GenParticle_y);
   //fChain->SetBranchAddress("LHEWeight", &LHEWeight_, &b_LHEWeight_);
   //fChain->SetBranchAddress("LHEWeight.id", LHEWeight_id, &b_LHEWeight_id);
   //fChain->SetBranchAddress("LHEWeight.weight", LHEWeight_weight, &b_LHEWeight_weight);
   fChain->SetBranchAddress("Electron", &Electron_, &b_Electron_);
   fChain->SetBranchAddress("Electron.pt", Electron_pt, &b_Electron_pt);
   fChain->SetBranchAddress("Electron.eta", Electron_eta, &b_Electron_eta);
   fChain->SetBranchAddress("Electron.phi", Electron_phi, &b_Electron_phi);
   fChain->SetBranchAddress("Electron.scEt", Electron_scEt, &b_Electron_scEt);
   fChain->SetBranchAddress("Electron.scEta", Electron_scEta, &b_Electron_scEta);
   fChain->SetBranchAddress("Electron.scPhi", Electron_scPhi, &b_Electron_scPhi);
   fChain->SetBranchAddress("Electron.ecalEnergy", Electron_ecalEnergy, &b_Electron_ecalEnergy);
   fChain->SetBranchAddress("Electron.pfPt", Electron_pfPt, &b_Electron_pfPt);
   fChain->SetBranchAddress("Electron.pfEta", Electron_pfEta, &b_Electron_pfEta);
   fChain->SetBranchAddress("Electron.pfPhi", Electron_pfPhi, &b_Electron_pfPhi);
   fChain->SetBranchAddress("Electron.trkIso", Electron_trkIso, &b_Electron_trkIso);
   fChain->SetBranchAddress("Electron.ecalIso", Electron_ecalIso, &b_Electron_ecalIso);
   fChain->SetBranchAddress("Electron.hcalIso", Electron_hcalIso, &b_Electron_hcalIso);
   fChain->SetBranchAddress("Electron.hcalDepth1Iso", Electron_hcalDepth1Iso, &b_Electron_hcalDepth1Iso);
   fChain->SetBranchAddress("Electron.chHadIso", Electron_chHadIso, &b_Electron_chHadIso);
   fChain->SetBranchAddress("Electron.gammaIso", Electron_gammaIso, &b_Electron_gammaIso);
   fChain->SetBranchAddress("Electron.neuHadIso", Electron_neuHadIso, &b_Electron_neuHadIso);
   fChain->SetBranchAddress("Electron.puIso", Electron_puIso, &b_Electron_puIso);
   fChain->SetBranchAddress("Electron.ecalPFClusIso", Electron_ecalPFClusIso, &b_Electron_ecalPFClusIso);
   fChain->SetBranchAddress("Electron.hcalPFClusIso", Electron_hcalPFClusIso, &b_Electron_hcalPFClusIso);
   fChain->SetBranchAddress("Electron.puppiChHadIso", Electron_puppiChHadIso, &b_Electron_puppiChHadIso);
   fChain->SetBranchAddress("Electron.puppiGammaIso", Electron_puppiGammaIso, &b_Electron_puppiGammaIso);
   fChain->SetBranchAddress("Electron.puppiNeuHadIso", Electron_puppiNeuHadIso, &b_Electron_puppiNeuHadIso);
   fChain->SetBranchAddress("Electron.puppiChHadIsoNoLep", Electron_puppiChHadIsoNoLep, &b_Electron_puppiChHadIsoNoLep);
   fChain->SetBranchAddress("Electron.puppiGammaIsoNoLep", Electron_puppiGammaIsoNoLep, &b_Electron_puppiGammaIsoNoLep);
   fChain->SetBranchAddress("Electron.puppiNeuHadIsoNoLep", Electron_puppiNeuHadIsoNoLep, &b_Electron_puppiNeuHadIsoNoLep);
   fChain->SetBranchAddress("Electron.d0", Electron_d0, &b_Electron_d0);
   fChain->SetBranchAddress("Electron.dz", Electron_dz, &b_Electron_dz);
   fChain->SetBranchAddress("Electron.sip3d", Electron_sip3d, &b_Electron_sip3d);
   fChain->SetBranchAddress("Electron.sieie", Electron_sieie, &b_Electron_sieie);
   fChain->SetBranchAddress("Electron.e1x5", Electron_e1x5, &b_Electron_e1x5);
   fChain->SetBranchAddress("Electron.e2x5", Electron_e2x5, &b_Electron_e2x5);
   fChain->SetBranchAddress("Electron.e5x5", Electron_e5x5, &b_Electron_e5x5);
   fChain->SetBranchAddress("Electron.r9", Electron_r9, &b_Electron_r9);
   fChain->SetBranchAddress("Electron.eoverp", Electron_eoverp, &b_Electron_eoverp);
   fChain->SetBranchAddress("Electron.hovere", Electron_hovere, &b_Electron_hovere);
   fChain->SetBranchAddress("Electron.fbrem", Electron_fbrem, &b_Electron_fbrem);
   fChain->SetBranchAddress("Electron.dEtaInSeed", Electron_dEtaInSeed, &b_Electron_dEtaInSeed);
   fChain->SetBranchAddress("Electron.dEtaIn", Electron_dEtaIn, &b_Electron_dEtaIn);
   fChain->SetBranchAddress("Electron.dPhiIn", Electron_dPhiIn, &b_Electron_dPhiIn);
   fChain->SetBranchAddress("Electron.mva", Electron_mva, &b_Electron_mva);
   fChain->SetBranchAddress("Electron.q", Electron_q, &b_Electron_q);
   fChain->SetBranchAddress("Electron.isConv", Electron_isConv, &b_Electron_isConv);
   fChain->SetBranchAddress("Electron.nMissingHits", Electron_nMissingHits, &b_Electron_nMissingHits);
   fChain->SetBranchAddress("Electron.typeBits", Electron_typeBits, &b_Electron_typeBits);
   fChain->SetBranchAddress("Electron.fiducialBits", Electron_fiducialBits, &b_Electron_fiducialBits);
   fChain->SetBranchAddress("Electron.classification", Electron_classification, &b_Electron_classification);
   fChain->SetBranchAddress("Electron.scID", Electron_scID, &b_Electron_scID);
   fChain->SetBranchAddress("Electron.trkID", Electron_trkID, &b_Electron_trkID);
   fChain->SetBranchAddress("Muon", &Muon_, &b_Muon_);
   fChain->SetBranchAddress("Muon.pt", Muon_pt, &b_Muon_pt);
   fChain->SetBranchAddress("Muon.eta", Muon_eta, &b_Muon_eta);
   fChain->SetBranchAddress("Muon.phi", Muon_phi, &b_Muon_phi);
   fChain->SetBranchAddress("Muon.ptErr", Muon_ptErr, &b_Muon_ptErr);
   fChain->SetBranchAddress("Muon.staPt", Muon_staPt, &b_Muon_staPt);
   fChain->SetBranchAddress("Muon.staEta", Muon_staEta, &b_Muon_staEta);
   fChain->SetBranchAddress("Muon.staPhi", Muon_staPhi, &b_Muon_staPhi);
   fChain->SetBranchAddress("Muon.pfPt", Muon_pfPt, &b_Muon_pfPt);
   fChain->SetBranchAddress("Muon.pfEta", Muon_pfEta, &b_Muon_pfEta);
   fChain->SetBranchAddress("Muon.pfPhi", Muon_pfPhi, &b_Muon_pfPhi);
   fChain->SetBranchAddress("Muon.trkIso", Muon_trkIso, &b_Muon_trkIso);
   fChain->SetBranchAddress("Muon.ecalIso", Muon_ecalIso, &b_Muon_ecalIso);
   fChain->SetBranchAddress("Muon.hcalIso", Muon_hcalIso, &b_Muon_hcalIso);
   fChain->SetBranchAddress("Muon.chHadIso", Muon_chHadIso, &b_Muon_chHadIso);
   fChain->SetBranchAddress("Muon.gammaIso", Muon_gammaIso, &b_Muon_gammaIso);
   fChain->SetBranchAddress("Muon.neuHadIso", Muon_neuHadIso, &b_Muon_neuHadIso);
   fChain->SetBranchAddress("Muon.puIso", Muon_puIso, &b_Muon_puIso);
   fChain->SetBranchAddress("Muon.puppiChHadIso", Muon_puppiChHadIso, &b_Muon_puppiChHadIso);
   fChain->SetBranchAddress("Muon.puppiGammaIso", Muon_puppiGammaIso, &b_Muon_puppiGammaIso);
   fChain->SetBranchAddress("Muon.puppiNeuHadIso", Muon_puppiNeuHadIso, &b_Muon_puppiNeuHadIso);
   fChain->SetBranchAddress("Muon.puppiChHadIsoNoLep", Muon_puppiChHadIsoNoLep, &b_Muon_puppiChHadIsoNoLep);
   fChain->SetBranchAddress("Muon.puppiGammaIsoNoLep", Muon_puppiGammaIsoNoLep, &b_Muon_puppiGammaIsoNoLep);
   fChain->SetBranchAddress("Muon.puppiNeuHadIsoNoLep", Muon_puppiNeuHadIsoNoLep, &b_Muon_puppiNeuHadIsoNoLep);
   fChain->SetBranchAddress("Muon.d0", Muon_d0, &b_Muon_d0);
   fChain->SetBranchAddress("Muon.dz", Muon_dz, &b_Muon_dz);
   fChain->SetBranchAddress("Muon.sip3d", Muon_sip3d, &b_Muon_sip3d);
   fChain->SetBranchAddress("Muon.tkNchi2", Muon_tkNchi2, &b_Muon_tkNchi2);
   fChain->SetBranchAddress("Muon.muNchi2", Muon_muNchi2, &b_Muon_muNchi2);
   fChain->SetBranchAddress("Muon.trkKink", Muon_trkKink, &b_Muon_trkKink);
   fChain->SetBranchAddress("Muon.glbKink", Muon_glbKink, &b_Muon_glbKink);
   fChain->SetBranchAddress("Muon.trkHitFrac", Muon_trkHitFrac, &b_Muon_trkHitFrac);
   fChain->SetBranchAddress("Muon.chi2LocPos", Muon_chi2LocPos, &b_Muon_chi2LocPos);
   fChain->SetBranchAddress("Muon.segComp", Muon_segComp, &b_Muon_segComp);
   fChain->SetBranchAddress("Muon.caloComp", Muon_caloComp, &b_Muon_caloComp);
   fChain->SetBranchAddress("Muon.q", Muon_q, &b_Muon_q);
   fChain->SetBranchAddress("Muon.nValidHits", Muon_nValidHits, &b_Muon_nValidHits);
   fChain->SetBranchAddress("Muon.typeBits", Muon_typeBits, &b_Muon_typeBits);
   fChain->SetBranchAddress("Muon.selectorBits", Muon_selectorBits, &b_Muon_selectorBits);
   fChain->SetBranchAddress("Muon.pogIDBits", Muon_pogIDBits, &b_Muon_pogIDBits);
   fChain->SetBranchAddress("Muon.nTkHits", Muon_nTkHits, &b_Muon_nTkHits);
   fChain->SetBranchAddress("Muon.nPixHits", Muon_nPixHits, &b_Muon_nPixHits);
   fChain->SetBranchAddress("Muon.nTkLayers", Muon_nTkLayers, &b_Muon_nTkLayers);
   fChain->SetBranchAddress("Muon.nPixLayers", Muon_nPixLayers, &b_Muon_nPixLayers);
   fChain->SetBranchAddress("Muon.nMatchStn", Muon_nMatchStn, &b_Muon_nMatchStn);
   fChain->SetBranchAddress("Muon.trkID", Muon_trkID, &b_Muon_trkID);
   fChain->SetBranchAddress("Tau", &Tau_, &b_Tau_);
   fChain->SetBranchAddress("Tau.pt", Tau_pt, &b_Tau_pt);
   fChain->SetBranchAddress("Tau.eta", Tau_eta, &b_Tau_eta);
   fChain->SetBranchAddress("Tau.phi", Tau_phi, &b_Tau_phi);
   fChain->SetBranchAddress("Tau.m", Tau_m, &b_Tau_m);
   fChain->SetBranchAddress("Tau.e", Tau_e, &b_Tau_e);
   fChain->SetBranchAddress("Tau.q", Tau_q, &b_Tau_q);
   fChain->SetBranchAddress("Tau.dzLeadChHad", Tau_dzLeadChHad, &b_Tau_dzLeadChHad);
   fChain->SetBranchAddress("Tau.d0LeadChHad", Tau_d0LeadChHad, &b_Tau_d0LeadChHad);
   fChain->SetBranchAddress("Tau.nSignalChHad", Tau_nSignalChHad, &b_Tau_nSignalChHad);
   fChain->SetBranchAddress("Tau.nSignalGamma", Tau_nSignalGamma, &b_Tau_nSignalGamma);
   fChain->SetBranchAddress("Tau.decaymode", Tau_decaymode, &b_Tau_decaymode);
   fChain->SetBranchAddress("Tau.antiEleMVA6", Tau_antiEleMVA6, &b_Tau_antiEleMVA6);
   fChain->SetBranchAddress("Tau.antiEleMVA6Cat", Tau_antiEleMVA6Cat, &b_Tau_antiEleMVA6Cat);
   fChain->SetBranchAddress("Tau.rawMuonRejection", Tau_rawMuonRejection, &b_Tau_rawMuonRejection);
   fChain->SetBranchAddress("Tau.rawIso3Hits", Tau_rawIso3Hits, &b_Tau_rawIso3Hits);
   fChain->SetBranchAddress("Tau.rawIsoMVA3oldDMwoLT", Tau_rawIsoMVA3oldDMwoLT, &b_Tau_rawIsoMVA3oldDMwoLT);
   fChain->SetBranchAddress("Tau.rawIsoMVA3oldDMwLT", Tau_rawIsoMVA3oldDMwLT, &b_Tau_rawIsoMVA3oldDMwLT);
   fChain->SetBranchAddress("Tau.rawIsoMVA3newDMwoLT", Tau_rawIsoMVA3newDMwoLT, &b_Tau_rawIsoMVA3newDMwoLT);
   fChain->SetBranchAddress("Tau.rawIsoMVA3newDMwLT", Tau_rawIsoMVA3newDMwLT, &b_Tau_rawIsoMVA3newDMwLT);
   fChain->SetBranchAddress("Tau.puppiChHadIso", Tau_puppiChHadIso, &b_Tau_puppiChHadIso);
   fChain->SetBranchAddress("Tau.puppiGammaIso", Tau_puppiGammaIso, &b_Tau_puppiGammaIso);
   fChain->SetBranchAddress("Tau.puppiNeuHadIso", Tau_puppiNeuHadIso, &b_Tau_puppiNeuHadIso);
   fChain->SetBranchAddress("Tau.puppiChHadIsoNoLep", Tau_puppiChHadIsoNoLep, &b_Tau_puppiChHadIsoNoLep);
   fChain->SetBranchAddress("Tau.puppiGammaIsoNoLep", Tau_puppiGammaIsoNoLep, &b_Tau_puppiGammaIsoNoLep);
   fChain->SetBranchAddress("Tau.puppiNeuHadIsoNoLep", Tau_puppiNeuHadIsoNoLep, &b_Tau_puppiNeuHadIsoNoLep);
   fChain->SetBranchAddress("Tau.hpsDisc", Tau_hpsDisc, &b_Tau_hpsDisc);
   fChain->SetBranchAddress("Photon", &Photon_, &b_Photon_);
   fChain->SetBranchAddress("Photon.pt", Photon_pt, &b_Photon_pt);
   fChain->SetBranchAddress("Photon.eta", Photon_eta, &b_Photon_eta);
   fChain->SetBranchAddress("Photon.phi", Photon_phi, &b_Photon_phi);
   fChain->SetBranchAddress("Photon.scEt", Photon_scEt, &b_Photon_scEt);
   fChain->SetBranchAddress("Photon.scEta", Photon_scEta, &b_Photon_scEta);
   fChain->SetBranchAddress("Photon.scPhi", Photon_scPhi, &b_Photon_scPhi);
   fChain->SetBranchAddress("Photon.trkIso", Photon_trkIso, &b_Photon_trkIso);
   fChain->SetBranchAddress("Photon.ecalIso", Photon_ecalIso, &b_Photon_ecalIso);
   fChain->SetBranchAddress("Photon.hcalIso", Photon_hcalIso, &b_Photon_hcalIso);
   fChain->SetBranchAddress("Photon.chHadIso", Photon_chHadIso, &b_Photon_chHadIso);
   fChain->SetBranchAddress("Photon.gammaIso", Photon_gammaIso, &b_Photon_gammaIso);
   fChain->SetBranchAddress("Photon.neuHadIso", Photon_neuHadIso, &b_Photon_neuHadIso);
   fChain->SetBranchAddress("Photon.mva", Photon_mva, &b_Photon_mva);
   fChain->SetBranchAddress("Photon.hovere", Photon_hovere, &b_Photon_hovere);
   fChain->SetBranchAddress("Photon.sthovere", Photon_sthovere, &b_Photon_sthovere);
   fChain->SetBranchAddress("Photon.sieie", Photon_sieie, &b_Photon_sieie);
   fChain->SetBranchAddress("Photon.sipip", Photon_sipip, &b_Photon_sipip);
   fChain->SetBranchAddress("Photon.r9", Photon_r9, &b_Photon_r9);
   fChain->SetBranchAddress("Photon.fiducialBits", Photon_fiducialBits, &b_Photon_fiducialBits);
   fChain->SetBranchAddress("Photon.typeBits", Photon_typeBits, &b_Photon_typeBits);
   fChain->SetBranchAddress("Photon.scID", Photon_scID, &b_Photon_scID);
   fChain->SetBranchAddress("Photon.hasPixelSeed", Photon_hasPixelSeed, &b_Photon_hasPixelSeed);
   fChain->SetBranchAddress("Photon.passElectronVeto", Photon_passElectronVeto, &b_Photon_passElectronVeto);
   fChain->SetBranchAddress("Photon.isConv", Photon_isConv, &b_Photon_isConv);
   fChain->SetBranchAddress("PV", &PV_, &b_PV_);
   fChain->SetBranchAddress("PV.nTracksFit", PV_nTracksFit, &b_PV_nTracksFit);
   fChain->SetBranchAddress("PV.ndof", PV_ndof, &b_PV_ndof);
   fChain->SetBranchAddress("PV.chi2", PV_chi2, &b_PV_chi2);
   fChain->SetBranchAddress("PV.x", PV_x, &b_PV_x);
   fChain->SetBranchAddress("PV.y", PV_y, &b_PV_y);
   fChain->SetBranchAddress("PV.z", PV_z, &b_PV_z);
   fChain->SetBranchAddress("AK4CHS", &AK4CHS_, &b_AK4CHS_);
   fChain->SetBranchAddress("AK4CHS.pt", AK4CHS_pt, &b_AK4CHS_pt);
   fChain->SetBranchAddress("AK4CHS.eta", AK4CHS_eta, &b_AK4CHS_eta);
   fChain->SetBranchAddress("AK4CHS.phi", AK4CHS_phi, &b_AK4CHS_phi);
   fChain->SetBranchAddress("AK4CHS.mass", AK4CHS_mass, &b_AK4CHS_mass);
   fChain->SetBranchAddress("AK4CHS.ptRaw", AK4CHS_ptRaw, &b_AK4CHS_ptRaw);
   fChain->SetBranchAddress("AK4CHS.unc", AK4CHS_unc, &b_AK4CHS_unc);
   fChain->SetBranchAddress("AK4CHS.area", AK4CHS_area, &b_AK4CHS_area);
   fChain->SetBranchAddress("AK4CHS.d0", AK4CHS_d0, &b_AK4CHS_d0);
   fChain->SetBranchAddress("AK4CHS.dz", AK4CHS_dz, &b_AK4CHS_dz);
   fChain->SetBranchAddress("AK4CHS.csv", AK4CHS_csv, &b_AK4CHS_csv);
   fChain->SetBranchAddress("AK4CHS.bmva", AK4CHS_bmva, &b_AK4CHS_bmva);
   fChain->SetBranchAddress("AK4CHS.cvb", AK4CHS_cvb, &b_AK4CHS_cvb);
   fChain->SetBranchAddress("AK4CHS.cvl", AK4CHS_cvl, &b_AK4CHS_cvl);
   fChain->SetBranchAddress("AK4CHS.qgid", AK4CHS_qgid, &b_AK4CHS_qgid);
   fChain->SetBranchAddress("AK4CHS.axis2", AK4CHS_axis2, &b_AK4CHS_axis2);
   fChain->SetBranchAddress("AK4CHS.ptD", AK4CHS_ptD, &b_AK4CHS_ptD);
   fChain->SetBranchAddress("AK4CHS.mult", AK4CHS_mult, &b_AK4CHS_mult);
   fChain->SetBranchAddress("AK4CHS.q", AK4CHS_q, &b_AK4CHS_q);
   fChain->SetBranchAddress("AK4CHS.mva", AK4CHS_mva, &b_AK4CHS_mva);
   fChain->SetBranchAddress("AK4CHS.beta", AK4CHS_beta, &b_AK4CHS_beta);
   fChain->SetBranchAddress("AK4CHS.betaStar", AK4CHS_betaStar, &b_AK4CHS_betaStar);
   fChain->SetBranchAddress("AK4CHS.dR2Mean", AK4CHS_dR2Mean, &b_AK4CHS_dR2Mean);
   fChain->SetBranchAddress("AK4CHS.pullY", AK4CHS_pullY, &b_AK4CHS_pullY);
   fChain->SetBranchAddress("AK4CHS.pullPhi", AK4CHS_pullPhi, &b_AK4CHS_pullPhi);
   fChain->SetBranchAddress("AK4CHS.chPullY", AK4CHS_chPullY, &b_AK4CHS_chPullY);
   fChain->SetBranchAddress("AK4CHS.chPullPhi", AK4CHS_chPullPhi, &b_AK4CHS_chPullPhi);
   fChain->SetBranchAddress("AK4CHS.neuPullY", AK4CHS_neuPullY, &b_AK4CHS_neuPullY);
   fChain->SetBranchAddress("AK4CHS.neuPullPhi", AK4CHS_neuPullPhi, &b_AK4CHS_neuPullPhi);
   fChain->SetBranchAddress("AK4CHS.chEmFrac", AK4CHS_chEmFrac, &b_AK4CHS_chEmFrac);
   fChain->SetBranchAddress("AK4CHS.neuEmFrac", AK4CHS_neuEmFrac, &b_AK4CHS_neuEmFrac);
   fChain->SetBranchAddress("AK4CHS.chHadFrac", AK4CHS_chHadFrac, &b_AK4CHS_chHadFrac);
   fChain->SetBranchAddress("AK4CHS.neuHadFrac", AK4CHS_neuHadFrac, &b_AK4CHS_neuHadFrac);
   fChain->SetBranchAddress("AK4CHS.muonFrac", AK4CHS_muonFrac, &b_AK4CHS_muonFrac);
   fChain->SetBranchAddress("AK4CHS.genpt", AK4CHS_genpt, &b_AK4CHS_genpt);
   fChain->SetBranchAddress("AK4CHS.geneta", AK4CHS_geneta, &b_AK4CHS_geneta);
   fChain->SetBranchAddress("AK4CHS.genphi", AK4CHS_genphi, &b_AK4CHS_genphi);
   fChain->SetBranchAddress("AK4CHS.genm", AK4CHS_genm, &b_AK4CHS_genm);
   fChain->SetBranchAddress("AK4CHS.partonFlavor", AK4CHS_partonFlavor, &b_AK4CHS_partonFlavor);
   fChain->SetBranchAddress("AK4CHS.hadronFlavor", AK4CHS_hadronFlavor, &b_AK4CHS_hadronFlavor);
   fChain->SetBranchAddress("AK4CHS.nCharged", AK4CHS_nCharged, &b_AK4CHS_nCharged);
   fChain->SetBranchAddress("AK4CHS.nNeutrals", AK4CHS_nNeutrals, &b_AK4CHS_nNeutrals);
   fChain->SetBranchAddress("AK4CHS.nParticles", AK4CHS_nParticles, &b_AK4CHS_nParticles);
   fChain->SetBranchAddress("AK8CHS", &AK8CHS_, &b_AK8CHS_);
   fChain->SetBranchAddress("AK8CHS.pt", AK8CHS_pt, &b_AK8CHS_pt);
   fChain->SetBranchAddress("AK8CHS.eta", AK8CHS_eta, &b_AK8CHS_eta);
   fChain->SetBranchAddress("AK8CHS.phi", AK8CHS_phi, &b_AK8CHS_phi);
   fChain->SetBranchAddress("AK8CHS.mass", AK8CHS_mass, &b_AK8CHS_mass);
   fChain->SetBranchAddress("AK8CHS.ptRaw", AK8CHS_ptRaw, &b_AK8CHS_ptRaw);
   fChain->SetBranchAddress("AK8CHS.unc", AK8CHS_unc, &b_AK8CHS_unc);
   fChain->SetBranchAddress("AK8CHS.area", AK8CHS_area, &b_AK8CHS_area);
   fChain->SetBranchAddress("AK8CHS.d0", AK8CHS_d0, &b_AK8CHS_d0);
   fChain->SetBranchAddress("AK8CHS.dz", AK8CHS_dz, &b_AK8CHS_dz);
   fChain->SetBranchAddress("AK8CHS.csv", AK8CHS_csv, &b_AK8CHS_csv);
   fChain->SetBranchAddress("AK8CHS.bmva", AK8CHS_bmva, &b_AK8CHS_bmva);
   fChain->SetBranchAddress("AK8CHS.cvb", AK8CHS_cvb, &b_AK8CHS_cvb);
   fChain->SetBranchAddress("AK8CHS.cvl", AK8CHS_cvl, &b_AK8CHS_cvl);
   fChain->SetBranchAddress("AK8CHS.qgid", AK8CHS_qgid, &b_AK8CHS_qgid);
   fChain->SetBranchAddress("AK8CHS.axis2", AK8CHS_axis2, &b_AK8CHS_axis2);
   fChain->SetBranchAddress("AK8CHS.ptD", AK8CHS_ptD, &b_AK8CHS_ptD);
   fChain->SetBranchAddress("AK8CHS.mult", AK8CHS_mult, &b_AK8CHS_mult);
   fChain->SetBranchAddress("AK8CHS.q", AK8CHS_q, &b_AK8CHS_q);
   fChain->SetBranchAddress("AK8CHS.mva", AK8CHS_mva, &b_AK8CHS_mva);
   fChain->SetBranchAddress("AK8CHS.beta", AK8CHS_beta, &b_AK8CHS_beta);
   fChain->SetBranchAddress("AK8CHS.betaStar", AK8CHS_betaStar, &b_AK8CHS_betaStar);
   fChain->SetBranchAddress("AK8CHS.dR2Mean", AK8CHS_dR2Mean, &b_AK8CHS_dR2Mean);
   fChain->SetBranchAddress("AK8CHS.pullY", AK8CHS_pullY, &b_AK8CHS_pullY);
   fChain->SetBranchAddress("AK8CHS.pullPhi", AK8CHS_pullPhi, &b_AK8CHS_pullPhi);
   fChain->SetBranchAddress("AK8CHS.chPullY", AK8CHS_chPullY, &b_AK8CHS_chPullY);
   fChain->SetBranchAddress("AK8CHS.chPullPhi", AK8CHS_chPullPhi, &b_AK8CHS_chPullPhi);
   fChain->SetBranchAddress("AK8CHS.neuPullY", AK8CHS_neuPullY, &b_AK8CHS_neuPullY);
   fChain->SetBranchAddress("AK8CHS.neuPullPhi", AK8CHS_neuPullPhi, &b_AK8CHS_neuPullPhi);
   fChain->SetBranchAddress("AK8CHS.chEmFrac", AK8CHS_chEmFrac, &b_AK8CHS_chEmFrac);
   fChain->SetBranchAddress("AK8CHS.neuEmFrac", AK8CHS_neuEmFrac, &b_AK8CHS_neuEmFrac);
   fChain->SetBranchAddress("AK8CHS.chHadFrac", AK8CHS_chHadFrac, &b_AK8CHS_chHadFrac);
   fChain->SetBranchAddress("AK8CHS.neuHadFrac", AK8CHS_neuHadFrac, &b_AK8CHS_neuHadFrac);
   fChain->SetBranchAddress("AK8CHS.muonFrac", AK8CHS_muonFrac, &b_AK8CHS_muonFrac);
   fChain->SetBranchAddress("AK8CHS.genpt", AK8CHS_genpt, &b_AK8CHS_genpt);
   fChain->SetBranchAddress("AK8CHS.geneta", AK8CHS_geneta, &b_AK8CHS_geneta);
   fChain->SetBranchAddress("AK8CHS.genphi", AK8CHS_genphi, &b_AK8CHS_genphi);
   fChain->SetBranchAddress("AK8CHS.genm", AK8CHS_genm, &b_AK8CHS_genm);
   fChain->SetBranchAddress("AK8CHS.partonFlavor", AK8CHS_partonFlavor, &b_AK8CHS_partonFlavor);
   fChain->SetBranchAddress("AK8CHS.hadronFlavor", AK8CHS_hadronFlavor, &b_AK8CHS_hadronFlavor);
   fChain->SetBranchAddress("AK8CHS.nCharged", AK8CHS_nCharged, &b_AK8CHS_nCharged);
   fChain->SetBranchAddress("AK8CHS.nNeutrals", AK8CHS_nNeutrals, &b_AK8CHS_nNeutrals);
   fChain->SetBranchAddress("AK8CHS.nParticles", AK8CHS_nParticles, &b_AK8CHS_nParticles);
   fChain->SetBranchAddress("AddAK8CHS", &AddAK8CHS_, &b_AddAK8CHS_);
   fChain->SetBranchAddress("AddAK8CHS.index", AddAK8CHS_index, &b_AddAK8CHS_index);
   fChain->SetBranchAddress("AddAK8CHS.mass_prun", AddAK8CHS_mass_prun, &b_AddAK8CHS_mass_prun);
   fChain->SetBranchAddress("AddAK8CHS.mass_trim", AddAK8CHS_mass_trim, &b_AddAK8CHS_mass_trim);
   fChain->SetBranchAddress("AddAK8CHS.mass_sd0", AddAK8CHS_mass_sd0, &b_AddAK8CHS_mass_sd0);
   fChain->SetBranchAddress("AddAK8CHS.c2_0", AddAK8CHS_c2_0, &b_AddAK8CHS_c2_0);
   fChain->SetBranchAddress("AddAK8CHS.c2_0P2", AddAK8CHS_c2_0P2, &b_AddAK8CHS_c2_0P2);
   fChain->SetBranchAddress("AddAK8CHS.c2_0P5", AddAK8CHS_c2_0P5, &b_AddAK8CHS_c2_0P5);
   fChain->SetBranchAddress("AddAK8CHS.c2_1P0", AddAK8CHS_c2_1P0, &b_AddAK8CHS_c2_1P0);
   fChain->SetBranchAddress("AddAK8CHS.c2_2P0", AddAK8CHS_c2_2P0, &b_AddAK8CHS_c2_2P0);
   fChain->SetBranchAddress("AddAK8CHS.e2_b1", AddAK8CHS_e2_b1, &b_AddAK8CHS_e2_b1);
   fChain->SetBranchAddress("AddAK8CHS.e3_b1", AddAK8CHS_e3_b1, &b_AddAK8CHS_e3_b1);
   fChain->SetBranchAddress("AddAK8CHS.e3_v1_b1", AddAK8CHS_e3_v1_b1, &b_AddAK8CHS_e3_v1_b1);
   fChain->SetBranchAddress("AddAK8CHS.e3_v2_b1", AddAK8CHS_e3_v2_b1, &b_AddAK8CHS_e3_v2_b1);
   fChain->SetBranchAddress("AddAK8CHS.e4_v1_b1", AddAK8CHS_e4_v1_b1, &b_AddAK8CHS_e4_v1_b1);
   fChain->SetBranchAddress("AddAK8CHS.e4_v2_b1", AddAK8CHS_e4_v2_b1, &b_AddAK8CHS_e4_v2_b1);
   fChain->SetBranchAddress("AddAK8CHS.e2_b2", AddAK8CHS_e2_b2, &b_AddAK8CHS_e2_b2);
   fChain->SetBranchAddress("AddAK8CHS.e3_b2", AddAK8CHS_e3_b2, &b_AddAK8CHS_e3_b2);
   fChain->SetBranchAddress("AddAK8CHS.e3_v1_b2", AddAK8CHS_e3_v1_b2, &b_AddAK8CHS_e3_v1_b2);
   fChain->SetBranchAddress("AddAK8CHS.e3_v2_b2", AddAK8CHS_e3_v2_b2, &b_AddAK8CHS_e3_v2_b2);
   fChain->SetBranchAddress("AddAK8CHS.e4_v1_b2", AddAK8CHS_e4_v1_b2, &b_AddAK8CHS_e4_v1_b2);
   fChain->SetBranchAddress("AddAK8CHS.e4_v2_b2", AddAK8CHS_e4_v2_b2, &b_AddAK8CHS_e4_v2_b2);
   fChain->SetBranchAddress("AddAK8CHS.e2_sdb1", AddAK8CHS_e2_sdb1, &b_AddAK8CHS_e2_sdb1);
   fChain->SetBranchAddress("AddAK8CHS.e3_sdb1", AddAK8CHS_e3_sdb1, &b_AddAK8CHS_e3_sdb1);
   fChain->SetBranchAddress("AddAK8CHS.e3_v1_sdb1", AddAK8CHS_e3_v1_sdb1, &b_AddAK8CHS_e3_v1_sdb1);
   fChain->SetBranchAddress("AddAK8CHS.e3_v2_sdb1", AddAK8CHS_e3_v2_sdb1, &b_AddAK8CHS_e3_v2_sdb1);
   fChain->SetBranchAddress("AddAK8CHS.e4_v1_sdb1", AddAK8CHS_e4_v1_sdb1, &b_AddAK8CHS_e4_v1_sdb1);
   fChain->SetBranchAddress("AddAK8CHS.e4_v2_sdb1", AddAK8CHS_e4_v2_sdb1, &b_AddAK8CHS_e4_v2_sdb1);
   fChain->SetBranchAddress("AddAK8CHS.e2_sdb2", AddAK8CHS_e2_sdb2, &b_AddAK8CHS_e2_sdb2);
   fChain->SetBranchAddress("AddAK8CHS.e3_sdb2", AddAK8CHS_e3_sdb2, &b_AddAK8CHS_e3_sdb2);
   fChain->SetBranchAddress("AddAK8CHS.e3_v1_sdb2", AddAK8CHS_e3_v1_sdb2, &b_AddAK8CHS_e3_v1_sdb2);
   fChain->SetBranchAddress("AddAK8CHS.e3_v2_sdb2", AddAK8CHS_e3_v2_sdb2, &b_AddAK8CHS_e3_v2_sdb2);
   fChain->SetBranchAddress("AddAK8CHS.e4_v1_sdb2", AddAK8CHS_e4_v1_sdb2, &b_AddAK8CHS_e4_v1_sdb2);
   fChain->SetBranchAddress("AddAK8CHS.e4_v2_sdb2", AddAK8CHS_e4_v2_sdb2, &b_AddAK8CHS_e4_v2_sdb2);
   fChain->SetBranchAddress("AddAK8CHS.qjet", AddAK8CHS_qjet, &b_AddAK8CHS_qjet);
   fChain->SetBranchAddress("AddAK8CHS.tau1", AddAK8CHS_tau1, &b_AddAK8CHS_tau1);
   fChain->SetBranchAddress("AddAK8CHS.tau2", AddAK8CHS_tau2, &b_AddAK8CHS_tau2);
   fChain->SetBranchAddress("AddAK8CHS.tau3", AddAK8CHS_tau3, &b_AddAK8CHS_tau3);
   fChain->SetBranchAddress("AddAK8CHS.tau4", AddAK8CHS_tau4, &b_AddAK8CHS_tau4);
   fChain->SetBranchAddress("AddAK8CHS.doublecsv", AddAK8CHS_doublecsv, &b_AddAK8CHS_doublecsv);
   fChain->SetBranchAddress("AddAK8CHS.Double_sub", AddAK8CHS_Double_sub, &b_AddAK8CHS_Double_sub);
   fChain->SetBranchAddress("AddAK8CHS.sj1_pt", AddAK8CHS_sj1_pt, &b_AddAK8CHS_sj1_pt);
   fChain->SetBranchAddress("AddAK8CHS.sj1_eta", AddAK8CHS_sj1_eta, &b_AddAK8CHS_sj1_eta);
   fChain->SetBranchAddress("AddAK8CHS.sj1_phi", AddAK8CHS_sj1_phi, &b_AddAK8CHS_sj1_phi);
   fChain->SetBranchAddress("AddAK8CHS.sj1_m", AddAK8CHS_sj1_m, &b_AddAK8CHS_sj1_m);
   fChain->SetBranchAddress("AddAK8CHS.sj1_csv", AddAK8CHS_sj1_csv, &b_AddAK8CHS_sj1_csv);
   fChain->SetBranchAddress("AddAK8CHS.sj1_qgid", AddAK8CHS_sj1_qgid, &b_AddAK8CHS_sj1_qgid);
   fChain->SetBranchAddress("AddAK8CHS.sj1_q", AddAK8CHS_sj1_q, &b_AddAK8CHS_sj1_q);
   fChain->SetBranchAddress("AddAK8CHS.sj2_pt", AddAK8CHS_sj2_pt, &b_AddAK8CHS_sj2_pt);
   fChain->SetBranchAddress("AddAK8CHS.sj2_eta", AddAK8CHS_sj2_eta, &b_AddAK8CHS_sj2_eta);
   fChain->SetBranchAddress("AddAK8CHS.sj2_phi", AddAK8CHS_sj2_phi, &b_AddAK8CHS_sj2_phi);
   fChain->SetBranchAddress("AddAK8CHS.sj2_m", AddAK8CHS_sj2_m, &b_AddAK8CHS_sj2_m);
   fChain->SetBranchAddress("AddAK8CHS.sj2_csv", AddAK8CHS_sj2_csv, &b_AddAK8CHS_sj2_csv);
   fChain->SetBranchAddress("AddAK8CHS.sj2_qgid", AddAK8CHS_sj2_qgid, &b_AddAK8CHS_sj2_qgid);
   fChain->SetBranchAddress("AddAK8CHS.sj2_q", AddAK8CHS_sj2_q, &b_AddAK8CHS_sj2_q);
   fChain->SetBranchAddress("AddAK8CHS.sj3_pt", AddAK8CHS_sj3_pt, &b_AddAK8CHS_sj3_pt);
   fChain->SetBranchAddress("AddAK8CHS.sj3_eta", AddAK8CHS_sj3_eta, &b_AddAK8CHS_sj3_eta);
   fChain->SetBranchAddress("AddAK8CHS.sj3_phi", AddAK8CHS_sj3_phi, &b_AddAK8CHS_sj3_phi);
   fChain->SetBranchAddress("AddAK8CHS.sj3_m", AddAK8CHS_sj3_m, &b_AddAK8CHS_sj3_m);
   fChain->SetBranchAddress("AddAK8CHS.sj3_csv", AddAK8CHS_sj3_csv, &b_AddAK8CHS_sj3_csv);
   fChain->SetBranchAddress("AddAK8CHS.sj3_qgid", AddAK8CHS_sj3_qgid, &b_AddAK8CHS_sj3_qgid);
   fChain->SetBranchAddress("AddAK8CHS.sj3_q", AddAK8CHS_sj3_q, &b_AddAK8CHS_sj3_q);
   fChain->SetBranchAddress("AddAK8CHS.sj4_pt", AddAK8CHS_sj4_pt, &b_AddAK8CHS_sj4_pt);
   fChain->SetBranchAddress("AddAK8CHS.sj4_eta", AddAK8CHS_sj4_eta, &b_AddAK8CHS_sj4_eta);
   fChain->SetBranchAddress("AddAK8CHS.sj4_phi", AddAK8CHS_sj4_phi, &b_AddAK8CHS_sj4_phi);
   fChain->SetBranchAddress("AddAK8CHS.sj4_m", AddAK8CHS_sj4_m, &b_AddAK8CHS_sj4_m);
   fChain->SetBranchAddress("AddAK8CHS.sj4_csv", AddAK8CHS_sj4_csv, &b_AddAK8CHS_sj4_csv);
   fChain->SetBranchAddress("AddAK8CHS.sj4_qgid", AddAK8CHS_sj4_qgid, &b_AddAK8CHS_sj4_qgid);
   fChain->SetBranchAddress("AddAK8CHS.sj4_q", AddAK8CHS_sj4_q, &b_AddAK8CHS_sj4_q);
   fChain->SetBranchAddress("AddAK8CHS.pullAngle", AddAK8CHS_pullAngle, &b_AddAK8CHS_pullAngle);
   fChain->SetBranchAddress("AddAK8CHS.topTagType", AddAK8CHS_topTagType, &b_AddAK8CHS_topTagType);
   fChain->SetBranchAddress("AddAK8CHS.top_n_subjets", AddAK8CHS_top_n_subjets, &b_AddAK8CHS_top_n_subjets);
   fChain->SetBranchAddress("AddAK8CHS.top_m_min", AddAK8CHS_top_m_min, &b_AddAK8CHS_top_m_min);
   fChain->SetBranchAddress("AddAK8CHS.top_m_123", AddAK8CHS_top_m_123, &b_AddAK8CHS_top_m_123);
   fChain->SetBranchAddress("AddAK8CHS.top_fRec", AddAK8CHS_top_fRec, &b_AddAK8CHS_top_fRec);
   fChain->SetBranchAddress("AddAK8CHS.topchi2", AddAK8CHS_topchi2, &b_AddAK8CHS_topchi2);
   fChain->SetBranchAddress("AK4Puppi", &AK4Puppi_, &b_AK4Puppi_);
   fChain->SetBranchAddress("AK4Puppi.pt", AK4Puppi_pt, &b_AK4Puppi_pt);
   fChain->SetBranchAddress("AK4Puppi.eta", AK4Puppi_eta, &b_AK4Puppi_eta);
   fChain->SetBranchAddress("AK4Puppi.phi", AK4Puppi_phi, &b_AK4Puppi_phi);
   fChain->SetBranchAddress("AK4Puppi.mass", AK4Puppi_mass, &b_AK4Puppi_mass);
   fChain->SetBranchAddress("AK4Puppi.ptRaw", AK4Puppi_ptRaw, &b_AK4Puppi_ptRaw);
   fChain->SetBranchAddress("AK4Puppi.unc", AK4Puppi_unc, &b_AK4Puppi_unc);
   fChain->SetBranchAddress("AK4Puppi.area", AK4Puppi_area, &b_AK4Puppi_area);
   fChain->SetBranchAddress("AK4Puppi.d0", AK4Puppi_d0, &b_AK4Puppi_d0);
   fChain->SetBranchAddress("AK4Puppi.dz", AK4Puppi_dz, &b_AK4Puppi_dz);
   fChain->SetBranchAddress("AK4Puppi.csv", AK4Puppi_csv, &b_AK4Puppi_csv);
   fChain->SetBranchAddress("AK4Puppi.bmva", AK4Puppi_bmva, &b_AK4Puppi_bmva);
   fChain->SetBranchAddress("AK4Puppi.cvb", AK4Puppi_cvb, &b_AK4Puppi_cvb);
   fChain->SetBranchAddress("AK4Puppi.cvl", AK4Puppi_cvl, &b_AK4Puppi_cvl);
   fChain->SetBranchAddress("AK4Puppi.qgid", AK4Puppi_qgid, &b_AK4Puppi_qgid);
   fChain->SetBranchAddress("AK4Puppi.axis2", AK4Puppi_axis2, &b_AK4Puppi_axis2);
   fChain->SetBranchAddress("AK4Puppi.ptD", AK4Puppi_ptD, &b_AK4Puppi_ptD);
   fChain->SetBranchAddress("AK4Puppi.mult", AK4Puppi_mult, &b_AK4Puppi_mult);
   fChain->SetBranchAddress("AK4Puppi.q", AK4Puppi_q, &b_AK4Puppi_q);
   fChain->SetBranchAddress("AK4Puppi.mva", AK4Puppi_mva, &b_AK4Puppi_mva);
   fChain->SetBranchAddress("AK4Puppi.beta", AK4Puppi_beta, &b_AK4Puppi_beta);
   fChain->SetBranchAddress("AK4Puppi.betaStar", AK4Puppi_betaStar, &b_AK4Puppi_betaStar);
   fChain->SetBranchAddress("AK4Puppi.dR2Mean", AK4Puppi_dR2Mean, &b_AK4Puppi_dR2Mean);
   fChain->SetBranchAddress("AK4Puppi.pullY", AK4Puppi_pullY, &b_AK4Puppi_pullY);
   fChain->SetBranchAddress("AK4Puppi.pullPhi", AK4Puppi_pullPhi, &b_AK4Puppi_pullPhi);
   fChain->SetBranchAddress("AK4Puppi.chPullY", AK4Puppi_chPullY, &b_AK4Puppi_chPullY);
   fChain->SetBranchAddress("AK4Puppi.chPullPhi", AK4Puppi_chPullPhi, &b_AK4Puppi_chPullPhi);
   fChain->SetBranchAddress("AK4Puppi.neuPullY", AK4Puppi_neuPullY, &b_AK4Puppi_neuPullY);
   fChain->SetBranchAddress("AK4Puppi.neuPullPhi", AK4Puppi_neuPullPhi, &b_AK4Puppi_neuPullPhi);
   fChain->SetBranchAddress("AK4Puppi.chEmFrac", AK4Puppi_chEmFrac, &b_AK4Puppi_chEmFrac);
   fChain->SetBranchAddress("AK4Puppi.neuEmFrac", AK4Puppi_neuEmFrac, &b_AK4Puppi_neuEmFrac);
   fChain->SetBranchAddress("AK4Puppi.chHadFrac", AK4Puppi_chHadFrac, &b_AK4Puppi_chHadFrac);
   fChain->SetBranchAddress("AK4Puppi.neuHadFrac", AK4Puppi_neuHadFrac, &b_AK4Puppi_neuHadFrac);
   fChain->SetBranchAddress("AK4Puppi.muonFrac", AK4Puppi_muonFrac, &b_AK4Puppi_muonFrac);
   fChain->SetBranchAddress("AK4Puppi.genpt", AK4Puppi_genpt, &b_AK4Puppi_genpt);
   fChain->SetBranchAddress("AK4Puppi.geneta", AK4Puppi_geneta, &b_AK4Puppi_geneta);
   fChain->SetBranchAddress("AK4Puppi.genphi", AK4Puppi_genphi, &b_AK4Puppi_genphi);
   fChain->SetBranchAddress("AK4Puppi.genm", AK4Puppi_genm, &b_AK4Puppi_genm);
   fChain->SetBranchAddress("AK4Puppi.partonFlavor", AK4Puppi_partonFlavor, &b_AK4Puppi_partonFlavor);
   fChain->SetBranchAddress("AK4Puppi.hadronFlavor", AK4Puppi_hadronFlavor, &b_AK4Puppi_hadronFlavor);
   fChain->SetBranchAddress("AK4Puppi.nCharged", AK4Puppi_nCharged, &b_AK4Puppi_nCharged);
   fChain->SetBranchAddress("AK4Puppi.nNeutrals", AK4Puppi_nNeutrals, &b_AK4Puppi_nNeutrals);
   fChain->SetBranchAddress("AK4Puppi.nParticles", AK4Puppi_nParticles, &b_AK4Puppi_nParticles);
   fChain->SetBranchAddress("AK8Puppi", &AK8Puppi_, &b_AK8Puppi_);
   fChain->SetBranchAddress("AK8Puppi.pt", AK8Puppi_pt, &b_AK8Puppi_pt);
   fChain->SetBranchAddress("AK8Puppi.eta", AK8Puppi_eta, &b_AK8Puppi_eta);
   fChain->SetBranchAddress("AK8Puppi.phi", AK8Puppi_phi, &b_AK8Puppi_phi);
   fChain->SetBranchAddress("AK8Puppi.mass", AK8Puppi_mass, &b_AK8Puppi_mass);
   fChain->SetBranchAddress("AK8Puppi.ptRaw", AK8Puppi_ptRaw, &b_AK8Puppi_ptRaw);
   fChain->SetBranchAddress("AK8Puppi.unc", AK8Puppi_unc, &b_AK8Puppi_unc);
   fChain->SetBranchAddress("AK8Puppi.area", AK8Puppi_area, &b_AK8Puppi_area);
   fChain->SetBranchAddress("AK8Puppi.d0", AK8Puppi_d0, &b_AK8Puppi_d0);
   fChain->SetBranchAddress("AK8Puppi.dz", AK8Puppi_dz, &b_AK8Puppi_dz);
   fChain->SetBranchAddress("AK8Puppi.csv", AK8Puppi_csv, &b_AK8Puppi_csv);
   fChain->SetBranchAddress("AK8Puppi.bmva", AK8Puppi_bmva, &b_AK8Puppi_bmva);
   fChain->SetBranchAddress("AK8Puppi.cvb", AK8Puppi_cvb, &b_AK8Puppi_cvb);
   fChain->SetBranchAddress("AK8Puppi.cvl", AK8Puppi_cvl, &b_AK8Puppi_cvl);
   fChain->SetBranchAddress("AK8Puppi.qgid", AK8Puppi_qgid, &b_AK8Puppi_qgid);
   fChain->SetBranchAddress("AK8Puppi.axis2", AK8Puppi_axis2, &b_AK8Puppi_axis2);
   fChain->SetBranchAddress("AK8Puppi.ptD", AK8Puppi_ptD, &b_AK8Puppi_ptD);
   fChain->SetBranchAddress("AK8Puppi.mult", AK8Puppi_mult, &b_AK8Puppi_mult);
   fChain->SetBranchAddress("AK8Puppi.q", AK8Puppi_q, &b_AK8Puppi_q);
   fChain->SetBranchAddress("AK8Puppi.mva", AK8Puppi_mva, &b_AK8Puppi_mva);
   fChain->SetBranchAddress("AK8Puppi.beta", AK8Puppi_beta, &b_AK8Puppi_beta);
   fChain->SetBranchAddress("AK8Puppi.betaStar", AK8Puppi_betaStar, &b_AK8Puppi_betaStar);
   fChain->SetBranchAddress("AK8Puppi.dR2Mean", AK8Puppi_dR2Mean, &b_AK8Puppi_dR2Mean);
   fChain->SetBranchAddress("AK8Puppi.pullY", AK8Puppi_pullY, &b_AK8Puppi_pullY);
   fChain->SetBranchAddress("AK8Puppi.pullPhi", AK8Puppi_pullPhi, &b_AK8Puppi_pullPhi);
   fChain->SetBranchAddress("AK8Puppi.chPullY", AK8Puppi_chPullY, &b_AK8Puppi_chPullY);
   fChain->SetBranchAddress("AK8Puppi.chPullPhi", AK8Puppi_chPullPhi, &b_AK8Puppi_chPullPhi);
   fChain->SetBranchAddress("AK8Puppi.neuPullY", AK8Puppi_neuPullY, &b_AK8Puppi_neuPullY);
   fChain->SetBranchAddress("AK8Puppi.neuPullPhi", AK8Puppi_neuPullPhi, &b_AK8Puppi_neuPullPhi);
   fChain->SetBranchAddress("AK8Puppi.chEmFrac", AK8Puppi_chEmFrac, &b_AK8Puppi_chEmFrac);
   fChain->SetBranchAddress("AK8Puppi.neuEmFrac", AK8Puppi_neuEmFrac, &b_AK8Puppi_neuEmFrac);
   fChain->SetBranchAddress("AK8Puppi.chHadFrac", AK8Puppi_chHadFrac, &b_AK8Puppi_chHadFrac);
   fChain->SetBranchAddress("AK8Puppi.neuHadFrac", AK8Puppi_neuHadFrac, &b_AK8Puppi_neuHadFrac);
   fChain->SetBranchAddress("AK8Puppi.muonFrac", AK8Puppi_muonFrac, &b_AK8Puppi_muonFrac);
   fChain->SetBranchAddress("AK8Puppi.genpt", AK8Puppi_genpt, &b_AK8Puppi_genpt);
   fChain->SetBranchAddress("AK8Puppi.geneta", AK8Puppi_geneta, &b_AK8Puppi_geneta);
   fChain->SetBranchAddress("AK8Puppi.genphi", AK8Puppi_genphi, &b_AK8Puppi_genphi);
   fChain->SetBranchAddress("AK8Puppi.genm", AK8Puppi_genm, &b_AK8Puppi_genm);
   fChain->SetBranchAddress("AK8Puppi.partonFlavor", AK8Puppi_partonFlavor, &b_AK8Puppi_partonFlavor);
   fChain->SetBranchAddress("AK8Puppi.hadronFlavor", AK8Puppi_hadronFlavor, &b_AK8Puppi_hadronFlavor);
   fChain->SetBranchAddress("AK8Puppi.nCharged", AK8Puppi_nCharged, &b_AK8Puppi_nCharged);
   fChain->SetBranchAddress("AK8Puppi.nNeutrals", AK8Puppi_nNeutrals, &b_AK8Puppi_nNeutrals);
   fChain->SetBranchAddress("AK8Puppi.nParticles", AK8Puppi_nParticles, &b_AK8Puppi_nParticles);
   fChain->SetBranchAddress("AddAK8Puppi", &AddAK8Puppi_, &b_AddAK8Puppi_);
   fChain->SetBranchAddress("AddAK8Puppi.index", AddAK8Puppi_index, &b_AddAK8Puppi_index);
   fChain->SetBranchAddress("AddAK8Puppi.mass_prun", AddAK8Puppi_mass_prun, &b_AddAK8Puppi_mass_prun);
   fChain->SetBranchAddress("AddAK8Puppi.mass_trim", AddAK8Puppi_mass_trim, &b_AddAK8Puppi_mass_trim);
   fChain->SetBranchAddress("AddAK8Puppi.mass_sd0", AddAK8Puppi_mass_sd0, &b_AddAK8Puppi_mass_sd0);
   fChain->SetBranchAddress("AddAK8Puppi.c2_0", AddAK8Puppi_c2_0, &b_AddAK8Puppi_c2_0);
   fChain->SetBranchAddress("AddAK8Puppi.c2_0P2", AddAK8Puppi_c2_0P2, &b_AddAK8Puppi_c2_0P2);
   fChain->SetBranchAddress("AddAK8Puppi.c2_0P5", AddAK8Puppi_c2_0P5, &b_AddAK8Puppi_c2_0P5);
   fChain->SetBranchAddress("AddAK8Puppi.c2_1P0", AddAK8Puppi_c2_1P0, &b_AddAK8Puppi_c2_1P0);
   fChain->SetBranchAddress("AddAK8Puppi.c2_2P0", AddAK8Puppi_c2_2P0, &b_AddAK8Puppi_c2_2P0);
   fChain->SetBranchAddress("AddAK8Puppi.e2_b1", AddAK8Puppi_e2_b1, &b_AddAK8Puppi_e2_b1);
   fChain->SetBranchAddress("AddAK8Puppi.e3_b1", AddAK8Puppi_e3_b1, &b_AddAK8Puppi_e3_b1);
   fChain->SetBranchAddress("AddAK8Puppi.e3_v1_b1", AddAK8Puppi_e3_v1_b1, &b_AddAK8Puppi_e3_v1_b1);
   fChain->SetBranchAddress("AddAK8Puppi.e3_v2_b1", AddAK8Puppi_e3_v2_b1, &b_AddAK8Puppi_e3_v2_b1);
   fChain->SetBranchAddress("AddAK8Puppi.e4_v1_b1", AddAK8Puppi_e4_v1_b1, &b_AddAK8Puppi_e4_v1_b1);
   fChain->SetBranchAddress("AddAK8Puppi.e4_v2_b1", AddAK8Puppi_e4_v2_b1, &b_AddAK8Puppi_e4_v2_b1);
   fChain->SetBranchAddress("AddAK8Puppi.e2_b2", AddAK8Puppi_e2_b2, &b_AddAK8Puppi_e2_b2);
   fChain->SetBranchAddress("AddAK8Puppi.e3_b2", AddAK8Puppi_e3_b2, &b_AddAK8Puppi_e3_b2);
   fChain->SetBranchAddress("AddAK8Puppi.e3_v1_b2", AddAK8Puppi_e3_v1_b2, &b_AddAK8Puppi_e3_v1_b2);
   fChain->SetBranchAddress("AddAK8Puppi.e3_v2_b2", AddAK8Puppi_e3_v2_b2, &b_AddAK8Puppi_e3_v2_b2);
   fChain->SetBranchAddress("AddAK8Puppi.e4_v1_b2", AddAK8Puppi_e4_v1_b2, &b_AddAK8Puppi_e4_v1_b2);
   fChain->SetBranchAddress("AddAK8Puppi.e4_v2_b2", AddAK8Puppi_e4_v2_b2, &b_AddAK8Puppi_e4_v2_b2);
   fChain->SetBranchAddress("AddAK8Puppi.e2_sdb1", AddAK8Puppi_e2_sdb1, &b_AddAK8Puppi_e2_sdb1);
   fChain->SetBranchAddress("AddAK8Puppi.e3_sdb1", AddAK8Puppi_e3_sdb1, &b_AddAK8Puppi_e3_sdb1);
   fChain->SetBranchAddress("AddAK8Puppi.e3_v1_sdb1", AddAK8Puppi_e3_v1_sdb1, &b_AddAK8Puppi_e3_v1_sdb1);
   fChain->SetBranchAddress("AddAK8Puppi.e3_v2_sdb1", AddAK8Puppi_e3_v2_sdb1, &b_AddAK8Puppi_e3_v2_sdb1);
   fChain->SetBranchAddress("AddAK8Puppi.e4_v1_sdb1", AddAK8Puppi_e4_v1_sdb1, &b_AddAK8Puppi_e4_v1_sdb1);
   fChain->SetBranchAddress("AddAK8Puppi.e4_v2_sdb1", AddAK8Puppi_e4_v2_sdb1, &b_AddAK8Puppi_e4_v2_sdb1);
   fChain->SetBranchAddress("AddAK8Puppi.e2_sdb2", AddAK8Puppi_e2_sdb2, &b_AddAK8Puppi_e2_sdb2);
   fChain->SetBranchAddress("AddAK8Puppi.e3_sdb2", AddAK8Puppi_e3_sdb2, &b_AddAK8Puppi_e3_sdb2);
   fChain->SetBranchAddress("AddAK8Puppi.e3_v1_sdb2", AddAK8Puppi_e3_v1_sdb2, &b_AddAK8Puppi_e3_v1_sdb2);
   fChain->SetBranchAddress("AddAK8Puppi.e3_v2_sdb2", AddAK8Puppi_e3_v2_sdb2, &b_AddAK8Puppi_e3_v2_sdb2);
   fChain->SetBranchAddress("AddAK8Puppi.e4_v1_sdb2", AddAK8Puppi_e4_v1_sdb2, &b_AddAK8Puppi_e4_v1_sdb2);
   fChain->SetBranchAddress("AddAK8Puppi.e4_v2_sdb2", AddAK8Puppi_e4_v2_sdb2, &b_AddAK8Puppi_e4_v2_sdb2);
   fChain->SetBranchAddress("AddAK8Puppi.qjet", AddAK8Puppi_qjet, &b_AddAK8Puppi_qjet);
   fChain->SetBranchAddress("AddAK8Puppi.tau1", AddAK8Puppi_tau1, &b_AddAK8Puppi_tau1);
   fChain->SetBranchAddress("AddAK8Puppi.tau2", AddAK8Puppi_tau2, &b_AddAK8Puppi_tau2);
   fChain->SetBranchAddress("AddAK8Puppi.tau3", AddAK8Puppi_tau3, &b_AddAK8Puppi_tau3);
   fChain->SetBranchAddress("AddAK8Puppi.tau4", AddAK8Puppi_tau4, &b_AddAK8Puppi_tau4);
   fChain->SetBranchAddress("AddAK8Puppi.doublecsv", AddAK8Puppi_doublecsv, &b_AddAK8Puppi_doublecsv);
   fChain->SetBranchAddress("AddAK8Puppi.Double_sub", AddAK8Puppi_Double_sub, &b_AddAK8Puppi_Double_sub);
   fChain->SetBranchAddress("AddAK8Puppi.sj1_pt", AddAK8Puppi_sj1_pt, &b_AddAK8Puppi_sj1_pt);
   fChain->SetBranchAddress("AddAK8Puppi.sj1_eta", AddAK8Puppi_sj1_eta, &b_AddAK8Puppi_sj1_eta);
   fChain->SetBranchAddress("AddAK8Puppi.sj1_phi", AddAK8Puppi_sj1_phi, &b_AddAK8Puppi_sj1_phi);
   fChain->SetBranchAddress("AddAK8Puppi.sj1_m", AddAK8Puppi_sj1_m, &b_AddAK8Puppi_sj1_m);
   fChain->SetBranchAddress("AddAK8Puppi.sj1_csv", AddAK8Puppi_sj1_csv, &b_AddAK8Puppi_sj1_csv);
   fChain->SetBranchAddress("AddAK8Puppi.sj1_qgid", AddAK8Puppi_sj1_qgid, &b_AddAK8Puppi_sj1_qgid);
   fChain->SetBranchAddress("AddAK8Puppi.sj1_q", AddAK8Puppi_sj1_q, &b_AddAK8Puppi_sj1_q);
   fChain->SetBranchAddress("AddAK8Puppi.sj2_pt", AddAK8Puppi_sj2_pt, &b_AddAK8Puppi_sj2_pt);
   fChain->SetBranchAddress("AddAK8Puppi.sj2_eta", AddAK8Puppi_sj2_eta, &b_AddAK8Puppi_sj2_eta);
   fChain->SetBranchAddress("AddAK8Puppi.sj2_phi", AddAK8Puppi_sj2_phi, &b_AddAK8Puppi_sj2_phi);
   fChain->SetBranchAddress("AddAK8Puppi.sj2_m", AddAK8Puppi_sj2_m, &b_AddAK8Puppi_sj2_m);
   fChain->SetBranchAddress("AddAK8Puppi.sj2_csv", AddAK8Puppi_sj2_csv, &b_AddAK8Puppi_sj2_csv);
   fChain->SetBranchAddress("AddAK8Puppi.sj2_qgid", AddAK8Puppi_sj2_qgid, &b_AddAK8Puppi_sj2_qgid);
   fChain->SetBranchAddress("AddAK8Puppi.sj2_q", AddAK8Puppi_sj2_q, &b_AddAK8Puppi_sj2_q);
   fChain->SetBranchAddress("AddAK8Puppi.sj3_pt", AddAK8Puppi_sj3_pt, &b_AddAK8Puppi_sj3_pt);
   fChain->SetBranchAddress("AddAK8Puppi.sj3_eta", AddAK8Puppi_sj3_eta, &b_AddAK8Puppi_sj3_eta);
   fChain->SetBranchAddress("AddAK8Puppi.sj3_phi", AddAK8Puppi_sj3_phi, &b_AddAK8Puppi_sj3_phi);
   fChain->SetBranchAddress("AddAK8Puppi.sj3_m", AddAK8Puppi_sj3_m, &b_AddAK8Puppi_sj3_m);
   fChain->SetBranchAddress("AddAK8Puppi.sj3_csv", AddAK8Puppi_sj3_csv, &b_AddAK8Puppi_sj3_csv);
   fChain->SetBranchAddress("AddAK8Puppi.sj3_qgid", AddAK8Puppi_sj3_qgid, &b_AddAK8Puppi_sj3_qgid);
   fChain->SetBranchAddress("AddAK8Puppi.sj3_q", AddAK8Puppi_sj3_q, &b_AddAK8Puppi_sj3_q);
   fChain->SetBranchAddress("AddAK8Puppi.sj4_pt", AddAK8Puppi_sj4_pt, &b_AddAK8Puppi_sj4_pt);
   fChain->SetBranchAddress("AddAK8Puppi.sj4_eta", AddAK8Puppi_sj4_eta, &b_AddAK8Puppi_sj4_eta);
   fChain->SetBranchAddress("AddAK8Puppi.sj4_phi", AddAK8Puppi_sj4_phi, &b_AddAK8Puppi_sj4_phi);
   fChain->SetBranchAddress("AddAK8Puppi.sj4_m", AddAK8Puppi_sj4_m, &b_AddAK8Puppi_sj4_m);
   fChain->SetBranchAddress("AddAK8Puppi.sj4_csv", AddAK8Puppi_sj4_csv, &b_AddAK8Puppi_sj4_csv);
   fChain->SetBranchAddress("AddAK8Puppi.sj4_qgid", AddAK8Puppi_sj4_qgid, &b_AddAK8Puppi_sj4_qgid);
   fChain->SetBranchAddress("AddAK8Puppi.sj4_q", AddAK8Puppi_sj4_q, &b_AddAK8Puppi_sj4_q);
   fChain->SetBranchAddress("AddAK8Puppi.pullAngle", AddAK8Puppi_pullAngle, &b_AddAK8Puppi_pullAngle);
   fChain->SetBranchAddress("AddAK8Puppi.topTagType", AddAK8Puppi_topTagType, &b_AddAK8Puppi_topTagType);
   fChain->SetBranchAddress("AddAK8Puppi.top_n_subjets", AddAK8Puppi_top_n_subjets, &b_AddAK8Puppi_top_n_subjets);
   fChain->SetBranchAddress("AddAK8Puppi.top_m_min", AddAK8Puppi_top_m_min, &b_AddAK8Puppi_top_m_min);
   fChain->SetBranchAddress("AddAK8Puppi.top_m_123", AddAK8Puppi_top_m_123, &b_AddAK8Puppi_top_m_123);
   fChain->SetBranchAddress("AddAK8Puppi.top_fRec", AddAK8Puppi_top_fRec, &b_AddAK8Puppi_top_fRec);
   fChain->SetBranchAddress("AddAK8Puppi.topchi2", AddAK8Puppi_topchi2, &b_AddAK8Puppi_topchi2);
   Notify();
}

Bool_t GetNegEvents::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void GetNegEvents::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t GetNegEvents::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef GetNegEvents_cxx
