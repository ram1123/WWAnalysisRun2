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
#include <TClonesArray.h>           

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

// HEADER FILE FOR JES 
#include "CondFormats/JetMETObjects/interface/JetCorrectorParameters.h"
#include "CondFormats/JetMETObjects/interface/FactorizedJetCorrector.h"
#include "CondFormats/JetMETObjects/interface/JetCorrectionUncertainty.h"

// HEADER FOR B-TAG
#include "CondFormats/BTauObjects/interface/BTagCalibration.h"
#include "CondTools/BTau/interface/BTagCalibrationReader.h"
#include "../BtagUnc.hh"

#include "../interface/setOutputTree.h"
#include "../interface/METzCalculator.h"
#include "../interface/METzCalculator_Run2.h"
#include "../interface/analysisUtils.h"
#include "../interface/readJSONFile.h"
#include "../interface/Utils.hh"
#include "../interface/GeneralizedEndpoint.h"

using namespace std;

//*******MAIN*******************************************************************

int main (int argc, char** argv)
{ 
  
  int t0 = time(NULL);
  
  std::string inputFolder = argv[1];
  std::string outputFile = argv[2];
  int isMC = atoi(argv[3]);
  std::string cluster = argv[4];
  std::string inputTreeName = argv[5];
  std::string inputFile = argv[6];
  std::string xSecWeight = argv[7];
  std::string TotalNumberOfEntries = argv[8];
  std::string TotalNumberOfNegativeEntries = argv[9];
  int applyTrigger = atoi(argv[10]);
  std::string jsonFileName = argv[11];
  int isLocal = atoi(argv[12]);
  int VBFSel  = atoi(argv[13]);
  
  std::string leptonName;

  if ( VBFSel==1)	cout<<"==> VBF selection method : Select two highest pT jets"<<endl;
  else if ( VBFSel==2)	cout<<"==> VBF selection method : Select pair with highest mjj..."<<endl;
  else if ( VBFSel==3)	cout<<"==> VBF selection method : Select pair with highest DeltaEta..."<<endl;
  else {	cout<<"\n\nERROR:	Enter valid vbf selection criteria....\n\n"<<endl;
  		exit(0);  
	}
  
  std::string iHLTFile="${CMSSW_BASE}/src/BaconAna/DataFormats/data/HLTFile_25ns";
  const std::string cmssw_base = getenv("CMSSW_BASE");
  std::string cmssw_base_env = "${CMSSW_BASE}";
  size_t start_pos = iHLTFile.find(cmssw_base_env);
  if(start_pos != std::string::npos) {
  	iHLTFile.replace(start_pos, cmssw_base_env.length(), cmssw_base);
  }

  const baconhep::TTrigger triggerMenu(iHLTFile);  
  std::cout<<"Apply trigger: "<<applyTrigger<<std::endl;

  TLorentzVector W_type0,W_type0_jes_up, W_type0_jes_dn, W_type0_jer_up, W_type0_jer_dn, W_type2, W_run2,W_puppi_type2, W_puppi_type0, W_puppi_run2, W_type0_LEP_Up, W_type0_LEP_Down;
  TLorentzVector LEP1, LEP2, LEP1_Up, LEP1_Down, LEP2_Up, LEP2_Down, SJ1_PuppiAK8, SJ2_PuppiAK8, SJ1, SJ2;
  TLorentzVector NU0,NU1,NU2,NU0_puppi,NU1_puppi,NU2_puppi;
  TLorentzVector NU0_jes_up, NU0_jes_dn;
  TLorentzVector NU0_puppi_jes_up, NU0_puppi_jes_dn;
  TLorentzVector NU0_jer_up, NU0_jer_dn;
  TLorentzVector JET, JET_PuppiAK8, AK4;
  TLorentzVector JET_jes_up, JET_jes_dn, JET_PuppiAK8_jes_up, JET_PuppiAK8_jes_dn;
  TLorentzVector AK4_JET1,AK4_JET2;
  TLorentzVector AK4_JET1_jes_up, AK4_JET1_jes_dn;
  TLorentzVector AK4_JET2_jes_up, AK4_JET2_jes_dn;
  TLorentzVector PuppiAK4_JET1,PuppiAK4_JET2;
  TLorentzVector PuppiAK4_JET1_jes_up, PuppiAK4_JET1_jes_dn;
  TLorentzVector PuppiAK4_JET2_jes_up, PuppiAK4_JET2_jes_dn;
  TLorentzVector VBF1,VBF2,TOT;
  TLorentzVector VBF1_jes_up, VBF1_jes_dn, VBF2_jes_up, VBF2_jes_dn;
  TLorentzVector ELE,MU;

  
  std::vector<TLorentzVector> tightMuon;
  std::vector<TLorentzVector> looseMuon;
  std::vector<TLorentzVector> tightEle;
  std::vector<TLorentzVector> looseEle;

  int ok=0, total=0;
  
  // Data structures to store info from TTrees
  baconhep::TEventInfo *info  	= new baconhep::TEventInfo();
  baconhep::TGenEventInfo *gen	= new baconhep::TGenEventInfo();
  TClonesArray *genPartArr 	= new TClonesArray("baconhep::TGenParticle");
  TClonesArray *muonArr    	= new TClonesArray("baconhep::TMuon");
  TClonesArray *electronArr	= new TClonesArray("baconhep::TElectron");
  TClonesArray *vertexArr	= new TClonesArray("baconhep::TVertex");
  TClonesArray *jetArr		= new TClonesArray("baconhep::TJet");
  TClonesArray *lheWgtArr	= new TClonesArray("baconhep::TLHEWeight");
  

  char command1[3000];
  char command2[3000];
  if ( cluster == "lxplus")
  	sprintf(command1, "eos find -f %s  | awk '!/log|fail/ {print $1}' | awk 'NF {print \"root://eoscms.cern.ch/\"$1}' > listTemp_%s.txt", (inputFolder).c_str(), outputFile.c_str());	// NF in awk command skips the blank line
  else 
	sprintf(command1,"xrdfs root://cmseos.fnal.gov ls %s | awk '{print \"root://cmseos.fnal.gov/\"$1}' > listTemp_%s.txt",(inputFolder).c_str(),  outputFile.c_str());
	//sprintf(command1,"eos root://cmseos.fnal.gov find -f %s | awk '!/log|fail/ {print $1}' | awk 'NF {print \"root://cmseos.fnal.gov/\"$1}' > listTemp_%s.txt",(inputFolder).c_str(),  outputFile.c_str());	// WORKS ONLY WITH INTERACTIVE NODE

  std::cout<<command1<<std::endl;
  sprintf(command2,"sed -i '/failed$/d' listTemp_%s.txt", outputFile.c_str());
  system(command1);
  system(command2);
  char list1[2000];
  sprintf (list1, "listTemp_%s.txt", outputFile.c_str());
  ifstream rootList (list1);
  char command3[300];
  sprintf(command3, "rm listTemp_%s.txt", outputFile.c_str());
  system(command3);

  int fileCounter=0;

  vector<TString>  sampleName; 

  while (!rootList.eof())
  {
	char iRun_tW[700];
	rootList >> iRun_tW;
	if(!rootList.good())break;
	sampleName.push_back(iRun_tW);
	fileCounter++;
  }

  TFile *infile=0;
  TTree *eventTree=0;

  std::ofstream file1;
  char TxtFileName[300];
  sprintf(TxtFileName, "%s.txt",outputFile.c_str());
  file1.open(TxtFileName);
  int cutEff[21]={0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};

  // Read and add pileup in histogram
  TFile* pileupFileMC = TFile::Open("puWeights_80x_37ifb.root");
  TH1D* puWeights = (TH1D*)pileupFileMC->Get("puWeights");
  TH1D* puWeightsUp = (TH1D*)pileupFileMC->Get("puWeightsUp");
  TH1D* puWeightsDown = (TH1D*)pileupFileMC->Get("puWeightsDown");
  puWeights->SetBins(75,0,75);
  puWeightsUp->SetBins(75,0,75);
  puWeightsDown->SetBins(75,0,75);

  //---------------- Root Files for ID, ISO, Trigger, GSF correctiosn for ELE and MU both: Starts -------------
  TFile* IDIsoEle = TFile::Open("egammaEffi_EGM2D_TightCutBasedIDSF.root","READ");
  TH1F *hIDIsoEle = (TH1F*)IDIsoEle->Get("EGamma_SF2D");

  TFile* GSFCorrEle = TFile::Open("egammaEffi_SF2D_GSF_tracking.root","READ");
  TH1F *hGSFCorrEle = (TH1F*)GSFCorrEle->Get("EGamma_SF2D");

  TFile* TriggerEle = TFile::Open("ElectronTrigger_SF.root","READ");
  TH1F* hTriggerEle = (TH1F*)TriggerEle->Get("HLT_Ele27");

  TFile* IDMuA = TFile::Open("MuonID_RunBCDEF_23SepReReco_19p72fb.root","READ");
  TH1F *hIDMuA = (TH1F*)IDMuA->Get("MC_NUM_TightID_DEN_genTracks_PAR_pt_eta/abseta_pt_ratio");

  TFile* IDMuB = TFile::Open("MuonID_RunGH_23SepReReco_16p146fb.root","READ");
  TH1F *hIDMuB = (TH1F*)IDMuB->Get("MC_NUM_TightID_DEN_genTracks_PAR_pt_eta/abseta_pt_ratio");

  TFile* IsoMuA = TFile::Open("MuonIso_RunBCDEF_23SepReReco_19p72fb.root","READ");
  TH1F *hIsoMuA = (TH1F*)IsoMuA->Get("TightISO_TightID_pt_eta/abseta_pt_ratio");

  TFile* IsoMuB = TFile::Open("MuonIso_RunGH_23SepReReco_16p146fb.root","READ");
  TH1F *hIsoMuB = (TH1F*)IsoMuB->Get("TightISO_TightID_pt_eta/abseta_pt_ratio");

  TFile* TriggerMuA = TFile::Open("MuonTrigger_RunBCDEF_23SepReReco_19p72fb.root","READ");
  TH1F* hTriggerMuA = (TH1F*)TriggerMuA->Get("IsoMu24_OR_IsoTkMu24_PtEtaBins/abseta_pt_ratio");

  TFile* TriggerMuB = TFile::Open("MuonTrigger_RunGH_23SepReReco_16p146fb.root","READ");
  TH1F* hTriggerMuB = (TH1F*)TriggerMuB->Get("IsoMu24_OR_IsoTkMu24_PtEtaBins/abseta_pt_ratio");

  //---------------- Root Files for ID, ISO, Trigger, GSF correctiosn for ELE and MU both: ENDS -------------
  

  //---------output tree----------------
  TFile* outROOT = TFile::Open((outputFile+(".root")).c_str(),"recreate");
  outROOT->cd();
  TTree* outTree = new TTree("otree", "otree");
  setOutputTree* WWTree = new setOutputTree(outTree);

  int nEvents=0;	
  Long64_t jentry2=0;
  int count_genEvents=0;

  int nInputFiles = sampleName.size();

  if (isLocal==1) nInputFiles = 2;
  cout<<"==> Total number of input files : "<<nInputFiles<<endl;

  TH1D *MCpu = new TH1D("MCpu","",75,0,75);
  TH1D *MCpu_up = new TH1D("MCpu_up","",75,0,75);
  TH1D *MCpu_down = new TH1D("MCpu_down","",75,0,75);
  
  Long64_t TotalNumberOfEvents = 0, nNegEvents = 0; 


  // Set up b-tag scale factor readers
  BTagCalibration calib("csvv2", "CSVv2_Moriond17_B_H.csv");
  BTagCalibrationReader bTagReader(BTagEntry::OP_LOOSE,  // working point: can be OP_LOOSE, OP_MEDIUM, OP_TIGHT 
                                   "central",             // label for the central value (see the scale factor file)
                                   {"up","down"});        // vector of labels for systematics
  bTagReader.load(calib, BTagEntry::FLAV_B, "comb");      // use the "comb" measurements for b-jets
  bTagReader.load(calib, BTagEntry::FLAV_C, "comb");      // use the "comb" measurements for c-jets
  bTagReader.load(calib, BTagEntry::FLAV_UDSG, "incl");   // use the "incl" measurements for light jets
  

  vector<const baconhep::TJet*> goodJetsv;
  
  // Loop on input files
  for(int i=0;i<nInputFiles;i++)
  {
     infile = TFile::Open(sampleName[i]);
     eventTree = (TTree*)infile->Get("Events");
     //std::cout << "\t File no. " << i << "\t"<< sampleName[i] <<std:: endl;
     
     TotalNumberOfEvents+=eventTree->GetEntries();

     if(isMC)
     { 
        eventTree->SetBranchAddress("Info", &info);    TBranch *infoBr = eventTree->GetBranch("Info");
  	TBranch *genBr=0;
     	eventTree->SetBranchAddress("GenEvtInfo", &gen); genBr = eventTree->GetBranch("GenEvtInfo");
	for (Long64_t jentry=0; jentry<eventTree->GetEntries();jentry++,jentry2++)
	{
	  //eventTree->GetEntry(jentry);
	    genBr->GetEntry(jentry);
    	    infoBr->GetEntry(jentry);	    
	    MCpu->Fill(info->nPUmean);
	    MCpu_up->Fill(info->nPUmeanp);
	    MCpu_down->Fill(info->nPUmeanm);
	    if (jentry2%50000 == 0) std::cout << "\t File no. " << i << "; Neg Event Count; read entry: " << jentry2 <<"/"<<TotalNumberOfEvents<<std:: endl;
	    if (gen->weight<0)	nNegEvents++;
	}
     }
     delete infile;
     infile=0, eventTree=0;
  }
  
  
  cout<<"==> Total number of events : "<<TotalNumberOfEvents<<endl;
  cout<<"==> Total number of negative events : "<<nNegEvents<<endl;

  float weight = std::atof(xSecWeight.c_str())/(std::atof(TotalNumberOfEntries.c_str()) - 2*std::atof(TotalNumberOfNegativeEntries.c_str()));
  cout<<"Weight of cross-sec/events = "<<weight<<endl;
  int totalEntries=0;


  JetCorrectorParameters paramAK4puppi("Summer16_23Sep2016V4_MC_JEC/Summer16_23Sep2016V4_MC_Uncertainty_AK4PFPuppi.txt");
  JetCorrectorParameters paramAK4chs("Summer16_23Sep2016V4_MC_JEC/Summer16_23Sep2016V4_MC_Uncertainty_AK4PFchs.txt");
  JetCorrectorParameters paramAK8chs("Summer16_23Sep2016V4_MC_JEC/Summer16_23Sep2016V4_MC_Uncertainty_AK8PFchs.txt");
  JetCorrectorParameters paramAK8puppi("Summer16_23Sep2016V4_MC_JEC/Summer16_23Sep2016V4_MC_Uncertainty_AK8PFPuppi.txt");
  JetCorrectionUncertainty *fJetUnc_AK4chs = new JetCorrectionUncertainty(paramAK4chs);

  //---------start loop on events------------
  std::cout << "---------start loop on events------------" << std::endl;
  jentry2=0;
  for(int i=0;i<nInputFiles;i++)
  {
  cout<<"\n\n=====	Processing File Number : "<<i<<"/"<<nInputFiles<<"\n\t"<<sampleName[i]<<"\n-------"<<endl;

  infile = TFile::Open(sampleName[i]);
  eventTree = (TTree*)infile->Get("Events");
  
  totalEntries+=eventTree->GetEntries();

  nEvents=eventTree->GetEntries();

  cout<<"\t==> Entries = "<<nEvents<<endl;



  eventTree->SetBranchAddress("Info", &info);    TBranch *infoBr = eventTree->GetBranch("Info");
  eventTree->SetBranchAddress("Muon", &muonArr); TBranch *muonBr = eventTree->GetBranch("Muon");
  eventTree->SetBranchAddress("Electron", &electronArr); TBranch *electronBr = eventTree->GetBranch("Electron");
  eventTree->SetBranchAddress("PV",   &vertexArr); TBranch *vertexBr = eventTree->GetBranch("PV");
  eventTree->SetBranchAddress("AK4CHS",   &jetArr); TBranch *jetBr = eventTree->GetBranch("AK4CHS");    
  TBranch *genBr=0, *genPartBr=0, *lhePartBr=0;
  if(isMC)
     { 
       eventTree->SetBranchAddress("GenEvtInfo", &gen); genBr = eventTree->GetBranch("GenEvtInfo");
       eventTree->SetBranchAddress("GenParticle",&genPartArr); genPartBr = eventTree->GetBranch("GenParticle");
       if(eventTree->GetListOfBranches()->FindObject("LHEWeight"))
       {
       eventTree->SetBranchAddress("LHEWeight",&lheWgtArr); lhePartBr = eventTree->GetBranch("LHEWeight");	       }
     }

  //for (Long64_t jentry=0; jentry<172; jentry++,jentry2++)
  for (Long64_t jentry=0; jentry<eventTree->GetEntries();jentry++,jentry2++)
  {
    //if (jentry2 != 87 && jentry2 != 113) continue;	// For debug
    infoBr->GetEntry(jentry);	    

    int GenPassCut = 0;

    tightMuon.clear();
    tightEle.clear();
    looseMuon.clear();
    looseEle.clear();
    

    if (jentry2%10000 == 0) std::cout << "\tread entry: " << jentry2 <<"/"<<TotalNumberOfEvents<<std:: endl;
    
    //*********************************
    WWTree->initializeVariables(); //initialize all variables

    WWTree->run   = info->runNum;
    WWTree->event = info->evtNum;
    WWTree->lumi  = info->lumiSec;


    /////////////////MC Info
    if (isMC==1)
    {
      lheWgtArr->Clear();
      if(lhePartBr)
	{
	  lhePartBr->GetEntry(jentry);
	}
      genPartArr->Clear();
      genBr->GetEntry(jentry);
      genPartBr->GetEntry(jentry);
      TLorentzVector hadW, lepW, VBFJ1, VBFJ2, VBFJ, temp;
      TLorentzVector genLep, genNeutrino, genWquarks, genVBFquarks;
      std::vector<TLorentzVector> v_GEN_hadW, v_GEN_lepW, v_GEN_VBFJ1, v_GEN_VBFJ2, v_GEN_VBFJ, v_GEN_temp;
      std::vector<TLorentzVector> v_genLep, v_genNeutrino, v_genWquarks, v_genVBFquarks;

      v_GEN_hadW.clear();	v_GEN_lepW.clear();	v_GEN_VBFJ1.clear();	v_GEN_VBFJ2.clear();
      v_GEN_VBFJ.clear();	v_GEN_temp.clear();	v_genLep.clear();	v_genNeutrino.clear();
      v_genWquarks.clear();	v_genVBFquarks.clear();
      
      for (int i = 0; i<genPartArr->GetEntries();i++)
	{
	  const baconhep::TGenParticle* genloop = (baconhep::TGenParticle*) ((*genPartArr)[i]);
	  Int_t parentPdg=dynamic_cast<baconhep::TGenParticle *>(genPartArr->At(genloop->parent>-1 ? genloop->parent : 0))->pdgId;
	  if( (abs(genloop->pdgId) == 11 || abs(genloop->pdgId) == 13 ) && abs(parentPdg) == 24)
	    {
	      genLep.SetPtEtaPhiM(genloop->pt, genloop->eta, genloop->phi, genloop->mass);
	      v_genLep.push_back(genLep);
	    }
	  if( (abs(genloop->pdgId) == 12 || abs(genloop->pdgId) == 14 ) && abs(parentPdg) == 24)
	    {
	      genNeutrino.SetPtEtaPhiM(genloop->pt, genloop->eta, genloop->phi, genloop->mass);
	      v_genNeutrino.push_back(genNeutrino);
	    }
	  if( (abs(genloop->pdgId) == 1 || abs(genloop->pdgId) == 3 || abs(genloop->pdgId) == 5 || abs(genloop->pdgId) == 2 || abs(genloop->pdgId) == 4 || abs(genloop->pdgId) == 6) && abs(parentPdg) == 24)
	    {
	      genWquarks.SetPtEtaPhiM(genloop->pt, genloop->eta, genloop->phi, genloop->mass);
	      v_genWquarks.push_back(genWquarks);
	    }
	  if( (abs(genloop->pdgId) == 1 || abs(genloop->pdgId) == 3 || abs(genloop->pdgId) == 5 || abs(genloop->pdgId) == 2 || abs(genloop->pdgId) == 4 || abs(genloop->pdgId) == 6) && genloop->status == 23 && abs(parentPdg) != 24)
	    {
	      genVBFquarks.SetPtEtaPhiM(genloop->pt, genloop->eta, genloop->phi, genloop->mass);
	      v_genVBFquarks.push_back(genVBFquarks);
	    }
	}
      if (v_genLep.size()==1 && v_genNeutrino.size()==1 && v_genWquarks.size()==2 && v_genVBFquarks.size()==2)
	{
	  WWTree->isGen           = 1;
	  WWTree->lep_pt_gen      = v_genLep[0].Pt();
	  WWTree->lep_eta_gen     = v_genLep[0].Eta();
	  WWTree->nu_pz_gen	  = v_genNeutrino[0].Pz();
	  WWTree->nu_pt_gen	  = v_genNeutrino[0].Pt();
	  WWTree->nu_eta_gen	  = v_genNeutrino[0].Eta();
	  WWTree->hadW_pt_gen     = (v_genWquarks[0]+v_genWquarks[1]).Pt();
	  WWTree->hadW_eta_gen    = (v_genWquarks[0]+v_genWquarks[1]).Eta();
	  WWTree->hadW_phi_gen    = (v_genWquarks[0]+v_genWquarks[1]).Phi();
	  WWTree->hadW_e_gen      = (v_genWquarks[0]+v_genWquarks[1]).E();
	  WWTree->hadW_m_gen      = (v_genWquarks[0]+v_genWquarks[1]).M();
	  
	  WWTree->lepW_pt_gen     = (v_genLep[0]+v_genNeutrino[0]).Pt();
	  WWTree->lepW_eta_gen    = (v_genLep[0]+v_genNeutrino[0]).Eta();
	  WWTree->lepW_phi_gen    = (v_genLep[0]+v_genNeutrino[0]).Phi();
	  WWTree->lepW_e_gen      = (v_genLep[0]+v_genNeutrino[0]).E();
	  WWTree->lepW_m_gen      = (v_genLep[0]+v_genNeutrino[0]).M();

	  WWTree->WW_mass_gen	= (v_genLep[0] + v_genNeutrino[0] + v_genWquarks[0] + v_genWquarks[1]).M();
	  WWTree->WW_mT_gen	= (v_genLep[0] + v_genNeutrino[0] + v_genWquarks[0] + v_genWquarks[1]).Mt();
	  WWTree->WW_pT_gen	= (v_genLep[0] + v_genNeutrino[0] + v_genWquarks[0] + v_genWquarks[1]).Pt();

	  WWTree->AK4_1_pt_gen	= v_genVBFquarks[0].Pt();
	  WWTree->AK4_1_eta_gen	= v_genVBFquarks[0].Eta();
	  WWTree->AK4_1_phi_gen	= v_genVBFquarks[0].Phi();
	  WWTree->AK4_1_e_gen	= v_genVBFquarks[0].E();
	  WWTree->AK4_1_mass_gen	= v_genVBFquarks[0].M();
	  
	  WWTree->AK4_2_pt_gen	= v_genVBFquarks[1].Pt();
	  WWTree->AK4_2_eta_gen	= v_genVBFquarks[1].Eta();
	  WWTree->AK4_2_phi_gen	= v_genVBFquarks[1].Phi();
	  WWTree->AK4_2_e_gen	= v_genVBFquarks[1].E();
	  WWTree->AK4_2_mass_gen	= v_genVBFquarks[1].M();
	  
	  WWTree->AK4_jj_DeltaEta_gen = abs(v_genVBFquarks[0].Eta() - v_genVBFquarks[1].Eta());
	  WWTree->AK4_jj_mass_gen = (v_genVBFquarks[0] + v_genVBFquarks[1]).M();
	  
	  count_genEvents++;
	}
	if ( WWTree->lep_pt_gen > 30 && abs(WWTree->lep_eta_gen) < 2.5 && WWTree->AK4_jj_DeltaEta_gen>2 && WWTree->hadW_pt_gen>200 && WWTree->AK4_1_pt_gen >30 && WWTree->AK4_2_pt_gen>30 && WWTree->AK4_jj_mass_gen > 500 ){
	        GenPassCut = 1;
	 }
    	for (int i = 0; i<lheWgtArr->GetEntries();i++)     // Note that i is starting from 446.
	{
		const baconhep::TLHEWeight *lhe = (baconhep::TLHEWeight*)((*lheWgtArr)[i]);
		//WWTree->LHEid[i] = lhe->id;
		WWTree->LHEWeight[i] = lhe->weight;
	}
    }

    
    WWTree->issignal = 0;
    WWTree->wSampleWeight = weight; //xsec/TotalNumberOfEntries
    WWTree->pu_Weight = 1.; //temporary value
    WWTree->pu_Weight_up = 1.; //temporary value
    WWTree->pu_Weight_down = 1.; //temporary value
    WWTree->top1_NNLO_Weight = 1.;
    WWTree->top2_NNLO_Weight = 1.;
    WWTree->id_eff_Weight = 1.;
    WWTree->id_eff_Weight2 = 1.;
    WWTree->trig_eff_Weight = 1.;
    WWTree->trig_eff_Weight2 = 1.;

    if (gen->weight>0)
      WWTree->genWeight=1.;
    else if (gen->weight<0) {
      WWTree->genWeight=-1.;
      //nNegEvents++;
    }
    cutEff[0]++;

    if (isMC==1)
    {
    if (GenPassCut == 1)   cutEff[1]++;
    }

    
    vertexArr->Clear();
    vertexBr->GetEntry(jentry);
    WWTree->nPV = vertexArr->GetEntries();
  
    //PILE-UP WEIGHT
    if (isMC==1) {
       if(int(info->nPUmean)<75){
           WWTree->pu_Weight = puWeights->GetBinContent(info->nPUmean); //our pu recipe
           WWTree->pu_Weight_up = puWeightsUp->GetBinContent(info->nPUmean); //our pu recipe
           WWTree->pu_Weight_down = puWeightsDown->GetBinContent(info->nPUmean); //our pu recipe
       }
       else{ //should not happen as we have a weight for all simulated n_pu multiplicities!
           std::cout<<"Warning! n_pu too big"<<std::endl;
	   // throw logic_error("n_pu too big");
	   WWTree->pu_Weight = 0.;
	   WWTree->pu_Weight_up = 0.;
	   WWTree->pu_Weight_down = 0.;
       } 
    }



    if(applyTrigger==1)
      if(!(triggerMenu.pass("HLT_IsoMu24_v*",info->triggerBits) || triggerMenu.pass("HLT_IsoTkMu24_v*",info->triggerBits) ||  triggerMenu.pass("HLT_Ele27_WPTight_Gsf_v*",info->triggerBits))) continue;
  
    /////////////////THE SELECTED LEPTON
    GeneralizedEndpoint ScaleSystematic;
    int nTightEle=0, nLooseEle=0;
    int nTightMu=0, nLooseMu=0;
    double pt_cut = 20;
    double leadelept_cut = 30;
    double leadmupt_cut = 27;
    electronArr->Clear();
    electronBr->GetEntry(jentry);
    const baconhep::TElectron *leadele = NULL;
    const baconhep::TElectron *subele = NULL;
    double leadeleE=-999, subeleE=-999;
    double iso = 1.5;
    for (int i=0; i<electronArr->GetEntries(); i++) {
      const baconhep::TElectron *ele = (baconhep::TElectron*)((*electronArr)[i]);
      if (ele->pt<=pt_cut) continue;
      if (fabs(ele->eta)>=2.5) continue;
      if(!passEleLooseSel(ele,info->rhoIso)) continue;
      nLooseEle++;
      if(!passEleTightSel(ele,info->rhoIso)) continue;
      ELE.SetPtEtaPhiM(ele->pt,ele->eta,ele->phi,0.0005109989461);
      tightEle.push_back(ELE);
      nTightEle++;
      iso = ele->chHadIso + TMath::Max( 0.0,(ele->gammaIso + ele->neuHadIso - info->rhoIso*eleEffArea(ele->eta)) );
      if(!leadele || ele->pt>leadele->pt)
	{
	  if(!(ele->pt>leadelept_cut)) continue;
	  subele = leadele;
	  leadele = ele;
	  leadeleE = ELE.E();
	  WWTree->l_iso1 = iso/ele->pt;
	}
      else if (!subele || ele->pt > subele->pt)
	{
	  subele = ele;
	  subeleE = ELE.E();
	  WWTree->l_iso2 = iso/ele->pt;
	}
    }
    if(leadele)
      {
	WWTree->l_pt1      = leadele->pt;
	WWTree->l_pt1_Up   = leadele->pt + 0.01*leadele->pt;
	WWTree->l_pt1_Down = leadele->pt - 0.01*leadele->pt;
	WWTree->l_eta1     = leadele->eta;
	WWTree->l_phi1     = leadele->phi;	
	WWTree->l_e1 	   = leadeleE;
	WWTree->l_charge1  = leadele->q;
      }
    if(subele)
      {
	WWTree->l_pt2      = subele->pt;
	WWTree->l_pt2_Up   = subele->pt + 0.01*subele->pt;
	WWTree->l_pt2_Down = subele->pt - 0.01*subele->pt;
	WWTree->l_eta2     = subele->eta;
	WWTree->l_phi2     = subele->phi;	
	WWTree->l_e2       = subeleE;
	WWTree->l_charge2  = subele->q;
      }
    muonArr->Clear();
    muonBr->GetEntry(jentry);
    const baconhep::TMuon *leadmu = NULL;
    const baconhep::TMuon *submu = NULL;
    double leadmue=-999, submue = -999;
    iso = 1.5;
    for(Int_t i=0; i<muonArr->GetEntries(); i++) {
      const baconhep::TMuon *mu = (baconhep::TMuon*)((*muonArr)[i]);
      if (mu->pt<pt_cut) continue;
      if (fabs(mu->eta)>=2.4) continue;
      if(!passMuonLooseSel(mu)) continue;
      nLooseMu++;
      if(!passMuonTightSel(mu)) continue;
      nTightMu++;
      MU.SetPtEtaPhiM(mu->pt,mu->eta,mu->phi,0.1056583745);
      tightMuon.push_back(MU);
      iso = mu->chHadIso + TMath::Max(mu->neuHadIso + mu->gammaIso - 0.5*(mu->puIso), double(0));
      if(!leadmu || mu->pt>leadmu->pt)
	{
	  if(!(mu->pt>leadmupt_cut)) continue;
	  submu = leadmu;
	  leadmu = mu;
	  leadmue = MU.E();
	  WWTree->l_iso1 = iso/mu->pt;
	}
      else if (!submu || mu->pt > submu->pt)
	{
	  submu = mu;
	  submue = MU.E();
	  WWTree->l_iso2 = iso/mu->pt;
	}
    }   
    if(leadmu)
      {
	if(leadmu->pt>500.0){
	WWTree->l_pt1  = ScaleSystematic.GeneralizedEndpointPt(leadmu->pt, leadmu->q, leadmu->eta, leadmu->phi, 0, 0);
	WWTree->l_pt1_Up  = ScaleSystematic.GeneralizedEndpointPt(leadmu->pt, leadmu->q, leadmu->eta, leadmu->phi, 1, 0);
	WWTree->l_pt1_Down  = ScaleSystematic.GeneralizedEndpointPt(leadmu->pt, leadmu->q, leadmu->eta, leadmu->phi, 2, 0);
	} else{
	WWTree->l_pt1      = leadmu->pt;
	WWTree->l_pt1_Up   = leadmu->pt + 0.01*leadmu->pt;
	WWTree->l_pt1_Down = leadmu->pt - 0.01*leadmu->pt;
	}
	WWTree->l_eta1 = leadmu->eta;
	WWTree->l_phi1 = leadmu->phi;	
	WWTree->l_e1 = leadmue;	
	WWTree->l_charge1 = leadmu->q;
	//if ( (WWTree->l_charge1 != 1) or (WWTree->l_charge1!=-1))
	//cout<<"Charge = "<< WWTree->l_charge1<<endl;
      }
    if(submu)
      {
	if(submu->pt>500.0){
	WWTree->l_pt2  = ScaleSystematic.GeneralizedEndpointPt(submu->pt, submu->q, submu->eta, submu->phi, 0, 0);
	WWTree->l_pt2_Up  = ScaleSystematic.GeneralizedEndpointPt(submu->pt, submu->q, submu->eta, submu->phi, 1, 0);
	WWTree->l_pt2_Down  = ScaleSystematic.GeneralizedEndpointPt(submu->pt, submu->q, submu->eta, submu->phi, 2, 0);
	} else{
	WWTree->l_pt2      = submu->pt;
	WWTree->l_pt2_Up   = submu->pt + 0.01*submu->pt;
	WWTree->l_pt2_Down = submu->pt - 0.01*submu->pt;
	}
	WWTree->l_eta2 = submu->eta;
	WWTree->l_phi2 = submu->phi;	
	WWTree->l_e2 = submue;	
	WWTree->l_charge2 = submu->q;
      }
    if(!(WWTree->l_pt1>0)) continue;
    if ((nTightMu+nTightEle)==0) continue; //no leptons with required ID
    if((nLooseEle+nLooseMu)>2) continue;
    if(nTightMu>0 && nLooseEle>0) continue;
    if(nTightEle>0 && nLooseMu>0) continue;
    if(nTightMu==1 && nLooseMu>1) continue;
    if(nTightEle==1 && nLooseEle>1) continue;
    if(nTightMu>0){
      WWTree->type=0;
      leptonName = "mu";	// Added this part for neutrino pz calculation in case there is w-boson.
    }else{
      WWTree->type=1;
      leptonName = "el";
    }
    
    if (WWTree->l_pt2<0)	cutEff[2]++;
    else cutEff[12]++;

    
    LEP1.SetPtEtaPhiM(WWTree->l_pt1,WWTree->l_eta1,WWTree->l_phi1,0.0005109989461);
    LEP2.SetPtEtaPhiM(WWTree->l_pt2,WWTree->l_eta2,WWTree->l_phi2,0.0005109989461);
    LEP1_Up.SetPtEtaPhiM(WWTree->l_pt1_Up,WWTree->l_eta1,WWTree->l_phi1,0.0005109989461);
    LEP2_Up.SetPtEtaPhiM(WWTree->l_pt2_Up,WWTree->l_eta2,WWTree->l_phi2,0.0005109989461);
    LEP1_Down.SetPtEtaPhiM(WWTree->l_pt1_Down,WWTree->l_eta1,WWTree->l_phi1,0.0005109989461);
    LEP2_Down.SetPtEtaPhiM(WWTree->l_pt2_Down,WWTree->l_eta2,WWTree->l_phi2,0.0005109989461);
    if(WWTree->l_pt2>0)
      {
	WWTree->dilep_pt  = (LEP1+LEP2).Pt();
	WWTree->dilep_eta = (LEP1+LEP2).Eta();
	WWTree->dilep_phi = (LEP1+LEP2).Phi();	
	WWTree->dilep_m = (LEP1+LEP2).M();	

	WWTree->dilep_pt_Up  = (LEP1_Up+LEP2_Up).Pt();
	WWTree->dilep_eta_Up = (LEP1_Up+LEP2_Up).Eta();
	WWTree->dilep_phi_Up = (LEP1_Up+LEP2_Up).Phi();	
	WWTree->dilep_m_Up = (LEP1_Up+LEP2_Up).M();	

	WWTree->dilep_pt_Down  = (LEP1_Down+LEP2_Down).Pt();
	WWTree->dilep_eta_Down = (LEP1_Down+LEP2_Down).Eta();
	WWTree->dilep_phi_Down = (LEP1_Down+LEP2_Down).Phi();	
	WWTree->dilep_m_Down = (LEP1_Down+LEP2_Down).M();	
      }

    //ID & GSF efficiency SF for electrons (https://twiki.cern.ch/twiki/bin/view/CMS/EgammaIDRecipesRun2#Electron_efficiencies_and_scale)
    if (strcmp(leptonName.c_str(),"el")==0 && isMC==1) {
	//  apply ID, ISO SF's
	WWTree->id_eff_Weight = GetSFs_Lepton(WWTree->l_pt1, WWTree->l_eta1, hIDIsoEle);	// Get Scale factor corresponding to the pt and eta.
	
	// apply GSF/RECO SF's for electrons
	WWTree->id_eff_Weight = WWTree->id_eff_Weight*GetSFs_Lepton(WWTree->l_pt1, WWTree->l_eta1, hGSFCorrEle);
	WWTree->trig_eff_Weight = 1.0/(GetSFs_Lepton(WWTree->l_pt1, WWTree->l_eta1, hTriggerEle));
    	
	if(WWTree->l_pt2>0){
	   //  apply ID, ISO SF's
	   WWTree->id_eff_Weight2 = GetSFs_Lepton(WWTree->l_pt2, WWTree->l_eta2, hIDIsoEle);	// Get Scale factor corresponding to the pt and eta.

	   // apply GSF/RECO SF's for electrons
	   WWTree->id_eff_Weight2 = WWTree->id_eff_Weight2*GetSFs_Lepton(WWTree->l_pt2, WWTree->l_eta2, hGSFCorrEle);
	   WWTree->trig_eff_Weight2 = 1.0/(GetSFs_Lepton(WWTree->l_pt2, WWTree->l_eta2, hTriggerEle));
	}
    }

    //ID&ISO efficiency SF for muons (https://twiki.cern.ch/twiki/bin/view/CMS/MuonWorkInProgressAndPagResults)
    if (strcmp(leptonName.c_str(),"mu")==0 && isMC==1) {
	//  apply ID SF's
	if (WWTree->run<278820){
		WWTree->id_eff_Weight = GetSFs_Lepton(WWTree->l_pt1, abs(WWTree->l_eta1), hIDMuA);
		WWTree->id_eff_Weight_Up = GetSFs_Lepton(WWTree->l_pt1_Up, abs(WWTree->l_eta1), hIDMuA);
		WWTree->id_eff_Weight_Down = GetSFs_Lepton(WWTree->l_pt1_Down, abs(WWTree->l_eta1), hIDMuA);
		if (WWTree->l_pt2>0) {
			WWTree->id_eff_Weight2 = GetSFs_Lepton(WWTree->l_pt2, abs(WWTree->l_eta2), hIDMuA);
			WWTree->id_eff_Weight2_Up = GetSFs_Lepton(WWTree->l_pt2_Up, abs(WWTree->l_eta2), hIDMuA);
			WWTree->id_eff_Weight2_Down = GetSFs_Lepton(WWTree->l_pt2_Down, abs(WWTree->l_eta2), hIDMuA);
		}
	}
	else{
		WWTree->id_eff_Weight = GetSFs_Lepton(WWTree->l_pt1, abs(WWTree->l_eta1), hIDMuB);
		WWTree->id_eff_Weight_Up = GetSFs_Lepton(WWTree->l_pt1_Up, abs(WWTree->l_eta1), hIDMuB);
		WWTree->id_eff_Weight_Down = GetSFs_Lepton(WWTree->l_pt1_Down, abs(WWTree->l_eta1), hIDMuB);
		if (WWTree->l_pt2>0) {
			WWTree->id_eff_Weight2 = GetSFs_Lepton(WWTree->l_pt2, abs(WWTree->l_eta2), hIDMuB);
			WWTree->id_eff_Weight2_Up = GetSFs_Lepton(WWTree->l_pt2_Up, abs(WWTree->l_eta2), hIDMuB);
			WWTree->id_eff_Weight2_Down = GetSFs_Lepton(WWTree->l_pt2_Down, abs(WWTree->l_eta2), hIDMuB);
		}
	}

	//  apply ISO SF's
	if (WWTree->run<278820){
		WWTree->id_eff_Weight = WWTree->id_eff_Weight*GetSFs_Lepton(WWTree->l_pt1, abs(WWTree->l_eta1), hIsoMuA);
		WWTree->id_eff_Weight_Up = WWTree->id_eff_Weight_Up*GetSFs_Lepton(WWTree->l_pt1_Up, abs(WWTree->l_eta1), hIsoMuA);
		WWTree->id_eff_Weight_Down = WWTree->id_eff_Weight_Down*GetSFs_Lepton(WWTree->l_pt1_Down, abs(WWTree->l_eta1), hIsoMuA);

		WWTree->id_eff_Weight2 = WWTree->id_eff_Weight2*GetSFs_Lepton(WWTree->l_pt2, abs(WWTree->l_eta2), hIsoMuA);
		WWTree->id_eff_Weight2_Up = WWTree->id_eff_Weight2_Up*GetSFs_Lepton(WWTree->l_pt2_Up, abs(WWTree->l_eta2), hIsoMuA);
		WWTree->id_eff_Weight2_Down = WWTree->id_eff_Weight2_Down*GetSFs_Lepton(WWTree->l_pt2_Down, abs(WWTree->l_eta2), hIsoMuA);
		}
	else{
		WWTree->id_eff_Weight = WWTree->id_eff_Weight*GetSFs_Lepton(WWTree->l_pt1, abs(WWTree->l_eta1), hIsoMuB);
		WWTree->id_eff_Weight_Up = WWTree->id_eff_Weight_Up*GetSFs_Lepton(WWTree->l_pt1_Up, abs(WWTree->l_eta1), hIsoMuB);
		WWTree->id_eff_Weight_Down = WWTree->id_eff_Weight_Down*GetSFs_Lepton(WWTree->l_pt1_Down, abs(WWTree->l_eta1), hIsoMuB);

		WWTree->id_eff_Weight2 = WWTree->id_eff_Weight2*GetSFs_Lepton(WWTree->l_pt2, abs(WWTree->l_eta2), hIsoMuB);
		WWTree->id_eff_Weight2_Up = WWTree->id_eff_Weight2_Up*GetSFs_Lepton(WWTree->l_pt2_Up, abs(WWTree->l_eta2), hIsoMuB);
		WWTree->id_eff_Weight2_Down = WWTree->id_eff_Weight2_Down*GetSFs_Lepton(WWTree->l_pt2_Down, abs(WWTree->l_eta2), hIsoMuB);
		}
	
	// apply Trigger SF's
	if (WWTree->run<278820){
		WWTree->trig_eff_Weight  = GetSFs_Lepton(WWTree->l_pt1, abs(WWTree->l_eta1), hTriggerMuA);
		WWTree->trig_eff_Weight_Up  = GetSFs_Lepton(WWTree->l_pt1_Up, abs(WWTree->l_eta1_Up), hTriggerMuA);
		WWTree->trig_eff_Weight_Down  = GetSFs_Lepton(WWTree->l_pt1_Down, abs(WWTree->l_eta1_Down), hTriggerMuA);
		if (WWTree->l_pt2>0) {
			WWTree->trig_eff_Weight2 = GetSFs_Lepton(WWTree->l_pt2, abs(WWTree->l_eta2), hTriggerMuA);
			WWTree->trig_eff_Weight2_Up = GetSFs_Lepton(WWTree->l_pt2_Up, abs(WWTree->l_eta2), hTriggerMuA);
			WWTree->trig_eff_Weight2_Down = GetSFs_Lepton(WWTree->l_pt2_Down, abs(WWTree->l_eta2), hTriggerMuA);
		}
	}
	else{
		WWTree->trig_eff_Weight  = GetSFs_Lepton(WWTree->l_pt1, abs(WWTree->l_eta1), hTriggerMuB);
		WWTree->trig_eff_Weight_Up  = GetSFs_Lepton(WWTree->l_pt1_Up, abs(WWTree->l_eta1), hTriggerMuB);
		WWTree->trig_eff_Weight_Down  = GetSFs_Lepton(WWTree->l_pt1_Down, abs(WWTree->l_eta1), hTriggerMuB);
		if (WWTree->l_pt2>0) {
			WWTree->trig_eff_Weight2 = GetSFs_Lepton(WWTree->l_pt2, abs(WWTree->l_eta2), hTriggerMuB);
			WWTree->trig_eff_Weight2_Up = GetSFs_Lepton(WWTree->l_pt2_Up, abs(WWTree->l_eta2), hTriggerMuB);
			WWTree->trig_eff_Weight2_Down = GetSFs_Lepton(WWTree->l_pt2_Down, abs(WWTree->l_eta2), hTriggerMuB);
		}
        }
	}
	
    //////////////THE MET
    
    // //preselection on met
     if (info->pfMETC < 0) continue;
     if (WWTree->l_pt2<0) cutEff[3]++;
     else cutEff[13]++;
    
    // Calculate Neutrino Pz using all the possible choices : 
    // type0 -> if real roots, pick the one nearest to the lepton Pz except when the Pz so chosen
    //               is greater than 300 GeV in which case pick the most central root.
    // type1 -> type = 1: if real roots, choose the one closest to the lepton Pz if complex roots, use only the real part.
    //          type = 2: if real roots, choose the most central solution. if complex roots, use only the real part. 
    //          type = 3: if real roots, pick the largest value of the cosine*
    
    TLorentzVector W_Met;
    WWTree->pfMET_Corr = info->pfMETC;
    WWTree->pfMET_Corr_phi = info->pfMETCphi;
    if (WWTree->l_pt2<0) {
    float Wmass = 80.385;
  
    TLorentzVector W_Met_jes_up, W_Met_jes_dn, W_Met_jer_up, W_Met_jer_dn, AK4Up, AK4Down, AK4Up_Puppi, AK4Down_Puppi;
  
    W_Met.SetPxPyPzE(info->pfMETC * TMath::Cos(info->pfMETCphi), info->pfMETC * TMath::Sin(info->pfMETCphi), 0., sqrt(info->pfMETC*info->pfMETC));
    W_Met_jer_up.SetPxPyPzE(info->pfMETCjerup * TMath::Cos(info->pfMETCphijerup), info->pfMETCjerup * TMath::Sin(info->pfMETCphijerup), 0., sqrt(info->pfMETCjerup*info->pfMETCjerup) );
    W_Met_jer_dn.SetPxPyPzE(info->pfMETCjerdn * TMath::Cos(info->pfMETCphijerdn), info->pfMETCjerdn * TMath::Sin(info->pfMETCphijerdn), 0., sqrt(info->pfMETCjerdn*info->pfMETCjerdn) );

    ////////////////////////////////////////////////////////////////
    //		
    //		MET JES Calculate
    //
    ////////////////////////////////////////////////////////////////
    jetArr->Clear();
    jetBr->GetEntry(jentry);
    for ( int i=0; i<jetArr->GetEntries(); i++) //loop on AK4 jet
    {
      const baconhep::TJet *jet = (baconhep::TJet*)((*jetArr)[i]);
      TLorentzVector AK4_LV_temp, AK4_LV_temp2;

      // Get uncertanity
      double unc = func(jet->pt, jet->eta, fJetUnc_AK4chs); 

      // Get AK4 LorentzVector 
      AK4_LV_temp.SetPtEtaPhiM(jet->pt, jet->eta, jet->phi, jet->mass);

      // calculate Up variation
      AK4_LV_temp2.SetPtEtaPhiM((1.+unc)*jet->pt, jet->eta, jet->phi, (1.+unc)*jet->mass);
      AK4Up += AK4_LV_temp2 - AK4_LV_temp;

      // calculate Down variation
      AK4_LV_temp2.SetPtEtaPhiM((1.-unc)*jet->pt, jet->eta, jet->phi, (1.-unc)*jet->mass);
      AK4Down += AK4_LV_temp2 - AK4_LV_temp;
    }
    W_Met_jes_up = W_Met + AK4Up;
    W_Met_jes_dn = W_Met + AK4Down;
    //////////////////////////////////////// END: MET JES Calculate
    if(LEP1.Pt()<=0 || W_Met.Pt() <= 0 ){ std::cerr<<" Negative Lepton - Neutrino Pt "<<std::endl; continue ; }
    if (WWTree->l_pt2<0) cutEff[4]++;	// There is no MET in two lepton case. So, cutEff[14] will not increase. So, added cutEff[14] after this met loop


    // type0 calculation of neutrino pZ
    METzCalculator NeutrinoPz_type0;
    METzCalculator NeutrinoPz_type0_jes_up;
    METzCalculator NeutrinoPz_type0_jes_dn;
    METzCalculator NeutrinoPz_type0_jer_up;
    METzCalculator NeutrinoPz_type0_jer_dn;

    METzCalculator_Run2 NeutrinoPz_run2;
    NeutrinoPz_type0.SetMET(W_Met);
    NeutrinoPz_type0.SetLepton(LEP1);
    NeutrinoPz_type0.SetLeptonType(leptonName.c_str());
    
    NeutrinoPz_type0_jes_up.SetMET(W_Met_jes_up);
    NeutrinoPz_type0_jes_up.SetLepton(LEP1);
    NeutrinoPz_type0_jes_up.SetLeptonType(leptonName.c_str());
    
    NeutrinoPz_type0_jes_dn.SetMET(W_Met_jes_dn);
    NeutrinoPz_type0_jes_dn.SetLepton(LEP1);
    NeutrinoPz_type0_jes_dn.SetLeptonType(leptonName.c_str());

    NeutrinoPz_type0_jer_up.SetMET(W_Met_jer_up);
    NeutrinoPz_type0_jer_up.SetLepton(LEP1);
    NeutrinoPz_type0_jer_up.SetLeptonType(leptonName.c_str());

    NeutrinoPz_type0_jer_dn.SetMET(W_Met_jer_dn);
    NeutrinoPz_type0_jer_dn.SetLepton(LEP1);
    NeutrinoPz_type0_jer_dn.SetLeptonType(leptonName.c_str());

    
    NeutrinoPz_run2.SetMET(W_Met);
    NeutrinoPz_run2.SetLepton(LEP1);
    NeutrinoPz_run2.SetLeptonType(leptonName.c_str());
  
    double pz1_type0 = NeutrinoPz_type0.Calculate(); // Default one -> according to type0
    //double pz2_type0 = NeutrinoPz_type0.getOther();  // Default one
    
    double pz1_run2 = NeutrinoPz_run2.Calculate();
  
    double pz1_type0_jes_up = NeutrinoPz_type0_jes_up.Calculate(); // Default one -> according to type0
    double pz1_type0_jes_dn = NeutrinoPz_type0_jes_dn.Calculate(); // Default one -> according to type0

    double pz1_type0_jer_up = NeutrinoPz_type0_jer_up.Calculate();
    double pz1_type0_jer_dn = NeutrinoPz_type0_jer_dn.Calculate();
  
    // don't touch the neutrino pT
    TLorentzVector W_neutrino_type0_met; 
    W_neutrino_type0_met.SetPxPyPzE(W_Met.Px(),W_Met.Py(),pz1_type0,sqrt(W_Met.Pt()*W_Met.Pt()+pz1_type0*pz1_type0));
    
    // change the neutrino pT in case of complex solution in order to make it real
    TLorentzVector W_neutrino_type0; 
    W_neutrino_type0.SetPxPyPzE(W_Met.Px(),W_Met.Py(),pz1_type0,sqrt(W_Met.Pt()*W_Met.Pt()+pz1_type0*pz1_type0));
    
    if(NeutrinoPz_type0.IsComplex()) {// if this is a complex, change MET
      double nu_pt1 = NeutrinoPz_type0.getPtneutrino(1);
      double nu_pt2 = NeutrinoPz_type0.getPtneutrino(2);
      TLorentzVector W_neutrino_1;
      W_neutrino_1.SetPxPyPzE(nu_pt1 * TMath::Cos(info->pfMETCphi), nu_pt1 * TMath::Sin(info->pfMETCphi), pz1_type0, sqrt(nu_pt1*nu_pt1 + pz1_type0*pz1_type0) );
      TLorentzVector W_neutrino_2;
      W_neutrino_2.SetPxPyPzE(nu_pt2 * TMath::Cos(info->pfMETCphi), nu_pt2 * TMath::Sin(info->pfMETCphi), pz1_type0, sqrt(nu_pt2*nu_pt2 + pz1_type0*pz1_type0) );
      
      if ( fabs((LEP1+W_neutrino_1).M()-Wmass) < fabs((LEP1+W_neutrino_2).M()-Wmass) ) W_neutrino_type0 = W_neutrino_1;
      else W_neutrino_type0 = W_neutrino_2;
    }
    
    // type2 calculation of neutrino pZ
    METzCalculator NeutrinoPz_type2;
    NeutrinoPz_type2.SetMET(W_Met);
    NeutrinoPz_type2.SetLepton(LEP1);
    NeutrinoPz_type2.SetLeptonType(leptonName.c_str());
    
    double pz1_type2 = NeutrinoPz_type2.Calculate(2); // Default one -> according to type2
    //double pz2_type2 = NeutrinoPz_type2.getOther();   // Default one
  
    // don't touch the neutrino pT
    TLorentzVector W_neutrino_type2_met; 
    W_neutrino_type2_met.SetPxPyPzE(W_Met.Px(),W_Met.Py(),pz1_type2,sqrt(W_Met.Pt()*W_Met.Pt()+pz1_type2*pz1_type2));
    
    // change the neutrino pT in case of complex solution in order to make it real
    TLorentzVector W_neutrino_type2; 
    W_neutrino_type2.SetPxPyPzE(W_Met.Px(),W_Met.Py(),pz1_type2,sqrt(W_Met.Pt()*W_Met.Pt()+pz1_type2*pz1_type2));
    
    if (NeutrinoPz_type2.IsComplex()) {// if this is a complex, change MET
      double nu_pt1 = NeutrinoPz_type2.getPtneutrino(1);
      double nu_pt2 = NeutrinoPz_type2.getPtneutrino(2);
      TLorentzVector W_neutrino_1;
      W_neutrino_1.SetPxPyPzE(nu_pt1 * TMath::Cos(info->pfMETCphi), nu_pt1 * TMath::Sin(info->pfMETCphi), pz1_type2, sqrt(nu_pt1*nu_pt1 + pz1_type2*pz1_type2) );
      TLorentzVector W_neutrino_2;
      W_neutrino_2.SetPxPyPzE(nu_pt2 * TMath::Cos(info->pfMETCphi), nu_pt2 * TMath::Sin(info->pfMETCphi), pz1_type2, sqrt(nu_pt2*nu_pt2 + pz1_type2*pz1_type2) );
      
      if ( fabs((LEP1+W_neutrino_1).M()-Wmass) < fabs((LEP1+W_neutrino_2).M()-Wmass) ) W_neutrino_type2 = W_neutrino_1;
      else W_neutrino_type2 = W_neutrino_2;
    }
    
    WWTree->pfMET = sqrt(info->pfMET*info->pfMET);
    WWTree->pfMET_jes_up = W_Met_jes_up.Pt();		// Calculated with corrected pfMET
    WWTree->pfMET_jes_dn = W_Met_jes_dn.Pt();		// Calculated with corrected pfMET
    WWTree->pfMET_Phi = info->pfMETphi;
    WWTree->pfMET_Corr = info->pfMETC;
    WWTree->pfMET_Corr_phi = info->pfMETCphi;
    WWTree->pfMET_Corr_Cov00 = info->pfMETCCov00;
    WWTree->pfMET_Corr_Cov01 = info->pfMETCCov01;
    WWTree->pfMET_Corr_Cov11 = info->pfMETCCov11;
    WWTree->pfMET_Corr_jerup = info->pfMETCjerup;
    WWTree->pfMET_Corr_jerdn = info->pfMETCjerdn;
    WWTree->pfMET_Corr_jenup = info->pfMETCjenup;
    WWTree->pfMET_Corr_jendn = info->pfMETCjendn;
    WWTree->pfMET_Corr_uncup = info->pfMETCuncup;
    WWTree->pfMET_Corr_uncdn = info->pfMETCuncdn;
    WWTree->pfMET_Corr_jrsup = info->pfMETCjrsup;
    WWTree->pfMET_Corr_jrsdn = info->pfMETCjrsdn;
    WWTree->pfMET_Corr_phijerup = info->pfMETCphijerup;
    WWTree->pfMET_Corr_phijerdn = info->pfMETCphijerdn;
    WWTree->pfMET_Corr_phijenup = info->pfMETCphijenup;
    WWTree->pfMET_Corr_phijendn = info->pfMETCphijendn;
    WWTree->pfMET_Corr_phiuncup = info->pfMETCphiuncup;
    WWTree->pfMET_Corr_phiuncdn = info->pfMETCphiuncdn;
    WWTree->pfMET_Corr_phijrsup = info->pfMETCphijrsup;
    WWTree->pfMET_Corr_phijrsdn = info->pfMETCphijrsdn;

    WWTree->nu_pz_type0 = pz1_type0;
    WWTree->nu_pz_type2 = pz1_type2;
    WWTree->nu_pz_run2 = pz1_run2;
    WWTree->nu_pz_isre = 1-NeutrinoPz_run2.IsComplex();
    WWTree->nu_pz_run2_oth = NeutrinoPz_run2.getOther();
    WWTree->nu_pz_run2_type = NeutrinoPz_run2.getType();
  
    
    /////////////////THE LEPTONIC W
  
    NU0.SetPtEtaPhiM( GetPt_MET(info->pfMETC, info->pfMETCphi, WWTree->nu_pz_type0) , GetEta_MET(info->pfMETC, info->pfMETCphi, WWTree->nu_pz_type0), info->pfMETCphi , 0.0 );

    NU0_jes_up.SetPtEtaPhiM( GetPt_MET(W_Met_jes_up.Pt(), W_Met_jes_up.Phi(), pz1_type0_jes_up), GetEta_MET(W_Met_jes_up.Pt(), W_Met_jes_up.Phi(), pz1_type0_jes_up), W_Met_jes_up.Phi() , 0.0 );
    NU0_jes_dn.SetPtEtaPhiM( GetPt_MET(W_Met_jes_dn.Pt(), W_Met_jes_dn.Phi(), pz1_type0_jes_dn), GetEta_MET(W_Met_jes_dn.Pt(), W_Met_jes_dn.Phi(), pz1_type0_jes_dn), W_Met_jes_dn.Phi() , 0.0 );
    NU0_jer_up.SetPtEtaPhiM( GetPt_MET(W_Met_jer_up.Pt(), W_Met_jer_up.Phi(), pz1_type0_jer_up), GetEta_MET(W_Met_jer_up.Pt(), W_Met_jer_up.Phi(), pz1_type0_jer_up), W_Met_jer_up.Phi() , 0.0);
    NU0_jer_dn.SetPtEtaPhiM( GetPt_MET(W_Met_jer_dn.Pt(), W_Met_jer_dn.Phi(), pz1_type0_jer_dn), GetEta_MET(W_Met_jer_dn.Pt(), W_Met_jer_dn.Phi(), pz1_type0_jer_dn), W_Met_jer_dn.Phi() , 0.0);


    
    NU2.SetPtEtaPhiM(GetPt_MET( info->pfMETC, info->pfMETCphi, WWTree->nu_pz_type2 ), GetEta_MET( info->pfMETC, info->pfMETCphi, WWTree->nu_pz_type2 ), info->pfMETCphi, 0.0);
    NU1.SetPtEtaPhiM(GetPt_MET( info->pfMETC, info->pfMETCphi, WWTree->nu_pz_run2  ), GetEta_MET( info->pfMETC, info->pfMETCphi, WWTree->nu_pz_run2  ), info->pfMETCphi, 0.0);

  
    W_type0 = LEP1 + NU0;
    W_type0_jes_up = LEP1 + NU0_jes_up;
    W_type0_jes_dn = LEP1 + NU0_jes_dn;
    W_type0_jer_up = LEP1 + NU0_jer_up;
    W_type0_jer_dn = LEP1 + NU0_jer_dn;
    W_type2 = LEP1 + NU2;
    W_run2 = LEP1 + NU1;

    W_type0_LEP_Up = LEP1_Up +NU0;
    W_type0_LEP_Down = LEP1_Down +NU0;
    

  
    WWTree->v_pt_type0 = W_type0.Pt();
    WWTree->v_pt_type0_LEP_Up = W_type0_LEP_Up.Pt();
    WWTree->v_pt_type0_LEP_Down = W_type0_LEP_Down.Pt();

    WWTree->v_eta_type0 = W_type0.Eta();
    WWTree->v_mt_type0 = TMath::Sqrt(2*LEP1.Et()*NU0.Et()*(1-TMath::Cos(LEP1.DeltaPhi(NU0))));
    WWTree->v_mass_type0 = W_type0.M();
    WWTree->v_mass_type0_LEP_Up = W_type0_LEP_Up.M();
    WWTree->v_mass_type0_LEP_Down = W_type0_LEP_Down.M();

    WWTree->v_pt_type0_jes_up = W_type0_jes_up.Pt();
    WWTree->v_eta_type0_jes_up = W_type0_jes_up.Eta();
    WWTree->v_mt_type0_jes_up = TMath::Sqrt(2*LEP1.Et()*NU0_jes_up.Et()*(1-TMath::Cos(LEP1.DeltaPhi(NU0_jes_up))));
    WWTree->v_mass_type0_jes_up = W_type0_jes_up.M();
    WWTree->v_pt_type0_jes_dn = W_type0_jes_dn.Pt();
    WWTree->v_eta_type0_jes_dn = W_type0_jes_dn.Eta();
    WWTree->v_mt_type0_jes_dn = TMath::Sqrt(2*LEP1.Et()*NU0_jes_dn.Et()*(1-TMath::Cos(LEP1.DeltaPhi(NU0_jes_dn))));
    WWTree->v_mass_type0_jes_dn = W_type0_jes_dn.M();
    WWTree->v_pt_type0_jer_up = W_type0_jer_up.Pt();
    WWTree->v_eta_type0_jer_up = W_type0_jer_up.Eta();
    WWTree->v_mt_type0_jer_up = TMath::Sqrt(2*LEP1.Et()*NU0_jer_up.Et()*(1-TMath::Cos(LEP1.DeltaPhi(NU0_jer_up))));
    WWTree->v_mass_type0_jer_up = W_type0_jer_up.M();
    WWTree->v_pt_type0_jer_dn = W_type0_jer_dn.Pt();
    WWTree->v_eta_type0_jer_dn = W_type0_jer_dn.Eta();
    WWTree->v_mt_type0_jer_dn = TMath::Sqrt(2*LEP1.Et()*NU0_jer_dn.Et()*(1-TMath::Cos(LEP1.DeltaPhi(NU0_jer_dn))));
    WWTree->v_mass_type0_jer_dn = W_type0_jer_dn.M();

    WWTree->v_pt_type2 = W_type2.Pt();
    WWTree->v_pt_run2 = W_run2.Pt();
    WWTree->v_eta_type2 = W_type2.Eta();
    WWTree->v_eta_run2 = W_run2.Eta();
    WWTree->v_phi = W_type2.Phi();
    WWTree->v_mt_type2 = TMath::Sqrt(2*LEP1.Et()*NU2.Et()*(1-TMath::Cos(LEP1.DeltaPhi(NU2))));
    WWTree->v_mt_run2 = TMath::Sqrt(2*LEP1.Et()*NU1.Et()*(1-TMath::Cos(LEP1.DeltaPhi(NU1))));
    WWTree->v_mass_type2 = W_type2.M();
    WWTree->v_mass_run2 = W_run2.M();
    }
    //if(LEP1.Pt()<=0 || LEP2.Pt() <= 0 ){ std::cerr<<" Negative Lepton - Neutrino Pt "<<std::endl; continue ; }
    if (WWTree->l_pt2>0) cutEff[14]++;  // There is no MET in two lepton case. So, cutEff[4] is placed in previous if condition.
  
      
    if (WWTree->l_pt2<0) 
    	WWTree->totalEventWeight = WWTree->genWeight*WWTree->pu_Weight*WWTree->top1_NNLO_Weight*WWTree->top2_NNLO_Weight*WWTree->trig_eff_Weight*WWTree->id_eff_Weight;
    else
        WWTree->totalEventWeight_2Lep = WWTree->genWeight*WWTree->pu_Weight*WWTree->top1_NNLO_Weight*WWTree->top2_NNLO_Weight*WWTree->trig_eff_Weight*WWTree->id_eff_Weight*WWTree->trig_eff_Weight2*WWTree->id_eff_Weight2;

    
    WWTree->nEvents = TotalNumberOfEvents;
    WWTree->nNegEvents = nNegEvents;
    WWTree->nTotEvents = std::atof(TotalNumberOfEntries.c_str());
    WWTree->nTotNegEvents = std::atof(TotalNumberOfNegativeEntries.c_str());

    outTree->Fill();

    goodJetsv.clear();
    }
    delete infile;
    infile=0, eventTree=0;
    /////////////////FILL THE TREE
  }
  //delete puWeight;	delete puWeight_up;	delete puWeight_down;
  delete MCpu;	delete MCpu_up;	delete MCpu_down;
  delete puWeightsDown;	delete puWeightsUp;	delete puWeights;
  //delete pileupHisto;
  //pileupFile->Close();
  pileupFileMC->Close();
  std::cout << "---------end loop on events------------" << std::endl;
  std::cout << std::endl;
  std::cout << "GEN events = " << count_genEvents << std::endl;


  
  std::cout << "----------------------" << std::endl;
  std::cout << " SUMMARY" << std::endl;
  std::cout << "----------------------" << std::endl;
  std::cout << std::endl;
  std::cout<<"MC matching: "<<(float)ok/(float)total<<std::endl;
  std::cout<<"negative events: "<<nNegEvents<<std::endl;
  std::cout << std::endl;
  std::cout<<"(0) all events:        "<<cutEff[0]<<"\t:\t"<<((float)cutEff[0]*100.0)/(float)cutEff[0]<<std::endl
  	   <<"(1) Gen Events:        "<<cutEff[1]<<"\t:\t"<<((float)cutEff[1]*100.0)/(float)cutEff[0]<<std::endl
	   <<"(2) tight lepton:      "<<cutEff[2]<<"\t:\t"<<((float)cutEff[2]*100.0)/(float)cutEff[0]<<std::endl
	   <<"(3) MET:               "<<cutEff[3]<<"\t:\t"<<((float)cutEff[3]*100.0)/(float)cutEff[2]<<std::endl
	   <<"(4) negative lep-MET:  "<<cutEff[4]<<"\t:\t"<<((float)cutEff[4]*100.0)/(float)cutEff[3]<<std::endl
	   <<"(5) 1 good AK8:        "<<cutEff[5]<<"\t:\t"<<((float)cutEff[5]*100.0)/(float)cutEff[4]<<std::endl
	   <<"(6) m(WV) > 0:       "<<cutEff[6]<<"\t:\t"<<((float)cutEff[6]*100.0)/(float)cutEff[5]<<std::endl
	   <<"(7) >=2 good VBF jets: "<<cutEff[7]<<"\t:\t"<<((float)cutEff[7]*100.0)/(float)cutEff[6]<<std::endl
	   <<"(8) Found VBF jets:  "<<cutEff[8]<<"\t:\t"<<((float)cutEff[8]*100.)/(float)cutEff[7]<<std::endl;
  
  std::cout<<"\n\n----------------------------"<< std::endl;
  std::cout<<"\tSUMMARY for 2 lepton case "<< std::endl;
  std::cout<<"---------------------------"<< std::endl;
  std::cout << std::endl;
  std::cout<<"(0) all events:        "<<cutEff[0]<<"\t:\t"<<((float)cutEff[0]*100.0)/(float)cutEff[0]<<std::endl
  	   <<"(1) Gen Events:        "<<cutEff[1]<<"\t:\t"<<((float)cutEff[1]*100.0)/(float)cutEff[0]<<std::endl
	   <<"(2) tight lepton:      "<<cutEff[12]<<"\t:\t"<<((float)cutEff[12]*100.0)/(float)cutEff[0]<<std::endl
	   <<"(3) MET:               "<<cutEff[13]<<"\t:\t"<<((float)cutEff[13]*100.0)/(float)cutEff[12]<<std::endl
	   <<"(4) negative lep-MET:  "<<cutEff[14]<<"\t:\t"<<((float)cutEff[14]*100.0)/(float)cutEff[13]<<std::endl
	   <<"(5) 1 good AK8:        "<<cutEff[15]<<"\t:\t"<<((float)cutEff[15]*100.0)/(float)cutEff[14]<<std::endl
	   <<"(6) m(WV) > 0:       "<<cutEff[16]<<"\t:\t"<<((float)cutEff[16]*100.0)/(float)cutEff[15]<<std::endl
	   <<"(7) >=2 good VBF jets: "<<cutEff[17]<<"\t:\t"<<((float)cutEff[17]*100.0)/(float)cutEff[16]<<std::endl
	   <<"(8) Found VBF jets:  "<<cutEff[18]<<"\t:\t"<<((float)cutEff[18]*100.)/(float)cutEff[17]<<std::endl;
  
 
 
  //--------close everything-------------
  delete info; delete gen;
  delete genPartArr; delete muonArr; delete electronArr; delete vertexArr;
  delete jetArr; delete lheWgtArr;
  outROOT->Write();
  outROOT->Close();
  int t1 = time(NULL);
  printf ("\n==> time to run this code = %0.3f min\n", (float)(t1 - t0)/60.0);
  return(0);
}
