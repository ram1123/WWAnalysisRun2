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
  
  string inputFolder = argv[1];
  string outputFile = argv[2];
  int isMC = atoi(argv[3]);
  string cluster = argv[4];
  string inputTreeName = argv[5];
  string inputFile = argv[6];
  string xSecWeight = argv[7];
  string TotalNumberOfEntries = argv[8];
  string TotalNumberOfNegativeEntries = argv[9];
  int applyTrigger = atoi(argv[10]);
  string jsonFileName = argv[11];
  int isLocal = atoi(argv[12]);
  int VBFSel  = atoi(argv[13]);
  
  string leptonName;

  if ( VBFSel==1)	cout<<"==> VBF selection method : Select two highest pT jets"<<endl;
  else if ( VBFSel==2)	cout<<"==> VBF selection method : Select pair with highest mjj..."<<endl;
  else if ( VBFSel==3)	cout<<"==> VBF selection method : Select pair with highest DeltaEta..."<<endl;
  else {	cout<<"\n\nERROR:	Enter valid vbf selection criteria....\n\n"<<endl;
  		exit(0);  
	}
  
  string iHLTFile="${CMSSW_BASE}/src/BaconAna/DataFormats/data/HLTFile_25ns";
  const string cmssw_base = getenv("CMSSW_BASE");
  string cmssw_base_env = "${CMSSW_BASE}";
  size_t start_pos = iHLTFile.find(cmssw_base_env);
  if(start_pos != string::npos) {
  	iHLTFile.replace(start_pos, cmssw_base_env.length(), cmssw_base);
  }

  const baconhep::TTrigger triggerMenu(iHLTFile);  
  cout<<"Apply trigger: "<<applyTrigger<<endl;

  TLorentzVector W_type0,W_type0_jes_up, W_type0_jes_dn, W_type0_jer_up, W_type0_jer_dn, W_type2, W_run2,W_puppi_type2, W_puppi_type0, W_puppi_run2, W_type0_LEP_Up, W_type0_LEP_Down;
  TLorentzVector LEP1, LEP2, LEP1_Up, LEP1_Down, LEP2_Up, LEP2_Down, SJ1_PuppiAK8, SJ2_PuppiAK8, SJ1, SJ2;
  TLorentzVector NU0,NU1,NU2,NU0_puppi,NU1_puppi,NU2_puppi;
  TLorentzVector JET, JET_PuppiAK8, AK4;
  TLorentzVector AK4_JET1,AK4_JET2;
  TLorentzVector PuppiAK4_JET1,PuppiAK4_JET2;
  TLorentzVector VBF1,VBF2,TOT;
  TLorentzVector ELE,MU;

  
  vector<TLorentzVector> tightMuon;
  vector<TLorentzVector> looseMuon;
  vector<TLorentzVector> tightEle;
  vector<TLorentzVector> looseEle;

  int ok=0, total=0;
  
  // Data structures to store info from TTrees
  baconhep::TEventInfo *info  	= new baconhep::TEventInfo();
  baconhep::TGenEventInfo *gen	= new baconhep::TGenEventInfo();
  TClonesArray *genPartArr 	= new TClonesArray("baconhep::TGenParticle");
  TClonesArray *lheWgtArr	= new TClonesArray("baconhep::TLHEWeight");
  

  char command1[3000];
  char command2[3000];
  if ( cluster == "lxplus")
  	sprintf(command1, "eos find -f %s  | awk '!/log|fail/ {print $1}' | awk 'NF {print \"root://eoscms.cern.ch/\"$1}' > listTemp_%s.txt", (inputFolder).c_str(), outputFile.c_str());	// NF in awk command skips the blank line
  else 
	sprintf(command1,"xrdfs root://cmseos.fnal.gov ls %s | awk '{print \"root://cmseos.fnal.gov/\"$1}' > listTemp_%s.txt",(inputFolder).c_str(),  outputFile.c_str());
	//sprintf(command1,"eos root://cmseos.fnal.gov find -f %s | awk '!/log|fail/ {print $1}' | awk 'NF {print \"root://cmseos.fnal.gov/\"$1}' > listTemp_%s.txt",(inputFolder).c_str(),  outputFile.c_str());	// WORKS ONLY WITH INTERACTIVE NODE

  cout<<command1<<endl;
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
  int cutEff[21]={0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};

  ofstream myfile;
  ofstream myfile2;
  myfile.open ("Neutrino_pt_type02.txt");
  myfile2.open ("Neutrino_pt_run2.txt");
  myfile << "IsComplex \t Truth \t Sol1 \t Sol2 \t type0 \t type1 \t type2"<< endl; 
  myfile2 << "IsComplex \t Truth \t Sol1 \t Sol2 \t run2"<< endl; 
  //myfile << "Truth \t Sol1 \t Sol2 \t type0 \t type1 \t type2 \t run2" << endl; 




  //---------output tree----------------
  TFile* outROOT = TFile::Open((outputFile+(".root")).c_str(),"recreate");
  outROOT->cd();
  TTree* outTree = new TTree("otree", "otree");
  setOutputTree* WWTree = new setOutputTree(outTree);

  int nEvents=0;	
  Long64_t jentry2=0;
  int count_genEvents=0;

  int nInputFiles = sampleName.size();

  if (isLocal==1) nInputFiles = 51;
  cout<<"==> Total number of input files : "<<nInputFiles<<endl;

  
  Long64_t TotalNumberOfEvents = 0, nNegEvents = 0; 

  

  
  // Loop on input files
  for(int i=0;i<nInputFiles;i++)
  {
     infile = TFile::Open(sampleName[i]);
     eventTree = (TTree*)infile->Get("Events");
     //cout << "\t File no. " << i << "\t"<< sampleName[i] << endl;
     
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
	    if (jentry2%50000 == 0) cout << "\t File no. " << i << "; Neg Event Count; read entry: " << jentry2 <<"/"<<TotalNumberOfEvents<< endl;
	    if (gen->weight<0)	nNegEvents++;
	}
     }
     delete infile;
     infile=0, eventTree=0;
  }
  
  
  cout<<"==> Total number of events : "<<TotalNumberOfEvents<<endl;
  cout<<"==> Total number of negative events : "<<nNegEvents<<endl;

  float weight = atof(xSecWeight.c_str())/(atof(TotalNumberOfEntries.c_str()) - 2*atof(TotalNumberOfNegativeEntries.c_str()));
  cout<<"Weight of cross-sec/events = "<<weight<<endl;
  int totalEntries=0;



  //---------start loop on events------------
  cout << "---------start loop on events------------" << endl;
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

    

    if (jentry2%10000 == 0) cout << "\tread entry: " << jentry2 <<"/"<<TotalNumberOfEvents<< endl;
    
    //*********************************
    WWTree->initializeVariables(); //initialize all variables

    WWTree->run   = info->runNum;
    WWTree->event = info->evtNum;
    WWTree->lumi  = info->lumiSec;


    /////////////////MC Info
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
      vector<TLorentzVector> v_GEN_hadW, v_GEN_lepW, v_GEN_VBFJ1, v_GEN_VBFJ2, v_GEN_VBFJ, v_GEN_temp;
      vector<TLorentzVector> v_genLep, v_genNeutrino, v_genWquarks, v_genVBFquarks;

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
	      if (abs(genloop->pdgId) == 11) leptonName = "el";	      
	      if (abs(genloop->pdgId) == 13) leptonName = "mu";
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
	  WWTree->lep_phi_gen     = v_genLep[0].Phi();
	  WWTree->lep_mass_gen     = v_genLep[0].M();
	  WWTree->nu_pz_gen	  = v_genNeutrino[0].Pz();
	  WWTree->nu_pt_gen	  = v_genNeutrino[0].Pt();
	  WWTree->nu_eta_gen	  = v_genNeutrino[0].Eta();
	  WWTree->nu_phi_gen	  = v_genNeutrino[0].Phi();
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

    //////////////THE MET
    if (GenPassCut == 0) continue;
    //cout << "GenPassCut = " << GenPassCut << endl;
    //cout << "GenPassCut Passed............"  << endl;
    
    // //preselection on met
     if (WWTree->nu_pt_gen < 30) continue;
     //cout<< "Neu pt = " << WWTree->nu_pt_gen << endl;
    
    // Calculate Neutrino Pz using all the possible choices : 
    // type0 -> if real roots, pick the one nearest to the lepton Pz except when the Pz so chosen
    //               is greater than 300 GeV in which case pick the most central root.
    // type1 -> type = 1: if real roots, choose the one closest to the lepton Pz if complex roots, use only the real part.
    //          type = 2: if real roots, choose the most central solution. if complex roots, use only the real part. 
    //          type = 3: if real roots, pick the largest value of the cosine*
    
    LEP1.SetPtEtaPhiM(WWTree->lep_pt_gen,WWTree->lep_eta_gen,WWTree->lep_phi_gen,WWTree->lep_mass_gen);
    TLorentzVector W_Met;
    WWTree->pfMET_Corr = WWTree->nu_pt_gen;
    WWTree->pfMET_Corr_phi = WWTree->nu_phi_gen;

    //cout << "===> " <<  WWTree->pfMET_Corr  << endl;
    float Wmass = 80.385;
  
    TLorentzVector W_Met_jes_up, W_Met_jes_dn, W_Met_jer_up, W_Met_jer_dn, AK4Up, AK4Down, AK4Up_Puppi, AK4Down_Puppi;
  
    W_Met.SetPxPyPzE(WWTree->nu_pt_gen * TMath::Cos(WWTree->nu_phi_gen), WWTree->nu_pt_gen * TMath::Sin(WWTree->nu_phi_gen), 0., sqrt(WWTree->nu_pt_gen*WWTree->nu_pt_gen));


    if(LEP1.Pt()<=0 || W_Met.Pt() <= 0 ){ cerr<<" Negative Lepton - Neutrino Pt "<<endl; continue ; }


    // type0 calculation of neutrino pZ
    METzCalculator NeutrinoPz_type0;

    NeutrinoPz_type0.SetMET(W_Met);
    NeutrinoPz_type0.SetLepton(LEP1);
    NeutrinoPz_type0.SetLeptonType(leptonName.c_str());
    
    double pz1_type0 = NeutrinoPz_type0.Calculate(0); // Default one -> according to type0
    double pz2_type0 = NeutrinoPz_type0.getOther();  // Default one
    
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
      W_neutrino_1.SetPxPyPzE(nu_pt1 * TMath::Cos(WWTree->nu_phi_gen), nu_pt1 * TMath::Sin(WWTree->nu_phi_gen), pz1_type0, sqrt(nu_pt1*nu_pt1 + pz1_type0*pz1_type0) );
      TLorentzVector W_neutrino_2;
      W_neutrino_2.SetPxPyPzE(nu_pt2 * TMath::Cos(WWTree->nu_phi_gen), nu_pt2 * TMath::Sin(WWTree->nu_phi_gen), pz1_type0, sqrt(nu_pt2*nu_pt2 + pz1_type0*pz1_type0) );
      
      if ( fabs((LEP1+W_neutrino_1).M()-Wmass) < fabs((LEP1+W_neutrino_2).M()-Wmass) ) W_neutrino_type0 = W_neutrino_1;
      else W_neutrino_type0 = W_neutrino_2;
    }
    
    // type1 calculation of neutrino pZ
    METzCalculator NeutrinoPz_type1;

    NeutrinoPz_type1.SetMET(W_Met);
    NeutrinoPz_type1.SetLepton(LEP1);
    NeutrinoPz_type1.SetLeptonType(leptonName.c_str());
    
    double pz1_type1 = NeutrinoPz_type1.Calculate(1); // Default one -> according to type1
    double pz2_type1 = NeutrinoPz_type1.getOther();  // Default one
    
    // don't touch the neutrino pT
    TLorentzVector W_neutrino_type1_met; 
    W_neutrino_type1_met.SetPxPyPzE(W_Met.Px(),W_Met.Py(),pz1_type1,sqrt(W_Met.Pt()*W_Met.Pt()+pz1_type1*pz1_type1));
    
    // change the neutrino pT in case of complex solution in order to make it real
    TLorentzVector W_neutrino_type1; 
    W_neutrino_type1.SetPxPyPzE(W_Met.Px(),W_Met.Py(),pz1_type1,sqrt(W_Met.Pt()*W_Met.Pt()+pz1_type1*pz1_type1));
    
    if(NeutrinoPz_type1.IsComplex()) {// if this is a complex, change MET
      double nu_pt1 = NeutrinoPz_type1.getPtneutrino(1);
      double nu_pt2 = NeutrinoPz_type1.getPtneutrino(2);
      TLorentzVector W_neutrino_1;
      W_neutrino_1.SetPxPyPzE(nu_pt1 * TMath::Cos(WWTree->nu_phi_gen), nu_pt1 * TMath::Sin(WWTree->nu_phi_gen), pz1_type1, sqrt(nu_pt1*nu_pt1 + pz1_type1*pz1_type1) );
      TLorentzVector W_neutrino_2;
      W_neutrino_2.SetPxPyPzE(nu_pt2 * TMath::Cos(WWTree->nu_phi_gen), nu_pt2 * TMath::Sin(WWTree->nu_phi_gen), pz1_type1, sqrt(nu_pt2*nu_pt2 + pz1_type1*pz1_type1) );
      
      if ( fabs((LEP1+W_neutrino_1).M()-Wmass) < fabs((LEP1+W_neutrino_2).M()-Wmass) ) W_neutrino_type1 = W_neutrino_1;
      else W_neutrino_type1 = W_neutrino_2;
    }
    
    // type2 calculation of neutrino pZ
    METzCalculator NeutrinoPz_type2;
    NeutrinoPz_type2.SetMET(W_Met);
    NeutrinoPz_type2.SetLepton(LEP1);
    NeutrinoPz_type2.SetLeptonType(leptonName.c_str());
    
    double pz1_type2 = NeutrinoPz_type2.Calculate(2); // Default one -> according to type2
    double pz2_type2 = NeutrinoPz_type2.getOther();   // Default one
  
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
      W_neutrino_1.SetPxPyPzE(nu_pt1 * TMath::Cos(WWTree->nu_phi_gen), nu_pt1 * TMath::Sin(WWTree->nu_phi_gen), pz1_type2, sqrt(nu_pt1*nu_pt1 + pz1_type2*pz1_type2) );
      TLorentzVector W_neutrino_2;
      W_neutrino_2.SetPxPyPzE(nu_pt2 * TMath::Cos(WWTree->nu_phi_gen), nu_pt2 * TMath::Sin(WWTree->nu_phi_gen), pz1_type2, sqrt(nu_pt2*nu_pt2 + pz1_type2*pz1_type2) );
      
      if ( fabs((LEP1+W_neutrino_1).M()-Wmass) < fabs((LEP1+W_neutrino_2).M()-Wmass) ) W_neutrino_type2 = W_neutrino_1;
      else W_neutrino_type2 = W_neutrino_2;
    }
    
    // run2 calculation of neutrino pZ
    METzCalculator_Run2 NeutrinoPz_run2;

    NeutrinoPz_run2.SetMET(W_Met);
    NeutrinoPz_run2.SetLepton(LEP1);
    NeutrinoPz_run2.SetLeptonType(leptonName.c_str());

    double pz1_run2 = NeutrinoPz_run2.Calculate();
    double pz2_run2 = NeutrinoPz_run2.getOther();

    
    /////////
    WWTree->pfMET = sqrt(info->pfMET*info->pfMET);
    WWTree->pfMET_Phi = info->pfMETphi;
    WWTree->pfMET_Corr = WWTree->nu_pt_gen;
    WWTree->pfMET_Corr_phi = WWTree->nu_phi_gen;
    WWTree->pfMET_Corr_Cov00 = info->pfMETCCov00;
    WWTree->pfMET_Corr_Cov01 = info->pfMETCCov01;
    WWTree->pfMET_Corr_Cov11 = info->pfMETCCov11;

    WWTree->nu_pz_type0 = pz1_type0;
    //WWTree->nu_pz_type1 = pz1_type1;
    WWTree->nu_pz_type2 = pz1_type2;


    //cout<< WWTree->nu_pz_gen << "\t" << WWTree->nu_pz_type0 << "\t" <<  WWTree->nu_pz_type2 << endl; 
    //cout<< WWTree->nu_pz_gen << "\t" << WWTree->nu_pz_type0 << "\t" << pz2_type0 << "\t" << WWTree->nu_pz_type2 << "\t" << pz2_type2 << endl;
    //myfile <<  WWTree->nu_pz_gen << "\t" << WWTree->nu_pz_type0 << "\t" << pz2_type0 << "\t" << WWTree->nu_pz_type0 << "\t" << pz1_type1 << "\t" << WWTree->nu_pz_type2 << endl;
    myfile <<  NeutrinoPz_type0.IsComplex() << "\t" << WWTree->nu_pz_gen << "\t" << WWTree->nu_pz_type0 << "\t" << pz2_type0 << "\t" << WWTree->nu_pz_type0 << "\t" << pz1_type1 << "\t" << WWTree->nu_pz_type2 << endl;
    myfile2 << NeutrinoPz_run2.IsComplex() << "\t" << WWTree->nu_pz_gen << "\t" << pz1_run2 << "\t" << pz2_run2 << "\t" << pz1_run2 << endl;
  
    
    /////////////////THE LEPTONIC W
  
    NU0.SetPtEtaPhiM( GetPt_MET(WWTree->nu_pt_gen, WWTree->nu_phi_gen, WWTree->nu_pz_type0) , GetEta_MET(WWTree->nu_pt_gen, WWTree->nu_phi_gen, WWTree->nu_pz_type0), WWTree->nu_phi_gen , 0.0 );



    
    NU2.SetPtEtaPhiM(GetPt_MET( WWTree->nu_pt_gen, WWTree->nu_phi_gen, WWTree->nu_pz_type2 ), GetEta_MET( WWTree->nu_pt_gen, WWTree->nu_phi_gen, WWTree->nu_pz_type2 ), WWTree->nu_phi_gen, 0.0);

  
    W_type0 = LEP1 + NU0;
    W_type2 = LEP1 + NU2;
    

  
    WWTree->v_pt_type0 = W_type0.Pt();
    WWTree->v_eta_type0 = W_type0.Eta();
    WWTree->v_mt_type0 = TMath::Sqrt(2*LEP1.Et()*NU0.Et()*(1-TMath::Cos(LEP1.DeltaPhi(NU0))));
    WWTree->v_mass_type0 = W_type0.M();

    WWTree->v_pt_type2 = W_type2.Pt();
    WWTree->v_eta_type2 = W_type2.Eta();
    WWTree->v_mass_type2 = W_type2.M();
    //if(LEP1.Pt()<=0 || LEP2.Pt() <= 0 ){ cerr<<" Negative Lepton - Neutrino Pt "<<endl; continue ; }
  

    outTree->Fill();
    }

    delete infile;
    infile=0, eventTree=0;
    /////////////////FILL THE TREE
  }
  myfile.close();
  myfile2.close();
  //delete pileupHisto;
  //pileupFile->Close();
  //file->Close();
  cout << "---------end loop on events------------" << endl;
  cout << endl;
  cout << "GEN events = " << count_genEvents << endl;


  
  cout << "----------------------" << endl;
  cout << " SUMMARY" << endl;
  cout << "----------------------" << endl;
  cout << endl;
  cout<<"MC matching: "<<(float)ok/(float)total<<endl;
  cout<<"negative events: "<<nNegEvents<<endl;
  cout << endl;
  cout<<"(0) all events:        "<<cutEff[0]<<"\t:\t"<<((float)cutEff[0]*100.0)/(float)cutEff[0]<<endl
  	   <<"(1) Gen Events:        "<<cutEff[1]<<"\t:\t"<<((float)cutEff[1]*100.0)/(float)cutEff[0]<<endl
	   <<"(2) tight lepton:      "<<cutEff[2]<<"\t:\t"<<((float)cutEff[2]*100.0)/(float)cutEff[0]<<endl
	   <<"(3) MET:               "<<cutEff[3]<<"\t:\t"<<((float)cutEff[3]*100.0)/(float)cutEff[2]<<endl
	   <<"(4) negative lep-MET:  "<<cutEff[4]<<"\t:\t"<<((float)cutEff[4]*100.0)/(float)cutEff[3]<<endl
	   <<"(5) 1 good AK8:        "<<cutEff[5]<<"\t:\t"<<((float)cutEff[5]*100.0)/(float)cutEff[4]<<endl
	   <<"(6) m(WV) > 0:       "<<cutEff[6]<<"\t:\t"<<((float)cutEff[6]*100.0)/(float)cutEff[5]<<endl
	   <<"(7) >=2 good VBF jets: "<<cutEff[7]<<"\t:\t"<<((float)cutEff[7]*100.0)/(float)cutEff[6]<<endl
	   <<"(8) Found VBF jets:  "<<cutEff[8]<<"\t:\t"<<((float)cutEff[8]*100.)/(float)cutEff[7]<<endl;
  
  cout<<"\n\n----------------------------"<< endl;
  cout<<"\tSUMMARY for 2 lepton case "<< endl;
  cout<<"---------------------------"<< endl;
  cout << endl;
  cout<<"(0) all events:        "<<cutEff[0]<<"\t:\t"<<((float)cutEff[0]*100.0)/(float)cutEff[0]<<endl
  	   <<"(1) Gen Events:        "<<cutEff[1]<<"\t:\t"<<((float)cutEff[1]*100.0)/(float)cutEff[0]<<endl
	   <<"(2) tight lepton:      "<<cutEff[12]<<"\t:\t"<<((float)cutEff[12]*100.0)/(float)cutEff[0]<<endl
	   <<"(3) MET:               "<<cutEff[13]<<"\t:\t"<<((float)cutEff[13]*100.0)/(float)cutEff[12]<<endl
	   <<"(4) negative lep-MET:  "<<cutEff[14]<<"\t:\t"<<((float)cutEff[14]*100.0)/(float)cutEff[13]<<endl
	   <<"(5) 1 good AK8:        "<<cutEff[15]<<"\t:\t"<<((float)cutEff[15]*100.0)/(float)cutEff[14]<<endl
	   <<"(6) m(WV) > 0:       "<<cutEff[16]<<"\t:\t"<<((float)cutEff[16]*100.0)/(float)cutEff[15]<<endl
	   <<"(7) >=2 good VBF jets: "<<cutEff[17]<<"\t:\t"<<((float)cutEff[17]*100.0)/(float)cutEff[16]<<endl
	   <<"(8) Found VBF jets:  "<<cutEff[18]<<"\t:\t"<<((float)cutEff[18]*100.)/(float)cutEff[17]<<endl;
  
 
 
  //--------close everything-------------
  delete info; delete gen;
  delete genPartArr; 
  delete lheWgtArr;
  outROOT->Write();
  outROOT->Close();
  int t1 = time(NULL);
  printf ("\n==> time to run this code = %0.3f min\n", (float)(t1 - t0)/60.0);
  return(0);
}
