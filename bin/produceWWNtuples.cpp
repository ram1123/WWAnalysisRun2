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
#include "TCanvas.h"
#include "TApplication.h"
#include "TLorentzVector.h"

#include "../interface/setInputTree.h"
#include "../interface/setOutputTree.h"
#include "../interface/METzCalculator.h"
#include "../interface/METzCalculator_Run2.h"
#include "../interface/analysisUtils.h"

using namespace std;

//*****PU WEIGHT***************

vector<double> generate_weights(TH1* data_npu_estimated){
  // see SimGeneral/MixingModule/python/mix_E7TeV_FlatDist10_2011EarlyData_inTimeOnly_cfi.py; copy and paste from there:
  const double npu_probs[60] = {
   2.560E-06,
 5.239E-06,
 1.420E-05,
 5.005E-05,
 1.001E-04,
 2.705E-04,
 1.999E-03,
 6.097E-03,
 1.046E-02,
 1.383E-02,
 1.685E-02,
 2.055E-02,
 2.572E-02,
 3.262E-02,
 4.121E-02,
 4.977E-02,
 5.539E-02,
 5.725E-02,
 5.607E-02,
 5.312E-02,
 5.008E-02,
 4.763E-02,
 4.558E-02,
 4.363E-02,
 4.159E-02,
 3.933E-02,
 3.681E-02,
 3.406E-02,
 3.116E-02,
 2.818E-02,
 2.519E-02,
 2.226E-02,
 1.946E-02,
 1.682E-02,
 1.437E-02,
 1.215E-02,
 1.016E-02,
 8.400E-03,
 6.873E-03,
 5.564E-03,
 4.457E-03,
 3.533E-03,
 2.772E-03,
 2.154E-03,
 1.656E-03,
 1.261E-03,
 9.513E-04,
 7.107E-04,
 5.259E-04,
 3.856E-04,
 2.801E-04,
 2.017E-04,
 1.439E-04,
 1.017E-04,
 7.126E-05,
 4.948E-05,
 3.405E-05,
 2.322E-05,
 1.570E-05,
 5.005E-06
};
  vector<double> result(60);
  /*  
  double s = 0.0;
  for(int npu=0; npu<60; ++npu){
    double npu_estimated = data_npu_estimated->GetBinContent(data_npu_estimated->GetXaxis()->FindBin(npu));                              
    result[npu] = npu_estimated / npu_probs[npu];
    s += npu_estimated;
  }
  // normalize weights such that the total sum of weights over thw whole sample is 1.0, i.e., sum_i  result[i] * npu_probs[i] should be 1.0 (!)
  for(int npu=0; npu<60; ++npu){
    result[npu] /= s;
    }*/
  
  for(int npu=0; npu<60; ++npu){
    if (data_npu_estimated->GetBinContent(data_npu_estimated->GetXaxis()->FindBin(npu))==NULL)
      result[npu] = 0.;
    else {
      double npu_estimated = data_npu_estimated->GetBinContent(data_npu_estimated->GetXaxis()->FindBin(npu));                              
      result[npu] = npu_estimated;
    }
  }
  
  return result;
}


//*******MAIN*******************************************************************

int main (int argc, char** argv)
{ 
  std::string inputFolder = argv[1];
  std::string outputFile = argv[2];
  bool isMC = argv[3];
  std::string leptonName = argv[4];
  std::string inputTreeName = argv[5];
  std::string inputFile = argv[6];
  std::string xSecWeight = argv[7];
  std::string numberOfEntries = argv[8];
  float weight = std::atof(xSecWeight.c_str())/std::atof(numberOfEntries.c_str());
  if (strcmp(leptonName.c_str(),"el")!=0 && strcmp(leptonName.c_str(),"mu")!=0) {
    std::cout<<"Error: wrong lepton category"<<std::endl;
    return(-1);
  }
  float genMass = atof(argv[9]);
  int applyTrigger = atoi(argv[10]);
  //applyTrigger=false;
  std::cout<<"apply trigger: "<<applyTrigger<<std::endl;

bool verbose = 0;
  TLorentzVector W,MET,LEP;
  TLorentzVector NU0,NU1,NU2;
  TLorentzVector JET, HADW, AK4, AK4_new;
  TLorentzVector VBF1,VBF2,TOT, VBF1_AK4, VBF2_AK4, Wjet1_AK4, Wjet2_AK4, TOT_Wjet;
  TLorentzVector ELE,MU;

  std::vector<TLorentzVector> tightMuon;
  std::vector<TLorentzVector> looseMuon;
  std::vector<TLorentzVector> tightEle;
  std::vector<TLorentzVector> looseEle;

  int ok=0, total=0;
  int evento=127270357;
  int count=0;
	
  setInputTree *ReducedTree = new setInputTree (inputTreeName.c_str());
  ReducedTree->Init();


  char command1[3000];
  //sprintf(command1, "ls  %s/%s/  | awk '{print \"/eos/uscms/store/user/rasharma/WWScattering/26August15/ReducedTrees/%s/\"$1}' > listTemp_%s.txt", (inputFolder).c_str(), (inputFile).c_str(), (inputFile).c_str(), outputFile.c_str());
  //sprintf(command1, "xrd eoscms dirlist %s/%s/  | awk '{print \"root://xrootd.unl.edu/\"$5}' > listTemp_%s.txt", (inputFolder).c_str(), (inputFile).c_str(), outputFile.c_str());
  //sprintf(command1, "xrd eoscms dirlist %s/%s/  | awk '{print \"root://eoscms.cern.ch/\"$5}' > listTemp_%s.txt", (inputFolder).c_str(), (inputFile).c_str(), outputFile.c_str());
  //std::cout<<command1<<std::endl;
  //system(command1);
  char list1[2000];
  sprintf (list1, "InputRootFiles/listTemp_%s.txt", inputFile.c_str());
  //sprintf (list1, "listTemp_%s.txt", outputFile.c_str());
  //ifstream rootList ("test.txt");
  //ifstream rootList ("listFiles_Wjet.txt");
  ifstream rootList (list1);

  int fileCounter=0;
  Long64_t totalEntries=0;
#if 0
  while (!rootList.eof())
    {
      char iRun_tW[700];
      rootList >> iRun_tW;
      ReducedTree->fChain->Add(iRun_tW);
      cout<<"file counter = "<<fileCounter<<endl;
      fileCounter++;
    }
#else
	//ReducedTree->fChain->Add("/afs/cern.ch/user/r/rasharma/work/WW_Scattering/AnalysisFrameWork/CMSSW_7_4_7_patch2/src/AllHadronicSUSY/ReducedSelection.root");
	//ReducedTree->fChain->Add("/afs/cern.ch/user/r/rasharma/work/public/Temp/W_LW_L_50k.root");
	while (rootList.good())
	{
	string line;
	getline(rootList,line);
	char iRun_tW[700];
	strcpy(iRun_tW, line.c_str());
	ReducedTree->fChain->Add(iRun_tW);
	cout<<"file counter = "<<fileCounter<<endl;
	fileCounter++;
	}
#endif
  std::cout<<"number of files found: "<<fileCounter-1<<std::endl;
  std::cout<<"total entries: "<<ReducedTree->fChain->GetEntries()<<std::endl;
  totalEntries=ReducedTree->fChain->GetEntries();

  char command3[300];
  //sprintf(command3, "rm listTemp_%s.txt", outputFile.c_str());
  //system(command3);

  int cutEff[20]={0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};

  //--------pile up file -----------------
  //    TFile* pileupFile = TFile::Open("190456-208686-13Julv2_Prompt_Moriond2013.69400.observed.root");  
  //TH1F *pileupHisto = (TH1F*)pileupFile->Get("pileup");
  TFile* pileupFile = TFile::Open("PU.root");  
  TH1F *pileupHisto = (TH1F*)pileupFile->Get("puweights");

  std::vector<double> weights;

  weights = generate_weights(pileupHisto);
  pileupFile->Close();

  //---------output tree----------------
  TFile* outROOT = TFile::Open((std::string("output/output_")+leptonName+std::string("/")+outputFile+(".root")).c_str(),"recreate");
  outROOT->cd();
  TTree* outTree = new TTree("otree", "otree");
  outTree->SetDirectory(0);

  setOutputTree *WWTree = new setOutputTree(outTree);

  //---------start loop on events------------
  Long64_t jentry2=0;

  int TotalEle= 0, TriggerPassEle = 0, MediumPassEle = 0, PtPassEle = 0, EtaPassEle = 0;
  int TotalMu= 0, TriggerPassMu = 0, TightPassMu = 0, PtPassMu = 0, EtaPassMu = 0, IsoPassMu = 0;
  int TotalAK4Jets = 0, TotalAK4Jets_MoreThan4 = 0;
  int IsJet = 0, JetPtEtaPass = 0, JetLoosePass = 0, JetLepCleanPass = 0;

  int BoolTotalEle = 0, BoolTriggerPassEle = 0, BoolMediumPassEle = 0, BoolPtPassEle = 0, BoolEtaPassEle = 0;
  int BoolTotalMu = 0, BoolTriggerPassMu = 0, BoolTightPassMu = 0, BoolPtPassMu = 0, BoolEtaPassMu = 0, BoolIsoPassMu = 0;
  int BoolTotalAK4Jets = 0, BoolTotalAK4Jets_MoreThan4 = 0;
  int BoolIsJet = 0, BoolJetPtEtaPass = 0, BoolJetLoosePass = 0, BoolJetLepCleanPass = 0;

TH1F * ptEle = new TH1F("ptEle","",100,0,300);
TH1F * ptMet = new TH1F("ptMet","",100,0,300);
TH1F * pt_W = new TH1F("pt_W","",100,0,300);
TH1F * mjj_01 = new TH1F("mjj_01","",100,0,2500);
TH1F * mjj_02 = new TH1F("mjj_02","",100,0,2500);
TH1F * mjj_03 = new TH1F("mjj_03","",100,0,2500);
TH1F * mjj_12 = new TH1F("mjj_12","",100,0,2500);
TH1F * mjj_13 = new TH1F("mjj_13","",100,0,2500);
TH1F * mjj_23 = new TH1F("mjj_23","",100,0,2500);
//TCanvas *c1 = new TCanvas("c1","",1);

  for (Long64_t jentry=0; jentry<ReducedTree->fChain->GetEntries();jentry++,jentry2++) {
  //for (Long64_t jentry=0; jentry<50000;jentry++,jentry2++) {

    Long64_t iEntry = ReducedTree->LoadTree(jentry);
    if (iEntry < 0) break;
    int nb = ReducedTree->fChain->GetEntry(jentry);   
    // if (Cut(ientry) < 0) continue;                                                                                                                           

    tightMuon.clear();
    tightEle.clear();
    looseMuon.clear();
    looseEle.clear();

    if(jentry2 % 1000 == 0)    
      std::cout << "read entry: " << jentry2 <<"/"<<totalEntries<<std:: endl;

    WWTree->initializeVariables(); //initialize all variables
    
    WWTree->issignal = 0;
    WWTree->wSampleWeight = weight; //xsec/numberOfEntries
    WWTree->totalEventWeight = 1.; //temporary value
    WWTree->eff_and_pu_Weight = 1.; //temporary value

    if (ReducedTree->genEventWeight>0)
      WWTree->genWeight=1.;
    else if (ReducedTree->genEventWeight<0)
      WWTree->genWeight=-1.;
    //    WWTree->genWeight = ReducedTree->genEventWeight;

    //PILE-UP WEIGHT
    if (isMC) {
      if(ReducedTree->NVtx<weights.size()){
	WWTree->eff_and_pu_Weight = weights[ReducedTree->NVtx];
	WWTree->totalEventWeight*=weights[ReducedTree->NVtx];
      }
      else{ //should not happen as we have a weight for all simulated n_pu multiplicities!
	std::cout<<"Warning! n_pu too big"<<std::endl;
	//	throw logic_error("n_pu too big");
        WWTree->eff_and_pu_Weight = 0.;
	WWTree->totalEventWeight*=0.;
      }    
    }    
    //require at least one lepton and one jet
    //    if ( strcmp(leptonName.c_str(),"el")==0 && ReducedTree->ElectronsNum==0) continue; 
    //    if ( strcmp(leptonName.c_str(),"mu")==0 && ReducedTree->MuonsNum==0) continue;      
        
    //save event variables
    WWTree->run   = ReducedTree->RunNum;
    WWTree->event = ReducedTree->EvtNum;
    WWTree->lumi = ReducedTree->LumiBlockNum;
   // WWTree->njets = ReducedTree->NJets;
    WWTree->nPV  = ReducedTree->NVtx;

    if(WWTree->event==evento) std::cout<<"debug: "<<count<<std::endl; count++;

    if (verbose)
    	cout<<"==================> debug 1 "<<endl;

    /////////////////THE SELECTED LEPTON
    int nTightLepton=0;
   BoolTotalEle = 0, BoolTriggerPassEle = 0, BoolMediumPassEle = 0, BoolPtPassEle = 0, BoolEtaPassEle = 0;
   BoolTotalMu = 0, BoolTriggerPassMu = 0, BoolTightPassMu = 0, BoolPtPassMu = 0, BoolEtaPassMu = 0, BoolIsoPassMu = 0;
   BoolTotalAK4Jets = 0, BoolTotalAK4Jets_MoreThan4 = 0;
    if (strcmp(leptonName.c_str(),"el")==0) {
    	if (verbose)
    	cout<<"==================> debug 2 "<<endl;
      float tempPt=0.;


      for (int i=0; i<ReducedTree->ElectronsNum; i++) {
	//if (applyTrigger==1 && ReducedTree->TriggerProducerTriggerPass->at(0)==0) continue; //trigger
	//if (ReducedTree->TriggerProducerTriggerPass->at(0)==0) continue; //trigger
	BoolTotalEle = 1;
	//if (ReducedTree->Electrons_isHEEP[i]==false) continue;       
	if (ReducedTree->Electrons_isMedium[i]==false) continue;       
	ptEle->Fill(ReducedTree->ElectronsPt[i]);
	//if (ReducedTree->Electrons_isLoose[i]==false) continue;       
	//if (ReducedTree->Electrons_isTight[i]==false) continue;       
	BoolMediumPassEle = 1;
        if (ReducedTree->ElectronsPt[i]<=25) continue;
	BoolPtPassEle = 1;
        if (fabs(ReducedTree->ElectronsEta[i])>=2.5) continue;
	BoolEtaPassEle = 1;
	if (ReducedTree->ElectronsPt[i]<tempPt) continue;
	ELE.SetPtEtaPhiE(ReducedTree->ElectronsPt[i],ReducedTree->ElectronsEta[i],ReducedTree->ElectronsPhi[i],ReducedTree->ElectronsE[i]);
	tightEle.push_back(ELE);
	WWTree->l_pt  = ReducedTree->ElectronsPt[i];
	WWTree->l_eta = ReducedTree->ElectronsEta[i];
	WWTree->l_phi = ReducedTree->ElectronsPhi[i];	
	WWTree->l_e= ReducedTree->ElectronsE[i];	
	tempPt = WWTree->l_pt;
	nTightLepton++;
      }
    }
    else if (strcmp(leptonName.c_str(),"mu")==0) {
    	if (verbose)
    	cout<<"==================> debug 3 "<<endl;
      float tempPt=0.;
      for (int i=0; i<ReducedTree->MuonsNum; i++) {
      BoolTotalMu = 1;
	//if (applyTrigger==1 && ReducedTree->TriggerProducerTriggerPass->at(0)==0) continue; //trigger
	//if (ReducedTree->TriggerProducerTriggerPass->at(1)==0) continue; //trigger
	//if (ReducedTree->Muons_isLoose[i]==false) continue;
	if (ReducedTree->Muons_isTight[i]==false) continue;
	//if (ReducedTree->Muons_isHighPt[i]==false) continue;
	ptEle->Fill(ReducedTree->ElectronsPt[i]);
	BoolTightPassMu = 1;
	//	if (ReducedTree->Muons_isPFMuon[i]==false) continue; //not in the synch ntuple!!
        if ((ReducedTree->Muons_trackIso[i]/ReducedTree->MuonsPt[i])>=0.1) continue;
	BoolIsoPassMu = 1;
        if (ReducedTree->MuonsPt[i]<25) continue;
	BoolPtPassMu = 1;
        if (fabs(ReducedTree->MuonsEta[i])>=2.4) continue;
	BoolEtaPassMu = 1;
	MU.SetPtEtaPhiE(ReducedTree->MuonsPt[i],ReducedTree->MuonsEta[i],ReducedTree->MuonsPhi[i],ReducedTree->MuonsE[i]);
	tightMuon.push_back(MU);
	if (ReducedTree->MuonsPt[i]<tempPt) continue;
	WWTree->l_pt  = ReducedTree->MuonsPt[i];
	WWTree->l_eta = ReducedTree->MuonsEta[i];
	WWTree->l_phi = ReducedTree->MuonsPhi[i];
	WWTree->l_e = ReducedTree->MuonsE[i];
	tempPt = WWTree->l_pt;
	nTightLepton++;
      }
    	if (verbose)
    	cout<<"==================> debug 9 "<<endl;
    }
    //======================= START::Counting events that passed above lep ids	=======================

	
	if(BoolEtaPassMu)
	EtaPassMu++;
	if(BoolPtPassMu)
	PtPassMu++;
	if(BoolIsoPassMu)
	IsoPassMu++;
	if(BoolTightPassMu)
	TightPassMu++;
//	if(BoolTriggerPassMu)
//	TriggerPassMu++;
        if(BoolTotalMu)
        TotalMu++;
	if(BoolEtaPassEle)
	EtaPassEle++;
	if(BoolPtPassEle)
	PtPassEle++;
	if(BoolMediumPassEle)
	MediumPassEle++;
        if(BoolTotalEle)
        TotalEle++;
//	if(BoolTriggerPassEle)
//	TriggerPassEle++;

    //======================= END::Counting events that passed above lep ids	=======================
    if (nTightLepton==0) continue; //no leptons with required ID
    if(WWTree->event==evento) std::cout<<"debug: "<<count<<std::endl; count++;

    //VETO ADDITIONAL LEPTONS
    int nLooseLepton=0;
    for (int i=0; i<ReducedTree->ElectronsNum; i++) {
    if(WWTree->event==evento) std::cout<<"debug: "<<i<<std::endl; count++;
      if (ReducedTree->Electrons_isLoose[i]==false) continue;       
      //if (ReducedTree->Electrons_isHEEP[i]==false) continue;       
    if(WWTree->event==evento) std::cout<<"debug: "<<i<<std::endl; count++;
      if (ReducedTree->ElectronsPt[i]<15) continue;       
      if (fabs(ReducedTree->ElectronsEta[i])>=2.5) continue;       
    if(WWTree->event==evento) std::cout<<"debug: "<<i<<std::endl; count++;
      ELE.SetPtEtaPhiE(ReducedTree->ElectronsPt[i],ReducedTree->ElectronsEta[i],ReducedTree->ElectronsPhi[i],ReducedTree->ElectronsE[i]);
      looseEle.push_back(ELE);      
      nLooseLepton++;
    }
    for (int i=0; i<ReducedTree->MuonsNum; i++) {
    if(WWTree->event==evento) std::cout<<"debug: "<<i<<std::endl; count++;
      if (ReducedTree->Muons_isLoose[i]==false) continue;
      //if (ReducedTree->Muons_isHighPt[i]==false) continue;
    if(WWTree->event==evento) std::cout<<"debug: "<<i<<std::endl; count++;
      if ((ReducedTree->Muons_trackIso[i]/ReducedTree->MuonsPt[i])>=0.1) continue;
    if(WWTree->event==evento) std::cout<<"debug: "<<i<<std::endl; count++;
      if (fabs(ReducedTree->MuonsEta[i])>=2.4) continue;
    if(WWTree->event==evento) std::cout<<"debug: "<<i<<std::endl; count++;
      if (ReducedTree->MuonsPt[i]<15) continue;
    if(WWTree->event==evento) std::cout<<"debug: "<<i<<std::endl; count++;
      MU.SetPtEtaPhiE(ReducedTree->MuonsPt[i],ReducedTree->MuonsEta[i],ReducedTree->MuonsPhi[i],ReducedTree->MuonsE[i]);
      looseMuon.push_back(MU);
      nLooseLepton++;
    }
    if(WWTree->event==evento)     std::cout<<nLooseLepton<<std::endl;
    if (nLooseLepton!=1) continue; //no additional leptons
    cutEff[0]++;
    ptMet->Fill(ReducedTree->METPt);
    if(WWTree->event==evento) std::cout<<"debug: "<<count<<std::endl; count++;

    //preselection on jet pt and met
    WWTree->Met_pt  = ReducedTree->METPt;
    if (ReducedTree->METPt < 20) continue; 
    cutEff[1]++;
    if(WWTree->event==evento) std::cout<<"debug: "<<count<<std::endl; count++;

    MET.SetPtEtaPhiE(ReducedTree->METPt,0.,ReducedTree->METPhi,0.);
    LEP.SetPtEtaPhiE(WWTree->l_pt,WWTree->l_eta,WWTree->l_phi,WWTree->l_e);

    //    if (ReducedTree->AK8Jets_PtCorr[0] < 150) continue;     
    
    //lepton Pt preselection
    //    if ( strcmp(leptonName.c_str(),"el")==0 && ReducedTree->ElectronsPt[0]<30) continue; 
    //    if ( strcmp(leptonName.c_str(),"mu")==0 && ReducedTree->MuonsPt[0]<30) continue; 

    //////////////THE MET

    // Calculate Neutrino Pz using all the possible choices : 
    // type0 -> if real roots, pick the one nearest to the lepton Pz except when the Pz so chosen
    //               is greater than 300 GeV in which case pick the most central root.
    // type1 -> type = 1: if real roots, choose the one closest to the lepton Pz if complex roots, use only the real part.
    //          type = 2: if real roots, choose the most central solution. if complex roots, use only the real part. 
    //          type = 3: if real roots, pick the largest value of the cosine*

    float Wmass = 80.385;

    TLorentzVector W_mu, W_Met;

    W_mu.SetPtEtaPhiE(WWTree->l_pt,WWTree->l_eta,WWTree->l_phi,WWTree->l_e);
    W_Met.SetPxPyPzE(ReducedTree->METPt * TMath::Cos(ReducedTree->METPhi), ReducedTree->METPt * TMath::Sin(ReducedTree->METPhi), 0., sqrt(ReducedTree->METPt*ReducedTree->METPt));

    if(W_mu.Pt()<=0 || W_Met.Pt() <= 0 ){ std::cerr<<" Negative Lepton - Neutrino Pt "<<std::endl; continue ; }
    if(WWTree->event==evento) std::cout<<"debug: "<<count<<std::endl; count++;

    // type0 calculation of neutrino pZ
    METzCalculator NeutrinoPz_type0;
    METzCalculator_Run2 NeutrinoPz_run2;
    NeutrinoPz_type0.SetMET(W_Met);
    NeutrinoPz_type0.SetLepton(W_mu);
    NeutrinoPz_type0.SetLeptonType(leptonName.c_str());

    NeutrinoPz_run2.SetMET(W_Met);
    NeutrinoPz_run2.SetLepton(W_mu);
    NeutrinoPz_run2.SetLeptonType(leptonName.c_str());

    double pz1_type0 = NeutrinoPz_type0.Calculate(); // Default one -> according to type0
    double pz2_type0 = NeutrinoPz_type0.getOther(); // Default one

    double pz1_run2 = NeutrinoPz_run2.Calculate(); 

    // don't touch the neutrino pT
    TLorentzVector W_neutrino_type0_met; 
    W_neutrino_type0_met.SetPxPyPzE(W_Met.Px(),W_Met.Py(),pz1_type0,sqrt(W_Met.Pt()*W_Met.Pt()+pz1_type0*pz1_type0));

    //    W_mass_type0_met = (W_neutrino_type0_met+W_mu).M();
    //    W_pz_type0_met = (W_neutrino_type0_met+W_mu).Pz();
    //    W_nu1_pz_type0_met = pz1_type0;
    //    W_nu2_pz_type0_met = pz2_type0;

    // chenge the neutrino pT in case of complex solution in order to make it real
    TLorentzVector W_neutrino_type0; 
    W_neutrino_type0.SetPxPyPzE(W_Met.Px(),W_Met.Py(),pz1_type0,sqrt(W_Met.Pt()*W_Met.Pt()+pz1_type0*pz1_type0));

    if (NeutrinoPz_type0.IsComplex()) {// if this is a complix, change MET
      double nu_pt1 = NeutrinoPz_type0.getPtneutrino(1);
      double nu_pt2 = NeutrinoPz_type0.getPtneutrino(2);
      TLorentzVector W_neutrino_1;
      W_neutrino_1.SetPxPyPzE(nu_pt1 * TMath::Cos(ReducedTree->METPhi),
			      nu_pt1 * TMath::Sin(ReducedTree->METPhi), pz1_type0, sqrt(nu_pt1*nu_pt1 + pz1_type0*pz1_type0) );
      TLorentzVector W_neutrino_2;
      W_neutrino_2.SetPxPyPzE(nu_pt2 * TMath::Cos(ReducedTree->METPhi),
			      nu_pt2 * TMath::Sin(ReducedTree->METPhi), pz1_type0, sqrt(nu_pt2*nu_pt2 + pz1_type0*pz1_type0) );

      if ( fabs((W_mu+W_neutrino_1).M()-Wmass) < fabs((W_mu+W_neutrino_2).M()-Wmass) ) W_neutrino_type0 = W_neutrino_1;
      else W_neutrino_type0 = W_neutrino_2;
    }

    //    W_mass_type0 = (W_mu+W_neutrino_type0).M();
    //    W_pz_type0 = (W_mu+W_neutrino_type0).Pz();
    //    W_nu1_pz_type0 = pz1_type0;
    //    W_nu2_pz_type0 = pz2_type0;

    // type2 calculation of neutrino pZ
    METzCalculator NeutrinoPz_type2;
    NeutrinoPz_type2.SetMET(W_Met);
    NeutrinoPz_type2.SetLepton(W_mu);
    NeutrinoPz_type2.SetLeptonType(leptonName.c_str());

    double pz1_type2 = NeutrinoPz_type2.Calculate(2); // Default one -> according to type2
    double pz2_type2 = NeutrinoPz_type2.getOther(); // Default one

    // don't touch the neutrino pT
    TLorentzVector W_neutrino_type2_met; 
    W_neutrino_type2_met.SetPxPyPzE(W_Met.Px(),W_Met.Py(),pz1_type2,sqrt(W_Met.Pt()*W_Met.Pt()+pz1_type2*pz1_type2));
    //    W_mass_type2_met = (W_neutrino_type2_met+W_mu).M();
    //    W_pz_type2_met = (W_neutrino_type2_met+W_mu).Pz();
    //    W_nu1_pz_type2_met = pz1_type2;
    //    W_nu2_pz_type2_met = pz2_type2;

    // chenge the neutrino pT in case of complex solution in order to make it real
    TLorentzVector W_neutrino_type2; 
    W_neutrino_type2.SetPxPyPzE(W_Met.Px(),W_Met.Py(),pz1_type2,sqrt(W_Met.Pt()*W_Met.Pt()+pz1_type2*pz1_type2));

    if (NeutrinoPz_type2.IsComplex()) {// if this is a complix, change MET
      double nu_pt1 = NeutrinoPz_type2.getPtneutrino(1);
      double nu_pt2 = NeutrinoPz_type2.getPtneutrino(2);
      TLorentzVector W_neutrino_1;
      W_neutrino_1.SetPxPyPzE(nu_pt1 * TMath::Cos(ReducedTree->METPhi),
			      nu_pt1 * TMath::Sin(ReducedTree->METPhi), pz1_type2, sqrt(nu_pt1*nu_pt1 + pz1_type2*pz1_type2) );
      TLorentzVector W_neutrino_2;
      W_neutrino_2.SetPxPyPzE(nu_pt2 * TMath::Cos(ReducedTree->METPhi),
			      nu_pt2 * TMath::Sin(ReducedTree->METPhi), pz1_type2, sqrt(nu_pt2*nu_pt2 + pz1_type2*pz1_type2) );

      if ( fabs((W_mu+W_neutrino_1).M()-Wmass) < fabs((W_mu+W_neutrino_2).M()-Wmass) ) W_neutrino_type2 = W_neutrino_1;
      else W_neutrino_type2 = W_neutrino_2;
    }

    //    W_mass_type2 = (W_mu+W_neutrino_type2).M();
    //    W_pz_type2 = (W_mu+W_neutrino_type2).Pz();
    //    W_nu1_pz_type2 = pz1_type2;
    //    W_nu2_pz_type2 = pz2_type2;

    WWTree->pfMET   = sqrt(ReducedTree->METPt*ReducedTree->METPt);
    WWTree->pfMET_Phi = ReducedTree->METPhi;
    WWTree->nu_pz_type0 = pz1_type0;
    WWTree->nu_pz_type2 = pz1_type2;
    WWTree->nu_pz_run2 = pz1_run2;
    WWTree->nu_pz_isre = 1-NeutrinoPz_run2.IsComplex();
    WWTree->nu_pz_run2_oth = NeutrinoPz_run2.getOther();
    WWTree->nu_pz_run2_type = NeutrinoPz_run2.getType();


    /////////////////THE LEPTONIC W
    
    LEP.SetPtEtaPhiE(WWTree->l_pt,WWTree->l_eta,WWTree->l_phi,WWTree->l_e);
    NU0.SetPxPyPzE(ReducedTree->METPt*TMath::Cos(ReducedTree->METPhi),ReducedTree->METPt*TMath::Sin(ReducedTree->METPhi),WWTree->nu_pz_type0,TMath::Sqrt(WWTree->pfMET*WWTree->pfMET+WWTree->nu_pz_type0*WWTree->nu_pz_type0));
    NU2.SetPxPyPzE(ReducedTree->METPt*TMath::Cos(ReducedTree->METPhi),ReducedTree->METPt*TMath::Sin(ReducedTree->METPhi),WWTree->nu_pz_type2,TMath::Sqrt(WWTree->pfMET*WWTree->pfMET+WWTree->nu_pz_type2*WWTree->nu_pz_type2));
    NU1.SetPxPyPzE(ReducedTree->METPt*TMath::Cos(ReducedTree->METPhi),ReducedTree->METPt*TMath::Sin(ReducedTree->METPhi),WWTree->nu_pz_run2,TMath::Sqrt(WWTree->pfMET*WWTree->pfMET+WWTree->nu_pz_run2*WWTree->nu_pz_run2));
    W = LEP + NU2;
    
    WWTree->v_pt = W.Pt();
    WWTree->v_eta = W.Eta();
    WWTree->v_phi = W.Phi();
    WWTree->v_mt = TMath::Sqrt(2*LEP.Et()*NU2.Et()*(1-TMath::Cos(LEP.DeltaPhi(NU2))));
    //    W_mt = W.Mt();


    //FOR THE SYNCHORNIZATION!!! REMOVE IT FOR THE REAL ANALYSIS!!!!
    //    NU2.SetPtEtaPhiE(ReducedTree->METPt,0.,ReducedTree->METPhi,0.);
    //    W = NU2+LEP; 
    ////

	pt_W->Fill(W.Pt());
    //if (W.Pt()<50) continue;
    if(WWTree->event==evento) std::cout<<"debug: "<<count<<std::endl; count++;


    //    if (WWTree->v_pt < 150) continue;
//    if (WWTree->deltaR_lak8jet < (TMath::Pi()/2.0))   continue;

//===================================== START::AK4 Jets Selection ============================================
//
    std::vector<int> indexGoodJets;
    indexGoodJets.clear();
    WWTree->njets=0;
    //WWTree->nBTagJet_loose=0;
    //WWTree->nBTagJet_medium=0;
    //WWTree->nBTagJet_tight=0;

    float oldDeltaR = 1000.;
    float oldDeltaRLep = 1000.;
    int indexCloserJet = -1;
    int indexCloserJetLep = -1;
  //int TotalAK4Jets = 0, TotalAK4Jets_MoreThan4 = 0;
  int BoolIsJet = 0, BoolJetPtEtaPass = 0, BoolJetLoosePass = 0, BoolJetLepCleanPass = 0;
    for (unsigned int i=0; i<ReducedTree->JetsNum; i++) //loop on AK4 jet
      {
      BoolIsJet = 1;
	bool isCleanedJet = true;
	if (ReducedTree->Jets_isLooseJetId[i]==false) continue;
	if (fabs(ReducedTree->JetsEta[i])>=3.0) continue;
	if (ReducedTree->Jets_PtCorr[i]<=25) continue;
	//if (ReducedTree->Jets_PtCorr[i]<=25 || fabs(ReducedTree->JetsEta[i])>=3.0)  continue;
	BoolJetPtEtaPass = 1;
	//if (ReducedTree->Jets_PtCorr[i]<=30 || ReducedTree->JetsPt[i]<=20 || fabs(ReducedTree->JetsEta[i])>=2.4)  continue;
	BoolJetLoosePass = 1;


	//CLEANING FROM LEPTONS
	//cout<<"\n\n\n=============================================\n\n\n"<<endl;
	//cout<<"tightEle size = "<<tightEle.size()<<endl;
	//if (tightEle.size()!=1) exit(EXIT_FAILURE);
	//cout<<"\n\n\n=============================================\n\n\n"<<endl;
	for (int j=0; j<tightEle.size(); j++) {
	  if (deltaR(tightEle.at(j).Eta(), tightEle.at(j).Phi(),
		     ReducedTree->JetsEta[i],   ReducedTree->JetsPhi[i]) <0.3) {
	    isCleanedJet = false;
	  }
	}
	for (int j=0; j<tightMuon.size(); j++) {
	  if (deltaR(tightMuon.at(j).Eta(), tightMuon.at(j).Phi(),
		     ReducedTree->JetsEta[i],   ReducedTree->JetsPhi[i]) <0.3) {
	    isCleanedJet = false;
	  }
	}

	if (isCleanedJet==false) continue;
	BoolJetLepCleanPass = 1;


	WWTree->njetsAK4++;
	TotalAK4Jets++;

	AK4_new.SetPtEtaPhiE(ReducedTree->Jets_PtCorr[i],ReducedTree->JetsEta[i],ReducedTree->JetsPhi[i],ReducedTree->Jets_ECorr[i]);

	float deltaRlep = W.DeltaR(AK4_new);
	if (deltaRlep<oldDeltaRLep) indexCloserJetLep = i;

	indexGoodJets.push_back(i); //save index of the "good" vbf jets candidate
//	cout<<"yes there is jet"<<endl;
      }
//      cout<<"=========================================================="<<endl;
  //int BoolTotalAK4Jets = 0, BoolTotalAK4Jets_MoreThan4 = 0;
//  int BoolIsJet = 0, BoolJetPtEtaPass = 0, BoolJetLoosePass = 0, BoolJetLepCleanPass = 0;

	if (BoolIsJet)	IsJet++;
	if (BoolJetPtEtaPass) JetPtEtaPass++;
	if (BoolJetLoosePass) JetLoosePass++;
	if (BoolJetLepCleanPass) JetLepCleanPass++;
  

	
      if (indexGoodJets.size()<4)  continue;
      TotalAK4Jets_MoreThan4++;

      // Assign first two highest Pt jet as VBF tagged jets, and
      // Next two jets as W-jets


	float DeltaEta = 0.;
	int nVBF1=-1, nVBF2=-1; //position of the two vbf jets

	int nGoodAK4VBFjets = 0;
	if(verbose)
	for(int i=0; i<indexGoodJets.size();i++)
	{
	cout<<"Event = "<<jentry<<"\tpt [ "<<i<<" ] = "<<ReducedTree->Jets_PtCorr[indexGoodJets.at(i)]<<endl;
	}

	for(int i=0; i<indexGoodJets.size()-1;i++)
	{
		for(int j=i+1; j<indexGoodJets.size();j++)
		{
			VBF1.SetPtEtaPhiE(ReducedTree->Jets_PtCorr[indexGoodJets.at(i)],ReducedTree->JetsEta[indexGoodJets.at(i)],ReducedTree->JetsPhi[indexGoodJets.at(i)],ReducedTree->Jets_ECorr[indexGoodJets.at(i)]);
			VBF2.SetPtEtaPhiE(ReducedTree->Jets_PtCorr[indexGoodJets.at(j)],ReducedTree->JetsEta[indexGoodJets.at(j)],ReducedTree->JetsPhi[indexGoodJets.at(j)],ReducedTree->Jets_ECorr[indexGoodJets.at(j)]);
			//cout<<"Found Before Check!!!!"<<endl;
	if(verbose)
			cout<<"Before if loop::DeltaEta = "<<abs(VBF1.Eta()-VBF2.Eta())<<"opp hemi = "<< VBF1.Eta()*VBF2.Eta()*cos(VBF1.Theta()-VBF2.Theta()) <<"\t mass of dijet = "<<(VBF1+VBF2).M()<<endl;
			if (i==0 && j==1)
				mjj_01->Fill((VBF1+VBF2).M());
			if (i==0 && j==2)
				mjj_02->Fill((VBF1+VBF2).M());
			if (i==0 && j==3)
				mjj_03->Fill((VBF1+VBF2).M());
			if (i==1 && j==2)
				mjj_12->Fill((VBF1+VBF2).M());
			if (i==1 && j==3)
				mjj_13->Fill((VBF1+VBF2).M());
			if (i==2 && j==3)
				mjj_23->Fill((VBF1+VBF2).M());
			//if (DeltaEta > abs(VBF1.Eta()-VBF2.Eta()) || VBF1.Eta()*VBF2.Eta()*cos(VBF1.Theta()-VBF2.Theta()) > 0) continue;

			if (DeltaEta > abs(VBF1.Eta()-VBF2.Eta()) || VBF1.Eta()*VBF2.Eta()*cos(VBF1.Theta()-VBF2.Theta()) > 0 || (VBF1+VBF2).M()<200) continue;
			if (abs(VBF1.Eta()-VBF2.Eta())<1.5) continue;
			//if (DeltaEta < abs(VBF1.Eta()-VBF2.Eta()) &&  VBF1.Eta()*VBF2.Eta()*cos(VBF1.Theta()-VBF2.Theta()) > 0 && (VBF1+VBF2).M()<300) continue;
			//if (DeltaEta < abs(VBF1.Eta()-VBF2.Eta())) continue;
			//if (VBF1.Eta()*VBF2.Eta()*cos(VBF1.Theta()-VBF2.Theta()) > 0 ) continue;
			//if ((VBF1+VBF2).M()<300) continue;
	if(verbose)

			cout<<"Found!!!!"<<endl;
			DeltaEta = abs(VBF1.Eta()-VBF2.Eta()); //take the jet pair with largest DeltaEta
			nVBF1 = indexGoodJets.at(i); //save position of the 1st vbf jet
			nVBF2 = indexGoodJets.at(j); //save position of the 2nd vbf jet
	if(verbose)
			cout<<"nVBF1 = "<<nVBF1<<"\tnVBF2 = "<<nVBF2<<"\tDeltaEta = "<<DeltaEta<<"\tEta1*Eta2 = "<<VBF1.Eta()*VBF2.Eta()*cos(VBF1.Theta()-VBF2.Theta())<<"\tMass = "<<(VBF1+VBF2).M()<<endl;
		}
	nGoodAK4VBFjets++;
	}

	if (nGoodAK4VBFjets == 0) continue;
	cutEff[3]++;

	if (nVBF1!=-1 && nVBF2!=-1) //save infos for vbf jet pair
	{
    //cutEff[2]++;
		VBF1.SetPtEtaPhiE(ReducedTree->Jets_PtCorr[nVBF1],ReducedTree->JetsEta[nVBF1],ReducedTree->JetsPhi[nVBF1],ReducedTree->Jets_ECorr[nVBF1]);
		VBF2.SetPtEtaPhiE(ReducedTree->Jets_PtCorr[nVBF2],ReducedTree->JetsEta[nVBF2],ReducedTree->JetsPhi[nVBF2],ReducedTree->Jets_ECorr[nVBF2]);
		TOT = VBF1 + VBF2 ;

	if(verbose)
	cout<<"nVBF1 = "<<nVBF1<<"\tnVBF2 = "<<nVBF2<<endl;
	if(verbose)
	cout<<"VBF Jet1 = "<<ReducedTree->Jets_PtCorr[nVBF1]<<"\tVBF Jet1 = "<<ReducedTree->Jets_PtCorr[nVBF2]<<endl;

	    WWTree->vbf_AK4_j1_pt = ReducedTree->Jets_PtCorr[nVBF1];
	    WWTree->vbf_AK4_j1_eta = ReducedTree->JetsEta[nVBF1];
	    WWTree->vbf_AK4_j1_phi = ReducedTree->JetsPhi[nVBF1];
	    WWTree->vbf_AK4_j1_e = ReducedTree->Jets_ECorr[nVBF1];
	    WWTree->vbf_AK4_j1_bDiscriminatorCSV = ReducedTree->Jets_bDiscriminatorICSV[nVBF1];
	    WWTree->vbf_AK4_j2_pt = ReducedTree->Jets_PtCorr[nVBF2];
	    WWTree->vbf_AK4_j2_eta = ReducedTree->JetsEta[nVBF2];
	    WWTree->vbf_AK4_j2_phi = ReducedTree->JetsPhi[nVBF2];
	    WWTree->vbf_AK4_j2_e = ReducedTree->Jets_ECorr[nVBF2];
	    WWTree->vbf_AK4_j2_bDiscriminatorCSV = ReducedTree->Jets_bDiscriminatorICSV[nVBF2];
	    WWTree->vbf_AK4_jj_pt = TOT.Pt();
	    WWTree->vbf_AK4_jj_eta = TOT.Eta();
	    WWTree->vbf_AK4_jj_phi = TOT.Phi();
	    WWTree->vbf_AK4_jj_m = TOT.M();	
	    WWTree->vbf_AK4_jj_DeltaEta = fabs(VBF1.Eta()-VBF2.Eta());	

	}

	int coutWjets = 0;	
	int nWjets1 = -1, nWjets2 = -1 ;
	for(int i=0; i<indexGoodJets.size();i++)
	{
	if(indexGoodJets.at(0) == nVBF1 || indexGoodJets.at(0) == nVBF2 ) continue;
	coutWjets++;
	if(coutWjets == 1)
	{
	Wjet1_AK4.SetPtEtaPhiE(ReducedTree->Jets_PtCorr[indexGoodJets.at(i)],ReducedTree->JetsEta[indexGoodJets.at(i)],ReducedTree->JetsPhi[indexGoodJets.at(i)],ReducedTree->Jets_ECorr[indexGoodJets.at(i)]);
	nWjets1 = indexGoodJets.at(i);	// Save position of first w-jets
	}
        if(coutWjets == 2)
	{
	Wjet2_AK4.SetPtEtaPhiE(ReducedTree->Jets_PtCorr[indexGoodJets.at(i)],ReducedTree->JetsEta[indexGoodJets.at(i)],ReducedTree->JetsPhi[indexGoodJets.at(i)],ReducedTree->Jets_ECorr[indexGoodJets.at(i)]);
	nWjets2 = indexGoodJets.at(i);  // Save position of second w-jets
	}
	if ( nWjets1 == -1 || nWjets2 == -1 ) continue;
	cout<<nWjets1<<"\t"<<nWjets2<<endl;
	cutEff[4]++;
	TOT_Wjet = Wjet1_AK4 + Wjet1_AK4 ;

	    WWTree->Wjets_AK4_j1_pt = ReducedTree->Jets_PtCorr[nWjets1];
	    WWTree->Wjets_AK4_j1_eta = ReducedTree->JetsEta[nWjets1];
	    WWTree->Wjets_AK4_j1_phi = ReducedTree->JetsPhi[nWjets1];
	    WWTree->Wjets_AK4_j1_e = ReducedTree->Jets_ECorr[nWjets1];
	    WWTree->Wjets_AK4_j1_bDiscriminatorCSV = ReducedTree->Jets_bDiscriminatorICSV[nWjets1];
	    WWTree->Wjets_AK4_j2_pt = ReducedTree->Jets_PtCorr[nWjets2];
	    WWTree->Wjets_AK4_j2_eta = ReducedTree->JetsEta[nWjets2];
	    WWTree->Wjets_AK4_j2_phi = ReducedTree->JetsPhi[nWjets2];
	    WWTree->Wjets_AK4_j2_e = ReducedTree->Jets_ECorr[nWjets2];
	    WWTree->Wjets_AK4_j2_bDiscriminatorCSV = ReducedTree->Jets_bDiscriminatorICSV[nWjets2];
	    WWTree->Wjets_AK4_jj_pt =  TOT_Wjet.Pt();
	    WWTree->Wjets_AK4_jj_eta = TOT_Wjet.Eta();
	    WWTree->Wjets_AK4_jj_phi = TOT_Wjet.Phi();
	    WWTree->Wjets_AK4_jj_m =   TOT_Wjet.M();	
	    WWTree->Wjets_AK4_jj_e =   TOT_Wjet.E();	

	    nWjets1 = nWjets2 = -1;

	    }
	
	

	
	

//
//===================================== END::AK4 W-Jets Selection ============================================

    //////////////////ANGULAR VARIABLES
    JET.SetPtEtaPhiE(WWTree->ungroomed_jet_pt,WWTree->ungroomed_jet_eta,WWTree->ungroomed_jet_phi,WWTree->ungroomed_jet_e);
    WWTree->deltaR_lak8jet = JET.DeltaR(LEP);
    WWTree->deltaphi_METak8jet = JET.DeltaPhi(NU2);
    WWTree->deltaphi_Vak8jet = JET.DeltaPhi(W);
    if (WWTree->deltaR_lak8jet>(TMath::Pi()/2.0) && fabs(WWTree->deltaphi_METak8jet)>2.0 && fabs(WWTree->deltaphi_Vak8jet)>2.0)
      WWTree->issignal=1;

    //FOUR-BODY INVARIANT MASS
    WWTree->mass_lvj_type0 = (LEP + NU0 + JET).M();
    WWTree->mass_lvj_type2 = (LEP + NU2 + JET).M();
    WWTree->mass_lvj_run2  = (LEP + NU1 + JET).M();

    /////////////////MC Infos
    if (isMC)
      {
	TLorentzVector hadW, lepW, temp;
	int posWhad =-1, posWlep =-1, posTemp=-1, posGenJet=-1;
	//	std::cout<<"entry: "<<iEntry<<" "<<GenNuNum<<std::endl;
	double deltaPhiOld=100.;
	WWTree->genGravMass=100.;	

	for (int i=0; i<ReducedTree->GenBosonNum; i++) {
	  for (int j=i+1; j<ReducedTree->GenBosonNum; j++) {

	    hadW.SetPtEtaPhiE(ReducedTree->GenBosonPt[i],ReducedTree->GenBosonEta[i],ReducedTree->GenBosonPhi[i],ReducedTree->GenBosonE[i]);
	    lepW.SetPtEtaPhiE(ReducedTree->GenBosonPt[j],ReducedTree->GenBosonEta[j],ReducedTree->GenBosonPhi[j],ReducedTree->GenBosonE[j]);


	    if (fabs((hadW+lepW).M()-genMass)< fabs(WWTree->genGravMass-genMass)) { //found the gen graviton
	      WWTree->genGravMass=(hadW+lepW).M();	
	      posWhad=i; //save positions of the two W's in random order, will fix them later in the code
	      posWlep=j;

	    }

	  }	
	}

	if (posWhad!=-1 && posWlep!=-1) {

	  float oldDR=100.;
	  bool isWhadOk=true;

	  if(WWTree->event==127270357) std::cout<<"debug: "<<std::endl;

	  for (int i=0; i<ReducedTree->GenJetsAK8Num; i++) {
	      if (deltaR(ReducedTree->GenBosonEta[posWhad], ReducedTree->GenBosonPhi[posWhad],
			 ReducedTree->GenJetsAK8Eta[i], ReducedTree->GenJetsAK8Phi[i])< oldDR ) 
		{
		  posGenJet=i;		
		  oldDR = deltaR(ReducedTree->GenBosonEta[posWhad], ReducedTree->GenBosonPhi[posWhad],ReducedTree->GenJetsAK8Eta[i], ReducedTree->GenJetsAK8Phi[i]);
		  isWhadOk = true;
		  if(WWTree->event==127270357) std::cout<<"debug: had "<<ReducedTree->GenJetsAK8Phi[i]<<" "<<ReducedTree->GenBosonPhi[posWhad]<<" "<<oldDR<<std::endl;
		}
	      if (deltaR(ReducedTree->GenBosonEta[posWlep], ReducedTree->GenBosonPhi[posWlep],
			 ReducedTree->GenJetsAK8Eta[i], ReducedTree->GenJetsAK8Phi[i])< oldDR ) 
		{
		  posGenJet=i;
		  oldDR = deltaR(ReducedTree->GenBosonEta[posWlep], ReducedTree->GenBosonPhi[posWlep],ReducedTree->GenJetsAK8Eta[i], ReducedTree->GenJetsAK8Phi[i]);
		  isWhadOk = false;
		  if(WWTree->event==127270357) std::cout<<"debug: lep "<<ReducedTree->GenJetsAK8Phi[i]<<" "<<ReducedTree->GenBosonPhi[posWlep]<<" "<<oldDR<<std::endl;
         	}      
	  }
	  if(WWTree->event==127270357) std::cout<<ReducedTree->GenJetsAK8Phi[posGenJet]<<std::endl;
	  if (isWhadOk==false) //wrong W's positions saved, switch them
	    {
	      posTemp = posWhad;
	      posWhad = posWlep;
	      posWlep = posTemp;
	    }

	  hadW.SetPtEtaPhiE(ReducedTree->GenBosonPt[posWhad],ReducedTree->GenBosonEta[posWhad],ReducedTree->GenBosonPhi[posWhad],ReducedTree->GenBosonE[posWhad]);
	  lepW.SetPtEtaPhiE(ReducedTree->GenBosonPt[posWlep],ReducedTree->GenBosonEta[posWlep],ReducedTree->GenBosonPhi[posWlep],ReducedTree->GenBosonE[posWlep]);

	  WWTree->W_pt_gen = ReducedTree->GenBosonPt[posWlep];
	  WWTree->W_pz_gen = lepW.Pz();
	  WWTree->W_rap_gen = lepW.Rapidity();
	  
	  WWTree->hadW_pt_gen = ReducedTree->GenBosonPt[posWhad];
	  WWTree->hadW_eta_gen = ReducedTree->GenBosonEta[posWhad];
	  WWTree->hadW_phi_gen = ReducedTree->GenBosonPhi[posWhad];
	  WWTree->hadW_e_gen = ReducedTree->GenBosonE[posWhad];
	  WWTree->hadW_m_gen = hadW.M();

	  WWTree->lepW_pt_gen = ReducedTree->GenBosonPt[posWlep];
	  WWTree->lepW_eta_gen = ReducedTree->GenBosonEta[posWlep];
	  WWTree->lepW_phi_gen = ReducedTree->GenBosonPhi[posWlep];
	  WWTree->lepW_e_gen = ReducedTree->GenBosonE[posWlep];
	  WWTree->lepW_m_gen = lepW.M();

	  WWTree->AK8_pt_gen = ReducedTree->GenJetsAK8Pt[posGenJet];
	  WWTree->AK8_eta_gen = ReducedTree->GenJetsAK8Eta[posGenJet];
	  WWTree->AK8_phi_gen = ReducedTree->GenJetsAK8Phi[posGenJet];
	  WWTree->AK8_e_gen = ReducedTree->GenJetsAK8E[posGenJet];
	  WWTree->AK8_pruned_mass_gen = ReducedTree->GenJetsAK8_prunedMass[posGenJet];
	  WWTree->AK8_softdrop_mass_gen = ReducedTree->GenJetsAK8_softdropMass[posGenJet];

          if (deltaR(ReducedTree->GenBosonEta[posWhad], ReducedTree->GenBosonPhi[posWhad],
                     WWTree->ungroomed_jet_eta, WWTree->ungroomed_jet_phi)<0.1)     ok++;
          total++;
	}

	deltaPhiOld=100.;
       	for (int i=0; i<ReducedTree->GenNuNum; i++) {
	  double deltaPhi = getDeltaPhi(ReducedTree->GenNuPhi[i],WWTree->v_phi);
	  if (abs(deltaPhi)>abs(deltaPhiOld))   continue;	  
	  temp.SetPtEtaPhiE(ReducedTree->GenNuPt[i],ReducedTree->GenNuEta[i],ReducedTree->GenNuPhi[i],ReducedTree->GenNuE[i]);
	  WWTree->nu_pz_gen=temp.Pz();	  
	  WWTree->nu_pt_gen=temp.Pt();	  
	  WWTree->nu_phi_gen=temp.Phi();	  
	  WWTree->nu_eta_gen=temp.Eta();
	  deltaPhiOld = deltaPhi;
	}		
      }
    
    //fill the tree
    if(WWTree->event==evento) std::cout<<"fill: "<<count<<std::endl; count++;
    outTree->Fill();
  }
  ptEle->SaveAs("ptEle.C");
  ptMet->SaveAs("ptMet.C");
  pt_W->SaveAs("pt_W.C");
  mjj_01->SaveAs("mjj_01.C");
  mjj_02->SaveAs("mjj_02.C");
  mjj_03->SaveAs("mjj_03.C");
  mjj_12->SaveAs("mjj_12.C");
  mjj_13->SaveAs("mjj_13.C");
  mjj_23->SaveAs("mjj_23.C");


//  sprintf("Total entries = %d \n Total Leptons = %d \n Medium Id = %d \n pt>25 = %d \n eta>2.5 = %d \n veto & 1Lep = %d \n met >30 = %d \n Wpt > 50 = %d \n HaveJet = %d \n JetPt Eta = %d \n JetLoose Id = %d \n Lepton Cleaning = %d \n Morethanor4Jet = %d \n",ReducedTree->fChain->GetEntries(), TotalEle, MediumPassEle, PtPassEle, EtaPassEle, cutEff[0] , cutEff[1], cutEff[2], IsJet, JetPtEtaPass, JetLoosePass, JetLepCleanPass, TotalAK4Jets_MoreThan4);
  cout<<"TotalEle = "<< TotalEle <<"\tTriggerPassEle = "<<  TriggerPassEle <<"\t MediumPassEle = "<<  MediumPassEle <<"\t PtPassEle = "<<  PtPassEle <<"\t EtaPassEle = "<<  EtaPassEle <<endl;
  cout<<" TotalMu= "<< TotalMu <<"\t TriggerPassMu = "<< TriggerPassMu <<"\t TightPassMu = "<< TightPassMu <<"\t PtPassMu = "<< PtPassMu <<"\t EtaPassMu = "<< EtaPassMu <<"\t IsoPassMu = "<<  IsoPassMu <<endl;;

  cout<<"TotalAK4Jets = "<<TotalAK4Jets<<"\tTotalAK4Jets_MoreThan4 = "<<TotalAK4Jets_MoreThan4<<endl;
  cout<<"IsJet = "<<IsJet<<"\tJetPtEtaPass = "<<JetPtEtaPass<<"\tJetLoosePass = "<<JetLoosePass<<"\tJetLepCleanPass = "<<JetLepCleanPass<<endl;

  std::cout<<"matching: "<<(float)ok/(float)total<<std::endl;

  std::cout<<"total entries: "<<ReducedTree->fChain->GetEntries()<<std::endl;
  std::cout<<"lepton eff: "<<cutEff[0]<<std::endl
	   <<"met eff:    "<<cutEff[1]<<std::endl
	   <<"W eff:      "<<cutEff[2]<<std::endl
	   <<"VBF Jet found:  "<<cutEff[3]<<std::endl
	   <<"WJets found:  "<<cutEff[4]<<std::endl;

  //--------close everything-------------
  ReducedTree->fChain->Delete();
  outTree->Write();
  outROOT->Close();

  return(0);
}
