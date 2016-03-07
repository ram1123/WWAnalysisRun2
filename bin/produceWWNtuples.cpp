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
void computeAngles(TLorentzVector thep4H, TLorentzVector thep4Z1, TLorentzVector thep4M11, TLorentzVector thep4M12, TLorentzVector thep4Z2, TLorentzVector thep4M21, TLorentzVector thep4M22, double& costheta1, double& costheta2, double& Phi, double& costhetastar, double& Phi1);
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
  TLorentzVector VBF1,VBF2,TOT, VBF1_AK4, VBF2_AK4, Wjet1_AK4, Wjet2_AK4, TOT_Wjet, TOT_Wjet_Final , Wjets1, Wjets2;
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
  char list1[2000];
  //sprintf (list1, "InputRootFiles/%s.txt", inputFile.c_str());
  sprintf (list1, "%s.txt", inputFile.c_str());
  ifstream rootList (list1);

  int fileCounter=0;
  Long64_t totalEntries=0;
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
  std::cout<<"number of files found: "<<fileCounter-1<<std::endl;
  std::cout<<"total entries: "<<ReducedTree->fChain->GetEntries()<<std::endl;
  totalEntries=ReducedTree->fChain->GetEntries();

  char command3[300];
  //sprintf(command3, "rm listTemp_%s.txt", outputFile.c_str());
  //system(command3);

  int cutEff[20]={0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};

  //--------pile up file -----------------
  TFile* pileupFile = TFile::Open("PU.root");  
  TH1F *pileupHisto = (TH1F*)pileupFile->Get("puweights");

  std::vector<double> weights;

  weights = generate_weights(pileupHisto);
  pileupFile->Close();

  //---------output tree----------------
  //TFile* outROOT = TFile::Open((std::string("/eos/uscms/store/user/rasharma/WWScattering/WWTrees_13Jan2016/output/output_")+leptonName+std::string("/")+outputFile+(".root")).c_str(),"recreate");	// Path here for lpc only because at lpc condor jobs output can be only saved at eos area
  TFile* outROOT = TFile::Open((leptonName+std::string("_")+outputFile+(".root")).c_str(),"recreate");
  outROOT->cd();
  TTree* outTree = new TTree("otree", "otree");
  outTree->SetDirectory(0);

  setOutputTree *WWTree = new setOutputTree(outTree);

  //---------start loop on events------------
  Long64_t jentry2=0;
/* 
 *
 *Patch for New Tree
 *
 */
 //TFile f("tree1.root","recreate");	// CutAnalysis
 TFile f((leptonName+std::string("_")+outputFile+("_1.root")).c_str(),"recreate");
 TTree t1("Leptopn","Tree for understanding Cut Flow"); // CutAnalysis
 TTree t2("Jet","Tree for understanding JET cut Flow: After Selecting Lepton and Good Jet"); // CutAnalysis
 TTree t3("MET","Tree for understanding MET: After Selecting Lepton and Good Jet"); // CutAnalysis

 Int_t Ele_num, Mu_num, Ele_nPV, Mu_nPV;
 Double_t Ele_pt, Ele_eta, Ele_phi, Ele_E;
 Double_t Mu_pt, Mu_eta, Mu_phi, Mu_E;

 t1.Branch("Ele_num"	,	&Ele_num	,	"Ele_num/I"	);
 t1.Branch("Ele_pt"	,	&Ele_pt		,	"Ele_pt/D"	);
 t1.Branch("Ele_eta"	,	&Ele_eta	,	"Ele_eta/D"	);
 t1.Branch("Ele_phi"	,	&Ele_phi	,	"Ele_phi/D"	);
 t1.Branch("Ele_nPV"	,	&Ele_nPV	,	"Ele_nPV/I"	);

 t1.Branch("Mu_num"	,	&Mu_num		,	"Mu_num/I"	);
 t1.Branch("Mu_pt"	,	&Mu_pt		,	"Mu_pt/D"	);
 t1.Branch("Mu_eta"	,	&Mu_eta		,	"Mu_eta/D"	);
 t1.Branch("Mu_phi"	,	&Mu_phi		,	"Mu_phi/D"	);
 t1.Branch("Mu_nPV"	,	&Mu_nPV		,	"Mu_nPV/I"	);

 Int_t Jet_num_1, Jet_num_2, Jet0_nPV;
 Double_t Jet0_pt, Jet0_eta, Jet0_phi, Jet0_E;

 t2.Branch("Jet_num_1"	,	&Jet_num_1	,	"Jet_num_1/I"	);
 t2.Branch("Jet_num_2"	,	&Jet_num_2	,	"Jet_num_2/I"	);
 t2.Branch("Jet0_pt"	,	&Jet0_pt	,	"Jet0_pt/D"	);
 t2.Branch("Jet0_eta"	,	&Jet0_eta	,	"Jet0_eta/D"	);
 t2.Branch("Jet0_phi"	,	&Jet0_phi	,	"Jet0_phi/D"	);
 t2.Branch("Jet0_E"	,	&Jet0_E		,	"Jet0_E/D"	);
 t2.Branch("Jet0_nPV"	,	&Jet0_nPV	,	"Jet0_nPV/I"	);

 Int_t pf_MET_nPV;
 Double_t pf_MET_pt , pf_MET_phi ;
 t3.Branch("pf_MET_pt"	,	&pf_MET_pt	,	"pf_MET_pt/D"	);
 t3.Branch("pf_MET_phi"	,	&pf_MET_phi	,	"pf_MET_phi/D"	);
 t3.Branch("pf_MET_nPV"	,	&pf_MET_nPV	,	"pf_MET_nPV/I"	);
/*
 * END
 */

  int MediumPassEle = 0, PtPassEle = 0, EtaPassEle = 0;
  int TightPassMu = 0, EtaPassMu = 0, IsoPassMu = 0;
  int TotalAK4Jets = 0, TotalAK4Jets_MoreThan4 = 0;
  int IsJet = 0, JetPtEtaPass = 0, JetLoosePass = 0, JetLepCleanPass = 0;

  int BoolTriggerPassEle = 0, BoolMediumPassEle = 0, BoolPtPassEle = 0, BoolEtaPassEle = 0;
  int BoolTriggerPassMu = 0, BoolTightPassMu = 0, BoolPtPassMu = 0, BoolEtaPassMu = 0, BoolIsoPassMu = 0;
  int BoolTotalAK4Jets = 0, BoolTotalAK4Jets_MoreThan4 = 0;
  int BoolIsJet = 0, BoolJetPtEtaPass = 0, BoolJetLoosePass = 0, BoolJetLepCleanPass = 0;



  for (Long64_t jentry=0; jentry<ReducedTree->fChain->GetEntries();jentry++,jentry2++) {
  //for (Long64_t jentry=0; jentry<50000;jentry++,jentry2++) {

    Long64_t iEntry = ReducedTree->LoadTree(jentry);
    if (iEntry < 0) break;
    int nb = ReducedTree->fChain->GetEntry(jentry);   

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

    WWTree->genWeight = ReducedTree->genEventWeight;
 //   cout<<"ReducedTree->genEventWeight = "<<ReducedTree->genEventWeight<<endl;

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
    WWTree->njets = ReducedTree->NJets;
    WWTree->nPV  = ReducedTree->NVtx;

    if(WWTree->event==evento) std::cout<<"debug: "<<count<<std::endl; count++;

    if (verbose)
    	cout<<"==================> debug 1 "<<endl;

    /////////////////THE SELECTED LEPTON
    int nTightLepton=0;
   BoolMediumPassEle = 0, BoolPtPassEle = 0, BoolEtaPassEle = 0;
   BoolTightPassMu = 0, BoolPtPassMu = 0, BoolEtaPassMu = 0, BoolIsoPassMu = 0;
   BoolTotalAK4Jets = 0, BoolTotalAK4Jets_MoreThan4 = 0;
    if (strcmp(leptonName.c_str(),"el")==0) {
      float tempPt=0.;
      int passTrigger=0;
      for (int i=0; i<ReducedTree->ElectronsNum; i++) {
	if (applyTrigger==1)
	  for (int t=0; t<ReducedTree->TriggerProducerTriggerNames->size(); t++)
	  {
	  if(TString(ReducedTree->TriggerProducerTriggerNames->at(t)).Contains("HLT_Ele23_WPLoose_Gsf_v"))
	  if (ReducedTree->TriggerProducerTriggerPass->at(t)==1) {passTrigger=1; //trigger
	  cout<<"Enter Ele"<<endl;}
	  }
	  if (passTrigger==0) continue;
	if (ReducedTree->Electrons_isMedium[i]==false) continue;       
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
      int passTrigger=0;
      float tempPt=0.;
      for (int i=0; i<ReducedTree->MuonsNum; i++) {
	if (applyTrigger==1)
	  for (int t=0; t<ReducedTree->TriggerProducerTriggerNames->size(); t++)
	  {
	  if(TString(ReducedTree->TriggerProducerTriggerNames->at(t)).Contains("HLT_IsoMu20_v"))
	  if (ReducedTree->TriggerProducerTriggerPass->at(t)==1) {passTrigger=1; //trigger
	  cout<<"Enter MU"<<endl;}
	  }
	  if (passTrigger==0) continue;

	//if (ReducedTree->TriggerProducerTriggerPass->at(0)==0) continue; //trigger
	//if (ReducedTree->Muons_isLoose[i]==false) continue;
	if (ReducedTree->Muons_isTight[i]==false) continue;
	//if (ReducedTree->Muons_isHighPt[i]==false) continue;
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
	if(BoolIsoPassMu)
	IsoPassMu++;
	if(BoolTightPassMu)
	TightPassMu++;
//	if(BoolTriggerPassMu)
//	TriggerPassMu++;
       	if(BoolEtaPassEle)
	EtaPassEle++;
	if(BoolPtPassEle)
	PtPassEle++;
	if(BoolMediumPassEle)
	MediumPassEle++;
//	if(BoolTriggerPassEle)

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
 
    if (strcmp(leptonName.c_str(),"el")==0) {
    Ele_num	= nLooseLepton; 	// CutAnalysis
    Ele_pt	= WWTree->l_pt;
    Ele_eta	= WWTree->l_eta;
    Ele_phi	= WWTree->l_phi;
    Ele_nPV	= WWTree->nPV;
    }
    if (strcmp(leptonName.c_str(),"mu")==0) {
    Mu_num	= nLooseLepton; 	// CutAnalysis
    Mu_pt	= WWTree->l_pt;
    Mu_eta	= WWTree->l_eta;
    Mu_phi	= WWTree->l_phi;
    Mu_nPV	= WWTree->nPV;
    }

    t1.Fill();
    if(WWTree->event==evento) std::cout<<"debug: "<<count<<std::endl; count++;

    //preselection on jet pt and met
    WWTree->Met_pt  = ReducedTree->METPt;
    //if (ReducedTree->METPt < 20) continue; 
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
  Jet_num_1 = ReducedTree->JetsNum;

    for (unsigned int i=0; i<ReducedTree->JetsNum; i++) //loop on AK4 jet
      {
      BoolIsJet = 1;
	bool isCleanedJet = true;
	if (fabs(ReducedTree->JetsEta[i])>=3.0) continue;
	if (ReducedTree->Jets_PtCorr[i]<=25) continue;
	//if (ReducedTree->Jets_PtCorr[i]<=25 || fabs(ReducedTree->JetsEta[i])>=3.0)  continue;
	BoolJetPtEtaPass = 1;
	if (ReducedTree->Jets_isLooseJetId[i]==false) continue;
	//if (ReducedTree->Jets_PtCorr[i]<=30 || ReducedTree->JetsPt[i]<=20 || fabs(ReducedTree->JetsEta[i])>=2.4)  continue;
	BoolJetLoosePass = 1;


	//CLEANING FROM LEPTONS
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

	if (BoolIsJet)	IsJet++;
	if (BoolJetPtEtaPass) JetPtEtaPass++;
	if (BoolJetLoosePass) JetLoosePass++;
	if (BoolJetLepCleanPass) JetLepCleanPass++;
 
     Jet_num_2	= indexGoodJets.size();
     if (indexGoodJets.size()<1) continue;
     Jet0_pt	= ReducedTree->Jets_PtCorr[indexGoodJets.at(0)];
     Jet0_eta	= ReducedTree->JetsEta[indexGoodJets.at(0)];
     Jet0_phi	= ReducedTree->JetsPhi[indexGoodJets.at(0)];
     Jet0_E	= ReducedTree->Jets_ECorr[indexGoodJets.at(0)];
     Jet0_nPV	= WWTree->nPV;

     t2.Fill();

    //MET.SetPtEtaPhiE(ReducedTree->METPt,0.,ReducedTree->METPhi,0.);
     pf_MET_pt = ReducedTree->METPt;
     pf_MET_phi= ReducedTree->METPhi;
     pf_MET_nPV= WWTree->nPV;
	
     t3.Fill();
      if (indexGoodJets.size()<4)  continue;
      TotalAK4Jets_MoreThan4++;


	float DeltaEta = 0.;
	int nVBF1=-1, nVBF2=-1; //position of the two vbf jets

	int nGoodAK4VBFjets = 0;
	if(verbose)
	for(int i=0; i<indexGoodJets.size();i++)
	{
	cout<<"index = "<<indexGoodJets.at(i)<<"  Event = "<<jentry<<"\tpt [ "<<i<<" ] = "<<ReducedTree->Jets_PtCorr[indexGoodJets.at(i)]<<endl;
	}

	//================================ Selection of VBF jets: BEGIN	=====================================
	for(int i=0; i<indexGoodJets.size()-1;i++)
	{
		for(int j=i+1; j<indexGoodJets.size();j++)
		{
			VBF1.SetPtEtaPhiE(ReducedTree->Jets_PtCorr[indexGoodJets.at(i)],ReducedTree->JetsEta[indexGoodJets.at(i)],ReducedTree->JetsPhi[indexGoodJets.at(i)],ReducedTree->Jets_ECorr[indexGoodJets.at(i)]);
			VBF2.SetPtEtaPhiE(ReducedTree->Jets_PtCorr[indexGoodJets.at(j)],ReducedTree->JetsEta[indexGoodJets.at(j)],ReducedTree->JetsPhi[indexGoodJets.at(j)],ReducedTree->Jets_ECorr[indexGoodJets.at(j)]);
			//cout<<"Found Before Check!!!!"<<endl;
	if(verbose)
			cout<<"Before if loop::DeltaEta = "<<abs(VBF1.Eta()-VBF2.Eta())<<"opp hemi = "<< VBF1.Eta()*VBF2.Eta()*cos(VBF1.Theta()-VBF2.Theta()) <<"\t mass of dijet = "<<(VBF1+VBF2).M()<<endl;
			if (DeltaEta > abs(VBF1.Eta()-VBF2.Eta()) || VBF1.Eta()*VBF2.Eta() > 0 || (VBF1+VBF2).M()<500) continue;
			if (abs(VBF1.Eta()-VBF2.Eta())<3.5) continue;

	if(verbose)
			cout<<"Found!!!!"<<endl;
			DeltaEta = abs(VBF1.Eta()-VBF2.Eta()); //take the jet pair with largest DeltaEta
			nVBF1 = indexGoodJets.at(i); //save position of the 1st vbf jet
			nVBF2 = indexGoodJets.at(j); //save position of the 2nd vbf jet
	if(verbose)
			cout<<"nVBF1 = "<<nVBF1<<"\tnVBF2 = "<<nVBF2<<"\tDeltaEta = "<<DeltaEta<<"\tEta1*Eta2 = "<<VBF1.Eta()*VBF2.Eta()*cos(VBF1.Theta()-VBF2.Theta())<<"\tMass = "<<(VBF1+VBF2).M()<<endl;
	nGoodAK4VBFjets++;
		}
	}

	if (nGoodAK4VBFjets == 0) continue;
	cutEff[3]++;

	if (ReducedTree->Jets_bDiscriminatorICSV[nVBF1]> 0.97 && ReducedTree->Jets_bDiscriminatorICSV[nVBF2]> 0.97) continue;
	cutEff[6]++;


	    WWTree->vbf_AK4_j_Num = nGoodAK4VBFjets;	

	if (nVBF1!=-1 && nVBF2!=-1) //save infos for vbf jet pair
	{
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
	    WWTree->vbf_AK4_jj_DeltaPhi = fabs(VBF1.Phi()-VBF2.Phi());	
    	    WWTree->vbf_AK4_jj_AssymPt = (ReducedTree->Jets_PtCorr[nVBF1]-ReducedTree->Jets_PtCorr[nVBF2])/(ReducedTree->Jets_PtCorr[nVBF1] + ReducedTree->Jets_PtCorr[nVBF2]);
    	    //WWTree->vbf_AK4_jj_AssymPt = fabs((ReducedTree->Jets_PtCorr[nVBF1]-ReducedTree->Jets_PtCorr[nVBF2])/(ReducedTree->Jets_PtCorr[nVBF1] + ReducedTree->Jets_PtCorr[nVBF2]));

	}

	//================================ Selection of VBF jets: END	=====================================
	//
	
	//================================ Selection of W-jets: STARTS	=====================================
	int nGoodAK4Wjets = 0;	
	int nWjets1 = -1, nWjets2 = -1 ;
	double DeltaMassWindow = 25.;
	for(int i=0; i<indexGoodJets.size()-1;i++)
	{
		for(int j=i+1; j<indexGoodJets.size();j++)
		{
	if(indexGoodJets.at(i) == nVBF1 || indexGoodJets.at(j) == nVBF2 ) continue;
	if(indexGoodJets.at(j) == nVBF1 || indexGoodJets.at(i) == nVBF2 ) continue;
	//coutWjets++;
			Wjet1_AK4.SetPtEtaPhiE(ReducedTree->Jets_PtCorr[indexGoodJets.at(i)],ReducedTree->JetsEta[indexGoodJets.at(i)],ReducedTree->JetsPhi[indexGoodJets.at(i)],ReducedTree->Jets_ECorr[indexGoodJets.at(i)]);
			Wjet2_AK4.SetPtEtaPhiE(ReducedTree->Jets_PtCorr[indexGoodJets.at(j)],ReducedTree->JetsEta[indexGoodJets.at(j)],ReducedTree->JetsPhi[indexGoodJets.at(j)],ReducedTree->Jets_ECorr[indexGoodJets.at(j)]);
			TOT_Wjet = Wjet1_AK4 + Wjet2_AK4 ;
			//cout<<"Found Before Check!!!!"<<endl;


			if (DeltaMassWindow < abs(TOT_Wjet.M() - 80.)) continue;
			//cout<<"(Wjet1_AK4+Wjet2_AK4).M()-80. = "<<(Wjet1_AK4+Wjet2_AK4).M()<<endl;
			
//			cout<< "====> " << Wjet1_AK4.Mag() << "\t" << Wjet2_AK4.Mag() << "\t" << TOT_Wjet.Mag() <<endl;

	if(verbose)
			cout<<"Found!!!!"<<endl;
			DeltaMassWindow = abs(TOT_Wjet.M()-80.); //take the jet pair with largest DeltaEta
			nWjets1 = indexGoodJets.at(i); //save position of the 1st vbf jet
			nWjets2 = indexGoodJets.at(j); //save position of the 2nd vbf jet
			nGoodAK4Wjets++;
		}
	}
	if (nGoodAK4Wjets == 0) continue;
	cutEff[2]++;
	if ( nWjets1 == -1 || nWjets2 == -1 ) continue;
	//cout<<nWjets1<<"\t"<<nWjets2<<endl;
	Wjets1.SetPtEtaPhiE(ReducedTree->Jets_PtCorr[nWjets1],ReducedTree->JetsEta[nWjets1],ReducedTree->JetsPhi[nWjets1],ReducedTree->Jets_ECorr[nWjets1]);
	Wjets2.SetPtEtaPhiE(ReducedTree->Jets_PtCorr[nWjets2],ReducedTree->JetsEta[nWjets2],ReducedTree->JetsPhi[nWjets2],ReducedTree->Jets_ECorr[nWjets2]);
	TOT_Wjet_Final = Wjets1 + Wjets2 ;

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
	    WWTree->Wjets_AK4_jj_pt =  TOT_Wjet_Final.Pt();
	    WWTree->Wjets_AK4_jj_eta = TOT_Wjet_Final.Eta();
	    WWTree->Wjets_AK4_jj_phi = TOT_Wjet_Final.Phi();
	    WWTree->Wjets_AK4_jj_m =   TOT_Wjet_Final.M();	
	    WWTree->Wjets_AK4_jj_e =   TOT_Wjet_Final.E();	
	    WWTree->WW_mass = (TOT + TOT_Wjet_Final).M();
    	    WWTree->Wjets_AK4_jj_AssymPt = (ReducedTree->Jets_PtCorr[nWjets1]-ReducedTree->Jets_PtCorr[nWjets2])/(ReducedTree->Jets_PtCorr[nWjets1] + ReducedTree->Jets_PtCorr[nWjets2]);
    	    //WWTree->Wjets_AK4_jj_AssymPt = fabs((ReducedTree->Jets_PtCorr[nWjets1]-ReducedTree->Jets_PtCorr[nWjets2])/(ReducedTree->Jets_PtCorr[nWjets1] + ReducedTree->Jets_PtCorr[nWjets2]));
	
	
	TLorentzVector p4_WHad = Wjets1  + Wjets2;
        TLorentzVector p4_WLep = W ;        // Leptonic W 4-vector
        TLorentzVector p4_WW = p4_WHad + p4_WLep;
        
        double a_costheta1, a_costheta2, a_costhetastar, a_Phi, a_Phi1;
        computeAngles( p4_WW, p4_WLep, LEP , NU2, p4_WHad, Wjets1, Wjets2, 
                      a_costheta1, a_costheta2, a_Phi, a_costhetastar, a_Phi1);

	WWTree-> costheta1 = a_costheta1;
	WWTree-> costheta2 = a_costheta2;
	WWTree-> costhetastar = a_costhetastar;
	WWTree-> Phi = a_Phi;
	WWTree-> Phi1 = a_Phi1;


	//================================ Selection of W-jets: END	=====================================
//
//===================================== END::AK4 W-Jets Selection ============================================
/* 
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
*/
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


  std::cout<<"matching: "<<(float)ok/(float)total<<std::endl;

  std::cout<<"===========================	OUTPUTS	=========================="<<std::endl;
  std::cout<<"total entries: \t"<<ReducedTree->fChain->GetEntries()<<std::endl;
  if( strcmp(leptonName.c_str(),"el")==0)
  {
  std::cout<<"MediumPassEle: \t"	<<MediumPassEle		<<"\t"<< setprecision(3) << fixed <<(float)MediumPassEle/(float)ReducedTree->fChain->GetEntries()<<std::endl
  	   <<"Pt and Eta Pass ele:\t "	<<EtaPassEle		<<"\t"<< setprecision(3) << fixed <<(float)EtaPassEle/(float)MediumPassEle<<std::endl
	   <<"Electron eff: \t"		<<cutEff[0]		<<"\t"<< setprecision(3) << fixed <<(float)cutEff[0]/(float)EtaPassEle<<std::endl
	   <<"met eff:   \t "		<<cutEff[1]		<<"\t"<< setprecision(3) << fixed <<(float)cutEff[1]/(float)cutEff[0]<<std::endl
	   <<"Jet pt & eta :\t "	<<JetPtEtaPass		<<"\t"<< setprecision(3) << fixed <<(float)JetPtEtaPass/(float)cutEff[1]<<std::endl
	   <<"Jet loose ID: \t"		<<JetLoosePass		<<"\t"<< setprecision(3) << fixed <<(float)JetLoosePass/(float)JetPtEtaPass<<std::endl
	   <<"Jet Cleaning: \t"		<<JetLepCleanPass	<<"\t"<< setprecision(3) << fixed <<(float)JetLepCleanPass/(float)JetLoosePass<<std::endl
	   <<"Events >=4Jet: \t"	<<TotalAK4Jets_MoreThan4<<"\t"<< setprecision(3) << fixed <<(float)TotalAK4Jets_MoreThan4/(float)JetLepCleanPass<<std::endl
	   <<"VBF Jet found:  \t"	<<cutEff[3]		<<"\t"<< setprecision(3) << fixed <<(float)cutEff[3]/(float)TotalAK4Jets_MoreThan4<<std::endl
	   <<"VBF Jet No B-tagged:\t  "	<<cutEff[6]		<<"\t"<< setprecision(3) << fixed <<(float)cutEff[6]/(float)cutEff[3]<<std::endl
	   <<"Wjet found: \t"		<<cutEff[2]		<<"\t"<< setprecision(3) << fixed <<(float)cutEff[2]/(float)cutEff[6]<<std::endl;
  }
  if( strcmp(leptonName.c_str(),"mu")==0)
  {
  std::cout<<"TightPassEle: \t"		<<TightPassMu		<<"\t"<< setprecision(3) << fixed <<(float)TightPassMu/(float)ReducedTree->fChain->GetEntries()<<std::endl
	   <<"Isolation passed: \t"	<<IsoPassMu		<<"\t"<< setprecision(3) << fixed <<(float)IsoPassMu/(float)TightPassMu<<std::endl
  	   <<"Pt and Eta Pass ele: \t"	<<EtaPassMu		<<"\t"<< setprecision(3) << fixed <<(float)EtaPassMu/(float)IsoPassMu<<std::endl
	   <<"Muon eff: \t"		<<cutEff[0]		<<"\t"<< setprecision(3) << fixed <<(float)cutEff[0]/(float)EtaPassMu<<std::endl
	   <<"met eff: \t"		<<cutEff[1]		<<"\t"<< setprecision(3) << fixed <<(float)cutEff[1]/(float)cutEff[0]<<std::endl
	   <<"Jet pt & eta: \t"		<<JetPtEtaPass		<<"\t"<< setprecision(3) << fixed <<(float)JetPtEtaPass/(float)cutEff[1]<<std::endl
	   <<"Jet loose ID: \t"		<<JetLoosePass		<<"\t"<< setprecision(3) << fixed <<(float)JetLoosePass/(float)JetPtEtaPass<<std::endl
	   <<"Jet Cleaning: \t"		<<JetLepCleanPass	<<"\t"<< setprecision(3) << fixed <<(float)JetLepCleanPass/(float)JetLoosePass<<std::endl
	   <<"Events >=4Jet: \t"	<<TotalAK4Jets_MoreThan4<<"\t"<< setprecision(3) << fixed <<(float)TotalAK4Jets_MoreThan4/(float)JetLepCleanPass<<std::endl
	   <<"VBF Jet found:  \t"	<<cutEff[3]		<<"\t"<< setprecision(3) << fixed <<(float)cutEff[3]/(float)TotalAK4Jets_MoreThan4<<std::endl
	   <<"VBF Jet No B-tagged:  \t"	<<cutEff[6]		<<"\t"<< setprecision(3) << fixed <<(float)cutEff[6]/(float)cutEff[3]<<std::endl
	   <<"Wjet found: \t"		<<cutEff[2]		<<"\t"<< setprecision(3) << fixed <<(float)cutEff[2]/(float)cutEff[6]<<std::endl;
  }
  std::cout<<"===========================	OUTPUTS	=========================="<<std::endl;
	 

  std::cout<<"===========================	Short OUTPUTS	=========================="<<std::endl;
  std::cout<<"total entries: \t"<<ReducedTree->fChain->GetEntries()<<std::endl;
  if( strcmp(leptonName.c_str(),"el")==0)
  {
  std::cout<<"Electron eff: \t"		<<cutEff[0]		<<"\t"<< setprecision(3) << fixed <<(float)cutEff[0]/(float)ReducedTree->fChain->GetEntries()<<std::endl
	   <<"met eff:   \t "		<<cutEff[1]		<<"\t"<< setprecision(3) << fixed <<(float)cutEff[1]/(float)cutEff[0]<<std::endl
	   <<"Jet eff: \t"		<<JetLepCleanPass	<<"\t"<< setprecision(3) << fixed <<(float)JetLepCleanPass/(float)cutEff[1]<<std::endl
	   <<"Events >=4Jet: \t"	<<TotalAK4Jets_MoreThan4<<"\t"<< setprecision(3) << fixed <<(float)TotalAK4Jets_MoreThan4/(float)JetLepCleanPass<<std::endl
	   <<"VBF Jet found:  \t"	<<cutEff[3]		<<"\t"<< setprecision(3) << fixed <<(float)cutEff[3]/(float)TotalAK4Jets_MoreThan4<<std::endl
	   <<"VBF Jet No B-tagged:\t  "	<<cutEff[6]		<<"\t"<< setprecision(3) << fixed <<(float)cutEff[6]/(float)cutEff[3]<<std::endl
	   <<"Wjet found: \t"		<<cutEff[2]		<<"\t"<< setprecision(3) << fixed <<(float)cutEff[2]/(float)cutEff[6]<<std::endl;
  }
  if( strcmp(leptonName.c_str(),"mu")==0)
  {
  std::cout<<"Muon eff: \t"		<<cutEff[0]		<<"\t"<<  setprecision(3) << fixed <<(float)cutEff[0]/(float)ReducedTree->fChain->GetEntries()<<std::endl
	   <<"met eff: \t"		<<cutEff[1]		<<"\t"<<  setprecision(3) << fixed <<(float)cutEff[1]/(float)cutEff[0]<<std::endl
	   <<"Jet eff: \t"		<<JetLepCleanPass	<<"\t"<<  setprecision(3) << fixed <<(float)JetLepCleanPass/(float)cutEff[1]<<std::endl
	   <<"Events >=4Jet: \t"	<<TotalAK4Jets_MoreThan4<<"\t"<<  setprecision(3) << fixed <<(float)TotalAK4Jets_MoreThan4/(float)JetLepCleanPass<<std::endl
	   <<"VBF Jet found:  \t"	<<cutEff[3]		<<"\t"<<  setprecision(3) << fixed <<(float)cutEff[3]/(float)TotalAK4Jets_MoreThan4<<std::endl
	   <<"VBF Jet No B-tagged:  \t"	<<cutEff[6]		<<"\t"<<  setprecision(3) << fixed <<(float)cutEff[6]/(float)cutEff[3]<<std::endl
	   <<"Wjet found: \t"		<<cutEff[2]		<<"\t"<<  setprecision(3) << fixed <<(float)cutEff[2]/(float)cutEff[6]<<std::endl;
  }
  std::cout<<"===========================       Short OUTPUTS	=========================="<<std::endl;



  //--------close everything-------------
  ReducedTree->fChain->Delete();
  outTree->Write();
  outROOT->Close();
  t1.Write();
  t2.Write();
  t3.Write();

  return(0);
}

//////////////////////////////////
//// P A P E R   4 - V E C T O R   D E F I N I T I O N   O F   P H I   A N D   P H I 1
//////////////////////////////////
void computeAngles(TLorentzVector thep4H, TLorentzVector thep4Z1, TLorentzVector thep4M11, TLorentzVector thep4M12, TLorentzVector thep4Z2, TLorentzVector thep4M21, TLorentzVector thep4M22, double& costheta1, double& costheta2, double& Phi, double& costhetastar, double& Phi1){
    
    ///////////////////////////////////////////////
    // check for z1/z2 convention, redefine all 4 vectors with convention
    ///////////////////////////////////////////////	
    TLorentzVector p4H, p4Z1, p4M11, p4M12, p4Z2, p4M21, p4M22;
    p4H = thep4H;
    
    p4Z1 = thep4Z1; p4M11 = thep4M11; p4M12 = thep4M12;
    p4Z2 = thep4Z2; p4M21 = thep4M21; p4M22 = thep4M22;
    //// costhetastar
	TVector3 boostX = -(thep4H.BoostVector());
	TLorentzVector thep4Z1inXFrame( p4Z1 );
	TLorentzVector thep4Z2inXFrame( p4Z2 );	
	thep4Z1inXFrame.Boost( boostX );
	thep4Z2inXFrame.Boost( boostX );
	TVector3 theZ1X_p3 = TVector3( thep4Z1inXFrame.X(), thep4Z1inXFrame.Y(), thep4Z1inXFrame.Z() );
	TVector3 theZ2X_p3 = TVector3( thep4Z2inXFrame.X(), thep4Z2inXFrame.Y(), thep4Z2inXFrame.Z() );    
    costhetastar = theZ1X_p3.CosTheta();
    
    //// --------------------------- costheta1
    TVector3 boostV1 = -(thep4Z1.BoostVector());
    TLorentzVector p4M11_BV1( p4M11 );
	TLorentzVector p4M12_BV1( p4M12 );	
    TLorentzVector p4M21_BV1( p4M21 );
	TLorentzVector p4M22_BV1( p4M22 );
    p4M11_BV1.Boost( boostV1 );
	p4M12_BV1.Boost( boostV1 );
	p4M21_BV1.Boost( boostV1 );
	p4M22_BV1.Boost( boostV1 );
    
    TLorentzVector p4V2_BV1 = p4M21_BV1 + p4M22_BV1;
    //// costheta1
    costheta1 = -p4V2_BV1.Vect().Dot( p4M11_BV1.Vect() )/p4V2_BV1.Vect().Mag()/p4M11_BV1.Vect().Mag();
    
    //// --------------------------- costheta2
    TVector3 boostV2 = -(thep4Z2.BoostVector());
    TLorentzVector p4M11_BV2( p4M11 );
	TLorentzVector p4M12_BV2( p4M12 );	
    TLorentzVector p4M21_BV2( p4M21 );
	TLorentzVector p4M22_BV2( p4M22 );
    p4M11_BV2.Boost( boostV2 );
	p4M12_BV2.Boost( boostV2 );
	p4M21_BV2.Boost( boostV2 );
	p4M22_BV2.Boost( boostV2 );
    
    TLorentzVector p4V1_BV2 = p4M11_BV2 + p4M12_BV2;
    //// costheta2
    costheta2 = -p4V1_BV2.Vect().Dot( p4M21_BV2.Vect() )/p4V1_BV2.Vect().Mag()/p4M21_BV2.Vect().Mag();
    
    //// --------------------------- Phi and Phi1
    //    TVector3 boostX = -(thep4H.BoostVector());
    TLorentzVector p4M11_BX( p4M11 );
	TLorentzVector p4M12_BX( p4M12 );	
    TLorentzVector p4M21_BX( p4M21 );
	TLorentzVector p4M22_BX( p4M22 );	
    
	p4M11_BX.Boost( boostX );
	p4M12_BX.Boost( boostX );
	p4M21_BX.Boost( boostX );
	p4M22_BX.Boost( boostX );
    
    TVector3 tmp1 = p4M11_BX.Vect().Cross( p4M12_BX.Vect() );
    TVector3 tmp2 = p4M21_BX.Vect().Cross( p4M22_BX.Vect() );    
    
    TVector3 normal1_BX( tmp1.X()/tmp1.Mag(), tmp1.Y()/tmp1.Mag(), tmp1.Z()/tmp1.Mag() ); 
    TVector3 normal2_BX( tmp2.X()/tmp2.Mag(), tmp2.Y()/tmp2.Mag(), tmp2.Z()/tmp2.Mag() ); 
    
    //// Phi
    TLorentzVector p4Z1_BX = p4M11_BX + p4M12_BX;    
    double tmpSgnPhi = p4Z1_BX.Vect().Dot( normal1_BX.Cross( normal2_BX) );
    double sgnPhi = tmpSgnPhi/fabs(tmpSgnPhi);
    Phi = sgnPhi * acos( -1.*normal1_BX.Dot( normal2_BX) );
    
    
    //////////////
    
    TVector3 beamAxis(0,0,1);
    TVector3 tmp3 = (p4M11_BX + p4M12_BX).Vect();
    
    TVector3 p3V1_BX( tmp3.X()/tmp3.Mag(), tmp3.Y()/tmp3.Mag(), tmp3.Z()/tmp3.Mag() );
    TVector3 tmp4 = beamAxis.Cross( p3V1_BX );
    TVector3 normalSC_BX( tmp4.X()/tmp4.Mag(), tmp4.Y()/tmp4.Mag(), tmp4.Z()/tmp4.Mag() );
    
    //// Phi1
    double tmpSgnPhi1 = p4Z1_BX.Vect().Dot( normal1_BX.Cross( normalSC_BX) );
    double sgnPhi1 = tmpSgnPhi1/fabs(tmpSgnPhi1);    
    Phi1 = sgnPhi1 * acos( normal1_BX.Dot( normalSC_BX) );    
    
    //    std::cout << "extractAngles: " << std::endl;
    //    std::cout << "costhetastar = " << costhetastar << ", costheta1 = " << costheta1 << ", costheta2 = " << costheta2 << std::endl;
    //    std::cout << "Phi = " << Phi << ", Phi1 = " << Phi1 << std::endl;    
    
}
