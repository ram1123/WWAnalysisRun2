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

#include "../interface/setInputTree.h"
#include "../interface/setOutputTree.h"
#include "../interface/METzCalculator.h"
#include "../interface/METzCalculator_Run2.h"
#include "../interface/analysisUtils.h"
#include "../interface/readJSONFile.h"

using namespace std;
void computeAngles(TLorentzVector thep4H, TLorentzVector thep4Z1, TLorentzVector thep4M11, TLorentzVector thep4M12, 
					  TLorentzVector thep4Z2, 
		  double& costheta1,  double& costhetastar, double& Phi1);

float getPUPPIweight(float puppipt, float puppieta, TF1* puppisd_corrGEN, TF1* puppisd_corrRECO_cen, TF1* puppisd_corrRECO_for ){

  float genCorr = 1.;
  float recoCorr = 1.;
  float totalWeight = 1.;

  genCorr = puppisd_corrGEN->Eval( puppipt );

  if( fabs(puppieta) <= 1.3 ) recoCorr = puppisd_corrRECO_cen->Eval( puppipt );
  else recoCorr = puppisd_corrRECO_for->Eval( puppipt );

  totalWeight = genCorr * recoCorr;
  return totalWeight;

}

//*****PU WEIGHT***************

vector<double> generate_weights(TH1* data_npu_estimated, int isForSynch){
  // see https://github.com/cms-sw/cmssw/blob/CMSSW_8_0_X/SimGeneral/MixingModule/python/mix_2016_25ns_Moriond17MC_PoissonOOTPU_cfi.py#L25 ; copy and paste from there: 
  //
  const double npu_probs[75] = {  1.78653e-05 ,2.56602e-05 ,5.27857e-05 ,8.88954e-05 ,0.000109362 ,0.000140973 ,0.000240998 ,0.00071209 ,0.00130121 ,0.00245255 ,0.00502589 ,0.00919534 ,0.0146697 ,0.0204126 ,0.0267586 ,0.0337697 ,0.0401478 ,0.0450159 ,0.0490577 ,0.0524855 ,0.0548159 ,0.0559937 ,0.0554468 ,0.0537687 ,0.0512055 ,0.0476713 ,0.0435312 ,0.0393107 ,0.0349812 ,0.0307413 ,0.0272425 ,0.0237115 ,0.0208329 ,0.0182459 ,0.0160712 ,0.0142498 ,0.012804 ,0.011571 ,0.010547 ,0.00959489 ,0.00891718 ,0.00829292 ,0.0076195 ,0.0069806 ,0.0062025 ,0.00546581 ,0.00484127 ,0.00407168 ,0.00337681 ,0.00269893 ,0.00212473 ,0.00160208 ,0.00117884 ,0.000859662 ,0.000569085 ,0.000365431 ,0.000243565 ,0.00015688 ,9.88128e-05 ,6.53783e-05 ,3.73924e-05 ,2.61382e-05 ,2.0307e-05 ,1.73032e-05 ,1.435e-05 ,1.36486e-05 ,1.35555e-05 ,1.37491e-05 ,1.34255e-05 ,1.33987e-05 ,1.34061e-05 ,1.34211e-05 ,1.34177e-05 ,1.32959e-05 ,1.33287e-05 };
  
  if (isForSynch==0) { //OFFICIAL RECIPE
    vector<double> result(75);
    double s = 0.0;
    for(int npu=0; npu<75; ++npu){
      double npu_estimated = data_npu_estimated->GetBinContent(data_npu_estimated->GetXaxis()->FindBin(npu));
      result[npu] = npu_estimated / npu_probs[npu];
      //cout<<"npu_estimated = "<<npu_estimated<<"\t result[npu] = "<<result[npu]<<endl;
      s += npu_estimated;
    }
    // normalize weights such that the total sum of weights over thw whole sample is 1.0, i.e., sum_i  result[i] * npu_probs[i] should be 1.0 (!)
    for(int npu=0; npu<50; ++npu){
      result[npu] /= s;
    }
    return result;
  }

  else { //THIS IS FOR THE SYNCH ONLY. THIS IS NOT THE OFFICIAL RECIPE!
    vector<double> result(50);
    for(int npu=0; npu<50; ++npu){
      if (data_npu_estimated->GetBinContent(data_npu_estimated->GetXaxis()->FindBin(npu))==0.)
	result[npu] = 0.;
      else {
	double npu_estimated = data_npu_estimated->GetBinContent(data_npu_estimated->GetXaxis()->FindBin(npu));                              
	result[npu] = npu_estimated;
      }
    }
    return result;
  }

}


//*******MAIN*******************************************************************

int main (int argc, char** argv)
{ 

  int t0 = time(NULL);

  std::string inputFolder = argv[1];
  std::string outputFile = argv[2];
  int isMC = atoi(argv[3]);
  std::string leptonName = argv[4];
  std::string inputTreeName = argv[5];
  std::string inputFile = argv[6];
  std::string xSecWeight = argv[7];
  std::string numberOfEntries = argv[8];
  float genMass = atof(argv[9]);
  int applyTrigger = atoi(argv[10]);
  std::string jsonFileName = argv[11];
  int isLocal = atoi(argv[12]);
  
//cout<<"\n\n>>>>>>>>>>>>>>>		Ramkrishna >>>>>>> "<<endl;
  float weight = std::atof(xSecWeight.c_str())/std::atof(numberOfEntries.c_str());
  if (strcmp(leptonName.c_str(),"el")!=0 && strcmp(leptonName.c_str(),"mu")!=0) {
    std::cout<<"Error: wrong lepton category"<<std::endl;
    return(-1);
  }
  
 //std::string jsonFileName="/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions15/13TeV/Cert_246908-260627_13TeV_PromptReco_Collisions15_25ns_JSON.txt";
 jsonFileName="Cert_271036-284044_13TeV_23Sep2016ReReco_Collisions16_JSON.txt";
  std::map<int, std::vector<std::pair<int, int> > > jsonMap;
  jsonMap = readJSONFile(jsonFileName);
  std::cout<<"JSON file: "<<jsonFileName<<std::endl;

  // define map with events                                                                                                                                     
  std::map<std::pair<int,std::pair<int,int> >,int> eventsMap;
//
  //applyTrigger=false;
  std::cout<<"apply trigger: "<<applyTrigger<<std::endl;

  TLorentzVector W,W_puppi,LEP;
  TLorentzVector NU0,NU1,NU2,NU0_puppi,NU1_puppi,NU2_puppi;
  TLorentzVector NU0_jes_up, NU0_jes_dn;
  TLorentzVector JET, JET_PuppiAK8, AK4;
  TLorentzVector JET_jes_up, JET_jes_dn, JET_PuppiAK8_jes_up, JET_PuppiAK8_jes_dn;
  TLorentzVector AK4_JET1,AK4_JET2;
  TLorentzVector AK4_JET1_jes_up, AK4_JET1_jes_dn;
  TLorentzVector AK4_JET2_jes_up, AK4_JET2_jes_dn;
  TLorentzVector PuppiAK4_JET1,PuppiAK4_JET2;
  TLorentzVector PuppiAK4_JET1_jes_up, PuppiAK4_JET1_jes_dn;
  TLorentzVector PuppiAK4_JET2_jes_up, PuppiAK4_JET2_jes_dn;
  TLorentzVector VBF1,VBF2,TOT;
  TLorentzVector ELE,MU;

  
  std::vector<TLorentzVector> tightMuon;
  std::vector<TLorentzVector> looseMuon;
  std::vector<TLorentzVector> tightEle;
  std::vector<TLorentzVector> looseEle;

  int ok=0, total=0;
	
  setInputTree *ReducedTree = new setInputTree (inputTreeName.c_str());
  ReducedTree->Init();


  char command1[3000];
  //cout<<"Ramkrishna "<<endl;
  //exit(0);
  sprintf(command1, "xrd eoscms dirlist %s/%s/  | awk '{print \"root://eoscms.cern.ch/\"$5}' > listTemp_%s.txt", (inputFolder).c_str(), (inputFile).c_str(), outputFile.c_str());
  std::cout<<command1<<std::endl;
  system(command1);
  char list1[2000];
  sprintf (list1, "listTemp_%s.txt", outputFile.c_str());
  ifstream rootList (list1);

  int fileCounter=0;
  Long64_t totalEntries=0;
  
  if (isLocal==1) {
    ReducedTree->fChain->Add("ReducedSelection.root");
    cout<<" ====> Found local file."<<endl;
  }
  else {
    while (!rootList.eof())
    {
      char iRun_tW[700];
      rootList >> iRun_tW;
      if(!rootList.good())break;
      //cout<<"===> "<<iRun_tW<<endl;
      ReducedTree->fChain->Add(iRun_tW);
      fileCounter++;
    }
  }
  
  std::cout<<"number of files found: "<<fileCounter<<std::endl;
  totalEntries=ReducedTree->fChain->GetEntries();
  std::cout<<"total entries: "<<totalEntries<<std::endl;
  //exit(0);	//Just to see total number of entries
  //cout<<"\n\n>>>>>>>>>>>>>>>		Ramkrishna >>>>>>> "<<endl;
  
  char command3[300];
  sprintf(command3, "rm listTemp_%s.txt", outputFile.c_str());
  system(command3);

  int cutEff[20]={0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
  
  
  //--------pile up file -----------------
  std::vector<double> weights_pu1; //these are made with our recipe
  TFile* pileupFile1 = TFile::Open("PileUpData2016_23Sep2016ReReco_69200ub.root");
  TH1F* pileupHisto1 = (TH1F*)pileupFile1->Get("pileup");  
  weights_pu1 = generate_weights(pileupHisto1,0);
  pileupFile1->Close();
  
  std::vector<double> weights_pu2; //these are made with the official recipe
  TFile* pileupFile2 = TFile::Open("PUxSynch.root");  
  TH1F* pileupHisto2 = (TH1F*)pileupFile2->Get("puweights");
  weights_pu2 = generate_weights(pileupHisto2,1);
  pileupFile2->Close();


  //puppi corrections
  TFile* file = TFile::Open( "puppiJecCorr.root" );
  TF1* Puppisd_corrGEN = (TF1*)file->Get("puppiJECcorr_gen");
  TF1* Puppisd_corrRECO_cen = (TF1*)file->Get("puppiJECcorr_reco_0eta1v3");
  TF1* Puppisd_corrRECO_for = (TF1*)file->Get("puppiJECcorr_reco_1v3eta2v5");

  //---------output tree----------------
  TFile* outROOT = TFile::Open((outputFile+(".root")).c_str(),"recreate");
  outROOT->cd();
  TTree* outTree = new TTree("otree", "otree");
  setOutputTree* WWTree = new setOutputTree(outTree);

  //float top_NNLO_weight[2];

  std::ifstream badEventsFile;
  std::multimap<int,int> badEventsList;
  
  int run, lumi, evt;
  if (isMC==0) {
    if (strcmp(leptonName.c_str(),"el")==0)
      badEventsFile.open("SingleElectron_csc2015.txt");
    else
      badEventsFile.open("SingleMuon_csc2015.txt");      
    while(!badEventsFile.eof()) 
      {
        badEventsFile >> run >> lumi >> evt;
        badEventsList.insert(std::pair<int,int>(run,evt));
      }      
  }
  badEventsFile.close();

  int nNegEvents=0; 
  int nEvents=0;	nEvents=ReducedTree->fChain->GetEntries();
  Long64_t jentry2=0;

  
  //---------start loop on events------------
  std::cout << "---------start loop on events------------" << std::endl;
  //for (Long64_t jentry=0; jentry<ReducedTree->fChain->GetEntries();jentry++,jentry2++)
  for (Long64_t jentry=0; jentry<1198;jentry++,jentry2++)
  {
    
    Long64_t iEntry = ReducedTree->LoadTree(jentry);
    if (iEntry < 0) break;
    ReducedTree->fChain->GetEntry(jentry);
    
    tightMuon.clear();
    tightEle.clear();
    looseMuon.clear();
    looseEle.clear();

    if (jentry2%100 == 0) std::cout << "read entry: " << jentry2 <<"/"<<totalEntries<<std:: endl;
    
    //*********************************
    // JSON FILE AND DUPLIACTES IN DATA
    WWTree->run   = ReducedTree->runNum;
    WWTree->event = ReducedTree->evtNum;
    WWTree->lumi  = ReducedTree->lumiSec;
    
    /*     
    bool skipEvent = false;
    if( isMC==0 ) //apply json file
    {
      if(AcceptEventByRunAndLumiSection(ReducedTree->runNum,ReducedTree->lumiSec,jsonMap) == false) skipEvent = true;                           
      std::pair<int,Long64_t> eventLSandID(ReducedTree->lumiSec,ReducedTree->evtNum);            
      std::pair<int,std::pair<int,Long64_t> > eventRUNandLSandID(ReducedTree->runNum,eventLSandID); 
      if( eventsMap[eventRUNandLSandID] == 1 ) skipEvent = true;                                             
      else eventsMap[eventRUNandLSandID] = 1;                                                          
    }
    if( skipEvent == true ) continue;    
   */
    
    WWTree->initializeVariables(); //initialize all variables
   /*  
    //    if (ReducedTree->passFilterHBHELooseRerun == 0) continue;
    if (ReducedTree->passFilterHBHE == 0) continue;
    if (ReducedTree->passFilterHBHEIso == 0) continue;
    //    if (ReducedTree->passFilterCSCHalo == 0) continue;
    if (ReducedTree->passFilterGoodVtx == 0) continue;
    if (ReducedTree->passFilterEEBadSC == 0) continue;
    if (ReducedTree->passFilterEcalDeadCellTriggerPrimitive == 0) continue;
    if (ReducedTree->passFilterGlobalTightHalo2016 == 0) continue;
    if (ReducedTree->passFilterBadChCand == 0) continue;
    if (ReducedTree->passFilterBadPFMuon == 0) continue;
    */

    WWTree->issignal = 0;
    WWTree->wSampleWeight = weight; //xsec/numberOfEntries
    WWTree->eff_and_pu_Weight = 1.; //temporary value
    WWTree->eff_and_pu_Weight_2 = 1.; //temporary value
    WWTree->eff_and_pu_Weight_3 = 1.; //temporary value
    WWTree->top1_NNLO_Weight = 1.;
    WWTree->top2_NNLO_Weight = 1.;
    WWTree->trig_eff_Weight = 1.;

    if (ReducedTree->weight>0)
      WWTree->genWeight=1.;
    else if (ReducedTree->weight<0) {
      WWTree->genWeight=-1.;
      nNegEvents++;
    }
    // std::cout<<ReducedTree->weight<<" "<<WWTree->genWeight<<std::endl;
    // WWTree->genWeight = ReducedTree->weight;
    
    //PILE-UP WEIGHT
    /*
    if (isMC==1) {
      if(ReducedTree->PV_<int(weights_pu1.size())){
	WWTree->eff_and_pu_Weight = weights_pu1[ReducedTree->npT]; //official pu recipe
	//cout<<"**** ReducedTree->npT = "<<ReducedTree->npT<<endl;
	//cout<<"===> weights_pu1 = "<<weights_pu1[ReducedTree->npT]<<endl;
      }
      else{ //should not happen as we have a weight for all simulated n_pu multiplicities!
	std::cout<<"Warning! n_pu too big"<<std::endl;
	// throw logic_error("n_pu too big");
        WWTree->eff_and_pu_Weight = 0.;
      }
      
      if(ReducedTree->PV_<int(weights_pu2.size())){
	WWTree->eff_and_pu_Weight_2 = weights_pu2[ReducedTree->PV_]; //our pu recipe
      }
      else{ //should not happen as we have a weight for all simulated n_pu multiplicities!
	std::cout<<"Warning! n_pu too big"<<std::endl;
	// throw logic_error("n_pu too big");
        WWTree->eff_and_pu_Weight_2 = 0.;
      }    
      WWTree->eff_and_pu_Weight_3 = ReducedTree->PUWeight; //our pu recipe
    }*/
    
    //save event variables
    WWTree->run   = ReducedTree->runNum;
    WWTree->event = ReducedTree->evtNum;
    WWTree->lumi  = ReducedTree->lumiSec;
    
    WWTree->nPV = ReducedTree->PV_;
    
    cutEff[0]++;
    
    /////////////////THE SELECTED LEPTON
    int nTightLepton=0;
    if (strcmp(leptonName.c_str(),"el")==0) { //electrons
      int passTrigger=0;
      float tempPt=0.;

      for (int i=0; i<ReducedTree->Electron_; i++) {
//	if (applyTrigger==1)
//	  for ( int t=0; t<ReducedTree->TriggerProducerTriggerNames->size(); t++)
//	    if(TString(ReducedTree->TriggerProducerTriggerNames->at(t)).Contains("HLT_Ele27_eta2p1_WPLoose_Gsf_v") ||
//	       TString(ReducedTree->TriggerProducerTriggerNames->at(t)).Contains("HLT_Ele105_CaloIdVT_GsfTrkIdT_v") )
//	      if (ReducedTree->TriggerProducerTriggerPass->at(t)==1) passTrigger=1; //trigger
//	if (passTrigger==0 && applyTrigger==1) continue;
	// if (ReducedTree->TriggerProducerTriggerPass->at(0)==0) continue; //trigger
	//if (ReducedTree->Electrons_isTight[i]==false) continue;       
        if (ReducedTree->Electron_pt[i]<=45) continue;
        if (fabs(ReducedTree->Electron_eta[i])>=2.1) continue;
	// if (fabs(ReducedTree->Electron_eta[i])>=2.5) continue; //this is already in the HEEP requirement
	if (ReducedTree->Electron_pt[i]<tempPt) continue;
	ELE.SetPtEtaPhiE(ReducedTree->Electron_pt[i],ReducedTree->Electron_eta[i],ReducedTree->Electron_phi[i],ReducedTree->Electron_ecalEnergy[i]);
	tightEle.push_back(ELE);
	WWTree->l_pt  = ReducedTree->Electron_pt[i];
	WWTree->l_eta = ReducedTree->Electron_eta[i];
	WWTree->l_phi = ReducedTree->Electron_phi[i];	
	WWTree->l_e= ReducedTree->Electron_ecalEnergy[i];	
	WWTree->l_charge= ReducedTree->Electron_q[i];
//	std::cout<<"ELE: pt = "<< WWTree->l_pt <<"\tEta = "<< WWTree->l_eta <<"\tPhi = "<< WWTree->l_phi <<"\tE = "<< WWTree->l_e <<"\tCharge = "<< WWTree->l_charge <<endl;
	tempPt = WWTree->l_pt;
	nTightLepton++;
      }
    } //end electrons
    else if (strcmp(leptonName.c_str(),"mu")==0) { //Muon
      int passTrigger=0;
      float tempPt=0.;
      for (int i=0; i<ReducedTree->Muon_; i++) {
	//if (applyTrigger==1)
	//  for ( int t=0; t<ReducedTree->TriggerProducerTriggerNames->size(); t++)
	//    if(TString(ReducedTree->TriggerProducerTriggerNames->at(t)).Contains("HLT_IsoMu27"))
	//      if (ReducedTree->TriggerProducerTriggerPass->at(t)==1) passTrigger=1; //trigger
	//if (passTrigger==0 && applyTrigger==1) continue;
	//// if (ReducedTree->TriggerProducerTriggerPass->at(1)==0) continue; //trigger
	//if (ReducedTree->Muon_isTight[i]==false) continue;
	// if (ReducedTree->Muon_isPFMuon[i]==false) continue; //not in the synch ntuple!!
        if ((ReducedTree->Muon_trkIso[i]/ReducedTree->Muon_pt[i])>=0.1) continue;
        if (ReducedTree->Muon_pt[i]<40) continue;
        if (fabs(ReducedTree->Muon_eta[i])>=2.1) continue;
	MU.SetPtEtaPhiE(ReducedTree->Muon_pt[i],ReducedTree->Muon_eta[i],ReducedTree->Muon_phi[i],ReducedTree->Muon_pt[i]);
	tightMuon.push_back(MU);
	if (ReducedTree->Muon_pt[i]<tempPt) continue;
	WWTree->l_pt  = ReducedTree->Muon_pt[i];
	WWTree->l_eta = ReducedTree->Muon_eta[i];
	WWTree->l_phi = ReducedTree->Muon_phi[i];
	WWTree->l_e = ReducedTree->Muon_pt[i];
	WWTree->l_charge= ReducedTree->Muon_q[i];
	std::cout<<"ELE: pt = "<< WWTree->l_pt <<"\tEta = "<< WWTree->l_eta <<"\tPhi = "<< WWTree->l_phi <<"\tE = "<< WWTree->l_e <<"\tCharge = "<< WWTree->l_charge <<endl;
	tempPt = WWTree->l_pt;
	nTightLepton++;
      }
    } //end Muon
    if (nTightLepton==0) continue; //no leptons with required ID
    cutEff[1]++;
    
    LEP.SetPtEtaPhiE(WWTree->l_pt,WWTree->l_eta,WWTree->l_phi,WWTree->l_e);
    
    
    //trigger SF for muon (from B2G-15-005, these are for the HLT_IsoMu20 trigger,
    //see https://indico.cern.ch/event/462268/contribution/9/attachments/1188638/1724574/2015.11.17_MuonPOG_SingleMuTrigEff_SF_KPLee_v2.pdf
    if (strcmp(leptonName.c_str(),"mu")==0 && isMC==1) {
      if (WWTree->run >= 257820) {
	if ( fabs(WWTree->l_eta) < 0.9) {
	  if      ( WWTree->l_pt<50)  WWTree->trig_eff_Weight = 1.0000;
	  else if ( WWTree->l_pt<60)  WWTree->trig_eff_Weight = 0.9981;
	  else                        WWTree->trig_eff_Weight = 0.9908;
	}	
	else if ( fabs(WWTree->l_eta) >= 0.9 && fabs(WWTree->l_eta)<1.2) {
	  if      ( WWTree->l_pt<50)  WWTree->trig_eff_Weight = 0.9992;
	  else if ( WWTree->l_pt<60)  WWTree->trig_eff_Weight = 1.0028;
	  else                        WWTree->trig_eff_Weight = 0.9847;
	}	
	else if ( fabs(WWTree->l_eta) >= 1.2 && fabs(WWTree->l_eta)<2.1) {
	  if      ( WWTree->l_pt<50)  WWTree->trig_eff_Weight = 0.9951;
	  else if ( WWTree->l_pt<60)  WWTree->trig_eff_Weight = 0.9940;
	  else                        WWTree->trig_eff_Weight = 0.9931;
	}
      }
      else {
	if ( fabs(WWTree->l_eta) < 0.9) {
	  if      ( WWTree->l_pt<50)  WWTree->trig_eff_Weight = 0.9770;
	  else if ( WWTree->l_pt<60)  WWTree->trig_eff_Weight = 0.9848;
	  else                        WWTree->trig_eff_Weight = 0.9763;
	}	
	else if ( fabs(WWTree->l_eta) >= 0.9 && fabs(WWTree->l_eta)<1.2) {
	  if      ( WWTree->l_pt<50)  WWTree->trig_eff_Weight = 0.9810;
	  else if ( WWTree->l_pt<60)  WWTree->trig_eff_Weight = 0.9890;
	  else                        WWTree->trig_eff_Weight = 0.9807;
	}	
	else if ( fabs(WWTree->l_eta) >= 1.2 && fabs(WWTree->l_eta)<2.1) {
	  if      ( WWTree->l_pt<50)  WWTree->trig_eff_Weight = 0.9767;
	  else if ( WWTree->l_pt<60)  WWTree->trig_eff_Weight = 0.9788;
	  else                        WWTree->trig_eff_Weight = 0.9808;
	}
      }
    }
    
    //trigger SF for electron (from B2G-15-005)
    if (strcmp(leptonName.c_str(),"el")==0 && isMC==1) {
      if ( fabs(WWTree->l_eta) < 0.4) {
	if      ( WWTree->l_pt<50)   WWTree->trig_eff_Weight = 0.980;
	else if ( WWTree->l_pt<200)  WWTree->trig_eff_Weight = 0.998;
	else                         WWTree->trig_eff_Weight = 1.011;
      }	
      else if ( fabs(WWTree->l_eta) >= 0.4 && fabs(WWTree->l_eta)<1.4 ) {
	if      ( WWTree->l_pt<50)   WWTree->trig_eff_Weight = 1.005;
	else if ( WWTree->l_pt<200)  WWTree->trig_eff_Weight = 1.012;
	else                         WWTree->trig_eff_Weight = 1.012;
      }
      else if ( fabs(WWTree->l_eta) >= 1.4 && fabs(WWTree->l_eta)<1.6 ) {
	if      ( WWTree->l_pt<50)   WWTree->trig_eff_Weight = 0.998;
	else if ( WWTree->l_pt<200)  WWTree->trig_eff_Weight = 1.004;
	else                         WWTree->trig_eff_Weight = 0.947;
      }
      else {
	if      ( WWTree->l_pt<50)   WWTree->trig_eff_Weight = 0.958;
	else if ( WWTree->l_pt<200)  WWTree->trig_eff_Weight = 0.944;
	else                         WWTree->trig_eff_Weight = 0.950;
      }
    }
    
    //ID efficiency SF for electrons (from B2G-15-005)
    if (strcmp(leptonName.c_str(),"el")==0 && isMC==1) {
      if ( fabs(WWTree->l_eta) < 1.4) {
	if      ( WWTree->l_pt<50)  WWTree->id_eff_Weight = 0.975;
	else if ( WWTree->l_pt<100) WWTree->id_eff_Weight = 0.976;
	else if ( WWTree->l_pt<150) WWTree->id_eff_Weight = 0.977;
	else                        WWTree->id_eff_Weight = 1.008;
      }	
      else if (fabs(WWTree->l_eta) >= 1.4 && fabs(WWTree->l_eta)<1.6 ) {
	if      ( WWTree->l_pt<50)  WWTree->id_eff_Weight = 0.968;
	else if ( WWTree->l_pt<100) WWTree->id_eff_Weight = 0.996;
	else if ( WWTree->l_pt<150) WWTree->id_eff_Weight = 0.983;
	else                        WWTree->id_eff_Weight = 1.045;
      }
      else if (fabs(WWTree->l_eta) >= 1.6 && fabs(WWTree->l_eta)<2.1 ) {
	if      ( WWTree->l_pt<50)  WWTree->id_eff_Weight = 0.987;
	else if ( WWTree->l_pt<100) WWTree->id_eff_Weight = 0.981;
	else if ( WWTree->l_pt<150) WWTree->id_eff_Weight = 0.996;
	else                        WWTree->id_eff_Weight = 0.977;
      }
    }
    
    //ID efficiency SF for Muon (from B2G-15-005)
    if (strcmp(leptonName.c_str(),"mu")==0 && isMC==1) {
      if ( fabs(WWTree->l_eta) < 1.2) {
	if      ( WWTree->l_pt<50)  WWTree->id_eff_Weight = 0.994;
	else if ( WWTree->l_pt<100) WWTree->id_eff_Weight = 0.992;
	else if ( WWTree->l_pt<150) WWTree->id_eff_Weight = 0.985;
	else                        WWTree->id_eff_Weight = 0.984;
      }	
      else if (fabs(WWTree->l_eta) >= 1.2 && fabs(WWTree->l_eta)<2.1 ) {
	if      ( WWTree->l_pt<50)  WWTree->id_eff_Weight = 0.995;
	else if ( WWTree->l_pt<100) WWTree->id_eff_Weight = 0.993;
	else if ( WWTree->l_pt<150) WWTree->id_eff_Weight = 0.986;
	else                        WWTree->id_eff_Weight = 0.960;
      }
    }
    
    //VETO ADDITIONAL LEPTONS
    int nLooseLepton=0;
    for (int i=0; i<ReducedTree->Electron_; i++) {
      //if (ReducedTree->Electrons_isLoose[i]==false) continue; //NB: CHANGE TO VETO!!!
      if (ReducedTree->Electron_pt[i]<30) continue;       
      // if (fabs(ReducedTree->Electron_eta[i])>=2.5) continue;
      ELE.SetPtEtaPhiE(ReducedTree->Electron_pt[i],ReducedTree->Electron_eta[i],ReducedTree->Electron_phi[i],ReducedTree->Electron_ecalEnergy[i]);
      looseEle.push_back(ELE);      
      nLooseLepton++;
    }
    for (int i=0; i<ReducedTree->Muon_; i++) {
      //if (ReducedTree->Muon_isLoose[i]==false) continue;
      if ((ReducedTree->Muon_trkIso[i]/ReducedTree->Muon_pt[i])>=0.1) continue;
      if (fabs(ReducedTree->Muon_eta[i])>=2.4) continue;
      if (ReducedTree->Muon_pt[i]<20) continue;
      MU.SetPtEtaPhiE(ReducedTree->Muon_pt[i],ReducedTree->Muon_eta[i],ReducedTree->Muon_phi[i],ReducedTree->Muon_pt[i]);
      looseMuon.push_back(MU);
      nLooseLepton++;
    }
    if (nLooseLepton!=1) continue; //no additional leptons
    cutEff[2]++;
    
    
    
    //////////////THE MET
    
    // //preselection on met
    // if (ReducedTree->pfMET < 30) continue;
    // cutEff[3]++;
    
    // Calculate Neutrino Pz using all the possible choices : 
    // type0 -> if real roots, pick the one nearest to the lepton Pz except when the Pz so chosen
    //               is greater than 300 GeV in which case pick the most central root.
    // type1 -> type = 1: if real roots, choose the one closest to the lepton Pz if complex roots, use only the real part.
    //          type = 2: if real roots, choose the most central solution. if complex roots, use only the real part. 
    //          type = 3: if real roots, pick the largest value of the cosine*

    float Wmass = 80.385;

    TLorentzVector W_Met, W_Met_jes_up, W_Met_jes_dn;
    
    W_Met.SetPxPyPzE(ReducedTree->pfMET * TMath::Cos(ReducedTree->pfMETphi), ReducedTree->pfMET * TMath::Sin(ReducedTree->pfMETphi), 0., sqrt(ReducedTree->pfMET*ReducedTree->pfMET));
    //W_Met_jes_up.SetPxPyPzE(ReducedTree->METPtUp * TMath::Cos(ReducedTree->METPhiUp), ReducedTree->METPtUp * TMath::Sin(ReducedTree->METPhiUp), 0., sqrt(ReducedTree->METPtUp*ReducedTree->METPtUp));
    //W_Met_jes_dn.SetPxPyPzE(ReducedTree->METPtDown * TMath::Cos(ReducedTree->METPhiDown), ReducedTree->METPtDown * TMath::Sin(ReducedTree->METPhiDown), 0., sqrt(ReducedTree->METPtDown*ReducedTree->METPtDown));
    
    // if(LEP.Pt()<=0 || W_Met.Pt() <= 0 ){ std::cerr<<" Negative Lepton - Neutrino Pt "<<std::endl; continue ; }
    
    // type0 calculation of neutrino pZ
    METzCalculator NeutrinoPz_type0;
    METzCalculator NeutrinoPz_type0_jes_up;
    METzCalculator NeutrinoPz_type0_jes_dn;
    METzCalculator_Run2 NeutrinoPz_run2;
    NeutrinoPz_type0.SetMET(W_Met);
    NeutrinoPz_type0.SetLepton(LEP);
    NeutrinoPz_type0.SetLeptonType(leptonName.c_str());
    
    NeutrinoPz_type0_jes_up.SetMET(W_Met_jes_up);
    NeutrinoPz_type0_jes_up.SetLepton(LEP);
    NeutrinoPz_type0_jes_up.SetLeptonType(leptonName.c_str());
    
    NeutrinoPz_type0_jes_dn.SetMET(W_Met_jes_dn);
    NeutrinoPz_type0_jes_dn.SetLepton(LEP);
    NeutrinoPz_type0_jes_dn.SetLeptonType(leptonName.c_str());
    
    NeutrinoPz_run2.SetMET(W_Met);
    NeutrinoPz_run2.SetLepton(LEP);
    NeutrinoPz_run2.SetLeptonType(leptonName.c_str());
    
    double pz1_type0 = NeutrinoPz_type0.Calculate(); // Default one -> according to type0
    //double pz2_type0 = NeutrinoPz_type0.getOther();  // Default one
    
    double pz1_run2 = NeutrinoPz_run2.Calculate();
    
    double pz1_type0_jes_up = NeutrinoPz_type0_jes_up.Calculate(); // Default one -> according to type0
    double pz1_type0_jes_dn = NeutrinoPz_type0_jes_dn.Calculate(); // Default one -> according to type0
    
    // don't touch the neutrino pT
    TLorentzVector W_neutrino_type0_met; 
    W_neutrino_type0_met.SetPxPyPzE(W_Met.Px(),W_Met.Py(),pz1_type0,sqrt(W_Met.Pt()*W_Met.Pt()+pz1_type0*pz1_type0));
    
    // change the neutrino pT in case of complex solution in order to make it real
    TLorentzVector W_neutrino_type0; 
    W_neutrino_type0.SetPxPyPzE(W_Met.Px(),W_Met.Py(),pz1_type0,sqrt(W_Met.Pt()*W_Met.Pt()+pz1_type0*pz1_type0));

    if (NeutrinoPz_type0.IsComplex()) {// if this is a complex, change MET
      double nu_pt1 = NeutrinoPz_type0.getPtneutrino(1);
      double nu_pt2 = NeutrinoPz_type0.getPtneutrino(2);
      TLorentzVector W_neutrino_1;
      W_neutrino_1.SetPxPyPzE(nu_pt1 * TMath::Cos(ReducedTree->pfMETphi), nu_pt1 * TMath::Sin(ReducedTree->pfMETphi), pz1_type0, sqrt(nu_pt1*nu_pt1 + pz1_type0*pz1_type0) );
      TLorentzVector W_neutrino_2;
      W_neutrino_2.SetPxPyPzE(nu_pt2 * TMath::Cos(ReducedTree->pfMETphi), nu_pt2 * TMath::Sin(ReducedTree->pfMETphi), pz1_type0, sqrt(nu_pt2*nu_pt2 + pz1_type0*pz1_type0) );
      
      if ( fabs((LEP+W_neutrino_1).M()-Wmass) < fabs((LEP+W_neutrino_2).M()-Wmass) ) W_neutrino_type0 = W_neutrino_1;
      else W_neutrino_type0 = W_neutrino_2;
    }
    
    // type2 calculation of neutrino pZ
    METzCalculator NeutrinoPz_type2;
    NeutrinoPz_type2.SetMET(W_Met);
    NeutrinoPz_type2.SetLepton(LEP);
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
      W_neutrino_1.SetPxPyPzE(nu_pt1 * TMath::Cos(ReducedTree->pfMETphi), nu_pt1 * TMath::Sin(ReducedTree->pfMETphi), pz1_type2, sqrt(nu_pt1*nu_pt1 + pz1_type2*pz1_type2) );
      TLorentzVector W_neutrino_2;
      W_neutrino_2.SetPxPyPzE(nu_pt2 * TMath::Cos(ReducedTree->pfMETphi), nu_pt2 * TMath::Sin(ReducedTree->pfMETphi), pz1_type2, sqrt(nu_pt2*nu_pt2 + pz1_type2*pz1_type2) );
      
      if ( fabs((LEP+W_neutrino_1).M()-Wmass) < fabs((LEP+W_neutrino_2).M()-Wmass) ) W_neutrino_type2 = W_neutrino_1;
      else W_neutrino_type2 = W_neutrino_2;
    }
    
    WWTree->pfMET = sqrt(ReducedTree->pfMET*ReducedTree->pfMET);
    //WWTree->pfMET_jes_up = sqrt(ReducedTree->METPtUp*ReducedTree->METPtUp);
    //WWTree->pfMET_jes_dn = sqrt(ReducedTree->METPtDown*ReducedTree->METPtDown);
    WWTree->pfMET_Phi = ReducedTree->pfMETphi;
    WWTree->nu_pz_type0 = pz1_type0;
    WWTree->nu_pz_type2 = pz1_type2;
    WWTree->nu_pz_run2 = pz1_run2;
    WWTree->nu_pz_isre = 1-NeutrinoPz_run2.IsComplex();
    WWTree->nu_pz_run2_oth = NeutrinoPz_run2.getOther();
    WWTree->nu_pz_run2_type = NeutrinoPz_run2.getType();
    
    
    /////////////////THE LEPTONIC W
    
    NU0.SetPxPyPzE(ReducedTree->pfMET*TMath::Cos(ReducedTree->pfMETphi),ReducedTree->pfMET*TMath::Sin(ReducedTree->pfMETphi),WWTree->nu_pz_type0,TMath::Sqrt(WWTree->pfMET*WWTree->pfMET+WWTree->nu_pz_type0*WWTree->nu_pz_type0));
    //NU0_jes_up.SetPxPyPzE(ReducedTree->METPtUp*TMath::Cos(ReducedTree->METPhiUp),ReducedTree->METPtUp*TMath::Sin(ReducedTree->METPhiUp),pz1_type0_jes_up,TMath::Sqrt(WWTree->pfMET_jes_up*WWTree->pfMET_jes_up+pz1_type0_jes_up*pz1_type0_jes_up));
    //NU0_jes_dn.SetPxPyPzE(ReducedTree->METPtDown*TMath::Cos(ReducedTree->METPhiDown),ReducedTree->METPtDown*TMath::Sin(ReducedTree->METPhiDown),pz1_type0_jes_dn,TMath::Sqrt(WWTree->pfMET_jes_dn*WWTree->pfMET_jes_dn+pz1_type0_jes_dn*pz1_type0_jes_dn));
    
    NU2.SetPxPyPzE(ReducedTree->pfMET*TMath::Cos(ReducedTree->pfMETphi),ReducedTree->pfMET*TMath::Sin(ReducedTree->pfMETphi),WWTree->nu_pz_type2,TMath::Sqrt(WWTree->pfMET*WWTree->pfMET+WWTree->nu_pz_type2*WWTree->nu_pz_type2));
    NU1.SetPxPyPzE(ReducedTree->pfMET*TMath::Cos(ReducedTree->pfMETphi),ReducedTree->pfMET*TMath::Sin(ReducedTree->pfMETphi),WWTree->nu_pz_run2,TMath::Sqrt(WWTree->pfMET*WWTree->pfMET+WWTree->nu_pz_run2*WWTree->nu_pz_run2));

    
    W = LEP + NU2;
    
    WWTree->v_pt = W.Pt();
    WWTree->v_eta = W.Eta();
    WWTree->v_phi = W.Phi();
    WWTree->v_mt = TMath::Sqrt(2*LEP.Et()*NU2.Et()*(1-TMath::Cos(LEP.DeltaPhi(NU2))));
    WWTree->v_mass = W.M();



    //////////////THE PUPPI MET
    
    //preselection on met
    if (ReducedTree->pfMET < 30 && ReducedTree->puppET < 30) continue;
    cutEff[3]++;
    
    // Calculate Neutrino Pz using all the possible choices : 
    // type0 -> if real roots, pick the one nearest to the lepton Pz except when the Pz so chosen
    //               is greater than 300 GeV in which case pick the most central root.
    // type1 -> type = 1: if real roots, choose the one closest to the lepton Pz if complex roots, use only the real part.
    //          type = 2: if real roots, choose the most central solution. if complex roots, use only the real part. 
    //          type = 3: if real roots, pick the largest value of the cosine*
    
    W_Met.SetPxPyPzE(ReducedTree->puppET * TMath::Cos(ReducedTree->puppETphi), ReducedTree->puppET * TMath::Sin(ReducedTree->puppETphi), 0., sqrt(ReducedTree->puppET*ReducedTree->puppET));
    //W_Met_jes_up.SetPxPyPzE(ReducedTree->METpuppiPtUp * TMath::Cos(ReducedTree->METpuppiPhiUp), ReducedTree->METpuppiPtUp * TMath::Sin(ReducedTree->METpuppiPhiUp), 0., sqrt(ReducedTree->METpuppiPtUp*ReducedTree->METpuppiPtUp));
    //W_Met_jes_dn.SetPxPyPzE(ReducedTree->METpuppiPtDown * TMath::Cos(ReducedTree->METpuppiPhiDown), ReducedTree->METpuppiPtDown * TMath::Sin(ReducedTree->METpuppiPhiDown), 0., sqrt(ReducedTree->METpuppiPtDown*ReducedTree->METpuppiPtDown));
    
    if(LEP.Pt()<=0 || W_Met.Pt() <= 0 ){ std::cerr<<" Negative Lepton - Neutrino Pt "<<std::endl; continue ; }
    cutEff[4]++;
    
    // type0 calculation of neutrino pZ
    NeutrinoPz_type0.SetMET(W_Met);
    NeutrinoPz_type0.SetLepton(LEP);
    NeutrinoPz_type0.SetLeptonType(leptonName.c_str());
    
    //NeutrinoPz_type0_jes_up.SetMET(W_Met_jes_up);
    //NeutrinoPz_type0_jes_up.SetLepton(LEP);
    //NeutrinoPz_type0_jes_up.SetLeptonType(leptonName.c_str());
    //
    //NeutrinoPz_type0_jes_dn.SetMET(W_Met_jes_dn);
    //NeutrinoPz_type0_jes_dn.SetLepton(LEP);
    //NeutrinoPz_type0_jes_dn.SetLeptonType(leptonName.c_str());
    
    NeutrinoPz_run2.SetMET(W_Met);
    NeutrinoPz_run2.SetLepton(LEP);
    NeutrinoPz_run2.SetLeptonType(leptonName.c_str());
    
     pz1_type0 = NeutrinoPz_type0.Calculate(); // Default one -> according to type0
    // pz2_type0 = NeutrinoPz_type0.getOther();  // Default one
    
     pz1_run2 = NeutrinoPz_run2.Calculate();
     //cout<<"pz1_run2 = "<<pz1_run2<<endl;
    
     //pz1_type0_jes_up = NeutrinoPz_type0_jes_up.Calculate(); // Default one -> according to type0
     //pz1_type0_jes_dn = NeutrinoPz_type0_jes_dn.Calculate(); // Default one -> according to type0
    
    // don't touch the neutrino pT
    W_neutrino_type0_met.SetPxPyPzE(W_Met.Px(),W_Met.Py(),pz1_type0,sqrt(W_Met.Pt()*W_Met.Pt()+pz1_type0*pz1_type0));
    
    // change the neutrino pT in case of complex solution in order to make it real
    W_neutrino_type0.SetPxPyPzE(W_Met.Px(),W_Met.Py(),pz1_type0,sqrt(W_Met.Pt()*W_Met.Pt()+pz1_type0*pz1_type0));

    if (NeutrinoPz_type0.IsComplex()) {// if this is a complex, change MET
      double nu_pt1 = NeutrinoPz_type0.getPtneutrino(1);
      double nu_pt2 = NeutrinoPz_type0.getPtneutrino(2);
      TLorentzVector W_neutrino_1;
      W_neutrino_1.SetPxPyPzE(nu_pt1 * TMath::Cos(ReducedTree->puppETphi), nu_pt1 * TMath::Sin(ReducedTree->puppETphi), pz1_type0, sqrt(nu_pt1*nu_pt1 + pz1_type0*pz1_type0) );
      TLorentzVector W_neutrino_2;
      W_neutrino_2.SetPxPyPzE(nu_pt2 * TMath::Cos(ReducedTree->puppETphi), nu_pt2 * TMath::Sin(ReducedTree->puppETphi), pz1_type0, sqrt(nu_pt2*nu_pt2 + pz1_type0*pz1_type0) );
      
      if ( fabs((LEP+W_neutrino_1).M()-Wmass) < fabs((LEP+W_neutrino_2).M()-Wmass) ) W_neutrino_type0 = W_neutrino_1;
      else W_neutrino_type0 = W_neutrino_2;
    }
    
    // type2 calculation of neutrino pZ
    NeutrinoPz_type2.SetMET(W_Met);
    NeutrinoPz_type2.SetLepton(LEP);
    NeutrinoPz_type2.SetLeptonType(leptonName.c_str());
    
    pz1_type2 = NeutrinoPz_type2.Calculate(2); // Default one -> according to type2
    //double pz2_type2 = NeutrinoPz_type2.getOther();   // Default one
    
    // don't touch the neutrino pT
    W_neutrino_type2_met.SetPxPyPzE(W_Met.Px(),W_Met.Py(),pz1_type2,sqrt(W_Met.Pt()*W_Met.Pt()+pz1_type2*pz1_type2));
    
    // change the neutrino pT in case of complex solution in order to make it real
    W_neutrino_type2.SetPxPyPzE(W_Met.Px(),W_Met.Py(),pz1_type2,sqrt(W_Met.Pt()*W_Met.Pt()+pz1_type2*pz1_type2));
    
    if (NeutrinoPz_type2.IsComplex()) {// if this is a complex, change MET
      double nu_pt1 = NeutrinoPz_type2.getPtneutrino(1);
      double nu_pt2 = NeutrinoPz_type2.getPtneutrino(2);
      TLorentzVector W_neutrino_1;
      W_neutrino_1.SetPxPyPzE(nu_pt1 * TMath::Cos(ReducedTree->puppETphi), nu_pt1 * TMath::Sin(ReducedTree->puppETphi), pz1_type2, sqrt(nu_pt1*nu_pt1 + pz1_type2*pz1_type2) );
      TLorentzVector W_neutrino_2;
      W_neutrino_2.SetPxPyPzE(nu_pt2 * TMath::Cos(ReducedTree->puppETphi), nu_pt2 * TMath::Sin(ReducedTree->puppETphi), pz1_type2, sqrt(nu_pt2*nu_pt2 + pz1_type2*pz1_type2) );
      
      if ( fabs((LEP+W_neutrino_1).M()-Wmass) < fabs((LEP+W_neutrino_2).M()-Wmass) ) W_neutrino_type2 = W_neutrino_1;
      else W_neutrino_type2 = W_neutrino_2;
    }
    
    WWTree->pfMETpuppi = sqrt(ReducedTree->puppET*ReducedTree->puppET);
    //WWTree->pfMETpuppi_jes_up = sqrt(ReducedTree->METpuppiPtUp*ReducedTree->METpuppiPtUp);
    //WWTree->pfMETpuppi_jes_dn = sqrt(ReducedTree->METpuppiPtDown*ReducedTree->METpuppiPtDown);
    WWTree->pfMETpuppi_Phi = ReducedTree->puppETphi;
    WWTree->nu_pz_type0 = pz1_type0;
    WWTree->nu_pz_type2 = pz1_type2;
    WWTree->nu_pz_run2 = pz1_run2;
    WWTree->nu_pz_isre = 1-NeutrinoPz_run2.IsComplex();
    WWTree->nu_pz_run2_oth = NeutrinoPz_run2.getOther();
    WWTree->nu_pz_run2_type = NeutrinoPz_run2.getType();

    
    /////////////////THE LEPTONIC W PUPPI
    
    NU0_puppi.SetPxPyPzE(ReducedTree->puppET*TMath::Cos(ReducedTree->puppETphi),ReducedTree->puppET*TMath::Sin(ReducedTree->puppETphi),WWTree->nu_pz_type0,TMath::Sqrt(WWTree->pfMETpuppi*WWTree->pfMETpuppi+WWTree->nu_pz_type0*WWTree->nu_pz_type0));
    //NU0_jes_up.SetPxPyPzE(ReducedTree->METpuppiPtUp*TMath::Cos(ReducedTree->METpuppiPhiUp),ReducedTree->METpuppiPtUp*TMath::Sin(ReducedTree->METpuppiPhiUp),pz1_type0_jes_up,TMath::Sqrt(WWTree->pfMETpuppi_jes_up*WWTree->pfMETpuppi_jes_up+pz1_type0_jes_up*pz1_type0_jes_up));
    //NU0_jes_dn.SetPxPyPzE(ReducedTree->METpuppiPtDown*TMath::Cos(ReducedTree->METpuppiPhiDown),ReducedTree->METpuppiPtDown*TMath::Sin(ReducedTree->METpuppiPhiDown),pz1_type0_jes_dn,TMath::Sqrt(WWTree->pfMETpuppi_jes_dn*WWTree->pfMETpuppi_jes_dn+pz1_type0_jes_dn*pz1_type0_jes_dn));
    
    NU2_puppi.SetPxPyPzE(ReducedTree->puppET*TMath::Cos(ReducedTree->puppETphi),ReducedTree->puppET*TMath::Sin(ReducedTree->puppETphi),WWTree->nu_pz_type2,TMath::Sqrt(WWTree->pfMETpuppi*WWTree->pfMETpuppi+WWTree->nu_pz_type2*WWTree->nu_pz_type2));
    NU1_puppi.SetPxPyPzE(ReducedTree->puppET*TMath::Cos(ReducedTree->puppETphi),ReducedTree->puppET*TMath::Sin(ReducedTree->puppETphi),WWTree->nu_pz_run2,TMath::Sqrt(WWTree->pfMETpuppi*WWTree->pfMETpuppi+WWTree->nu_pz_run2*WWTree->nu_pz_run2));
 //   	cout<<"\tNUmass = "<<NU1.M()<<"\t"<<NU1.Px()<<"\t"<<NU1.Py()<<"\t"<<NU1.Pz()<<"\t"<<NU1.E()<<endl;
	//cout<<"WWTree: "<<WWTree->nu_pz_run2<<endl;
    
    W_puppi = LEP + NU2_puppi;
    
    WWTree->v_puppi_pt = W_puppi.Pt();
    WWTree->v_puppi_eta = W_puppi.Eta();
    WWTree->v_puppi_phi = W_puppi.Phi();
    WWTree->v_puppi_mt = TMath::Sqrt(2*LEP.Et()*NU2_puppi.Et()*(1-TMath::Cos(LEP.DeltaPhi(NU2_puppi))));
    WWTree->v_puppi_mass = W_puppi.M();
    
    
    
    ///////////THE FAT JET - AK8
    float tempPt=0., tempMass=0.;
    int nGoodAK8jets=0;
    int ttb_jet_position=-1; //position of AK8 jet in ttbar-topology
    // if (ReducedTree->AK8CHS_ < 1 ) continue; 
    
    for ( int i=0; i<ReducedTree->AK8CHS_; i++)
    {
      bool isCleanedJet = true;
      if (ReducedTree->AK8CHS_pt[i]<100 || fabs(ReducedTree->AK8CHS_eta[i])>2.4)  continue; //be careful: this is not inside the synchntuple code
      if (ReducedTree->AddAK8CHS_mass_prun[i]>tempMass) {
        if ( (ReducedTree->AK8CHS_eta[i]>0 && WWTree->l_eta<0) || 
             (ReducedTree->AK8CHS_eta[i]<0 && WWTree->l_eta>0)) { //jet and lepton in opposite hemisphere for ttb
          ttb_jet_position=i; //save AK8 jet in ttbar topology
          tempMass=ReducedTree->AddAK8CHS_mass_prun[i];
        }
      }
      if (ReducedTree->AK8CHS_pt[i]<=tempPt) continue; //save the jet with the largest pt
      //if (ReducedTree->AK8Jets_AK8isLooseJetId[i]==false) continue; //fat jet must satisfy loose ID
      
      //CLEANING FROM LEPTONS
      for (std::size_t j=0; j<tightEle.size(); j++) {
        if (deltaR(tightEle.at(j).Eta(), tightEle.at(j).Phi(),
                   ReducedTree->AK8CHS_eta[i], ReducedTree->AK8CHS_phi[i]) < 1.0)
          isCleanedJet = false;
      }
      for ( std::size_t j=0; j<tightMuon.size(); j++) {
        if (deltaR(tightMuon.at(j).Eta(), tightMuon.at(j).Phi(),
                   ReducedTree->AK8CHS_eta[i], ReducedTree->AK8CHS_phi[i]) < 1.0)
          isCleanedJet = false;
      }
      
      if (isCleanedJet==false) continue; //jet is overlapped with a lepton
      
      WWTree->ungroomed_jet_pt  = ReducedTree->AK8CHS_pt[i];
      WWTree->ungroomed_jet_eta = ReducedTree->AK8CHS_eta[i];
      WWTree->ungroomed_jet_phi = ReducedTree->AK8CHS_phi[i];
      WWTree->ungroomed_jet_e   = ReducedTree->AK8CHS_ptRaw[i];
      //WWTree->ungroomed_jet_pt_jes_up = (ReducedTree->AK8CHS_pt[i]/ReducedTree->AK8Jets_AK8correction[i])*ReducedTree->AK8Jets_AK8correctionUp[i];
      //WWTree->ungroomed_jet_pt_jes_dn = (ReducedTree->AK8CHS_pt[i]/ReducedTree->AK8Jets_AK8correction[i])*ReducedTree->AK8Jets_AK8correctionDown[i];
      
      //WWTree->jet_pt_so    = ReducedTree->AK8Jets_softDropPt[i];
      WWTree->jet_mass_pr  = ReducedTree->AddAK8CHS_mass_prun[i];
      WWTree->jet_mass_so  = ReducedTree->AddAK8CHS_mass_sd0[i];
      WWTree->jet_mass_tr  = ReducedTree->AddAK8CHS_mass_trim[i];
      ////WWTree->jet_mass_fi  = ReducedTree->AK8Jets_filteredMass[i];
      WWTree->jet_tau2tau1 = ReducedTree->AddAK8CHS_tau2[i]/ReducedTree->AddAK8CHS_tau1[i];
      //WWTree->jet_mass_pr_jes_up = (ReducedTree->AddAK8CHS_mass_prun[i]/ReducedTree->AK8Jets_AK8massCorrection[i])*ReducedTree->AK8Jets_AK8massCorrectionUp[i];
      //WWTree->jet_mass_pr_jes_dn = (ReducedTree->AddAK8CHS_mass_prun[i]/ReducedTree->AK8Jets_AK8massCorrection[i])*ReducedTree->AK8Jets_AK8massCorrectionDown[i];
      
      tempPt = WWTree->ungroomed_jet_pt;
      nGoodAK8jets++;
    }
    if (WWTree->ungroomed_jet_pt > 0)
    {
      JET.SetPtEtaPhiE(WWTree->ungroomed_jet_pt,WWTree->ungroomed_jet_eta,WWTree->ungroomed_jet_phi,WWTree->ungroomed_jet_e);
      /* JET_jes_up.SetPtEtaPhiE(WWTree->ungroomed_jet_pt*(ReducedTree->AK8Jets_AK8correctionUp[hadWpos]/ReducedTree->AK8Jets_AK8correction[hadWpos]),
                              WWTree->ungroomed_jet_eta,
                              WWTree->ungroomed_jet_phi,
                              WWTree->ungroomed_jet_e*(ReducedTree->AK8Jets_AK8correctionUp[hadWpos]/ReducedTree->AK8Jets_AK8correction[hadWpos]));
      JET_jes_dn.SetPtEtaPhiE(WWTree->ungroomed_jet_pt*(ReducedTree->AK8Jets_AK8correctionDown[hadWpos]/ReducedTree->AK8Jets_AK8correction[hadWpos]),
                              WWTree->ungroomed_jet_eta,
                              WWTree->ungroomed_jet_phi,
                              WWTree->ungroomed_jet_e*(ReducedTree->AK8Jets_AK8correctionDown[hadWpos]/ReducedTree->AK8Jets_AK8correction[hadWpos]));
	*/
    }
    
    
    ///////////THE FAT JET - PuppiAK8
    tempPt=0., tempMass=0.;
    int nGoodPuppiAK8jets=0;
    int hadWPuppiAK8pos = -1;
    int ttb_PuppiAK8_jet_position=-1; //position of AK8 jet in ttbar-topology
    
    for ( int i=0; i<ReducedTree->AddAK8Puppi_; i++)
    {
      bool isCleanedJet = true;
      if (ReducedTree->AK8Puppi_pt[i]<100 || fabs(ReducedTree->AK8Puppi_eta[i])>2.4)  continue; //be careful: this is not inside the synchntuple code
      if (ReducedTree->AddAK8Puppi_mass_prun[i]>tempMass) {
        if ( (ReducedTree->AK8Puppi_eta[i]>0 && WWTree->l_eta<0) || 
             (ReducedTree->AK8Puppi_eta[i]<0 && WWTree->l_eta>0)) { //jet and lepton in opposite hemisphere for ttb
          ttb_PuppiAK8_jet_position=i; //save AK8 jet in ttbar topology
          tempMass=ReducedTree->AddAK8Puppi_mass_prun[i];
        }
      }
      if (ReducedTree->AK8Puppi_pt[i]<=tempPt) continue; //save the jet with the largest pt
      //if (ReducedTree->PuppiAK8Jets_PuppiAK8isLooseJetId[i]==false) continue; //fat jet must satisfy loose ID
      
      //CLEANING FROM LEPTONS
      for ( std::size_t j=0; j<tightEle.size(); j++) {
        if (deltaR(tightEle.at(j).Eta(), tightEle.at(j).Phi(),
                   ReducedTree->AK8Puppi_eta[i], ReducedTree->AK8Puppi_phi[i]) < 1.0)
          isCleanedJet = false;
      }
      for ( std::size_t j=0; j<tightMuon.size(); j++) {
        if (deltaR(tightMuon.at(j).Eta(), tightMuon.at(j).Phi(),
                   ReducedTree->AK8Puppi_eta[i], ReducedTree->AK8Puppi_phi[i]) < 1.0)
          isCleanedJet = false;
      }
      
      if (isCleanedJet==false) continue; //jet is overlapped with a lepton
      
      WWTree->ungroomed_PuppiAK8_jet_pt  = ReducedTree->AK8Puppi_pt[i];
      WWTree->ungroomed_PuppiAK8_jet_eta = ReducedTree->AK8Puppi_eta[i];
      WWTree->ungroomed_PuppiAK8_jet_phi = ReducedTree->AK8Puppi_phi[i];
      WWTree->ungroomed_PuppiAK8_jet_e   = ReducedTree->AK8Puppi_ptRaw[i];
      //WWTree->ungroomed_PuppiAK8_jet_pt_jes_up = (ReducedTree->AK8Puppi_pt[i]/ReducedTree->PuppiAK8Jets_PuppiAK8correction[i])*ReducedTree->PuppiAK8Jets_PuppiAK8correctionUp[i];
      //WWTree->ungroomed_PuppiAK8_jet_pt_jes_dn = (ReducedTree->AK8Puppi_pt[i]/ReducedTree->PuppiAK8Jets_PuppiAK8correction[i])*ReducedTree->PuppiAK8Jets_PuppiAK8correctionDown[i];
      
//      WWTree->PuppiAK8_jet_pt_so    = ReducedTree->PuppiAK8Jets_softDropPt[i];
      WWTree->PuppiAK8_jet_mass_pr  = ReducedTree->AddAK8Puppi_mass_prun[i];
//      WWTree->PuppiAK8_jet_mass_so  = (ReducedTree->AddAK8Puppi_mass_sd0[i]/ReducedTree->PuppiAK8Jets_PuppiAK8massCorrection[i])*getPUPPIweight(WWTree->ungroomed_PuppiAK8_jet_pt,WWTree->ungroomed_PuppiAK8_jet_eta, Puppisd_corrGEN, Puppisd_corrRECO_cen, Puppisd_corrRECO_for);
      WWTree->PuppiAK8_jet_mass_so  = ReducedTree->AddAK8Puppi_mass_sd0[i];
      WWTree->PuppiAK8_jet_mass_tr  = ReducedTree->AddAK8Puppi_mass_trim[i];
//     WWTree->PuppiAK8_jet_mass_fi  = ReducedTree->PuppiAK8Jets_filteredMass[i];
      WWTree->PuppiAK8_jet_tau2tau1 = ReducedTree->AddAK8Puppi_tau2[i]/ReducedTree->AddAK8Puppi_tau1[i];
 //     WWTree->PuppiAK8_jet_mass_pr_jes_up = (ReducedTree->AddAK8Puppi_mass_prun[i]/ReducedTree->PuppiAK8Jets_PuppiAK8massCorrection[i])*ReducedTree->PuppiAK8Jets_PuppiAK8massCorrectionUp[i];
   //   WWTree->PuppiAK8_jet_mass_pr_jes_dn = (ReducedTree->AddAK8Puppi_mass_prun[i]/ReducedTree->PuppiAK8Jets_PuppiAK8massCorrection[i])*ReducedTree->PuppiAK8Jets_PuppiAK8massCorrectionDown[i];
      
      tempPt = WWTree->ungroomed_PuppiAK8_jet_pt;
      nGoodPuppiAK8jets++;
      hadWPuppiAK8pos = i;
    }
    if (WWTree->ungroomed_PuppiAK8_jet_pt > 0.)
    {
      JET_PuppiAK8.SetPtEtaPhiE(WWTree->ungroomed_PuppiAK8_jet_pt,WWTree->ungroomed_PuppiAK8_jet_eta,WWTree->ungroomed_PuppiAK8_jet_phi,WWTree->ungroomed_PuppiAK8_jet_e);
      /* 
      JET_PuppiAK8_jes_up.SetPtEtaPhiE(WWTree->ungroomed_PuppiAK8_jet_pt*(ReducedTree->PuppiAK8Jets_PuppiAK8correctionUp[hadWPuppiAK8pos]/ReducedTree->PuppiAK8Jets_PuppiAK8correction[hadWPuppiAK8pos]),
                                       WWTree->ungroomed_PuppiAK8_jet_eta,
                                       WWTree->ungroomed_PuppiAK8_jet_phi,
                                       WWTree->ungroomed_PuppiAK8_jet_e*(ReducedTree->PuppiAK8Jets_PuppiAK8correctionUp[hadWPuppiAK8pos]/ReducedTree->PuppiAK8Jets_PuppiAK8correction[hadWPuppiAK8pos]));
      JET_PuppiAK8_jes_dn.SetPtEtaPhiE(WWTree->ungroomed_PuppiAK8_jet_pt*(ReducedTree->PuppiAK8Jets_PuppiAK8correctionDown[hadWPuppiAK8pos]/ReducedTree->PuppiAK8Jets_PuppiAK8correction[hadWPuppiAK8pos]),
                                       WWTree->ungroomed_PuppiAK8_jet_eta,
                                       WWTree->ungroomed_PuppiAK8_jet_phi,
                                       WWTree->ungroomed_PuppiAK8_jet_e*(ReducedTree->PuppiAK8Jets_PuppiAK8correctionDown[hadWPuppiAK8pos]/ReducedTree->PuppiAK8Jets_PuppiAK8correction[hadWPuppiAK8pos]));
	*/
    }
    
    
    // FAT JET SELECTION
    bool isGoodFatJet = true;
    if (nGoodAK8jets==0 && nGoodPuppiAK8jets==0) isGoodFatJet = false; //not found a good hadronic W candidate
    if (WWTree->ungroomed_jet_pt<150 && WWTree->ungroomed_PuppiAK8_jet_pt<150) isGoodFatJet = false;
    if (isGoodFatJet) cutEff[5]++;
    
    
    
    //////////////////UNMERGED HADRONIC W
    float tempPt1 = 0.; int pos1 = -1;
    float tempPt2 = 0.; int pos2 = -1;
    int nGoodAK4jets=0;
    for ( int i=0; i<ReducedTree->AK4CHS_; i++) //loop on AK4 jet
    {
      bool isCleanedJet = true;
      if (ReducedTree->AK4CHS_pt[i]<=30  || fabs(ReducedTree->AK4CHS_eta[i])>=2.4)  continue;
      //if (ReducedTree->Jets_isLooseJetId[i]==false) continue;
      if (ReducedTree->AK4CHS_csv[i]>0.890) continue;
      
      //CLEANING FROM LEPTONS
      for ( std::size_t j=0; j<tightEle.size(); j++) {
        if (deltaR(tightEle.at(j).Eta(), tightEle.at(j).Phi(),
                   ReducedTree->AK4CHS_eta[i], ReducedTree->AK4CHS_phi[i]) < 0.3) {
          isCleanedJet = false;
        }
      }
      for ( std::size_t j=0; j<tightMuon.size(); j++) {
        if (deltaR(tightMuon.at(j).Eta(), tightMuon.at(j).Phi(),
                   ReducedTree->AK4CHS_eta[i], ReducedTree->AK4CHS_phi[i]) < 0.3) {
          isCleanedJet = false;
        }
      }
      
      if (isCleanedJet==false) continue;
      
      if (ReducedTree->AK4CHS_pt[i]<tempPt1 && ReducedTree->AK4CHS_pt[i]<tempPt2) continue;
      
      if (ReducedTree->AK4CHS_pt[i]>tempPt1)
      {
        WWTree->AK4_jet1_pt  = ReducedTree->AK4CHS_pt[i];
        WWTree->AK4_jet1_eta = ReducedTree->AK4CHS_eta[i];
        WWTree->AK4_jet1_phi = ReducedTree->AK4CHS_phi[i];
        WWTree->AK4_jet1_e   = ReducedTree->AK4CHS_ptRaw[i];
        //WWTree->AK4_jet1_pt_jes_up = (ReducedTree->AK4CHS_pt[i]/ReducedTree->Jets_AK4correction[i])*ReducedTree->Jets_AK4correctionUp[i];
        //WWTree->AK4_jet1_pt_jes_dn = (ReducedTree->AK4CHS_pt[i]/ReducedTree->Jets_AK4correction[i])*ReducedTree->Jets_AK4correctionDown[i];
        
        tempPt1 = WWTree->AK4_jet1_pt;
        
        if (pos1!=-1)
        {
          WWTree->AK4_jet2_pt  = ReducedTree->AK4CHS_pt[pos1];
          WWTree->AK4_jet2_eta = ReducedTree->AK4CHS_eta[pos1];
          WWTree->AK4_jet2_phi = ReducedTree->AK4CHS_phi[pos1];
          WWTree->AK4_jet2_e   = ReducedTree->AK4CHS_ptRaw[pos1];
         // WWTree->AK4_jet2_pt_jes_up = (ReducedTree->AK4CHS_pt[pos1]/ReducedTree->Jets_AK4correction[pos1])*ReducedTree->Jets_AK4correctionUp[pos1];
         // WWTree->AK4_jet2_pt_jes_dn = (ReducedTree->AK4CHS_pt[pos1]/ReducedTree->Jets_AK4correction[pos1])*ReducedTree->Jets_AK4correctionDown[pos1];
          
          tempPt2 = WWTree->AK4_jet2_pt;
        }
        pos1 = i;
        nGoodAK4jets++;
      }
      else if (ReducedTree->AK4CHS_pt[i]>tempPt2)
      {
        WWTree->AK4_jet2_pt  = ReducedTree->AK4CHS_pt[i];
        WWTree->AK4_jet2_eta = ReducedTree->AK4CHS_eta[i];
        WWTree->AK4_jet2_phi = ReducedTree->AK4CHS_phi[i];
        WWTree->AK4_jet2_e   = ReducedTree->AK4CHS_ptRaw[i];
       // WWTree->AK4_jet2_pt_jes_up = (ReducedTree->AK4CHS_pt[i]/ReducedTree->Jets_AK4correction[i])*ReducedTree->Jets_AK4correctionUp[i];
       // WWTree->AK4_jet2_pt_jes_dn = (ReducedTree->AK4CHS_pt[i]/ReducedTree->Jets_AK4correction[i])*ReducedTree->Jets_AK4correctionDown[i];
        
        tempPt2 = WWTree->AK4_jet2_pt;
        pos2 = i;
        nGoodAK4jets++;
      }
    }
    
    if (WWTree->AK4_jet1_pt > 0.)
    {
      AK4_JET1.SetPtEtaPhiE(WWTree->AK4_jet1_pt,WWTree->AK4_jet1_eta,WWTree->AK4_jet1_phi,WWTree->AK4_jet1_e);
      //AK4_JET1_jes_up.SetPtEtaPhiE(WWTree->AK4_jet1_pt*(ReducedTree->Jets_AK4correctionUp[pos1]/ReducedTree->Jets_AK4correction[pos1]),
      //                             WWTree->AK4_jet1_eta,
      //                             WWTree->AK4_jet1_phi,
      //                             WWTree->AK4_jet1_e*(ReducedTree->Jets_AK4correctionUp[pos1]/ReducedTree->Jets_AK4correction[pos1]));
      //AK4_JET1_jes_dn.SetPtEtaPhiE(WWTree->AK4_jet1_pt*(ReducedTree->Jets_AK4correctionDown[pos1]/ReducedTree->Jets_AK4correction[pos1]),
      //                             WWTree->AK4_jet1_eta,
      //                             WWTree->AK4_jet1_phi,
      //                             WWTree->AK4_jet1_e*(ReducedTree->Jets_AK4correctionDown[pos1]/ReducedTree->Jets_AK4correction[pos1]));
    }
    if (WWTree->AK4_jet2_pt > 0.)
    {
      AK4_JET2.SetPtEtaPhiE(WWTree->AK4_jet2_pt,WWTree->AK4_jet2_eta,WWTree->AK4_jet2_phi,WWTree->AK4_jet2_e);
      //AK4_JET2_jes_up.SetPtEtaPhiE(WWTree->AK4_jet2_pt*(ReducedTree->Jets_AK4correctionUp[pos2]/ReducedTree->Jets_AK4correction[pos2]),
      //                             WWTree->AK4_jet2_eta,
      //                             WWTree->AK4_jet2_phi,
      //                             WWTree->AK4_jet2_e*(ReducedTree->Jets_AK4correctionUp[pos2]/ReducedTree->Jets_AK4correction[pos2]));
      //AK4_JET2_jes_dn.SetPtEtaPhiE(WWTree->AK4_jet2_pt*(ReducedTree->Jets_AK4correctionDown[pos2]/ReducedTree->Jets_AK4correction[pos2]),
      //                             WWTree->AK4_jet2_eta,
      //                             WWTree->AK4_jet2_phi,
      //                             WWTree->AK4_jet2_e*(ReducedTree->Jets_AK4correctionDown[pos2]/ReducedTree->Jets_AK4correction[pos2]));
    }
    
    if (WWTree->AK4_jet2_pt>0) {
      WWTree->AK4_jetjet_pt = (AK4_JET1+AK4_JET2).Pt();
      WWTree->AK4_jetjet_mass = (AK4_JET1+AK4_JET2).M();
      WWTree->AK4_jetjet_deltaeta = deltaEta(AK4_JET1.Eta(),AK4_JET2.Eta());
      WWTree->AK4_jetjet_deltaphi = deltaPhi(AK4_JET1.Phi(),AK4_JET2.Phi());
      WWTree->AK4_jetjet_deltar = deltaR(AK4_JET1.Eta(),AK4_JET1.Phi(),AK4_JET2.Eta(),AK4_JET2.Phi());
    }

    
    tempPt1 = 0.; int pos1Puppi = -1;
    tempPt2 = 0.; int pos2Puppi = -1;
    int nGoodPuppiAK4jets=0;
    for ( int i=0; i<ReducedTree->AK4Puppi_; i++) //loop on PuppiAK4 jet
    {
      bool isCleanedJet = true;
      if ( ReducedTree->AK4Puppi_pt[i]<=20 || fabs(ReducedTree->AK4Puppi_eta[i])>=2.4)  continue;
      //if (ReducedTree->JetsPuppi_isLooseJetId[i]==false) continue;
      if (ReducedTree->AK4Puppi_csv[i]>0.890) continue;
      
      //CLEANING FROM LEPTONS
      for ( std::size_t j=0; j<tightEle.size(); j++) {
        if (deltaR(tightEle.at(j).Eta(), tightEle.at(j).Phi(),
                   ReducedTree->AK4Puppi_eta[i], ReducedTree->AK4Puppi_phi[i]) < 0.3) {
          isCleanedJet = false;
        }
      }
      for ( std::size_t j=0; j<tightMuon.size(); j++) {
        if (deltaR(tightMuon.at(j).Eta(), tightMuon.at(j).Phi(),
                   ReducedTree->AK4Puppi_eta[i], ReducedTree->AK4Puppi_phi[i]) < 0.3) {
          isCleanedJet = false;
        }
      }
      
      if (isCleanedJet==false) continue;
      
      if (ReducedTree->AK4Puppi_pt[i]<tempPt1 && ReducedTree->AK4Puppi_pt[i]<tempPt2) continue;
      
      if (ReducedTree->AK4Puppi_pt[i]>tempPt1)
      {
        WWTree->PuppiAK4_jet1_pt  = ReducedTree->AK4Puppi_pt[i];
        WWTree->PuppiAK4_jet1_eta = ReducedTree->AK4Puppi_eta[i];
        WWTree->PuppiAK4_jet1_phi = ReducedTree->AK4Puppi_phi[i];
        WWTree->PuppiAK4_jet1_e   = ReducedTree->AK4Puppi_ptRaw[i];
      //  WWTree->PuppiAK4_jet1_pt_jes_up = (ReducedTree->AK4Puppi_pt[i]/ReducedTree->JetsPuppi_AK4correction[i])*ReducedTree->JetsPuppi_AK4correctionUp[i];
      //  WWTree->PuppiAK4_jet1_pt_jes_dn = (ReducedTree->AK4Puppi_pt[i]/ReducedTree->JetsPuppi_AK4correction[i])*ReducedTree->JetsPuppi_AK4correctionDown[i];
        
        tempPt1 = WWTree->PuppiAK4_jet1_pt;
        
        if (pos1Puppi!=-1)
        {
          WWTree->PuppiAK4_jet2_pt  = ReducedTree->AK4Puppi_pt[pos1Puppi];
          WWTree->PuppiAK4_jet2_eta = ReducedTree->AK4Puppi_eta[pos1Puppi];
          WWTree->PuppiAK4_jet2_phi = ReducedTree->AK4Puppi_phi[pos1Puppi];
          WWTree->PuppiAK4_jet2_e   = ReducedTree->AK4Puppi_ptRaw[pos1Puppi];
        //  WWTree->PuppiAK4_jet2_pt_jes_up = (ReducedTree->AK4Puppi_pt[pos1Puppi]/ReducedTree->JetsPuppi_AK4correction[pos1Puppi])*ReducedTree->JetsPuppi_AK4correctionUp[pos1Puppi];
        //  WWTree->PuppiAK4_jet2_pt_jes_dn = (ReducedTree->AK4Puppi_pt[pos1Puppi]/ReducedTree->JetsPuppi_AK4correction[pos1Puppi])*ReducedTree->JetsPuppi_AK4correctionDown[pos1Puppi];
          
          tempPt2 = WWTree->PuppiAK4_jet2_pt;
        }
        pos1Puppi = i;
        nGoodPuppiAK4jets++;
      }
      else if (ReducedTree->AK4Puppi_pt[i]>tempPt2)
      {
        WWTree->PuppiAK4_jet2_pt  = ReducedTree->AK4Puppi_pt[i];
        WWTree->PuppiAK4_jet2_eta = ReducedTree->AK4Puppi_eta[i];
        WWTree->PuppiAK4_jet2_phi = ReducedTree->AK4Puppi_phi[i];
        WWTree->PuppiAK4_jet2_e   = ReducedTree->AK4Puppi_ptRaw[i];
        //WWTree->PuppiAK4_jet2_pt_jes_up = (ReducedTree->AK4Puppi_pt[i]/ReducedTree->JetsPuppi_AK4correction[i])*ReducedTree->JetsPuppi_AK4correctionUp[i];
        //WWTree->PuppiAK4_jet2_pt_jes_dn = (ReducedTree->AK4Puppi_pt[i]/ReducedTree->JetsPuppi_AK4correction[i])*ReducedTree->JetsPuppi_AK4correctionDown[i];
        
        tempPt2 = WWTree->PuppiAK4_jet2_pt;
        pos2Puppi = i;
        nGoodPuppiAK4jets++;
      }
    }
    
    if (WWTree->PuppiAK4_jet1_pt > 0.)
    {
      PuppiAK4_JET1.SetPtEtaPhiE(WWTree->PuppiAK4_jet1_pt,WWTree->PuppiAK4_jet1_eta,WWTree->PuppiAK4_jet1_phi,WWTree->PuppiAK4_jet1_e);
      //PuppiAK4_JET1_jes_up.SetPtEtaPhiE(WWTree->PuppiAK4_jet1_pt*(ReducedTree->JetsPuppi_AK4correctionUp[pos1Puppi]/ReducedTree->JetsPuppi_AK4correction[pos1Puppi]),
      //                                  WWTree->PuppiAK4_jet1_eta,
      //                                  WWTree->PuppiAK4_jet1_phi,
      //                                  WWTree->PuppiAK4_jet1_e*(ReducedTree->JetsPuppi_AK4correctionUp[pos1Puppi]/ReducedTree->JetsPuppi_AK4correction[pos1Puppi]));
      //PuppiAK4_JET1_jes_dn.SetPtEtaPhiE(WWTree->PuppiAK4_jet1_pt*(ReducedTree->JetsPuppi_AK4correctionDown[pos1Puppi]/ReducedTree->JetsPuppi_AK4correction[pos1Puppi]),
      //                                  WWTree->PuppiAK4_jet1_eta,
      //                                  WWTree->PuppiAK4_jet1_phi,
      //                                  WWTree->PuppiAK4_jet1_e*(ReducedTree->JetsPuppi_AK4correctionDown[pos1Puppi]/ReducedTree->JetsPuppi_AK4correction[pos1Puppi]));
    }
    if (WWTree->PuppiAK4_jet2_pt > 0.)
    {
      PuppiAK4_JET2.SetPtEtaPhiE(WWTree->PuppiAK4_jet2_pt,WWTree->PuppiAK4_jet2_eta,WWTree->PuppiAK4_jet2_phi,WWTree->PuppiAK4_jet2_e);
      //PuppiAK4_JET2_jes_up.SetPtEtaPhiE(WWTree->PuppiAK4_jet2_pt*(ReducedTree->Jets_AK4correctionUp[pos2Puppi]/ReducedTree->Jets_AK4correction[pos2Puppi]),
      //                                  WWTree->PuppiAK4_jet2_eta,
      //                                  WWTree->PuppiAK4_jet2_phi,
      //                                  WWTree->PuppiAK4_jet2_e*(ReducedTree->JetsPuppi_AK4correctionUp[pos2Puppi]/ReducedTree->JetsPuppi_AK4correction[pos2Puppi]));
      //PuppiAK4_JET2_jes_dn.SetPtEtaPhiE(WWTree->PuppiAK4_jet2_pt*(ReducedTree->JetsPuppi_AK4correctionDown[pos2Puppi]/ReducedTree->JetsPuppi_AK4correction[pos2Puppi]),
      //                                  WWTree->PuppiAK4_jet2_eta,
      //                                  WWTree->PuppiAK4_jet2_phi,
      //                                  WWTree->PuppiAK4_jet2_e*(ReducedTree->JetsPuppi_AK4correctionDown[pos2Puppi]/ReducedTree->JetsPuppi_AK4correction[pos2Puppi]));
    }
    
    if (WWTree->PuppiAK4_jet2_pt>0) {
      WWTree->PuppiAK4_jetjet_pt = (PuppiAK4_JET1+PuppiAK4_JET2).Pt();
      WWTree->PuppiAK4_jetjet_mass = (PuppiAK4_JET1+PuppiAK4_JET2).M();
      WWTree->PuppiAK4_jetjet_deltaeta = deltaEta(PuppiAK4_JET1.Eta(),PuppiAK4_JET2.Eta());
      WWTree->PuppiAK4_jetjet_deltaphi = deltaPhi(PuppiAK4_JET1.Phi(),PuppiAK4_JET2.Phi());
      WWTree->PuppiAK4_jetjet_deltar = deltaR(PuppiAK4_JET1.Eta(),PuppiAK4_JET1.Phi(),PuppiAK4_JET2.Eta(),PuppiAK4_JET2.Phi());
    }
    
    
    // UNMERGED JETS SELECTION
    bool isGoodUnmergedJets = true;
    if (nGoodAK4jets<2 && nGoodPuppiAK4jets<2) isGoodUnmergedJets = false; //not found a good hadronic W candidate
    if (WWTree->AK4_jetjet_pt<150 && WWTree->PuppiAK4_jetjet_pt<150) isGoodUnmergedJets = false;
    if (isGoodUnmergedJets) cutEff[6]++;
    
    if (!isGoodFatJet && !isGoodUnmergedJets) continue;
    cutEff[7]++;
    
    
    //////////////////ANGULAR VARIABLES
    WWTree->deltaR_lak8jet = deltaR(JET.Eta(),JET.Phi(),LEP.Eta(),LEP.Phi());
    WWTree->deltaphi_METak8jet = deltaPhi(JET.Phi(),NU2.Phi());
    WWTree->deltaphi_Vak8jet = deltaPhi(JET.Phi(),W.Phi());
    WWTree->deltaR_lPuppiak8jet = deltaR(JET_PuppiAK8.Eta(),JET_PuppiAK8.Phi(),LEP.Eta(),LEP.Phi());
    WWTree->deltaphi_METPuppiak8jet = deltaPhi(JET_PuppiAK8.Phi(),NU2_puppi.Phi());
    WWTree->deltaphi_VPuppiak8jet = deltaPhi(JET_PuppiAK8.Phi(),W_puppi.Phi());
    if (WWTree->deltaR_lak8jet>(TMath::Pi()/2.0) && fabs(WWTree->deltaphi_METak8jet)>2.0 && fabs(WWTree->deltaphi_Vak8jet)>2.0 && nGoodAK8jets>0)
      WWTree->issignal=1;
    if (WWTree->deltaR_lPuppiak8jet>(TMath::Pi()/2.0) && fabs(WWTree->deltaphi_METPuppiak8jet)>2.0 && fabs(WWTree->deltaphi_VPuppiak8jet)>2.0 && nGoodPuppiAK8jets>0)
      WWTree->issignal_PuppiAK8=1;

    if (WWTree->AK4_jet2_pt>0) {
      WWTree->deltaR_lak4jetjet = deltaR((AK4_JET1+AK4_JET2).Eta(),(AK4_JET1+AK4_JET2).Phi(),LEP.Eta(),LEP.Phi());
      WWTree->deltaphi_METak4jetjet = deltaPhi((AK4_JET1+AK4_JET2).Phi(),NU2.Phi());
      WWTree->deltaphi_Vak4jetjet = deltaPhi((AK4_JET1+AK4_JET2).Phi(),W.Phi());
      if (WWTree->deltaR_lak4jetjet>(TMath::Pi()/2.0) && fabs(WWTree->deltaphi_METak4jetjet)>2.0 && fabs(WWTree->deltaphi_Vak4jetjet)>2.0)
        WWTree->issignal_AK4jetjet=1;
    }
    if (WWTree->PuppiAK4_jet2_pt>0) {
      WWTree->deltaR_lPuppiak4jetjet = deltaR((PuppiAK4_JET1+PuppiAK4_JET2).Eta(),(PuppiAK4_JET1+PuppiAK4_JET2).Phi(),LEP.Eta(),LEP.Phi());
      WWTree->deltaphi_METPuppiak4jetjet = deltaPhi((PuppiAK4_JET1+PuppiAK4_JET2).Phi(),NU2.Phi());
      WWTree->deltaphi_VPuppiak4jetjet = deltaPhi((PuppiAK4_JET1+PuppiAK4_JET2).Phi(),W.Phi());
      if (WWTree->deltaR_lPuppiak4jetjet>(TMath::Pi()/2.0) && fabs(WWTree->deltaphi_METPuppiak4jetjet)>2.0 && fabs(WWTree->deltaphi_VPuppiak4jetjet)>2.0)
        WWTree->issignal_PuppiAK4jetjet=1;
    }

    
    
    
    //////////////////FOUR-BODY INVARIANT MASS
    WWTree->mass_lvj_type0 = (LEP + NU0 + JET).M();
    WWTree->mass_lvj_type2 = (LEP + NU2 + JET).M();
    WWTree->mass_lvj_run2  = (LEP + NU1 + JET).M();
    if (isnan(WWTree->mass_lvj_run2) == 1)
    	cout<<"==============> Run2 mass is NAN"<<"\t LEP mass = "<<LEP.M()<<"\tNUmass = "<<NU1.M()<<"\t"<<NU1.Px()<<"\t"<<NU1.Py()<<"\t"<<NU1.Pz()<<"\t"<<NU1.E()<<endl;
    WWTree->mass_lvj_type0_met_jes_up = (LEP + NU0_jes_up + JET_jes_up).M();
    WWTree->mass_lvj_type0_met_jes_dn = (LEP + NU0_jes_dn + JET_jes_dn).M();
    
    //cout<<"\n=================="<<endl;
    //cout<<"Hadronic W:	Eta = "<<( JET_PuppiAK8).Eta()<<"\tRapidity = "<<( JET_PuppiAK8).Rapidity()<<endl;
    //cout<<"Leptonic W:	Eta = "<<(LEP + NU0_puppi ).Eta()<<"\tRapidity = "<<(LEP + NU0_puppi).Rapidity()<<endl;
    //cout<<"WW System:	Eta = "<<(LEP + NU0_puppi + JET_PuppiAK8).Eta()<<"\tRapidity = "<<(LEP + NU0_puppi + JET_PuppiAK8).Rapidity()<<endl;
    WWTree->TempLepWEta = (LEP + NU0_puppi ).Eta();
    WWTree->TempLepWRapidity = (LEP + NU0_puppi ).Rapidity();
    WWTree->TempHadWEta = (JET_PuppiAK8 ).Eta();
    WWTree->TempHadWRapidity = (JET_PuppiAK8 ).Rapidity();
    WWTree->TempWWEta = (LEP + NU0_puppi + JET_PuppiAK8 ).Eta();
    WWTree->TempWWRapidity = (LEP + NU0_puppi + JET_PuppiAK8 ).Rapidity();
    WWTree->pt_lvj_type0_PuppiAK8 = (LEP + NU0_puppi + JET_PuppiAK8).Pt();
    WWTree->pt_lvj_type2_PuppiAK8 = (LEP + NU2_puppi + JET_PuppiAK8).Pt();
    WWTree->pt_lvj_run2_PuppiAK8  = (LEP + NU1_puppi + JET_PuppiAK8).Pt();
    WWTree->eta_lvj_type0_PuppiAK8 = (LEP + NU0_puppi + JET_PuppiAK8).Eta();
    WWTree->eta_lvj_type2_PuppiAK8 = (LEP + NU2_puppi + JET_PuppiAK8).Eta();
    WWTree->eta_lvj_run2_PuppiAK8  = (LEP + NU1_puppi + JET_PuppiAK8).Eta();
    WWTree->rapidity_lvj_type0_PuppiAK8 = (LEP + NU0_puppi + JET_PuppiAK8).Rapidity();
    WWTree->rapidity_lvj_type2_PuppiAK8 = (LEP + NU2_puppi + JET_PuppiAK8).Rapidity();
    WWTree->rapidity_lvj_run2_PuppiAK8  = (LEP + NU1_puppi + JET_PuppiAK8).Rapidity();
    WWTree->phi_lvj_type0_PuppiAK8 = (LEP + NU0_puppi + JET_PuppiAK8).Phi();
    WWTree->phi_lvj_type2_PuppiAK8 = (LEP + NU2_puppi + JET_PuppiAK8).Phi();
    WWTree->phi_lvj_run2_PuppiAK8  = (LEP + NU1_puppi + JET_PuppiAK8).Phi();
    WWTree->energy_lvj_type0_PuppiAK8 = (LEP + NU0_puppi + JET_PuppiAK8).E();
    WWTree->energy_lvj_type2_PuppiAK8 = (LEP + NU2_puppi + JET_PuppiAK8).E();
    WWTree->energy_lvj_run2_PuppiAK8  = (LEP + NU1_puppi + JET_PuppiAK8).E();
    WWTree->mass_lvj_type0_PuppiAK8 = (LEP + NU0_puppi + JET_PuppiAK8).M();
    WWTree->mass_lvj_type2_PuppiAK8 = (LEP + NU2_puppi + JET_PuppiAK8).M();
    WWTree->mass_lvj_run2_PuppiAK8  = (LEP + NU1_puppi + JET_PuppiAK8).M();
    WWTree->mt_lvj_type0_PuppiAK8 = (LEP + NU0_puppi + JET_PuppiAK8).Mt();
    WWTree->mt_lvj_type2_PuppiAK8 = (LEP + NU2_puppi + JET_PuppiAK8).Mt();
    WWTree->mt_lvj_run2_PuppiAK8  = (LEP + NU1_puppi + JET_PuppiAK8).Mt();
    WWTree->mass_lvj_type0_met_PuppiAK8_jes_up = (LEP + NU0_jes_up + JET_PuppiAK8_jes_up).M();
    WWTree->mass_lvj_type0_met_PuppiAK8_jes_dn = (LEP + NU0_jes_dn + JET_PuppiAK8_jes_dn).M();
    
    WWTree->mass_lvjj_type0_AK4 = (LEP + NU0 + AK4_JET1 + AK4_JET2).M();
    WWTree->mass_lvjj_type2_AK4 = (LEP + NU2 + AK4_JET1 + AK4_JET2).M();
    WWTree->mass_lvjj_run2_AK4  = (LEP + NU1 + AK4_JET1 + AK4_JET2).M();
    WWTree->mass_lvjj_type0_met_jes_up_AK4 = (LEP + NU0_jes_up + AK4_JET1_jes_up + AK4_JET2_jes_up).M();
    WWTree->mass_lvjj_type0_met_jes_dn_AK4 = (LEP + NU0_jes_dn + AK4_JET1_jes_dn + AK4_JET2_jes_dn).M();
    
    WWTree->mass_lvjj_type0_PuppiAK4 = (LEP + NU0 + PuppiAK4_JET1 + PuppiAK4_JET2).M();
    WWTree->mass_lvjj_type2_PuppiAK4 = (LEP + NU2 + PuppiAK4_JET1 + PuppiAK4_JET2).M();
    WWTree->mass_lvjj_run2_PuppiAK4  = (LEP + NU1 + PuppiAK4_JET1 + PuppiAK4_JET2).M();
    WWTree->mass_lvjj_type0_met_jes_up_PuppiAK4 = (LEP + NU0_jes_up + PuppiAK4_JET1_jes_up + PuppiAK4_JET2_jes_up).M();
    WWTree->mass_lvjj_type0_met_jes_dn_PuppiAK4 = (LEP + NU0_jes_dn + PuppiAK4_JET1_jes_dn + PuppiAK4_JET2_jes_dn).M();
    

    
    //--- ttbar topology ------
    if (ttb_PuppiAK8_jet_position>=0)
    {
      WWTree->ttb_ungroomed_jet_pt  = ReducedTree->AK8CHS_pt[ttb_PuppiAK8_jet_position];
      WWTree->ttb_ungroomed_jet_eta = ReducedTree->AK8CHS_eta[ttb_PuppiAK8_jet_position];
      WWTree->ttb_ungroomed_jet_phi = ReducedTree->AK8CHS_phi[ttb_PuppiAK8_jet_position];
      WWTree->ttb_ungroomed_jet_e   = ReducedTree->AK8CHS_ptRaw[ttb_PuppiAK8_jet_position];
      WWTree->ttb_jet_mass_pr   = ReducedTree->AddAK8CHS_mass_prun[ttb_PuppiAK8_jet_position];
      WWTree->ttb_jet_mass_so   = ReducedTree->AddAK8CHS_mass_sd0[ttb_PuppiAK8_jet_position];
//      WWTree->ttb_jet_pt_so   = ReducedTree->AK8Jets_softDropPt[ttb_PuppiAK8_jet_position];
      WWTree->ttb_jet_mass_tr   = ReducedTree->AddAK8CHS_mass_trim[ttb_PuppiAK8_jet_position];
      //WWTree->ttb_jet_mass_fi   = ReducedTree->AK8Jets_filteredMass[ttb_PuppiAK8_jet_position];
      WWTree->ttb_jet_tau2tau1   = ReducedTree->AddAK8CHS_tau2[ttb_PuppiAK8_jet_position]/ReducedTree->AddAK8CHS_tau1[ttb_PuppiAK8_jet_position];
      
      WWTree->ttb_deltaeta_lak8jet = deltaEta(WWTree->ttb_ungroomed_jet_eta,WWTree->l_eta);
    }
    
    
    /////////VBF and b-tag section
    WWTree->njets=0;
    WWTree->nBTagJet_loose=0;
    WWTree->nBTagJet_medium=0;
    WWTree->nBTagJet_tight=0;
    
    WWTree->njets_unmerged=0;
    WWTree->nBTagJet_loose_unmerged=0;
    WWTree->nBTagJet_medium_unmerged=0;
    WWTree->nBTagJet_tight_unmerged=0;

    float oldDeltaR = 1000.;
    float oldDeltaRLep = 1000.;
    int indexCloserJet = -1;
    int indexCloserJetLep = -1;
    float deltaRbtag_prev=100.;
    float deltaRbtag_prev_loose=100.;
    
    std::vector<int> indexGoodVBFJets;
    
    for ( int i=0; i<ReducedTree->AK4CHS_; i++) //loop on AK4 jet
    {
      bool isCleaned = true;
      bool isCleanedFromFatJet = true;
      bool isCleanedFromUnmergedJets = true;
      
      if (ReducedTree->AK4CHS_pt[i]<=30 ) continue;
      //if (ReducedTree->Jets_isLooseJetId[i]==false) continue;
      
      //CLEANING FROM FAT JET
      if (nGoodAK8jets > 0) {
        if (deltaR(WWTree->ungroomed_jet_eta, WWTree->ungroomed_jet_phi,
                   ReducedTree->AK4CHS_eta[i],ReducedTree->AK4CHS_phi[i]) < 0.8 )
          isCleanedFromFatJet = false;
      }
      
      //CLEANING FROM UNMERGED JETS
      if (nGoodAK4jets>0) {
        if (deltaR(WWTree->AK4_jet1_eta, WWTree->AK4_jet1_phi,
                   ReducedTree->AK4CHS_eta[i],ReducedTree->AK4CHS_phi[i]) < 0.4 )
          isCleanedFromUnmergedJets = false;
      }
      if (nGoodAK4jets>1) {
        if (deltaR(WWTree->AK4_jet2_eta, WWTree->AK4_jet2_phi,
                   ReducedTree->AK4CHS_eta[i],ReducedTree->AK4CHS_phi[i]) < 0.4 )
          isCleanedFromUnmergedJets = false;
      }      
      
      //CLEANING FROM LEPTONS
      for ( std::size_t j=0; j<tightEle.size(); j++) {
        if (deltaR(tightEle.at(j).Eta(), tightEle.at(j).Phi(),
                   ReducedTree->AK4CHS_eta[i],   ReducedTree->AK4CHS_phi[i]) < 0.3) {
          isCleaned = false;
        }
      }
      for ( std::size_t j=0; j<tightMuon.size(); j++) {
        if (deltaR(tightMuon.at(j).Eta(), tightMuon.at(j).Phi(),
                   ReducedTree->AK4CHS_eta[i],   ReducedTree->AK4CHS_phi[i]) < 0.3) {
          isCleaned = false;
        }
      }
      
      if (isCleaned==false) continue;
      
      
      if (isCleanedFromUnmergedJets==true && fabs(ReducedTree->AK4CHS_eta[i])<2.4)
      {
        WWTree->njets_unmerged++;
        if (ReducedTree->AK4CHS_csv[i]>0.605) WWTree->nBTagJet_loose_unmerged++;
        if (ReducedTree->AK4CHS_csv[i]>0.890) WWTree->nBTagJet_medium_unmerged++;
        if (ReducedTree->AK4CHS_csv[i]>0.970) WWTree->nBTagJet_tight_unmerged++;
      }
      
      if (isCleanedFromFatJet==false) continue;
      
      indexGoodVBFJets.push_back(i); //save index of the "good" vbf jets candidates
      
      if (fabs(ReducedTree->AK4CHS_eta[i])>=2.4) continue;
      
      WWTree->njets++;
      AK4.SetPtEtaPhiE(ReducedTree->AK4CHS_pt[i],ReducedTree->AK4CHS_eta[i],ReducedTree->AK4CHS_phi[i],ReducedTree->AK4CHS_ptRaw[i]);
      
      
      //fill B-Tag info
      if (ReducedTree->AK4CHS_csv[i]>0.605) {
        WWTree->nBTagJet_loose++;
      }
      
      if (ReducedTree->AK4CHS_csv[i]>0.890) {  
        WWTree->nBTagJet_medium++;
      }
      
      if (ReducedTree->AK4CHS_csv[i]>0.970) {
        WWTree->nBTagJet_tight++;
      }
      
      
      //------------------------------
      // !!! VBF non-Puppi missing !!!
      //------------------------------
    }
    
    
    
    /////////VBF and b-tag section Puppi
    WWTree->njetsPuppi=0;
    WWTree->nBTagJetPuppi_loose=0;
    WWTree->nBTagJetPuppi_medium=0;
    WWTree->nBTagJetPuppi_tight=0;
    
    WWTree->njetsPuppi_unmerged=0;
    WWTree->nBTagJetPuppi_loose_unmerged=0;
    WWTree->nBTagJetPuppi_medium_unmerged=0;
    WWTree->nBTagJetPuppi_tight_unmerged=0;

    oldDeltaR = 1000.;
    oldDeltaRLep = 1000.;
    indexCloserJet = -1;
    indexCloserJetLep = -1;
    deltaRbtag_prev=100.;
    deltaRbtag_prev_loose=100.;
    
    std::vector<int> indexGoodVBFJetsPuppi;
    
    for ( int i=0; i<ReducedTree->AK4Puppi_; i++) //loop on PuppiAK4 jet
    {
      bool isCleaned = true;
      bool isCleanedFromFatJet = true;
      bool isCleanedFromUnmergedJets = true;
      
      if (ReducedTree->AK4Puppi_pt[i]<=30 || ReducedTree->AK4Puppi_pt[i]<=20) continue;
      //if (ReducedTree->JetsPuppi_isLooseJetId[i]==false) continue;
      
      //CLEANING FROM FAT JET
      if (nGoodPuppiAK8jets > 0) {
        if (deltaR(WWTree->ungroomed_PuppiAK8_jet_eta, WWTree->ungroomed_PuppiAK8_jet_phi,
                   ReducedTree->AK4Puppi_eta[i],ReducedTree->AK4Puppi_phi[i]) < 0.8 )
          isCleanedFromFatJet = false;
      }
      
      //CLEANING FROM UNMERGED JETS
      if (nGoodPuppiAK4jets>0) {
        if (deltaR(WWTree->PuppiAK4_jet1_eta, WWTree->PuppiAK4_jet1_phi,
                   ReducedTree->AK4Puppi_eta[i],ReducedTree->AK4Puppi_phi[i]) < 0.4 )
          isCleanedFromUnmergedJets = false;
      }
      if (nGoodPuppiAK4jets>1) {
        if (deltaR(WWTree->PuppiAK4_jet2_eta, WWTree->PuppiAK4_jet2_phi,
                   ReducedTree->AK4Puppi_eta[i],ReducedTree->AK4Puppi_phi[i]) < 0.4 )
          isCleanedFromUnmergedJets = false;
      }      
      
      //CLEANING FROM LEPTONS
      for ( std::size_t j=0; j<tightEle.size(); j++) {
        if (deltaR(tightEle.at(j).Eta(), tightEle.at(j).Phi(),
                   ReducedTree->AK4Puppi_eta[i],   ReducedTree->AK4Puppi_phi[i]) < 0.3) {
          isCleaned = false;
        }
      }
      for ( std::size_t j=0; j<tightMuon.size(); j++) {
        if (deltaR(tightMuon.at(j).Eta(), tightMuon.at(j).Phi(),
                   ReducedTree->AK4Puppi_eta[i],   ReducedTree->AK4Puppi_phi[i]) < 0.3) {
          isCleaned = false;
        }
      }
      
      if (isCleaned==false) continue;
      
      
      if (isCleanedFromUnmergedJets==true && fabs(ReducedTree->AK4Puppi_eta[i])<2.4)
      {
        WWTree->njetsPuppi_unmerged++;
        if (ReducedTree->AK4Puppi_csv[i]>0.605) WWTree->nBTagJetPuppi_loose_unmerged++;
        if (ReducedTree->AK4Puppi_csv[i]>0.890) WWTree->nBTagJetPuppi_medium_unmerged++;
        if (ReducedTree->AK4Puppi_csv[i]>0.970) WWTree->nBTagJetPuppi_tight_unmerged++;
      }
      
      if (isCleanedFromFatJet==false) continue;
      
      indexGoodVBFJetsPuppi.push_back(i); //save index of the "good" vbf jets candidates
      
      if (fabs(ReducedTree->AK4Puppi_eta[i])>=2.4) continue;
      
      WWTree->njetsPuppi++;
      AK4.SetPtEtaPhiE(ReducedTree->AK4Puppi_pt[i],ReducedTree->AK4Puppi_eta[i],ReducedTree->AK4Puppi_phi[i],ReducedTree->AK4Puppi_ptRaw[i]);
      
      
      //fill B-Tag info
      if (ReducedTree->AK4Puppi_csv[i]>0.605) { 
        WWTree->nBTagJetPuppi_loose++;
        float deltaRbtag = JET_PuppiAK8.DeltaR(AK4);
        if (deltaRbtag>0.8 && deltaRbtag<deltaRbtag_prev_loose) {
          WWTree->deltaR_AK8_closestBtagJet_loose = deltaRbtag;
          deltaRbtag_prev_loose = deltaRbtag;
        }	  
      }
      
      if (ReducedTree->AK4Puppi_csv[i]>0.890) {  
        WWTree->nBTagJetPuppi_medium++;
        float deltaRbtag = JET_PuppiAK8.DeltaR(AK4);
        if (deltaRbtag>0.8 && deltaRbtag<deltaRbtag_prev) {
          WWTree->deltaR_AK8_closestBtagJet = deltaRbtag;
          deltaRbtag_prev = deltaRbtag;
        }	  
      }
      
      if (ReducedTree->AK4Puppi_csv[i]>0.970) {
        WWTree->nBTagJetPuppi_tight++;
      }
      
      float deltaRlep = W.DeltaR(AK4);
      if (deltaRlep<oldDeltaRLep) indexCloserJetLep = i;
      
      float deltaR = JET_PuppiAK8.DeltaR(AK4);
      if (deltaR<0.8) continue; //the vbf jets must be outside the had W cone
      
      if (WWTree->njets!=0) {
        if (WWTree->jet2_pt!=0) {
          WWTree->jet3_pt=ReducedTree->AK4Puppi_pt[i];
          WWTree->jet3_eta=ReducedTree->AK4Puppi_eta[i];
          WWTree->jet3_phi=ReducedTree->AK4Puppi_phi[i];
          WWTree->jet3_e=ReducedTree->AK4Puppi_ptRaw[i];
          WWTree->jet3_btag=ReducedTree->AK4Puppi_csv[i];
        }
        else {
          WWTree->jet2_pt=ReducedTree->AK4Puppi_pt[i];
          WWTree->jet2_eta=ReducedTree->AK4Puppi_eta[i];
          WWTree->jet2_phi=ReducedTree->AK4Puppi_phi[i];
          WWTree->jet2_e=ReducedTree->AK4Puppi_ptRaw[i];
          WWTree->jet2_btag=ReducedTree->AK4Puppi_csv[i];
        }
      }	
      
      if (deltaR<oldDeltaR)  indexCloserJet = i; //index of the closest jet to the AK8
    }
    
    
    if (indexCloserJet>=0) { //fill hadronic top mass
      AK4.SetPtEtaPhiE(ReducedTree->AK4Puppi_pt[indexCloserJet],ReducedTree->AK4Puppi_eta[indexCloserJet],ReducedTree->AK4Puppi_phi[indexCloserJet],ReducedTree->AK4Puppi_ptRaw[indexCloserJet]);
      WWTree->mass_ungroomedjet_closerjet  = (JET_PuppiAK8 + AK4).M();
      WWTree->AK8_closerjet_pt = AK4.Pt();
      WWTree->AK8_closerjet_eta = AK4.Eta();
      WWTree->AK8_closerjet_phi = AK4.Phi();
      WWTree->AK8_closerjet_e = AK4.E();
    }
    if (indexCloserJetLep>=0) { //fill leptonic top mass
      AK4.SetPtEtaPhiE(ReducedTree->AK4Puppi_pt[indexCloserJetLep],ReducedTree->AK4Puppi_eta[indexCloserJetLep],ReducedTree->AK4Puppi_phi[indexCloserJetLep],ReducedTree->AK4Puppi_ptRaw[indexCloserJetLep]);
      WWTree->mass_leptonic_closerjet  = (W + AK4).M();
    }
    
    
    if (indexGoodVBFJetsPuppi.size()>=2) 
    {
      float tempPtMax=0.;
      int nVBF1=-1, nVBF2=-1; //position of the two vbf jets
      
      for (std::size_t i=0; i<indexGoodVBFJetsPuppi.size()-1; i++) {
        for ( std::size_t ii=i+1; ii<indexGoodVBFJetsPuppi.size(); ii++) {
          VBF1.SetPtEtaPhiE(ReducedTree->AK4Puppi_pt[indexGoodVBFJetsPuppi.at(i)],ReducedTree->AK4Puppi_eta[indexGoodVBFJetsPuppi.at(i)],ReducedTree->AK4Puppi_phi[indexGoodVBFJetsPuppi.at(i)],ReducedTree->AK4Puppi_ptRaw[indexGoodVBFJetsPuppi.at(i)]);
          VBF2.SetPtEtaPhiE(ReducedTree->AK4Puppi_pt[indexGoodVBFJetsPuppi.at(ii)],ReducedTree->AK4Puppi_eta[indexGoodVBFJetsPuppi.at(ii)],ReducedTree->AK4Puppi_phi[indexGoodVBFJetsPuppi.at(ii)],ReducedTree->AK4Puppi_ptRaw[indexGoodVBFJetsPuppi.at(ii)]);
          TOT = VBF1 + VBF2;
          if (TOT.Pt() < tempPtMax) continue;
          tempPtMax = TOT.Pt(); //take the jet pair with largest Pt
          nVBF1 = indexGoodVBFJetsPuppi.at(i); //save position of the 1st vbf jet
          nVBF2 = indexGoodVBFJetsPuppi.at(ii); //save position of the 2nd vbf jet
        }
      }
      
      if (nVBF1!=-1 && nVBF2!=-1) //save infos for vbf jet pair
      {
        nVBF1 = indexGoodVBFJetsPuppi.at(0); //save position of the 1st vbf jet
        nVBF2 = indexGoodVBFJetsPuppi.at(1); //save position of the 2nd vbf jet
        // nVBF1=0; nVBF2=1;
        
        VBF1.SetPtEtaPhiE(ReducedTree->AK4CHS_pt[nVBF1],ReducedTree->AK4CHS_eta[nVBF1],ReducedTree->AK4CHS_phi[nVBF1],ReducedTree->AK4CHS_ptRaw[nVBF1]);
        VBF2.SetPtEtaPhiE(ReducedTree->AK4CHS_pt[nVBF2],ReducedTree->AK4CHS_eta[nVBF2],ReducedTree->AK4CHS_phi[nVBF2],ReducedTree->AK4CHS_ptRaw[nVBF2]);
        TOT = VBF1 + VBF2;
	
        WWTree->vbf_maxpt_j1_pt = ReducedTree->AK4CHS_pt[nVBF1];
        WWTree->vbf_maxpt_j1_eta = ReducedTree->AK4CHS_eta[nVBF1];
        WWTree->vbf_maxpt_j1_phi = ReducedTree->AK4CHS_phi[nVBF1];
        WWTree->vbf_maxpt_j1_e = ReducedTree->AK4CHS_ptRaw[nVBF1];
        WWTree->vbf_maxpt_j1_bDiscriminatorCSV = ReducedTree->AK4CHS_csv[nVBF1];
        WWTree->vbf_maxpt_j2_pt = ReducedTree->AK4CHS_pt[nVBF2];
        WWTree->vbf_maxpt_j2_eta = ReducedTree->AK4CHS_eta[nVBF2];
        WWTree->vbf_maxpt_j2_phi = ReducedTree->AK4CHS_phi[nVBF2];
        WWTree->vbf_maxpt_j2_e = ReducedTree->AK4CHS_ptRaw[nVBF2];
        WWTree->vbf_maxpt_j2_bDiscriminatorCSV = ReducedTree->AK4CHS_csv[nVBF2];
        WWTree->vbf_maxpt_jj_pt = TOT.Pt();
        WWTree->vbf_maxpt_jj_eta = TOT.Eta();
        WWTree->vbf_maxpt_jj_phi = TOT.Phi();
        WWTree->vbf_maxpt_jj_m = TOT.M();	
      }
    }
    
    
    /////////////////MC Infos
    if (isMC==1)
    {
    	TLorentzVector hadW, lepW, VBFJ1, VBFJ2, VBFJ, temp;

	for (int i = 0; i<ReducedTree->GenParticle_;i++)
	{
	   if( (abs(ReducedTree->GenParticle_pdgId[i]) == 11 || abs(ReducedTree->GenParticle_pdgId[i]) == 13 || abs(ReducedTree->GenParticle_pdgId[i]) == 15) && ReducedTree->GenParticle_status[i] == 23)
	   if( (abs(ReducedTree->GenParticle_parent[i]) == 24))
	   	cout<<"Lep W = "<<ReducedTree->GenParticle_parent[i]<<endl;


	   if( (abs(ReducedTree->GenParticle_pdgId[i]) == 1 || abs(ReducedTree->GenParticle_pdgId[i]) == 2 || abs(ReducedTree->GenParticle_pdgId[i]) == 3 || abs(ReducedTree->GenParticle_pdgId[i]) == 4 || abs(ReducedTree->GenParticle_pdgId[i]) == 5 || abs(ReducedTree->GenParticle_pdgId[i]) == 6) && ReducedTree->GenParticle_status[i] == 23 )
	   if( (abs(ReducedTree->GenParticle_parent[i]) == 24))
	   	cout<<"Had W = "<<ReducedTree->GenParticle_parent[i]<<endl;
	}
    }
    

    WWTree->totalEventWeight = WWTree->genWeight*WWTree->eff_and_pu_Weight*WWTree->top1_NNLO_Weight*WWTree->top2_NNLO_Weight*WWTree->trig_eff_Weight;
    WWTree->totalEventWeight_2 = WWTree->genWeight*WWTree->eff_and_pu_Weight_2*WWTree->top1_NNLO_Weight*WWTree->top2_NNLO_Weight*WWTree->trig_eff_Weight;
    WWTree->totalEventWeight_3 = WWTree->genWeight*WWTree->eff_and_pu_Weight_3*WWTree->top1_NNLO_Weight*WWTree->top2_NNLO_Weight*WWTree->trig_eff_Weight;

    
    bool isBadEvent=false;
    if (isMC==0) {
      std::multimap<int,int>::iterator it = badEventsList.begin(); //filter bad events
      for (it=badEventsList.begin(); it!=badEventsList.end(); ++it) {
	if (it->first == WWTree->run && it->second == WWTree->event)
	  isBadEvent = true;
      }
    }
    if (isBadEvent) continue;
    
    //    WWTree->wSampleWeight = std::atof(xSecWeight.c_str())/nEvents; //xsec/numberOfEntries
    WWTree->nEvents = nEvents;
    WWTree->nNegEvents = nNegEvents;

    double a_costheta1, a_costheta2, a_costhetastar, a_Phi, a_Phi1;
    computeAngles( LEP + NU0_puppi + JET_PuppiAK8, LEP + NU0_puppi, LEP, NU0_puppi, JET_PuppiAK8,  a_costheta1, a_costhetastar, a_Phi1);
    WWTree->costheta1_type0 = (float) a_costheta1;                
    WWTree->costhetastar_type0 = (float) a_costhetastar;
    WWTree->phi1_type0 = (float) a_Phi1;

    computeAngles( LEP + NU2_puppi + JET_PuppiAK8, LEP + NU2_puppi, LEP, NU2_puppi, JET_PuppiAK8,  a_costheta1, a_costhetastar, a_Phi1);
    WWTree->costheta1_type2 = (float) a_costheta1;                
    WWTree->costhetastar_type2 = (float) a_costhetastar;
    WWTree->phi1_type2 = (float) a_Phi1;

    computeAngles( LEP + NU1_puppi + JET_PuppiAK8, LEP + NU1_puppi, LEP, NU1_puppi, JET_PuppiAK8,  a_costheta1, a_costhetastar, a_Phi1);
    WWTree->costheta1_run2 = (float) a_costheta1;                
    WWTree->costhetastar_run2 = (float) a_costhetastar;
    WWTree->phi1_run2 = (float) a_Phi1;

    if (fabs(VBF1.Eta() - VBF2.Eta()) == 0.0)
    {
    	 WWTree->VBSCentrality_type0 = -999.0;
    	 WWTree->VBSCentrality_type2 = -999.0;
    	 WWTree->VBSCentrality_run2 = -999.0;
    }
    else
    {
    	WWTree->VBSCentrality_type0 = (fabs(VBF1.Eta()- (((LEP + NU0_puppi).Rapidity()+JET_PuppiAK8.Rapidity()))- VBF2.Eta() ))/fabs(VBF1.Eta() - VBF2.Eta());
    	WWTree->VBSCentrality_type2 = (fabs(VBF1.Eta()- (((LEP + NU2_puppi).Rapidity()+JET_PuppiAK8.Rapidity()))- VBF2.Eta() ))/fabs(VBF1.Eta() - VBF2.Eta());
    	WWTree->VBSCentrality_run2 = (fabs(VBF1.Eta()- (((LEP + NU1_puppi).Rapidity()+JET_PuppiAK8.Rapidity()))- VBF2.Eta() ))/fabs(VBF1.Eta() - VBF2.Eta());
    }
     WWTree->RpT_type0 = (JET_PuppiAK8.Pt()*(LEP + NU0_puppi).Pt())/(VBF1.Pt()*VBF2.Pt());
     WWTree->RpT_type2 = (JET_PuppiAK8.Pt()*(LEP + NU2_puppi).Pt())/(VBF1.Pt()*VBF2.Pt());
     WWTree->RpT_run2 = (JET_PuppiAK8.Pt()*(LEP + NU1_puppi).Pt())/(VBF1.Pt()*VBF2.Pt());
     WWTree->ZeppenfeldWH = JET_PuppiAK8.Rapidity() - (VBF1.Eta() + VBF2.Eta())/2.0;
     WWTree->ZeppenfeldWL_type0 = (LEP + NU0_puppi).Rapidity() - (VBF1.Eta() + VBF2.Eta())/2.0;
     WWTree->ZeppenfeldWL_type2 = (LEP + NU2_puppi).Rapidity() - (VBF1.Eta() + VBF2.Eta())/2.0;
     WWTree->ZeppenfeldWL_run2 = (LEP + NU1_puppi).Rapidity() - (VBF1.Eta() + VBF2.Eta())/2.0;
     WWTree->LeptonProjection_type0 = (LEP.Pt()*cos(LEP.Theta()-(LEP + NU0_puppi).Theta()))/(LEP + NU0_puppi).Pt();
     WWTree->LeptonProjection_type2 = (LEP.Pt()*cos(LEP.Theta()-(LEP + NU2_puppi).Theta()))/(LEP + NU2_puppi).Pt();
     WWTree->LeptonProjection_run2 = (LEP.Pt()*cos(LEP.Theta()-(LEP + NU1_puppi).Theta()))/(LEP + NU1_puppi).Pt();

    
    /////////////////FILL THE TREE
    outTree->Fill();
  }
  std::cout << "---------end loop on events------------" << std::endl;
  std::cout << std::endl;
  
  std::cout << "----------------------" << std::endl;
  std::cout << " SUMMARY" << std::endl;
  std::cout << "----------------------" << std::endl;
  std::cout << std::endl;
  std::cout<<"MC matching: "<<(float)ok/(float)total<<std::endl;
  std::cout<<"negative events: "<<nNegEvents<<std::endl;
  std::cout << std::endl;
  std::cout<<"(0) all events:        "<<cutEff[0]<<std::endl
	   <<"(1) tight lepton:      "<<cutEff[1]<<std::endl
	   <<"(2) loose lepton veto: "<<cutEff[2]<<std::endl
	   <<"(3) MET:               "<<cutEff[3]<<std::endl
	   <<"(4) negative lep-MET:  "<<cutEff[4]<<std::endl
	   <<"(5) 1 good AK8:        "<<cutEff[5]<<std::endl
	   <<"(6) 2 good AK4:        "<<cutEff[6]<<std::endl
	   <<"(7) 1 AK8 or 2 AK4:    "<<cutEff[7]<<std::endl;
  
  //--------close everything-------------
  ReducedTree->fChain->Delete();
  outTree->Write();
  outROOT->Close();
  int t1 = time(NULL);
  printf ("time to run this code = %d secs\n", t1 - t0);
  return(0);
}


//////////////////////////////////
//Ref: https://github.com/ram1123/LHEAnalyzer/blob/LHEanalyzer/LHEanalyzer.cpp
//////////////////////////////////

void computeAngles(TLorentzVector thep4H, TLorentzVector thep4Z1, TLorentzVector thep4M11, TLorentzVector thep4M12, 
					  TLorentzVector thep4Z2, 
		  double& costheta1,  double& costhetastar, double& Phi1)
{
    ///////////////////////////////////////////////
    // check for z1/z2 convention, redefine all 4 vectors with convention
    ///////////////////////////////////////////////	
    TLorentzVector p4H, p4Z1, p4M11, p4M12, p4Z2;
    p4H = thep4H;
    
    p4Z1 = thep4Z1; p4M11 = thep4M11; p4M12 = thep4M12;
    p4Z2 = thep4Z2;
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
    	p4M11_BV1.Boost( boostV1 );
	p4M12_BV1.Boost( boostV1 );
    	TLorentzVector p4V2_BV1 (p4Z2);
	p4V2_BV1.Boost(boostV1);
	//cout<<"Boost = "<<p4V2_BV1.Vect().Mag()<<endl;
    //// costheta1
    
    	costheta1 = -p4V2_BV1.Vect().Dot( p4M11_BV1.Vect() )/p4V2_BV1.Vect().Mag()/p4M11_BV1.Vect().Mag();
	//cout<<costheta1<<endl;

    /* 
    //// --------------------------- costheta2
    TVector3 boostV2 = -(thep4Z2.BoostVector());
    TLorentzVector p4M11_BV2( p4M11 );
	TLorentzVector p4M12_BV2( p4M12 );	
    //TLorentzVector p4M21_BV2( p4M21 );
//	TLorentzVector p4M22_BV2( p4M22 );
    p4M11_BV2.Boost( boostV2 );
	p4M12_BV2.Boost( boostV2 );
//	p4M21_BV2.Boost( boostV2 );
//	p4M22_BV2.Boost( boostV2 );
    
    TLorentzVector p4V1_BV2 = p4M11_BV2 + p4M12_BV2;
    //// costheta2
//    costheta2 = -p4V1_BV2.Vect().Dot( p4M21_BV2.Vect() )/p4V1_BV2.Vect().Mag()/p4M21_BV2.Vect().Mag();
    */


    //// --------------------------- Phi and Phi1
    //    TVector3 boostX = -(thep4H.BoostVector());
    TLorentzVector p4M11_BX( p4M11 );
	TLorentzVector p4M12_BX( p4M12 );	
 //   TLorentzVector p4M21_BX( p4M21 );
//	TLorentzVector p4M22_BX( p4M22 );	
    
	p4M11_BX.Boost( boostX );
	p4M12_BX.Boost( boostX );
//	p4M21_BX.Boost( boostX );
//	p4M22_BX.Boost( boostX );
    
    TVector3 tmp1 = p4M11_BX.Vect().Cross( p4M12_BX.Vect() );
    //TVector3 tmp2 = tmp1;
  //  TVector3 tmp2 = p4M21_BX.Vect().Cross( p4M22_BX.Vect() );    
    
    TVector3 normal1_BX( tmp1.X()/tmp1.Mag(), tmp1.Y()/tmp1.Mag(), tmp1.Z()/tmp1.Mag() ); 
    //TVector3 normal2_BX( tmp2.X()/tmp2.Mag(), tmp2.Y()/tmp2.Mag(), tmp2.Z()/tmp2.Mag() ); 
    
    //// Phi
    TLorentzVector p4Z1_BX = p4M11_BX + p4M12_BX;    
   // double tmpSgnPhi = p4Z1_BX.Vect().Dot( normal1_BX.Cross( normal2_BX) );
    //double sgnPhi = tmpSgnPhi/fabs(tmpSgnPhi);
    //Phi = sgnPhi * acos( -1.*normal1_BX.Dot( normal2_BX) );
    
    
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
