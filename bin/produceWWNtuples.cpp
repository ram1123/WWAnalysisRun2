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
   //int VBFSel  = atoi(argv[13]);

   std::string leptonName;

   std::string iHLTFile="${CMSSW_BASE}/src/BaconAna/DataFormats/data/HLTFile_25ns";
   const std::string cmssw_base = getenv("CMSSW_BASE");
   std::string cmssw_base_env = "${CMSSW_BASE}";
   size_t start_pos = iHLTFile.find(cmssw_base_env);
   if(start_pos != std::string::npos) {
      iHLTFile.replace(start_pos, cmssw_base_env.length(), cmssw_base);
   }

   const baconhep::TTrigger triggerMenu(iHLTFile);  
   std::cout<<"Apply trigger: "<<applyTrigger<<std::endl;

   TLorentzVector LEP1, LEP2, LEP1_Up, LEP1_Down, LEP2_Up, LEP2_Down;
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
   TClonesArray *lheWgtArr	= new TClonesArray("baconhep::TLHEWeight");

   char list1[2000];
   sprintf (list1, "%s", inputFile.c_str());
   ifstream rootList (list1);

   int fileCounter=0;
   // vector to store root file names
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

   // Array to store the number of events left after each cuts applied in this macro
   int cutEff[21]={0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};

   // Read and add pileup in histogram
   TFile* pileupFileMC = TFile::Open("inputfiles/puWeights_90x_41ifb.root");
   TH1D* puWeights = (TH1D*)pileupFileMC->Get("puWeights");
   TH1D* puWeightsUp = (TH1D*)pileupFileMC->Get("puWeightsUp");
   TH1D* puWeightsDown = (TH1D*)pileupFileMC->Get("puWeightsDown");
   puWeights->SetBins(100,0,100);
   puWeightsUp->SetBins(100,0,100);
   puWeightsDown->SetBins(100,0,100);

   //---------------- Root Files for ID, ISO, Trigger, GSF correctiosn for ELE and MU both: Starts -------------
   TFile* IDIsoEle   = TFile::Open("inputfiles/egammaEffi.txt_EGM2D_runBCDEF_passingTight94X.root","READ");
   TH1F *hIDIsoEle   = (TH1F*)IDIsoEle->Get("EGamma_SF2D");

   TFile* GSFCorrEle = TFile::Open("inputfiles/egammaEffi.txt_EGM2D_runBCDEF_passingRECO.root","READ");
   TH1F *hGSFCorrEle = (TH1F*)GSFCorrEle->Get("EGamma_SF2D");

   //TFile* TriggerEle = TFile::Open("ElectronTrigger_SF.root","READ");
   //TH1F* hTriggerEle = (TH1F*)TriggerEle->Get("HLT_Ele27");

   TFile* IDMu = TFile::Open("inputfiles/RunBCDEF_SF_ID.root","READ");
   TH1F *hIDMu = (TH1F*)IDMu->Get("NUM_TightID_DEN_genTracks_pt_abseta");

   TFile* IsoMu = TFile::Open("inputfiles/RunBCDEF_SF_ISO.root","READ");
   TH1F *hIsoMu = (TH1F*)IsoMu->Get("NUM_TightRelIso_DEN_TightIDandIPCut_pt_abseta");

   TFile* TriggerMu = TFile::Open("inputfiles/EfficienciesAndSF_RunBtoF_Nov17Nov2017.root","READ");
   TH1F* hTriggerMu = (TH1F*)TriggerMu->Get("IsoMu27_PtEtaBins/pt_abseta_ratio");
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

   if (isLocal==1) nInputFiles = 11;
   cout<<"==> Total number of input files : "<<nInputFiles<<endl;

   Long64_t TotalNumberOfEvents = 0, nNegEvents = 0; 

   /*! 
    *  This for loop on the number of input file is only for calculating
    *  the total number of events and the total number of negative events
    *  in the root files.
    */
   for(int InputFiles=0; InputFiles<nInputFiles; InputFiles++)
   {
      infile = TFile::Open(sampleName[InputFiles]);
      if (infile->IsZombie()) continue;
      if (!(infile->GetListOfKeys()->Contains("Events"))) continue;
      eventTree = (TTree*)infile->Get("Events");
      //std::cout << "\t File no. " << i << "\t"<< sampleName[i] <<std:: endl;

      TotalNumberOfEvents+=eventTree->GetEntries();
      if(isMC)
      {
	 eventTree->SetBranchAddress("Info", &info);    
	 TBranch *infoBr = eventTree->GetBranch("Info");

	 TBranch *genBr=0;
	 eventTree->SetBranchAddress("GenEvtInfo", &gen); 
	 genBr = eventTree->GetBranch("GenEvtInfo");

	 for (Long64_t jentry=0; jentry<eventTree->GetEntries();jentry++,jentry2++)
	 {
	    genBr->GetEntry(jentry);
	    infoBr->GetEntry(jentry);	    
	    //if (jentry2%10000 == 0) std::cout << "\t File no. " << InputFiles << "; Neg Event Count; read entry: " << jentry2 <<"/"<<TotalNumberOfEvents<<std:: endl;
	    if (jentry2%10000 == 0) std::cout << "\t File no. " << InputFiles << " :" << sampleName[InputFiles] << "; Neg Event Count; read entry: " << jentry2 <<"/"<<TotalNumberOfEvents<<std:: endl;
	    if (gen->weight<0)	nNegEvents++;
	 }
      }
      delete infile;
      infile=0, eventTree=0;
   }
   cout<<"==> Total number of events : "<<TotalNumberOfEvents<<endl;
   cout<<"==> Total number of negative events : "<<nNegEvents<<endl;

   // Get the weight
   float weight = std::atof(xSecWeight.c_str())/(std::atof(TotalNumberOfEntries.c_str()) - 2*std::atof(TotalNumberOfNegativeEntries.c_str()));
   cout<<"Weight of cross-sec/events = "<<weight<<endl;
   int totalEntries=0;

   //---------start loop on events------------
   std::cout << "---------start loop on events------------" << std::endl;
   jentry2=0;
   for(int i=0;i<nInputFiles;i++)
   {
      cout<<"\n\n=====	Processing File Number : "<<i<<"/"<<nInputFiles<<"\n\t"<<sampleName[i]<<"\n-------"<<endl;

      infile = TFile::Open(sampleName[i]);
      if (infile->IsZombie()) continue;
      if (!(infile->GetListOfKeys()->Contains("Events"))) continue;
      eventTree = (TTree*)infile->Get("Events");

      totalEntries	+= eventTree->GetEntries();

      nEvents=eventTree->GetEntries();
      cout<<"\t==> Entries = "<<nEvents<<endl;

      eventTree->SetBranchAddress("Info", &info);    TBranch *infoBr = eventTree->GetBranch("Info");
      eventTree->SetBranchAddress("Muon", &muonArr); TBranch *muonBr = eventTree->GetBranch("Muon");
      eventTree->SetBranchAddress("Electron", &electronArr); TBranch *electronBr = eventTree->GetBranch("Electron");
      eventTree->SetBranchAddress("PV",   &vertexArr); TBranch *vertexBr = eventTree->GetBranch("PV");
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
	    if(lhePartBr)	lhePartBr->GetEntry(jentry);
	    genPartArr->Clear();
	    genBr->GetEntry(jentry);
	    genPartBr->GetEntry(jentry);

	    //TLorentzVector lepW, temp;
	    TLorentzVector genLep;//, genNeutrino;
	    //std::vector<TLorentzVector> v_GEN_lepW, v_GEN_temp;
	    std::vector<TLorentzVector> v_genLep;	//,v_genNeutrino;

	    //v_GEN_lepW.clear();	
	    //v_GEN_temp.clear();	
	    v_genLep.clear();	

	    for (int i = 0; i<genPartArr->GetEntries();i++)
	    {
	       const baconhep::TGenParticle* genloop = (baconhep::TGenParticle*) ((*genPartArr)[i]);
	       Int_t parentPdg=dynamic_cast<baconhep::TGenParticle *>(genPartArr->At(genloop->parent>-1 ? genloop->parent : 0))->pdgId;
	       if( (abs(genloop->pdgId) == 11 || abs(genloop->pdgId) == 13 ) && abs(parentPdg) == 24)
	       {
		  genLep.SetPtEtaPhiM(genloop->pt, genloop->eta, genloop->phi, genloop->mass);
		  v_genLep.push_back(genLep);
	       }
	    }

	    if ((v_genLep.size()==1 || v_genLep.size()==2) )
	    {
	       WWTree->isGen           = 1;
	       WWTree->lep_pt_gen      = v_genLep[0].Pt();
	       WWTree->lep_eta_gen     = v_genLep[0].Eta();

	       count_genEvents++;
	    }

	    if ( WWTree->lep_pt_gen > 30 && abs(WWTree->lep_eta_gen) < 2.5 )
	    {
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
	 else if (gen->weight<0) 
	 {
	    WWTree->genWeight=-1.;
	    //nNegEvents++;
	 }
	 cutEff[0]++;

	 if (isMC==1)
	    if (GenPassCut == 1)   cutEff[1]++;


	 vertexArr->Clear();
	 vertexBr->GetEntry(jentry);
	 WWTree->nPV = vertexArr->GetEntries();

	 //PILE-UP WEIGHT
	 if (isMC==1) 
	 {
	    if(int(info->nPUmean)<100)
	    {
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
	    if(!(triggerMenu.pass("HLT_IsoMu27_v*",info->triggerBits) || triggerMenu.pass("HLT_IsoTkMu27_v*",info->triggerBits) ||  triggerMenu.pass("HLT_Ele35_WPTight_Gsf_v*",info->triggerBits))) continue;

	 /////////////////THE SELECTED LEPTON
	 GeneralizedEndpoint ScaleSystematic;
	 int nTightEle=0, nLooseEle=0;
	 int nTightMu=0, nLooseMu=0;
	 double pt_cut = 20;
	 double leadelept_cut = 40;
	 double leadmupt_cut = 30;
	 electronArr->Clear();
	 electronBr->GetEntry(jentry);
	 const baconhep::TElectron *leadele = NULL;
	 const baconhep::TElectron *subele = NULL;
	 double leadeleE=-999, subeleE=-999;
	 double iso = 1.5;
	 for (int i=0; i<electronArr->GetEntries(); i++) 
	 {
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
	 for(Int_t i=0; i<muonArr->GetEntries(); i++) 
	 {
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
	    } else
	    {
	       WWTree->l_pt1      = leadmu->pt;
	       WWTree->l_pt1_Up   = leadmu->pt + 0.01*leadmu->pt;
	       WWTree->l_pt1_Down = leadmu->pt - 0.01*leadmu->pt;
	    }
	    WWTree->l_eta1 = leadmu->eta;
	    WWTree->l_phi1 = leadmu->phi;	
	    WWTree->l_e1 = leadmue;	
	    WWTree->l_charge1 = leadmu->q;
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
	 }	else{
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
	    //std::cout<< "Start Reading..." << endl;
	    WWTree->id_eff_Weight = GetSFs_Lepton(WWTree->l_eta1, WWTree->l_pt1, hIDIsoEle, "eta", "pt");	// Get Scale factor corresponding to the pt and eta.
	    //std::cout << "Ele1: " << WWTree->l_eta1 << "\t" << WWTree->l_pt1 << "\t" << WWTree->id_eff_Weight << std::endl;

	    // apply GSF/RECO SF's for electrons
	    WWTree->id_eff_Weight = WWTree->id_eff_Weight*GetSFs_Lepton(WWTree->l_eta1, WWTree->l_pt1, hGSFCorrEle, "eta", "pt");
	    //WWTree->trig_eff_Weight = (GetSFs_Lepton(WWTree->l_pt1, WWTree->l_eta1, hTriggerEle, "eta", "pt"));

	    if(WWTree->l_pt2>0){
	       //apply ID, ISO SF's
	       WWTree->id_eff_Weight2 = GetSFs_Lepton(WWTree->l_eta2, WWTree->l_pt2, hIDIsoEle, "eta", "pt");	// Get Scale factor corresponding to the pt and eta.

	       // apply GSF/RECO SF's for electrons
	       WWTree->id_eff_Weight2 = WWTree->id_eff_Weight2*GetSFs_Lepton(WWTree->l_eta2, WWTree->l_pt2, hGSFCorrEle, "eta", "pt");
	       //WWTree->trig_eff_Weight2 = 1.0/(GetSFs_Lepton(WWTree->l_pt2, WWTree->l_eta2, hTriggerEle, "eta", "pt"));
	    }
	 }

	 //ID&ISO efficiency SF for muons (https://twiki.cern.ch/twiki/bin/view/CMS/MuonWorkInProgressAndPagResults)
	 if (strcmp(leptonName.c_str(),"mu")==0 && isMC==1) {
	    //apply ID SF's
	    WWTree->id_eff_Weight = GetSFs_Lepton(WWTree->l_pt1, abs(WWTree->l_eta1), hIDMu, "pt", "eta");
	    WWTree->id_eff_Weight_Up = GetSFs_Lepton(WWTree->l_pt1_Up, abs(WWTree->l_eta1), hIDMu, "pt", "eta");
	    WWTree->id_eff_Weight_Down = GetSFs_Lepton(WWTree->l_pt1_Down, abs(WWTree->l_eta1), hIDMu, "pt", "eta");
	    if (WWTree->l_pt2>0) {
	       WWTree->id_eff_Weight2 = GetSFs_Lepton(WWTree->l_pt2, abs(WWTree->l_eta2), hIDMu, "pt", "eta");
	       WWTree->id_eff_Weight2_Up = GetSFs_Lepton(WWTree->l_pt2_Up, abs(WWTree->l_eta2), hIDMu, "pt", "eta");
	       WWTree->id_eff_Weight2_Down = GetSFs_Lepton(WWTree->l_pt2_Down, abs(WWTree->l_eta2), hIDMu, "pt", "eta");
	    }
	    //  apply ISO SF's
	    WWTree->id_eff_Weight = WWTree->id_eff_Weight*GetSFs_Lepton(WWTree->l_pt1, abs(WWTree->l_eta1), hIsoMu, "pt", "eta");
	    WWTree->id_eff_Weight_Up = WWTree->id_eff_Weight_Up*GetSFs_Lepton(WWTree->l_pt1_Up, abs(WWTree->l_eta1), hIsoMu, "pt", "eta");
	    WWTree->id_eff_Weight_Down = WWTree->id_eff_Weight_Down*GetSFs_Lepton(WWTree->l_pt1_Down, abs(WWTree->l_eta1), hIsoMu, "pt", "eta");
	    if (WWTree->l_pt2>0) {
	       WWTree->id_eff_Weight2 = WWTree->id_eff_Weight2*GetSFs_Lepton(WWTree->l_pt2, abs(WWTree->l_eta2), hIsoMu, "pt", "eta");
	       WWTree->id_eff_Weight2_Up = WWTree->id_eff_Weight2_Up*GetSFs_Lepton(WWTree->l_pt2_Up, abs(WWTree->l_eta2), hIsoMu, "pt", "eta");
	       WWTree->id_eff_Weight2_Down = WWTree->id_eff_Weight2_Down*GetSFs_Lepton(WWTree->l_pt2_Down, abs(WWTree->l_eta2), hIsoMu, "pt", "eta");
	    }

	    // apply Trigger SF's
	    WWTree->trig_eff_Weight  = GetSFs_Lepton(WWTree->l_pt1, abs(WWTree->l_eta1), hTriggerMu, "pt", "eta");
	    WWTree->trig_eff_Weight_Up  = GetSFs_Lepton(WWTree->l_pt1_Up, abs(WWTree->l_eta1_Up), hTriggerMu, "pt", "eta");
	    WWTree->trig_eff_Weight_Down  = GetSFs_Lepton(WWTree->l_pt1_Down, abs(WWTree->l_eta1_Down), hTriggerMu, "pt", "eta");
	    if (WWTree->l_pt2>0) {
	       WWTree->trig_eff_Weight2 = GetSFs_Lepton(WWTree->l_pt2, abs(WWTree->l_eta2), hTriggerMu, "pt", "eta");
	       WWTree->trig_eff_Weight2_Up = GetSFs_Lepton(WWTree->l_pt2_Up, abs(WWTree->l_eta2), hTriggerMu, "pt", "eta");
	       WWTree->trig_eff_Weight2_Down = GetSFs_Lepton(WWTree->l_pt2_Down, abs(WWTree->l_eta2), hTriggerMu, "pt", "eta");
	    }
	 }

	 //if(LEP1.Pt()<=0 || LEP2.Pt() <= 0 ){ std::cerr<<" Negative Lepton - Neutrino Pt "<<std::endl; continue ; }
	 if (WWTree->l_pt2>0) cutEff[14]++;  // There is no MET in two lepton case. So, cutEff[4] is placed in previous if condition.
 /////////////////////////////////////////////////////////////////////////////////////////	END btag weight calculation
	 if (WWTree->l_pt2<0) 
	    WWTree->totalEventWeight = WWTree->genWeight*WWTree->pu_Weight*WWTree->top1_NNLO_Weight*WWTree->top2_NNLO_Weight*WWTree->trig_eff_Weight*WWTree->id_eff_Weight;
	 else
	    WWTree->totalEventWeight_2Lep = WWTree->genWeight*WWTree->pu_Weight*WWTree->top1_NNLO_Weight*WWTree->top2_NNLO_Weight*WWTree->trig_eff_Weight*WWTree->id_eff_Weight*WWTree->trig_eff_Weight2*WWTree->id_eff_Weight2;
	 WWTree->nEvents = TotalNumberOfEvents;
	 WWTree->nNegEvents = nNegEvents;
	 WWTree->nTotEvents = std::atof(TotalNumberOfEntries.c_str());
	 WWTree->nTotNegEvents = std::atof(TotalNumberOfNegativeEntries.c_str());


	 outTree->Fill();
	 //goodJetsv.clear();
      }
      delete infile;
      infile=0, eventTree=0;
      /////////////////FILL THE TREE
   }
   delete puWeightsDown;	delete puWeightsUp;	delete puWeights;
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
   std::cout	<<"(0) all events:        "<<cutEff[0]<<"\t:\t"<<((float)cutEff[0]*100.0)/(float)cutEff[0]<<std::endl
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
   std::cout	<<"(0) all events:        "<<cutEff[0]<<"\t:\t"<<((float)cutEff[0]*100.0)/(float)cutEff[0]<<std::endl
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
   delete lheWgtArr;
   outROOT->Write();
   outROOT->Close();
   int t1 = time(NULL);
   printf ("\n==> time to run this code = %0.3f min\n", (float)(t1 - t0)/60.0);
   return(0);
}
