//============================================================================================
// B-tagging Uncertainty functions:
//-----------------------------
//  - jet b-tagging SFs
//  - b-tag MC efficiencies for b-jets, c-jets, and light flavor(u/d/s/g) jets
//  - b-tag event weight depending on # of selected b-tags in event
//============================================================================================

#ifndef BTAGUNC_HH
#define BTAGUNC_HH


// bacon object headers
#include "BaconAna/DataFormats/interface/TJet.hh"

//ROOT headers
#include <TLorentzVector.h>
#include <TMath.h>

//C++ headers
#include <vector>
#include <iostream>
#include <cmath>
#include <string>
#include <cassert>

// B-tag calibration and SF headers
#include "CondFormats/BTauObjects/interface/BTagEntry.h"
#include "CondFormats/BTauObjects/interface/BTagCalibration.h"
#include "CondTools/BTau/interface/BTagCalibrationReader.h"

using namespace std;

//=== FUNCTION DECLARATIONS ======================================================================================
vector<float> getJetSFs(std::vector<const baconhep::TJet*> vGoodJet, BTagCalibrationReader &reader, string sysTypeHF, string sysTypeLF);
float getBtagEventReweight(int NminBjets, std::vector <const baconhep::TJet*> &vJet, std::vector <float> SF);
float getBtagEventReweightEtaBin(int NminBjets, std::vector <const baconhep::TJet*> &vJet, std::vector <float> SF);
//=== FUNCTION DEFINITIONS ======================================================================================

vector<float> getJetSFs(std::vector<const baconhep::TJet*> vGoodJet, BTagCalibrationReader &reader, string sysTypeHF, string sysTypeLF)
{

  std::vector<float> vSF;
  // std::vector<float> vSFerr;
  
  for(unsigned int ijet=0; ijet<vGoodJet.size(); ++ijet){
    
    bool isB=false; bool isC=false; bool isLF=false;

    if(fabs(vGoodJet[ijet]->hadronFlavor)==5){      isB=true;}//  std::cout << "isB" << std::endl;}
    else if(fabs(vGoodJet[ijet]->hadronFlavor)==4){ isC=true;}//  std::cout << "isC" << std::endl;}
    else if(fabs(vGoodJet[ijet]->hadronFlavor)==0){ isLF=true;}// std::cout << "isLF" << std::endl;}
     
     float SF=0;
     if(isB){
       SF = reader.eval_auto_bounds(sysTypeHF, BTagEntry::FLAV_B, vGoodJet[ijet]->eta, vGoodJet[ijet]->pt);
     } else if (isC){
       SF = reader.eval_auto_bounds(sysTypeHF, BTagEntry::FLAV_C,vGoodJet[ijet]->eta, vGoodJet[ijet]->pt);
     } else if (isLF){
       SF = reader.eval_auto_bounds(sysTypeLF, BTagEntry::FLAV_UDSG,vGoodJet[ijet]->eta, vGoodJet[ijet]->pt);
     }

     //  SFerr = fabs(SF_UP - SF) > fabs(SF_DN - SF) ? fabs(SF_DN - SF):fabs(SF_UP - SF);
     vSF.push_back(SF);
     //  vSFerr.push_back(SFerr);
  }
  return vSF;
}
//-------------------------------------------------------------------------------------------------------------------------
float getBtagEventReweight(int NminBjets, std::vector <const baconhep::TJet*> &vJet, std::vector <float> SF){
 

  if(vJet.size()==0) return 1;
  
  
  // float err1 = 0; 
  // float err2 = 0;
  
  float wtbtag;
  // const float wtbtagErr;
  
  std::vector<float> bjet_ptbin_eff, cjet_ptbin_eff, ljet_ptbin_eff;
  // calculation done in all dilepton ttbar (powheg)
  bjet_ptbin_eff.push_back(0.680421);   cjet_ptbin_eff.push_back(0.156886);   ljet_ptbin_eff.push_back(0.0253288); 
  bjet_ptbin_eff.push_back(0.713573);   cjet_ptbin_eff.push_back(0.153847);   ljet_ptbin_eff.push_back(0.0196695);
  bjet_ptbin_eff.push_back(0.728557);   cjet_ptbin_eff.push_back(0.160366);   ljet_ptbin_eff.push_back(0.0190403);
  bjet_ptbin_eff.push_back(0.727522);   cjet_ptbin_eff.push_back(0.165694);   ljet_ptbin_eff.push_back(0.0196844);
  bjet_ptbin_eff.push_back(0.718537);   cjet_ptbin_eff.push_back(0.170418);   ljet_ptbin_eff.push_back(0.0247658);
  bjet_ptbin_eff.push_back(0.682109);   cjet_ptbin_eff.push_back(0.170722);   ljet_ptbin_eff.push_back(0.0266844);
  bjet_ptbin_eff.push_back(0.621643);   cjet_ptbin_eff.push_back(0.179256);   ljet_ptbin_eff.push_back(0.03663);
  
  std::vector<float> ptbinlow, ptbinhigh;
  ptbinlow.push_back(30);  ptbinhigh.push_back(50);
  ptbinlow.push_back(50);  ptbinhigh.push_back(70);
  ptbinlow.push_back(70);  ptbinhigh.push_back(100);
  ptbinlow.push_back(100); ptbinhigh.push_back(140);
  ptbinlow.push_back(140); ptbinhigh.push_back(200);
  ptbinlow.push_back(200); ptbinhigh.push_back(300);
  ptbinlow.push_back(300); ptbinhigh.push_back(670);
  
  
  //if jet pt < 30 things are messed up...   
  
  float mcProd_0 = 1.;
  float dataProd_0 = 1.;
  
  float mcProd_1 = 1.;
  float dataProd_1 = 1.;
  
  float mcSum = 0.;
  float dataSum = 0.;
  
  for(unsigned int ij=0; ij<vJet.size(); ++ij){
    float mcTag = 1.;
    float mcNoTag = 0.;
    float dataTag = 1.;
    float dataNoTag = 0.;
    mcProd_1 = 1.; dataProd_1 = 1.;
    
     for(unsigned int ipt=0; ipt<ptbinhigh.size(); ++ipt){
       if(vJet[ij]->pt > 670){
	 if(vJet[ij]->hadronFlavor == 5)       { mcTag  = bjet_ptbin_eff[ipt]; mcNoTag  = (1 - bjet_ptbin_eff[ipt]); dataTag  = bjet_ptbin_eff[ipt]*SF[ij]; dataNoTag  = (1 - bjet_ptbin_eff[ipt]*SF[ij]); 
	 } else if(vJet[ij]->hadronFlavor == 4){ mcTag  = cjet_ptbin_eff[ipt]; mcNoTag  = (1 - cjet_ptbin_eff[ipt]); dataTag  = cjet_ptbin_eff[ipt]*SF[ij]; dataNoTag  = (1 - cjet_ptbin_eff[ipt]*SF[ij]);
	 } else if(vJet[ij]->hadronFlavor == 0){ mcTag  = ljet_ptbin_eff[ipt]; mcNoTag  = (1 - ljet_ptbin_eff[ipt]); dataTag  = ljet_ptbin_eff[ipt]*SF[ij]; dataNoTag  = (1 - ljet_ptbin_eff[ipt]*SF[ij]);
	 }
      }
       if(vJet[ij]->pt > ptbinlow[ipt] && vJet[ij]->pt <= ptbinhigh[ipt]){
	if(vJet[ij]->hadronFlavor == 5)       { mcTag  = bjet_ptbin_eff[ipt]; mcNoTag  = (1 - bjet_ptbin_eff[ipt]); dataTag  = bjet_ptbin_eff[ipt]*SF[ij]; dataNoTag  = (1 - bjet_ptbin_eff[ipt]*SF[ij]);
	} else if(vJet[ij]->hadronFlavor == 4){ mcTag  = cjet_ptbin_eff[ipt]; mcNoTag  = (1 - cjet_ptbin_eff[ipt]); dataTag  = cjet_ptbin_eff[ipt]*SF[ij]; dataNoTag  = (1 - cjet_ptbin_eff[ipt]*SF[ij]);
	} else if(vJet[ij]->hadronFlavor == 0){ mcTag  = ljet_ptbin_eff[ipt]; mcNoTag  = (1 - ljet_ptbin_eff[ipt]); dataTag  = ljet_ptbin_eff[ipt]*SF[ij]; dataNoTag  = (1 - ljet_ptbin_eff[ipt]*SF[ij]);
	}
      }
    }
    
     mcProd_0   *= mcNoTag;
     dataProd_0 *= dataNoTag;

    for(unsigned int ik=0; ik < vJet.size(); ++ik){
      if(ik == ij) continue;
      for(unsigned int jpt=0; jpt<ptbinhigh.size(); ++jpt){
	if(vJet[ik]->pt > 670){
	  if(vJet[ik]->hadronFlavor == 5)    {    mcNoTag = (1 - bjet_ptbin_eff[jpt]); dataNoTag  = (1 - bjet_ptbin_eff[jpt]*SF[ik]); 
	  } else if(vJet[ik]->hadronFlavor == 4){ mcNoTag = (1 - cjet_ptbin_eff[jpt]); dataNoTag  = (1 - cjet_ptbin_eff[jpt]*SF[ik]);
	  } else if(vJet[ik]->hadronFlavor == 0){ mcNoTag = (1 - ljet_ptbin_eff[jpt]); dataNoTag  = (1 - ljet_ptbin_eff[jpt]*SF[ik]);
	  }
	}
	if(vJet[ik]->pt > ptbinlow[jpt] && vJet[ik]->pt <= ptbinhigh[jpt]){
	  if(vJet[ik]->hadronFlavor == 5)     {    mcNoTag  = (1 - bjet_ptbin_eff[jpt]);  dataNoTag  = (1 - bjet_ptbin_eff[jpt]*SF[ik]);
	  } else if(vJet[ik]->hadronFlavor == 4){  mcNoTag  = (1 - cjet_ptbin_eff[jpt]);  dataNoTag  = (1 - cjet_ptbin_eff[jpt]*SF[ik]);
	  } else if(vJet[ik]->hadronFlavor == 0){  mcNoTag  = (1 - ljet_ptbin_eff[jpt]);  dataNoTag  = (1 - ljet_ptbin_eff[jpt]*SF[ik]);
	  }
	}
      }
      mcProd_1   *= mcNoTag;
      dataProd_1 *= dataNoTag;
    }
    
    mcSum   += mcTag*mcProd_1;
    dataSum += dataTag*dataProd_1;
    
    
  }

  
  // re-weighting for events with 0 jets
  if(NminBjets==0){
    /*
    std::cout<<"0 btag"<<std::endl;
    std::cout<<"dataProd_0: " << dataProd_0 << std::endl; 
    std::cout<<"mcProd_0:   " << mcProd_0 << std::endl; 
    */

    return wtbtag = dataProd_0/mcProd_0;
    
    //re-weighting for events with exactly 1 b-tag jet
  } else if(NminBjets==-1){
    /*
    std::cout<<"1 btag" <<std::endl;
    std::cout<<"dataSum: " << dataSum << std::endl;
    std::cout<<"mcSum:   " << mcSum << std::endl;
    */
    return wtbtag = dataSum/mcSum;
    
    // re-weighting for events with 1 or more b-tag jets
  } else if(NminBjets==1){
    /*
    std::cout<<">=1 btag"<<std::endl;
    std::cout<<"dataProd: "<<(1 - dataProd_0) <<std::endl;
    std::cout<<"mcProd:   "<<(1 - mcProd_0) <<std::endl;
    */
    return wtbtag = (1 - dataProd_0)/(1 - mcProd_0);
    // re-weighting for events with 2 or more b-tag jets
  } else if(NminBjets==2){
    /*
    std::cout<<">=2 btag"<<std::endl;
    std::cout<<"1 - dataProd_0 - dataSum: " << 1 - dataProd_0 - dataSum <<std::endl;
    std::cout<<"1 - mcProd_0 - mcSum:     " << 1 - mcProd_0 - mcSum <<std::endl;
    */
    return wtbtag = (1 - dataProd_0 - dataSum)/(1 - mcProd_0 - mcSum);
    
  } 
  
  else return 1;
   
  
  
  //return wtbtagErr = sqrt(2*(pow(err1+err2,2)))*wtbtag;
  
}
//-------------------------------------------------------------------------------------------------------------------------
float getBtagEventReweightEtaBin(int NminBjets, std::vector <const baconhep::TJet*> &vJet, std::vector <float> SF){
  
  if(vJet.size()==0) return 1;
  
  // float err1 = 0; 
  // float err2 = 0;
  
  float wtbtag;
  // const float wtbtagErr;
  
  std::vector<float> bjet_ptbin_eff[3], cjet_ptbin_eff[3], ljet_ptbin_eff[3];
  // calculation done in all dilepton ttbar (powheg)
  bjet_ptbin_eff[0].push_back(0.707814);   cjet_ptbin_eff[0].push_back(0.15622);    ljet_ptbin_eff[0].push_back(0.0169301); 
  bjet_ptbin_eff[0].push_back(0.738506);   cjet_ptbin_eff[0].push_back(0.150696);   ljet_ptbin_eff[0].push_back(0.0120694);
  bjet_ptbin_eff[0].push_back(0.751684);   cjet_ptbin_eff[0].push_back(0.158773);   ljet_ptbin_eff[0].push_back(0.011405);
  bjet_ptbin_eff[0].push_back(0.750919);   cjet_ptbin_eff[0].push_back(0.164317);   ljet_ptbin_eff[0].push_back(0.0108449);
  bjet_ptbin_eff[0].push_back(0.74029);    cjet_ptbin_eff[0].push_back(0.155843);   ljet_ptbin_eff[0].push_back(0.0118652);
  bjet_ptbin_eff[0].push_back(0.70986);    cjet_ptbin_eff[0].push_back(0.157589);   ljet_ptbin_eff[0].push_back(0.0141109);
  bjet_ptbin_eff[0].push_back(0.645886);   cjet_ptbin_eff[0].push_back(0.157835);   ljet_ptbin_eff[0].push_back(0.0200094);

  bjet_ptbin_eff[1].push_back(0.697665);   cjet_ptbin_eff[1].push_back(0.166117);   ljet_ptbin_eff[1].push_back(0.0237088); 
  bjet_ptbin_eff[1].push_back(0.728531);   cjet_ptbin_eff[1].push_back(0.163026);   ljet_ptbin_eff[1].push_back(0.0182443);
  bjet_ptbin_eff[1].push_back(0.742855);   cjet_ptbin_eff[1].push_back(0.167246);   ljet_ptbin_eff[1].push_back(0.0175251);
  bjet_ptbin_eff[1].push_back(0.740688);   cjet_ptbin_eff[1].push_back(0.17496);    ljet_ptbin_eff[1].push_back(0.0176429);
  bjet_ptbin_eff[1].push_back(0.729423);   cjet_ptbin_eff[1].push_back(0.18088);    ljet_ptbin_eff[1].push_back(0.0209487);
  bjet_ptbin_eff[1].push_back(0.69149);    cjet_ptbin_eff[1].push_back(0.183375);   ljet_ptbin_eff[1].push_back(0.0236471);
  bjet_ptbin_eff[1].push_back(0.625445);   cjet_ptbin_eff[1].push_back(0.187591);   ljet_ptbin_eff[1].push_back(0.0359101);

  bjet_ptbin_eff[2].push_back(0.607004);   cjet_ptbin_eff[2].push_back(0.144066);   ljet_ptbin_eff[2].push_back(0.033486); 
  bjet_ptbin_eff[2].push_back(0.636906);   cjet_ptbin_eff[2].push_back(0.142641);   ljet_ptbin_eff[2].push_back(0.0272722);
  bjet_ptbin_eff[2].push_back(0.646279);   cjet_ptbin_eff[2].push_back(0.150943);   ljet_ptbin_eff[2].push_back(0.0266951);
  bjet_ptbin_eff[2].push_back(0.636208);   cjet_ptbin_eff[2].push_back(0.151915);   ljet_ptbin_eff[2].push_back(0.0287519);
  bjet_ptbin_eff[2].push_back(0.623899);   cjet_ptbin_eff[2].push_back(0.166242);   ljet_ptbin_eff[2].push_back(0.0392322);
  bjet_ptbin_eff[2].push_back(0.55471);    cjet_ptbin_eff[2].push_back(0.161047);   ljet_ptbin_eff[2].push_back(0.0404986);
  bjet_ptbin_eff[2].push_back(0.496788);   cjet_ptbin_eff[2].push_back(0.189958);   ljet_ptbin_eff[2].push_back(0.0535281);
  
  std::vector<float> ptbinlow, ptbinhigh;
  ptbinlow.push_back(30);  ptbinhigh.push_back(50);
  ptbinlow.push_back(50);  ptbinhigh.push_back(70);
  ptbinlow.push_back(70);  ptbinhigh.push_back(100);
  ptbinlow.push_back(100); ptbinhigh.push_back(140);
  ptbinlow.push_back(140); ptbinhigh.push_back(200);
  ptbinlow.push_back(200); ptbinhigh.push_back(300);
  ptbinlow.push_back(300); ptbinhigh.push_back(670);
  
  std::vector<float> etabinlow, etabinhigh;
  etabinlow.push_back(0.0); etabinhigh.push_back(0.5);
  etabinlow.push_back(0.5); etabinhigh.push_back(1.5);
  etabinlow.push_back(1.5); etabinhigh.push_back(2.4);
  
  //if jet pt < 30 things are messed up...   
  
  float mcProd_0 = 1.;
  float dataProd_0 = 1.;
  
  float mcProd_1 = 1.;
  float dataProd_1 = 1.;
  
  float mcSum = 0.;
  float dataSum = 0.;
  
  for(unsigned int ij=0; ij<vJet.size(); ++ij){
    float mcTag = 1.;
    float mcNoTag = 1.;
    float dataTag = 1.;
    float dataNoTag = 0.;
    mcProd_1 = 1.; dataProd_1 = 1.;
    for(unsigned int ieta=0; ieta<etabinhigh.size(); ++ieta){ 
      for(unsigned int ipt=0; ipt<ptbinhigh.size(); ++ipt){
	if(vJet[ij]->pt > 670){
	  if(fabs(vJet[ij]->eta) > etabinlow[ieta] && fabs(vJet[ij]->eta) < etabinhigh[ieta]){ 
	    if(vJet[ij]->hadronFlavor == 5)       { mcTag  = bjet_ptbin_eff[ieta][ipt]; mcNoTag  = (1 - bjet_ptbin_eff[ieta][ipt]); dataTag  = bjet_ptbin_eff[ieta][ipt]*SF[ij]; dataNoTag  = (1 - bjet_ptbin_eff[ieta][ipt]*SF[ij]); 
	    } else if(vJet[ij]->hadronFlavor == 4){ mcTag  = cjet_ptbin_eff[ieta][ipt]; mcNoTag  = (1 - cjet_ptbin_eff[ieta][ipt]); dataTag  = cjet_ptbin_eff[ieta][ipt]*SF[ij]; dataNoTag  = (1 - cjet_ptbin_eff[ieta][ipt]*SF[ij]);
	    } else if(vJet[ij]->hadronFlavor == 0){ mcTag  = ljet_ptbin_eff[ieta][ipt]; mcNoTag  = (1 - ljet_ptbin_eff[ieta][ipt]); dataTag  = ljet_ptbin_eff[ieta][ipt]*SF[ij]; dataNoTag  = (1 - ljet_ptbin_eff[ieta][ipt]*SF[ij]);
	    }
	  }
	}
	if(fabs(vJet[ij]->eta) > etabinlow[ieta] && fabs(vJet[ij]->eta) < etabinhigh[ieta]){ 
	  if(vJet[ij]->pt > ptbinlow[ipt] && vJet[ij]->pt <= ptbinhigh[ipt]){
	    if(vJet[ij]->hadronFlavor == 5)       { mcTag  = bjet_ptbin_eff[ieta][ipt]; mcNoTag  = (1 - bjet_ptbin_eff[ieta][ipt]); dataTag  = bjet_ptbin_eff[ieta][ipt]*SF[ij]; dataNoTag  = (1 - bjet_ptbin_eff[ieta][ipt]*SF[ij]);
	    } else if(vJet[ij]->hadronFlavor == 4){ mcTag  = cjet_ptbin_eff[ieta][ipt]; mcNoTag  = (1 - cjet_ptbin_eff[ieta][ipt]); dataTag  = cjet_ptbin_eff[ieta][ipt]*SF[ij]; dataNoTag  = (1 - cjet_ptbin_eff[ieta][ipt]*SF[ij]);
	    } else if(vJet[ij]->hadronFlavor == 0){ mcTag  = ljet_ptbin_eff[ieta][ipt]; mcNoTag  = (1 - ljet_ptbin_eff[ieta][ipt]); dataTag  = ljet_ptbin_eff[ieta][ipt]*SF[ij]; dataNoTag  = (1 - ljet_ptbin_eff[ieta][ipt]*SF[ij]);
	    }
	  }
	}
      }
    }
     mcProd_0   *= mcNoTag;
     dataProd_0 *= dataNoTag;

    for(unsigned int ik=0; ik < vJet.size(); ++ik){
      if(ik == ij) continue;
      for(unsigned int jeta=0; jeta<etabinhigh.size(); ++jeta){
	for(unsigned int jpt=0; jpt<ptbinhigh.size(); ++jpt){
	  if(fabs(vJet[ik]->eta) > etabinlow[jeta] && fabs(vJet[ik]->eta) < etabinhigh[jeta]){
	    if(vJet[ik]->pt > 670){
	      if(vJet[ik]->hadronFlavor == 5)    {    mcNoTag = (1 - bjet_ptbin_eff[jeta][jpt]); dataNoTag  = (1 - bjet_ptbin_eff[jeta][jpt]*SF[ik]); 
	      } else if(vJet[ik]->hadronFlavor == 4){ mcNoTag = (1 - cjet_ptbin_eff[jeta][jpt]); dataNoTag  = (1 - cjet_ptbin_eff[jeta][jpt]*SF[ik]);
	      } else if(vJet[ik]->hadronFlavor == 0){ mcNoTag = (1 - ljet_ptbin_eff[jeta][jpt]); dataNoTag  = (1 - ljet_ptbin_eff[jeta][jpt]*SF[ik]);
	      }
	    }
	  }
	  if(fabs(vJet[ik]->eta) > etabinlow[jeta] && fabs(vJet[ik]->eta) < etabinhigh[jeta]){
	    if(vJet[ik]->pt > ptbinlow[jpt] && vJet[ik]->pt <= ptbinhigh[jpt]){
	      if(vJet[ik]->hadronFlavor == 5)     {    mcNoTag  = (1 - bjet_ptbin_eff[jeta][jpt]);  dataNoTag  = (1 - bjet_ptbin_eff[jeta][jpt]*SF[ik]);
	      } else if(vJet[ik]->hadronFlavor == 4){  mcNoTag  = (1 - cjet_ptbin_eff[jeta][jpt]);  dataNoTag  = (1 - cjet_ptbin_eff[jeta][jpt]*SF[ik]);
	      } else if(vJet[ik]->hadronFlavor == 0){  mcNoTag  = (1 - ljet_ptbin_eff[jeta][jpt]);  dataNoTag  = (1 - ljet_ptbin_eff[jeta][jpt]*SF[ik]);
	      }
	    }
	  }
	}
      }
      mcProd_1   *= mcNoTag;
      dataProd_1 *= dataNoTag;
    }
    
    mcSum   += mcTag*mcProd_1;
    dataSum += dataTag*dataProd_1;
    
    
  }

  
  // re-weighting for events with 0 jets
  if(NminBjets==0){
    if (isnan(dataProd_0/mcProd_0)){
    std::cout<<"0 btag"<<"\t";
    std::cout<<"dataProd_0: " << dataProd_0 << "\t"; 
    std::cout<<"mcProd_0:   " << mcProd_0 << std::endl; 
    }

    return wtbtag = (mcProd_0 == 0.0 || dataProd_0 == 0.0) ? 1.0 : dataProd_0/mcProd_0;
    
    //re-weighting for events with exactly 1 b-tag jet
  } else if(NminBjets==-1){
    /*
    std::cout<<"1 btag" <<std::endl;
    std::cout<<"dataSum: " << dataSum << std::endl;
    std::cout<<"mcSum:   " << mcSum << std::endl;
    */
    return wtbtag = (mcSum ==0.0 || dataSum ==0.0 ) ? 1.0 : dataSum/mcSum;
    
    // re-weighting for events with 1 or more b-tag jets
  } else if(NminBjets==1){
    /*
    std::cout<<">=1 btag"<<std::endl;
    std::cout<<"dataProd: "<<(1 - dataProd_0) <<std::endl;
    std::cout<<"mcProd:   "<<(1 - mcProd_0) <<std::endl;
    */
    return wtbtag = (dataProd_0 == 1.0 || mcProd_0 ==1.0) ? 1.0 : (1 - dataProd_0)/(1 - mcProd_0);
    // re-weighting for events with 2 or more b-tag jets
  } else if(NminBjets==2){
    /*
    std::cout<<">=2 btag"<<std::endl;
    std::cout<<"1 - dataProd_0 - dataSum: " << 1 - dataProd_0 - dataSum <<std::endl;
    std::cout<<"1 - mcProd_0 - mcSum:     " << 1 - mcProd_0 - mcSum <<std::endl;
    */
    return wtbtag = (dataProd_0 + dataSum == 1.0 || mcProd_0 + mcSum == 1.0) ? 1.0 : (1 - dataProd_0 - dataSum)/(1 - mcProd_0 - mcSum);
    
  } 
  
  else return 1;
   
  
  
  //return wtbtagErr = sqrt(2*(pow(err1+err2,2)))*wtbtag;
  
}


#endif
