#ifndef JET_ID_HH
#define JET_ID_HH

#include <iostream>
#include <TF1.h>
#include <TH1F.h>
#include "BaconAna/DataFormats/interface/TJet.hh"
#include "BaconAna/DataFormats/interface/TAddJet.hh"

bool passAK4JetLoose(const baconhep::TJet *jet, int era);

bool passAK4JetLoose(const baconhep::TJet *jet, int era) {
  if (era!=2016 && era!=2017) {
    std::cout << "Invalid Era!" << std::endl;
    return false;
  }
  else if (era==2016) {
    //2016
    //https://twiki.cern.ch/twiki/bin/view/CMS/JetID13TeVRun2016
    if(fabs(jet->eta)<= 2.7){
      if(jet->neuHadFrac >= 0.99) return false;
      if(jet->neuEmFrac  >= 0.99) return false;
      if(jet->nParticles <= 1)    return false;
    }
    if(fabs(jet->eta)<= 2.4) {
      if(jet->chHadFrac == 0)     return false;
      if(jet->nCharged  == 0)     return false;
      if(jet->chEmFrac  >= 0.99)  return false;
    }
    if(fabs(jet->eta) > 2.7 && fabs(jet->eta) <= 3.0) {
      if(jet->neuEmFrac <= 0.01)  return false;
      if(jet->neuHadFrac >= 0.98) return false;
      if(jet->nNeutrals <= 2)     return false;
    }
    if(fabs(jet->eta) > 3.0) {
      if(jet->neuEmFrac >= 0.90)  return false;
      if(jet->nNeutrals <= 10)    return false;
    }
  }
  else if (era==2017) {
    //2017
    //https://twiki.cern.ch/twiki/bin/view/CMS/JetID13TeVRun2017
    //(Tight)
    if(fabs(jet->eta)<= 2.7){
      if(jet->neuHadFrac >= 0.90) return false;
      if(jet->neuEmFrac  >= 0.90) return false;
      if(jet->nParticles <= 1)    return false;
    }
    if(fabs(jet->eta)<= 2.4) {
      if(jet->chHadFrac == 0)     return false;
      if(jet->nCharged  == 0)     return false;
    }
    if(fabs(jet->eta) > 2.7 && fabs(jet->eta) <= 3.0) {
      if(jet->neuEmFrac <= 0.02)  return false;
      if(jet->neuEmFrac >= 0.99)  return false;
      if(jet->nNeutrals <= 2)     return false;
    }
    if(fabs(jet->eta) > 3.0) {
      if(jet->neuHadFrac <= 0.02) return false;
      if(jet->neuEmFrac  >= 0.90) return false;
      if(jet->nNeutrals  <= 10)   return false;
    }
  }
  return true;
}

float getPUPPIweight(TF1* puppisd_corrGEN, TF1* puppisd_corrRECO_cen, TF1* puppisd_corrRECO_for, float puppipt, float puppieta ){
  float genCorr  = 1.;
  float recoCorr = 1.;
  float totalWeight = 1.;
        
  genCorr =  puppisd_corrGEN->Eval( puppipt );
  if( fabs(puppieta)  <= 1.3 ){
    recoCorr = puppisd_corrRECO_cen->Eval( puppipt );
  }
  else{
    recoCorr = puppisd_corrRECO_for->Eval( puppipt );
  }
  
  totalWeight = genCorr * recoCorr;


  return totalWeight;
}

#endif
