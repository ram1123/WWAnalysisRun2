#ifndef ELECTRON_IDS_HH
#define ELECTRON_IDS_HH

#include <iostream>
#include "BaconAna/DataFormats/interface/TElectron.hh"

bool passEleLoose(const baconhep::TElectron *electron, float iso, int era);
bool passEleTight(const baconhep::TElectron *electron, float iso, int era);
float eleEffArea(const float eta, int era);

bool passEleLoose(const baconhep::TElectron *electron, float iso, int era) {
  if(electron->isConv) return false;

  if (era!=2016 && era!=2017) {
    std::cout << "Invalid Era!" << std::endl;
    return false;
  }
  else if (era==2016) {
  //2016
  //https://twiki.cern.ch/twiki/bin/view/CMS/CutBasedElectronIdentificationRun2#Working_points_for_2016_data_for
    if(fabs(electron->scEta)<1.479) {
      if(iso                          >= 0.0994*(electron->pt)) 
	return false;
      if(electron->sieie              >= 0.01100)
	return false;
      if(fabs(electron->dEtaInSeed)   >= 0.00477)
	return false;
      if(fabs(electron->dPhiIn)       >= 0.22200)
	return false;
      if(electron->hovere             >= 0.29800)
	return false;
      if(fabs(1.0 - electron->eoverp) >= 0.24100*(electron->ecalEnergy))
	return false;
      if(electron->nMissingHits       >  1)
	return false;
    } else {
      if(iso                          >= 0.107*(electron->pt))
	return false;
      if(electron->sieie              >= 0.03140)
	return false;
      if(fabs(electron->dEtaInSeed)   >= 0.00868)
	return false;
      if(fabs(electron->dPhiIn)       >= 0.21300)
	return false;
      if(electron->hovere             >= 0.10100)
	return false;
      if(fabs(1.0 - electron->eoverp) >= 0.14000*(electron->ecalEnergy))
	return false;
      if(electron->nMissingHits       >  1)
	return false;
    }
  }
  else if (era==2017) {
    //2017
    //https://twiki.cern.ch/twiki/bin/view/CMS/CutBasedElectronIdentificationRun2#Offline_selection_criteria_for_V
    float rho = eleEffArea(electron->eta, era);
    
    if(fabs(electron->scEta)<1.479) {
      if(iso >= 0.112 + 0.506/(electron->pt))
	return false;
      if(electron->sieie              >= 0.0112)
	return false;
      if(fabs(electron->dEtaInSeed)   >= 0.00377)
	return false;
      if(fabs(electron->dPhiIn)       >= 0.0884)
	return false;
      if(electron->hovere             >= 0.05+(1.16+0.0324*rho)/(electron->ecalEnergy))
	return false;
      if(fabs(1.0 - electron->eoverp) >= 0.193*(electron->ecalEnergy)) 
	return false;
      if(electron->nMissingHits       >  1)
	return false;
    } else {
      if(iso >= 0.108 + 0.963/(electron->pt))
	return false;
      if(electron->sieie              >= 0.0425)
	return false;
      if(fabs(electron->dEtaInSeed)   >= 0.00674)
	return false;
      if(fabs(electron->dPhiIn)       >= 0.169)
	return false;
      if(electron->hovere             >= 0.0441+(2.54+0.183*rho)/(electron->ecalEnergy))
	return false;
      if(fabs(1.0 - electron->eoverp) >= 0.111*(electron->ecalEnergy)) 
	return false;
      if(electron->nMissingHits       >  1)
	return false;
    }
  }

  return true;
}

bool passEleTight(const baconhep::TElectron *electron, float iso, int era) {

  if(electron->isConv) return false;

  if (era!=2016 && era!=2017) {
    std::cout << "Invalid Era!" << std::endl;
    return false;
  }
  else if (era==2016) {
    //2016
    //https://twiki.cern.ch/twiki/bin/view/CMS/CutBasedElectronIdentificationRun2#Working_points_for_2016_data_for
    if(fabs(electron->scEta)<1.479) {
      if(iso                        >= 0.0588*(electron->pt))
	return false;
      if(electron->sieie              >= 0.00998)
	return false;
      if(fabs(electron->dEtaInSeed)   >= 0.00308)
	return false;
      if(fabs(electron->dPhiIn)       >= 0.08160)
	return false;
      if(electron->hovere             >= 0.04140)
	return false;
      if(fabs(1.0 - electron->eoverp) >= 0.12900*(electron->ecalEnergy))
	return false;
      if(electron->nMissingHits       >  1)
	return false;
    } else {
      if(iso                          >= 0.0571*(electron->pt)) 
	return false;
      if(electron->sieie              >= 0.02920)
	return false;
      if(fabs(electron->dEtaInSeed)   >= 0.00605)
	return false;
      if(fabs(electron->dPhiIn)       >= 0.03940)
	return false;
      if(electron->hovere             >= 0.06410)
	return false;
      if(fabs(1.0 - electron->eoverp) >= 0.12900*(electron->ecalEnergy)) 
	return false;
      if(electron->nMissingHits       >  1)                              
	return false;
    }
  }
  else if (era==2017) {
    //2017
    //https://twiki.cern.ch/twiki/bin/view/CMS/CutBasedElectronIdentificationRun2#Offline_selection_criteria_for_V
    float rho = eleEffArea(electron->eta,era);
    
    if(fabs(electron->scEta)<1.479) {
      if(iso                          >= 0.0287 + 0.506/(electron->pt))
	return false;
      if(electron->sieie              >= 0.0104)
	return false;
      if(fabs(electron->dEtaInSeed)   >= 0.00255)
	return false;
      if(fabs(electron->dPhiIn)       >= 0.022)
	return false;
      if(electron->hovere             >= 0.026+(1.15+0.0324*rho)/(electron->ecalEnergy))
	return false;
      if(fabs(1.0 - electron->eoverp) >= 0.159*(electron->ecalEnergy)) 
	return false;
      if(electron->nMissingHits       >  1)
	return false;
    } else {
      if(iso                          >= 0.0445 + 0.963/(electron->pt))
	return false;
      if(electron->sieie              >= 0.0353)
	return false;
      if(fabs(electron->dEtaInSeed)   >= 0.00501)
	return false;
      if(fabs(electron->dPhiIn)       >= 0.0236)
	return false;
      if(electron->hovere             >= 0.0188+(2.06+0.183*rho)/(electron->ecalEnergy))
	return false;
      if(fabs(1.0 - electron->eoverp) >= 0.0197*(electron->ecalEnergy)) 
	return false;
      if(electron->nMissingHits       >  1)
	return false;
    }
  }
  return true;
}

float eleEffArea(const float eta, int era) {

  if (era==2016) {
    //2016
    // https://github.com/ikrav/cmssw/blob/egm_id_80X_v1/RecoEgamma/ElectronIdentification/data/Summer16/effAreaElectrons_cone03_pfNeuHadronsAndPhotons_80X.txt
    if     (fabs(eta) >= 0.0000 && fabs(eta) < 1.0000) { return 0.1703; }
    else if(fabs(eta) >= 1.0000 && fabs(eta) < 1.4790) { return 0.1715; }
    else if(fabs(eta) >= 1.4790 && fabs(eta) < 2.0000) { return 0.1213; }
    else if(fabs(eta) >= 2.0000 && fabs(eta) < 2.2000) { return 0.1230; }
    else if(fabs(eta) >= 2.2000 && fabs(eta) < 2.3000) { return 0.1635; }
    else if(fabs(eta) >= 2.3000 && fabs(eta) < 2.4000) { return 0.1937; }
    else                                               { return 0.2393; }
  }
  else if (era==2017) {
    //2017
    //https://github.com/cms-sw/cmssw/blob/CMSSW_10_2_X/RecoEgamma/ElectronIdentification/data/Fall17/effAreaElectrons_cone03_pfNeuHadronsAndPhotons_94X.txt
    if     (fabs(eta) >= 0.0000 && fabs(eta) < 1.0000) { return 0.1440; }
    else if(fabs(eta) >= 1.0000 && fabs(eta) < 1.4790) { return 0.1562; }
    else if(fabs(eta) >= 1.4790 && fabs(eta) < 2.0000) { return 0.1032; }
    else if(fabs(eta) >= 2.0000 && fabs(eta) < 2.2000) { return 0.0859; }
    else if(fabs(eta) >= 2.2000 && fabs(eta) < 2.3000) { return 0.1116; }
    else if(fabs(eta) >= 2.3000 && fabs(eta) < 2.4000) { return 0.1321; }
    else                                               { return 0.1654; }  
  }

  std::cout << "Invalid Era!" << std::endl;
  return 0.0;

}
  
#endif

