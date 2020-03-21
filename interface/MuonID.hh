#ifndef MUON_IDS_HH
#define MUON_IDS_HH

#include "BaconAna/DataFormats/interface/TMuon.hh"

bool passMuonLoose(const baconhep::TMuon *muon);
bool passMuonTight(const baconhep::TMuon *muon);

//https://twiki.cern.ch/twiki/bin/view/CMS/SWGuideMuonIdRun2 
bool passMuonLoose(const baconhep::TMuon *muon)
{
  if(!(muon->pogIDBits & baconhep::kPOGLooseMuon)) return false;
  
  // PF-isolation with Delta-beta correction
  double iso = muon->chHadIso + TMath::Max(muon->neuHadIso + muon->gammaIso - 0.5*(muon->puIso), double(0));
  if(iso >= 0.25*(muon->pt)) return false;
  
  return true;
}

//https://twiki.cern.ch/twiki/bin/view/CMS/SWGuideMuonIdRun2 
bool passMuonTight(const baconhep::TMuon *muon)
{
  if(!(muon->pogIDBits & baconhep::kPOGTightMuon)) return false;
  
  // PF-isolation with Delta-beta correction
  double iso = muon->chHadIso + TMath::Max(muon->neuHadIso + muon->gammaIso - 0.5*(muon->puIso), double(0));
  if(iso >= 0.15*(muon->pt)) return false;
  return true;
}

#endif
