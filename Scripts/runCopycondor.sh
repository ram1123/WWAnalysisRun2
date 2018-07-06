#!/bin/bash
cd /uscms_data/d3/rasharma/aQGC_analysis/SecondStep/CMSSW_8_0_26_patch1/src/WWAnalysis/WWAnalysisRun2/Scripts/../
eval `scram runtime -sh`
cd -
echo "sample Name : "$1
xrdfs root://cmseos.fnal.gov/ mkdir /store/user/lnujj/WpWm_aQGC_Ntuples_Ram/FirstStepOutput/BaconNtuples//$1/
xrdcp -r root://cmseos.fnal.gov//eos/uscms/store/user/aapyan2/Run2//$1/  root://cmseos.fnal.gov//store/user/lnujj/WpWm_aQGC_Ntuples_Ram/FirstStepOutput/BaconNtuples//$1/
