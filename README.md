#####Table Of Content

* [Instructions](#instructions)
* [POG Recipes for Moriond 2017](#pog-recipes-for-moriond-2017)
* [How To Make PU Distribution for data](#how-to-make-pu-distribution-for-data)
* [To Do List](#to-do-list)
---
The package contains a code to produce ntuple for WW semileptonic final state.
It takes in input ntuples produced from miniAOD with the Bacon (https://github.com/ksung25/BaconProd )


## Instructions

	cmsrel CMSSW_8_0_11
	cd CMSSW_8_0_11/src
	cmsenv
	git clone git@github.com:ksung25/BaconAna.git
	mkdir WWAnalysis
	cd WWAnalysis
	git clone https://github.com/ram1123/WWAnalysisRun2.git;
	cd WWAnalysisRun2
	git checkout bacon_80x
	scramv1 b -j 8
	python python/produceWWNtuples.py -l el

* To see the various options availabe with **produceWWNtuples.py** do,

		python python/produceWWNtuples.py --help

## To Do List
- [ ] Apply latest ID, ISO, & Trigger Efficiency for electron and muons
- [ ] Apply ID for jets
- [ ] Apply PU
- [ ] Clean the code

## POG Recipes for Moriond 2017

* Pile-up reweighting xsec = 69.2mb

Ref: https://twiki.cern.ch/twiki/bin/viewauth/CMS/POGRecipesICHEP2016

## How To Make PU Distribution for data

	pileupCalc.py -i MyAnalysisJSON.txt --inputLumiJSON pileup_latest.txt  --calcMode true --minBiasXsec 69200 --maxPileupBin 50 --numPileupBins 50  MyDataPileupHistogram.root

where,
* MyAnalysisJSON.txt is the JSON file we are using.	
* pileup_latest.txt : this is input json file. It can be found at link: /afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions15/13TeV/PileUp/pileup_latest.txt

Reference: https://twiki.cern.ch/twiki/bin/viewauth/CMS/PileupJSONFileforData#Pileup_JSON_Files_For_Run_II



