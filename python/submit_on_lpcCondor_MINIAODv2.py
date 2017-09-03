#! /usr/bin/env python
import os
import glob
import math
from array import array
import sys
import time
import subprocess

currentDir = os.getcwd();
CMSSWDir =  currentDir+"/../";

#inputFolder = "/store/user/lnujj/WpWm_aQGC_Ntuples_Ram/FirstStepOutput/Feb142016";
inputFolder = "/store/user/lnujj/WpWm_aQGC_Ntuples_Ram/FirstStepOutput/BaconNtuples/";
outputFolder = currentDir+"/output/";
exePathName = currentDir+"/WWAnalysisRun2/produceWWNtuples"

lumi = 35900.0

dryRun = False;
doMC = False;
doData = True;

#category = ["el"];
category = ["el","mu"];
#category = ["EleMu"];

samples = [
    ( 0.9114,	"WplusToLNuWminusTo2JJJ_EWK_LO_SM_MJJ100PTJ10_TuneCUETP8M1_13TeV-madgraph-pythia8",	0),
    ( 0.9107,	"WplusTo2JWminusToLNuJJ_EWK_LO_SM_MJJ100PTJ10_TuneCUETP8M1_13TeV-madgraph-pythia8",	0),
    ( 542.00,	"TTToSemilepton_powheg_1",	0),
    ( 542.00,	"TTToSemilepton_powheg_2",	0),
    ( 542.00,	"TTToSemilepton_powheg_3",	0),
    ( 542.00,	"TTToSemilepton_powheg_4",	0),
    ( 542.00,	"TTToSemilepton_powheg_5",	0),
    ( 542.00,	"TTToSemilepton_powheg_6",	0),
    ( 0.01398,	"ZZZ_13TeV_amcatnlo_pythia8",	0),
    ( 0.1651,	"WWZ_13TeV_amcatnlo_pythia8",	0),
    ( 49.997,	"WWTo1L1Nu2Q_13TeV_amcatnloFXFX_madspin_pythia8",	0),
    ( 3.22,	"ZZTo2L2Q_13TeV_amcatnloFXFX_madspin_pythia8",	0),
    ( 5.595,	"WZTo2L2Q_13TeV_amcatnloFXFX_madspin_pythia8",	0),
    ( 10.71,	"WZTo1L1Nu2Q_13TeV_amcatnloFXFX_madspin_pythia8",	0),
    ( 5765.40,	"DYJetsToLL_M-50_amcatnlo_1",	0),
    ( 5765.40,	"DYJetsToLL_M-50_amcatnlo_2",	0),
    ( 5765.40,	"DYJetsToLL_M-50_amcatnlo_3",	0),
    ( 5765.40,	"DYJetsToLL_M-50_amcatnlo_4",	0),
    ( 5765.40,	"DYJetsToLL_M-50_amcatnlo_5",	0),
    ( 1.68,	"DYJetsToLL_M-50_HT-600to800",	0),
    #( 			"WJetsToLNu_HT_70To100_13TeV"
    ( 1506.4,		"WJetsToLNu_HT_100To200_13TeV", 		0),
    ( 1506.4,		"WJetsToLNu_HT_100To200_13TeV_ext1",		0),
    ( 1506.4,		"WJetsToLNu_HT_100To200_13TeV_ext2",		0),
    ( 435.237,		"WJetsToLNu_HT_200To400_13TeV",			0),
    ( 435.237,		"WJetsToLNu_HT_200To400_13TeV_ext1",		0),
    ( 435.237,		"WJetsToLNu_HT_200To400_13TeV_ext2",		0),
    ( 59.1811,		"WJetsToLNu_HT_400To600_13TeV",			0),
    ( 59.1811,		"WJetsToLNu_HT_400To600_13TeV_ext1",		0),
    ( 14.5805,		"WJetsToLNu_HT_600To800_13TeV",			0),
    ( 14.5805,		"WJetsToLNu_HT_600To800_13TeV_ext1",		0),
    ( 6.65621,		"WJetsToLNu_HT_800To1200_13TeV",		0),
    ( 6.65621,		"WJetsToLNu_HT_800To1200_13TeV_ext1",		0),
    ( 1.60809,		"WJetsToLNu_HT_1200To2500_13TeV",		0),
    ( 1.60809,		"WJetsToLNu_HT_1200To2500_13TeV_ext1",		0),
    ( 0.0389136,	"WJetsToLNu_HT_2500ToInf_13TeV",		0),
    ( 0.0389136,	"WJetsToLNu_HT_2500ToInf_13TeV_ext1",		0),
    ( 19.47,		"ST_tW_top_5f_NoFullyHadronicDecays_13TeV_powheg_TuneCUETP8M1",	0),
    ( 19.47,		"ST_tW_antitop_5f_NoFullyHadronicDecays_13TeV_powheg_TuneCUETP8M1",	0)
    ]

nameDataMu = [
    "SingleMuonRun2016B_03Feb2017_ver1_v1",
    "SingleMuonRun2016B_03Feb2017_ver2_v2",
    "SingleMuonRun2016C_03Feb2017_v1",
    "SingleMuonRun2016D_03Feb2017_v1",
    "SingleMuonRun2016E_03Feb2017_v1",
    "SingleMuonRun2016F_03Feb2017_v1",
    "SingleMuonRun2016G_03Feb2017_v1",
    "SingleMuonRun2016H_03Feb2017_ver2_v1",
    "SingleMuonRun2016H_03Feb2017_ver3_v1"
    ];

nameDataEl = [
"SingleElectron_Run2016B-03Feb2017_ver1-v1",
"SingleElectron_Run2016B-03Feb2017_ver2-v2",
"SingleElectron_Run2016C-03Feb2017-v1",
"SingleElectron_Run2016D-03Feb2017-v1",
"SingleElectron_Run2016E-03Feb2017-v1",
"SingleElectronRun2016F_03Feb2017_v1",
"SingleElectron_Run2016G-03Feb2017-v1",
"SingleElectronRun2016H_03Feb2017_ver2_v1",
"SingleElectronRun2016H_03Feb2017_ver3_v1"
];


inputlist="python/produceWWNtuples.py, ElectronTrigger_SF.root,  MuonTrigger_RunGH_23SepReReco_16p146fb.root, MuonID_RunBCDEF_23SepReReco_19p72fb.root, PileUpData2016_23Sep2016ReReco_69200ub.root, MuonID_RunGH_23SepReReco_16p146fb.root, egammaEffi_EGM2D_TightCutBasedIDSF.root, MuonIso_RunBCDEF_23SepReReco_19p72fb.root, egammaEffi_SF2D_GSF_tracking.root, MuonIso_RunGH_23SepReReco_16p146fb.root, puWeights_80x_37ifb.root, MuonTrigger_RunBCDEF_23SepReReco_19p72fb.root"

nameData = {"el": nameDataEl, "mu":nameDataMu};

#command = "python "+currentDir+"/WWAnalysisRun2/python/produceWWNtuples.py -i "+inputFolder+" $*";
command = "python produceWWNtuples.py -i "+inputFolder+" $*";

outScript = open("runstep2condor.sh","w");
outScript.write('#!/bin/bash');
outScript.write("\n"+'cd '+CMSSWDir);
outScript.write("\n"+'eval `scram runtime -sh`');
outScript.write("\n"+"cd -");
outScript.write("\n"+command);
outScript.write("\n");
outScript.close();
os.system("chmod 777 runstep2condor.sh");

outJDL = open("runstep2condor.jdl","w");
outJDL.write("Executable = runstep2condor.sh\n");
outJDL.write("Universe = vanilla\n");
outJDL.write("Requirements =FileSystemDomain==\"fnal.gov\" && Arch==\"X86_64\"");
outJDL.write("\n");
outJDL.write("Notification = ERROR\n");
outJDL.write("Should_Transfer_Files = YES\n");
outJDL.write("WhenToTransferOutput = ON_EXIT\n");
#outJDL.write("include : list-infiles.sh |\n");
outJDL.write("transfer_input_files = "+inputlist+"\n");
outJDL.write("x509userproxy = $ENV(X509_USER_PROXY)\n");

for a in range(len(category)):

    #MC
    if( doMC ):
        for i in range(len(samples)):
            outJDL.write("Output = "+str(samples[i][1])+"_"+category[a]+"_New.stdout\n");
            outJDL.write("Error = "+str(samples[i][1])+"_"+category[a]+"_New.stdout\n");
            outJDL.write("Arguments = -n "+str(samples[i][1])+" -o WWTree_"+str(samples[i][1])+"_"+category[a]+" -w "+str(samples[i][0])+" -lumi "+str(lumi)+" --ismc 1 -trig 1 -c lpc\n");
            outJDL.write("Queue\n");
    
    #data
    if( doData ):
        for i in range(len(nameData[category[a]])):
            outJDL.write("Output = "+(nameData[category[a]])[i]+"_New.stdout\n");
            outJDL.write("Error = "+(nameData[category[a]])[i]+"_New.stdout\n");
            outJDL.write("Arguments = -n "+(nameData[category[a]])[i]+" -o WWTree_"+(nameData[category[a]])[i]+"_"+category[a]+"_New"+" -w 1. -no 1. --ismc 0 -trig 1 -c lpc\n");
            outJDL.write("Queue\n");

outJDL.close();
print "===> Set Proxy Using:";
print "\tvoms-proxy-init --voms cms --valid 168:00";
print "\"condor_submit runstep2condor.jdl\" to submit";

