#! /usr/bin/env python
import os
import glob
import math
from array import array
import sys
import time
import subprocess

currentDir = os.getcwd();
CMSSWDir = currentDir+"/../";

inputFolder = "/store/cmst3/group/monojet/production/12/";
#inputFolder = "/eos/uscms/store/user/lnujj/WpWm_aQGC_Ntuples_Ram/FirstStepOutput/Jan102016";
#outputFolder = currentDir+"/output/";
outputFolder = currentDir+"/WWAnalysisRun2/output_test/";
exeName = currentDir+"/WWAnalysisRun2/produceWWNtuples"

dryRun = False;
doMC = True;
doData = False;

category = ["mu","el"];
#category = ["el"];
#category = ["mu"];



samples = [
    ( 49.997,		"WWToLNuQQ_13TeV_powheg",			0),
    ( 49.997,		"WWToLNuQQ_13TeV_powheg_ext",			0),
    ( 61526.7,		"WJetsToLNu_13TeV",				0),
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
    #( 6.65621,		"WJetsToLNu_HT_800To1200_13TeV_ext1",	
    ( 1.60809,		"WJetsToLNu_HT_1200To2500_13TeV",		0),
    ( 1.60809,		"WJetsToLNu_HT_1200To2500_13TeV_ext1",		0),
    ( 0.0389136,	"WJetsToLNu_HT_2500ToInf_13TeV",		0),
    ( 0.0389136,	"WJetsToLNu_HT_2500ToInf_13TeV_ext1",		0),
    ( 49.997,		"WWTo1L1Nu2Q_13TeV_amcatnloFXFX_madspin_pythia8",0),
    ( 49.997,		"WWToLNuQQ_13TeV_powheg",			0),
    ( 49.997,		"WWToLNuQQ_13TeV_powheg_ext",			0),
    ( 10.71,		"WZTo1L1Nu2Q_13TeV_amcatnloFXFX_madspin_pythia8",0),
    ( 3.22,		"ZZTo2L2Q_13TeV_amcatnloFXFX_madspin_pythia8",	0),
    ( 831.76,		"TTToSemilepton_powheg",			0),
    ( 11.36,		"ST_s_channel_4f_leptonDecays_13TeV_amcatnlo_pythia8_TuneCUETP8M1",	0),
    ( 80.95,		"ST_t_channel_antitop_4f_inclusiveDecays_13TeV_powhegV2_madspin_pythia8_TuneCUETP8M1",	0),
    ( 136.02,		"ST_t_channel_top_4f_inclusiveDecays_13TeV_powhegV2_madspin_pythia8_TuneCUETP8M1",	0),
    ( 19.46,		"ST_tW_antitop_5f_NoFullyHadronicDecays_13TeV_powheg_TuneCUETP8M1",	0),
    ( 19.46,		"ST_tW_top_5f_NoFullyHadronicDecays_13TeV_powheg_TuneCUETP8M1",		0)
    ]

nameDataMu = [
"data_mu_2016_runB_v3_1",
"data_mu_2016_runB_v3_10",
"data_mu_2016_runB_v3_11",
"data_mu_2016_runB_v3_12",
"data_mu_2016_runB_v3_13",
"data_mu_2016_runB_v3_14",
"data_mu_2016_runB_v3_15",
"data_mu_2016_runB_v3_16",
"data_mu_2016_runB_v3_17",
"data_mu_2016_runB_v3_18",
"data_mu_2016_runB_v3_19",
"data_mu_2016_runB_v3_2",
"data_mu_2016_runB_v3_3",
"data_mu_2016_runB_v3_4",
"data_mu_2016_runB_v3_5",
"data_mu_2016_runB_v3_6",
"data_mu_2016_runB_v3_7",
"data_mu_2016_runB_v3_8",
"data_mu_2016_runB_v3_9",
"data_mu_2016_runC_v1",
"data_mu_2016_runC_v1_1",
"data_mu_2016_runC_v1_2",
"data_mu_2016_runC_v1_3",
"data_mu_2016_runC_v1_4",
"data_mu_2016_runC_v1_5",
"data_mu_2016_runC_v1_6",
"data_mu_2016_runC_v1_7",
"data_mu_2016_runD_v1_1",
"data_mu_2016_runD_v1_10",
"data_mu_2016_runD_v1_11",
"data_mu_2016_runD_v1_12",
"data_mu_2016_runD_v1_2",
"data_mu_2016_runD_v1_3",
"data_mu_2016_runD_v1_4",
"data_mu_2016_runD_v1_5",
"data_mu_2016_runD_v1_6",
"data_mu_2016_runD_v1_7",
"data_mu_2016_runD_v1_8",
"data_mu_2016_runD_v1_9",
"data_mu_2016_runE_v1_1",
"data_mu_2016_runE_v1_10",
"data_mu_2016_runE_v1_2",
"data_mu_2016_runE_v1_3",
"data_mu_2016_runE_v1_4",
"data_mu_2016_runE_v1_5",
"data_mu_2016_runE_v1_6",
"data_mu_2016_runE_v1_7",
"data_mu_2016_runE_v1_8",
"data_mu_2016_runE_v1_9",
"data_mu_2016_runF_v1_1",
"data_mu_2016_runF_v1_2",
"data_mu_2016_runF_v1_3",
"data_mu_2016_runF_v1_4",
"data_mu_2016_runF_v1_5",
"data_mu_2016_runF_v1_6",
"data_mu_2016_runF_v1_7",
"data_mu_2016_runF_v1_8",
"data_mu_2016_runG_v1_1",
"data_mu_2016_runG_v1_10",
"data_mu_2016_runG_v1_11",
"data_mu_2016_runG_v1_12",
"data_mu_2016_runG_v1_13",
"data_mu_2016_runG_v1_14",
"data_mu_2016_runG_v1_15",
"data_mu_2016_runG_v1_16",
"data_mu_2016_runG_v1_17",
"data_mu_2016_runG_v1_2",
"data_mu_2016_runG_v1_3",
"data_mu_2016_runG_v1_4",
"data_mu_2016_runG_v1_5",
"data_mu_2016_runG_v1_6",
"data_mu_2016_runG_v1_7",
"data_mu_2016_runG_v1_8",
"data_mu_2016_runG_v1_9",
"data_mu_2016_runH_v3"
];

nameDataEl = [
"data_el_2016_runC_v1_1",
"data_el_2016_runC_v1_10",
"data_el_2016_runC_v1_11",
"data_el_2016_runC_v1_12",
"data_el_2016_runC_v1_13",
"data_el_2016_runC_v1_14",
"data_el_2016_runC_v1_2",
"data_el_2016_runC_v1_3",
"data_el_2016_runC_v1_4",
"data_el_2016_runC_v1_5",
"data_el_2016_runC_v1_6",
"data_el_2016_runC_v1_7",
"data_el_2016_runC_v1_8",
"data_el_2016_runC_v1_9",
"data_el_2016_runD_v1_1",
"data_el_2016_runD_v1_10",
"data_el_2016_runD_v1_11",
"data_el_2016_runD_v1_12",
"data_el_2016_runD_v1_2",
"data_el_2016_runD_v1_3",
"data_el_2016_runD_v1_4",
"data_el_2016_runD_v1_5",
"data_el_2016_runD_v1_6",
"data_el_2016_runD_v1_7",
"data_el_2016_runD_v1_8",
"data_el_2016_runD_v1_9",
"data_el_2016_runF_v1_1",
"data_el_2016_runF_v1_10",
"data_el_2016_runF_v1_11",
"data_el_2016_runF_v1_12",
"data_el_2016_runF_v1_13",
"data_el_2016_runF_v1_14",
"data_el_2016_runF_v1_15",
"data_el_2016_runF_v1_16",
"data_el_2016_runF_v1_17",
"data_el_2016_runF_v1_2",
"data_el_2016_runF_v1_3",
"data_el_2016_runF_v1_4",
"data_el_2016_runF_v1_5",
"data_el_2016_runF_v1_6",
"data_el_2016_runF_v1_7",
"data_el_2016_runF_v1_8",
"data_el_2016_runF_v1_9",
"data_el_2016_runG_v1_1",
"data_el_2016_runG_v1_10",
"data_el_2016_runG_v1_11",
"data_el_2016_runG_v1_12",
"data_el_2016_runG_v1_13",
"data_el_2016_runG_v1_14",
"data_el_2016_runG_v1_15",
"data_el_2016_runG_v1_16",
"data_el_2016_runG_v1_17",
"data_el_2016_runG_v1_18",
"data_el_2016_runG_v1_19",
"data_el_2016_runG_v1_2",
"data_el_2016_runG_v1_3",
"data_el_2016_runG_v1_4",
"data_el_2016_runG_v1_5",
"data_el_2016_runG_v1_6",
"data_el_2016_runG_v1_7",
"data_el_2016_runG_v1_8",
"data_el_2016_runG_v1_9",
"data_el_2016_runH_v2_1",
"data_el_2016_runH_v2_10",
"data_el_2016_runH_v2_11",
"data_el_2016_runH_v2_12",
"data_el_2016_runH_v2_13",
"data_el_2016_runH_v2_14",
"data_el_2016_runH_v2_15",
"data_el_2016_runH_v2_16",
"data_el_2016_runH_v2_17",
"data_el_2016_runH_v2_18",
"data_el_2016_runH_v2_19",
"data_el_2016_runH_v2_2",
"data_el_2016_runH_v2_3",
"data_el_2016_runH_v2_4",
"data_el_2016_runH_v2_5",
"data_el_2016_runH_v2_6",
"data_el_2016_runH_v2_7",
"data_el_2016_runH_v2_8",
"data_el_2016_runH_v2_9",
"data_el_2016_runH_v3"
];


nameData = {"el": nameDataEl, "mu":nameDataMu};



for a in range(len(category)):
    
    #MC
    if( doMC ):
        for i in range(len(samples)):
            fn = "WWAnalysisRun2/Job/job_"+samples[i][1]+"_"+category[a];
            outScript = open(fn+".sh","w");
            command = "python "+currentDir+"/WWAnalysisRun2/python/produceWWNtuples.py -i "+inputFolder+" -n "+str(samples[i][1])+" -o WWTree_"+str(samples[i][1])+" -l "+category[a]+" -w "+str(samples[i][0])+" -mass "+str(samples[i][2])+" --ismc 1 -trig 0";
            print command;
            outScript.write('#!/bin/bash');
            outScript.write("\n"+'cd '+CMSSWDir);
            outScript.write("\n"+'eval `scram runtime -sh`');
            #outScript.write("\n"+'cd '+currentDir);
            #outScript.write("\n"+'source scripts/setup.sh');
            outScript.write("\n"+"ls");
            outScript.write("\n"+"cd WWAnalysis/WWAnalysisRun2");
            outScript.write("\n"+"cp "+currentDir+"/WWAnalysisRun2/*.txt ./");
            outScript.write("\n"+"cp "+currentDir+"/WWAnalysisRun2/*.root ./");
            outScript.write("\n"+command);
            outScript.write("\n"+"mv WWTree_"+(samples[i][1])+".root "+outputFolder+"/output_"+category[a]);
            outScript.write("\n");
            outScript.close();
            os.system("chmod 777 "+currentDir+"/"+fn+".sh");
            command2 = "bsub -q 1nd -cwd "+currentDir+" "+currentDir+"/"+fn+".sh";
            print command2
            if( dryRun != True ):
                os.system(command2);
                time.sleep(3)
    
    #data
    if( doData ):
        for i in range(len(nameData[category[a]])):
            fn = "Job/job_"+(nameData[category[a]])[i]+"_"+category[a];
            outScript = open(fn+".sh","w");
            command = "python "+currentDir+"/WWAnalysisRun2/python/produceWWNtuples.py --exe "+exeName+" -i "+inputFolder+" -n "+(nameData[category[a]])[i]+" -o WWTree_"+(nameData[category[a]])[i]+" -l "+category[a]+" -w 1. -no 1. -mass 0 --ismc 0 -trig 1";
            print command;
            outScript.write('#!/bin/bash');
            outScript.write("\n"+'cd '+CMSSWDir);
            outScript.write("\n"+'eval `scram runtime -sh`');
            #outScript.write("\n"+'cd '+currentDir);
            #outScript.write("\n"+'source scripts/setup.sh');
            outScript.write("\n"+"cd -");
            outScript.write("\n"+"cp "+currentDir+"/*.txt ./");
            outScript.write("\n"+"cp "+currentDir+"/*.root ./");
            outScript.write("\n"+command);
            outScript.write("\n"+"mv WWTree_"+(nameData[category[a]])[i]+".root "+outputFolder+"/output_"+category[a]);
            outScript.write("\n");
            outScript.close();
            os.system("chmod 777 "+currentDir+"/"+fn+".sh");
            command2 = "bsub -q 1nd -cwd "+currentDir+" "+currentDir+"/"+fn+".sh";
            print command2
            if( dryRun != True ):
                os.system(command2);
                time.sleep(3)
