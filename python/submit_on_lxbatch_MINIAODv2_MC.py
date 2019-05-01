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

lumi = 35900.0

inputFolder = "/store/cmst3/group/monojet/production/12a/";	# Path of data files
outputFolder = currentDir+"/WWAnalysisRun2/output/";
exeName = currentDir+"/WWAnalysisRun2/produceWWNtuples"

MCs = 1 # 1 for signal and 2 for bkg
dryRun = False;
doMC = True;
doData = False;	# For data run another script

category = ["EleMu_new"];
#category = ["el"];
#category = ["mu"];



if MCs==1:
	inputFolder = "/store/user/arapyan/Run2/";	# Path of signal sample
	samples = [
	( 0.9114,	"WplusToLNuWminusTo2JJJ_EWK_LO_SM_MJJ100PTJ10_TuneCUETP8M1_13TeV-madgraph-pythia8",	0),
	( 0.9107,	"WplusTo2JWminusToLNuJJ_EWK_LO_SM_MJJ100PTJ10_TuneCUETP8M1_13TeV-madgraph-pythia8",	0),
	( 0.0879,	"WplusToLNuWplusTo2JJJ_EWK_LO_SM_MJJ100PTJ10_TuneCUETP8M1_13TeV-madgraph-pythia8",	0),
	( 0.0326,	"WminusToLNuWminusTo2JJJ_EWK_LO_SM_MJJ100PTJ10_TuneCUETP8M1_13TeV-madgraph-pythia8",	0),
	( 0.1825,	"WplusToLNuZTo2JJJ_EWK_LO_SM_MJJ100PTJ10_TuneCUETP8M1_13TeV-madgraph-pythia8",		0),
	( 0.0540,	"WplusTo2JZTo2LJJ_EWK_LO_SM_MJJ100PTJ10_TuneCUETP8M1_13TeV-madgraph-pythia8",		0),
	( 0.1000,	"WminusToLNuZTo2JJJ_EWK_LO_SM_MJJ100PTJ10_TuneCUETP8M1_13TeV-madgraph-pythia8",		0),
	( 0.0298,	"WminusTo2JZTo2LJJ_EWK_LO_SM_MJJ100PTJ10_TuneCUETP8M1_13TeV-madgraph-pythia8",		0),
	( 0.0159,	"ZTo2LZTo2JJJ_EWK_LO_SM_MJJ100PTJ10_TuneCUETP8M1_13TeV-madgraph-pythia8",		0),
	
	( 5.546,	"WplusTo2JWminusToLNuJJ_QCD_LO_SM_MJJ100PTJ10_TuneCUETP8M1_13TeV-madgraph-pythia8",	0),
	( 5.558,	"WplusToLNuWminusTo2JJJ_QCD_LO_SM_MJJ100PTJ10_TuneCUETP8M1_13TeV-madgraph-pythia8",	0),
	( 0.086,	"WplusToLNuWplusTo2JJJ_QCD_LO_SM_MJJ100PTJ10_TuneCUETP8M1_13TeV-madgraph-pythia8",	0),
	( 0.038,	"WminusToLNuWminusTo2JJJ_QCD_LO_SM_MJJ100PTJ10_TuneCUETP8M1_13TeV-madgraph-pythia8",	0),
	( 2.159,	"WplusToLNuZTo2JJJ_QCD_LO_SM_MJJ100PTJ10_TuneCUETP8M1_13TeV-madgraph-pythia8",		0),
	( 0.640,	"WplusTo2JZTo2LJJ_QCD_LO_SM_MJJ100PTJ10_TuneCUETP8M1_13TeV-madgraph-pythia8",		0),
	( 1.302,	"WminusToLNuZTo2JJJ_QCD_LO_SM_MJJ100PTJ10_TuneCUETP8M1_13TeV-madgraph-pythia8",		0),
	( 0.387,	"WminusTo2JZTo2LJJ_QCD_LO_SM_MJJ100PTJ10_TuneCUETP8M1_13TeV-madgraph-pythia8",		0),
	( 0.376,	"ZTo2LZTo2JJJ_QCD_LO_SM_MJJ100PTJ10_TuneCUETP8M1_13TeV-madgraph-pythia8",		0),
	
	# Cross-Section To be added
	( 17.94,	"WplusToLNuWminusTo2JJJ_EWK_LO_aQGC_MJJ100PTJ10_TuneCUETP8M1_13TeV-madgraph-pythia8",	0),
	( 17.92,	"WplusTo2JWminusToLNuJJ_EWK_LO_aQGC_MJJ100PTJ10_TuneCUETP8M1_13TeV-madgraph-pythia8",	0),
	( 3.451,	"WplusToLNuWplusTo2JJJ_EWK_LO_aQGC_MJJ100PTJ10_TuneCUETP8M1_13TeV-madgraph-pythia8",	0),
	( 0.5067,	"WminusToLNuWminusTo2JJJ_EWK_LO_aQGC_MJJ100PTJ10_TuneCUETP8M1_13TeV-madgraph-pythia8",	0),
	( 1.895,	"WplusToLNuZTo2JJJ_EWK_LO_aQGC_MJJ100PTJ10_TuneCUETP8M1_13TeV-madgraph-pythia8",	0),
	( 0.5686,	"WplusTo2JZTo2LJJ_EWK_LO_aQGC_MJJ100PTJ10_TuneCUETP8M1_13TeV-madgraph-pythia8",		0),
	( 0.7414,	"WminusToLNuZTo2JJJ_EWK_LO_aQGC_MJJ100PTJ10_TuneCUETP8M1_13TeV-madgraph-pythia8",	0),
	( 0.2223,	"WminusTo2JZTo2LJJ_EWK_LO_aQGC_MJJ100PTJ10_TuneCUETP8M1_13TeV-madgraph-pythia8",	0),
	( 3.361,	"ZTo2LZTo2JJJ_EWK_LO_aQGC_MJJ100PTJ10_TuneCUETP8M1_13TeV-madgraph-pythia8",		0)
	]

if MCs==2:
	inputFolder = "/store/cmst3/group/monojet/production/12/";	# Path of MC files
	samples = [
	    #( , "QCD_HT50to100_13TeV",		0),
	    ( 208.98,	"DYJetsToLL_M-50_HT-70to100",	0),
	    ( 181.30,	"DYJetsToLL_M-50_HT-100to200",	0),
	    ( 181.30,	"DYJetsToLL_M-50_HT-100to200_ext1",	0),
	    ( 50.42,	"DYJetsToLL_M-50_HT-200to400",	0),
	    ( 50.42,	"DYJetsToLL_M-50_HT-200to400_ext1",	0),
	    ( 6.98,	"DYJetsToLL_M-50_HT-400to600",	0),
	    ( 6.98,	"DYJetsToLL_M-50_HT-400to600_ext1",	0),
	    ( 1.68,	"DYJetsToLL_M-50_HT-600to800",	0),
	    ( 0.78,	"DYJetsToLL_M-50_HT-800to1200",	0),
	    ( 0.19,	"DYJetsToLL_M-50_HT-1200to2500",	0),
	    ( 0.0044,	"DYJetsToLL_M-50_HT-2500toInf",	0),
	    ( 27990000,	"QCD_HT100to200_13TeV",		0),
	    ( 1712000,	"QCD_HT200to300_13TeV",		0),
	    ( 1712000,	"QCD_HT200to300_13TeV_ext",	0),
	    ( 347700,	"QCD_HT300to500_13TeV",		0),
	    ( 347700,	"QCD_HT300to500_13TeV_ext",	0),
	    ( 32100,	"QCD_HT500to700_13TeV",		0),
	    ( 32100,	"QCD_HT500to700_13TeV_ext",	0),
	    ( 6831,	"QCD_HT700to1000_13TeV",	0),
	    ( 6831,	"QCD_HT700to1000_13TeV_ext",	0),
	    ( 1207,	"QCD_HT1000to1500_13TeV",	0),
	    ( 1207,	"QCD_HT1000to1500_13TeV_ext",	0),
	    ( 119.9,	"QCD_HT1500to2000_13TeV",	0),
	    ( 119.9,	"QCD_HT1500to2000_13TeV_ext",	0),
	    ( 25.24,	"QCD_HT2000toInf_13TeV",	0),
	    ( 25.24,	"QCD_HT2000toInf_13TeV_ext",	0),
	    ( 49.997,		"WWToLNuQQ_13TeV_powheg",			0),
	    ( 49.997,		"WWToLNuQQ_13TeV_powheg_ext",			0),
	    ( 61526.7,		"WJetsToLNu_13TeV",				0),
	    ##( 			"WJetsToLNu_HT_70To100_13TeV"
	    #( 1506.4,		"WJetsToLNu_HT_100To200_13TeV", 		0),
	    #( 1506.4,		"WJetsToLNu_HT_100To200_13TeV_ext1",		0),
	    #( 1506.4,		"WJetsToLNu_HT_100To200_13TeV_ext2",		0),
	    #( 435.237,		"WJetsToLNu_HT_200To400_13TeV",			0),
	    #( 435.237,		"WJetsToLNu_HT_200To400_13TeV_ext1",		0),
	    #( 435.237,		"WJetsToLNu_HT_200To400_13TeV_ext2",		0),
	    #( 59.1811,		"WJetsToLNu_HT_400To600_13TeV",			0),
	    #( 59.1811,		"WJetsToLNu_HT_400To600_13TeV_ext1",		0),
	    #( 14.5805,		"WJetsToLNu_HT_600To800_13TeV",			0),
	    #( 14.5805,		"WJetsToLNu_HT_600To800_13TeV_ext1",		0),
	    #( 6.65621,		"WJetsToLNu_HT_800To1200_13TeV",		0),
	    ##( 6.65621,		"WJetsToLNu_HT_800To1200_13TeV_ext1",	
	    #( 1.60809,		"WJetsToLNu_HT_1200To2500_13TeV",		0),
	    #( 1.60809,		"WJetsToLNu_HT_1200To2500_13TeV_ext1",		0),
	    #( 0.0389136,	"WJetsToLNu_HT_2500ToInf_13TeV",		0),
	    #( 0.0389136,	"WJetsToLNu_HT_2500ToInf_13TeV_ext1",		0),
	    #( 49.997,		"WWTo1L1Nu2Q_13TeV_amcatnloFXFX_madspin_pythia8",0),
	    #( 10.71,		"WZTo1L1Nu2Q_13TeV_amcatnloFXFX_madspin_pythia8",0),
	    #( 3.22,		"ZZTo2L2Q_13TeV_amcatnloFXFX_madspin_pythia8",	0),
	    #( 542.00,		"TTToSemilepton_powheg",			0),
	    ( 11.36,		"ST_s_channel_4f_leptonDecays_13TeV_amcatnlo_pythia8_TuneCUETP8M1",	0),
	    ( 80.95,		"ST_t_channel_antitop_4f_inclusiveDecays_13TeV_powhegV2_madspin_pythia8_TuneCUETP8M1",	0),
	    ( 136.02,		"ST_t_channel_top_4f_inclusiveDecays_TuneCUETP8M2T4_13TeV_powhegV2_madspin",	0),
	    ( 19.46,		"ST_tW_antitop_5f_NoFullyHadronicDecays_13TeV_powheg_TuneCUETP8M1",	0),
	    ( 19.46,		"ST_tW_top_5f_NoFullyHadronicDecays_13TeV_powheg_TuneCUETP8M1",		0)
	    ]

nameDataMu = [
#"SingleMuonRun2016B_23Sep2016_v1",
#"SingleMuonRun2016B_23Sep2016_v3",
#"SingleMuonRun2016C_23Sep2016_v1",
#"SingleMuonRun2016D_23Sep2016_v1",
#"SingleMuonRun2016E_23Sep2016_v1",
#"SingleMuonRun2016F_23Sep2016_v1",
#"SingleMuonRun2016G_23Sep2016_v1",
#"SingleMuonRun2016H_PromptReco_v1",
#"SingleMuonRun2016H_PromptReco_v2",
"SingleMuonRun2016H_PromptReco_v3_test"
];

nameDataEl = [
#"SingleElectronRun2016B_23Sep2016_v2",
#"SingleElectronRun2016C_23Sep2016_v1",
#"SingleElectronRun2016H_PromptReco_v1",
"SingleElectronRun2016H_PromptReco_v3_test"
];


nameData = {"el": nameDataEl, "mu":nameDataMu};



for a in range(len(category)):
    
    #MC
    if( doMC ):
        for i in range(len(samples)):
            fn = "WWAnalysisRun2/Job/job_"+samples[i][1]+"_"+category[a];
            outScript = open(fn+".sh","w");
            command = "python "+currentDir+"/WWAnalysisRun2/python/produceWWNtuples.py -i "+inputFolder+" -n "+str(samples[i][1])+" -o WWTree_"+str(samples[i][1])+"_"+category[a]+" -w "+str(samples[i][0])+" -lumi "+str(lumi)+" --ismc 1 -trig 1";
            print command;
            outScript.write('#!/bin/bash');
	    outScript.write("\n"+"workDir=`pwd`");
	    outScript.write("\n"+"echo `hostname`");
	    outScript.write("\n"+"echo $workDir");
            outScript.write("\n"+'cd '+CMSSWDir);
            outScript.write("\n"+'eval `scram runtime -sh`');
            #outScript.write("\n"+"cd WWAnalysis/WWAnalysisRun2");
	    outScript.write("\n"+"cp WWAnalysis/WWAnalysisRun2/python/produceWWNtuples.py WWAnalysis/WWAnalysisRun2/*.root ${workDir}");
	    outScript.write("\n"+"cd ${workDir}");
            outScript.write("\n"+command);
	    outScript.write("\n"+"echo \"====> LISTING ALL FILES..... \"");
	    #outScript.write("\n"+"ls");
	    outScript.write("\n"+"echo \"====> MOVING OUTPUT ROOT FILES\"");
	    outScript.write("\n"+"echo \"mv WWTree_"+(samples[i][1])+"_"+category[a]+".root "+outputFolder+"/output_"+category[a]+"\"");
	    outScript.write("\n"+"mv WWTree_"+(samples[i][1])+"_"+category[a]+".root "+outputFolder+"/output_"+category[a]);
	    outScript.write("\n"+"echo \"#########i########### Finished COPYING...\"");
            outScript.write("\n");
            outScript.close();
            os.system("chmod 777 "+currentDir+"/"+fn+".sh");
            command2 = "bsub -o out.%J -q 2nd -cwd "+currentDir+" "+currentDir+"/"+fn+".sh  -J "+str(samples[i][1][18:27]);
            print command2
            if( dryRun != True ):
                os.system(command2);
                time.sleep(3)
    
    #data
    if( doData ):
        for i in range(len(nameData[category[a]])):
            fn = "WWAnalysisRun2/Job/job_"+(nameData[category[a]])[i]+"_"+category[a];
            outScript = open(fn+".sh","w");
            command = "python "+currentDir+"/WWAnalysisRun2/python/produceWWNtuples.py -i "+inputFolder+" -n "+(nameData[category[a]])[i]+" -o WWTree_"+(nameData[category[a]])[i]+" -l "+category[a]+" -w 1. -no 1. -lumi "+lumi+ " --ismc 0 -trig 1";
            print command;
            outScript.write('#!/bin/bash');
            outScript.write("\n"+'cd '+CMSSWDir);
            outScript.write("\n"+'eval `scram runtime -sh`');
	    outScript.write("\n"+"cd WWAnalysis/WWAnalysisRun2");
            outScript.write("\n"+command);
	    outScript.write("\n"+"echo \"====> LISTING ALL FILES..... \"");
	    #outScript.write("\n"+"ls");
	    outScript.write("\n"+"echo \"====> MOVING OUTPUT ROOT FILES\"");
            outScript.write("\n"+"echo \"mv WWTree_"+(nameData[category[a]])[i]+".root "+outputFolder+"/output_"+category[a]+"\"");
            outScript.write("\n"+"mv WWTree_"+(nameData[category[a]])[i]+".root "+outputFolder+"/output_"+category[a]);
	    outScript.write("\n"+"echo \"#########i########### Finishedi COPYING...\"");
            outScript.write("\n");
            outScript.close();
            os.system("chmod 777 "+currentDir+"/"+fn+".sh");
            command2 = "bsub -q 2nd -cwd "+currentDir+" "+currentDir+"/"+fn+".sh";
            print command2
            if( dryRun != True ):
                os.system(command2);
                time.sleep(1)
