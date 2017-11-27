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

#inputFolder = "/store/cmst3/group/monojet/production/12";
inputFolder = "/store/user/arapyan/Run2/";
#inputFolder = "/eos/cms/store/user/ksung/production/12/"
inputFolder = "/store/user/lnujj/WpWm_aQGC_Ntuples_Ram/FirstStepOutput/BaconNtuples_New";
outputFolder = "/store/user/lnujj/WpWm_aQGC_Ntuples_Ram/FirstStepOutput/BaconNtuples";

dryRun = False;
doMC = False;
doData = True;

sampleName = [
"DY1JetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8",
"DY2JetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8",
"DY3JetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8",
"DY4JetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8",
"DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8",
#"DYToLL_0J_13TeV-amcatnloFXFX-pythia8",
#"DYToLL_1J_13TeV-amcatnloFXFX-pythia8",
#"DYToLL_2J_13TeV-amcatnloFXFX-pythia8"
#"Summer16_ST_s_channel_4f_leptonDecays",
#"Summer16_ST_t_channel_antitop_4f_inclusiveDecays_TuneCUETP8M2T4",
#"Summer16_ST_t_channel_top_4f_inclusiveDecays_TuneCUETP8M2T4"
#"WplusToLNuWminusTo2JJJ_EWK_LO_SM_MJJ100PTJ10_TuneCUETP8M1_13TeV-madgraph-pythia8",	
#"WplusTo2JWminusToLNuJJ_EWK_LO_SM_MJJ100PTJ10_TuneCUETP8M1_13TeV-madgraph-pythia8",	
#"WplusToLNuWplusTo2JJJ_EWK_LO_SM_MJJ100PTJ10_TuneCUETP8M1_13TeV-madgraph-pythia8",	
#"WminusToLNuWminusTo2JJJ_EWK_LO_SM_MJJ100PTJ10_TuneCUETP8M1_13TeV-madgraph-pythia8",	
#"WplusToLNuZTo2JJJ_EWK_LO_SM_MJJ100PTJ10_TuneCUETP8M1_13TeV-madgraph-pythia8",		
#"WplusTo2JZTo2LJJ_EWK_LO_SM_MJJ100PTJ10_TuneCUETP8M1_13TeV-madgraph-pythia8",		
#"WminusToLNuZTo2JJJ_EWK_LO_SM_MJJ100PTJ10_TuneCUETP8M1_13TeV-madgraph-pythia8",		
#"WminusTo2JZTo2LJJ_EWK_LO_SM_MJJ100PTJ10_TuneCUETP8M1_13TeV-madgraph-pythia8",		
#"ZTo2LZTo2JJJ_EWK_LO_SM_MJJ100PTJ10_TuneCUETP8M1_13TeV-madgraph-pythia8",		
#
#"WplusTo2JWminusToLNuJJ_QCD_LO_SM_MJJ100PTJ10_TuneCUETP8M1_13TeV-madgraph-pythia8",	
#"WplusToLNuWminusTo2JJJ_QCD_LO_SM_MJJ100PTJ10_TuneCUETP8M1_13TeV-madgraph-pythia8",	
#"WplusToLNuWplusTo2JJJ_QCD_LO_SM_MJJ100PTJ10_TuneCUETP8M1_13TeV-madgraph-pythia8",	
#"WminusToLNuWminusTo2JJJ_QCD_LO_SM_MJJ100PTJ10_TuneCUETP8M1_13TeV-madgraph-pythia8",	
#"WplusToLNuZTo2JJJ_QCD_LO_SM_MJJ100PTJ10_TuneCUETP8M1_13TeV-madgraph-pythia8",		
#"WplusTo2JZTo2LJJ_QCD_LO_SM_MJJ100PTJ10_TuneCUETP8M1_13TeV-madgraph-pythia8",		
#"WminusToLNuZTo2JJJ_QCD_LO_SM_MJJ100PTJ10_TuneCUETP8M1_13TeV-madgraph-pythia8",		
#"WminusTo2JZTo2LJJ_QCD_LO_SM_MJJ100PTJ10_TuneCUETP8M1_13TeV-madgraph-pythia8",		
#"ZTo2LZTo2JJJ_QCD_LO_SM_MJJ100PTJ10_TuneCUETP8M1_13TeV-madgraph-pythia8"	
#
#"WplusToLNuWminusTo2JJJ_EWK_LO_aQGC_MJJ100PTJ10_TuneCUETP8M1_13TeV-madgraph-pythia8",	
#"WplusTo2JWminusToLNuJJ_EWK_LO_aQGC_MJJ100PTJ10_TuneCUETP8M1_13TeV-madgraph-pythia8",	
#"WplusToLNuWplusTo2JJJ_EWK_LO_aQGC_MJJ100PTJ10_TuneCUETP8M1_13TeV-madgraph-pythia8",	
#"WminusToLNuWminusTo2JJJ_EWK_LO_aQGC_MJJ100PTJ10_TuneCUETP8M1_13TeV-madgraph-pythia8",	
#"WplusToLNuZTo2JJJ_EWK_LO_aQGC_MJJ100PTJ10_TuneCUETP8M1_13TeV-madgraph-pythia8",	
#"WplusTo2JZTo2LJJ_EWK_LO_aQGC_MJJ100PTJ10_TuneCUETP8M1_13TeV-madgraph-pythia8",		
#"WminusToLNuZTo2JJJ_EWK_LO_aQGC_MJJ100PTJ10_TuneCUETP8M1_13TeV-madgraph-pythia8",	
#"WminusTo2JZTo2LJJ_EWK_LO_aQGC_MJJ100PTJ10_TuneCUETP8M1_13TeV-madgraph-pythia8",	
#"ZTo2LZTo2JJJ_EWK_LO_aQGC_MJJ100PTJ10_TuneCUETP8M1_13TeV-madgraph-pythia8",	
###"DY1JetsToLL_M_50_13TeV",
###"DY2JetsToLL_M_50_13TeV",
###"DY3JetsToLL_M_50_13TeV",
###"DY4JetsToLL_M_50_13TeV",
###"DYJetsToLL_M-10to50_amcatnlo",
###"DYJetsToLL_M-10to50_amcatnlo_ext1",
###"DYJetsToLL_M_50_13TeV_ext",
###"DYJetsToLL_M_50_13TeV_Fast",
###"DYJetsToLL_M-50_amcatnlo",
###"DYJetsToLL_M-50_HT-100to200",
###"DYJetsToLL_M-50_HT-100to200_ext1",
###"DYJetsToLL_M-50_HT-1200to2500",
###"DYJetsToLL_M-50_HT-200to400",
###"DYJetsToLL_M-50_HT-200to400_ext1",
###"DYJetsToLL_M-50_HT-2500toInf",
###"DYJetsToLL_M-50_HT-400to600",
###"DYJetsToLL_M-50_HT-400to600_ext1",
###"DYJetsToLL_M-50_HT-600to800",
###"DYJetsToLL_M-50_HT-70to100",
###"DYJetsToLL_M-50_HT-800to1200",
###"DYToLL_0J_amcatnlo",
###"DYToLL_0J_amcatnlo_backup",
###"DYToLL_1J_amcatnlo",
###"DYToLL_1J_amcatnlo_backup",
###"DYToLL_2J_amcatnlo",
###"DYToLL_2J_amcatnlo_backup",
##"QCD_HT1000to1500_13TeV",
##"QCD_HT1000to1500_13TeV_ext",
##"QCD_HT100to200_13TeV",
##"QCD_HT1500to2000_13TeV",
##"QCD_HT1500to2000_13TeV_ext",
##"QCD_HT2000toInf_13TeV",
##"QCD_HT2000toInf_13TeV_ext",
##"QCD_HT200to300_13TeV",
##"QCD_HT200to300_13TeV_ext",
##"QCD_HT300to500_13TeV",
##"QCD_HT300to500_13TeV_ext",
##"QCD_HT500to700_13TeV",
##"QCD_HT500to700_13TeV_ext",
##"QCD_HT50to100_13TeV",
##"QCD_HT700to1000_13TeV",
##"QCD_HT700to1000_13TeV_ext",
##"ST_t_channel_top_4f_scaleup_inclusiveDecays_13TeV_powhegV2_madspin_pythia8",
##"ST_tW_antitop_5f_inclusiveDecays_13TeV_powheg_pythia8_TuneCUETP8M2T4",
##"ST_tW_antitop_5f_NoFullyHadronicDecays_13TeV_powheg_TuneCUETP8M1",
##"ST_tW_top_5f_inclusiveDecays_13TeV_powheg_pythia8_TuneCUETP8M2T4",
##"ST_tW_top_5f_NoFullyHadronicDecays_13TeV_powheg_TuneCUETP8M1",
##"TT_powheg",
###"TT_powheg_backup",
##"TTTo2L2Nu_powheg",
##"TTTo2L2Nu_powheg_pythia8",
###"TTToSemilepton_powheg",
##"TTWJetsToLNu_ext1",
##"TTWJetsToLNu_ext2",
##"TTWJetsToQQ",
##"TTZToLLNuNu_M-10",
##"TTZToLLNuNu_M-10_ext2",
##"TTZToQQ",
###"W1JetsToLNu_13TeV_v1",
###"W2JetsToLNu_13TeV_v1",
###"W3JetsToLNu_13TeV_v1",
###"W4JetsToLNu_13TeV_v1",
###"WJetsToLNu_13TeV",
###"WJetsToLNu_HT_100To200_13TeV",
###"WJetsToLNu_HT_100To200_13TeV_ext1",
###"WJetsToLNu_HT_100To200_13TeV_ext2",
###"WJetsToLNu_HT_1200To2500_13TeV",
###"WJetsToLNu_HT_1200To2500_13TeV_ext1",
###"WJetsToLNu_HT_200To400_13TeV",
###"WJetsToLNu_HT_200To400_13TeV_ext1",
###"WJetsToLNu_HT_200To400_13TeV_ext2",
###"WJetsToLNu_HT_2500ToInf_13TeV",
###"WJetsToLNu_HT_2500ToInf_13TeV_ext1",
###"WJetsToLNu_HT_400To600_13TeV",
###"WJetsToLNu_HT_400To600_13TeV_ext1",
###"WJetsToLNu_HT_600To800_13TeV",
###"WJetsToLNu_HT_600To800_13TeV_ext1",
###"WJetsToLNu_HT_70To100_13TeV",
###"WJetsToLNu_HT_800To1200_13TeV",
###"WJetsToLNu_HT_800To1200_13TeV_ext1",
##"WW_13TeV_pythia8",
##"WWTo1L1Nu2Q_13TeV_amcatnloFXFX_madspin_pythia8",
##"WWTo2L2Nu_13TeV_powheg",
##"WWTo4Q_13TeV_powheg",
##"WWToLNuQQ_13TeV_powheg",
##"WWToLNuQQ_13TeV_powheg_ext",
##"WWW_4F_13TeV_amcatnlo_pythia8",
##"WWZ_13TeV_amcatnlo_pythia8",
##"WZ_13TeV_pythia8",
##"WZTo1L1Nu2Q_13TeV_amcatnloFXFX_madspin_pythia8",
##"WZTo2L2Q_13TeV_amcatnloFXFX_madspin_pythia8",
##"WZToLNu2Q_13TeV_powheg_pythia8",
##"ZZ_13TeV_pythia8",
##"ZZTo2L2Nu_13TeV_powheg_pythia8",
##"ZZTo2L2Q_13TeV_amcatnloFXFX_madspin_pythia8",
##"ZZTo2Q2Nu_13TeV_amcatnloFXFX_madspin_pythia8",
##"ZZto4L_13TeV_amcatnloFXFX_pythia8",
##"ZZTo4Q_13TeV_amcatnloFXFX_madspin_pythia8",
##"ZZZ_13TeV_amcatnlo_pythia8"
]

command1 = "xrdfs root://cmseos.fnal.gov/ mkdir "+outputFolder+"/$1/"
command2 = "xrdcp -r root://cmseos.fnal.gov/"+inputFolder+"/$1/  root://cmseos.fnal.gov/"+outputFolder+"/$1/"
#command2 = "xrdcp -r -f -s root://eoscms.cern.ch/"+inputFolder+"/$1/  root://cmseos.fnal.gov/"+outputFolder+"/$1/"

outScript = open("runCopycondor.sh","w");
outScript.write('#!/bin/bash');
outScript.write("\n"+'cd '+CMSSWDir);
outScript.write("\n"+'eval `scram runtime -sh`');
outScript.write("\n"+"cd -");
outScript.write("\n"+"echo \"sample Name : \"$1");
outScript.write("\n"+command1);
outScript.write("\n"+command2);
outScript.write("\n");
outScript.close();
os.system("chmod 777 runCopycondor.sh");

outJDL = open("runCopycondor.jdl","w");
outJDL.write("Executable = runCopycondor.sh\n");
outJDL.write("Universe = vanilla\n");
outJDL.write("Requirements =FileSystemDomain==\"fnal.gov\" && Arch==\"X86_64\"");
outJDL.write("\n");
outJDL.write("Notification = ERROR\n");
outJDL.write("Should_Transfer_Files = YES\n");
outJDL.write("WhenToTransferOutput = ON_EXIT\n");
#outJDL.write("include : list-infiles.sh |\n");
outJDL.write("x509userproxy = $ENV(X509_USER_PROXY)\n");

for i in range(len(sampleName)):
    outJDL.write("Output = "+str(sampleName[i])+"_Copy.stdout\n");
    outJDL.write("Error = "+str(sampleName[i])+"_Copy.stdout\n");
    outJDL.write("Arguments = "+str(sampleName[i])+"\n");
    outJDL.write("Queue\n");

outJDL.close();
print "===> Set Proxy Using:";
print "\tvoms-proxy-init --voms cms --valid 168:00";
print "\"condor_submit runCopycondor.jdl\" to submit";

