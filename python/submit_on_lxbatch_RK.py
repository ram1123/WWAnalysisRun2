#! /usr/bin/env python
import os
import glob
import math
from array import array
import sys
import time
import subprocess

currentDir = os.getcwd();
CMSSWDir = currentDir;
ReducedTreeDir = "";

ExtraComment = "";

name = [ "WJets", "TTbar_powheg","ST_t-channel_top_4f_leptonDecays_13TeV-powheg-pythia8_TuneCUETP8M1", "ST_t-channel_antitop_4f_leptonDecays_13TeV-powheg-pythia8_TuneCUETP8M1", "tWch", "tWch_bar", "WW", "ZZ", "WZ",
	"WWJJToLNuQQ_LL_13TeV-madgraph-pythia8", "WWJJToLNuQQ_LT_13TeV-madgraph-pythia8", "WWJJToLNuQQ_TT_13TeV-madgraph-pythia8"
	
#        "WJets","TTbar_amcatnlo", , "tch", "tch_bar", "tWch", "tWch_bar", "WW", "ZZ", "WZ",
#        "WW_excl", "WZ_excl", "ZZ_excl", "TTbar_powheg", "TTbar_madgraph",
#        "WJets_50ns", "WW_50ns", "WZ_50ns", "ZZ_50ns", "TTbar_amcatnlo_50ns", "TTbar_powheg_50ns", "TTbar_madgraph_50ns",
#        "tch_50ns", "tch_bar_50ns", "tWch_bar_50ns", "tWch_50ns",
];

category = ["mu","el"];
xSecWeight = [              "61526.7", "831.76",  "43.8", "26.07", "35.6", "35.6", "118.7", "15.4", "66.1",
		"0.041300", "0.25973", "0.44750"
			
#              "43.53", "10.96", "3.38", "831.76", "831.76",
#              "61526.7", "118.7", "15.4", "66.1", "831.76", "831.76", "831.76",
#              "43.8", "26.07", "35.6", "35.6",
];

N = [     "24151270.", "42730273.", "2966200.", "1695400.", "995600.", "1000000.", "994416.", "996168.", "991232.",
	"499600", "2496000", "4998400"
#     "1969600.", "24711046.", "18898680.", "19899500.", "11339232.",
#     "24089991.", "989608.", "998848.", "996920.", "4994250.", "19665194.", "4992231.",
#     "1273800.", "681900.", "1000000.", "998400.",
];

mass = [        "0","0", "0", "0", "0", "0", "0", "0", "0",
		"0","0","0"
#        "0", "0", "0", "0", "0",
#        "0", "0", "0", "0", "0", "0", "0",
#        "0", "0", "0", "0",
];

nameData = ["data_mu_prompt_v3_25ns_runD_ReMiniAOD","data_el_prompt_v3_25ns_runD_ReMiniAOD","data_mu_prompt_v4_25ns_runD","data_el_prompt_v4_25ns_runD"];

#MC
for a in range(len(category)):
    for i in range(len(name)):
        fn = "Job/Job_"+name[i]+"_"+ExtraComment+"_"+category[a];
        outScript = open(fn+".csh","w");
        command = "python produceWWNtuples.py -mc True -l "+category[a]+" -n "+name[i]+" -o WWTree_"+name[i]+" -w "+xSecWeight[i]+" -no "+N[i]+" -mass "+mass[i];
        #command = "python python/produceWWNtuples.py -mc True -l "+category[a]+" -n "+name[i]+" -o WWTree_"+name[i]+" -w "+xSecWeight[i]+" -no "+N[i]+" -mass "+mass[i];
        #command = "python python/produceWWNtuples.py -n "+name[i]+" -o WWTree_"+name[i]+" -l "+category[a]+" -w "+xSecWeight[i]+" -no "+N[i]+" -mass "+mass[i]+" --ismc True -trig 1";
        print command;
        outScript.write('#!/bin/tcsh');
	outScript.write(os.linesep+'source /cvmfs/cms.cern.ch/cmsset_default.csh');
        outScript.write(os.linesep+'cd '+CMSSWDir);
	outScript.write(os.linesep+'setenv SCRAM_ARCH slc6_amd64_gcc472');
        #outScript.write(os.linesep+'eval `scram runtime -sh`');
	outScript.write(os.linesep+'cmsenv');
        #outScript.write(os.linesep+'cd '+CMSSWDir);
	outScript.write(os.linesep+'cd ${_CONDOR_SCRATCH_DIR}');
	outScript.write(os.linesep+'setenv LD_LIBRARY_PATH /uscms_data/d3/rasharma/WWScattering/CMSSW_7_1_5/biglib/slc6_amd64_gcc481:/uscms_data/d3/rasharma/WWScattering/CMSSW_7_1_5/lib/slc6_amd64_gcc481:/uscms_data/d3/rasharma/WWScattering/CMSSW_7_1_5/external/slc6_amd64_gcc481/lib:/cvmfs/cms.cern.ch/slc6_amd64_gcc481/cms/cmssw/CMSSW_7_1_5/biglib/slc6_amd64_gcc481:/cvmfs/cms.cern.ch/slc6_amd64_gcc481/cms/cmssw/CMSSW_7_1_5/lib/slc6_amd64_gcc481:/cvmfs/cms.cern.ch/slc6_amd64_gcc481/cms/cmssw/CMSSW_7_1_5/external/slc6_amd64_gcc481/lib:/cvmfs/cms.cern.ch/slc6_amd64_gcc481/external/llvm/3.4-cms2/lib:/cvmfs/cms.cern.ch/slc6_amd64_gcc481/external/gcc/4.8.1/lib64:/cvmfs/cms.cern.ch/slc6_amd64_gcc481/external/gcc/4.8.1/lib\:${_CONDOR_SCRATCH_DIR}');
        outScript.write(os.linesep+command+os.linesep);
        outScript.close();
        os.system("chmod 777 "+currentDir+"/"+fn+".csh");
        #command2 = "bsub -q 1nd -cwd "+currentDir+" "+currentDir+"/"+fn+".sh";
        #os.system(command2);
        #print command2

for a in range(len(category)):
    for i in range(len(name)):
		fn = "Job/Condor_"+name[i]+"_"+ExtraComment+"_"+category[a];
		outScript = open(fn,"w");
		outScript.write('universe = vanilla');
		outScript.write(os.linesep+'Executable = Job/Job_'+name[i]+"_"+ExtraComment+'_'+category[a]+'.csh');
		outScript.write(os.linesep+'Should_Transfer_Files = YES');
		outScript.write(os.linesep+'WhenToTransferOutput = ON_EXIT');
		outScript.write(os.linesep+'transfer_input_files = '+CMSSWDir+'/produceWWNtuples.py, '+CMSSWDir+'/produceWWNtuples.exe, '+CMSSWDir+'/PU.root, '+CMSSWDir+'/'+name[i]+'.txt');
		outScript.write(os.linesep+'Output		= condor_logs/'+name[i]+"_"+category[a]+"_"+ExtraComment+'_$(Cluster)_$(Process).stdout')
		outScript.write(os.linesep+'Error		= condor_logs/'+name[i]+"_"+category[a]+"_"+ExtraComment+'_$(Cluster)_$(Process).stderr');
		outScript.write(os.linesep+'Log		= condor_logs/'+name[i]+"_"+category[a]+"_"+ExtraComment+'_$(Cluster)_$(Process).log');
		outScript.write(os.linesep+'Queue 1'+os.linesep);
		outScript.close();
		command2 = "condor_submit "+currentDir+"/"+fn;
		print command2
		os.system(command2);


#data
for a in range(len(category)):
    for i in range(len(nameData)):
        fn = "Job/Job_"+nameData[i]+"_"+ExtraComment+"_"+category[a];
        outScript = open(fn+".csh","w");
        command = "python produceWWNtuples.py -n "+nameData[i]+" -o WWTree_"+nameData[i]+" -l "+category[a]+" -w 1. -no 1. -mass 0 --ismc False -trig 1";
        #command = "python python/produceWWNtuples.py -n "+nameData[i]+" -o WWTree_"+nameData[i]+" -l "+category[a]+" -w 1. -no 1. -mass 0 --ismc False -trig 1";
        print command;
        outScript.write('#!/bin/tcsh');
	outScript.write(os.linesep+'source /cvmfs/cms.cern.ch/cmsset_default.csh');
        outScript.write(os.linesep+'cd '+CMSSWDir);
	outScript.write(os.linesep+'setenv SCRAM_ARCH slc6_amd64_gcc472');
        outScript.write(os.linesep+'eval `scram runtime -sh`');
        #outScript.write(os.linesep+'cd '+CMSSWDir);
	outScript.write(os.linesep+'cmsenv');
	outScript.write(os.linesep+'cd ${_CONDOR_SCRATCH_DIR}');
	outScript.write(os.linesep+'setenv LD_LIBRARY_PATH /uscms_data/d3/rasharma/WWScattering/CMSSW_7_1_5/biglib/slc6_amd64_gcc481:/uscms_data/d3/rasharma/WWScattering/CMSSW_7_1_5/lib/slc6_amd64_gcc481:/uscms_data/d3/rasharma/WWScattering/CMSSW_7_1_5/external/slc6_amd64_gcc481/lib:/cvmfs/cms.cern.ch/slc6_amd64_gcc481/cms/cmssw/CMSSW_7_1_5/biglib/slc6_amd64_gcc481:/cvmfs/cms.cern.ch/slc6_amd64_gcc481/cms/cmssw/CMSSW_7_1_5/lib/slc6_amd64_gcc481:/cvmfs/cms.cern.ch/slc6_amd64_gcc481/cms/cmssw/CMSSW_7_1_5/external/slc6_amd64_gcc481/lib:/cvmfs/cms.cern.ch/slc6_amd64_gcc481/external/llvm/3.4-cms2/lib:/cvmfs/cms.cern.ch/slc6_amd64_gcc481/external/gcc/4.8.1/lib64:/cvmfs/cms.cern.ch/slc6_amd64_gcc481/external/gcc/4.8.1/lib\:${_CONDOR_SCRATCH_DIR}');
        outScript.write(os.linesep+command+os.linesep);
        outScript.close();
        os.system("chmod 777 "+currentDir+"/"+fn+".csh");
        #command2 = "bsub -q 1nd -cwd "+currentDir+" "+currentDir+"/"+fn+".sh";
        #os.system(command2);
        #print command2

for a in range(len(category)):
	for i in range(len(nameData)):
		fn = "Job/Condor_"+nameData[i]+"_"+ExtraComment+"_"+category[a];
		outScript = open(fn,"w");
		outScript.write('universe = vanilla');
		outScript.write(os.linesep+'Executable = Job/Job_'+nameData[i]+"_"+ExtraComment+'_'+category[a]+'.csh');
		outScript.write(os.linesep+'Should_Transfer_Files = YES');
		outScript.write(os.linesep+'WhenToTransferOutput = ON_EXIT');
		outScript.write(os.linesep+'transfer_input_files = '+CMSSWDir+'/produceWWNtuples.py, '+CMSSWDir+'/produceWWNtuples.exe, '+CMSSWDir+'/PU.root, '+CMSSWDir+'/'+nameData[i]+'.txt');
		#outScript.write(os.linesep+'transfer_input_files = $pwd/python/produceWWNtuples.py, $pwd/produceWWNtuples.exe, $pwd/PU.root');
		outScript.write(os.linesep+'Output		= condor_logs/'+nameData[i]+"_"+category[a]+"_"+ExtraComment+'_$(Cluster)_$(Process).stdout');
		outScript.write(os.linesep+'Error		= condor_logs/'+nameData[i]+"_"+category[a]+"_"+ExtraComment+'_$(Cluster)_$(Process).stderr');
		outScript.write(os.linesep+'Log		= condor_logs/'+nameData[i]+"_"+category[a]+"_"+ExtraComment+'_$(Cluster)_$(Process).log');
		outScript.write(os.linesep+'Queue 1'+os.linesep);
		outScript.close();
		command2 = "condor_submit "+currentDir+"/"+fn;
		print command2
		os.system(command2);
