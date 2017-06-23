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

#inputFolder = "/eos/cms/store/caf/user/lbrianza/WWReducedTree_run2";
inputFolder = "/eos/uscms/store/user/lnujj/WpWm_aQGC_Ntuples_Ram/FirstStepOutput/Jan102016";
outputFolder = currentDir+"/output/";
exeName = currentDir+"/produceWWNtuples.exe"

dryRun = False;
doMC = True;
doData = False;

category = ["mu","el"];
#category = ["el"];
#category = ["mu"];



samples = [
    ( 0.25973,		"Signal_LT",		2471400,	0.),
    ( 0.04130,		"Signal_LL",		499800,		0.),
    ( 0.44750,		"Signal_TT",		4988000,	0.),
    ( 49.9970,		"WW_excl",		1999200,	0.),
    ( 49.9970,		"WW_excl_ext1",		6998600,	0.),
    ( 49.9970,		"WW_excl_amcatnlo",	5176114,	0.),
    ( 10.7100,		"WZ_excl_amcatnlo_2",	24221923,	0.),
    ( 10.7100,		"WZ_excl_amcatnlo_1",	24221923,	0.),
    ( 3.22000,		"ZZ_excl_amcatnlo",	15345572,	0.),
    ( 3.36000,		"sch",			1000000,	0.),
    ( 43.8000,		"tch_5",		67240808,	0.),
    ( 43.8000,		"tch_4",		67240808,	0.),
    ( 43.8000,		"tch_3",		67240808,	0.),
    ( 43.8000,		"tch_2",		67240808,	0.),
    ( 43.8000,		"tch_1",		67240808,	0.),
    ( 26.0700,          "tch_bar_3",  		38811017,	0.),
    ( 26.0700,          "tch_bar_2",  		38811017,	0.),
    ( 26.0700,          "tch_bar_1",  		38811017,	0.),
    ( 35.6000,		"tWch_ext1",		6952830,	0.),
    ( 35.6000,		"tWch_bar_ext1",	6933094,	0.),
    ( 831.760,		"TTbar_amcatnlo_4", 	23561608,    	0.),
    ( 831.760,		"TTbar_amcatnlo_3", 	23561608,    	0.),
    ( 831.760,		"TTbar_amcatnlo_2", 	23561608,    	0.),
    ( 831.760,		"TTbar_amcatnlo_1", 	23561608,    	0.),
    ( 831.760,		"TTbar_powheg", 	77229341,    	0.),
    ( 1627.45,       	"WJets100", 		10235198,    	0.),
    ( 1627.45,       	"WJets100ext1", 	29503700,    	0.),
    ( 1627.45,       	"WJets100ext2_2", 	39617787,    	0.),
    ( 1627.45,       	"WJets100ext2_1", 	39617787,    	0.),
    ( 435.24,       	"WJets200", 		4950373,    	0.),
    ( 435.24,       	"WJets200ext1", 	14815928,    	0.),
    ( 435.24,       	"WJets200ext2", 	19914590,    	0.),
    ( 59.18,       	"WJets400", 		1963464,    	0.),
    ( 59.18,       	"WJets400ext1", 	5796237,    	0.),
    ( 14.58,       	"WJets600", 		3779141,    	0.),
    ( 14.58,       	"WJets600ext1", 	14908339,    	0.),
    ( 6.655,       	"WJets800ext1", 	6200954,    	0.),
    ( 1.60809,       	"WJets1200", 		244532,    	0.),
    ( 1.60809,       	"WJets1200ext1", 	6627909,    	0.),
    ( 0.0389136,       	"WJets2500", 		253561,    	0.),
    ( 0.0389136,       	"WJets2500ext1", 	2384260,    	0.)
    #(    0.3363900000000,         "Higgs650",   399600,  650.),
    ##(    0.63980000000,         "Higgs750",    96200,  750.),
    #(    0.033540000000,        "Higgs1000",   400000, 1000.),
    #(    0.067650000000,      "VBFHiggs650",   400000,  650.),
    ##(    0.19150000000,      "VBFHiggs750",   100000,  750.),
    #(    0.023750000000,     "VBFHiggs1000",   399634, 1000.),
    #(    1.77400363645440040e-01,  "BulkGraviton600",    50000,  600.),
    #(    7.22813367744000179e-02, "BulkGraviton700",     48000,  700.),
    #(    4.79197963636363647e-02, "BulkGraviton750",     46000,  750.),
    #(    3.31548423242400067e-02,  "BulkGraviton800",    48000,  800.),
    #(    1.67515749198720032e-02, "BulkGraviton900",     49000,  900.),
    #(    8.99252925667200220e-03, "BulkGraviton1000",    50000, 1000.),
    #(    2.99686333939200040e-03, "BulkGraviton1200",    50000, 1200.),
    #(    1.15022191555440019e-03, "BulkGraviton1400",    50000, 1400.),
    #(    4.83220517270400096e-04, "BulkGraviton1600",    50000, 1600.),
    #(    2.18919227679360021e-04, "BulkGraviton1800",    50000, 1800.),
    #(    2.23940299210069489e-05, "BulkGraviton2000",    50000, 2000.),
    #(    4.19924014292041532e-06, "BulkGraviton2500",    50000, 2500.),
    #(    9.19110485489025090e-07, "BulkGraviton3000",    50000, 3000.),
    #(    2.25464182856213949e-07, "BulkGraviton3500",    50000, 3500.),
    #(    5.82648753152116075e-08, "BulkGraviton4000",    49800, 4000.),
    #(    0.000000000000000000000, "BulkGraviton4500",    50000, 4500.)
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
            fn = "Job/job_"+samples[i][1]+"_"+category[a];
            outScript = open(fn+".sh","w");
            command = "python "+currentDir+"/python/produceWWNtuples.py --exe "+exeName+" -i "+inputFolder+" -n "+str(samples[i][1])+" -o WWTree_"+str(samples[i][1])+" -l "+category[a]+" -w "+str(samples[i][0])+" -no "+str(samples[i][2])+" -mass "+str(samples[i][3])+" --ismc 1 -trig 0";
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
            outScript.write("\n"+"mv WWTree_"+(samples[i][1])+".root "+outputFolder+"/output_"+category[a]);
            outScript.write("\n");
            outScript.close();
            os.system("chmod 777 "+currentDir+"/"+fn+".sh");
            command2 = "bsub -q cmscaf1nd -cwd "+currentDir+" "+currentDir+"/"+fn+".sh";
            print command2
            if( dryRun != True ):
                os.system(command2);
                time.sleep(3)
    
    #data
    if( doData ):
        for i in range(len(nameData[category[a]])):
            fn = "Job/job_"+(nameData[category[a]])[i]+"_"+category[a];
            outScript = open(fn+".sh","w");
            command = "python "+currentDir+"/python/produceWWNtuples.py --exe "+exeName+" -i "+inputFolder+" -n "+(nameData[category[a]])[i]+" -o WWTree_"+(nameData[category[a]])[i]+" -l "+category[a]+" -w 1. -no 1. -mass 0 --ismc 0 -trig 1";
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
            command2 = "bsub -q cmscaf1nd -cwd "+currentDir+" "+currentDir+"/"+fn+".sh";
            print command2
            if( dryRun != True ):
                os.system(command2);
                time.sleep(3)
