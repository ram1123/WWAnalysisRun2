import os

inPath = "/eos/uscms/store/user/rasharma/SecondStep/WWTree_2018_01_07_12h02/HaddedFiles/"
outPath = inPath+"Hadds_for_BkgEstimation/"

os.system("hadd "+outPath+"WWTree_data_golden.root " + inPath + "Data.root")
os.system("hadd "+outPath+"WWTree_VV.root  " + inPath + "*_QCD_LO_SM.root")
os.system("hadd "+outPath+"QCD_HTbin.root  " + inPath + "QCD_HT*.root")
os.system("hadd "+outPath+"WWTree_STop.root  " + inPath + "ST_s_channel.root  " + inPath + "ST_tW_antitop_5f_NoFullyHadronicDecays.root  " + inPath + "ST_tW_top_5f_NoFullyHadronicDecays.root  " + inPath + "ST_t_channel_antitop.root ST_t_channel_top_4f.root ")
os.system("hadd "+outPath+"WWTree_TTbar.root  " + inPath + "TTToSemilepton.root")
os.system("hadd "+outPath+"WWTree_WJets.root  " + inPath + "WJetsToLNu_HT_*.root")
os.system("hadd "+outPath+"WWTree_Signal_aQGC.root  " + inPath + "*_EWK_LO_aQGC.root")
os.system("hadd "+outPath+"WWTree_Signal_SM.root   " + inPath + "*_EWK_LO_SM.root")
os.system("hadd "+outPath+"PseudoData.root " + outPath + "WWTree_VV.root  " + outPath + "QCD_HTbin.root  " + outPath + "WWTree_STop.root  " + outPath + "WWTree_TTbar.root  " + outPath + "WWTree_WJets.root  " + outPath + "WWTree_Signal_SM.root")
