import os

#inPath = "/eos/uscms/store/user/rasharma/SecondStep/WWTree_CleanedCode_Isolated_NaNFixed_Btag30GeV_AlphaRatioBkgEst_2018_03_26_03h12/HaddedFiles/"
#inPath = "/eos/uscms/store/user/rasharma/SecondStep/WWTree_CleanedCode_Isolated_NaNFixed_Btag30GeV_AlphaRatioBkgEst_2018_03_27_02h28/HaddedFiles/"
inPath = "/eos/uscms/store/user/rasharma/SecondStep/WWTree_CommonNtuple_For1and2Lepton_2018_05_15_04h15/HaddedFiles/"
outPath = inPath+"Hadds_for_BkgEstimation/"

#os.system("mkdir outPath")

os.system("hadd -f "+outPath+"WWTree_data_golden.root " + inPath + "Data.root")
os.system("hadd -f "+outPath+"WWTree_VV.root  " + inPath + "*_QCD_LO_SM.root")
os.system("hadd -f "+outPath+"QCD_HTbin.root  " + inPath + "QCD_HT*.root")
os.system("hadd -f "+outPath+"WWTree_STop.root  " + inPath + "ST_s_channel.root  " + inPath + "ST_tW_antitop_5f_NoFullyHadronicDecays.root  " + inPath + "ST_tW_top_5f_NoFullyHadronicDecays.root  " + inPath + "ST_t_channel_antitop.root " + inPath + "ST_t_channel_top_4f.root ")
os.system("hadd -f "+outPath+"WWTree_TTbar.root  " + inPath + "TTToSemilepton.root")
os.system("hadd -f "+outPath+"WWTree_WJets.root  " + inPath + "WJetsToLNu_HT_*.root")
os.system("hadd -f "+outPath+"WWTree_ZJets.root  " + inPath + "DY*JetsToLL.root")
os.system("hadd -f "+outPath+"WWTree_VJets.root  " + inPath + "WJetsToLNu_HT_*.root "+ inPath + "DY*JetsToLL.root")
os.system("hadd -f "+outPath+"WWTree_Signal_aQGC.root  " + inPath + "*_EWK_LO_aQGC.root")
os.system("hadd -f "+outPath+"WWTree_Signal_SM.root   " + inPath + "*_EWK_LO_SM.root")
os.system("hadd -f "+outPath+"PseudoData.root " + outPath + "WWTree_VV.root  " + outPath + "QCD_HTbin.root  " + outPath + "WWTree_STop.root  " + outPath + "WWTree_TTbar.root  " + outPath + "WWTree_WJets.root  " + outPath + "WWTree_Signal_SM.root")
