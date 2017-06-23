hadd -f WWTree_data_mu_2016_runB_v3_mu.root WWTree_data_mu_2016_runB_v3_*_mu.root
hadd -f WWTree_data_mu_2016_runC_v1_mu.root WWTree_data_mu_2016_runC_v1_*_mu.root
hadd -f WWTree_data_mu_2016_runD_v1_mu.root WWTree_data_mu_2016_runD_v1_*_mu.root
hadd -f WWTree_data_mu_2016_runE_v1_mu.root WWTree_data_mu_2016_runE_v1_*_mu.root
hadd -f WWTree_data_mu_2016_runF_v1_mu.root WWTree_data_mu_2016_runF_v1_*_mu.root
hadd -f WWTree_data_mu_2016_runG_v1_mu.root WWTree_data_mu_2016_runG_v1_*_mu.root
hadd -f WWTree_data_mu_2016_runH_v2_mu.root WWTree_data_mu_2016_runH_v2_*_mu.root
hadd -f WWTree_data_golden_B3C1D1E1F1G1H3.root WWTree_data_mu_2016_runB_v3_mu.root WWTree_data_mu_2016_runC_v1_mu.root WWTree_data_mu_2016_runD_v1_mu.root WWTree_data_mu_2016_runE_v1_mu.root WWTree_data_mu_2016_runF_v1_mu.root WWTree_data_mu_2016_runG_v1_mu.root WWTree_data_mu_2016_runH_v3_mu.root WWTree_data_mu_2016_runH_v2_mu.root
mv WWTree_data_mu_2016_runB_v3_*mu.root others/
mv WWTree_data_mu_2016_runC_v1_*mu.root others/
mv WWTree_data_mu_2016_runD_v1_*mu.root others/
mv WWTree_data_mu_2016_runE_v1_*mu.root others/
mv WWTree_data_mu_2016_runF_v1_*mu.root others/
mv WWTree_data_mu_2016_runG_v1_*mu.root others/
mv WWTree_data_mu_2016_runH_v2_*mu.root others/
mv WWTree_data_mu_2016_runH_v3_mu.root others/

hadd -f WWTree_DYJetsToLL_amcatnlo_ext1_mu.root WWTree_DYJetsToLL_amcatnlo_ext1_1_mu.root WWTree_DYJetsToLL_amcatnlo_ext1_2_mu.root
hadd -f WWTree_Signal_LLpLTpTT_mu.root	WWTree_Signal_LL_mu.root WWTree_Signal_LT_mu.root WWTree_Signal_TT_mu.root
hadd -f WWTree_tch_mu.root	WWTree_tch_1_mu.root	WWTree_tch_2_mu.root	WWTree_tch_3_mu.root	WWTree_tch_4_mu.root	WWTree_tch_5_mu.root	
hadd -f WWTree_tch_bar_mu.root	WWTree_tch_bar_1_mu.root	WWTree_tch_bar_2_mu.root	WWTree_tch_bar_3_mu.root	
hadd -f WWTree_TTbar_amcatnlo_mu.root	WWTree_TTbar_amcatnlo_1_mu.root	WWTree_TTbar_amcatnlo_2_mu.root	WWTree_TTbar_amcatnlo_3_mu.root	WWTree_TTbar_amcatnlo_4_mu.root	
hadd -f WWTree_TTbar_powheg_mu.root	WWTree_TTbar_powheg_1_mu.root	WWTree_TTbar_powheg_2_mu.root	WWTree_TTbar_powheg_3_mu.root	WWTree_TTbar_powheg_4_mu.root	
hadd -f WWTree_tWch_bara_mu.root	WWTree_tWch_bar_mu.root	WWTree_tWch_bar_ext1_mu.root	
hadd -f WWTree_tWcha_mu.root	WWTree_tWch_mu.root	WWTree_tWch_ext1_mu.root	
hadd -f WWTree_WJets100a_mu.root	WWTree_WJets100_mu.root	WWTree_WJets100ext1_mu.root	WWTree_WJets100ext2_1_mu.root	WWTree_WJets100ext2_2_mu.root	
hadd -f WWTree_WJets1200a_mu.root	WWTree_WJets1200_mu.root	WWTree_WJets1200ext1_mu.root	
hadd -f WWTree_WJets200a_mu.root	WWTree_WJets200_mu.root	WWTree_WJets200ext1_mu.root	WWTree_WJets200ext2_mu.root	
hadd -f WWTree_WJets2500a_mu.root	WWTree_WJets2500_mu.root	WWTree_WJets2500ext1_mu.root	
hadd -f WWTree_WJets400a_mu.root	WWTree_WJets400_mu.root	WWTree_WJets400ext1_mu.root	
hadd -f WWTree_WJets600a_mu.root	WWTree_WJets600_mu.root	WWTree_WJets600ext1_mu.root	
hadd -f WWTree_WW_excla_mu.root	WWTree_WW_excl_mu.root	WWTree_WW_excl_ext1_mu.root	
hadd -f WWTree_WZ_excl_amcatnlo_mu.root	WWTree_WZ_excl_amcatnlo_1_mu.root	WWTree_WZ_excl_amcatnlo_2_mu.root	

mv WWTree_DYJetsToLL_amcatnlo_ext1_1_mu.root WWTree_DYJetsToLL_amcatnlo_ext1_2_mu.root others/
#mv WWTree_Signal_LL_mu.root WWTree_Signal_LT_mu.root WWTree_Signal_TT_mu.root others/
mv WWTree_tch_1_mu.root	WWTree_tch_2_mu.root	WWTree_tch_3_mu.root	WWTree_tch_4_mu.root	WWTree_tch_5_mu.root	 others/
mv WWTree_tch_bar_1_mu.root	WWTree_tch_bar_2_mu.root	WWTree_tch_bar_3_mu.root	 others/
mv WWTree_TTbar_amcatnlo_1_mu.root	WWTree_TTbar_amcatnlo_2_mu.root	WWTree_TTbar_amcatnlo_3_mu.root	WWTree_TTbar_amcatnlo_4_mu.root	 others/
mv WWTree_TTbar_powheg_1_mu.root	WWTree_TTbar_powheg_2_mu.root	WWTree_TTbar_powheg_3_mu.root	WWTree_TTbar_powheg_4_mu.root	 others/
mv WWTree_tWch_bar_mu.root	WWTree_tWch_bar_ext1_mu.root	 others/
mv WWTree_tWch_mu.root	WWTree_tWch_ext1_mu.root	 others/
mv WWTree_WJets100_mu.root	WWTree_WJets100ext1_mu.root	WWTree_WJets100ext2_1_mu.root	WWTree_WJets100ext2_2_mu.root	 others/
mv WWTree_WJets1200_mu.root	WWTree_WJets1200ext1_mu.root	 others/
mv WWTree_WJets200_mu.root	WWTree_WJets200ext1_mu.root	WWTree_WJets200ext2_mu.root	 others/
mv WWTree_WJets2500_mu.root	WWTree_WJets2500ext1_mu.root	 others/
mv WWTree_WJets400_mu.root	WWTree_WJets400ext1_mu.root	 others/
mv WWTree_WJets600_mu.root	WWTree_WJets600ext1_mu.root	 others/
mv WWTree_WW_excl_mu.root	WWTree_WW_excl_ext1_mu.root	 others/
mv WWTree_WZ_excl_amcatnlo_1_mu.root	WWTree_WZ_excl_amcatnlo_2_mu.root	 others/
