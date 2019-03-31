today=`date +%d%m%Y_%H%M%S`
rm *.dat
for d in *.stdout;do 
	#echo $d
	tempVar=`grep "TreeMaker2/PreSelection" $d | awk '{print $12}'`
	#echo $tempVar
	tempFile=`grep "TreeMaker2/PreSelection" $d | awk '{print $10}'`
	#echo $tempFile
	if [[ "$tempVar" == "el" ]];then
		#echo $tempFile
		#awk '{ if (lines > 0) {print; --lines; }} /SUMMARY/ {lines = 13}' < $d | awk 'FNR > 5' | awk '{print $(NF)}'
		awk '{ if (lines > 0) {print; --lines; }} /SUMMARY/ {lines = 13}' < $d | awk 'FNR > 5' | awk '{print $(NF)}' >> ${tempFile}.dat
	fi
	if [[ "$tempVar" == "mu" ]];then
		#echo $tempFile
		awk '{ if (lines > 0) {print; --lines; }} /SUMMARY/ {lines = 13}' < $d | awk 'FNR > 5' | awk '{print $(NF)}' >> ${tempFile}.dat
	fi
	#echo -e "\n==============================================\n"
done


echo "Wjet_HTbin_inclusive" > Wjet_HTbin_inclusive_el.dat
awk '{a[FNR]+=$1;} END{for (i=1; i<=FNR; i++) print a[i]}' WWTree_WJets100*el.dat WWTree_WJets200*el.dat WWTree_WJets400*el.dat WWTree_WJets600*el.dat WWTree_WJets800ext1_el.dat WWTree_WJets1200*el.dat WWTree_WJets2500*el.dat >> Wjet_HTbin_inclusive_el.dat
echo "DYJetsToLL" > DYJetsToLL_el.dat
awk '{a[FNR]+=$1;} END{for (i=1; i<=FNR; i++) print a[i]}' WWTree_DYJetsToLL_*el.dat >> DYJetsToLL_el.dat
echo "WJets100_el" > WJets100_el.dat
awk '{a[FNR]+=$1;} END{for (i=1; i<=FNR; i++) print a[i]}' WWTree_WJets100*el.dat  >> WJets100_el.dat
echo "WJets1200_el" > WJets1200_el.dat
awk '{a[FNR]+=$1;} END{for (i=1; i<=FNR; i++) print a[i]}' WWTree_WJets1200*el.dat >> WJets1200_el.dat
echo "WJets200_el" > WJets200_el.dat
awk '{a[FNR]+=$1;} END{for (i=1; i<=FNR; i++) print a[i]}' WWTree_WJets200*el.dat >> WJets200_el.dat
echo "WJets2500_el" > WJets2500_el.dat
awk '{a[FNR]+=$1;} END{for (i=1; i<=FNR; i++) print a[i]}' WWTree_WJets2500*el.dat >> WJets2500_el.dat
echo "WJets400_el" > WJets400_el.dat
awk '{a[FNR]+=$1;} END{for (i=1; i<=FNR; i++) print a[i]}' WWTree_WJets400*el.dat >> WJets400_el.dat
echo "WJets600_el" > WJets600_el.dat
awk '{a[FNR]+=$1;} END{for (i=1; i<=FNR; i++) print a[i]}' WWTree_WJets600*el.dat >> WJets600_el.dat
(echo WJets800ext1_el; cat WWTree_WJets800ext1_el.dat; )	> WJets800ext1_el.dat
(echo WJets_amcatnlo; cat WWTree_WJets_amcatnlo_el.dat; ) > WJets_amcatnlo_el.dat
(echo WJets_madgraph; cat WWTree_WJets_madgraph_el.dat; ) > WJets_madgraph_el.dat
echo "WW_excl_el" > WW_excl_el.dat
awk '{a[FNR]+=$1;} END{for (i=1; i<=FNR; i++) print a[i]}' WWTree_WW_excl*el.dat >> WW_excl_el.dat
echo "tch_el" > tch_el.dat
awk '{a[FNR]+=$1;} END{for (i=1; i<=FNR; i++) print a[i]}' WWTree_tch_*_el.dat >> tch_el.dat
echo "tch_bar_el" > tch_bar_el.dat
awk '{a[FNR]+=$1;} END{for (i=1; i<=FNR; i++) print a[i]}' WWTree_tch_bar_*_el.dat >> tch_bar_el.dat
echo "TTbar_amcatnlo_el" > TTbar_amcatnlo_el.dat
awk '{a[FNR]+=$1;} END{for (i=1; i<=FNR; i++) print a[i]}' WWTree_TTbar_amcatnlo_*_el.dat >> TTbar_amcatnlo_el.dat
echo "TTbar_powheg_el" > TTbar_powheg_el.dat
awk '{a[FNR]+=$1;} END{for (i=1; i<=FNR; i++) print a[i]}' WWTree_TTbar_powheg_*_el.dat >> TTbar_powheg_el.dat
(echo Signal_LL_el; cat WWTree_Signal_LL_el.dat; )	> Signal_LL_el.dat
(echo Signal_LT_el; cat WWTree_Signal_LT_el.dat; )	> Signal_LT_el.dat
(echo Signal_TT_el; cat WWTree_Signal_TT_el.dat; )	> Signal_TT_el.dat
(echo WW_excl_amcatnlo_el; cat WWTree_WW_excl_amcatnlo_el.dat; )	> WW_excl_amcatnlo_el.dat
echo "WZ_excl_amcatnlo_el" > WZ_excl_amcatnlo_el.dat
awk '{a[FNR]+=$1;} END{for (i=1; i<=FNR; i++) print a[i]}' WWTree_WZ_excl_amcatnlo_*_el.dat >> WZ_excl_amcatnlo_el.dat
(echo ZZ_excl_amcatnlo_el; cat WWTree_ZZ_excl_amcatnlo_el.dat; )	> ZZ_excl_amcatnlo_el.dat
(echo sch_el; cat WWTree_sch_el.dat; )	> sch_el.dat
(echo tWch_bar_ext1_el; cat WWTree_tWch_bar_ext1_el.dat; )	> tWch_bar_ext1_el.dat
(echo tWch_ext1_el; cat WWTree_tWch_ext1_el.dat; )	> tWch_ext1_el.dat

paste CutNames.txt Signal_LL_el.dat Signal_LT_el.dat Signal_TT_el.dat WJets_amcatnlo_el.dat WJets_madgraph_el.dat Wjet_HTbin_inclusive_el.dat TTbar_amcatnlo_el.dat TTbar_powheg_el.dat WW_excl_el.dat WW_excl_amcatnlo_el.dat WZ_excl_amcatnlo_el.dat ZZ_excl_amcatnlo_el.dat sch_el.dat tWch_ext1_el.dat tWch_bar_ext1_el.dat tch_bar_el.dat tch_el.dat > Ele_Event_Selection.dat 
awk 'BEGIN{print "<html>	\n<head>	\n<style>	\ntable, th, td { \n     border: 1px solid black; \n     border-collapse: collapse; \n} \n</style> \n</head> \n<body>	\n<h1>Cut Flow Table For Electron Channel</h1>	\n<table>"} {print "<tr>";for(i=1;i<=NF;i++)print "<td>" $i"</td>";print "</tr>"} END{print "</table>\n"}' Ele_Event_Selection.dat > CutFlowTable_${today}_Electrons.htm
awk 'BEGIN{OFS="\t&\t"}{print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11}' Ele_Event_Selection.dat >  CutFlowTable_${today}_Electrons.tex
##echo  "Adding all data for Electrons...."
##echo  "data_el_2016" > el_data_evt.dat
##awk '{a[FNR]+=$1;} END{for (i=1; i<=FNR; i++) print a[i]}'	el_data_el_2016_run*.dat	>>	el_data_evt.dat
##
#
#

echo "Wjet_HTbin_inclusive" > Wjet_HTbin_inclusive_mu.dat
awk '{a[FNR]+=$1;} END{for (i=1; i<=FNR; i++) print a[i]}' WWTree_WJets100*mu.dat WWTree_WJets200*mu.dat WWTree_WJets400*mu.dat WWTree_WJets600*mu.dat WWTree_WJets800ext1_mu.dat WWTree_WJets1200*mu.dat WWTree_WJets2500*mu.dat >> Wjet_HTbin_inclusive_mu.dat
echo "DYJetsToLL" > DYJetsToLL_mu.dat
awk '{a[FNR]+=$1;} END{for (i=1; i<=FNR; i++) print a[i]}' WWTree_DYJetsToLL_*mu.dat >> DYJetsToLL_mu.dat
echo "WJets100_mu" > WJets100_mu.dat
awk '{a[FNR]+=$1;} END{for (i=1; i<=FNR; i++) print a[i]}' WWTree_WJets100*mu.dat  >> WJets100_mu.dat
echo "WJets1200_mu" > WJets1200_mu.dat
awk '{a[FNR]+=$1;} END{for (i=1; i<=FNR; i++) print a[i]}' WWTree_WJets1200*mu.dat >> WJets1200_mu.dat
echo "WJets200_mu" > WJets200_mu.dat
awk '{a[FNR]+=$1;} END{for (i=1; i<=FNR; i++) print a[i]}' WWTree_WJets200*mu.dat >> WJets200_mu.dat
echo "WJets2500_mu" > WJets2500_mu.dat
awk '{a[FNR]+=$1;} END{for (i=1; i<=FNR; i++) print a[i]}' WWTree_WJets2500*mu.dat >> WJets2500_mu.dat
echo "WJets400_mu" > WJets400_mu.dat
awk '{a[FNR]+=$1;} END{for (i=1; i<=FNR; i++) print a[i]}' WWTree_WJets400*mu.dat >> WJets400_mu.dat
echo "WJets600_mu" > WJets600_mu.dat
awk '{a[FNR]+=$1;} END{for (i=1; i<=FNR; i++) print a[i]}' WWTree_WJets600*mu.dat >> WJets600_mu.dat
(echo WJets800ext1_mu; cat WWTree_WJets800ext1_mu.dat; )	> WJets800ext1_mu.dat
(echo WJets_amcatnlo; cat WWTree_WJets_amcatnlo_mu.dat; ) > WJets_amcatnlo_mu.dat
(echo WJets_madgraph; cat WWTree_WJets_madgraph_mu.dat; ) > WJets_madgraph_mu.dat
echo "WW_excl_mu" > WW_excl_mu.dat
awk '{a[FNR]+=$1;} END{for (i=1; i<=FNR; i++) print a[i]}' WWTree_WW_excl*mu.dat >> WW_excl_mu.dat
echo "tch_mu" > tch_mu.dat
awk '{a[FNR]+=$1;} END{for (i=1; i<=FNR; i++) print a[i]}' WWTree_tch_*_mu.dat >> tch_mu.dat
echo "tch_bar_mu" > tch_bar_mu.dat
awk '{a[FNR]+=$1;} END{for (i=1; i<=FNR; i++) print a[i]}' WWTree_tch_bar_*_mu.dat >> tch_bar_mu.dat
echo "TTbar_amcatnlo_mu" > TTbar_amcatnlo_mu.dat
awk '{a[FNR]+=$1;} END{for (i=1; i<=FNR; i++) print a[i]}' WWTree_TTbar_amcatnlo_*_mu.dat >> TTbar_amcatnlo_mu.dat
echo "TTbar_powheg_mu" > TTbar_powheg_mu.dat
awk '{a[FNR]+=$1;} END{for (i=1; i<=FNR; i++) print a[i]}' WWTree_TTbar_powheg_*_mu.dat >> TTbar_powheg_mu.dat
(echo Signal_LL_mu; cat WWTree_Signal_LL_mu.dat; )	> Signal_LL_mu.dat
(echo Signal_LT_mu; cat WWTree_Signal_LT_mu.dat; )	> Signal_LT_mu.dat
(echo Signal_TT_mu; cat WWTree_Signal_TT_mu.dat; )	> Signal_TT_mu.dat
(echo WW_excl_amcatnlo_mu; cat WWTree_WW_excl_amcatnlo_mu.dat; )	> WW_excl_amcatnlo_mu.dat
echo "WZ_excl_amcatnlo_mu" > WZ_excl_amcatnlo_mu.dat
awk '{a[FNR]+=$1;} END{for (i=1; i<=FNR; i++) print a[i]}' WWTree_WZ_excl_amcatnlo_*_mu.dat >> WZ_excl_amcatnlo_mu.dat
(echo ZZ_excl_amcatnlo_mu; cat WWTree_ZZ_excl_amcatnlo_mu.dat; )	> ZZ_excl_amcatnlo_mu.dat
(echo sch_mu; cat WWTree_sch_mu.dat; )	> sch_mu.dat
(echo tWch_bar_ext1_mu; cat WWTree_tWch_bar_ext1_mu.dat; )	> tWch_bar_ext1_mu.dat
(echo tWch_ext1_mu; cat WWTree_tWch_ext1_mu.dat; )	> tWch_ext1_mu.dat

paste CutNames.txt Signal_LL_mu.dat Signal_LT_mu.dat Signal_TT_mu.dat WJets_amcatnlo_mu.dat WJets_madgraph_mu.dat Wjet_HTbin_inclusive_mu.dat TTbar_amcatnlo_mu.dat TTbar_powheg_mu.dat WW_excl_mu.dat WW_excl_amcatnlo_mu.dat WZ_excl_amcatnlo_mu.dat ZZ_excl_amcatnlo_mu.dat sch_mu.dat tWch_ext1_mu.dat tWch_bar_ext1_mu.dat tch_bar_mu.dat tch_mu.dat > Mu_Event_Selection.dat 
awk 'BEGIN{print "<html>	\n<head>	\n<style>	\ntable, th, td { \n     border: 1px solid black; \n     border-collapse: collapse; \n} \n</style> \n</head> \n<body>	\n<h1>Cut Flow Table For Muon Channel</h1>	\n<table>"} {print "<tr>";for(i=1;i<=NF;i++)print "<td>" $i"</td>";print "</tr>"} END{print "</table>\n"}' Mu_Event_Selection.dat > CutFlowTable_${today}_Muons.htm
awk 'BEGIN{OFS="\t&\t"}{print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11}' Mu_Event_Selection.dat >  CutFlowTable_${today}_Muons.tex
##echo  "Adding all data for Electrons...."
##echo  "data_el_2016" > el_data_evt.dat
##awk '{a[FNR]+=$1;} END{for (i=1; i<=FNR; i++) print a[i]}'	el_data_el_2016_run*.dat	>>	el_data_evt.dat
##
#
#
