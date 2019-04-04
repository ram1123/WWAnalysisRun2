for files in WWTree_*_el.root
do
	xrdcp -f $files root://cmseos.fnal.gov//store/user/lnujj/WpWm_aQGC_Ntuples_Ram/SecondStepOutput_21March2017/output_el/
done	
#for files in *.stdout
#do
#	xrdcp -f $files root://cmseos.fnal.gov//store/user/lnujj/WpWm_aQGC_Ntuples_Ram/SecondStepOutput_21March2017/Job_Logs/
#done	
#for files in WWTree_*_mu.root
#do
#	xrdcp -f $files root://cmseos.fnal.gov//store/user/lnujj/WpWm_aQGC_Ntuples_Ram/SecondStepOutput_21March2017/output_mu/
#done	
