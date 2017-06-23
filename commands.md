
	grep "time to run this code =" *.stdout | awk '{print $10,$7/60}'

	ls | grep data_el | awk '{print "\""$1"\""","}'
	
	voms-proxy-init

	condor_submit runstep2condor.jdl

	grep -r --exclude=\*.{root,o,exe,swp,bcup} genGravMass *
