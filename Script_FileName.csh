set tmp=`echo $PWD | sed 's/eos\/uscms\///'`
foreach file ( * )
	find ${file} -name "*.root" | awk '{print "root://xrootd.unl.edu/'${tmp}'/"$1}' > ${file}.txt
end 
