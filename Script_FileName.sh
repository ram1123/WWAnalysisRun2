tmp=$(echo $PWD | sed 's/eos\/uscms\///')
echo $tmp
for file in *; do
	find ${file} -name "*.root" | awk '{print "root://xrootd.unl.edu/'${tmp}'/"$1}' > ${file}.txt
done 
