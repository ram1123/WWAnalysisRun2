#!/bin/tcsh

set files = `ls *.stdout`

#foreach file (*.stdout)
#	echo $file
#end

#foreach file (`grep "(0) all events:        0        :       -nan" *.stdout`)
#foreach file ( `grep -E "Error|all events:        0" *.stdout | awk -F ":" '{print $1,$3}'` )
foreach file ( $files )
	grep -E "Error|all events:        0" $file
end
echo ""
echo "==========="
echo ""
#foreach file (`grep -L "time to run this code" *.stdout | awk -F ":" '{print $1,$3}'`)
foreach file ( $files )
	grep -L "time to run this code" $file
end
