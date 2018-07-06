#!/bin/tcsh

foreach file ( `grep -E "Error|all events:        0" *.stdout | awk -F ":" '{print $1,$3}'` )
	echo "$file"
end
echo ""
echo "==========="
echo "	not found time line... "
echo "==========="
echo ""
foreach file (`grep -L "time to run this code" *.stdout | awk -F ":" '{print $1,$3}'`)
	echo $file
end
echo ""
