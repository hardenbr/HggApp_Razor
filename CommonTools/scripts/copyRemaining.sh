#!/bin/bash -f
# ex.: ./copyRemaining.csh /castor/cern.ch/user/e/emanuele/VecBos21X/zeePYTHIA/ /castor/cern.ch/user/e/emanuele/CMST3/Vecbos2.1.X/ZeePYTHIA zeePYTHIA.csh

castordir=$1
cmst3dir=$2
script=$3

echo "CASTORDIR = $castordir"
echo "CMST3DIR = $cmst3dir"

rfdir $cmst3dir | awk '{print $9}' > filesinCMST3.txt
rfdir $castordir | awk '{print $9}' > filesinCASTOR.txt

cat filesinCMST3.txt > merged.txt
cat filesinCASTOR.txt >> merged.txt

sort merged.txt > merged_sorted.txt
uniq -u merged_sorted.txt > missing.txt

echo "#!/bin/csh -f" > $3
echo "setenv STAGE_HOST castorcms" >> $3
echo "setenv STAGE_SVCCLASS cmst3" >> $3

cat missing.txt | awk -v castordir="$castordir" -v cmst3dir="$cmst3dir" '{print "rfcp " castordir $1 " " cmst3dir}' >> $3

rm -f filesinCMST3.txt filesinCASTOR.txt merged.txt merged_sorted.txt missing.txt 

echo "Done. Run the script: source " $3 " now."

