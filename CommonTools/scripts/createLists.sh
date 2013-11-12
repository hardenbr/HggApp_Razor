#!/bin/bash -f
# usage: ./createLists.sh /castor/cern.ch/user/e/emanuele/Higgs32X/V02 cmst3_32X

maindir=$1
listdir=$2

echo "scanning $maindir"

rfdir $maindir | awk '{print $9}' > datasets.txt
rfdir $maindir | awk '{print "'"$maindir"'" "/" $9}' | xargs -i echo "rfdir " {} " | grep -v \" 0 \" | awk '{print \"{}/\" \$9}'" > commands.txt 


N=0
while read LINE ; do
    N=$((N+1))
    echo "Processing $LINE"
    names[${N}]=$LINE
done < datasets.txt


N=0
while read LINE ; do
    N=$((N+1))
    namescommand[${N}]=$LINE
    namesnum=${#namescommand}
done < commands.txt

rm -f finalcommand.sh

for ((i=1;i<$N+1;i++)); do
    echo ${namescommand[${i}]} ">" $listdir"/"${names[${i}]}".list" >> finalcommand.sh
done

echo "NOW reading from castor. It may take time..."

source finalcommand.sh

rm -f datasets.txt
rm -f commands.txt
rm -f finalcommand.sh

echo "LISTS are done in dir $listdir."
