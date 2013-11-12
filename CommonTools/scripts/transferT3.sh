#!/bin/bash -f
# usage: ./transferT3.sh ../../HiggsAnalysisTools/cmst3_32X /castor/cern.ch/user/e/emanuele/CMST3/Higgs3.2.X

listdir=$1
cmst3dir=$2

echo "scanning $listdir"

\ls $listdir | grep ".list" | awk -F "." '{print $1}' | xargs -i rfmkdir "$cmst3dir/{}"

\ls $listdir | grep ".list" | awk -F "." '{print $1}' > datasets.txt

ls $listdir | grep list | xargs -i head -n 1 "$listdir/{}" | awk -F "default" '{print $1}' > castordirs.txt

N=0
while read LINE ; do
    N=$((N+1))
    echo "Processing $LINE"
    names[${N}]=$LINE
done < datasets.txt

N=0
while read LINE ; do
    N=$((N+1))
    castornames[${N}]=$LINE
done < castordirs.txt

rm -f createscripts.txt

mkdir scriptscopy

for ((i=1;i<$N+1;i++)); do
    echo ./copyRemaining.sh ${castornames[${i}]} $cmst3dir/${names[${i}]} scriptscopy/${names[${i}]}.csh >> scriptscopy/createscripts.txt
done

rm -f datasets.txt
rm -f castordirs.txt

source scriptscopy/createscripts.txt

\ls scriptscopy | grep csh | grep -v copyAll | awk -F "/" '{print $1}' | awk -F "." '{print "source scriptscopy/" $1 ".csh >&! scriptscopy/" $1 ".log &"}' > scriptscopy/copyAll.csh

source ./scriptscopy/copyAll.csh
