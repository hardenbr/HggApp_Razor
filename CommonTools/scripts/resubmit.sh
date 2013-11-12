#! /bin/sh

outputdir=$1

ls $outputdir/*Counters.root | awk -F _ '{print $13}' > done.txt
ls | awk -F _ '{print $2}' | awk -F . '{print $1}' > all.txt
cat all.txt > merged.txt
cat done.txt >> merged.txt
sort merged.txt > sorted.txt
chmod a+x *src
uniq -u sorted.txt | awk '{print "bsub -q 8nh -o ../log/resubmit_" $1 ".log source submit_" $1 ".src"}' > resubmit.sh

echo "done. Now run source resubmit.sh"
