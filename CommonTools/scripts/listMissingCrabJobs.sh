#! /bin/bash
# usage: listMissingCrabJobs.sh listfile numtotjob
# eg: listMissingCrabJobs.sh doubleelectron.list 325

file=$1
numjobs=$2

cat $file | awk -F "default_data" '{print $2}' | awk -F "_" '{print $2}' | sort > done
seq 1 $numjobs | sort > todo
diff todo done | grep "<" | awk '{print $2}' > resub
cat resub | tr '\n' ','
echo
rm done todo resub
