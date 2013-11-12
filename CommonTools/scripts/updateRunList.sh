#!/bin/bash -f

echo "updating the 7 TeV data runlist in cmst3_35X/Data7TeV/data_minbias.list:"
mkdir tmp_runlist
CommonTools/scripts/createLists.sh /castor/cern.ch/user/m/meridian/VecBos/MinBias-v3 tmp_runlist

cd tmp_runlist
mkdir finalList
touch finalList/data_minbias.list
ls | grep "list" | xargs -i cat {} >> finalList/data_minbias.list
diff finalList/data_minbias.list ../cmst3_35X/Data7TeV/data_minbias.list > finalList/runsToProcess.list
mv  finalList/data_minbias.list ../cmst3_35X/Data7TeV/data_minbias.list
mv  finalList/runsToProcess.list ../cmst3_35X/Data7TeV
cd ..
#rm -r tmp_runlist
echo "DONE. The updated run list is: cmst3_35X/Data7TeV/data_minbias.list"
echo "and the new runs to process are in: cmst3_35X/Data7TeV/runsToProcess.list"
