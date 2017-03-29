from scratch: 

cmsrel CMSSW_7_4_14
cd CMSSW_7_4_14/src
git clone git://github.com/garrasmelos/UserCode test
scram b -j 16
cd test/HSCPRecHits/test/
cmsRun hscp_rpcrechits.py

To produce the plots:

cd test/HSCPRecHits/test/plots
mkdir histograms images rootfiles
cp <your rootfile> rootfile/.
sed -i -e 's/m1599_Gate1_BX_ToF/<your root file name w/o '.root'>/g' plotHistosRecHits.h
root -q -l -b plotHistosRecHits.cc

