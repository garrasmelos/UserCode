from scratch: 

cmsrel CMSSW_7_4_14
cd CMSSW_7_4_14/src
git clone git://github.com/garrasmelos/UserCode test
scram b -j 16
cd test/HSCPRecHits/test/
cmsRun hscp_rpcrechits.py

