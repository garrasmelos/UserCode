from scratch: 

cmsrel CMSSW_8_0_21
cd CMSSW_8_0_21/src
git clone git://github.com/garrasmelos/UserCode test
scram b -j 16
cd test/HSCPRecHits/test/
cmsRun hscp_rpcrechits.py
