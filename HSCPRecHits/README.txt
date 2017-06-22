from scratch: 

cmsrel CMSSW_8_0_21
cd CMSSW_8_0_21/src
git clone git://github.com/garrasmelos/UserCode HSCPAnalysis
scram b -j 16
cd HSCPAnalysis/HSCPRecHits/test/
cmsRun hscp_rpcrechits.py
