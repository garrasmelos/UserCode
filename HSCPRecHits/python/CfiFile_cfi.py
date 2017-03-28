import FWCore.ParameterSet.Config as cms

demo2 = cms.EDAnalyzer('HSCPRecHits',
							digisLabel = cms.InputTag("simMuonRPCDigis","RPCDigiSimLink"), #edm::DetSetVector<RPCDigiSimLink>     "simMuonRPCDigis"           "RPCDigiSimLink"   "DIGI2RAW"
							digis2Label = cms.InputTag("simMuonRPCDigis",""),
							recHitLabel = cms.InputTag("rpcRecHits",""),
							genEventInfo = cms.InputTag("generator"),
                     hepEventInfo = cms.InputTag("generator")
)	
